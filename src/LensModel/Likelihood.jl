module Likelihood


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..LensModelIO
using ..Constants
using ..Lenses


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export logL_sourceplane
export logL_sourceplane_flux
export logL_sourceplane_timedelay
export logL_imageplane


# --------------------------------------------------------------------------------------------------
# Constants
# --------------------------------------------------------------------------------------------------
# Minimum allowed |det A| at a flux-measured image (|mu| <= 1e10). An image with
# a finite observed magnitude cannot lie on a critical curve.
const DETA_MIN = 1.0e-10

# Sloped finite penalty per violating image. Finite (rather than -Inf) so the
# NM simplex still sees a downhill direction (cf. glafic's chi2pen_* design);
# effectively a rejection for the MCMC stage.
const CHI2_PEN_DETA = 1.0e8

# Floor on |mu_mod| after the first-order correction, protecting log10 against
# overshoot of the linear extrapolation very close to critical curves.
const MU_MIN = 1.0e-6


# --------------------------------------------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------------------------------------------
@inline function _inverse(a11, a12, a21, a22)
   det = a11*a22 - a12*a21
   if det == 0.0
      @warn "Singular 2x2 matrix."
      return 0.0, 0.0, 0.0, 0.0
   end

   ia11 =  a22 / det
   ia12 = -a12 / det
   ia21 = -a21 / det
   ia22 =  a11 / det

   return ia11, ia12, ia21, ia22
end

# Inverse and rotated covariance matrix (θ from +x axis in CCW direction)
@inline function _inverse_covariance(σx, σy, θ)
   # Cosine and sine of the angle
   cθ = cos(θ)
   sθ = sin(θ)

   # S = Rᵀ * diag(σx², σy²) * R
   S11 = σx^2 * cθ^2 + σy^2 * sθ^2
   S12 = (σx^2 - σy^2) * sθ * cθ
   S21 = S12
   S22 = σx^2 * sθ^2 + σy^2 * cθ^2

   return _inverse(S11, S12, S21, S22)
end

@inline function _weighted_position(β_ind::Vector{NTuple{2, Float64}},
                                    A::NTuple{4, Vector{Float64}}, 
                                    σx::Vector{Float64}, σy::Vector{Float64}, 
                                    θ::Vector{Float64}, n::Int64)
   # Weight matrix: W = Σᵢ μᵢᵀ * Sᵢ⁻¹ * μᵢ
   sumW11 = 0.0
   sumW12 = 0.0
   sumW21 = 0.0
   sumW22 = 0.0

   # Individual weight matrices
   W_all = Vector{NTuple{4,Float64}}(undef, n)
   iS_all = Vector{NTuple{4,Float64}}(undef, n)

   # Position vector: b = Σᵢ μᵢᵀ * Sᵢ⁻¹ * μᵢ * βᵢ
   b1 = 0.0
   b2 = 0.0
   
   @inbounds for i in 1:n
      # Individual source position
      βx_i, βy_i = β_ind[i]

      # μ = A⁻¹
      μ11, μ12, μ21, μ22 = _inverse(A[1][i], A[2][i], A[3][i], A[4][i])

      # Inverse covariance matrix
      iS11, iS12, iS21, iS22 = _inverse_covariance(σx[i], σy[i], θ[i])
      iS_all[i] = (iS11, iS12, iS21, iS22)

      # T = S⁻¹ * μ
      T11 = iS11 * μ11 + iS12 * μ21
      T12 = iS11 * μ12 + iS12 * μ22
      T21 = iS21 * μ11 + iS22 * μ21
      T22 = iS21 * μ12 + iS22 * μ22

      # W = μᵀ * T 
      W11 = μ11 * T11 + μ21 * T21
      W12 = μ11 * T12 + μ21 * T22
      W21 = μ12 * T11 + μ22 * T21
      W22 = μ12 * T12 + μ22 * T22

      # Store individual weight matrix
      W_all[i] = (W11, W12, W21, W22)

      # Accumulate total W and Vector b
      sumW11 = sumW11 + W11
      sumW12 = sumW12 + W12
      sumW21 = sumW21 + W21
      sumW22 = sumW22 + W22

      # b = μᵀ * T * β
      b1 = b1 + (W11 * βx_i + W21 * βy_i)
      b2 = b2 + (W12 * βx_i + W22 * βy_i)
   end

   # Solve for the source position: β_model = W⁻¹ * b
   iW11, iW12, iW21, iW22 = _inverse(sumW11, sumW12, sumW21, sumW22)
   βx_model = iW11 * b1 + iW12 * b2
   βy_model = iW21 * b1 + iW22 * b2
   β_model = (βx_model, βy_model)
   
   return β_model, W_all, iS_all
end


# --------------------------------------------------------------------------------------------------
# Source plane likelihood functions
# --------------------------------------------------------------------------------------------------
function logL_sourceplane(model::ModelConfig, 
                          adis::Vector{Float64}, 
                          αx_all::Vector{Vector{Float64}}, 
                          αy_all::Vector{Vector{Float64}}, 
                          A_all::Vector{NTuple{4, Vector{Float64}}})
   # Initialize chi2 for position
   χ2_total = 0.0

   # Identity tuple
   I4 = (1.0, 0.0, 0.0, 1.0)

# Storage for the shared per-knot source positions (kid-indexed)
   n_knots = sum(length(src.knots) for src in model.source_config.sources)
   β_ind_s = Vector{Vector{NTuple{2, Float64}}}(undef, n_knots)
   β_mod_s = Vector{NTuple{2, Float64}}(undef, n_knots)

   # Calculate chi2 for position
   sid = 1
   kid = 1
   for src in model.source_config.sources
      # Distance ratio for this source
      adis_value = adis[sid]
      
      for knot in src.knots
         # Knot positions and measurement errors
         x  = knot.x
         y  = knot.y
         σx = knot.σx
         σy = knot.σy
         σθ = knot.σθ         
         n  = length(x)
         
         # Deflection vector at the knot positions
         αx = @. adis_value * αx_all[kid]
         αy = @. adis_value * αy_all[kid]

         # Deformation tensor at the knot positions
         A = @. adis_value * A_all[kid]
         for i in eachindex(A)
            @. A[i] = I4[i] - A[i]
         end

         # Individual source positions (paired)
         β_ind = @. tuple(x - αx, y - αy)

         # Get weighted source position (Section 4.1 in https://arxiv.org/pdf/astro-ph/0102340)
         β_model, W_all, _ = _weighted_position(β_ind, A, σx, σy, σθ, n)
         βx_model, βy_model = β_model

         # Store shared per-knot source positions for flux / time-delay chi2 (kid-indexed)
         β_ind_s[kid] = β_ind
         β_mod_s[kid] = β_model

         # Calculate knot χ2
         χ2_knot = 0.0
         for i in 1:n
            βx_i, βy_i = β_ind[i]
            δβx = βx_i - βx_model
            δβy = βy_i - βy_model
            
            # Read weight matrix components
            w11, w12, w21, w22 = W_all[i]

            # χ² = δβᵀ * W * δβ
            χ2_knot = χ2_knot + δβx * (w11 * δβx + w12 * δβy) + δβy * (w21 * δβx + w22 * δβy)
         end
         χ2_total = χ2_total + χ2_knot
         
         kid = kid + 1
      end
      sid = sid + 1
   end
   return -0.5 * χ2_total, β_ind_s, β_mod_s
end


function logL_sourceplane_flux(model::ModelConfig, 
                               adis::Vector{Float64}, 
                               β_ind_all::Vector{Vector{NTuple{2, Float64}}},
                               β_mod_all::Vector{NTuple{2, Float64}},
                               A_all::Vector{NTuple{4, Vector{Float64}}},
                               A_xp_all::Vector{NTuple{4, Vector{Float64}}},
                               A_xm_all::Vector{NTuple{4, Vector{Float64}}},
                               A_yp_all::Vector{NTuple{4, Vector{Float64}}},
                               A_ym_all::Vector{NTuple{4, Vector{Float64}}},
                               h_fd::Float64)
   # Initialize chi2 for parity
   χ2_total = 0.0

   # Identity tuple
   I4 = (1.0, 0.0, 0.0, 1.0)

   # Calculate chi2 for position
   sid = 1
   kid = 1
   for src in model.source_config.sources
      # Distance ratio for this source
      adis_value = adis[sid]
      
      for knot in src.knots         
         # Knot magnitude values and errors
         m  = knot.m
         σm = knot.σm
         n = length(m)

         # Skip knots with magnitude measurement (kid must still advance)
         if n == 0
            kid += 1
            continue
         end

         # Precomputed source positions for this knot
         β_ind = β_ind_all[kid]
         βx_model, βy_model = β_mod_all[kid]

         # Deformation tensor A = I − adis * ψ_ij
         A = @. adis_value * A_all[kid]
         for i in eachindex(A)
            @. A[i] = I4[i] - A[i]
         end
         detA = @. A[1] * A[4] - A[2] * A[3]

         # Signed determinants at the four shifted positions,
         # det A(θ) = det[I − adis·ψ_ij(θ)], assembled per source from raw tensors
         detA_xp = @. (1.0 - adis_value * A_xp_all[kid][1]) * (1.0 - adis_value * A_xp_all[kid][4]) - adis_value^2 * A_xp_all[kid][2] * A_xp_all[kid][3]
         detA_xm = @. (1.0 - adis_value * A_xm_all[kid][1]) * (1.0 - adis_value * A_xm_all[kid][4]) - adis_value^2 * A_xm_all[kid][2] * A_xm_all[kid][3]
         detA_yp = @. (1.0 - adis_value * A_yp_all[kid][1]) * (1.0 - adis_value * A_yp_all[kid][4]) - adis_value^2 * A_yp_all[kid][2] * A_yp_all[kid][3]
         detA_ym = @. (1.0 - adis_value * A_ym_all[kid][1]) * (1.0 - adis_value * A_ym_all[kid][4]) - adis_value^2 * A_ym_all[kid][2] * A_ym_all[kid][3]

         # Best-fit source magnitude (profiled) and chi2, with first-order
         # magnification correction:
         #   dμ/dx = −(1/detA²)·d(detA)/dx   (signed μ = 1/detA)
         m_src = 0.0
         wsum  = 0.0
         dm = Vector{Float64}(undef, n)
         for i in 1:n
            d0 = detA[i]

            if abs(d0) < DETA_MIN
               # Image on a critical curve: sloped finite penalty, then floor
               χ2_total += CHI2_PEN_DETA * (1.0 + log10(DETA_MIN / max(abs(d0), 1.0e-300)))
               if d0 < 0.0
                  d0 = -DETA_MIN
               else
                  d0 = DETA_MIN
               end
            end

            # Signed magnification at observed position
            μ0 = 1.0 / d0

            # Model-image displacement δθ = A⁻¹ δβ
            βx_i, βy_i = β_ind[i]
            δβx = βx_model - βx_i
            δβy = βy_model - βy_i
            δθx = ( A[4][i] * δβx - A[2][i] * δβy) / d0
            δθy = (-A[3][i] * δβx + A[1][i] * δβy) / d0

            # ∇μ from central differences of detA
            ddet_dx = (detA_xp[i] - detA_xm[i]) / (2.0 * h_fd)
            ddet_dy = (detA_yp[i] - detA_ym[i]) / (2.0 * h_fd)
            dμ_dx = -ddet_dx / d0^2
            dμ_dy = -ddet_dy / d0^2

            # Signed magnification with first order correction
            μ_mod = μ0 + dμ_dx * δθx + dμ_dy * δθy

            # Guard against overshoot of the linear extrapolation
            abs_μ = max(abs(μ_mod), MU_MIN)

            dm[i] = m[i] + 2.5 * log10(abs_μ)
            w = 1.0 / σm[i]^2

            m_src = m_src + dm[i] * w
            wsum  = wsum + w
         end
         m_src = m_src / wsum

         # Chi2 for this knot, with the source magnitude m_src profiled out
         χ2_knot = 0.0
         for i in 1:n
            χ2_knot = χ2_knot + (dm[i] - m_src)^2 / σm[i]^2
         end
         χ2_total = χ2_total + χ2_knot

         kid = kid + 1
      end
      sid = sid + 1
   end
   return -0.5 * χ2_total
end


function logL_sourceplane_timedelay(model::ModelConfig, 
                                    adis::Vector{Float64}, 
                                    z_d::Float64, 
                                    D_d::Float64, 
                                    β_ind_all::Vector{Vector{NTuple{2, Float64}}},
                                    β_mod_all::Vector{NTuple{2, Float64}},
                                    ψ_all::Vector{Vector{Float64}},
                                    αx_all::Vector{Vector{Float64}}, 
                                    αy_all::Vector{Vector{Float64}})
   # Initialize chi2
   χ2_total = 0.0

   # Multipicative constant
   constant_factor = (1.0 + z_d) * D_d * ANGLE_ARCSEC^2 / CONST_C / DAY2SECOND

   # Calculate chi2
   sid = 1
   kid = 1
   for src in model.source_config.sources
      # Distance ratio for this source
      adis_value = adis[sid]
   
      for knot in src.knots
         # Knot time delay values and errors
         Δt_obs  = knot.td
         σ_Δt    = knot.σ_td
         n = length(Δt_obs)

         # Skip knots without time-delay measurements (kid must still advance)
         if n == 0
            kid = kid + 1
            continue
         end

         # Precomputed source positions for this knot
         β_ind = β_ind_all[kid]
         βx_model, βy_model = β_mod_all[kid]

         # Scaled deflections at the observed positions
         αx = @. adis_value * αx_all[kid]
         αy = @. adis_value * αy_all[kid]
         ψ  = ψ_all[kid]

         # Best-fit time delay value
         Δt_0 = 0.0
         wsum = 0.0
         td_fit_approx = Vector{Float64}(undef, n)
         pref = constant_factor / adis_value
         for i in 1:n
            βx_i, βy_i = β_ind[i]
            δβx = βx_model - βx_i
            δβy = βy_model - βy_i

            ϕ_i = 0.5 * (αx[i]^2 + αy[i]^2) - adis_value * ψ[i] - (αx[i] * δβx + αy[i] * δβy)
            td_fit_approx[i] = pref * ϕ_i

            w = 1.0 / σ_Δt[i]^2
            Δt_0 = Δt_0 + (Δt_obs[i] - td_fit_approx[i]) * w
            wsum = wsum + w
         end
         Δt_0 = Δt_0 / wsum

         # Chi2 for this knot, with the zero-point Δt_0 profiled out
         χ2_knot = 0.0
         for i in 1:n
            χ2_knot = χ2_knot + (Δt_obs[i] - td_fit_approx[i] - Δt_0)^2 / σ_Δt[i]^2
         end
         χ2_total = χ2_total + χ2_knot

         kid = kid + 1
      end
      sid = sid + 1
   end
   return -0.5 * χ2_total
end


# --------------------------------------------------------------------------------------------------
# Image plane likelihood functions
# --------------------------------------------------------------------------------------------------
function logL_imageplane(model::ModelConfig, 
                         adis::Vector{Float64}, 
                         αx_all::Vector{Vector{Float64}}, 
                         αy_all::Vector{Vector{Float64}}, 
                         A_all::Vector{NTuple{4, Vector{Float64}}})
   # Initialize chi2 for parity
   χ2_total = 0.0

   # Identity tuple
   I4 = (1.0, 0.0, 0.0, 1.0)

   # Calculate chi2 for position
   sid = 1
   kid = 1
   for src in model.source_config.sources
      # Distance ratio for this source
      adis_value = adis[sid]
      
      for knot in src.knots
         # Knot positions and measurement errors
         x  = knot.x
         y  = knot.y
         σx = knot.σx
         σy = knot.σy
         σθ = knot.σθ         
         n  = length(x)
         
         # Deflection vector at the knot positions
         αx = @. adis_value * αx_all[kid]
         αy = @. adis_value * αy_all[kid]

         # Deformation tensor at the knot positions
         A = @. adis_value * A_all[kid]
         for i in eachindex(A)
            @. A[i] = I4[i] - A[i]
         end

         # Individual source positions using broadcasting
         βx_ind = @. x - αx
         βy_ind = @. y - αy

         # Get weighted source position (Section 4.1 in https://arxiv.org/pdf/astro-ph/0102340)
         βx_model, βy_model, _, iS_all = _weighted_position(βx_ind, βy_ind, A, σx, σy, σθ, n)
         
         # Predict image positions
         # TODO: Calculate chi2 for position
         
         # Calculate knot χ2
         χ2_knot = 0.0
         for i in 1:n
            δθx = x[i] - θx_model[i]
            δθy = y[i] - θy_model[i]
            
            # Read inverse covariance matrix components
            iS11, iS12, iS21, iS22 = iS_all[i]

            # χ² = δθᵀ * S⁻¹ * δθ
            χ2_knot = χ2_knot + δθx * (iS11 * δθx + iS12 * δθy) + δθy * (iS21 * δθx + iS22 * δθy)
         end
         χ2_total = χ2_total + χ2_knot
         
         kid = kid + 1
      end
      sid = sid + 1
   end
   return -0.5 * χ2_total
end

end