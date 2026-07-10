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
export logL_imageplane_fast
export logL_imageplane_fast_flux
export logL_imageplane_fast_timedelay



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

# Image-plane Newton solver settings (used by logL_imageplane_fast).
# The lens equation β = θ − a_dis·α(θ) is solved for θ near each observed image.
const IMG_SOLVE_MAX_ITER = 50        # Maximum Newton iterations per image
const IMG_SOLVE_TOL      = 1.0e-8    # Convergence tol on |β − θ + a_dis·α| (arcsec); far
                                     # below mas astrometry, above the Float64 round-off floor
const IMG_DET_MIN        = 1.0e-12   # |det J| below which the Newton step is abandoned

# Finite penalty added per image when the Newton solve fails to converge (e.g. the
# trial image sits on/near a critical curve). Finite (rather than −Inf) so the NM
# simplex still sees a gradient; effectively a rejection for the MCMC stage.
const CHI2_PEN_NOCONV = 1.0e8


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
                                    σx::Vector{Float64}, 
                                    σy::Vector{Float64}, 
                                    θ::Vector{Float64}, 
                                    n::Int64)
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

# Solve the lens equation β_target = θ − a_dis·α(θ) for the image position θ using
# Newton's method, starting from the observed image position (θx0, θy0). The Jacobian
# of the map θ ↦ β is J = I − a_dis·ψ_ij, evaluated on the fly from the lens model.
#
# Returns (θx, θy, converged). When converged == false the last iterate is returned
# (the caller applies a finite penalty); this happens if |det J| collapses (image on a
# critical curve) or the iteration does not reach tol within max_iter.
@inline function _solve_image_position(lens::AbstractLens, adis_value::Float64,
                                       βx_target::Float64, βy_target::Float64,
                                       θx0::Float64, θy0::Float64)
   θx = θx0
   θy = θy0

   @inbounds for _ in 1:IMG_SOLVE_MAX_ITER
      # Residual of the lens equation at the current guess
      ψx, ψy = Lenses.get_deflection(lens, θx, θy)
      fx = θx - adis_value * ψx - βx_target
      fy = θy - adis_value * ψy - βy_target

      # Converged?
      if abs(fx) < IMG_SOLVE_TOL && abs(fy) < IMG_SOLVE_TOL
         return θx, θy, true
      end

      # Jacobian J = I − a_dis·ψ_ij (symmetric)
      ψxx, ψyy, ψxy = Lenses.get_jacobian(lens, θx, θy)
      J11 = 1.0 - adis_value * ψxx
      J22 = 1.0 - adis_value * ψyy
      J12 = -adis_value * ψxy

      det = J11 * J22 - J12 * J12
      if abs(det) < IMG_DET_MIN
         return θx, θy, false
      end

      # Newton step Δθ = J⁻¹ f
      Δx = ( J22 * fx - J12 * fy) / det
      Δy = (-J12 * fx + J11 * fy) / det
      θx = θx - Δx
      θy = θy - Δy

      if !isfinite(θx) || !isfinite(θy)
         return θx0, θy0, false
      end
   end
   return θx, θy, false
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
            # Skip unmeasured images (marked with sigma <= 0)
            if σm[i] <= 0.0
               dm[i] = 0.0   # Placeholder, never used
               continue
            end

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
            # Skip unmeasured images
            if σm[i] <= 0.0
               continue
            end
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
            # Skip unmeasured images (marked with sigma <= 0)
            if σ_Δt[i] <= 0.0
               td_fit_approx[i] = 0.0   # Placeholder, never used
               continue
            end

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
            # Skip unmeasured images
            if σ_Δt[i] <= 0.0
               continue
            end
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
function logL_imageplane_fast(model::ModelConfig, 
                              lens::AbstractLens,
                              adis::Vector{Float64}, 
                              αx_all::Vector{Vector{Float64}}, 
                              αy_all::Vector{Vector{Float64}}, 
                              A_all::Vector{NTuple{4, Vector{Float64}}})
   # Initialize chi2 for position
   χ2_total = 0.0

   # Identity tuple
   I4 = (1.0, 0.0, 0.0, 1.0)

   # Storage for the per-knot model source positions (kid-indexed)
   n_knots = sum(length(src.knots) for src in model.source_config.sources)
   β_mod_s = Vector{NTuple{2, Float64}}(undef, n_knots)

   # Storage for the recovered image positions (kid-indexed), reused by flux / time-delay
   θ_mod_s = Vector{Vector{NTuple{2, Float64}}}(undef, n_knots)

   # Track whether every observed image was recovered.
   all_converged = true

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

         # Deflection vector at the observed image positions
         αx = @. adis_value * αx_all[kid]
         αy = @. adis_value * αy_all[kid]

         # Deformation tensor A = I − a_dis·ψ_ij at the observed image positions
         A = @. adis_value * A_all[kid]
         for i in eachindex(A)
            @. A[i] = I4[i] - A[i]
         end

         # Individual source positions (paired)
         β_ind = @. tuple(x - αx, y - αy)

         # Get weighted model source position (same as source-plane scheme) and the
         # per-image inverse covariance matrices S⁻¹ used for the image-plane χ².
         β_model, _, iS_all = _weighted_position(β_ind, A, σx, σy, σθ, n)
         βx_model, βy_model = β_model

         # Store model source position for flux / time-delay chi2 (kid-indexed)
         β_mod_s[kid] = β_model

         # Map the model source back to the image plane and accumulate the χ²
         θ_knot = Vector{NTuple{2, Float64}}(undef, n)
         χ2_knot = 0.0
         for i in 1:n
            # Solve β_model = θ − a_dis·α(θ) starting from the observed image position
            θx_m, θy_m, converged = _solve_image_position(lens, adis_value, βx_model, βy_model, x[i], y[i])

            # Penalize images whose predicted position could not be recovered
            if !converged
               all_converged = false
               θ_knot[i] = (x[i], y[i])          # placeholder, unused when !all_converged
               χ2_knot = χ2_knot + CHI2_PEN_NOCONV
               continue
            end

            # Store the recovered image position for flux / time-delay reuse
            θ_knot[i] = (θx_m, θy_m)

            # Image-plane residual δθ = θ_obs − θ_model
            δθx = x[i] - θx_m
            δθy = y[i] - θy_m

            # χ² = δθᵀ * S⁻¹ * δθ
            iS11, iS12, iS21, iS22 = iS_all[i]
            χ2_knot = χ2_knot + δθx * (iS11 * δθx + iS12 * δθy) + δθy * (iS21 * δθx + iS22 * δθy)
         end
         θ_mod_s[kid] = θ_knot
         χ2_total = χ2_total + χ2_knot

         kid = kid + 1
      end
      sid = sid + 1
   end
   return -0.5 * χ2_total, β_mod_s, θ_mod_s, all_converged
end


function logL_imageplane_fast_flux(model::ModelConfig,
                              lens::AbstractLens,
                              adis::Vector{Float64},
                              θ_mod_all::Vector{Vector{NTuple{2, Float64}}})
   # Initialize chi2 for flux
   χ2_total = 0.0

   sid = 1
   kid = 1
   for src in model.source_config.sources
      # Distance ratio for this source
      adis_value = adis[sid]

      for knot in src.knots
         # Knot magnitude values and errors
         m  = knot.m
         σm = knot.σm
         n  = length(m)

         # Skip knots without magnitude measurements (kid must still advance)
         if n == 0
            kid += 1
            continue
         end

         # Recovered image positions for this knot
         θ_mod = θ_mod_all[kid]

         # Best-fit (profiled) source magnitude from exact magnifications at θ_model
         m_src = 0.0
         wsum  = 0.0
         dm = Vector{Float64}(undef, n)
         for i in 1:n
            # Skip unmeasured images (marked with sigma <= 0)
            if σm[i] <= 0.0
               dm[i] = 0.0   # Placeholder, never used
               continue
            end

            # Exact signed determinant det[I − a_dis·ψ_ij] at the predicted image
            θx, θy = θ_mod[i]
            ψxx, ψyy, ψxy = Lenses.get_jacobian(lens, θx, θy)
            detA = (1.0 - adis_value * ψxx) * (1.0 - adis_value * ψyy) - (adis_value * ψxy)^2

            # Guard against an image sitting on a critical curve
            if abs(detA) < DETA_MIN
               χ2_total += CHI2_PEN_DETA
               detA = detA < 0.0 ? -DETA_MIN : DETA_MIN
            end

            # Signed magnification, floored against overflow
            abs_μ = max(abs(1.0 / detA), MU_MIN)

            dm[i] = m[i] + 2.5 * log10(abs_μ)
            w = 1.0 / σm[i]^2

            m_src = m_src + dm[i] * w
            wsum  = wsum + w
         end
         m_src = m_src / wsum

         # Chi2 for this knot, with the source magnitude m_src profiled out
         χ2_knot = 0.0
         for i in 1:n
            if σm[i] <= 0.0
               continue
            end
            χ2_knot = χ2_knot + (dm[i] - m_src)^2 / σm[i]^2
         end
         χ2_total = χ2_total + χ2_knot

         kid = kid + 1
      end
      sid = sid + 1
   end
   return -0.5 * χ2_total
end


function logL_imageplane_fast_timedelay(model::ModelConfig,
                                   lens::AbstractLens,
                                   adis::Vector{Float64},
                                   z_d::Float64,
                                   D_d::Float64,
                                   β_mod_all::Vector{NTuple{2, Float64}},
                                   θ_mod_all::Vector{Vector{NTuple{2, Float64}}})
   # Initialize chi2
   χ2_total = 0.0

   # Multiplicative constant (arrival time in days)
   constant_factor = (1.0 + z_d) * D_d * ANGLE_ARCSEC^2 / CONST_C / DAY2SECOND

   sid = 1
   kid = 1
   for src in model.source_config.sources
      # Distance ratio for this source
      adis_value = adis[sid]

      for knot in src.knots
         # Knot time delay values and errors
         Δt_obs = knot.td
         σ_Δt   = knot.σ_td
         n = length(Δt_obs)

         # Skip knots without time-delay measurements (kid must still advance)
         if n == 0
            kid = kid + 1
            continue
         end

         # Recovered image positions and model source position for this knot
         θ_mod = θ_mod_all[kid]
         βx_model, βy_model = β_mod_all[kid]
         pref = constant_factor / adis_value

         # Best-fit arrival-time zero-point
         Δt_0 = 0.0
         wsum = 0.0
         td_fit = Vector{Float64}(undef, n)
         for i in 1:n
            # Skip unmeasured images (marked with sigma <= 0)
            if σ_Δt[i] <= 0.0
               td_fit[i] = 0.0   # Placeholder, never used
               continue
            end

            # Exact Fermat potential at the predicted image position
            θx, θy = θ_mod[i]
            ψ = Lenses.get_potential(lens, θx, θy)
            fermat = 0.5 * ((θx - βx_model)^2 + (θy - βy_model)^2) - adis_value * ψ
            td_fit[i] = pref * fermat

            w = 1.0 / σ_Δt[i]^2
            Δt_0 = Δt_0 + (Δt_obs[i] - td_fit[i]) * w
            wsum = wsum + w
         end
         Δt_0 = Δt_0 / wsum

         # Chi2 for this knot, with the zero-point Δt_0 profiled out
         χ2_knot = 0.0
         for i in 1:n
            if σ_Δt[i] <= 0.0
               continue
            end
            χ2_knot = χ2_knot + (Δt_obs[i] - td_fit[i] - Δt_0)^2 / σ_Δt[i]^2
         end
         χ2_total = χ2_total + χ2_knot

         kid = kid + 1
      end
      sid = sid + 1
   end
   return -0.5 * χ2_total
end

end