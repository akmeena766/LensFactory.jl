module Likelihood

# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..LensModelIO

# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export chi2_position
export chi2_fluxratio
export chi2_timedelay
export chi2_parity

# --------------------------------------------------------------------------------------------------
# Functions
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

@inline function _weighted_position(βx::Vector{Float64}, βy::Vector{Float64}, A::NTuple{4, Vector{Float64}}, σx::Vector{Float64}, σy::Vector{Float64}, θ::Vector{Float64}, n::Int64)
   # Weight matrix: W = Σᵢ μᵢᵀ * Sᵢ⁻¹ * μᵢ
   sumW11 = 0.0
   sumW12 = 0.0
   sumW21 = 0.0
   sumW22 = 0.0

   # Individual weight matrices
   W_all = Vector{NTuple{4,Float64}}(undef, n)

   # Position vector: b = Σᵢ μᵢᵀ * Sᵢ⁻¹ * μᵢ * βᵢ
   b1 = 0.0
   b2 = 0.0
   
   @inbounds for i in 1:n
      # Cosine and sine of the angle
      cθ = cos(θ[i])
      sθ = sin(θ[i])

      # μ = A⁻¹
      μ11, μ12, μ21, μ22 = _inverse(A[1][i], A[2][i], A[3][i], A[4][i])

      # Inverse covariance matrix
      iS11, iS12, iS21, iS22 = _inverse_covariance(σx[i], σy[i], θ[i])

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
      sumW11 += W11
      sumW12 += W12
      sumW21 += W21
      sumW22 += W22

      # b = μᵀ * T * β
      b1 = b1 + (W11 * βx[i] + W21 * βy[i])
      b2 = b2 + (W12 * βx[i] + W22 * βy[i])
   end

   # Solve for the source position: β_model = W⁻¹ * b
   iW11, iW12, iW21, iW22 = _inverse(sumW11, sumW12, sumW21, sumW22)
   βx_model = iW11 * b1 + iW12 * b2
   βy_model = iW21 * b1 + iW22 * b2
   
   return βx_model, βy_model, W_all
end


function chi2_position(model::ModelConfig, adis::Vector{Float64}, αx_all::Vector{Vector{Float64}}, αy_all::Vector{Vector{Float64}}, A_all::Vector{NTuple{4, Vector{Float64}}})
   # Initialize chi2 for position
   chi2 = 0.0

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
         
         # Number of images for this knot
         n = length(x)
         
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
         βx_model, βy_model, W_all = _weighted_position(βx_ind, βy_ind, A, σx, σy, σθ, n)

         # Calculate knot log-likelihood
         χ² = 0.0
         for i in 1:n
            δβx = βx_ind[i] - βx_model
            δβy = βy_ind[i] - βy_model
            
            w11, w12, w21, w22 = W_all[i]

            # χ² = δβᵀ * W * δβ
            χ² = χ² + δβx * (w11 * δβx + w12 * δβy) + δβy * (w21 * δβx + w22 * δβy)
         end
         chi2 = chi2 - χ²
         
         kid = kid + 1
      end
      sid = sid + 1
   end
   return chi2   
end


function chi2_parity(model::ModelConfig, adis::Vector{Float64}, A_all::Vector{NTuple{4, Vector{Float64}}})
   # Initialize chi2 for parity
   chi2 = 0.0
   
   # Penalty for wrong parity
   penalty = model.source_config.parity_force

   # Identity tuple
   I4 = (1.0, 0.0, 0.0, 1.0)

   # Calculate log-likelihood for each source
   sid = 1
   kid = 1
   for src in model.source_config.sources
      # Distance ratio for this source
      adis_value = adis[sid]

      for knot in src.knots
         # Input knot image parities
         parity_input = knot.parity
         
         # Deformation tensor at the knot positions
         A = @. adis_value * A_all[kid]
         for i in eachindex(A)
            @. A[i] = I4[i] - A[i]
         end
         
         # Model parity of knot images
         parity_model = @. Int64(sign(A[1] * A[4] - A[2] * A[3]))

         # Parity log-likelihood
         for i in eachindex(parity_input)
            if parity_input[i] != parity_model[i]
               chi2 = chi2 - penalty
            end
         end   
         kid = kid + 1
      end      
      sid = sid + 1
   end
   return chi2
end


function chi2_fluxratio()
   return 0.0
end

function chi2_timedelay()
   return 0.0
end

end