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
export LogL_position
export LogL_fluxratio
export LogL_timedelay
export LogL_parity

# --------------------------------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------------------------------
@inline function _inverse(a11, a12, a21, a22)
   det = a11*a22 - a12*a21
   det == 0.0 && error("Singular 2x2 matrix.")

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
      μ11, μ12, μ21, μ22 = _inverse(A[i][1], A[i][2], A[i][3], A[i][4])

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


function LogL_position(model::ModelConfig, adis::Vector{Float64}, αx_all::Vector{Vector{Float64}}, αy_all::Vector{Vector{Float64}}, A_all::Vector{NTuple{4, Vector{Float64}}})
   # Initialize log-likelihood
   logL = 0.0

   # Identity tuple
   I4 = (1.0, 0.0, 0.0, 1.0)

   # Calculate log-likelihood
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
         θ  = knot.θ
         
         # Number of images for this knot
         n = length(x)

         # Deflection vector at the knot positions
         αx = @. adis_value * αx_all[kid]
         αy = @. adis_value * αy_all[kid]

         # Deformation tensor at the knot positions
         for i in eachindex(A_all[kid])
            @. A_all[kid][i] = I4[i] - (adis_value * A_all[kid][i])
         end

         # Individual source positions using broadcasting
         βx_ind = @. x - αx
         βy_ind = @. y - αy

         # Get weighted source position (Section 4.1 in https://arxiv.org/pdf/astro-ph/0102340)
         βx_model, βy_model, W_all = _weighted_position(βx_ind, βy_ind, A_all[kid], σx, σy, θ, n)

         # Calculate knot log-likelihood
         χ² = 0.0
         for i in 1:n
            δβx = βx_ind[i] - βx_model
            δβy = βy_ind[i] - βy_model
            
            w11, w12, w21, w22 = W_all[i]

            # χ² = δβᵀ * W * δβ
            χ² = χ² + δβx * (w11 * δβx + w12 * δβy) + δβy * (w21 * δβx + w22 * δβy)
         end
         logL = logL - 0.5 * χ²
         
         # Knot increment
         kid = kid + 1
      end
      sid = sid + 1
   end
   error("Log(L): ",LogL)
   return logL
end


function LogL_parity()
   # Penalty for wrong parity
   
   
   return 0.0
end


function LogL_fluxratio()
   return 0.0
end

function LogL_timedelay()
   return 0.0
end

end