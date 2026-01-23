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
@inline function _rotated_covariance(σx, σy, θ)
   cθ = cos(θ)
   sθ = sin(θ)

   a = σx^2 * cθ^2 + σy^2 * sθ^2
   b = (σx^2 - σy^2) * sθ * cθ
   c = σx^2 * sθ^2 + σy^2 * cθ^2

   return a, b, c
end

@inline function _inverse_covariance(a, b, c)
   det = a*c - b*b
   det <= 0.0 && error("Invalid covariance matrix")

   ia =  c / det
   ib = -b / det
   ic =  a / det

   return ia, ib, ic
end

function LogL_position(model::ModelConfig, adis::Vector{Float64}, αx_all::Vector{Vector{Float64}}, αy_all::Vector{Vector{Float64}}, A_all::Vector{Vector{NTuple{4,Float64}}})
   # Initialize log-likelihood
   logL = 0.0

   # Calculate log-likelihood
   sid = 1
   kid = 1
   for src in model.source_config.sources
      # Distance ratio for this source
      adis = adis[sid]
      
      for knot in src.knots
         # Knot positions and measurement errors
         x  = knot.x
         y  = knot.y
         σx = knot.σx
         σy = knot.σy
         
         # Deflection and deformation tensor at the knot positions
         αx = αx_all[kid]
         αy = αy_all[kid]
         A  = A_all[kid]

         # Number of images for this knot
         n = length(x)

         # Individual source positions using broadcasting
         βx_individual = @. x - adis * αx
         βy_individual = @. y - adis * αy

         # Get weighted source position
                  


         kid = kid + 1
      end
      sid = sid + 1
   end

   return logL
end

function LogL_fluxratio()
   return 0.0
end

function LogL_parity()
   return 0.0
end

function LogL_timedelay()
   return 0.0
end

end