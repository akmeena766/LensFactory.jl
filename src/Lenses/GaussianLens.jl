module GaussianLens

# Inbuilt packages to use
using SpecialFunctions

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = dx^2 + dy^2

   ψ_up = ψ + κs * θs^2 * (2.0 * log(sqrt(dr)) - expinti(-0.5 * dr))
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)

   dx = 0.0
   dy = 0.0
   dr = 0.0

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = dx^2 + dy^2
         ψ[i, j] = ψ[i, j] + κs * θs^2 * (2.0 * log(sqrt(dr)) - expinti(-0.5 * dr))
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = dx^2 + dy^2

   ψx_up = ψx + 2.0 * κs * θs * (1.0 - exp(-0.5 * dr)) * dx / dr
   ψy_up = ψy + 2.0 * κs * θs * (1.0 - exp(-0.5 * dr)) * dy / dr
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)

   dx = 0.0
   dy = 0.0
   dr = 0.0

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = dx^2 + dy^2

         ψx[i, j] = ψx[i, j] + 2.0 * κs * θs * (1.0 - exp(-0.5 * dr)) * dx / dr
         ψy[i, j] = ψy[i, j] + 2.0 * κs * θs * (1.0 - exp(-0.5 * dr)) * dy / dr
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   
   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = dx^2 + dy^2

   α_r = 2.0 * κs * (1.0 - exp(-0.5 * dr)) / sqrt(dr)
   κ_r = κs * exp(-0.5 * dr)

   ψxx_up = ψxx + 2.0 * κ_r * dx^2 / dr - α_r * (dx^2 - dy^2) / dr^(3/2)
   ψyy_up = ψyy + 2.0 * κ_r * dy^2 / dr + α_r * (dx^2 - dy^2) / dr^(3/2)
   ψxy_up = ψxy + 2.0 * (κ_r - α_r / sqrt(dr)) * dx * dy / dr
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   
   dx = 0.0
   dy = 0.0
   dr = 0.0
   α_r = 0.0
   κ_r = 0.0  

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = dx^2 + dy^2

         α_r = 2.0 * κs * (1.0 - exp(-0.5 * dr)) / sqrt(dr)
         κ_r = κs * exp(-0.5 * dr)

         ψxx[i, j] = ψxx[i, j] + 2.0 * κ_r * dx^2 / dr - α_r * (dx^2 - dy^2) / dr^(3/2)
         ψyy[i, j] = ψyy[i, j] + 2.0 * κ_r * dy^2 / dr + α_r * (dx^2 - dy^2) / dr^(3/2)
         ψxy[i, j] = ψxy[i, j] + 2.0 * (κ_r - α_r / sqrt(dr)) * dx * dy / dr
      end
   end
end

end