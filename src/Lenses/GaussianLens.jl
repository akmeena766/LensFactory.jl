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
   κs = κs * θs^2

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr2 = dx^2 + dy^2

   ψ_up = ψ + κs * (log(dr2) - expinti(-0.5 * dr2))
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   κs = κs * θs^2

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds @simd for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr2 = dx^2 + dy^2

         ψ[i, j] = ψ[i, j] + κs * (log(dr2) - expinti(-0.5 * dr2))
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   κs = 2.0 * κs * θs

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr2 = dx^2 + dy^2

   ψx_up = ψx + κs * (1.0 - exp(-0.5 * dr2)) * dx / dr2
   ψy_up = ψy + κs * (1.0 - exp(-0.5 * dr2)) * dy / dr2
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   κs = 2.0 * κs * θs

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds @simd for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr2 = dx^2 + dy^2

         ψx[i, j] = ψx[i, j] + κs * (1.0 - exp(-0.5 * dr2)) * dx / dr2
         ψy[i, j] = ψy[i, j] + κs * (1.0 - exp(-0.5 * dr2)) * dy / dr2
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
   dx2 = dx^2
   dy2 = dy^2
   dr = sqrt(dx2 + dy2)
   dr2 = dr^2
   dr3 = dr * dr2

   α_r = κs * 2.0 *  (1.0 - exp(-0.5 * dr2)) / dr
   κ_r = κs * exp(-0.5 * dr2)

   ψxx_up = ψxx + 2.0 * κ_r * dx2 / dr2 - α_r * (dx2 - dy2) / dr3
   ψyy_up = ψyy + 2.0 * κ_r * dy2 / dr2 + α_r * (dx2 - dy2) / dr3
   ψxy_up = ψxy + 2.0 * (κ_r - α_r / dr) * dx * dy / dr2
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds @simd for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dx2 = dx^2
         dy2 = dy^2
         dr = sqrt(dx2 + dy2)
         dr2 = dr^2
         dr3 = dr * dr2

         α_r = κs * 2.0 * (1.0 - exp(-0.5 * dr2)) / dr
         κ_r = κs * exp(-0.5 * dr2)

         ψxx[i, j] = ψxx[i, j] + 2.0 * κ_r * dx2 / dr2 - α_r * (dx2 - dy2) / dr3
         ψyy[i, j] = ψyy[i, j] + 2.0 * κ_r * dy2 / dr2 + α_r * (dx2 - dy2) / dr3
         ψxy[i, j] = ψxy[i, j] + 2.0 * (κ_r - α_r / dr) * dx * dy / dr2
      end
   end
end

end