module GaussianLens

# Inbuilt packages to use
using SpecialFunctions

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!


function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   κ_s = (4.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)

   dx = θx - θxc
   dy = θy - θyc
   dr = (dx^2 + dy^2) / θs^2

   ψ_up = ψ + κ_s * θs^2 * (log(sqrt(dr)) - 0.5 * expinti(-0.5 * dr))
   return ψ_up
end

function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κ_s = (4.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)

   dx = 0.0
   dy = 0.0
   dr = 0.0

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         dr = (dx^2 + dy^2) / θs^2
         ψ[i, j] = ψ[i, j] + κ_s * θs^2 * (log(sqrt(dr)) - 0.5 * expinti(-0.5 * dr))
      end
   end
end


function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   κ_s = (4.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)

   dx = θx - θxc
   dy = θy - θyc
   dr = (dx^2 + dy^2) / θs^2

   ψx_up = ψx + κ_s * θs * (1.0 - exp(-0.5 * dr)) * dx / dr
   ψy_up = ψy + κ_s * θs * (1.0 - exp(-0.5 * dr)) * dy / dr
   return ψx_up, ψy_up
end

function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κ_s = (4.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)

   dx = 0.0
   dy = 0.0
   dr = 0.0

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         dr = (dx^2 + dy^2) / θs^2

         ψx[i, j] = ψx[i, j] + κ_s * θs * (1.0 - exp(-0.5 * dr)) * dx / dr
         ψy[i, j] = ψy[i, j] + κ_s * θs * (1.0 - exp(-0.5 * dr)) * dy / dr
      end
   end
end


function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   κ_s = (4.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   
   dx = θx - θxc
   dy = θy - θyc
   dr = (dx^2 + dy^2) / θs^2

   α_r = κ_s * (1.0 - exp(-0.5 * dr)) / sqrt(dr)
   κ_r = κ_s * exp(-0.5 * dr)

   ψxx_up = ψxx + κ_r * dx^2 / dr - α_r * (dx^2 - dy^2) / dr^(3/2)
   ψyy_up = ψyy + κ_r * dy^2 / dr + α_r * (dx^2 - dy^2) / dr^(3/2)
   ψxy_up = ψxy + (0.5 * κ_r - α_r / sqrt(dr)) * 2.0 * dx * dy / dr
   return ψxx_up, ψyy_up, ψxy_up
end

function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κ_s = (4.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   
   dx = 0.0
   dy = 0.0
   dr = 0.0
   α_r = 0.0
   κ_r = 0.0   

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         dr = (dx^2 + dy^2) / θs^2

         α_r = κ_s * (1.0 - exp(-0.5 * dr)) / sqrt(dr)
         κ_r = κ_s * exp(-0.5 * dr)

         ψxx[i, j] = ψxx[i, j] + κ_r * dx^2 / dr - α_r * (dx^2 - dy^2) / dr^(3/2)
         ψyy[i, j] = ψyy[i, j] + κ_r * dy^2 / dr + α_r * (dx^2 - dy^2) / dr^(3/2)
         ψxy[i, j] = ψxy[i, j] + (0.5 * κ_r - α_r / sqrt(dr)) * 2.0 * dx * dy / dr
      end
   end
end

end