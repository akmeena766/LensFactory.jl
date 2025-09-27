module MultiGaussianLens

# Inbuilt packages to use
using SpecialFunctions

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!


function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::S, θyc::S, mass::S, θs::S, nl::Int64) where {T <: RV, S <: Vector{RV}}
   ψ_up = ψ

   for k in 1:nl
      κs = (2.0 * CONST_G * mass[k] / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)

      dx = (θx - θxc[k]) / θs[k]
      dy = (θy - θyc[k]) / θs[k]
      dr = sqrt(dx^2 + dy^2)

      ψ_up = ψ_up + κs * θs[k]^2 * (2.0 * log(dr) - expinti(-0.5 * dr^2))
   end
   return ψ_up
end

function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::S, θyc::S, mass::S, θs::S, nl::Int64) where {T <: ROA, S <: Vector{RV}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      κs = (2.0 * CONST_G * mass[k] / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)

      @inbounds for j in ax2
         @inbounds for i in ax1
            dx = (θx[i, j] - θxc[k]) / θs[k]
            dy = (θy[i, j] - θyc[k]) / θs[k]
            dr = sqrt(dx^2 + dy^2)

            ψ[i, j] = ψ[i, j] + κs * θs[k]^2 * (2.0 * log(dr) - expinti(-0.5 * dr^2))
         end
      end
   end
end


function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::S, θyc::S, mass::S, θs::S, nl::Int64) where {T <: RV, S <: Vector{RV}}
   ψx_up = ψx
   ψy_up = ψy

   for k in 1:nl
      κs = (2.0 * CONST_G * mass[k] / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)
      
      dx = (θx[i, j] - θxc[k]) / θs[k]
      dy = (θy[i, j] - θyc[k]) / θs[k]
      dr = sqrt(dx^2 + dy^2)

      ψx_up = ψx_up + 2.0 * κs * θs[k] * (1.0 - exp(-0.5 * dr^2)) * dx / dr^2
      ψy_up = ψy_up + 2.0 * κs * θs[k] * (1.0 - exp(-0.5 * dr^2)) * dy / dr^2
   end
   return ψx_up, ψy_up
end

function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::S, θyc::S, mass::S, θs::S, nl::Int64) where {T <: ROA, S <: Vector{RV}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      κs = (2.0 * CONST_G * mass[k] / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)
      
      @inbounds for j in ax2
         @inbounds for i in ax1
            dx = (θx[i, j] - θxc[k]) / θs[k]
            dy = (θy[i, j] - θyc[k]) / θs[k]
            dr = sqrt(dx^2 + dy^2)

            ψx[i, j] = ψx[i, j] + 2.0 * κs * θs[k] * (1.0 - exp(-0.5 * dr^2)) * dx / dr^2
            ψy[i, j] = ψy[i, j] + 2.0 * κs * θs[k] * (1.0 - exp(-0.5 * dr^2)) * dy / dr^2
         end
      end
   end
end


function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::S, θyc::S, mass::S, θs::S, nl::Int64) where {T <: RV, S <: Vector{RV}}
   ψxx_up = ψxx
   ψyy_up = ψyy
   ψxy_up = ψxy

   for k in 1:nl
      κs = (2.0 * CONST_G * mass[k] / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)
      
      dx = (θx - θxc[k]) / θs[k]
      dy = (θy - θyc[k]) / θs[k]
      dr = sqrt(dx^2 + dy^2)

      α_r = κs * θs * 2.0 * (1.0 - exp(-0.5 * dr^2)) / dr
      κ_r = κs * exp(-0.5 * dr^2)

      ψxx_up = ψxx_up + 2.0 * κ_r * dx^2 / dr^2 - α_r * (dx^2 - dy^2) / dr^3
      ψyy_up = ψyy_up + 2.0 * κ_r * dy^2 / dr^2 + α_r * (dx^2 - dy^2) / dr^3
      ψxy_up = ψxy_up + 2.0 * (κ_r - α_r / dr) * dx * dy / dr^2
   end
end

function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::S, θyc::S, mass::S, θs::S, nl::Int64) where {T <: ROA, S <: Vector{RV}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      κs = (2.0 * CONST_G * mass[k] / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)
      
      @inbounds for j in ax2
         @inbounds for i in ax1
            dx = (θx[i, j] - θxc[k]) / θs[k]
            dy = (θy[i, j] - θyc[k]) / θs[k]
            dr = sqrt(dx^2 + dy^2)

            α_r = κs * θs * 2.0 * (1.0 - exp(-0.5 * dr^2)) / dr
            κ_r = κs * exp(-0.5 * dr^2)

            ψxx[i, j] = ψxx[i, j] + 2.0 * κ_r * dx^2 / dr^2 - α_r * (dx^2 - dy^2) / dr^3
            ψyy[i, j] = ψyy[i, j] + 2.0 * κ_r * dy^2 / dr^2 + α_r * (dx^2 - dy^2) / dr^3
            ψxy[i, j] = ψxy[i, j] + 2.0 * (κ_r - α_r / dr) * dx * dy / dr^2
         end
      end
   end
end

end