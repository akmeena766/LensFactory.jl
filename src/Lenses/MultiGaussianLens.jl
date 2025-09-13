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
      κ_s = (4.0 * CONST_G * mass[k] / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)

      dx = θx - θxc[k]
      dy = θy - θyc[k]
      dr = (dx^2 + dy^2) / θs[k]^2

      ψ_up = ψ_up + κ_s * θs[i]^2 * (log(sqrt(dr)) - 0.5 * expinti(-0.5 * dr))
   end
   return ψ_up
end

function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::S, θyc::S, mass::S, θs::S, nl::Int64) where {T <: ROA, S <: Vector{RV}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      κ_s = (4.0 * CONST_G * mass[k] / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)

      @inbounds for j in ax2
         @inbounds for i in ax1
            dx = θx[i, j] - θxc[k]
            dy = θy[i, j] - θyc[k]
            dr = (dx^2 + dy^2) / θs^2
            ψ[i, j] = ψ[i, j] + κ_s * θs[k]^2 * (log(sqrt(dr)) - 0.5 * expinti(-0.5 * dr))
         end
      end
   end
end


function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::S, θyc::S, mass::S, θs::S, nl::Int64) where {T <: RV, S <: Vector{RV}}
   ψx_up = ψx
   ψy_up = ψy

   for k in 1:nl
      κ_s = (4.0 * CONST_G * mass[k] / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)
      
      dx = θx - θxc[k]
      dy = θy - θyc[k]
      dr = (dx^2 + dy^2) / θs[k]^2

      ψx_up = ψx_up + κ_s * θs[i] * (1.0 - exp(-0.5 * dr)) * dx / dr
      ψy_up = ψy_up + κ_s * θs[i] * (1.0 - exp(-0.5 * dr)) * dy / dr
   end
   return ψx_up, ψy_up
end

function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::S, θyc::S, mass::S, θs::S, nl::Int64) where {T <: ROA, S <: Vector{RV}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      κ_s = (4.0 * CONST_G * mass[k] / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)
      
      @inbounds for j in ax2
         @inbounds for i in ax1
            dx = θx[i, j] - θxc[k]
            dy = θy[i, j] - θyc[k]
            dr = (dx^2 + dy^2) / θs[k]^2

            ψx[i, j] = ψx[i, j] + κ_s * θs[i] * (1.0 - exp(-0.5 * dr)) * dx / dr
            ψy[i, j] = ψy[i, j] + κ_s * θs[i] * (1.0 - exp(-0.5 * dr)) * dy / dr
         end
      end
   end
end


function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::S, θyc::S, mass::S, θs::S, nl::Int64) where {T <: RV, S <: Vector{RV}}
   ψxx_up = ψxx
   ψyy_up = ψyy
   ψxy_up = ψxy

   for k in 1:nl
      κ_s = (4.0 * CONST_G * mass[k] / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)
      
      dx = θx - θxc[k]
      dy = θy - θyc[k]
      dr = (dx^2 + dy^2) / θs[k]^2

      α_r = κ_s * (1.0 - exp(-0.5 * dr)) / sqrt(dr)
      κ_r = κ_s * exp(-0.5 * dr)

      ψxx_up = ψxx_up + κ_r * dx^2 / dr - α_r * (dx^2 - dy^2) / dr^(3/2)
      ψyy_up = ψyy_up + κ_r * dy^2 / dr + α_r * (dx^2 - dy^2) / dr^(3/2)
      ψxy_up = ψxy_up + (0.5 * κ_r - α_r / sqrt(dr)) * 2.0 * dx * dy / dr
   end
end

function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::S, θyc::S, mass::S, θs::S, nl::Int64) where {T <: ROA, S <: Vector{RV}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      κ_s = (4.0 * CONST_G * mass[k] / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)
      
      @inbounds for j in ax2
         @inbounds for i in ax1
            dx = θx[i, j] - θxc[k]
            dy = θy[i, j] - θyc[k]
            dr = (dx^2 + dy^2) / θs[k]^2

            α_r = κ_s * (1.0 - exp(-0.5 * dr)) / sqrt(dr)
            κ_r = κ_s * exp(-0.5 * dr)

            ψxx[i, j] = ψxx[i, j] + κ_r * dx^2 / dr - α_r * (dx^2 - dy^2) / dr^(3/2)
            ψyy[i, j] = ψyy[i, j] + κ_r * dy^2 / dr + α_r * (dx^2 - dy^2) / dr^(3/2)
            ψxy[i, j] = ψxy[i, j] + (0.5 * κ_r - α_r / sqrt(dr)) * 2.0 * dx * dy / dr
         end
      end
   end
end

end