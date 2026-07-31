module MultiGaussianLens


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using SpecialFunctions


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..Constants


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export potential!
export deflection!
export jacobian!


# --------------------------------------------------------------------------------------------------
# Main functions
# --------------------------------------------------------------------------------------------------
function potential!(ψ::U, θx::S, θy::S, D_d::Real, θxc::T, θyc::T, mass::T, θs::T, nl::Int64) where {U<:Real, S<:Real, T<:Vector{<:Real}}
   ψ_up = ψ
   for k in 1:nl
      κs = (2.0 * CONST_G * mass[k] * MASS_SUN / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)
      κs = κs * θs[k]^2

      dx = (θx - θxc[k]) / θs[k]
      dy = (θy - θyc[k]) / θs[k]
      dr = sqrt(dx^2 + dy^2)

      ψ_up = ψ_up + κs * (2.0 * log(dr) - expinti(-0.5 * dr^2))
   end
   return ψ_up
end

function potential!(ψ::U, θx::S, θy::S, D_d::Real, θxc::T, θyc::T, mass::T, θs::T, nl::Int64) where {U<:ROA, S<:ROA, T<:Vector{<:Real}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      κs = (2.0 * CONST_G * mass[k] * MASS_SUN / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)
      κs = κs * θs[k]^2

      @inbounds for j in ax2
         @inbounds for i in ax1
            dx = (θx[i, j] - θxc[k]) / θs[k]
            dy = (θy[i, j] - θyc[k]) / θs[k]
            dr = sqrt(dx^2 + dy^2)

            ψ[i, j] = ψ[i, j] + κs * (2.0 * log(dr) - expinti(-0.5 * dr^2))
         end
      end
   end
end


function deflection!(ψx::U, ψy::U, θx::S, θy::S, D_d::Real, θxc::T, θyc::T, mass::T, θs::T, nl::Int64) where {U<:Real, S<:Real, T<:Vector{<:Real}}
   ψx_up = ψx
   ψy_up = ψy
   for k in 1:nl
      κs = (2.0 * CONST_G * mass[k] * MASS_SUN / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)
      κs = 2.0 * κs * θs[k]

      dx = (θx - θxc[k]) / θs[k]
      dy = (θy - θyc[k]) / θs[k]
      dr = sqrt(dx^2 + dy^2)

      ψx_up = ψx_up + κs * (1.0 - exp(-0.5 * dr^2)) * dx / dr^2
      ψy_up = ψy_up + κs * (1.0 - exp(-0.5 * dr^2)) * dy / dr^2
   end
   return ψx_up, ψy_up
end

function deflection!(ψx::U, ψy::U, θx::S, θy::S, D_d::Real, θxc::T, θyc::T, mass::T, θs::T, nl::Int64) where {U<:ROA, S<:ROA, T<:Vector{<:Real}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      κs = (2.0 * CONST_G * mass[k] * MASS_SUN / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)
      κs = 2.0 * κs * θs[k]

      @inbounds for j in ax2
         @inbounds for i in ax1
            dx = (θx[i, j] - θxc[k]) / θs[k]
            dy = (θy[i, j] - θyc[k]) / θs[k]
            dr = sqrt(dx^2 + dy^2)

            ψx[i, j] = ψx[i, j] + κs * (1.0 - exp(-0.5 * dr^2)) * dx / dr^2
            ψy[i, j] = ψy[i, j] + κs * (1.0 - exp(-0.5 * dr^2)) * dy / dr^2
         end
      end
   end
end


function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, D_d::Real, θxc::T, θyc::T, mass::T, θs::T, nl::Int64) where {U<:Real, S<:Real, T<:Vector{<:Real}}
   ψxx_up = ψxx
   ψyy_up = ψyy
   ψxy_up = ψxy
   for k in 1:nl
      κs = (2.0 * CONST_G * mass[k] * MASS_SUN / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)
      
      dx = (θx - θxc[k]) / θs[k]
      dy = (θy - θyc[k]) / θs[k]
      dr = sqrt(dx^2 + dy^2)

      exp_term = exp(-0.5 * dr^2)
      α_r = κs * 2.0 * (1.0 - exp_term) / dr
      κ_r = κs * exp_term

      ψxx_up = ψxx_up + 2.0 * κ_r * dx^2 / dr^2 - α_r * (dx^2 - dy^2) / dr^3
      ψyy_up = ψyy_up + 2.0 * κ_r * dy^2 / dr^2 + α_r * (dx^2 - dy^2) / dr^3
      ψxy_up = ψxy_up + 2.0 * (κ_r - α_r / dr) * dx * dy / dr^2
   end
   return ψxx_up, ψyy_up, ψxy_up
end


function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, D_d::Real, θxc::T, θyc::T, mass::T, θs::T, nl::Int64) where {U<:ROA, S<:ROA, T<:Vector{<:Real}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      κs = (2.0 * CONST_G * mass[k] * MASS_SUN / CONST_C^2) / (D_d * θs[k]^2 * ANGLE_ARCSEC^2)
      
      @inbounds for j in ax2
         @inbounds for i in ax1
            dx = (θx[i, j] - θxc[k]) / θs[k]
            dy = (θy[i, j] - θyc[k]) / θs[k]
            dr = sqrt(dx^2 + dy^2)

            exp_term = exp(-0.5 * dr^2)
            α_r = κs * 2.0 * (1.0 - exp_term) / dr
            κ_r = κs * exp_term

            ψxx[i, j] = ψxx[i, j] + 2.0 * κ_r * dx^2 / dr^2 - α_r * (dx^2 - dy^2) / dr^3
            ψyy[i, j] = ψyy[i, j] + 2.0 * κ_r * dy^2 / dr^2 + α_r * (dx^2 - dy^2) / dr^3
            ψxy[i, j] = ψxy[i, j] + 2.0 * (κ_r - α_r / dr) * dx * dy / dr^2
         end
      end
   end
end

end