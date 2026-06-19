module MultiPlummerLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!


function potential!(ψ::Real, θx::S, θy::S, D_d::Real, θxc::T, θyc::T, mass::T, θs::T, nl::Int64) where {S<:Real, T<:Vector{<:Real}}
   ψ_up = ψ
   for k in 1:nl
      θE2 = (2.0 * CONST_G * mass[k] * MASS_SUN / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
      θs2 = θs[k]^2

      dx = (θx - θxc[k])
      dy = (θy - θyc[k])

      ψ_up = ψ_up + θE2 * log(θs2 + dx^2 + dy^2)
   end
   return ψ_up
end

function potential!(ψ::ROA, θx::S, θy::S, D_d::Real, θxc::T, θyc::T, mass::T, θs::T, nl::Int64) where {S<:ROA, T<:Vector{<:Real}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      θE2 = (2.0 * CONST_G * mass[k] * MASS_SUN / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
      θs2 = θs[k]^2

      @inbounds for j in ax2
         @inbounds for i in ax1
            dx = θx[i, j] - θxc[k]
            dy = θy[i, j] - θyc[k]

            ψ[i, j] = ψ[i, j] + θE2 * log(θs2 + dx^2 + dy^2)
         end
      end
   end
end

function deflection!(ψx::Real, ψy::Real, θx::S, θy::S, D_d::Real, θxc::T, θyc::T, mass::T, θs::T, nl::Int64) where {S<:Real, T<:Vector{<:Real}}
   ψx_up = ψx
   ψy_up = ψy
   for k in 1:nl
      θE2 = (4.0 * CONST_G * mass[k] * MASS_SUN / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
      θs2 = θs[k]^2

      dx = θx - θxc[k]
      dy = θy - θyc[k]
      θr = θs2 + dx^2 + dy^2

      ψx_up = ψx_up + θE2 * dx / θr
      ψy_up = ψy_up + θE2 * dy / θr
   end
   return ψx_up, ψy_up
end

function deflection!(ψx::ROA, ψy::ROA, θx::S, θy::S, D_d::Real, θxc::T, θyc::T, mass::T, θs::T, nl::Int64) where {S<:ROA, T<:Vector{<:Real}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      θE2 = (4.0 * CONST_G * mass[k] * MASS_SUN / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
      θs2 = θs[k]^2

      @inbounds for j in ax2
         @inbounds for i in ax1
            dx = θx[i, j] - θxc[k]
            dy = θy[i, j] - θyc[k]
            θr = θs2 + dx^2 + dy^2

            ψx[i, j] = ψx[i, j] + θE2 * dx / θr
            ψy[i, j] = ψy[i, j] + θE2 * dy / θr
         end
      end
   end
end

function jacobian!(ψxx::Real, ψyy::Real, ψxy::Real, θx::S, θy::S, D_d::Real, θxc::T, θyc::T, mass::T, θs::T, nl::Int64) where {S<:Real, T<:Vector{<:Real}}
   ψxx_up = ψxx
   ψyy_up = ψyy
   ψxy_up = ψxy
   for k in 1:nl
      θE2 = (4.0 * CONST_G * mass[k] * MASS_SUN / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
      θs2 = θs[k]^2

      dx = θx - θxc[k]
      dy = θy - θyc[k]
      θr = (θs2 + dx^2 + dy^2)^2

      ψxx_up = ψxx_up + θE2 * (θs2 - dx^2 + dy^2) / θr
      ψyy_up = ψyy_up + θE2 * (θs2 + dx^2 - dy^2) / θr
      ψxy_up = ψxy_up - θE2 * 2.0 * dx * dy / θr
   end
   return ψxx_up, ψyy_up, ψxy_up
end

function jacobian!(ψxx::ROA, ψyy::ROA, ψxy::ROA, θx::S, θy::S, D_d::Real, θxc::T, θyc::T, mass::T, θs::T, nl::Int64) where {S<:ROA, T<:Vector{<:Real}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      θE2 = (4.0 * CONST_G * mass[k] * MASS_SUN / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
      θs2 = θs[k]^2

      @inbounds for j in ax2
         @inbounds for i in ax1
            dx = θx[i, j] - θxc[k]
            dy = θy[i, j] - θyc[k]
            θr = (θs2 + dx^2 + dy^2)^2

            ψxx[i, j] = ψxx[i, j] + θE2 * (θs2 - dx^2 + dy^2) / θr
            ψyy[i, j] = ψyy[i, j] + θE2 * (θs2 + dx^2 - dy^2) / θr
            ψxy[i, j] = ψxy[i, j] - θE2 * 2.0 * dx * dy / θr
         end
      end
   end
end

end