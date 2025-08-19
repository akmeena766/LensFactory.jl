module PlummerLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!
export einstein_angle


"""
    potential!(ψ::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where T <: ROA

```math
ψ(\\pmb{θ}) = \\frac{4{\\rm G}M}{{\\rm c}^2} \\frac{1}{D_d} \\ln |\\pmb{θ} - \\pmb{θ}_c|
```
"""
function potential!(Dol::RV, θxc::RV, θyc::RV, mass::RV, θs::RV, θx::T, θy::T, ψ::T) where {T<:ROA}
   θE2::Float64 = 2.0 * CONST_G * mass / CONST_C^2 / Dol

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψ[i, j] = θE2 * log(θs^2 + (θx[i, j] - θxc)^2 + (θy[i, j] - θyc)^2)
      end
   end
end

function deflection!(Dol::RV, θxc::RV, θyc::RV, mass::RV, θs::RV, θx::T, θy::T, ψx::T, ψy::T) where {T<:ROA}
   θE2::RV = 4.0 * CONST_G * mass / CONST_C^2 / Dol
   θr::Float64 = 0.0

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         θr = θs^2 + (θx[i, j] - θxc)^2 + (θy[i, j] - θyc)^2
         ψx[i, j] = θE2 * (θx[i, j] - θxc) / θr
         ψy[i, j] = θE2 * (θy[i, j] - θyc) / θr
      end
   end
end

function jacobian!(Dol::RV, θxc::RV, θyc::RV, mass::RV, θs::RV, θx::T, θy::T, ψxx::T, ψyy::T, ψxy::T) where {T<:ROA}
   θE2::RV = 4.0 * CONST_G * mass / CONST_C^2 / Dol
   θr::Float64 = 0.0
   θ1::Float64 = 0.0
   θ2::Float64 = 0.0

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         θ1 = θx[i, j] - θxc
         θ2 = θy[i, j] - θyc
         θr = (θs^2 + θ1^2 + θ2^2)^2
         ψxx[i, j] = +θE2 * (θs^2 - θ1^2 + θ2^2) / θr
         ψyy[i, j] = +θE2 * (θs^2 + θ1^2 - θ2^2) / θr
         ψxy[i, j] = -θE2 * 2.0 * θ1 * θ2 / θr
      end
   end
end

function einstein_angle(Dol::Real, Dls::Real, Dos::Real, mass::Real, θs::Real)::Real
   return sqrt((4.0 * CONST_G * mass / CONST_C^2) * (Dls / Dol / Dos) - θs^2)
end

end