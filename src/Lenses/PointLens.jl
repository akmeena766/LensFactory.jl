module PointLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!
export einstein_angle


"""
    potential!(ψ::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where {T <: ROA}

```math
\\psi(\\pmb{\\theta}) = \\frac{4{\\rm G}M}{{\\rm c}^2} \\frac{1}{D_d} \\ln |\\pmb{\\theta} - \\pmb{\\theta}_c|
```
"""
function potential!(ψ::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where {T <: ROA}
   θE2::Float64 = 2.0 * CONST_G * mass / CONST_C^2 / Dol
   
   ax2, ax1 = axes(θx, 1), axes(θx, 2)
   @inbounds for i in ax1
      @inbounds for j in ax2
         ψ[j, i] = ψ[j, i] + θE2 * log((θx[j, i] - θxc)^2 + (θy[j, i] - θyc)^2)
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where {T <: ROA}
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where {T <: ROA}
   θE2::Float64 = 4.0 * CONST_G * mass / CONST_C^2 / Dol
   θr::Float64 = 0

   ax2, ax1 = axes(θx, 1), axes(θx, 2)
   @inbounds for i in ax1
      @inbounds for j in ax2
         θr = (θx[j, i] - θxc)^2 + (θy[j, i] - θyc)^2
         ψx[j, i] = ψx[j, i] + θE2 * (θx[j, i] - θxc) / θr
         ψy[j, i] = ψy[j, i] + θE2 * (θy[j, i] - θyc) / θr
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where {T <: ROA}
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where {T <: ROA}
   θE2::RV = 4.0 * CONST_G * mass / CONST_C^2 / Dol
   θr::Float64 = 0
   θ1::Float64 = 0
   θ2::Float64 = 0

   ax2, ax1 = axes(θx, 1), axes(θx, 2)
   @inbounds for i in ax1
      @inbounds for j in ax2
         θ1 = θx[j, i] - θxc
         θ2 = θy[j, i] - θyc
         θr = (θ1^2 + θ2^2)^2
         ψxx[j, i] = ψxx[j, i] - θE2 * (θ1^2 - θ2^2) / θr
         ψyy[j, i] = ψyy[j, i] + θE2 * (θ1^2 - θ2^2) / θr
         ψxy[j, i] = ψxy[j, i] - θE2 * 2.0 * θ1 * θ2 / θr
      end
   end
end

"""
    einstein_angle(Dol::RV, Dls::RV, Dos::RV, mass::RV)::RV
"""
function einstein_angle(Dol::RV, Dls::RV, Dos::RV, mass::RV)::RV
  return  √( (4.0 * CONST_G * mass / CONST_C^2) * (Dls / Dol / Dos) )
end

end

