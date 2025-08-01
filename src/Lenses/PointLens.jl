module PointLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!
export einstein_angle

function potential!(ψ::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where {T <: ROA}
   θE2::Float64 = 2.0 * CONST_G * mass / CONST_C^2 / Dol
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψ[i, j] = ψ[i, j] + θE2 * log((θx[i, j] - θxc)^2 + (θy[i, j] - θyc)^2)
      end
   end
end

function deflection!(ψx::T, ψy::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where {T <: ROA}
   θE2::Float64 = 4.0 * CONST_G * mass / CONST_C^2 / Dol
   θr::Float64 = 0

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         θr = (θx[i, j] - θxc)^2 + (θy[i, j] - θyc)^2
         ψx[i, j] = ψx[i, j] + θE2 * (θx[i, j] - θxc) / θr
         ψy[i, j] = ψy[i, j] + θE2 * (θy[i, j] - θyc) / θr
      end
   end
end

function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where {T <: ROA}
   θE2::RV = 4.0 * CONST_G * mass / CONST_C^2 / Dol
   θr::Float64 = 0
   θ1::Float64 = 0
   θ2::Float64 = 0

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         θ1 = θx[i,j] - θxc
         θ2 = θy[i,j] - θyc
         θr = (θ1^2 + θ2^2)^2
         ψxx[i, j] = ψxx[i, j] - θE2 * (θ1^2 - θ2^2) / θr
         ψyy[i, j] = ψyy[i, j] + θE2 * (θ1^2 - θ2^2) / θr
         ψxy[i, j] = ψxy[i, j] - θE2 * 2.0 * θ1 * θ2 / θr
      end
   end
end

function einstein_angle(Dol::RV, Dls::RV, Dos::RV, mass::RV)::RV
  return  √( (4.0 * CONST_G * mass / CONST_C^2) * (Dls / Dol / Dos) )
end

end

