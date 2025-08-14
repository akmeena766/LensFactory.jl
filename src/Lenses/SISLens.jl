
module SISLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!
export einstein_angle


"""
    potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: ROA

```math
ψ(\\pmb{θ}) = 4 π \\left(\\frac{v_d}{c} \\right)^2 |\\pmb{θ} - \\pmb{θ}_c|
```
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: ROA
   θE::Float64 = 4.0 * pi * (vd / CONST_C)^2
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψ[i, j] = ψ[i, j] + θE * sqrt( (θx[i, j] - θxc)^2 + (θy[i, j] - θyc)^2 ) + 1.0E-15
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: ROA
   θE::Float64 = 4.0 * pi * (vd / CONST_C)^2 
   θr::Float64 = 0.0

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         θr = sqrt( (θx[i, j] - θxc)^2 + (θy[i, j] - θyc)^2 ) + 1.0E-15

         ψx[i, j] = ψx[i, j] + θE * (θx[i, j] - θxc) / θr
         ψy[i, j] = ψy[i, j] + θE * (θy[i, j] - θyc) / θr
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: ROA
   θE::Float64 = 4.0 * pi * (vd / CONST_C)^2 
   
   dx::Float64 = 0.0
   dy::Float64 = 0.0
   θr::Float64 = 0.0
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dy = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         θr = (dx^2 + dy^2)^(3/2) + 1.0E-15

         ψxx[i, j] = ψxx[i, j] + θE * dy^2 / θr
         ψyy[i, j] = ψyy[i, j] + θE * dx^2 / θr
         ψxy[i, j] = ψxy[i, j] - θE * dx * dy / θr
      end
   end
end


"""
    einstein_angle(Dds::RV, Ds::RV, vd::RV)::RV

```math
θ_E = 4 π \\frac{D_{ds}}{D_s} \\left(\\frac{v_d}{c} \\right)^2
```
"""
function einstein_angle(Dds::RV, Ds::RV, vd::RV)::RV
   return 4.0 * pi * ( vd / CONST_C )^2 * (Dds / Ds)
end

end