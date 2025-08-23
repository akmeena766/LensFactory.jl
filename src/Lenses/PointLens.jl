module PointLens

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
function potential!(ψ::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where T <: ROA
   θE2::Float64 = 2.0 * CONST_G * mass / CONST_C^2 / Dol
   dx::Float64 = 0.0
   dy::Float64 = 0.0

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc + 1.0E-20
         dy = θy[i, j] - θyc + 1.0E-20
         ψ[i, j] = ψ[i, j] + θE2 * log(dx^2 + dy^2)
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where T <: ROA

```math
\\pmb{\\hat{α}} (\\pmb{θ}) = \\frac{4{\\rm G}M}{{\\rm c}^2} \\frac{1}{D_d} \\frac{\\pmb{θ} - \\pmb{θ}_c}{|\\pmb{θ} - \\pmb{θ}_c|^2}
```
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where T <: ROA
   θE2::Float64 = 4.0 * CONST_G * mass / CONST_C^2 / Dol
   
   dx::Float64 = 0.0
   dy::Float64 = 0.0
   θr::Float64 = 0.0

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc + 1.0E-20
         dy = θy[i, j] - θyc + 1.0E-20
         θr = dx^2 + dy^2
         ψx[i, j] = ψx[i, j] + θE2 * dx / θr
         ψy[i, j] = ψy[i, j] + θE2 * dy / θr
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where T <: ROA

```math
ψ_{xx} (\\pmb{θ}) = \\frac{4{\\rm G}M}{{\\rm c}^2} \\frac{1}{D_d} \\frac{(\\pmb{θ}_y - \\pmb{θ}_{yc})^2 - (\\pmb{θ}_x - \\pmb{θ}_{xc})^2}{|\\pmb{θ} - \\pmb{θ}_c|^4} \\\\
ψ_{yy} (\\pmb{θ}) = \\frac{4{\\rm G}M}{{\\rm c}^2} \\frac{1}{D_d} \\frac{(\\pmb{θ}_x - \\pmb{θ}_{xc})^2 - (\\pmb{θ}_y - \\pmb{θ}_{yc})^2}{|\\pmb{θ} - \\pmb{θ}_c|^4} \\\\
ψ_{xy} (\\pmb{θ}) = \\frac{4{\\rm G}M}{{\\rm c}^2} \\frac{1}{D_d} \\frac{-2 \\: (\\pmb{θ}_x - \\pmb{θ}_{xc}) \\: (\\pmb{θ}_y - \\pmb{θ}_{yc})}{|\\pmb{θ} - \\pmb{θ}_c|^4}
```
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, Dol::RV, θxc::RV, θyc::RV, mass::RV) where T <: ROA
   θE2::RV = 4.0 * CONST_G * mass / CONST_C^2 / Dol
   
   dx::Float64 = 0
   dy::Float64 = 0
   θr::Float64 = 0

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc + 1.0E-20
         dy = θy[i, j] - θyc + 1.0E-20
         θr = (dx^2 + dy^2)^2
         ψxx[i, j] = ψxx[i, j] - θE2 * (dx^2 - dy^2) / θr
         ψyy[i, j] = ψyy[i, j] + θE2 * (dx^2 - dy^2) / θr
         ψxy[i, j] = ψxy[i, j] - θE2 * 2.0 * dx * dy / θr
      end
   end
end

"""
    einstein_angle(Dol::RV, Dls::RV, Dos::RV, mass::RV)::RV

```math
θ_E = \\sqrt{\\frac{4{\\rm G} M}{c^2} \\frac{D_{ds}}{D_{d}D_{s}}}
```
"""
function einstein_angle(Dd::RV, Dds::RV, Ds::RV, mass::RV)::RV
  return  √( (4.0 * CONST_G * mass / CONST_C^2) * (Dds / Dd / Ds) )
end

end

