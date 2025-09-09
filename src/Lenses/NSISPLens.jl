module NSISPLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!
export einstein_angle


"""
    potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
```math
ψ(\\pmb{θ}) = 4 π \\left(\\frac{v_d}{c} \\right)^2 \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2}
```
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: RV
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   
   dx = θx - θxc
   dy = θy - θyc

   ψ_up = ψ + θE * sqrt(θs^2 + dx^2 + dy^2)
   return ψ_up
end

function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   
   dx = 0.0
   dy = 0.0

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         ψ[i, j] = ψ[i, j] + θE * sqrt(θs^2 + dx^2 + dy^2)
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
```math
\\pmb{\\hat{α}} (\\pmb{θ}) = 4 π \\left(\\frac{v_d}{c} \\right)^2 
                              \\frac{\\pmb{θ} - \\pmb{θ}_c}{\\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2}}
```
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: RV
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC

   dx = θx - θxc
   dy = θy - θyc
   θr = sqrt(θs^2 + dx^2 + dy^2)

   ψx_up = ψx + θE * dx / θr
   ψy_up = ψy + θE * dy / θr
   return ψx_up, ψy_up
end

function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   
   dx = 0.0
   dy = 0.0
   θr = 0.0

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         θr = sqrt(θs^2 + dx^2 + dy^2)
         ψx[i, j] = ψx[i, j] + θE * dx / θr
         ψy[i, j] = ψy[i, j] + θE * dy / θr
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
```math
\\begin{align*}
ψ_{xx} (\\pmb{θ}) &= 4 π \\left(\\frac{v_d}{c} \\right)^2 
                     \\frac{θ_s^2 + (\\pmb{θ}_y - \\pmb{θ}_{yc})^2}{\\left( θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2 \\right)^{3/2}} \\\\[5pt]
ψ_{yy} (\\pmb{θ}) &= 4 π \\left(\\frac{v_d}{c} \\right)^2 
                     \\frac{θ_s^2 + (\\pmb{θ}_x - \\pmb{θ}_{xc})^2}{\\left( θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2 \\right)^{3/2}} \\\\[5pt]
ψ_{xy} (\\pmb{θ}) &= 4 π \\left(\\frac{v_d}{c} \\right)^2 
                     \\frac{- \\: (\\pmb{θ}_x - \\pmb{θ}_{xc}) \\: (\\pmb{θ}_y - \\pmb{θ}_{yc})}{\\left( θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2 \\right)^{3/2}}
\\end{align*}
```
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: RV
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   
   dx = θx - θxc
   dy = θy - θyc
   θr = (θs^2 + dx^2 + dy^2)^(3/2)

   ψxx_up = ψxx + θE * dy^2 / θr
   ψyy_up = ψyy + θE * dx^2 / θr
   ψxy_up = ψxy - θE * dx * dy / θr
   return ψxx_up, ψyy_up, ψxy_up
end

function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   
   dx = 0.0
   dy = 0.0
   θr = 0.0
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         θr = (θs^2 + dx^2 + dy^2)^(3/2)
         ψxx[i, j] = ψxx[i, j] + θE * dy^2 / θr
         ψyy[i, j] = ψyy[i, j] + θE * dx^2 / θr
         ψxy[i, j] = ψxy[i, j] - θE * dx * dy / θr
      end
   end
end


"""
    einstein_angle(;D_ds::Float64=NaN, D_s::Float64=NaN, v_d::RV=NaN, x_s::Float64=NaN)
```math
θ_E = \\sqrt{ 4 π \\frac{D_{ds}}{D_s} \\left(\\frac{v_d}{c} \\right)^2 - x_s^2 }
```
"""
function einstein_angle(;D_ds::Float64=NaN, D_s::Float64=NaN, v_d::RV=NaN, x_s::Float64=NaN)
   return sqrt((4π * (v_d / CONST_C)^2 * (D_ds/D_s))^2 - x_s^2)
end

end