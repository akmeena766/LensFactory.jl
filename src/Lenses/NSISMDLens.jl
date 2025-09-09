module NSISMDLens

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
ψ(\\pmb{θ}) = 4 π \\left(\\frac{v_d}{c} \\right)^2 
\\left[ \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2} - θ_s \\ln \\left( \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2} + θ_s \\right) \\right]
```
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: RV
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC

   dx = θx - θxc
   dy = θy - θyc
   θr = sqrt(θs^2 + dx^2 + dy^2)

   ψ_up = ψ + θE * (θr - θs * log(θr + θs))
   return ψ_up
end

function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
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
         ψ[i, j] = ψ[i, j] + θE * (θr - θs * log(θr + θs))
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
```math
\\pmb{\\hat{α}} (\\pmb{θ}) = 4 π \\left(\\frac{v_d}{c} \\right)^2 
                              \\frac{\\pmb{θ} - \\pmb{θ}_c}{θ_s + \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2}}
```
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: RV
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC

   dx = θx - θxc
   dy = θy - θyc
   θr = sqrt(θs^2 + dx^2 + dy^2)

   ψx_up = ψx + θE * dx / (θs + θr)
   ψy_up = ψy + θE * dy / (θs + θr)
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
         ψx[i, j] = ψx[i, j] + θE * dx / (θs + θr)
         ψy[i, j] = ψy[i, j] + θE * dy / (θs + θr)
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
```math
\\begin{align*}
ψ_{xx} (\\pmb{θ}) &= 4 π \\left(\\frac{v_d}{c} \\right)^2 
\\left[ \\frac{1}{θ_s + \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2}} - 
\\frac{(\\pmb{θ}_x - \\pmb{θ}_{xc})^2}{\\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2} \\left( θ_s + \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2} \\right)^2 } \\right] \\\\[5pt]
ψ_{yy} (\\pmb{θ}) &= 4 π \\left(\\frac{v_d}{c} \\right)^2 
\\left[ \\frac{1}{θ_s + \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2}} - 
\\frac{(\\pmb{θ}_y - \\pmb{θ}_{yc})^2}{\\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2} \\left( θ_s + \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2} \\right)^2 } \\right] \\\\[5pt]
ψ_{xy} (\\pmb{θ}) &= 4 π \\left(\\frac{v_d}{c} \\right)^2 
\\left[ \\frac{-(\\pmb{θ}_x - \\pmb{θ}_{xc}) \\, (\\pmb{θ}_y - \\pmb{θ}_{yc})}{\\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2} \\left( θ_s + \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2} \\right)^2 } \\right] \\\\[5pt]
\\end{align*}
```
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: RV
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC

   dx = θx - θxc
   dy = θy - θyc
   θr = sqrt(θs^2 + dx^2 + dy^2)

   ψxx_up = ψxx + θE * ((1.0 / (θs + θr)) - (dx^2 / θr / (θs + θr)^2))
   ψyy_up = ψyy + θE * ((1.0 / (θs + θr)) - (dy^2 / θr / (θs + θr)^2))
   ψxy_up = ψxy - θE * dx * dy / θr / (θs + θr)^2
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
         θr = sqrt(θs^2 + dx^2 + dy^2)
         ψxx[i, j] = ψxx[i, j] + θE * ( ( 1.0 / (θs + θr) ) - ( dx^2 / θr / (θs + θr)^2 ) )
         ψyy[i, j] = ψyy[i, j] + θE * ( ( 1.0 / (θs + θr) ) - ( dy^2 / θr / (θs + θr)^2 ) )
         ψxy[i, j] = ψxy[i, j] - θE * dx * dy / θr / (θs + θr)^2
      end
   end
end


"""
    einstein_angle(;D_ds::Float64=NaN, D_s::Float64=NaN, v_d::RV=NaN, x_s::Float64=NaN)
```math
θ_E = \\sqrt{ 4 π \\frac{D_{ds}}{D_s} \\left(\\frac{v_d}{c} \\right)^2 }
\\sqrt{4 π \\frac{D_{ds}}{D_s} \\left(\\frac{v_d}{c} \\right)^2 - 2\\,x_s}
```
"""
function einstein_angle(;D_ds::Float64=NaN, D_s::Float64=NaN, v_d::RV=NaN, x_s::Float64=NaN)
   θE = 4π * (D_ds / D_s) * (v_d / CONST_C)^2
   return sqrt(θE^2 - 2.0 * x_s * θE)
end

end