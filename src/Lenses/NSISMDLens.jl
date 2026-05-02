module NSISMDLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!
export einstein_angle


"""
    potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: RV
   θE = 4.0 * pi * (vd * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   θs2 = θs^2

   dx = θx - θxc
   dy = θy - θyc
   dr = sqrt(θs2 + dx^2 + dy^2)

   ψ_up = ψ + θE * (dr - θs * log(dr + θs))
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
Calculate potential at given coordinates for NSISMD lens and update the potential (ψ) in place.
The lensing potential is given as,

```math
ψ(θ_x, θ_y) = 4 π \\left(\\frac{v_d}{{\\rm c}} \\right)^2 
   \\left[ \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2} - 
   θ_s \\, \\ln \\left(θ_s + \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2} \\right) \\right]. 
```

# Arguments
- `ψ`: Potential at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d`: Velocity dispersion (in ``\\rm \\mathbf{km/s}``).
- `θs`: Core radius of the lens (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `nothing`: Updates the potential (ψ) in place.
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
   θE = 4.0 * pi * (vd * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   θs2 = θs^2
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         dr = sqrt(θs2 + dx^2 + dy^2)

         ψ[i, j] = ψ[i, j] + θE * (dr - θs * log(dr + θs))
      end
   end
   return nothing
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: RV
   θE = 4.0 * pi * (vd * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   θs2 = θs^2

   dx = θx - θxc
   dy = θy - θyc
   dr = sqrt(θs2 + dx^2 + dy^2)

   ψx_up = ψx + θE * dx / (θs + dr)
   ψy_up = ψy + θE * dy / (θs + dr)
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
Calculate deflection at given coordinates for NSISMD lens and update the deflection components
(ψx, ψy) in place. The deflection angle components are given as,
```math
\\begin{align*} 
ψ_x(θ_x, θ_y) &= 4 π \\left(\\frac{v_d}{{\\rm c}} \\right)^2 \\frac{θ_x - θ_{x,c}}{θ_s + \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2}}, \\\\
ψ_y(θ_x, θ_y) &= 4 π \\left(\\frac{v_d}{{\\rm c}} \\right)^2 \\frac{θ_y - θ_{y,c}}{θ_s + \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2}}.
\\end{align*}
```

# Arguments
- `ψx`: x-component of the deflection at given coordinates
- `ψy`: y-component of the deflection at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d`: Velocity dispersion (in ``\\rm \\mathbf{km/s}``).
- `θs`: Core radius of the lens (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `nothing`: Updates the deflection (ψx, ψy) in place.
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
   θE = 4.0 * pi * (vd * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   θs2 = θs^2

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         dr = sqrt(θs2 + dx^2 + dy^2)
         ψx[i, j] = ψx[i, j] + θE * dx / (θs + dr)
         ψy[i, j] = ψy[i, j] + θE * dy / (θs + dr)
      end
   end
   return nothing
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: RV
   θE = 4.0 * pi * (vd * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   θs2 = θs^2

   dx = θx - θxc
   dy = θy - θyc
   dr = sqrt(θs2 + dx^2 + dy^2)
   inv_dr_dθ = 1.0 /  (θs + dr)

   ψxx_up = ψxx + θE * (inv_dr_dθ - (dx^2 / dr) * inv_dr_dθ^2)
   ψyy_up = ψyy + θE * (inv_dr_dθ - (dy^2 / dr) * inv_dr_dθ^2)
   ψxy_up = ψxy - θE * dx * dy / dr / (θs + dr)^2
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
Calculate jacobian at given coordinates for NSISMD lens and update the jacobian components 
(ψxx, ψyy, ψxy) in place. The jacobian components are given as,
```math
\\begin{align*}
ψ_{xx} &= 4 π \\left(\\frac{v_d}{{\\rm c}} \\right)^2 \\left[ \\frac{1}{θ_s + Θ} - \\frac{(θ_x - θ_{x,c})^2}{Θ \\, (θ_s + Θ)^2} \\right], \\\\
ψ_{yy} &= 4 π \\left(\\frac{v_d}{{\\rm c}} \\right)^2 \\left[ \\frac{1}{θ_s + Θ} - \\frac{(θ_y - θ_{y,c})^2}{Θ \\, (θ_s + Θ)^2} \\right], \\\\
ψ_{xy} &= 4 π \\left(\\frac{v_d}{{\\rm c}} \\right)^2 \\frac{- (θ_x - θ_{x,c}) (θ_y - θ_{y,c})}{Θ \\, (θ_s + Θ)^2},
\\end{align*}
```
where
```math
Θ = \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2}.
```

# Arguments
- `ψxx`: xx-component of the jacobian at given coordinates
- `ψyy`: yy-component of the jacobian at given coordinates
- `ψxy`: xy-component of the jacobian at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d`: Velocity dispersion (in ``\\rm \\mathbf{km/s}``).
- `θs`: Core radius of the lens (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `nothing`: Updates the jacobian (ψxx, ψyy, ψxy) in place.
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV) where T <: ROA
   θE = 4.0 * pi * (vd * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   θs2 = θs^2

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         dr = sqrt(θs2 + dx^2 + dy^2)
         inv_dr_dθ = 1.0 /  (θs + dr)

         ψxx[i, j] = ψxx[i, j] + θE * (inv_dr_dθ - (dx^2 / dr) * inv_dr_dθ^2)
         ψyy[i, j] = ψyy[i, j] + θE * (inv_dr_dθ - (dy^2 / dr) * inv_dr_dθ^2)
         ψxy[i, j] = ψxy[i, j] - θE * dx * dy / dr / (θs + dr)^2
      end
   end
   return nothing
end


"""
    einstein_angle(;D_ds::Float64=NaN, D_s::Float64=NaN, v_d::RV=NaN, x_s::Float64=NaN)
Calculate the Einstein angle for NSISMD lens,
```math
\\theta_E = \\sqrt{\\theta_E^2 - 2 \\, x_s \\, \\theta_E},
```
where,
```math
\\theta_E = 4 \\pi \\frac{D_{ds}}{D_s} \\left( \\frac{v_d}{{\\rm c}} \\right)^2.
```

# Keyword Arguments
- `D_ds = NaN`: ADD from the observer to the lens (in ``\\rm \\mathbf{meters}``).
- `D_s = NaN`: ADD from the observer to the source (in ``\\rm \\mathbf{meters}``).
- `v_d = NaN`: Velocity dispersion (in ``\\rm \\mathbf{km/s}``).
- `x_s = NaN`: Core radius of the lens (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `θE`: Einstein angle (in ``\\rm \\mathbf{arcseconds}``).
"""
function einstein_angle(; D_ds::Float64=NaN, D_s::Float64=NaN, v_d::RV=NaN, x_s::Float64=NaN)
   θE = 4π * (D_ds / D_s) * (v_d * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   return sqrt(θE^2 - 2.0 * x_s * θE)
end

end