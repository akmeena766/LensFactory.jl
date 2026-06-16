module PointLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!
export einstein_angle


"""
    potential!(ψ::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, mass::Real) where T <: Real
"""
function potential!(ψ::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, mass::Real) where T <: Real
   θE2 = (2.0 * CONST_G * mass * MASS_SUN / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
   
   dx = θx - θxc
   dy = θy - θyc

   ψ_up = ψ + θE2 * log(dx^2 + dy^2)
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, mass::Real) where T <: ROA
Calculate potential at given coordinates for a point mass lens and update the potential (ψ) in place.
The lensing potential is given as,
```math
\\psi(θ_x, θ_y) = \\frac{4 \\, \\rm{G} \\, M}{\\rm{c}^2} \\frac{1}{D_d} \\ln |\\pmb{θ} - \\pmb{θ}_c|.
```

# Arguments
- `ψ`: Potential at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `D_d::Real`: ADD from the observer to the lens (in ``\\rm \\mathbf{meters}``).
- `θxc::Real`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc::Real`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).

# Returns
- `nothing`: Updates the potential (ψ) in place.
"""
function potential!(ψ::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, mass::Real) where T <: ROA
   θE2 = (2.0 * CONST_G * mass * MASS_SUN / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         ψ[i, j] = ψ[i, j] + θE2 * log(dx^2 + dy^2)
      end
   end
   return nothing
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, mass::Real) where T <: Real
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, mass::Real) where T <: Real
   θE2 = (4.0 * CONST_G * mass * MASS_SUN / CONST_C^2 / D_d) / ANGLE_ARCSEC^2

   dx = θx - θxc
   dy = θy - θyc
   dr2 = dx^2 + dy^2

   ψx_up = ψx + θE2 * dx / dr2
   ψy_up = ψy + θE2 * dy / dr2
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, mass::Real) where T <: ROA
Calculate deflection at given coordinates for a point mass lens and update the deflection components
(ψx, ψy) in place. The deflection angle components are given as,
```math
\\begin{align*} 
ψ_x(θ_x, θ_y) &= \\frac{4 \\, \\rm{G} \\, M}{\\rm{c}^2} \\frac{1}{D_d} \\frac{θ_x - θ_{x,c}}{|\\pmb{θ} - \\pmb{θ}_c|^2}, \\\\
ψ_y(θ_x, θ_y) &= \\frac{4 \\, \\rm{G} \\, M}{\\rm{c}^2} \\frac{1}{D_d} \\frac{θ_y - θ_{y,c}}{|\\pmb{θ} - \\pmb{θ}_c|^2}.
\\end{align*}
```

# Arguments
- `ψx`: x-component of the deflection at given coordinates
- `ψy`: y-component of the deflection at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `D_d::Real`: ADD from the observer to the lens (in ``\\rm \\mathbf{meters}``).
- `θxc::Real`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc::Real`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).

# Returns
- `nothing`: Updates the deflection (ψx, ψy) in place.
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, mass::Real) where T <: ROA
   θE2 = (4.0 * CONST_G * mass * MASS_SUN / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         dr2 = dx^2 + dy^2
         
         ψx[i, j] = ψx[i, j] + θE2 * dx / dr2
         ψy[i, j] = ψy[i, j] + θE2 * dy / dr2
      end
   end
   return nothing
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, mass::Real) where T <: Real
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, mass::Real) where T <: Real
   θE2 = (4.0 * CONST_G * mass * MASS_SUN / CONST_C^2 / D_d) / ANGLE_ARCSEC^2

   dx = θx - θxc
   dy = θy - θyc

   dx2 = dx^2
   dy2 = dy^2
   dr4 = (dx2 + dy2)^2

   ψxx_up = ψxx - θE2 * (dx2 - dy2) / dr4
   ψyy_up = ψyy + θE2 * (dx2 - dy2) / dr4
   ψxy_up = ψxy - θE2 * 2.0 * dx * dy / dr4
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, mass::Real) where T <: ROA
Calculate jacobian at given coordinates for a point mass lens and update the jacobian components 
(ψxx, ψyy, ψxy) in place. The jacobian components are given as,
```math
\\begin{align*}
\\psi_{xx} &= \\frac{4 \\, \\rm{G} \\, M}{\\rm{c}^2} \\frac{1}{D_d} \\frac{(θ_y - θ_{y,c})^2 - (θ_x - θ_{x,c})^2}{|\\pmb{θ} - \\pmb{θ}_c|^4}, \\\\
\\psi_{yy} &= \\frac{4 \\, \\rm{G} \\, M}{\\rm{c}^2} \\frac{1}{D_d} \\frac{(θ_x - θ_{x,c})^2 - (θ_y - θ_{y,c})^2}{|\\pmb{θ} - \\pmb{θ}_c|^4}, \\\\
\\psi_{xy} &= \\frac{4 \\, \\rm{G} \\, M}{\\rm{c}^2} \\frac{1}{D_d} \\frac{-2 \\, (θ_x - θ_{x,c}) (θ_y - θ_{y,c})}{|\\pmb{θ} - \\pmb{θ}_c|^4}.
\\end{align*}
```

# Arguments
- `ψxx`: x-component of the jacobian at given coordinates
- `ψyy`: y-component of the jacobian at given coordinates
- `ψxy`: xy-component of the jacobian at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `D_d::Real`: ADD from the observer to the lens (in ``\\rm \\mathbf{meters}``).
- `θxc::Real`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc::Real`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).

# Returns
- `nothing`: Updates the jacobian (ψxx, ψyy, ψxy) in place.
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, mass::Real) where T <: ROA
   θE2 = (4.0 * CONST_G * mass * MASS_SUN / CONST_C^2 / D_d) / ANGLE_ARCSEC^2

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         
         dx2 = dx^2
         dy2 = dy^2
         dr4 = (dx2 + dy2)^2

         ψxx[i, j] = ψxx[i, j] - θE2 * (dx2 - dy2) / dr4
         ψyy[i, j] = ψyy[i, j] + θE2 * (dx2 - dy2) / dr4
         ψxy[i, j] = ψxy[i, j] - θE2 * 2.0 * dx * dy / dr4
      end
   end
   return nothing
end

"""
    einstein_angle(;D_d::Float64=NaN, D_ds::Float64=NaN, D_s::Float64=NaN, mass::Float64=NaN)
Calculate the Einstein angle for a point mass lens,
```math
\\theta_E = \\sqrt{\\frac{4 \\, \\rm{G} \\, M}{\\rm{c}^2} \\frac{D_{ds}}{D_d D_s}}.
```

# Keyword Arguments
- `D_d::Float64`=NaN: ADD from observer to lens (in ``\\rm \\mathbf{meters}``).
- `D_ds::Float64`=NaN: ADD from lens to source (in ``\\rm \\mathbf{meters}``).
- `D_s::Float64`=NaN: ADD from observer to source (in ``\\rm \\mathbf{meters}``).
- `mass::Float64`=NaN: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).

# Returns
- `θE`: Einstein angle (in ``\\rm \\mathbf{arcseconds}``)
"""
function einstein_angle(;D_d::Float64=NaN, D_ds::Float64=NaN, D_s::Float64=NaN, mass::Float64=NaN)
  return  sqrt((4.0 * CONST_G * mass * MASS_SUN / CONST_C^2) * (D_ds / D_d / D_s)) / ANGLE_ARCSEC
end

end

