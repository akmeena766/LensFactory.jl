module NSISPLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!
export einstein_angle


"""
    potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T) where {U<:Real, S<:Real, T<:Real}
"""
function potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T) where {U<:Real, S<:Real, T<:Real}
   θE = 4π * (v_d * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   θs2 = θs^2

   dx = θx - θxc
   dy = θy - θyc

   ψ_up = ψ + θE * sqrt(θs2 + dx^2 + dy^2)
   return ψ_up
end

"""
    potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T) where {U<:ROA, S<:ROA, T<:Real}
Calculate potential at given coordinates for NSISP lens and update the potential (ψ) in place.
The lensing potential is given as,

```math
\\begin{align*} 
ψ(θ_x, θ_y) = 4 π \\left(\\frac{v_d}{{\\rm c}} \\right)^2 \\sqrt{θ_s^2 + |\\pmb{θ} - \\pmb{θ}_c|^2}. 
\\end{align*}
```

# Arguments
- `ψ`  : Potential at given coordinates
- `θx` : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy` : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d`: Velocity dispersion (in ``\\rm \\mathbf{km/s}``).
- `θs` : Core radius of the lens (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `nothing`: Updates the potential (ψ) in place.
"""
function potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T) where {U<:ROA, S<:ROA, T<:Real}
   θE = 4π * (v_d * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   θs2 = θs^2

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         ψ[i, j] = ψ[i, j] + θE * sqrt(θs2 + dx^2 + dy^2)
      end
   end
   return nothing
end


"""
    deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T) where {U<:Real, S<:Real, T<:Real}
"""
function deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T) where {U<:Real, S<:Real, T<:Real}
   θE = 4π * (v_d * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   θs2 = θs^2

   dx = θx - θxc
   dy = θy - θyc
   dr = sqrt(θs2 + dx^2 + dy^2)

   ψx_up = ψx + θE * dx / dr
   ψy_up = ψy + θE * dy / dr
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T) where {U<:ROA, S<:ROA, T<:Real}
Calculate deflection at given coordinates for NSISP lens and update the deflection components
(ψx, ψy) in place.

# Arguments
- `ψx` : x-component of the deflection at given coordinates
- `ψy` : y-component of the deflection at given coordinates
- `θx` : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy` : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d`: Velocity dispersion (in ``\\rm \\mathbf{km/s}``).
- `θs` : Core radius of the lens (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `nothing`: Updates the deflection (ψx, ψy) in place.
"""
function deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T) where {U<:ROA, S<:ROA, T<:Real}
   θE = 4π * (v_d * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   θs2 = θs^2

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         dr = sqrt(θs2 + dx^2 + dy^2)

         ψx[i, j] = ψx[i, j] + θE * dx / dr
         ψy[i, j] = ψy[i, j] + θE * dy / dr
      end
   end
   return nothing
end


"""
    jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T) where {U<:Real, S<:Real, T<:Real}
"""
function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T) where {U<:Real, S<:Real, T<:Real}
   θE = 4π * (v_d * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   θs2 = θs^2
   
   dx = θx - θxc
   dy = θy - θyc
   dr2 = θs2 + dx^2 + dy^2
   dr3 = dr2 * sqrt(dr2)

   ψxx_up = ψxx + θE * (θs2 + dy^2) / dr3
   ψyy_up = ψyy + θE * (θs2 + dx^2) / dr3
   ψxy_up = ψxy - θE * dx * dy / dr3
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T) where {U<:ROA, S<:ROA, T<:Real}
Calculate jacobian at given coordinates for NSISP lens and and update the jacobian components 
(ψxx, ψyy, ψxy) in place.

# Arguments
- `ψxx`: xx-component of the jacobian at given coordinates
- `ψyy`: yy-component of the jacobian at given coordinates
- `ψxy`: xy-component of the jacobian at given coordinates
- `θx` : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy` : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d`: Velocity dispersion (in ``\\rm \\mathbf{km/s}``).
- `θs` : Core radius of the lens (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `nothing`: Updates the jacobian (ψxx, ψyy, ψxy) in place.
"""
function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T) where {U<:ROA, S<:ROA, T<:Real}
   θE = 4π * (v_d * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   θs2 = θs^2
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         dr2 = θs2 + dx^2 + dy^2
         dr3 = dr2 * sqrt(dr2)

         ψxx[i, j] = ψxx[i, j] + θE * (θs2 + dy^2) / dr3
         ψyy[i, j] = ψyy[i, j] + θE * (θs2 + dx^2) / dr3
         ψxy[i, j] = ψxy[i, j] - θE * dx * dy / dr3
      end
   end
   return nothing
end


"""
    einstein_angle(; D_ds::Real = NaN, 
                     D_s::Real  = NaN, 
                     v_d::Real  = NaN, 
                     x_s::Real  = NaN)
Calculate the Einstein angle for NSIS lens,
```math
\\theta_E = \\sqrt{\\left[4 \\pi \\frac{D_{ds}}{D_s} \\left( \\frac{v_d}{{\\rm c}} \\right)^2 \\right]^2 - x_s^2}.
```

# Keyword Arguments
- `D_ds`: ADD from the observer to the lens (in ``\\rm \\mathbf{meters}``).
- `D_s`: ADD from the observer to the source (in ``\\rm \\mathbf{meters}``).
- `v_d`: Velocity dispersion (in ``\\rm \\mathbf{km/s}``).
- `x_s`: Core radius of the lens (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `θE`: Einstein angle (in ``\\rm \\mathbf{arcseconds}``).
"""
function einstein_angle(; D_ds::Real=NaN, D_s::Real=NaN, v_d::Real=NaN, x_s::Real=NaN)
   return sqrt((4π * (v_d * 1.0E3 / CONST_C)^2 * (D_ds/D_s) / ANGLE_ARCSEC)^2 - x_s^2)
end

end