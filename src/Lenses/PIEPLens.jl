module PIEPLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!

"""
    potential!(ψ::Real, θx::Real, θy::Real, θxc::T, θyc::T, v_d::T, θs::T, ϵ::T, pa::T) where T <: Real
"""
function potential!(ψ::Real, θx::Real, θy::Real, θxc::T, θyc::T, v_d::T, θs::T, ϵ::T, pa::T) where T <: Real
   θE = 4.0 * pi * (v_d * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Pre-compute angles
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos_pa + (θy - θyc) * sin_pa
   dy_r = - (θx - θxc) * sin_pa + (θy - θyc) * cos_pa

   ψ_up = ψ + θE * sqrt(θs^2 + dx_r^2 + dy_r^2 / q^2)
   return ψ_up
end

"""
    potential!(ψ::ROA, θx::ROA, θy::ROA, θxc::T, θyc::T, v_d::T, θs::T, ϵ::T, pa::T) where T <: Real
Calculate potential at given coordinates for PIEP lens and update the potential values in-place.

# Arguments
- `ψ`: Potential at given coordinates
- `θx` : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy` : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d`: Velocity dispersion of the lens (in ``\\rm \\mathbf{km/s}``).
- `θs` : Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
- `ϵ`  : Ellipticity of the lens.
- `pa` : Position angle of the lens (in ``\\rm \\mathbf{degrees}``).
"""
function potential!(ψ::ROA, θx::ROA, θy::ROA, θxc::T, θyc::T, v_d::T, θs::T, ϵ::T, pa::T) where T <: Real
   θE = 4.0 * pi * (v_d * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Pre-compute angles
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx_r = + (θx[i, j] - θxc) * cos_pa + (θy[i, j] - θyc) * sin_pa
         dy_r = - (θx[i, j] - θxc) * sin_pa + (θy[i, j] - θyc) * cos_pa

         ψ[i, j] = ψ[i, j] + θE * sqrt(θs^2 + dx_r^2 + dy_r^2 / q^2)
      end
   end
end


"""
    deflection!(ψx::Real, ψy::Real, θx::Real, θy::Real, θxc::T, θyc::T, v_d::T, θs::T, ϵ::T, pa::T) where T <: Real
"""
function deflection!(ψx::Real, ψy::Real, θx::Real, θy::Real, θxc::T, θyc::T, v_d::T, θs::T, ϵ::T, pa::T) where T <: Real
   θE = 4.0 * pi * (v_d * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Pre-compute angles
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos_pa + (θy - θyc) * sin_pa
   dy_r = - (θx - θxc) * sin_pa + (θy - θyc) * cos_pa
   dr_r = sqrt(θs^2 + dx_r^2 + dy_r^2 / q^2)

   # Deflection in the rotated frame
   ψx_r = θE * dx_r / dr_r
   ψy_r = θE * dy_r / dr_r / q^2

   # Rotate back to original frame
   ψx_up = ψx + ψx_r * cos_pa - ψy_r * sin_pa
   ψy_up = ψy + ψx_r * sin_pa + ψy_r * cos_pa

   return ψx_up, ψy_up
end

"""
    deflection!(ψx::ROA, ψy::ROA, θx::ROA, θy::ROA, θxc::T, θyc::T, v_d::T, θs::T, ϵ::T, pa::T) where T <: Real
Calculate deflection at given coordinates for PIEP lens and update the deflection values in-place.

# Arguments
- `ψx` : x-component of deflection at given coordinates
- `ψy` : y-component of deflection at given coordinates
- `θx` : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy` : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d`: Velocity dispersion of the lens (in ``\\rm \\mathbf{km/s}``).
- `θs` : Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
- `ϵ`  : Ellipticity of the lens.
- `pa` : Position angle of the lens (in ``\\rm \\mathbf{degrees}``).
"""
function deflection!(ψx::ROA, ψy::ROA, θx::ROA, θy::ROA, θxc::T, θyc::T, v_d::T, θs::T, ϵ::T, pa::T) where T <: Real
   θE = 4.0 * pi * (v_d * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Pre-compute angles
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Coordinate in the rotated frame
         dx_r = + (θx[i, j] - θxc) * cos_pa + (θy[i, j] - θyc) * sin_pa
         dy_r = - (θx[i, j] - θxc) * sin_pa + (θy[i, j] - θyc) * cos_pa
         dr_r = sqrt(θs^2 + dx_r^2 + dy_r^2 / q^2)

         # Deflection in the rotated frame
         ψx_r = θE * dx_r / dr_r
         ψy_r = θE * dy_r / dr_r / q^2

         # Rotate back to original frame
         ψx[i, j] = ψx[i, j] + ψx_r * cos_pa - ψy_r * sin_pa
         ψy[i, j] = ψy[i, j] + ψx_r * sin_pa + ψy_r * cos_pa
      end
   end
end


"""
    jacobian!(ψxx::Real, ψyy::Real, ψxy::Real, θx::Real, θy::Real, θxc::T, θyc::T, v_d::T, θs::T, ϵ::T, pa::T) where T <: Real
"""
function jacobian!(ψxx::Real, ψyy::Real, ψxy::Real, θx::Real, θy::Real, θxc::T, θyc::T, v_d::T, θs::T, ϵ::T, pa::T) where T <: Real
   θE = 4.0 * pi * (v_d * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Pre-compute angles
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)
   sin_2pa = sin(2.0 * pa_rad)
   cos_2pa = cos(2.0 * pa_rad)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos_pa + (θy - θyc) * sin_pa
   dy_r = - (θx - θxc) * sin_pa + (θy - θyc) * cos_pa
   dr_r = sqrt(θs^2 + dx_r^2 + dy_r^2 / q^2)

   # Deformation tensor components in rotated frame
   ψxx_r = + θE * (θs^2 + dy_r^2 / q^2) / dr_r^3
   ψyy_r = + θE * (θs^2 + dx_r^2) / dr_r^3 / q^2
   ψxy_r = - θE * dx_r * dy_r / dr_r^3 / q^2

   # Rotate back to the original frame
   ψxx_up = ψxx + ψxx_r * cos_pa^2 - ψxy_r * sin_2pa + ψyy_r * sin_pa^2
   ψyy_up = ψyy + ψxx_r * sin_pa^2 + ψxy_r * sin_2pa + ψyy_r * cos_pa^2
   ψxy_up = ψxy + 0.5 * sin_2pa * (ψxx_r - ψyy_r) + cos_2pa * ψxy_r

   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::ROA, ψyy::ROA, ψxy::ROA, θx::ROA, θy::ROA, θxc::T, θyc::T, v_d::T, θs::T, ϵ::T, pa::T) where T <: Real
Calculate Jacobian at given coordinates for PIEP lens and update the Jacobian values in-place.

# Arguments
- `ψxx`: x-component of Jacobian at given coordinates
- `ψyy`: y-component of Jacobian at given coordinates
- `ψxy`: xy-component of Jacobian at given coordinates
- `θx` : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy` : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d`: Velocity dispersion of the lens (in ``\\rm \\mathbf{km/s}``).
- `θs` : Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
- `ϵ`  : Ellipticity of the lens.
- `pa` : Position angle of the lens (in ``\\rm \\mathbf{degrees}``).
"""
function jacobian!(ψxx::ROA, ψyy::ROA, ψxy::ROA, θx::ROA, θy::ROA, θxc::T, θyc::T, v_d::T, θs::T, ϵ::T, pa::T) where T <: Real
   θE = 4.0 * pi * (v_d * 1.0E3 / CONST_C)^2 / ANGLE_ARCSEC
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Pre-compute angles
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)
   sin_2pa = sin(2.0 * pa_rad)
   cos_2pa = cos(2.0 * pa_rad)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Coordinate in the rotated frame
         dx_r = + (θx[i, j] - θxc) * cos_pa + (θy[i, j] - θyc) * sin_pa
         dy_r = - (θx[i, j] - θxc) * sin_pa + (θy[i, j] - θyc) * cos_pa
         dr_r = sqrt(θs^2 + dx_r^2 + dy_r^2 / q^2)

         # Deformation tensor components in rotated frame
         ψxx_r = + θE * (θs^2 + dy_r^2 / q^2) / dr_r^3
         ψyy_r = + θE * (θs^2 + dx_r^2) / dr_r^3 / q^2
         ψxy_r = - θE * dx_r * dy_r / dr_r^3 / q^2

         # Rotate back to the original frame
         ψxx[i, j] = ψxx[i, j] + ψxx_r * cos_pa^2 - ψxy_r * sin_2pa + ψyy_r * sin_pa^2
         ψyy[i, j] = ψyy[i, j] + ψxx_r * sin_pa^2 + ψxy_r * sin_2pa + ψyy_r * cos_pa^2
         ψxy[i, j] = ψxy[i, j] + 0.5 * sin_2pa * (ψxx_r - ψyy_r) + cos_2pa * ψxy_r
      end
   end
end

end