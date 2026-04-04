module SIELens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!


"""
    potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, v_d::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, v_d::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Get b_sie(q)
   bq = (4.0 * pi * (v_d * 1.0E3 / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q)
   sq = θs / sqrt(q)

   # Pre-compute angles
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos_pa + (θy - θyc) * sin_pa
   dy_r = - (θx - θxc) * sin_pa + (θy - θyc) * cos_pa
   dr_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)

   # Get deflection vector in rotated frame
   ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (dr_r + sq))
   ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (dr_r + q^2 * sq))

   # Get potential
   ψ_up = ψ + dx_r * ψx_r + dy_r * ψy_r + bq * q * sq * log((1.0 + q) * sq / sqrt((dr_r + sq)^2 + (1.0 - q^2) * dx_r^2) + 1E-12)
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, v_d::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
Calculate potential at given coordinates for SIE lens and update the potential values in-place.

# Arguments
- `ψ`: Potential at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc::RV`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc::RV`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d::RV`: Velocity dispersion of the lens (in ``\\rm \\mathbf{km/s}``).
- `θs::RV`: Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
- `ϵ::RV`: Ellipticity of the lens.
- `pa::RV`: Position angle of the lens (in ``\\rm \\mathbf{degrees}``).
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, v_d::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Get b_sie(q)
   bq = (4.0 * pi * (v_d * 1.0E3 / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q)
   sq = θs / sqrt(q)

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
         dr_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)

         # Get deflection vector in rotated frame
         ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (dr_r + sq))
         ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (dr_r + q^2 * sq))

         ψ[i, j] = ψ[i, j] + dx_r * ψx_r + dy_r * ψy_r + bq * q * sq * log((1.0 + q) * sq / sqrt((dr_r + sq)^2 + (1.0 - q^2) * dx_r^2) + 1E-12)
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, v_d::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, v_d::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Get b_sie(q)
   bq = (4.0 * pi * (v_d * 1.0E3 / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q)
   sq = θs / sqrt(q)

   # Pre-compute angles
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos_pa + (θy - θyc) * sin_pa
   dy_r = - (θx - θxc) * sin_pa + (θy - θyc) * cos_pa
   dr_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)

   # Get deflection vector in rotated frame
   ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (dr_r + sq))
   ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (dr_r + q^2 * sq))

   # Rotate back to original frame
   ψx_up = ψx + ψx_r * cos_pa - ψy_r * sin_pa
   ψy_up = ψy + ψx_r * sin_pa + ψy_r * cos_pa

   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, v_d::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
Calculate deflection at given coordinates for SIE lens and update the deflection values in-place.

# Arguments
- `ψx`: x-component of deflection at given coordinates
- `ψy`: y-component of deflection at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc::RV`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc::RV`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d::RV`: Velocity dispersion of the lens (in ``\\rm \\mathbf{km/s}``).
- `θs::RV`: Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
- `ϵ::RV`: Ellipticity of the lens.
- `pa::RV`: Position angle of the lens (in ``\\rm \\mathbf{degrees}``).
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, v_d::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Get b_sie(q)
   bq = (4.0 * pi * (v_d * 1.0E3 / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q)
   sq = θs / sqrt(q)

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
         dr_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)

         # Get deflection vector in rotated frame
         ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (dr_r + sq))
         ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (dr_r + q^2 * sq))

         # Rotate back to original frame
         ψx[i, j] = ψx[i, j] + ψx_r * cos_pa - ψy_r * sin_pa
         ψy[i, j] = ψy[i, j] + ψx_r * sin_pa + ψy_r * cos_pa
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, v_d::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, v_d::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Get b_sie(q)
   bq = (4.0 * pi * (v_d * 1.0E3 / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q)
   sq = θs / sqrt(q)

   # Pre-compute angles
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)
   cos_2pa = cos(2 * pa_rad)
   sin_2pa = sin(2 * pa_rad)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos_pa + (θy - θyc) * sin_pa
   dy_r = - (θx - θxc) * sin_pa + (θy - θyc) * cos_pa
   dr_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)

   # Common factor
   common_factor = (1+q^2) * sq^2 + 2 * dr_r * sq + dx_r^2 + dy_r^2

   # Get deformation tensor in rotated frame
   ψxx_r = + bq * q * (q^2 * sq^2 + dy_r^2 + sq * dr_r) / dr_r / common_factor
   ψyy_r = + bq * q * (sq^2 + dx_r^2 + sq * dr_r) / dr_r / common_factor
   ψxy_r = - bq * q * dx_r * dy_r / dr_r / common_factor

   # Rotate back to the original frame
   ψxx_up = ψxx + ψxx_r * cos_pa^2 - ψxy_r * sin_2pa + ψyy_r * sin_pa^2
   ψyy_up = ψyy + ψxx_r * sin_pa^2 + ψxy_r * sin_2pa + ψyy_r * cos_pa^2
   ψxy_up = ψxy + 0.5 * sin_2pa * (ψxx_r - ψyy_r) + cos_2pa * ψxy_r

   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, v_d::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
Calculate Jacobian at given coordinates for SIE lens and update the Jacobian values in-place.

# Arguments
- `ψxx`: xx-component of Jacobian at given coordinates
- `ψyy`: yy-component of Jacobian at given coordinates
- `ψxy`: xy-component of Jacobian at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc::RV`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc::RV`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d::RV`: Velocity dispersion of the lens (in ``\\rm \\mathbf{km/s}``).
- `θs::RV`: Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
- `ϵ::RV`: Ellipticity of the lens.
- `pa::RV`: Position angle of the lens (in ``\\rm \\mathbf{degrees}``).
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, v_d::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Get b_sie(q)
   bq = (4.0 * pi * (v_d * 1.0E3 / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q)
   sq = θs / sqrt(q)

   # Pre-compute angles
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)
   cos_2pa = cos(2 * pa_rad)
   sin_2pa = sin(2 * pa_rad)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Coordinate in the rotated frame
         dx_r = + (θx[i, j] - θxc) * cos_pa + (θy[i, j] - θyc) * sin_pa
         dy_r = - (θx[i, j] - θxc) * sin_pa + (θy[i, j] - θyc) * cos_pa
         dr_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)

         # Common factor
         common_factor = (1+q^2) * sq^2 + 2 * dr_r * sq + dx_r^2 + dy_r^2

         # Get deformation tensor in rotated frame
         ψxx_r = + bq * q * (q^2 * sq^2 + dy_r^2 + sq * dr_r) / dr_r / common_factor
         ψyy_r = + bq * q * (sq^2 + dx_r^2 + sq * dr_r) / dr_r / common_factor
         ψxy_r = - bq * q * dx_r * dy_r / dr_r / common_factor

         # Rotate back to the original frame
         ψxx[i, j] = ψxx[i, j] + ψxx_r * cos_pa^2 - ψxy_r * sin_2pa + ψyy_r * sin_pa^2
         ψyy[i, j] = ψyy[i, j] + ψxx_r * sin_pa^2 + ψxy_r * sin_2pa + ψyy_r * cos_pa^2
         ψxy[i, j] = ψxy[i, j] + 0.5 * sin_2pa * (ψxx_r - ψyy_r) + cos_2pa * ψxy_r
      end
   end
end

end