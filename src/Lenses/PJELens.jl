module PJELens


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..Constants


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export potential!
export deflection!
export jacobian!


# --------------------------------------------------------------------------------------------------
# Main functions
# --------------------------------------------------------------------------------------------------
"""
    potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T, θt::T, ϵ::T, pa::T) where {U<:Real, S<:Real, T<:Real}
"""
function potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T, θt::T, ϵ::T, pa::T) where {U<:Real, S<:Real, T<:Real}
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Get b_sie(q)
   bq = (4.0 * pi * (v_d * 1.0E3 / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q) and a(q)
   sq = θs / sqrt(q)
   aq = θt / sqrt(q)

   # Pre-compute angles
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos_pa + (θy - θyc) * sin_pa
   dy_r = - (θx - θxc) * sin_pa + (θy - θyc) * cos_pa
   ds_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)
   da_r = sqrt(q^2 * (aq^2 + dx_r^2) + dy_r^2)

   # Get deflection vector corresponding to θs in rotated frame
   ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (ds_r + sq))
   ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (ds_r + q^2 * sq))
   ψ1 = dx_r * ψx_r + dy_r * ψy_r + bq * q * sq * log((1.0 + q) * sq / sqrt((ds_r + sq)^2 + (1.0 - q^2) * dx_r^2))

   # Get the deflection vector corresponding to θt in rotated frame
   ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (da_r + aq))
   ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (da_r + q^2 * aq))
   ψ2 = dx_r * ψx_r + dy_r * ψy_r + bq * q * aq * log((1.0 + q) * aq / sqrt((da_r + aq)^2 + (1.0 - q^2) * dx_r^2))

   ψ_up = ψ + (ψ1 - ψ2)
   return ψ_up
end

"""
    potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T, θt::T, ϵ::T, pa::T) where {U<:ROA, S<:ROA, T<:Real}
Calculate potential at given coordinates for PJE lens and update the potential values in-place.

# Arguments
- `ψ`: Potential at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc::Real`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc::Real`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d::Real`: Velocity dispersion of the lens (in ``\\rm \\mathbf{km/s}``).
- `θs::Real`: Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
- `θt::Real`: Truncation radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
- `ϵ::Real`: Ellipticity of the lens.
- `pa::Real`: Position angle of the lens (in ``\\rm \\mathbf{degrees}``).
"""
function potential!(ψ::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T, θt::T, ϵ::T, pa::T) where {U<:ROA, S<:ROA, T<:Real}
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Get b_sie(q)
   bq = (4.0 * pi * (v_d * 1.0E3 / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q) and a(q)
   sq = θs / sqrt(q)
   aq = θt / sqrt(q)

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
         ds_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)
         da_r = sqrt(q^2 * (aq^2 + dx_r^2) + dy_r^2)

         # Get deflection vector corresponding to θs in rotated frame
         ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (ds_r + sq))
         ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (ds_r + q^2 * sq))
         ψ1 = dx_r * ψx_r + dy_r * ψy_r + bq * q * sq * log((1.0 + q) * sq / sqrt((ds_r + sq)^2 + (1.0 - q^2) * dx_r^2))

         # Get the deflection vector corresponding to θt in rotated frame
         ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (da_r + aq))
         ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (da_r + q^2 * aq))
         ψ2 = dx_r * ψx_r + dy_r * ψy_r + bq * q * aq * log((1.0 + q) * aq / sqrt((da_r + aq)^2 + (1.0 - q^2) * dx_r^2))

         ψ[i, j] = ψ[i, j] + (ψ1 - ψ2)
      end
   end
end


"""
    deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T, θt::T, ϵ::T, pa::T) where {U<:Real, S<:Real, T<:Real}
"""
function deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T, θt::T, ϵ::T, pa::T) where {U<:Real, S<:Real, T<:Real}
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Get b_sie(q)
   bq = (4.0 * pi * (v_d * 1.0E3 / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q) and a(q)
   sq = θs / sqrt(q)
   aq = θt / sqrt(q)

   # Pre-compute angles
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos_pa + (θy - θyc) * sin_pa
   dy_r = - (θx - θxc) * sin_pa + (θy - θyc) * cos_pa
   ds_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)
   da_r = sqrt(q^2 * (aq^2 + dx_r^2) + dy_r^2)

   # Get deflection vector corresponding to θs in rotated frame
   ψx_r1 =  atan(sqrt(1 - q^2) * dx_r / (ds_r + sq))
   ψy_r1 = atanh(sqrt(1 - q^2) * dy_r / (ds_r + q^2 * sq))
   
   # Get the deflection vector corresponding to θt in rotated frame
   ψx_r2 =  atan(sqrt(1 - q^2) * dx_r / (da_r + aq))
   ψy_r2 = atanh(sqrt(1 - q^2) * dy_r / (da_r + q^2 * aq))

   # Add the two components
   ψx_r = (bq * q / sqrt(1 - q^2)) * (ψx_r1 - ψx_r2)
   ψy_r = (bq * q / sqrt(1 - q^2)) * (ψy_r1 - ψy_r2)

   # Rotate back to original frame and update the values
   ψx_up = ψx + ψx_r * cos_pa - ψy_r * sin_pa
   ψy_up = ψy + ψx_r * sin_pa + ψy_r * cos_pa
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T, θt::T, ϵ::T, pa::T) where {U<:ROA, S<:ROA, T<:Real}
Calculate deflection at given coordinates for PJE lens and update the deflection values in-place.

# Arguments
- `ψx` : x-component of deflection at given coordinates
- `ψy` : y-component of deflection at given coordinates
- `θx` : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy` : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d`: Velocity dispersion of the lens (in ``\\rm \\mathbf{km/s}``).
- `θs` : Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
- `θt` : Truncation radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
- `ϵ`  : Ellipticity of the lens.
- `pa` : Position angle of the lens (in ``\\rm \\mathbf{degrees}``).
"""
function deflection!(ψx::U, ψy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T, θt::T, ϵ::T, pa::T) where {U<:ROA, S<:ROA, T<:Real}
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Get b_sie(q)
   bq = (4.0 * pi * (v_d * 1.0E3 / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q) and a(q)
   sq = θs / sqrt(q)
   aq = θt / sqrt(q)

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
         ds_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)
         da_r = sqrt(q^2 * (aq^2 + dx_r^2) + dy_r^2)

         # Get deflection vector corresponding to θs in rotated frame
         ψx_r1 =  atan(sqrt(1 - q^2) * dx_r / (ds_r + sq))
         ψy_r1 = atanh(sqrt(1 - q^2) * dy_r / (ds_r + q^2 * sq))
   
         # Get the deflection vector corresponding to θt in rotated frame
         ψx_r2 =  atan(sqrt(1 - q^2) * dx_r / (da_r + aq))
         ψy_r2 = atanh(sqrt(1 - q^2) * dy_r / (da_r + q^2 * aq))

         # Add the two components
         ψx_r = (bq * q / sqrt(1 - q^2)) * (ψx_r1 - ψx_r2)
         ψy_r = (bq * q / sqrt(1 - q^2)) * (ψy_r1 - ψy_r2)

         # Rotate back to original frame and update the values
         ψx[i, j] = ψx[i, j] + ψx_r * cos_pa - ψy_r * sin_pa
         ψy[i, j] = ψy[i, j] + ψx_r * sin_pa + ψy_r * cos_pa
      end
   end
end


"""
    jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T, θt::T, ϵ::T, pa::T) where {U<:Real, S<:Real, T<:Real}
"""
function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T, θt::T, ϵ::T, pa::T) where {U<:Real, S<:Real, T<:Real}
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Get b_sie(q)
   bq = (4.0 * pi * (v_d * 1.0E3 / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q) and a(q)
   sq = θs / sqrt(q)
   aq = θt / sqrt(q)

   # Pre-compute angles
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)
   sin_2pa = sin(2.0 * pa_rad)
   cos_2pa = cos(2.0 * pa_rad)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos_pa + (θy - θyc) * sin_pa
   dy_r = - (θx - θxc) * sin_pa + (θy - θyc) * cos_pa
   ds_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)
   da_r = sqrt(q^2 * (aq^2 + dx_r^2) + dy_r^2)

   # Get deformation tensor in rotated frame
   common_factor = (1+q^2) * sq^2 + 2 * ds_r * sq + dx_r^2 + dy_r^2
   ψxx_r1 = + (q^2 * sq^2 + dy_r^2 + sq * ds_r) / ds_r / common_factor
   ψyy_r1 = + (sq^2 + dx_r^2 + sq * ds_r) / ds_r / common_factor
   ψxy_r1 = - dx_r * dy_r / ds_r / common_factor

   # Get deformation tensor in rotated frame
   common_factor = (1+q^2) * aq^2 + 2 * da_r * aq + dx_r^2 + dy_r^2
   ψxx_r2 = + (q^2 * aq^2 + dy_r^2 + aq * da_r) / da_r / common_factor
   ψyy_r2 = + (aq^2 + dx_r^2 + aq * da_r) / da_r / common_factor
   ψxy_r2 = - dx_r * dy_r / da_r / common_factor

   # Add the two components
   ψxx_r = bq * q * (ψxx_r1 - ψxx_r2)
   ψyy_r = bq * q * (ψyy_r1 - ψyy_r2)
   ψxy_r = bq * q * (ψxy_r1 - ψxy_r2)

   # Rotate back to original frame and update the values
   ψxx_up = ψxx + ψxx_r * cos_pa^2 - ψxy_r * sin_2pa + ψyy_r * sin_pa^2
   ψyy_up = ψyy + ψxx_r * sin_pa^2 + ψxy_r * sin_2pa + ψyy_r * cos_pa^2
   ψxy_up = ψxy + 0.5 * sin_2pa * (ψxx_r - ψyy_r) + cos_2pa * ψxy_r
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T, θt::T, ϵ::T, pa::T) where {U<:ROA, S<:ROA, T<:Real}
Calculate Jacobian at given coordinates for PJE lens and update the Jacobian values in-place.

# Arguments
- `ψxx`: xx-component of Jacobian at given coordinates
- `ψyy`: yy-component of Jacobian at given coordinates
- `ψxy`: xy-component of Jacobian at given coordinates
- `θx` : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy` : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θxc`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d`: Velocity dispersion of the lens (in ``\\rm \\mathbf{km/s}``).
- `θs` : Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
- `θt` : Truncation radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
- `ϵ`  : Ellipticity of the lens.
- `pa` : Position angle of the lens (in ``\\rm \\mathbf{degrees}``).
"""
function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, θxc::T, θyc::T, v_d::T, θs::T, θt::T, ϵ::T, pa::T) where {U<:ROA, S<:ROA, T<:Real}
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)

   # Get b_sie(q)
   bq = (4.0 * pi * (v_d * 1.0E3 / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q) and a(q)
   sq = θs / sqrt(q)
   aq = θt / sqrt(q)

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
         ds_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)
         da_r = sqrt(q^2 * (aq^2 + dx_r^2) + dy_r^2)

         # Get deformation tensor in rotated frame
         common_factor = (1+q^2) * sq^2 + 2 * ds_r * sq + dx_r^2 + dy_r^2
         ψxx_r1 = + (q^2 * sq^2 + dy_r^2 + sq * ds_r) / ds_r / common_factor
         ψyy_r1 = + (sq^2 + dx_r^2 + sq * ds_r) / ds_r / common_factor
         ψxy_r1 = - dx_r * dy_r / ds_r / common_factor

         # Get deformation tensor in rotated frame
         common_factor = (1+q^2) * aq^2 + 2 * da_r * aq + dx_r^2 + dy_r^2
         ψxx_r2 = + (q^2 * aq^2 + dy_r^2 + aq * da_r) / da_r / common_factor
         ψyy_r2 = + (aq^2 + dx_r^2 + aq * da_r) / da_r / common_factor
         ψxy_r2 = - dx_r * dy_r / da_r / common_factor

         # Add the two components
         ψxx_r = bq * q * (ψxx_r1 - ψxx_r2)
         ψyy_r = bq * q * (ψyy_r1 - ψyy_r2)
         ψxy_r = bq * q * (ψxy_r1 - ψxy_r2)

         # Rotate back to original frame and update the values
         ψxx[i, j] = ψxx[i, j] + ψxx_r * cos_pa^2 - ψxy_r * sin_2pa + ψyy_r * sin_pa^2
         ψyy[i, j] = ψyy[i, j] + ψxx_r * sin_pa^2 + ψxy_r * sin_2pa + ψyy_r * cos_pa^2
         ψxy[i, j] = ψxy[i, j] + 0.5 * sin_2pa * (ψxx_r - ψyy_r) + cos_2pa * ψxy_r
      end
   end
end

end