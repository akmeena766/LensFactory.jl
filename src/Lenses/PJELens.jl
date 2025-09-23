module PJELens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!

"""
    potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, θt::RV, ϵ::RV, pa::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, θt::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = 1.0 - ϵ

   # Get b_sie(q)
   bq = (4.0 * pi * (vd / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q) and a(q)
   sq = θs / sqrt(q)
   aq = θt / sqrt(q)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos(deg2rad(pa)) + (θy - θyc) * sin(deg2rad(pa))
   dy_r = - (θx - θxc) * sin(deg2rad(pa)) + (θy - θyc) * cos(deg2rad(pa))
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
    potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, θt::RV, ϵ::RV, pa::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, θt::RV, ϵ::RV, pa::RV) where T <: ROA
   # Get axis-ratio
   q = 1.0 - ϵ

   # Get b_sie(q)
   bq = (4.0 * pi * (vd / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q) and a(q)
   sq = θs / sqrt(q)
   aq = θt / sqrt(q)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Coordinate in the rotated frame
         dx_r = + (θx[i, j] - θxc) * cos(deg2rad(pa)) + (θy[i, j] - θyc) * sin(deg2rad(pa))
         dy_r = - (θx[i, j] - θxc) * sin(deg2rad(pa)) + (θy[i, j] - θyc) * cos(deg2rad(pa))
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
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, θt::RV, ϵ::RV, pa::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, θt::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = 1.0 - ϵ

   # Get b_sie(q)
   bq = (4.0 * pi * (vd / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q) and a(q)
   sq = θs / sqrt(q)
   aq = θt / sqrt(q)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos(deg2rad(pa)) + (θy - θyc) * sin(deg2rad(pa))
   dy_r = - (θx - θxc) * sin(deg2rad(pa)) + (θy - θyc) * cos(deg2rad(pa))
   ds_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)
   da_r = sqrt(q^2 * (aq^2 + dx_r^2) + dy_r^2)

   # Get deflection vector corresponding to θs in rotated frame
   ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (ds_r + sq))
   ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (ds_r + q^2 * sq))
   
   # Rotate back to original frame
   ψx_1 = ψx_r * cos(deg2rad(pa)) - ψy_r * sin(deg2rad(pa))
   ψy_1 = ψx_r * sin(deg2rad(pa)) + ψy_r * cos(deg2rad(pa))

   # Get the deflection vector corresponding to θt in rotated frame
   ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (da_r + aq))
   ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (da_r + q^2 * aq))

   # Rotate back to original frame
   ψx_2 = ψx_r * cos(deg2rad(pa)) - ψy_r * sin(deg2rad(pa))
   ψy_2 = ψx_r * sin(deg2rad(pa)) + ψy_r * cos(deg2rad(pa))

   ψx_up = ψx + (ψx_1 - ψx_2)
   ψy_up = ψy + (ψy_1 - ψy_2)
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, θt::RV, ϵ::RV, pa::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, θt::RV, ϵ::RV, pa::RV) where T <: ROA
   # Get axis-ratio
   q = 1.0 - ϵ

   # Get b_sie(q)
   bq = (4.0 * pi * (vd / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q) and a(q)
   sq = θs / sqrt(q)
   aq = θt / sqrt(q)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Coordinate in the rotated frame
         dx_r = + (θx[i, j] - θxc) * cos(deg2rad(pa)) + (θy[i, j] - θyc) * sin(deg2rad(pa))
         dy_r = - (θx[i, j] - θxc) * sin(deg2rad(pa)) + (θy[i, j] - θyc) * cos(deg2rad(pa))
         ds_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)
         da_r = sqrt(q^2 * (aq^2 + dx_r^2) + dy_r^2)

         # Get deflection vector corresponding to θs in rotated frame
         ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (ds_r + sq))
         ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (ds_r + q^2 * sq))
   
         # Rotate back to original frame
         ψx_1 = ψx_r * cos(deg2rad(pa)) - ψy_r * sin(deg2rad(pa))
         ψy_1 = ψx_r * sin(deg2rad(pa)) + ψy_r * cos(deg2rad(pa))

         # Get the deflection vector corresponding to θt in rotated frame
         ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (da_r + aq))
         ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (da_r + q^2 * aq))

         # Rotate back to original frame
         ψx_2 = ψx_r * cos(deg2rad(pa)) - ψy_r * sin(deg2rad(pa))
         ψy_2 = ψx_r * sin(deg2rad(pa)) + ψy_r * cos(deg2rad(pa))

         ψx[i, j] = ψx[i, j] + (ψx_1 - ψx_2)
         ψy[i, j] = ψy[i, j] + (ψy_1 - ψy_2)
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, θt::RV, ϵ::RV, pa::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, θt::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = 1.0 - ϵ

   # Get b_sie(q)
   bq = (4.0 * pi * (vd / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q) and a(q)
   sq = θs / sqrt(q)
   aq = θt / sqrt(q)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos(deg2rad(pa)) + (θy - θyc) * sin(deg2rad(pa))
   dy_r = - (θx - θxc) * sin(deg2rad(pa)) + (θy - θyc) * cos(deg2rad(pa))
   ds_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)
   da_r = sqrt(q^2 * (aq^2 + dx_r^2) + dy_r^2)

   # Get deformation tensor in rotated frame
   common_factor = (1+q^2) * sq^2 + 2 * ds_r * sq + dx_r^2 + dy_r^2
   ψxx_r = + bq * q * (q^2 * sq^2 + dy_r^2 + sq * ds_r) / ds_r / common_factor
   ψyy_r = + bq * q * (sq^2 + dx_r^2 + sq * ds_r) / ds_r / common_factor
   ψxy_r = - bq * q * dx_r * dy_r / ds_r / common_factor

   ψxx_1 = ψxx_r * cos(deg2rad(pa))^2 - ψxy_r * sin(deg2rad(2*pa)) + ψyy_r * sin(deg2rad(pa))^2
   ψyy_1 = ψxx_r * sin(deg2rad(pa))^2 + ψxy_r * sin(deg2rad(2*pa)) + ψyy_r * cos(deg2rad(pa))^2
   ψxy_1 = 0.5 * sin(deg2rad(2*pa)) * (ψxx_r - ψyy_r) + cos(deg2rad(2*pa)) * ψxy_r

   # Get deformation tensor in rotated frame
   common_factor = (1+q^2) * aq^2 + 2 * da_r * aq + dx_r^2 + dy_r^2
   ψxx_r = + bq * q * (q^2 * aq^2 + dy_r^2 + aq * da_r) / da_r / common_factor
   ψyy_r = + bq * q * (aq^2 + dx_r^2 + aq * da_r) / da_r / common_factor
   ψxy_r = - bq * q * dx_r * dy_r / da_r / common_factor

   ψxx_2 = ψxx_r * cos(deg2rad(pa))^2 - ψxy_r * sin(deg2rad(2*pa)) + ψyy_r * sin(deg2rad(pa))^2
   ψyy_2 = ψxx_r * sin(deg2rad(pa))^2 + ψxy_r * sin(deg2rad(2*pa)) + ψyy_r * cos(deg2rad(pa))^2
   ψxy_2 = 0.5 * sin(deg2rad(2*pa)) * (ψxx_r - ψyy_r) + cos(deg2rad(2*pa)) * ψxy_r
   
   ψxx_up = ψxx + (ψxx_1 - ψxx_2)
   ψyy_up = ψyy + (ψyy_1 - ψyy_2)
   ψxy_up = ψxy + (ψxy_1 - ψxy_2)
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, θt::RV, ϵ::RV, pa::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, θt::RV, ϵ::RV, pa::RV) where T <: ROA
   # Get axis-ratio
   q = 1.0 - ϵ

   # Get b_sie(q)
   bq = (4.0 * pi * (vd / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q) and a(q)
   sq = θs / sqrt(q)
   aq = θt / sqrt(q)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Coordinate in the rotated frame
         dx_r = + (θx[i, j] - θxc) * cos(deg2rad(pa)) + (θy[i, j] - θyc) * sin(deg2rad(pa))
         dy_r = - (θx[i, j] - θxc) * sin(deg2rad(pa)) + (θy[i, j] - θyc) * cos(deg2rad(pa))
         ds_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)
         da_r = sqrt(q^2 * (aq^2 + dx_r^2) + dy_r^2)

         # Get deformation tensor in rotated frame
         common_factor = (1+q^2) * sq^2 + 2 * ds_r * sq + dx_r^2 + dy_r^2
         ψxx_r = + bq * q * (q^2 * sq^2 + dy_r^2 + sq * ds_r) / ds_r / common_factor
         ψyy_r = + bq * q * (sq^2 + dx_r^2 + sq * ds_r) / ds_r / common_factor
         ψxy_r = - bq * q * dx_r * dy_r / ds_r / common_factor

         ψxx_1 = ψxx_r * cos(deg2rad(pa))^2 - ψxy_r * sin(deg2rad(2*pa)) + ψyy_r * sin(deg2rad(pa))^2
         ψyy_1 = ψxx_r * sin(deg2rad(pa))^2 + ψxy_r * sin(deg2rad(2*pa)) + ψyy_r * cos(deg2rad(pa))^2
         ψxy_1 = 0.5 * sin(deg2rad(2*pa)) * (ψxx_r - ψyy_r) + cos(deg2rad(2*pa)) * ψxy_r

         # Get deformation tensor in rotated frame
         common_factor = (1+q^2) * aq^2 + 2 * da_r * aq + dx_r^2 + dy_r^2
         ψxx_r = + bq * q * (q^2 * aq^2 + dy_r^2 + aq * da_r) / da_r / common_factor
         ψyy_r = + bq * q * (aq^2 + dx_r^2 + aq * da_r) / da_r / common_factor
         ψxy_r = - bq * q * dx_r * dy_r / da_r / common_factor

         ψxx_2 = ψxx_r * cos(deg2rad(pa))^2 - ψxy_r * sin(deg2rad(2*pa)) + ψyy_r * sin(deg2rad(pa))^2
         ψyy_2 = ψxx_r * sin(deg2rad(pa))^2 + ψxy_r * sin(deg2rad(2*pa)) + ψyy_r * cos(deg2rad(pa))^2
         ψxy_2 = 0.5 * sin(deg2rad(2*pa)) * (ψxx_r - ψyy_r) + cos(deg2rad(2*pa)) * ψxy_r

         ψxx[i, j] = ψxx[i, j] + (ψxx_1 - ψxx_2)
         ψyy[i, j] = ψyy[i, j] + (ψyy_1 - ψyy_2)
         ψxy[i, j] = ψxy[i, j] + (ψxy_1 - ψxy_2)
      end
   end
end

end