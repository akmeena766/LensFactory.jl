module SIELens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!


"""
    potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = 1.0 - ϵ

   # Get b_sie(q)
   bq = (4.0 * pi * (vd / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q)
   sq = θs / sqrt(q)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos(deg2rad(pa)) + (θy - θyc) * sin(deg2rad(pa))
   dy_r = - (θx - θxc) * sin(deg2rad(pa)) + (θy - θyc) * cos(deg2rad(pa))
   dr_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)

   # Get deflection vector in rotated frame
   ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (dr_r + sq))
   ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (dr_r + q^2 * sq))

   # Get potential
   ψ_up = ψ + dx_r * ψx_r + dy_r * ψy_r + bq * q * sq * log((1.0 + q) * sq / sqrt((dr_r + sq)^2 + (1.0 - q^2) * dx_r^2))
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   # Get axis-ratio
   q = 1.0 - ϵ

   # Get b_sie(q)
   bq = (4.0 * pi * (vd / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q)
   sq = θs / sqrt(q)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Coordinate in the rotated frame
         dx_r = + (θx[i, j] - θxc) * cos(deg2rad(pa)) + (θy[i, j] - θyc) * sin(deg2rad(pa))
         dy_r = - (θx[i, j] - θxc) * sin(deg2rad(pa)) + (θy[i, j] - θyc) * cos(deg2rad(pa))
         dr_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)

         # Get deflection vector in rotated frame
         ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (dr_r + sq))
         ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (dr_r + q^2 * sq))

         ψ[i, j] = ψ[i, j] + dx_r * ψx_r + dy_r * ψy_r + bq * q * sq * log((1.0 + q) * sq / sqrt((dr_r + sq)^2 + (1.0 - q^2) * dx_r^2))
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = 1.0 - ϵ

   # Get b_sie(q)
   bq = (4.0 * pi * (vd / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q)
   sq = θs / sqrt(q)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos(deg2rad(pa)) + (θy - θyc) * sin(deg2rad(pa))
   dy_r = - (θx - θxc) * sin(deg2rad(pa)) + (θy - θyc) * cos(deg2rad(pa))
   dr_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)

   # Get deflection vector in rotated frame
   ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (dr_r + sq))
   ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (dr_r + q^2 * sq))

   # Rotate back to original frame
   ψx_up = ψx + ψx_r * cos(deg2rad(pa)) - ψy_r * sin(deg2rad(pa))
   ψy_up = ψy + ψx_r * sin(deg2rad(pa)) + ψy_r * cos(deg2rad(pa))

   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   # Get axis-ratio
   q = 1.0 - ϵ

   # Get b_sie(q)
   bq = (4.0 * pi * (vd / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q)
   sq = θs / sqrt(q)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Coordinate in the rotated frame
         dx_r = + (θx[i, j] - θxc) * cos(deg2rad(pa)) + (θy[i, j] - θyc) * sin(deg2rad(pa))
         dy_r = - (θx[i, j] - θxc) * sin(deg2rad(pa)) + (θy[i, j] - θyc) * cos(deg2rad(pa))
         dr_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)

         # Get deflection vector in rotated frame
         ψx_r = (bq * q / sqrt(1 - q^2)) *  atan(sqrt(1 - q^2) * dx_r / (dr_r + sq))
         ψy_r = (bq * q / sqrt(1 - q^2)) * atanh(sqrt(1 - q^2) * dy_r / (dr_r + q^2 * sq))

         # Rotate back to original frame
         ψx[i, j] = ψx[i, j] + ψx_r * cos(deg2rad(pa)) - ψy_r * sin(deg2rad(pa))
         ψy[i, j] = ψy[i, j] + ψx_r * sin(deg2rad(pa)) + ψy_r * cos(deg2rad(pa))
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = 1.0 - ϵ

   # Get b_sie(q)
   bq = (4.0 * pi * (vd / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q)
   sq = θs / sqrt(q)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos(deg2rad(pa)) + (θy - θyc) * sin(deg2rad(pa))
   dy_r = - (θx - θxc) * sin(deg2rad(pa)) + (θy - θyc) * cos(deg2rad(pa))
   dr_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)

   # Common factor
   common_factor = (1+q^2) * sq^2 + 2 * dr_r * sq + dx_r^2 + dy_r^2

   # Get deformation tensor in rotated frame
   ψxx_r = + bq * q * (q^2 * sq^2 + dy_r^2 + sq * dr_r) / dr_r / common_factor
   ψyy_r = + bq * q * (sq^2 + dx_r^2 + sq * dr_r) / dr_r / common_factor
   ψxy_r = - bq * q * dx_r * dy_r / dr_r / common_factor

   # Rotate back to the original frame
   ψxx_up = ψxx + ψxx_r * cos(deg2rad(pa))^2 - ψxy_r * sin(deg2rad(2*pa)) + ψyy_r * sin(deg2rad(pa))^2
   ψyy_up = ψyy + ψxx_r * sin(deg2rad(pa))^2 + ψxy_r * sin(deg2rad(2*pa)) + ψyy_r * cos(deg2rad(pa))^2
   ψxy_up = ψxy + 0.5 * sin(deg2rad(2*pa)) * (ψxx_r - ψyy_r) + cos(deg2rad(2*pa)) * ψxy_r

   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   # Get axis-ratio
   q = 1.0 - ϵ

   # Get b_sie(q)
   bq = (4.0 * pi * (vd / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

   # Get s(q)
   sq = θs / sqrt(q)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Coordinate in the rotated frame
         dx_r = + (θx[i, j] - θxc) * cos(deg2rad(pa)) + (θy[i, j] - θyc) * sin(deg2rad(pa))
         dy_r = - (θx[i, j] - θxc) * sin(deg2rad(pa)) + (θy[i, j] - θyc) * cos(deg2rad(pa))
         dr_r = sqrt(q^2 * (sq^2 + dx_r^2) + dy_r^2)

         # Common factor
         common_factor = (1+q^2) * sq^2 + 2 * dr_r * sq + dx_r^2 + dy_r^2

         # Get deformation tensor in rotated frame
         ψxx_r = + bq * q * (q^2 * sq^2 + dy_r^2 + sq * dr_r) / dr_r / common_factor
         ψyy_r = + bq * q * (sq^2 + dx_r^2 + sq * dr_r) / dr_r / common_factor
         ψxy_r = - bq * q * dx_r * dy_r / dr_r / common_factor

         # Rotate back to the original frame
         ψxx[i, j] = ψxx[i, j] + ψxx_r * cos(deg2rad(pa))^2 - ψxy_r * sin(deg2rad(2*pa)) + ψyy_r * sin(deg2rad(pa))^2
         ψyy[i, j] = ψyy[i, j] + ψxx_r * sin(deg2rad(pa))^2 + ψxy_r * sin(deg2rad(2*pa)) + ψyy_r * cos(deg2rad(pa))^2
         ψxy[i, j] = ψxy[i, j] + 0.5 * sin(deg2rad(2*pa)) * (ψxx_r - ψyy_r) + cos(deg2rad(2*pa)) * ψxy_r
      end
   end
end

end