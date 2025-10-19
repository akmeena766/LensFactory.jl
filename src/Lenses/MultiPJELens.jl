module MultiPJELens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!


function potential!(ψ::T, θx::T, θy::T, θxc::S, θyc::S, vd::S, θs::S, θt::S, ϵ::S, pa::S, nl::Int64) where {T <: RV, S <: Vector{<:RV}}
   ψ_up = ψ
   for k in 1:nl
      # Get axis-ratio
      q = 1.0 - ϵ[k]

      # Get b_sie(q)
      bq = (4.0 * pi * (vd[k] / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

      # Get s(q) and a(q)
      sq = θs[k] / sqrt(q)
      aq = θt[k] / sqrt(q)

      # Coordinate in the rotated frame
      dx_r = + (θx - θxc[k]) * cos(deg2rad(pa[k])) + (θy - θyc[k]) * sin(deg2rad(pa[k]))
      dy_r = - (θx - θxc[k]) * sin(deg2rad(pa[k])) + (θy - θyc[k]) * cos(deg2rad(pa[k]))
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

      ψ_up = ψ_up + (ψ1 - ψ2)
   end
   return ψ_up
end

function potential!(ψ::T, θx::T, θy::T, θxc::S, θyc::S, vd::S, θs::S, θt::S, ϵ::S, pa::S, nl::Int64) where {T <: ROA, S <: Vector{<:RV}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      # Get axis-ratio
      q = 1.0 - ϵ[k]

      # Get b_sie(q)
      bq = (4.0 * pi * (vd[k] / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

      # Get s(q) and a(q)
      sq = θs[k] / sqrt(q)
      aq = θt[k] / sqrt(q)

      @inbounds for j in ax2
         @inbounds for i in ax1
            # Coordinate in the rotated frame
            dx_r = + (θx[i, j] - θxc[k]) * cos(deg2rad(pa[k])) + (θy[i, j] - θyc[k]) * sin(deg2rad(pa[k]))
            dy_r = - (θx[i, j] - θxc[k]) * sin(deg2rad(pa[k])) + (θy[i, j] - θyc[k]) * cos(deg2rad(pa[k]))
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
end


function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::S, θyc::S, vd::S, θs::S, θt::S, ϵ::S, pa::S, nl::Int64) where {T <: RV, S <: Vector{<:RV}}
   ψx_up = ψx
   ψy_up = ψy
   for k in 1:nl
      # Get axis-ratio
      q = 1.0 - ϵ[k]

      # Get b_sie(q)
      bq = (4.0 * pi * (vd[k] / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

      # Get s(q) and a(q)
      sq = θs[k] / sqrt(q)
      aq = θt[k] / sqrt(q)

      # Coordinate in the rotated frame
      dx_r = + (θx - θxc[k]) * cos(deg2rad(pa[k])) + (θy - θyc[k]) * sin(deg2rad(pa[k]))
      dy_r = - (θx - θxc[k]) * sin(deg2rad(pa[k])) + (θy - θyc[k]) * cos(deg2rad(pa[k]))
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
      ψx_up = ψx_up + ψx_r * cos(deg2rad(pa[k])) - ψy_r * sin(deg2rad(pa[k]))
      ψy_up = ψy_up + ψx_r * sin(deg2rad(pa[k])) + ψy_r * cos(deg2rad(pa[k]))
   end
   return ψx_up, ψy_up
end


function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::S, θyc::S, vd::S, θs::S, θt::S, ϵ::S, pa::S, nl::Int64) where {T <: ROA, S <: Vector{<:RV}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      # Get axis-ratio
      q = 1.0 - ϵ[k]

      # Get b_sie(q)
      bq = (4.0 * pi * (vd[k] / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

      # Get s(q) and a(q)
      sq = θs[k] / sqrt(q)
      aq = θt[k] / sqrt(q)

      @inbounds for j in ax2
         @inbounds for i in ax1
            # Coordinate in the rotated frame
            dx_r = + (θx[i, j] - θxc[k]) * cos(deg2rad(pa[k])) + (θy[i, j] - θyc[k]) * sin(deg2rad(pa[k]))
            dy_r = - (θx[i, j] - θxc[k]) * sin(deg2rad(pa[k])) + (θy[i, j] - θyc[k]) * cos(deg2rad(pa[k]))
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
            ψx[i, j] = ψx[i, j] + ψx_r * cos(deg2rad(pa[k])) - ψy_r * sin(deg2rad(pa[k]))
            ψy[i, j] = ψy[i, j] + ψx_r * sin(deg2rad(pa[k])) + ψy_r * cos(deg2rad(pa[k]))
         end
      end
   end
end


function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::S, θyc::S, vd::S, θs::S, θt::S, ϵ::S, pa::S, nl::Int64) where {T <: RV, S <: Vector{<:RV}}
   ψxx_up = ψxx
   ψyy_up = ψyy
   ψxy_up = ψxy
   for k in 1:nl
      # Get axis-ratio
      q = 1.0 - ϵ[k]

      # Get b_sie(q)
      bq = (4.0 * pi * (vd[k] / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

      # Get s(q) and a(q)
      sq = θs[k] / sqrt(q)
      aq = θt[k] / sqrt(q)

      # Coordinate in the rotated frame
      dx_r = + (θx - θxc[k]) * cos(deg2rad(pa[k])) + (θy - θyc[k]) * sin(deg2rad(pa[k]))
      dy_r = - (θx - θxc[k]) * sin(deg2rad(pa[k])) + (θy - θyc[k]) * cos(deg2rad(pa[k]))
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
      ψxx_up = ψxx_up + ψxx_r * cos(deg2rad(pa[k]))^2 - ψxy_r * sin(deg2rad(2*pa[k])) + ψyy_r * sin(deg2rad(pa[k]))^2
      ψyy_up = ψyy_up + ψxx_r * sin(deg2rad(pa[k]))^2 + ψxy_r * sin(deg2rad(2*pa[k])) + ψyy_r * cos(deg2rad(pa[k]))^2
      ψxy_up = ψxy_up  + 0.5 * sin(deg2rad(2*pa[k])) * (ψxx_r - ψyy_r) + cos(deg2rad(2*pa[k])) * ψxy_r
   end
   return ψxx_up, ψyy_up, ψxy_up
end

function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::S, θyc::S, vd::S, θs::S, θt::S, ϵ::S, pa::S, nl::Int64) where {T <: ROA, S <: Vector{<:RV}}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   for k in 1:nl
      # Get axis-ratio
      q = 1.0 - ϵ[k]

      # Get b_sie(q)
      bq = (4.0 * pi * (vd[k] / CONST_C)^2) / sqrt(q) / ANGLE_ARCSEC

      # Get s(q) and a(q)
      sq = θs[k] / sqrt(q)
      aq = θt[k] / sqrt(q)

      @inbounds for j in ax2
         @inbounds for i in ax1
            # Coordinate in the rotated frame
            dx_r = + (θx[i, j] - θxc[k]) * cos(deg2rad(pa[k])) + (θy[i, j] - θyc[k]) * sin(deg2rad(pa[k]))
            dy_r = - (θx[i, j] - θxc[k]) * sin(deg2rad(pa[k])) + (θy[i, j] - θyc[k]) * cos(deg2rad(pa[k]))
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
            ψxx[i, j] = ψxx[i, j] + ψxx_r * cos(deg2rad(pa[k]))^2 - ψxy_r * sin(deg2rad(2*pa[k])) + ψyy_r * sin(deg2rad(pa[k]))^2
            ψyy[i, j] = ψyy[i, j] + ψxx_r * sin(deg2rad(pa[k]))^2 + ψxy_r * sin(deg2rad(2*pa[k])) + ψyy_r * cos(deg2rad(pa[k]))^2
            ψxy[i, j] = ψxy[i, j] + 0.5 * sin(deg2rad(2*pa[k])) * (ψxx_r - ψyy_r) + cos(deg2rad(2*pa[k])) * ψxy_r
         end
      end
   end
end

end