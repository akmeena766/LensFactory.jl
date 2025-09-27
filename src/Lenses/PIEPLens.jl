module PIEPLens

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
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   q = 1.0 - ϵ

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos(deg2rad(pa)) + (θy - θyc) * sin(deg2rad(pa))
   dy_r = - (θx - θxc) * sin(deg2rad(pa)) + (θy - θyc) * cos(deg2rad(pa))

   ψ_up = ψ + θE * sqrt(θs^2 + dx_r^2 + dy_r^2 / q^2)
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   q = 1.0 - ϵ

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx_r = + (θx[i, j] - θxc) * cos(deg2rad(pa)) + (θy[i, j] - θyc) * sin(deg2rad(pa))
         dy_r = - (θx[i, j] - θxc) * sin(deg2rad(pa)) + (θy[i, j] - θyc) * cos(deg2rad(pa))

         ψ[i, j] = ψ[i, j] + θE * sqrt(θs^2 + dx_r^2 + dy_r^2 / q^2)
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   q = 1.0 - ϵ

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos(deg2rad(pa)) + (θy - θyc) * sin(deg2rad(pa))
   dy_r = - (θx - θxc) * sin(deg2rad(pa)) + (θy - θyc) * cos(deg2rad(pa)) 
   dr_r = sqrt(θs^2 + dx_r^2 + dy_r^2 / q^2)

   # Deflection in the rotated frame
   ψx_r = θE * dx_r / dr_r
   ψy_r = θE * dy_r / dr_r / q^2

   # Rotate back to original frame
   ψx_up = ψx + ψx_r * cos(deg2rad(pa)) - ψy_r * sin(deg2rad(pa))
   ψy_up = ψy + ψx_r * sin(deg2rad(pa)) + ψy_r * cos(deg2rad(pa))

   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   q = 1.0 - ϵ

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Coordinate in the rotated frame
         dx_r = + (θx[i, j] - θxc) * cos(deg2rad(pa)) + (θy[i, j] - θyc) * sin(deg2rad(pa))
         dy_r = - (θx[i, j] - θxc) * sin(deg2rad(pa)) + (θy[i, j] - θyc) * cos(deg2rad(pa))
         dr_r = sqrt(θs^2 + dx_r^2 + dy_r^2 / q^2)

         # Deflection in the rotated frame
         ψx_r = θE * dx_r / dr_r
         ψy_r = θE * dy_r / dr_r / q^2

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
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   q = 1.0 - ϵ

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos(deg2rad(pa)) + (θy - θyc) * sin(deg2rad(pa))
   dy_r = - (θx - θxc) * sin(deg2rad(pa)) + (θy - θyc) * cos(deg2rad(pa))
   dr_r = sqrt(θs^2 + dx_r^2 + dy_r^2 / q^2)

   # Deformation tensor components in rotated frame
   ψxx_r = + θE * (θs^2 + dy_r^2 / q^2) / dr_r^3
   ψyy_r = + θE * (θs^2 + dx_r^2) / dr_r^3 / q^2
   ψxy_r = - θE * dx_r * dy_r / dr_r^3 / q^2

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
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   q = 1.0 - ϵ

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Coordinate in the rotated frame
         dx_r = (θx[i, j] - θxc) * cos(deg2rad(pa)) + (θy[i, j] - θyc) * sin(deg2rad(pa))
         dy_r = (θy[i, j] - θyc) * cos(deg2rad(pa)) - (θx[i, j] - θxc) * sin(deg2rad(pa))
         dr_r = sqrt(θs^2 + dx_r^2 + dy_r^2 / q^2)

         # Deformation tensor components in rotated frame
         ψxx_r = + θE * (θs^2 + dy_r^2 / q^2) / dr_r^3
         ψyy_r = + θE * (θs^2 + dx_r^2) / dr_r^3 / q^2
         ψxy_r = - θE * dx_r * dy_r / dr_r^3 / q^2

         # Rotate back to the original frame
         ψxx[i, j] = ψxx[i, j] + ψxx_r * cos(deg2rad(pa))^2 - ψxy_r * sin(deg2rad(2*pa)) + ψyy_r * sin(deg2rad(pa))^2
         ψyy[i, j] = ψyy[i, j] + ψxx_r * sin(deg2rad(pa))^2 + ψxy_r * sin(deg2rad(2*pa)) + ψyy_r * cos(deg2rad(pa))^2
         ψxy[i, j] = ψxy[i, j] + 0.5 * sin(deg2rad(2*pa)) * (ψxx_r - ψyy_r) + cos(deg2rad(2*pa)) * ψxy_r
      end
   end
end

end