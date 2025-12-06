module eHernquistMDLens

# Inbuilt Julia functions to import
using QuadGK

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!

@inline function F_x(x::RV)
   if x < 1.0 
      arg = sqrt( 1.0 - x^2 )
      f_x = atanh( arg ) / arg
   else 
      arg = sqrt( x^2 - 1.0 )
      f_x =  atan( arg ) / arg
   end
   return f_x
end


function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
      # Get axis-ratio
   q = 1.0 - ϵ
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   κs = κs * θs_p^2

   # Precompute trigonometric functions
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos_pa + (θy - θyc) * sin_pa
   dy_r = - (θx - θxc) * sin_pa + (θy - θyc) * cos_pa

   # Scaled coordinates
   x = dx_r / θs_p
   y = dy_r / θs_p

   
end

function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   
end

function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: RV

end

function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA

end

function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
   
end

function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   
end

end