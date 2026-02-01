module eNFWMDLens

# Inbuilt Julia functions to import
using QuadGK

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!

@inline function F_x(x::RV)
   if x < 0.999999999
      if x > 1E-6 
         arg = sqrt( 1.0 - x^2 )
         return atanh( arg ) / arg
      else
         return log(2.0 / x)
      end
   elseif x > 1.000000001 
      arg = sqrt( x^2 - 1.0 )
      return atan( arg ) / arg
   else
      return 1.0
   end
end

@inline function α_r(x::RV)
   return (F_x(x) + log(0.5 * x)) / x
end

@inline function κ_r(x::RV)
   return 0.5 * (1.0 - F_x(x)) / (x^2 - 1.0)
end

@inline function κ_dr(x::RV)
   return 0.5 * (3.0 * x^2 * F_x(x) - 2.0 * x^2 - 1.0) / (x * (x^2 - 1.0)^2)
end

function I_integrand(u::RV, x::RV, y::RV, q::RV)
   ξ_u = sqrt(u * (x^2 + y^2 / (1.0 - (1.0 - q^2) * u)))
   return ξ_u * α_r(ξ_u) / (u * sqrt(1.0 - (1.0 - q^2) * u))
end

function I_integral(x::RV, y::RV, q::RV)
   I, _ = quadgk(u -> I_integrand(u, x, y, q), 0, 1)
   return I
end

function J_integrand(u::RV, x::RV, y::RV, q::RV, n::Int64)
   ξ_u = sqrt(u * (x^2 + y^2 / (1.0 - (1.0 - q^2) * u)))
   return κ_r(ξ_u) / (1.0 - (1.0 - q^2) * u)^(n + 0.5)
end

function J_integral(x::RV, y::RV, q::RV, n::Int64)
   J, _ = quadgk(u -> J_integrand(u, x, y, q, n), 0, 1)
   return J
end

function K_integrand(u::RV, x::RV, y::RV, q::RV, n::Int64)
   ξ_u = sqrt(u * (x^2 + y^2 / (1.0 - (1.0 - q^2) * u)))
   return 0.5 * u * κ_dr(ξ_u) / ξ_u / (1.0 - (1.0 - q^2) * u)^(n + 0.5)
end

function K_integral(x::RV, y::RV, q::RV, n::Int64)
   K, _ = quadgk(u -> K_integrand(u, x, y, q, n), 0, 1)
   return K
end

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)
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

   # Calculate integral
   I = I_integral(x, y, q)

   # Calculate potential
   ψ_up = ψ + κs * 0.5 * q * I

   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)
   κs = κs * θs_p^2

   # Precompute trigonometric functions
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Coordinate in the rotated frame
         dx_r = + (θx[i, j] - θxc) * cos_pa + (θy[i, j] - θyc) * sin_pa
         dy_r = - (θx[i, j] - θxc) * sin_pa + (θy[i, j] - θyc) * cos_pa

         # Scaled coordinates
         x = dx_r / θs_p
         y = dy_r / θs_p

         # Calculate integral
         I = I_integral(x, y, q)

         # Calculate potential
         ψ[i, j] += κs * 0.5 * q * I
      end
   end
   return nothing
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)
   κs = κs * θs_p

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

      # Calculate integrals
   J_0 = J_integral(x, y, q, 0)
   J_1 = J_integral(x, y, q, 1)

   # Calculate deflection vector in rotated frame
   ψx_r = κs * q * x * J_0
   ψy_r = κs * q * y * J_1

   # Get deflection vector in original frame
   ψx_up = ψx + ψx_r * cos_pa - ψy_r * sin_pa
   ψy_up = ψy + ψx_r * sin_pa + ψy_r * cos_pa

   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)
   κs = κs * θs_p

   # Precompute trigonometric functions
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Coordinate in the rotated frame
         dx_r = + (θx[i, j] - θxc) * cos_pa + (θy[i, j] - θyc) * sin_pa
         dy_r = - (θx[i, j] - θxc) * sin_pa + (θy[i, j] - θyc) * cos_pa

         # Scaled coordinates
         x = dx_r / θs_p
         y = dy_r / θs_p

         # Calculate integrals
         J_0 = J_integral(x, y, q, 0)
         J_1 = J_integral(x, y, q, 1)

         # Calculate deflection vector in rotated frame
         ψx_r = κs * q * x * J_0
         ψy_r = κs * q * y * J_1

         # Get deflection vector in original frame
         ψx[i, j] += ψx_r * cos_pa - ψy_r * sin_pa
         ψy[i, j] += ψx_r * sin_pa + ψy_r * cos_pa
      end
   end
   return nothing
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)

   # Precompute trigonometric functions
   pa_rad = deg2rad(pa)
   cos_pa = cos(pa_rad)
   sin_pa = sin(pa_rad)
   sin_2pa = sin(2.0 * pa_rad)
   cos_2pa = cos(2.0 * pa_rad)

   # Coordinate in the rotated frame
   dx_r = + (θx - θxc) * cos_pa + (θy - θyc) * sin_pa
   dy_r = - (θx - θxc) * sin_pa + (θy - θyc) * cos_pa

   # Scaled coordinates
   x = dx_r / θs_p
   y = dy_r / θs_p
   
   # Calculate integrals
   J_0 = J_integral(x, y, q, 0)
   J_1 = J_integral(x, y, q, 1)
   K_0 = K_integral(x, y, q, 0)
   K_1 = K_integral(x, y, q, 1)
   K_2 = K_integral(x, y, q, 2)

   # Calculate Jacobian in rotated frame
   ψxx_r = κs * q * (J_0 + 2.0 * x^2 * K_0)
   ψyy_r = κs * q * (J_1 + 2.0 * y^2 * K_2)
   ψxy_r = κs * q * (2.0 * x * y * K_1)

   # Get jacobian in original frame
   ψxx_up = ψxx + ψxx_r * cos_pa^2 - ψxy_r * sin_2pa + ψyy_r * sin_pa^2
   ψyy_up = ψyy + ψxx_r * sin_pa^2 + ψxy_r * sin_2pa + ψyy_r * cos_pa^2
   ψxy_up = ψxy + 0.5 * sin_2pa * (ψxx_r - ψyy_r) + cos_2pa * ψxy_r

   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)

   # Precompute trigonometric functions
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

         # Scaled coordinates
         x = dx_r / θs_p
         y = dy_r / θs_p

         # Calculate integrals
         J_0 = J_integral(x, y, q, 0)
         J_1 = J_integral(x, y, q, 1)
         K_0 = K_integral(x, y, q, 0)
         K_1 = K_integral(x, y, q, 1)
         K_2 = K_integral(x, y, q, 2)

         # Calculate Jacobian in rotated frame
         ψxx_r = κs * q * (J_0 + 2.0 * x^2 * K_0)
         ψyy_r = κs * q * (J_1 + 2.0 * y^2 * K_2)
         ψxy_r = κs * q * (2.0 * x * y * K_1)

         # Get jacobian in original frame
         ψxx[i, j] += ψxx_r * cos_pa^2 - ψxy_r * sin_2pa + ψyy_r * sin_pa^2
         ψyy[i, j] += ψxx_r * sin_pa^2 + ψxy_r * sin_2pa + ψyy_r * cos_pa^2
         ψxy[i, j] += 0.5 * sin_2pa * (ψxx_r - ψyy_r) + cos_2pa * ψxy_r
      end
   end
end

end