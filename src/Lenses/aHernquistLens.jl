module aHernquistLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!


@inline function ϕ_CSE(x::RV, y::RV, x_s::RV, q::RV)
   # Axis ratio squared
   q2 = q^2

   # Small psi
   ψ = sqrt(q2 * (x_s^2 + x^2) + y^2)

   # Capital Psi
   Ψ = (ψ + x_s)^2 + (1.0 - q2) * x^2

   return (0.5 * q / x_s) * log(Ψ) + (q / x_s) * log((1.0 + q) * x_s)
end

@inline function ϕx_CSE(x::RV, y::RV, x_s::RV, q::RV)
   # Axis ratio squared
   q2 = q^2

   # Small psi
   ψ = sqrt(q2 * (x_s^2 + x^2) + y^2)

   # Capital Psi
   Ψ = (ψ + x_s)^2 + (1.0 - q2) * x^2

   return (q / x_s) * (x / ψ) * (ψ + q^2 * x_s) / Ψ
end

@inline function ϕy_CSE(x::RV, y::RV, x_s::RV, q::RV)
   # Axis ratio squared
   q2 = q^2

   # Small psi
   ψ = sqrt(q2 * (x_s^2 + x^2) + y^2)

   # Capital Psi
   Ψ = (ψ + x_s)^2 + (1.0 - q2) * x^2

   return (q / x_s) * (y / ψ) * (ψ + x_s) / Ψ
end

@inline function ϕxx_CSE(x::RV, y::RV, x_s::RV, q::RV)
   # Axis ratio squared
   q2 = q^2

   # Small psi
   ψ = sqrt(q2 * (x_s^2 + x^2) + y^2)

   # Capital Psi
   Ψ = (ψ + x_s)^2 + (1.0 - q2) * x^2

   inv_Ψ = 1.0 / Ψ
   inv_ψ = 1.0 / ψ

   return (q / x_s) * inv_Ψ * (1.0 + q2 * x_s * (q2 * x_s^2 + y^2) * inv_ψ^3 - 2.0 * x^2 * (ψ + q2 * x_s)^2 * inv_ψ^2 * inv_Ψ)
end

@inline function ϕyy_CSE(x::RV, y::RV, x_s::RV, q::RV)
   # Axis ratio squared
   q2 = q^2

   # Small psi
   ψ = sqrt(q2 * (x_s^2 + x^2) + y^2)

   # Capital Psi
   Ψ = (ψ + x_s)^2 + (1.0 - q2) * x^2

   inv_Ψ = 1.0 / Ψ
   inv_ψ = 1.0 / ψ

   return (q / x_s) * inv_Ψ * (1.0 + q2 * x_s * (x_s^2 + x^2) * inv_ψ^3 - 2.0 * y^2 * (ψ + x_s)^2 * inv_ψ^2 * inv_Ψ)
end

@inline function ϕxy_CSE(x::RV, y::RV, x_s::RV, q::RV)
   # Axis ratio squared
   q2 = q^2

   # Small psi
   ψ = sqrt(q2 * (x_s^2 + x^2) + y^2)

   # Capital Psi
   Ψ = (ψ + x_s)^2 + (1.0 - q2) * x^2

   inv_Ψ = 1.0 / Ψ
   inv_ψ = 1.0 / ψ

   return - (q * x * y * inv_Ψ / x_s) * (q2 * x_s * inv_ψ^3 + 2.0 * (ψ + q2 * x_s) * (ψ + x_s) * inv_ψ^2 * inv_Ψ)
end


"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
"""
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

   ψ_r = 0.0
   for k in 1:41
      ψ_r = ψ_r + Ai[k] * ϕ_CSE(x, y, Si[k], q)
   end
   ψ_up = ψ + κs * ψ_r

   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
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

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Coordinate in the rotated frame
         dx_r = + (θx[i, j] - θxc) * cos_pa + (θy[i, j] - θyc) * sin_pa
         dy_r = - (θx[i, j] - θxc) * sin_pa + (θy[i, j] - θyc) * cos_pa

         # Scaled coordinates
         x = dx_r / θs_p
         y = dy_r / θs_p

         ψ_r = 0.0
         for k in 1:41
            ψ_r = ψ_r + Ai[k] * ϕ_CSE(x, y, Si[k], q)
         end
         ψ[i, j] = ψ[i, j] + κs * ψ_r
      end
   end
   return nothing
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = 1.0 - ϵ
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
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
   
   ψx_r = 0.0
   ψy_r = 0.0
   for k in 1:41
      ψx_r = ψx_r + Ai[k] * ϕx_CSE(x, y, Si[k], q)
      ψy_r = ψy_r + Ai[k] * ϕy_CSE(x, y, Si[k], q)
   end
   ψx_r = κs * ψx_r
   ψy_r = κs * ψy_r

   # Get deflection vector in original frame
   ψx_up = ψx + ψx_r * cos_pa - ψy_r * sin_pa
   ψy_up = ψy + ψx_r * sin_pa + ψy_r * cos_pa

   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   # Get axis-ratio
   q = 1.0 - ϵ
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
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

         ψx_r = 0.0
         ψy_r = 0.0
         for k in 1:41
            ψx_r = ψx_r + Ai[k] * ϕx_CSE(x, y, Si[k], q)
            ψy_r = ψy_r + Ai[k] * ϕy_CSE(x, y, Si[k], q)
         end
         ψx_r = κs * ψx_r
         ψy_r = κs * ψy_r

         # Get deflection vector in original frame
         ψx[i, j] = ψx[i, j] + ψx_r * cos_pa - ψy_r * sin_pa
         ψy[i, j] = ψy[i, j] + ψx_r * sin_pa + ψy_r * cos_pa
      end
   end
   return nothing
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: RV
   # Get axis-ratio
   q = 1.0 - ϵ
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)

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

   ψxx_r = 0.0
   ψyy_r = 0.0
   ψxy_r = 0.0
   for k in 1:41
      ψxx_r = ψxx_r + Ai[k] * ϕxx_CSE(x, y, Si[k], q)
      ψyy_r = ψyy_r + Ai[k] * ϕyy_CSE(x, y, Si[k], q)
      ψxy_r = ψxy_r + Ai[k] * ϕxy_CSE(x, y, Si[k], q)
   end
   ψxx_r = κs * ψxx_r
   ψyy_r = κs * ψyy_r
   ψxy_r = κs * ψxy_r

   # Get jacobian in original frame
   ψxx_up = ψxx + ψxx_r * cos_pa^2 - ψxy_r * sin_2pa + ψyy_r * sin_pa^2
   ψyy_up = ψyy + ψxx_r * sin_pa^2 + ψxy_r * sin_2pa + ψyy_r * cos_pa^2
   ψxy_up = ψxy + 0.5 * sin_2pa * (ψxx_r - ψyy_r) + cos_2pa * ψxy_r

   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass:: RV, θs::RV, ϵ::RV, pa::RV) where T <: ROA
   # Get axis-ratio
   q = 1.0 - ϵ
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)

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

         ψxx_r = 0.0
         ψyy_r = 0.0
         ψxy_r = 0.0
         for k in 1:41
            ψxx_r = ψxx_r + Ai[k] * ϕxx_CSE(x, y, Si[k], q)
            ψyy_r = ψyy_r + Ai[k] * ϕyy_CSE(x, y, Si[k], q)
            ψxy_r = ψxy_r + Ai[k] * ϕxy_CSE(x, y, Si[k], q)
         end
         ψxx_r = κs * ψxx_r
         ψyy_r = κs * ψyy_r
         ψxy_r = κs * ψxy_r

         # Get jacobian in original frame
         ψxx[i, j] = ψxx[i, j] + ψxx_r * cos_pa^2 - ψxy_r * sin_2pa + ψyy_r * sin_pa^2
         ψyy[i, j] = ψyy[i, j] + ψxx_r * sin_pa^2 + ψxy_r * sin_2pa + ψyy_r * cos_pa^2
         ψxy[i, j] = ψxy[i, j] + 0.5 * sin_2pa * (ψxx_r - ψyy_r) + cos_2pa * ψxy_r
      end
   end
   return nothing
end


const Ai = [9.200445e-18,
            2.184724e-16,
            3.548079e-15,
            2.823716e-14,
            1.091876e-13,
            6.998697e-13,
            3.142264e-12,
            1.457280e-11,
            4.472783e-11,
            2.042079e-10,
            8.708137e-10,
            2.423649e-09,
            7.353440e-09,
            5.470738e-08,
            2.445878e-07,
            4.541672e-07,
            3.227611e-06,
            1.110690e-05,
            3.725101e-05,
            1.056271e-04,
            6.531501e-04,
            2.121330e-03,
            8.285518e-03,
            4.084190e-02,
            5.760942e-02,
            1.788945e-01,
            2.092774e-01,
            3.697750e-01,
            3.440555e-01,
            5.792737e-01,
            2.325935e-01,
            5.227961e-01,
            3.079968e-01,
            1.633456e-01,
            7.410900e-02,
            3.123329e-02,
            1.292488e-02,
            2.156527e00,
            1.652553e-02,
            2.314934e-02,
            3.992313e-01]


const Si = [1.199110e-06,
            3.751762e-06,
            9.927207e-06,
            2.206076e-05,
            3.781528e-05,
            6.659808e-05,
            1.154366e-04,
            1.924150e-04,
            3.040440e-04,
            4.683051e-04,
            7.745084e-04,
            1.175953e-03,
            1.675459e-03,
            2.801948e-03,
            9.712807e-03,
            5.469589e-03,
            1.104654e-02,
            1.893893e-02,
            2.792864e-02,
            4.152834e-02,
            6.640398e-02,
            1.107083e-01,
            1.648028e-01,
            2.839601e-01,
            4.129439e-01,
            8.239115e-01,
            6.031726e-01,
            1.145604e00,
            1.401895e00,
            2.512223e00,
            2.038025e00,
            4.644014e00,
            9.301590e00,
            2.039273e01,
            4.896534e01,
            1.252311e02,
            3.576766e02,
            2.579464e04,
            2.944679e04,
            2.834717e03,
            5.931328e04]
end