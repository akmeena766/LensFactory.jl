module aNFWLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!


@inline function ϕ_CSE(x::Real, y::Real, x_s::Real, q::Real)
   # Axis ratio squared
   q2 = q^2

   # Small psi
   ψ = sqrt(q2 * (x_s^2 + x^2) + y^2)

   # Capital Psi
   Ψ = (ψ + x_s)^2 + (1.0 - q2) * x^2

   return 0.5 * (q / x_s) * log(Ψ) - (q / x_s) * log((1.0 + q) * x_s)
end

@inline function ϕx_CSE(x::Real, y::Real, x_s::Real, q::Real)
   # Axis ratio squared
   q2 = q^2

   # Small psi
   ψ = sqrt(q2 * (x_s^2 + x^2) + y^2)

   # Capital Psi
   Ψ = (ψ + x_s)^2 + (1.0 - q2) * x^2   

   return (q / x_s) * (x / ψ) * (ψ + q^2 * x_s) / Ψ
end

@inline function ϕy_CSE(x::Real, y::Real, x_s::Real, q::Real)
   # Axis ratio squared
   q2 = q^2

   # Small psi
   ψ = sqrt(q2 * (x_s^2 + x^2) + y^2)

   # Capital Psi
   Ψ = (ψ + x_s)^2 + (1.0 - q2) * x^2   

   return (q / x_s) * (y / ψ) * (ψ + x_s) / Ψ
end

@inline function ϕxx_CSE(x::Real, y::Real, x_s::Real, q::Real)
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

@inline function ϕyy_CSE(x::Real, y::Real, x_s::Real, q::Real)
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

@inline function ϕxy_CSE(x::Real, y::Real, x_s::Real, q::Real)
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
    potential!(ψ::Real, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, ϵ::T, pa::T) where {S<:Real, T<:Real}
"""
function potential!(ψ::Real, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, ϵ::T, pa::T) where {S<:Real, T<:Real}
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = 4.0 * k_s
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
   for k in 1:44
      ψ_r = ψ_r + Ai[k] * ϕ_CSE(x, y, Si[k], q)
   end
   ψ_up = ψ + κs * ψ_r
   
   return ψ_up
end

"""
    potential!(ψ::ROA, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, ϵ::T, pa::T) where {S<:ROA, T<:Real}
"""
function potential!(ψ::ROA, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, ϵ::T, pa::T) where {S<:ROA, T<:Real}
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = 4.0 * k_s
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
         for k in 1:44
            ψ_r = ψ_r + Ai[k] * ϕ_CSE(x, y, Si[k], q)
         end
         ψ[i, j] = ψ[i, j] + κs * ψ_r
      end
   end
   return nothing
end

"""
    deflection!(ψx::Real, ψy::Real, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, ϵ::T, pa::T) where {S<:Real, T<:Real}
"""
function deflection!(ψx::Real, ψy::Real, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, ϵ::T, pa::T) where {S<:Real, T<:Real}
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = 4.0 * k_s
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
   for k in 1:44
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
    deflection!(ψx::ROA, ψy::ROA, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, ϵ::T, pa::T) where {S<:ROA, T<:Real}
"""
function deflection!(ψx::ROA, ψy::ROA, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, ϵ::T, pa::T) where {S<:ROA, T<:Real}
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = 4.0 * k_s
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
         for k in 1:44
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
    jacobian!(ψxx::Real, ψyy::Real, ψxy::Real, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, ϵ::T, pa::T) where {S<:Real, T<:Real}
"""
function jacobian!(ψxx::Real, ψyy::Real, ψxy::Real, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, ϵ::T, pa::T) where {S<:Real, T<:Real}
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = 4.0 * k_s

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
   for k in 1:44
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
    jacobian!(ψxx::ROA, ψyy::ROA, ψxy::ROA, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, ϵ::T, pa::T) where {S<:ROA, T<:Real}
"""
function jacobian!(ψxx::ROA, ψyy::ROA, ψxy::ROA, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, ϵ::T, pa::T) where {S<:ROA, T<:Real}
   # Get axis-ratio
   q = (1.0 - ϵ) / (1.0 + ϵ)
   θs_p = θs / sqrt(q)

   # Get normalization constant κs
   κs = 4.0 * k_s

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
         for k in 1:44
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


const Ai = [1.648988e-18,
            6.274458e-16,
            3.646620e-17,
            3.459206e-15,
            2.457389e-14,
            1.059319e-13,
            4.211597e-13,
            1.142832e-12,
            4.391215e-12,
            1.556500e-11,
            6.951271e-11,
            3.147466e-10,
            1.379109e-09,
            3.829778e-09,
            1.384858e-08,
            5.370951e-08,
            1.804384e-07,
            5.788608e-07,
            3.205256e-06,
            1.102422e-05,
            4.093971e-05,
            1.282206e-04,
            4.575541e-04,
            7.995270e-04,
            5.013701e-03,
            1.403508e-02,
            5.230727e-02,
            1.898907e-01,
            3.643448e-01,
            7.203734e-01,
            1.717667e00,
            2.217566e00,
            3.187447e00,
            8.194898e00,
            1.765210e01,
            1.974319e01,
            2.783688e01,
            4.482311e01,
            5.598897e01,
            1.426485e02,
            2.279833e02,
            5.401335e02,
            9.743682e02,
            1.775124e03]


const Si = [1.082411e-06,
            8.786566e-06,
            3.292868e-06,
            1.860019e-05,
            3.274231e-05,
            6.232485e-05,
            9.256333e-05,
            1.546762e-04,
            2.097321e-04,
            3.391140e-04,
            5.178790e-04,
            8.636736e-04,
            1.405152e-03,
            2.193855e-03,
            3.179572e-03,
            4.970987e-03,
            7.631970e-03,
            1.119413e-02,
            1.827267e-02,
            2.945251e-02,
            4.562723e-02,
            6.782509e-02,
            1.596987e-01,
            1.127751e-01,
            2.169469e-01,
            3.423835e-01,
            5.194527e-01,
            8.623185e-01,
            1.382737e00,
            2.034929e00,
            3.402979e00,
            5.594276e00,
            8.052345e00,
            1.349045e01,
            2.603825e01,
            4.736823e01,
            6.559320e01,
            1.087932e02,
            1.477673e02,
            2.495341e02,
            4.305999e02,
            7.760206e02,
            2.143057e03,
            1.935749e03]
end