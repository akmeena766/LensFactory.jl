module Multipole


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

"""
    potential!(ψ::U, θx::S, θy::S, ϵ::T, θϵ::T, m::Int64, n::T) where {U<:Real, S<:Real, T<:Real}
"""
function potential!(ψ::U, θx::S, θy::S, ϵ::T, θϵ::T, m::Int64, n::T) where {U<:Real, S<:Real, T<:Real}
   # Get radial coordiate
   dr = sqrt(θx^2 + θy^2)

   # Handle origin
   if dr == 0.0
      return 0.0
   end

   # Calculate θ
   θ = atan(θy, θx)

   # Overall phase
   ϕ = m * (θ - deg2rad(θϵ))

   ψ_up = ψ + (ϵ / m) * dr^n * cos(ϕ)
   return ψ_up
end

"""
    potential!(ψ::U, θx::S, θy::S, ϵ::T, θϵ::T, m::Int64, n::T) where {U<:ROA, S<:ROA, T<:Real}
"""
function potential!(ψ::U, θx::S, θy::S, ϵ::T, θϵ::T, m::Int64, n::T) where {U<:ROA, S<:ROA, T<:Real}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j]
         dy = θy[i, j]
         dr = sqrt(θx[i, j]^2 + θy[i, j]^2)

         # Handle origin
         if dr == 0.0
            ψ[i, j] = ψ[i, j] + 0.0
            continue
         end

         # Calculate θ
         θ = atan(dy, dx)

         # Overall phase
         ϕ = m * (θ - deg2rad(θϵ))

         ψ[i, j] = ψ[i, j] + (ϵ / m) * dr^n * cos(ϕ)
      end
   end
end

"""
    deflection!(ψx::U, ψy::U, θx::S, θy::S, ϵ::T, θϵ::T, m::Int64, n::T) where {U<:Real, S<:Real, T<:Real}
"""
function deflection!(ψx::U, ψy::U, θx::S, θy::S, ϵ::T, θϵ::T, m::Int64, n::T) where {U<:Real, S<:Real, T<:Real}
   # Get radial coordiate
   dr = sqrt(θx^2 + θy^2)

   # Handle origin
   if dr == 0.0
      return 0.0, 0.0
   end

   # Calculate θ
   θ = atan(θy, θx)

   # Overall phase
   ϕ = m * (θ - deg2rad(θϵ))

   pre_factor = (ϵ / m) * dr^(n - 2.0)
   ψx_up = ψx + pre_factor * (n * θx * cos(ϕ) + m * θy * sin(ϕ))
   ψy_up = ψy + pre_factor * (n * θy * cos(ϕ) - m * θx * sin(ϕ))
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::U, ψy::U, θx::S, θy::S, ϵ::T, θϵ::T, m::Int64, n::T) where {U<:ROA, S<:ROA, T<:Real}
"""
function deflection!(ψx::U, ψy::U, θx::S, θy::S, ϵ::T, θϵ::T, m::Int64, n::T) where {U<:ROA, S<:ROA, T<:Real}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j]
         dy = θy[i, j]
         dr = sqrt(θx[i, j]^2 + θy[i, j]^2)

         # Handle origin
         if dr == 0.0
            ψx[i, j] = ψx[i, j] + 0.0
            ψy[i, j] = ψy[i, j] + 0.0
            continue
         end

         # Calculate θ
         θ = atan(dy, dx)

         # Overall phase
         ϕ = m * (θ - deg2rad(θϵ))
         cos_phi = cos(ϕ)
         sin_phi = sin(ϕ)

         pre_factor = (ϵ / m) * dr^(n - 2.0)
         ψx[i, j] = ψx[i, j] + pre_factor * (n * dx * cos_phi + m * dy * sin_phi)
         ψy[i, j] = ψy[i, j] + pre_factor * (n * dy * cos_phi - m * dx * sin_phi)
      end
   end
end


"""
    jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, ϵ::T, θϵ::T, m::Int64, n::T) where {U<:Real, S<:Real, T<:Real}
"""
function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, ϵ::T, θϵ::T, m::Int64, n::T) where {U<:Real, S<:Real, T<:Real}
   # Get radial coordiate
   dr = sqrt(θx^2 + θy^2)

   # Handle origin
   if dr == 0.0
      return 0.0, 0.0, 0.0
   end

   # Calculate θ
   θ = atan(θy, θx)

   # Overall phase
   ϕ = m * (θ - deg2rad(θϵ))
   cos_phi = cos(ϕ)
   sin_phi = sin(ϕ)
   
   ψxx_tmp = (n - 2.0) * θx * (n * θx * cos_phi + m * θy * sin_phi) + (n * dr^2 - m^2 * θy^2) * cos_phi + m * n * θx * θy * sin_phi
   ψyy_tmp = (n - 2.0) * θy * (n * θy * cos_phi - m * θx * sin_phi) + (n * dr^2 - m^2 * θx^2) * cos_phi - m * n * θx * θy * sin_phi
   ψxy_tmp = ((n - 2.0) * n + m^2) * θx * θy * cos_phi + m * (n - 1.0) * (θy^2 - θx^2) * sin_phi

   pre_factor = (ϵ / m) * dr^(n - 4.0)
   ψxx_up = ψxx + pre_factor * ψxx_tmp
   ψyy_up = ψyy + pre_factor * ψyy_tmp
   ψxy_up = ψxy + pre_factor * ψxy_tmp
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, ϵ::T, θϵ::T, m::Int64, n::T) where {U<:ROA, S<:ROA, T<:Real}
"""
function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, ϵ::T, θϵ::T, m::Int64, n::T) where {U<:ROA, S<:ROA, T<:Real}
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j]
         dy = θy[i, j]
         dr = sqrt(θx[i, j]^2 + θy[i, j]^2)

         # Handle origin
         if dr == 0.0
            ψxx[i, j] = ψxx[i, j] + 0.0
            ψyy[i, j] = ψyy[i, j] + 0.0
            ψxy[i, j] = ψxy[i, j] + 0.0
            continue
         end

         # Calculate θ
         θ = atan(dy, dx)

         # Overall phase
         ϕ = m * (θ - deg2rad(θϵ))
         cos_phi = cos(ϕ)
         sin_phi = sin(ϕ)

         ψxx_tmp = (n - 2.0) * dx * (n * dx * cos_phi + m * dy * sin_phi) + (n * dr^2 - m^2 * dy^2) * cos_phi + m * n * dx * dy * sin_phi
         ψyy_tmp = (n - 2.0) * dy * (n * dy * cos_phi - m * dx * sin_phi) + (n * dr^2 - m^2 * dx^2) * cos_phi - m * n * dx * dy * sin_phi
         ψxy_tmp = ((n - 2.0) * n + m^2) * dx * dy * cos_phi + m * (n - 1.0) * (dy^2 - dx^2) * sin_phi

         pre_factor = (ϵ / m) * dr^(n - 4.0)
         ψxx[i, j] = ψxx[i, j] + pre_factor * ψxx_tmp
         ψyy[i, j] = ψyy[i, j] + pre_factor * ψyy_tmp
         ψxy[i, j] = ψxy[i, j] + pre_factor * ψxy_tmp
      end
   end
end

end