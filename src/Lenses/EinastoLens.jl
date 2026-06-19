module EinastoLens

# Inbuilt Julia functions to import
using SpecialFunctions
using QuadGK

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!


function m_Ein(θ::Real, n::Real)
   Pax, _ = gamma_inc(3.0 / n, (2.0 / n) * θ^n)
   return (1.0 / n) * (0.5 * n)^(3.0 / n) * gamma(3.0 / n) * Pax
end

@inline function I_κ(z::Real, θ::Real, n::Real)
   return exp(-(2.0 / n) * (θ^2 + z^2)^(0.5 * n))
end

function κ(θ::Real, n::Real)
   i_value, _  = quadgk(x -> I_κ(x, θ, n), 0, Inf)
   return 0.5 * i_value
end

@inline function I_α(z::Real, θ::Real, n::Real)
   Pax, _ = gamma_inc(3.0 / n, (2.0 / n) * θ^n * (1.0 + z^2)^(0.5 * n))
   return  Pax / (1.0 + z^2)^1.5
end

function α(θ::Real, n::Real)
   i_value, _ = quadgk(x -> I_α(x, θ, n), 0, Inf)
   return gamma(3.0 / n) * i_value * (0.5 * n)^(3.0 / n) / θ / n
end

function ϕ(θ::Real, n::Real)
   i_value, _ = quadgk(x -> α(x, n), 0, θ)
   return i_value
end

"""
    potential!(ψ::Real, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, n::T) where {S<:Real, T<:Real}
"""
function potential!(ψ::Real, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, n::T) where {S<:Real, T<:Real}
   κs = 4.0 * k_s

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)

   ψ_up = ψ + κs * θs^2 * ϕ(dr, n)
   return ψ_up
end

"""
    potential!(ψ::ROA, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, n::T) where {S<:ROA, T<:Real}
"""
function potential!(ψ::ROA, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, n::T) where {S<:ROA, T<:Real}
   κs = 4.0 * k_s

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = sqrt(dx^2 + dy^2)
         ψ[i, j] = ψ[i, j] + κs * θs^2 * ϕ(dr, n)
      end
   end
end


"""
    deflection!(ψx::Real, ψy::Real, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, n::T) where {S<:Real, T<:Real}
"""
function deflection!(ψx::Real, ψy::Real, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, n::T) where {S<:Real, T<:Real}
   κs = 4.0 * k_s

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)

   ψx_up = ψx + κs * θs * α(dr, n) * dx / dr
   ψy_up = ψy + κs * θs * α(dr, n) * dy / dr
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::ROA, ψy::ROA, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, n::T) where {S<:ROA, T<:Real}
"""
function deflection!(ψx::ROA, ψy::ROA, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, n::T) where {S<:ROA, T<:Real}
   κs = 4.0 * k_s

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = sqrt(dx^2 + dy^2)

         ψx[i, j] = ψx[i, j] + κs * θs * α(dr, n) * dx / dr
         ψy[i, j] = ψy[i, j] + κs * θs * α(dr, n) * dy / dr
      end
   end
end


"""
    jacobian!(ψxx::Real, ψyy::Real, ψxy::Real, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, n::T) where {S<:Real, T<:Real}
"""
function jacobian!(ψxx::Real, ψyy::Real, ψxy::Real, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, n::T) where {S<:Real, T<:Real}
   κs = 4.0 * k_s

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)

   α_r = κs * α(dr, n)
   κ_r = κs * κ(dr, n)

   ψxx_up = ψxx + 2.0 * κ_r * dx^2 / dr^2 - α_r * (dx^2 - dy^2) / dr^3
   ψyy_up = ψyy + 2.0 * κ_r * dy^2 / dr^2 + α_r * (dx^2 - dy^2) / dr^3
   ψxy_up = ψxy + 2.0 * (κ_r - α_r / dr) * dx * dy / dr^2
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::ROA, ψyy::ROA, ψxy::ROA, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, n::T) where {S<:ROA, T<:Real}
"""
function jacobian!(ψxx::ROA, ψyy::ROA, ψxy::ROA, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T, n::T) where {S<:ROA, T<:Real}
   κs = 4.0 * k_s

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = sqrt(dx^2 + dy^2)

         α_r = κs * α(dr, n)
         κ_r = κs * κ(dr, n)

         ψxx[i, j] = ψxx[i, j] + 2.0 * κ_r * dx^2 / dr^2 - α_r * (dx^2 - dy^2) / dr^3
         ψyy[i, j] = ψyy[i, j] + 2.0 * κ_r * dy^2 / dr^2 + α_r * (dx^2 - dy^2) / dr^3
         ψxy[i, j] = ψxy[i, j] + 2.0 * (κ_r - α_r / dr) * dx * dy / dr^2
      end
   end
end

end