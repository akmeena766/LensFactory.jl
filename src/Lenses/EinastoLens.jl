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


function m_Ein(θ::RV, n::RV)
   Pax, _ = gamma_inc(3.0 / n, (2.0 / n) * θ^n)
   return (1.0 / n) * (n / 2.0)^(3.0 / n) * gamma(3.0 / n) * Pax
end

@inline function I_κ(z::RV, θ::RV, n::RV)
   return exp(-(2.0 / n) * (θ^2 + z^2)^(0.5 * n))
end

function κ(θ::RV, n::RV)
   i_value, _  = quadgk(x -> I_κ(x, θ, n), 0, Inf)
   return 0.5 * i_value
end

@inline function I_α(z::RV, θ::RV, n::RV)
   return (m_Ein(θ * sqrt(1.0 + z^2), n) + θ * m_Ein(θ * sqrt(1.0 + z^2) / z, n)) / (1.0 + z)^1.5
end

function α(θ::RV, n::RV)
   i_value, _ = quadgk(x -> I_α(x, θ, n), 0, Inf)
   return i_value / θ
end

function ϕ(θ::RV, n::RV)
   i_value, _ = quadgk(x -> α(θ, n), 0, θ)
   return i_value
end

function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, n::RV) where T <: RV
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)

   ψ_up = ψ + κs * θs^2 * ϕ(dr, n)
   return ψ_up
end

function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, n::RV) where T <: ROA
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)

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

function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, n::RV) where T <: RV
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)

   ψx_up = ψx + κs * θs * α(dr, n) * dx / dr
   ψy_up = ψy + κs * θs * α(dr, n) * dy / dr
   return ψx_up, ψy_up
end

function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, n::RV) where T <: ROA
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)

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


function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, n::RV) where T <: RV
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)

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

function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, n::RV) where T <: ROA
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)

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