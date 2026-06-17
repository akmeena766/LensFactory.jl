module gNFWLens

# Inbuilt Julia functions to import
using QuadGK

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!

@inline function I_κ(z::Real, θ::Real, n::Real)
   return 1.0 / (θ^2 + z^2)^(0.5 * n) / (1.0 + (θ^2 + z^2)^0.5)^(3.0 - n)
end

function κ(θ::Real, n::Real)
   i_value, _  = quadgk(x -> I_κ(x, θ, n), 0, Inf)
   return 0.5 * i_value
end

@inline function I_α(z::Real, n::Real)
   return z * κ(z, n)
end

function α(θ::Real, n::Real)
   i_value, _ = quadgk(x -> I_α(x, n), 0, θ)
   return 2.0 * i_value / θ
end

function ϕ(θ::Real, n::Real)
   i_value, _ = quadgk(x -> α(x, n), 0, θ)
   return i_value
end


"""
    potential!(ψ::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, k_s:: Real, θs::Real, n::Real) where T <: Real
"""
function potential!(ψ::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, k_s:: Real, θs::Real, n::Real) where T <: Real
   κs = 4.0 * k_s

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)

   ψ_up = ψ + κs * θs^2 * ϕ(dr, n)
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, k_s:: Real, θs::Real, n::Real) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, k_s:: Real, θs::Real, n::Real) where T <: ROA
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
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, k_s:: Real, θs::Real, n::Real) where T <: Real
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, k_s:: Real, θs::Real, n::Real) where T <: Real
   κs = 4.0 * k_s

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)

   ψx_up = ψx + κs * θs * α(dr, n) * dx / dr
   ψy_up = ψy + κs * θs * α(dr, n) * dy / dr
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, k_s:: Real, θs::Real, n::Real) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, k_s:: Real, θs::Real, n::Real) where T <: ROA
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
    deflection!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, k_s:: Real, θs::Real, n::Real) where T <: Real
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, k_s:: Real, θs::Real, n::Real) where T <: Real
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
    deflection!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, k_s:: Real, θs::Real, n::Real) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::Real, θxc::Real, θyc::Real, k_s:: Real, θs::Real, n::Real) where T <: ROA
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