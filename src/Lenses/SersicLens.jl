module SersicLens

# Julia packages to import
using SpecialFunctions
using HypergeometricFunctions

# LensFactory modules to use
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!
export scale_to_halflight


@inline function b_n(n::RV)
   if 0.06 < n < 0.36
      return 0.01945 - n * (0.8902 - n * (10.95 - n * (19.67 - 13.47 * n)))
   else
      return 2.0*n - (1.0/3.0) + (4.0/405.0/n) + (46.0/25515/n^2) + (131.0/1148175.0/n^3) - (2194697.0/30690717750.0/n^4)
   end
end


"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θe::RV, n::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θe::RV, n::RV) where T <: RV
   bn = b_n(n)
   θs = θe / bn^n
   κs = (4.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * gamma(2*n+1) * ANGLE_ARCSEC^2)

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)

   ψ_up = ψ + 0.5 * κs * θs^2 * dr^2 * pFq( (2*n, 2*n), (2*n+1, 2*n+1), -dr^(1/n))
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θe::RV, n::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θe::RV, n::RV) where T <: ROA
   bn = b_n(n)
   θs = θe / bn^n
   κs = (4.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * gamma(2*n+1) * ANGLE_ARCSEC^2)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = sqrt(dx^2 + dy^2)
         ψ[i, j] = ψ[i, j] + 0.5 * κs * θs^2 * dr^2 * pFq( (2*n, 2*n), (2*n+1, 2*n+1), -dr^(1/n))
      end
   end
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θe::RV, n::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θe::RV, n::RV) where T <: RV
   bn = b_n(n)
   θs = θe / bn^n
   κs = (4.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * gamma(2*n+1) * ANGLE_ARCSEC^2)

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)
   
   P, _ = gamma_inc(2*n, dr^(1/n))
   ψx_up = ψx + κs * θs * gamma(2*n+1) * P * dx / dr^2
   ψy_up = ψy + κs * θs * gamma(2*n+1) * P * dy / dr^2
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θe::RV, n::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θe::RV, n::RV) where T <: ROA
   bn = b_n(n)
   θs = θe / bn^n
   κs = (4.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * gamma(2*n+1) * ANGLE_ARCSEC^2)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = sqrt(dx^2 + dy^2)

         P, _ = gamma_inc(2*n, dr^(1/n))
         ψx[i, j] = ψx[i, j] + κs * θs * gamma(2*n+1) * P * dx / dr^2
         ψy[i, j] = ψy[i, j] + κs * θs * gamma(2*n+1) * P * dy / dr^2
      end
   end
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θe::RV, n::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θe::RV, n::RV) where T <: RV
   bn = b_n(n)
   θs = θe / bn^n
   κs = (4.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * gamma(2*n+1) * ANGLE_ARCSEC^2)

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)

   P, _ = gamma_inc(2*n, dr^(1/n))
   α_r = κs * gamma(2*n+1) * P / dr
   κ_r = κs * exp(-dr^(1/n))

   ψxx_up = ψxx + 2.0 * κ_r * dx^2 / dr^2 - α_r * (dx^2 - dy^2) / dr^3
   ψyy_up = ψyy + 2.0 * κ_r * dy^2 / dr^2 + α_r * (dx^2 - dy^2) / dr^3
   ψxy_up = ψxy + 2.0 * (κ_r - α_r / dr) * dx * dy / dr^2
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θe::RV, n::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θe::RV, n::RV) where T <: ROA
   bn = b_n(n)
   θs = θe / bn^n
   κs = (4.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * gamma(2*n+1) * ANGLE_ARCSEC^2)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = sqrt(dx^2 + dy^2)

         P, _ = gamma_inc(2*n, dr^(1/n))
         α_r = κs * gamma(2*n+1) * P / dr
         κ_r = κs * exp(-dr^(1.0/n))

         ψxx[i, j] = ψxx[i, j] + 2.0 * κ_r * dx^2 / dr^2 - α_r * (dx^2 - dy^2) / dr^3
         ψyy[i, j] = ψyy[i, j] + 2.0 * κ_r * dy^2 / dr^2 + α_r * (dx^2 - dy^2) / dr^3
         ψxy[i, j] = ψxy[i, j] + 2.0 * (κ_r - α_r / dr) * dx * dy / dr^2
      end
   end
end


"""
    scale_to_halflight(;n::RV=NaN, x_s::RV=NaN)
"""
function scale_to_halflight(;n::RV=NaN, x_s::RV=NaN)
   return x_s * b_n(n)^n
end

end