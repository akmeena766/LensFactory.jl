module NFWLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!

@inline function F_x(x::Real)
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

"""
    potential!(ψ::U, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T) where {U<:Real, S<:Real, T<:Real}
"""
function potential!(ψ::U, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T) where {U<:Real, S<:Real, T<:Real}
   κs = 4.0 * k_s
   κs = κs * θs^2 * 0.5

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)

   ψ_up = ψ + κs * ((dr^2 - 1.0) * F_x(dr)^2 + log(0.5 * dr)^2)
   return ψ_up
end

"""
    potential!(ψ::U, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T) where {U<:ROA, S<:ROA, T<:Real}
"""
function potential!(ψ::U, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T) where {U<:ROA, S<:ROA, T<:Real}
   κs = 4.0 * k_s
   κs = κs * θs^2 * 0.5

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = sqrt(dx^2 + dy^2)
         ψ[i, j] = ψ[i, j] + κs * ((dr^2 - 1.0) * F_x(dr)^2 + log(0.5 * dr)^2)
      end
   end
end


"""
    deflection!(ψx::U, ψy::U, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T) where {U<:Real, S<:Real, T<:Real}
"""
function deflection!(ψx::U, ψy::U, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T) where {U<:Real, S<:Real, T<:Real}
   κs = 4.0 * k_s
   κs = κs * θs

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)

   α_r = (F_x(dr) + log(0.5 * dr)) / dr

   ψx_up = ψx + κs * α_r * dx / dr
   ψy_up = ψy + κs * α_r * dy / dr
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::U, ψy::U, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T) where {U<:ROA, S<:ROA, T<:Real}
"""
function deflection!(ψx::U, ψy::U, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T) where {U<:ROA, S<:ROA, T<:Real}
   κs = 4.0 * k_s
   κs = κs * θs

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = sqrt(dx^2 + dy^2)
         
         α_r = (F_x(dr) + log(0.5 * dr)) / dr

         ψx[i, j] = ψx[i, j] + κs * α_r * dx / dr
         ψy[i, j] = ψy[i, j] + κs * α_r * dy / dr
      end
   end
end


"""
    jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s:: T, θs::T) where {U<:Real, S<:Real, T<:Real}
"""
function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s:: T, θs::T) where {U<:Real, S<:Real, T<:Real}
   κs = 4.0 * k_s
   
   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)

   α_r = κs * (F_x(dr) + log(0.5 * dr)) / dr
   κ_r = κs * 0.5 * (1.0 - F_x(dr)) / (dr^2 - 1.0)

   ψxx_up = ψxx + 2.0 * κ_r * dx^2 / dr^2 - α_r * (dx^2 - dy^2) / dr^3
   ψyy_up = ψyy + 2.0 * κ_r * dy^2 / dr^2 + α_r * (dx^2 - dy^2) / dr^3
   ψxy_up = ψxy + 2.0 * (κ_r - α_r / dr) * dx * dy / dr^2
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T) where {U<:ROA, S<:ROA, T<:Real}
"""
function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, D_d::T, θxc::T, θyc::T, k_s::T, θs::T) where {U<:ROA, S<:ROA, T<:Real}
   κs = 4.0 * k_s

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = sqrt(dx^2 + dy^2)

         α_r = κs * (F_x(dr) + log(0.5 * dr)) / dr
         κ_r = κs * 0.5 * (1.0 - F_x(dr)) / (dr^2 - 1.0)

         ψxx[i, j] = ψxx[i, j] + 2.0 * κ_r * dx^2 / dr^2 - α_r * (dx^2 - dy^2) / dr^3
         ψyy[i, j] = ψyy[i, j] + 2.0 * κ_r * dy^2 / dr^2 + α_r * (dx^2 - dy^2) / dr^3
         ψxy[i, j] = ψxy[i, j] + 2.0 * (κ_r - α_r / dr) * dx * dy / dr^2
      end
   end
end

end