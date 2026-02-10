module tNFWLens

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

@inline function L_x(x::RV, τ::RV)
   return log(x / (sqrt(τ^2 + x^2) + τ))
end


"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, k_s:: RV, θs::RV, θt::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, k_s:: RV, θs::RV, θt::RV) where T <: RV
   κs = 4.0 * k_s
   τ = θt / θs

   t2 = τ * τ
   t4 = t2 * t2
   term0 = 0.5 / (t2 + 1.0)^3

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)
   dr2 = dx^2 + dy^2
   
   term1 = τ^3 * pi * ((3.0 * t2 - 1.0) * log(τ + sqrt(t2 + dr2)) - 4.0 * τ * sqrt(t2 + dr2))
   term2 = (3.0 * t4 - 6.0 * t2 - 1.0) * τ * sqrt(t2 + dr2) * L_x(dr, τ) + t4 * (t2 - 3.0) * L_x(dr, τ)^2
   term3 = 8.0 * t4 * (dr2 - 1.0) * F_x(dr) + t4 * (t2 - 3.0) * (dr2 - 1.0) * F_x(dr)^2
   term4 = t2 * (2.0 * t2 * (t2 - 3.0) * log(τ) - 3.0 * t4 - 2.0 * t2 + 1.0) * log(dr)
   term5 = t2 * (t2 * (4.0 * τ * pi + (t2 - 3.0) * log(2)^2 + 8.0 * log(2)) 
         - log(2.0 * τ) * (1.0 + 6.0 * t2 - 3.0 * t4 + t2 * (t2 - 3.0) * log(2.0 * τ) + τ * pi * (3.0 * t2 - 1.0)))
   ψ_up = ψ + κs * θs^2 * term0 * (term1 + term2 + term3 + term4 + term5)
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, k_s:: RV, θs::RV, θt::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, k_s:: RV, θs::RV, θt::RV) where T <: ROA
   κs = 4.0 * k_s
   τ = θt / θs

   t2 = τ * τ
   t4 = t2 * t2
   term0 = 0.5 / (t2 + 1.0)^3

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = sqrt(dx^2 + dy^2)
         dr2 = dx^2 + dy^2
         
         term1 = τ^3 * pi * ((3.0 * t2 - 1.0) * log(τ + sqrt(t2 + dr2)) - 4.0 * τ * sqrt(t2 + dr2))
         term2 = (3.0 * t4 - 6.0 * t2 - 1.0) * τ * sqrt(t2 + dr2) * L_x(dr, τ) + t4 * (t2 - 3.0) * L_x(dr, τ)^2
         term3 = 8.0 * t4 * (dr2 - 1.0) * F_x(dr) + t4 * (t2 - 3.0) * (dr2 - 1.0) * F_x(dr)^2
         term4 = t2 * (2.0 * t2 * (t2 - 3.0) * log(τ) - 3.0 * t4 - 2.0 * t2 + 1.0) * log(dr)
         term5 = t2 * (t2 * (4.0 * τ * pi + (t2 - 3.0) * log(2)^2 + 8.0 * log(2)) 
               - log(2.0 * τ) * (1.0 + 6.0 * t2 - 3.0 * t4 + t2 * (t2 - 3.0) * log(2.0 * τ) + τ * pi * (3.0 * t2 - 1.0)))

         ψ[i, j] = ψ[i, j] + κs * θs^2 * term0 * (term1 + term2 + term3 + term4 + term5)
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, k_s:: RV, θs::RV, θt::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, k_s:: RV, θs::RV, θt::RV) where T <: RV
   κs = 4.0 * k_s
   τ = θt / θs
   
   t2 = τ * τ
   t4 = t2 * t2
   term0 = 0.5 * t4 / (t2 + 1.0)^3

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)
   dr2 = dx^2 + dy^2

   term1 = 2.0 * F_x(dr) * (t2 + 1.0 + 4.0 * (dr2 - 1.0))   
   term2 = (pi * (3.0 * t2 - 1.0) + 2.0 * τ * (t2 - 3.0) * log(τ)) / τ
   term3 = - t2 * τ * pi * (4.0 * (t2 + dr2) - t2 - 1.0) + (t2 * (1.0 - t4) + (t2 + dr2) * (3.0 * t4 - 6.0 * t2 - 1.0)) * L_x(dr, τ)
   α_r = term0 * (term1 + term2 + term3 / t2 / τ / sqrt(t2 + dr2)) / dr

   ψx_up = ψx + κs * θs * α_r * dx / dr
   ψy_up = ψy + κs * θs * α_r * dy / dr
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, k_s:: RV, θs::RV, θt::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, k_s:: RV, θs::RV, θt::RV) where T <: ROA
   κs = 4.0 * k_s
   τ = θt / θs
   
   t2 = τ * τ
   t4 = t2 * t2
   term0 = 0.5 * t4 / (t2 + 1.0)^3

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = sqrt(dx^2 + dy^2)
         dr2 = dx^2 + dy^2

         term1 = 2.0 * F_x(dr) * (t2 + 1.0 + 4.0 * (dr2 - 1.0))   
         term2 = (pi * (3.0 * t2 - 1.0) + 2.0 * τ * (t2 - 3.0) * log(τ)) / τ
         term3 = - t2 * τ * pi * (4.0 * (t2 + dr2) - t2 - 1.0) + (t2 * (1.0 - t4) + (t2 + dr2) * (3.0 * t4 - 6.0 * t2 - 1.0)) * L_x(dr, τ)
         α_r = term0 * (term1 + term2 + term3 / t2 / τ / sqrt(t2 + dr2)) / dr
         
         ψx[i, j] = ψx[i, j] + κs * θs * α_r * dx / dr
         ψy[i, j] = ψy[i, j] + κs * θs * α_r * dy / dr
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, k_s:: RV, θs::RV, θt::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, k_s:: RV, θs::RV, θt::RV) where T <: RV
   κs = 4.0 * k_s
   τ = θt / θs
   
   t2 = τ * τ
   t4 = t2 * t2
   term0_α = 0.50 * t4 / (t2 + 1.0)^3
   term0_κ = 0.25 * t4 / (t2 + 1.0)^3

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)
   dr2 = dx^2 + dy^2

   term1 = 2.0 * F_x(dr) * (t2 + 1.0 + 4.0 * (dr2 - 1.0))   
   term2 = (pi * (3.0 * t2 - 1.0) + 2.0 * τ * (t2 - 3.0) * log(τ)) / τ
   term3 = - t2 * τ * pi * (4.0 * (t2 + dr2) - t2 - 1.0) + (t2 * (1.0 - t4) + (t2 + dr2) * (3.0 * t4 - 6.0 * t2 - 1.0)) * L_x(dr, τ)
   α_r = κs * term0_α * (term1 + term2 + term3 / (t2 * τ * sqrt(t2 + dr2))) / dr

   term1 = 2.0 * (t2 + 1.0) * (1.0 - F_x(dr)) / (dr2 - 1.0) + 8.0 * F_x(dr)
   term2 = (t4 - 1.0) / t2 / (t2 + dr2) - pi * (4.0 * (t2 + dr2) + t2 + 1.0) / (t2 + dr2)^(3/2)
   term3 = (t2 * (t4 - 1.0) + (t2 + dr2) * (3.0 * t4 - 6.0 * t2 - 1.0)) * L_x(dr, τ) / (t2 * τ * (t2 + dr2)^(3/2))
   κ_r = κs * term0_κ * (term1 + term2 + term3)

   ψxx_up = ψxx + 2.0 * κ_r * dx^2 / dr^2 - α_r * (dx^2 - dy^2) / dr^3
   ψyy_up = ψyy + 2.0 * κ_r * dy^2 / dr^2 + α_r * (dx^2 - dy^2) / dr^3
   ψxy_up = ψxy + 2.0 * (κ_r - α_r / dr) * dx * dy / dr^2
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, k_s:: RV, θs::RV, θt::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, k_s:: RV, θs::RV, θt::RV) where T <: ROA
   κs = 4.0 * k_s
   τ = θt / θs
   
   t2 = τ * τ
   t4 = t2 * t2
   term0_α = 0.50 * t4 / (t2 + 1.0)^3
   term0_κ = 0.25 * t4 / (t2 + 1.0)^3

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = sqrt(dx^2 + dy^2)
         dr2 = dx^2 + dy^2

         term1 = 2.0 * F_x(dr) * (t2 + 1.0 + 4.0 * (dr2 - 1.0))   
         term2 = (pi * (3.0 * t2 - 1.0) + 2.0 * τ * (t2 - 3.0) * log(τ)) / τ
         term3 = - t2 * τ * pi * (4.0 * (t2 + dr2) - t2 - 1.0) + (t2 * (1.0 - t4) + (t2 + dr2) * (3.0 * t4 - 6.0 * t2 - 1.0)) * L_x(dr, τ)
         α_r = κs * term0_α * (term1 + term2 + term3 / (t2 * τ * sqrt(t2 + dr2))) / dr

         term1 = 2.0 * (t2 + 1.0) * (1.0 - F_x(dr)) / (dr2 - 1.0) + 8.0 * F_x(dr)
         term2 = (t4 - 1.0) / t2 / (t2 + dr2) - pi * (4.0 * (t2 + dr2) + t2 + 1.0) / (t2 + dr2)^(3/2)
         term3 = (t2 * (t4 - 1.0) + (t2 + dr2) * (3.0 * t4 - 6.0 * t2 - 1.0)) * L_x(dr, τ) / (t2 * τ * (t2 + dr2)^(3/2))
         κ_r = κs * term0_κ * (term1 + term2 + term3)

         ψxx[i, j] = ψxx[i, j] + 2.0 * κ_r * dx^2 / dr^2 - α_r * (dx^2 - dy^2) / dr^3
         ψyy[i, j] = ψyy[i, j] + 2.0 * κ_r * dy^2 / dr^2 + α_r * (dx^2 - dy^2) / dr^3
         ψxy[i, j] = ψxy[i, j] + 2.0 * (κ_r - α_r / dr) * dx * dy / dr^2
      end
   end
end

end