module tNFWLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!

@inline function F_x(x::RV)
   if x < 1.0 
      arg = sqrt(1.0 - x^2)
      return atanh(arg) / arg
   else 
      arg = sqrt(x^2 - 1.0)
      return atan(arg) / arg
   end
end

@inline function L_x(x::RV, τ::RV)
   return log((sqrt(τ^2 + x^2) - τ) / x)
end

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, θt::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, θt::RV) where T <: RV
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)
   τ = θt / θs

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr = sqrt(dx^2 + dy^2)

   term1 = τ^3 * pi * ((3.0 * τ^2 - 1.0) * log(τ + sqrt(τ^2 + dr^2)) - 4.0 * τ * sqrt(τ^2 + dr^2))
   term2 = (3.0 * τ^4 - 6.0 * τ^2 - 1.0) * τ * sqrt(τ^2 + dr^2) * L_x(dr, τ) + τ^4 * (τ^2 - 3.0) * L_x(dr, τ)^2
   term3 = 8.0 * τ^4 * (dr^2 - 1.0) * F_x(dr) + τ^4 * (τ^2 - 3.0) * (dr^2 - 1.0) * F_x(dr)^2
   term4 = τ^2 * (2.0 * τ^2 * (τ^2 - 3.0) * log(τ) - 3.0 * τ^4 - 2.0 * τ^2 + 1.0) * log(dr)
   term5 = τ^2 * (τ^2 * (4.0 * τ * pi + (τ^2 - 3.0) * log(2)^2 + 8.0 * log(2)) 
         - log(2.0 * τ) * (1.0 + 6.0 * τ^2 - 3.0 * τ^4 + τ^2 * (τ^2 - 3.0) * log(2.0 * τ) + τ * pi * (3.0 * τ^2 - 1.0)))
   
   ψ_up = ψ + κs * θs^2 * (term1 + term2 + term3 + term4 + term5) * 0.5 / (τ^2 + 1.0)^3
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, θt::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, θt::RV) where T <: ROA
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)
   τ = θt / θs

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr = sqrt(dx^2 + dy^2)
         
         term1 = τ^3 * pi * ((3.0 * τ^2 - 1.0) * log(τ + sqrt(τ^2 + dr^2)) - 4.0 * τ * sqrt(τ^2 + dr^2))
         term2 = (3.0 * τ^4 - 6.0 * τ^2 - 1.0) * τ * sqrt(τ^2 + dr^2) * L_x(dr, τ) + τ^4 * (τ^2 - 3.0) * L_x(dr, τ)^2
         term3 = 8.0 * τ^4 * (dr^2 - 1.0) * F_x(dr) + τ^4 * (τ^2 - 3.0) * (dr^2 - 1.0) * F_x(dr)^2
         term4 = τ^2 * (2.0 * τ^2 * (τ^2 - 3.0) * log(τ) - 3.0 * τ^4 - 2.0 * τ^2 + 1.0) * log(dr)
         term5 = τ^2 * (τ^2 * (4.0 * τ * pi + (τ^2 - 3.0) * log(2)^2 + 8.0 * log(2)) 
               - log(2.0 * τ) * (1.0 + 6.0 * τ^2 - 3.0 * τ^4 + τ^2 * (τ^2 - 3.0) * log(2.0 * τ) + τ * pi * (3.0 * τ^2 - 1.0)))

         ψ[i, j] = ψ[i, j] + κs * θs^2 * (term1 + term2 + term3 + term4 + term5) * 0.5 / (τ^2 + 1.0)^3
      end
   end
end


function deflection!(ψx::T, ψT::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, θt::RV) where T <: RV
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)
   τ = θt / θs
end

function deflection!(ψx::T, ψT::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, ρs:: RV, θs::RV, θt::RV) where T <: ROA
   κs = 4.0 * ρs * D_d * θs * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)
   τ = θt / θs
end

end