module PlummerLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!
export einstein_angle


"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   θE2 = (2.0 * CONST_G * mass / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
   θs2 = θs^2

   dx = θx - θxc
   dy = θy - θyc
   
   ψ_up = ψ + θE2 * log(θs2 + dx^2 + dy^2)
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   θE2 = (2.0 * CONST_G * mass / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
   θs2 = θs^2

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds @simd for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         dr2 = θs2 + dx^2 + dy^2

         ψ[i, j] = ψ[i, j] + θE2 * log(dr2)
      end
   end
   return nothing
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   θE2 = (4.0 * CONST_G * mass / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
   θs2 = θs^2

   dx = θx - θxc
   dy = θy - θyc
   dr2 = θs2 + dx^2 + dy^2

   ψx_up = ψx + θE2 * dx / dr2
   ψy_up = ψy + θE2 * dy / dr2
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   θE2 = (4.0 * CONST_G * mass / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
   θs2 = θs^2

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds @simd for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         dr2 = θs2 + dx^2 + dy^2
         
         ψx[i, j] = ψx[i, j] + θE2 * dx / dr2
         ψy[i, j] = ψy[i, j] + θE2 * dy / dr2
      end
   end
   return nothing
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   θE2 = (4.0 * CONST_G * mass / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
   θs2 = θs^2

   dx = θx - θxc
   dy = θy - θyc
   dr4 = (θs2 + dx^2 + dy^2)^2

   ψxx_up = ψxx + θE2 * (θs2 - dx^2 + dy^2) / dr4
   ψyy_up = ψyy + θE2 * (θs2 + dx^2 - dy^2) / dr4
   ψxy_up = ψxy - θE2 * 2.0 * dx * dy / dr4
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   θE2 = (4.0 * CONST_G * mass / CONST_C^2 / D_d) / ANGLE_ARCSEC^2
   θs2 = θs^2

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds @simd for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         dr4 = (θs2 + dx^2 + dy^2)^2
         ψxx[i, j] = ψxx[i, j] + θE2 * (θs2 - dx^2 + dy^2) / dr4
         ψyy[i, j] = ψyy[i, j] + θE2 * (θs2 + dx^2 - dy^2) / dr4
         ψxy[i, j] = ψxy[i, j] - θE2 * 2.0 * dx * dy / dr4
      end
   end
   return nothing
end


"""
    einstein_angle(;D_d::Float64=NaN, D_ds::Float64=NaN, D_s::Float64=NaN, mass::Float64=NaN, x_s::Float64=NaN)
"""
function einstein_angle(;D_d::Float64=NaN, D_ds::Float64=NaN, D_s::Float64=NaN, mass::Float64=NaN, x_s::Float64=NaN)
   return sqrt((4.0 * CONST_G * mass / CONST_C^2) * (D_ds / D_d / D_s) / ANGLE_ARCSEC^2 - x_s^2)
end

end