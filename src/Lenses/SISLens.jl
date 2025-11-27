module SISLens

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!
export einstein_angle


"""
    potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: RV
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   
   dx = θx - θxc
   dy = θy - θyc
   
   ψ_up = ψ + θE * sqrt(dx^2 + dy^2)
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: ROA
"""
function potential!(ψ::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: ROA
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc

         ψ[i, j] = ψ[i, j] + θE * sqrt(dx^2 + dy^2)
      end
   end
   return nothing
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: RV
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC

   dx = θx - θxc
   dy = θy - θyc
   dr = sqrt(dx^2 + dy^2)
   
   ψx_up = ψx + θE * dx / dr
   ψy_up = ψy + θE * dy / dr
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: ROA
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: ROA
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         dr = sqrt(dx^2 + dy^2)

         ψx[i, j] = ψx[i, j] + θE * dx / dr
         ψy[i, j] = ψy[i, j] + θE * dy / dr
      end
   end
   return nothing
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: RV
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   
   dx = θx - θxc
   dy = θy - θyc
   dr2 = dx^2 + dy^2
   dr3 = dr2 * sqrt(dr2)

   ψxx_up = ψxx + θE * dy^2 / dr3
   ψyy_up = ψyy + θE * dx^2 / dr3
   ψxy_up = ψxy - θE * dx * dy / dr3
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, θxc::RV, θyc::RV, vd::RV) where T <: ROA
   θE = 4.0 * pi * (vd / CONST_C)^2 / ANGLE_ARCSEC
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = θx[i, j] - θxc
         dy = θy[i, j] - θyc
         dr2 = dx^2 + dy^2
         dr3 = dr2 * sqrt(dr2)
         
         ψxx[i, j] = ψxx[i, j] + θE * dy^2 / dr3
         ψyy[i, j] = ψyy[i, j] + θE * dx^2 / dr3
         ψxy[i, j] = ψxy[i, j] - θE * dx * dy / dr3
      end
   end
   return nothing
end


"""
    einstein_angle(;D_ds::Float64=NaN, D_s::Float64=NaN, v_d::RV=NaN)
"""
function einstein_angle(;D_ds::Float64=NaN, D_s::Float64=NaN, v_d::RV=NaN)
   return 4π * (v_d / CONST_C)^2 * (D_ds / D_s) / ANGLE_ARCSEC
end

end