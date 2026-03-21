module HernquistLens

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

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)

   dx = (θx - θxc) / θs + eps()
   dy = (θy - θyc) / θs + eps()
   dr = sqrt(dx^2 + dy^2)

   ψ_up = ψ + κs * θs^2 * (2.0 * F_x(dr) + log(0.25 * dr^2))
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
Calculate potential at given coordinates for Hernquist lens and update the potential values in-place.

# Arguments
- `ψ`: Potential at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `D_d::RV`: Angular diameter distance to the lens (in ``\\rm \\mathbf{Mpc}``).
- `θxc::RV`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc::RV`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::RV`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `θs::RV`: Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs + eps()
         dy = (θy[i, j] - θyc) / θs + eps()
         dr = sqrt(dx^2 + dy^2)
         ψ[i, j] = ψ[i, j] + κs * θs^2 * (2.0 * F_x(dr) + log(0.25 * dr^2))
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)

   dx = (θx - θxc) / θs + eps()
   dy = (θy - θyc) / θs + eps()
   dr = sqrt(dx^2 + dy^2)

   α_r = 2.0 * dr * (1.0 - F_x(dr)) / (dr^2 - 1.0)

   ψx_up = ψx + κs * θs * α_r * dx / dr
   ψy_up = ψy + κs * θs * α_r * dy / dr
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
Calculate deflection at given coordinates for Hernquist lens and update the deflection values in-place.

# Arguments
- `ψx`: x-component of deflection at given coordinates
- `ψy`: y-component of deflection at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `D_d::RV`: Angular diameter distance to the lens (in ``\\rm \\mathbf{Mpc}``).
- `θxc::RV`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc::RV`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::RV`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `θs::RV`: Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs + eps()
         dy = (θy[i, j] - θyc) / θs + eps()
         dr = sqrt(dx^2 + dy^2)

         α_r = 2.0 * dr * (1.0 - F_x(dr)) / (dr^2 - 1.0)
         
         ψx[i, j] = ψx[i, j] + κs * θs * α_r * dx / dr
         ψy[i, j] = ψy[i, j] + κs * θs * α_r * dy / dr
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   
   dx = (θx - θxc) / θs + eps()
   dy = (θy - θyc) / θs + eps()
   dr = sqrt(dx^2 + dy^2)

   α_r = κs * 2.0 * dr * (1.0 - F_x(dr)) / (dr^2 - 1.0)
   κ_r = κs * (-3.0 + (2.0 + dr^2) * F_x(dr)) / (dr^2 - 1.0)^2

   ψxx_up = ψxx + 2.0 * κ_r * dx^2 / dr^2 - α_r * (dx^2 - dy^2) / dr^3
   ψyy_up = ψyy + 2.0 * κ_r * dy^2 / dr^2 + α_r * (dx^2 - dy^2) / dr^3
   ψxy_up = ψxy + 2.0 * (κ_r - α_r / dr) * dx * dy / dr^2

   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
Calculate Jacobian at given coordinates for Hernquist lens and update the Jacobian values in-place.

# Arguments
- `ψxx`: xx-component of Jacobian at given coordinates
- `ψyy`: yy-component of Jacobian at given coordinates
- `ψxy`: xy-component of Jacobian at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `D_d::RV`: Angular diameter distance to the lens (in ``\\rm \\mathbf{Mpc}``).
- `θxc::RV`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc::RV`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::RV`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `θs::RV`: Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κs = (2.0 * CONST_G * mass / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs + eps()
         dy = (θy[i, j] - θyc) / θs + eps()
         dr = sqrt(dx^2 + dy^2)

         α_r = κs * 2.0 * dr * (1.0 - F_x(dr)) / (dr^2 - 1.0)
         κ_r = κs * (-3.0 + (2.0 + dr^2) * F_x(dr)) / (dr^2 - 1.0)^2

         ψxx[i, j] = ψxx[i, j] + 2.0 * κ_r * dx^2 / dr^2 - α_r * (dx^2 - dy^2) / dr^3
         ψyy[i, j] = ψyy[i, j] + 2.0 * κ_r * dy^2 / dr^2 + α_r * (dx^2 - dy^2) / dr^3
         ψxy[i, j] = ψxy[i, j] + 2.0 * (κ_r - α_r / dr) * dx * dy / dr^2
      end
   end
end

end