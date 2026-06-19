module ExternalEffects


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..Constants


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export potential! 
export deflection!
export jacobian!

"""
    potential!(ψ::Real, θx::Real, θy::Real, kappa::T, gamma::T, angle::T) where T <: Real
"""
function potential!(ψ::Real, θx::Real, θy::Real, kappa::T, gamma::T, angle::T) where T <: Real
   angle = deg2rad(angle)
   gamma1 = gamma * cos(2.0 * angle)
   gamma2 = gamma * sin(2.0 * angle)

   f1 = 0.5 * (kappa + gamma1)
   f2 = 0.5 * (kappa - gamma1)

   ψ_up = ψ + f1 * θx^2 + f2 * θy^2 + gamma2 * θx * θy
   return ψ_up
end

"""
    potential!(ψ::ROA, θx::ROA, θy::ROA, kappa::T, gamma::T, angle::T) where T <: Real
Calculate potential at given coordinates for constant external convergence and shear and update the 
potential in place.

# Arguments
- `ψ` : Potential at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `κ` : External convergence.
- `γ` : External shear value.
- `ϕ` : External Shear angle (in ``\\rm \\mathbf{degrees}``).
"""
function potential!(ψ::ROA, θx::ROA, θy::ROA, kappa::T, gamma::T, angle::T) where T <: Real
   angle = deg2rad(angle)
   gamma1 = gamma * cos(2.0 * angle)
   gamma2 = gamma * sin(2.0 * angle)

   f1 = 0.5 * (kappa + gamma1)
   f2 = 0.5 * (kappa - gamma1)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψ[i, j] = ψ[i, j] + f1 * θx[i, j]^2 + f2 * θy[i, j]^2 + gamma2 * θx[i, j] * θy[i, j]
      end
   end
end


"""
    deflection!(ψx::Real, ψy::Real, θx::Real, θy::Real, kappa::T, gamma::T, angle::T) where T <: Real
"""
function deflection!(ψx::Real, ψy::Real, θx::Real, θy::Real, kappa::T, gamma::T, angle::T) where T <: Real
   angle = deg2rad(angle)
   gamma1 = gamma * cos(2.0 * angle)
   gamma2 = gamma * sin(2.0 * angle)

   f1 = (kappa + gamma1)
   f2 = (kappa - gamma1)

   ψx_up = ψx + f1 * θx + gamma2 * θy
   ψy_up = ψy + f2 * θy + gamma2 * θx
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::ROA, ψy::ROA, θx::ROA, θy::ROA, kappa::T, gamma::T, angle::T) where T <: Real
Calculate deflection at given coordinates for constant external convergence and shear and update 
the deflection in place.

# Arguments
- `ψx`: x-component of deflection at given coordinates
- `ψy`: y-component of deflection at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `κ` : External convergence.
- `γ` : External shear value.
- `ϕ` : External Shear angle (in ``\\rm \\mathbf{degrees}``).
"""
function deflection!(ψx::ROA, ψy::ROA, θx::ROA, θy::ROA, kappa::T, gamma::T, angle::T) where T <: Real
   angle = deg2rad(angle)
   gamma1 = gamma * cos(2.0 * angle)
   gamma2 = gamma * sin(2.0 * angle)

   f1 = (kappa + gamma1)
   f2 = (kappa - gamma1)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
      ψx[i, j] = ψx[i, j] + f1 * θx[i, j] + gamma2 * θy[i, j]
      ψy[i, j] = ψy[i, j] + f2 * θy[i, j] + gamma2 * θx[i, j]
      end
   end
end


"""
    jacobian!(ψxx::Real, ψyy::Real, ψxy::Real, θx::Real, θy::Real, kappa::T, gamma::T, angle::T) where T <: Real
"""
function jacobian!(ψxx::Real, ψyy::Real, ψxy::Real, θx::Real, θy::Real, kappa::T, gamma::T, angle::T) where T <: Real
   angle = deg2rad(angle)
   gamma1 = gamma * cos(2.0 * angle)
   gamma2 = gamma * sin(2.0 * angle)

   f1 = (kappa + gamma1)
   f2 = (kappa - gamma1)

   ψxx_up = ψxx + f1
   ψyy_up = ψyy + f2
   ψxy_up = ψxy + gamma2
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::ROA, ψyy::ROA, ψxy::ROA, θx::ROA, θy::ROA, kappa::T, gamma::T, angle::T) where T <: Real
Calculate Jacobian at given coordinates for constant external convergence and shear and update the 
Jacobian in place.

# Arguments
- `ψxx`: xx-component of Jacobian at given coordinates
- `ψyy`: yy-component of Jacobian at given coordinates
- `ψxy`: xy-component of Jacobian at given coordinates
- `θx` : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy` : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `κ`  : External convergence.
- `γ`  : External shear value.
- `ϕ`  : External Shear angle (in ``\\rm \\mathbf{degrees}``).
"""
function jacobian!(ψxx::ROA, ψyy::ROA, ψxy::ROA, θx::ROA, θy::ROA, kappa::T, gamma::T, angle::T) where T <: Real
   angle = deg2rad(angle)
   gamma1 = gamma * cos(2.0 * angle)
   gamma2 = gamma * sin(2.0 * angle)
   
   f1 = (kappa + gamma1)
   f2 = (kappa - gamma1)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ψxx[i, j] = ψxx[i, j] + f1
         ψyy[i, j] = ψyy[i, j] + f2
         ψxy[i, j] = ψxy[i, j] + gamma2
      end
   end
end

end