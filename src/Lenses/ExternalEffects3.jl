module ExternalEffects3


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
    potential!(ψ::T, θx::T, θy::T, δ::RV, ϕ::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, δ::RV, ϕ::RV) where T <: RV
   # Precompute angles
   ϕ = deg2rad(ϕ)
   cos_phi = cos(ϕ)
   sin_phi = sin(ϕ)
   
   # Get radial coordinate
   r = sqrt(θx^2 + θy^2)

   # cos(θ - ϕ) and sin(θ - ϕ)
   cos_value = (θx * cos_phi + θy * sin_phi) / r
   sin_value = (θy * cos_phi - θx * sin_phi) / r

   ψ_up = ψ + δ * r^3 * cos_value * sin_value^2
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, δ::RV, ϕ::RV) where T <: ROA
Calculate potential at given coordinates for "restricted" third order perturbations corresponding
to SIS lens model and update the potential in place.

# Arguments
- `ψ`: Potential at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `δ::RV`: Amplitude of third order perturbations.
- `ϕ::RV`: External Shear angle (in ``\\rm \\mathbf{degrees}``).
"""
function potential!(ψ::T, θx::T, θy::T, δ::RV, ϕ::RV) where T <: ROA
   # Precompute angles
   ϕ = deg2rad(ϕ)
   cos_phi = cos(ϕ)
   sin_phi = sin(ϕ)
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Get radial coordinate
         dx = θx[i, j]
         dy = θy[i, j]
         r = sqrt(dx^2 + dy^2)

         # cos(θ - ϕ) and sin(θ - ϕ)
         cos_value = (dx * cos_phi + dy * sin_phi) / r
         sin_value = (dy * cos_phi - dx * sin_phi) / r

         ψ[i, j] = ψ[i, j] + δ * r^3 * cos_value * sin_value^2
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, δ::RV, ϕ::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, δ::RV, ϕ::RV) where T <: RV
   # Precompute angles
   ϕ = deg2rad(ϕ)
   cos_phi = cos(ϕ)
   sin_phi = sin(ϕ)
   
   # Get radial coordinate
   r = sqrt(θx^2 + θy^2)

   # cos(θ - ϕ) and sin(θ - ϕ)
   cos_value = (θx * cos_phi + θy * sin_phi) / r
   sin_value = (θy * cos_phi - θx * sin_phi) / r

   # Deflection
   term1 = 3.0 * cos_value * sin_value^2
   term2 = sin_value * (2.0 * cos_value^2 - sin_value^2)
   ψx_up = ψx + δ * r * (term1 * θx - term2 * θy)
   ψy_up = ψy + δ * r * (term1 * θy + term2 * θx)
   return ψx_up, ψy_up
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, δ::RV, ϕ::RV) where T <: ROA
Calculate deflection at given coordinates for "restricted" third order perturbations corresponding
to SIS lens model and update the deflection in place.

# Arguments
- `ψx`: x-component of deflection at given coordinates
- `ψy`: y-component of deflection at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `δ::RV`: Amplitude of third order perturbations.
- `ϕ::RV`: External Shear angle (in ``\\rm \\mathbf{degrees}``).
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, δ::RV, ϕ::RV) where T <: ROA
   # Precompute angles
   ϕ = deg2rad(ϕ)
   cos_phi = cos(ϕ)
   sin_phi = sin(ϕ)
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Get radial coordinate
         dx = θx[i, j]
         dy = θy[i, j]
         r = sqrt(dx^2 + dy^2)

         # cos(θ - ϕ) and sin(θ - ϕ)
         cos_value = (dx * cos_phi + dy * sin_phi) / r
         sin_value = (dy * cos_phi - dx * sin_phi) / r

         # Deflection
         term1 = 3.0 * cos_value * sin_value^2
         term2 = sin_value * (2.0 * cos_value^2 - sin_value^2)
         ψx[i, j] = ψx[i, j] + δ * r * (term1 * dx - term2 * dy)
         ψy[i, j] = ψy[i, j] + δ * r * (term1 * dy + term2 * dx)
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, δ::RV, ϕ::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, δ::RV, ϕ::RV) where T <: RV

end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, δ::RV, ϕ::RV) where T <: ROA
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, δ::RV, ϕ::RV) where T <: ROA
   
end

end