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
    potential!(ψ::Real, θx::Real, θy::Real, δ::T, ϕ::T) where T <: Real
"""
function potential!(ψ::Real, θx::Real, θy::Real, δ::T, ϕ::T) where T <: Real
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
    potential!(ψ::ROA, θx::ROA, θy::ROA, δ::T, ϕ::T) where T <: Real
Calculate potential at given coordinates for "restricted" third order perturbations corresponding
to SIS lens model and update the potential in place.

# Arguments
- `ψ` : Potential at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `δ` : Amplitude of third order perturbations.
- `ϕ` : External Shear angle (in ``\\rm \\mathbf{degrees}``).
"""
function potential!(ψ::ROA, θx::ROA, θy::ROA, δ::T, ϕ::T) where T <: Real
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
    deflection!(ψx::Real, ψy::Real, θx::Real, θy::Real, δ::T, ϕ::T) where T <: Real
"""
function deflection!(ψx::Real, ψy::Real, θx::Real, θy::Real, δ::T, ϕ::T) where T <: Real
   # Precompute angles
   ϕ = deg2rad(ϕ)
   cos_phi = cos(ϕ)
   sin_phi = sin(ϕ)
   
   # Get radial coordinate
   r = sqrt(θx^2 + θy^2)

   # cos(θ - ϕ) and sin(θ - ϕ)
   cos_value = (θx * cos_phi + θy * sin_phi) / r
   sin_value = (θy * cos_phi - θx * sin_phi) / r

   # cos3(θ - ϕ) and sin3(θ - ϕ)
   cos_3value = 4.0 * cos_value^3 - 3.0 * cos_value
   sin_3value = 3.0 * sin_value - 4.0 * sin_value^3

   # Deflection
   term1 = 3.0 * (cos_value - cos_3value)
   term2 = (sin_value - 3.0 * sin_3value)
   ψx_up = ψx + 0.25 * δ * r * (term1 * θx + term2 * θy)
   ψy_up = ψy + 0.25 * δ * r * (term1 * θy - term2 * θx)
   return ψx_up, ψy_up
end


"""
    deflection!(ψx::ROA, ψy::ROA, θx::ROA, θy::ROA, δ::T, ϕ::T) where T <: Real
Calculate deflection at given coordinates for "restricted" third order perturbations corresponding
to SIS lens model and update the deflection in place.

# Arguments
- `ψx`: x-component of deflection at given coordinates
- `ψy`: y-component of deflection at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `δ` : Amplitude of third order perturbations.
- `ϕ` : External Shear angle (in ``\\rm \\mathbf{degrees}``).
"""
function deflection!(ψx::ROA, ψy::ROA, θx::ROA, θy::ROA, δ::T, ϕ::T) where T <: Real
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

         # cos3(θ - ϕ) and sin3(θ - ϕ)
         cos_3value = 4.0 * cos_value^3 - 3.0 * cos_value
         sin_3value = 3.0 * sin_value - 4.0 * sin_value^3

         # Deflection
         term1 = 3.0 * (cos_value - cos_3value)
         term2 = (sin_value - 3.0 * sin_3value)
         ψx[i, j] = ψx[i, j] + 0.25 * δ * r * (term1 * dx + term2 * dy)
         ψy[i, j] = ψy[i, j] + 0.25 * δ * r * (term1 * dy - term2 * dx)
      end
   end
end


"""
    jacobian!(ψxx::Real, ψyy::Real, ψxy::Real, θx::Real, θy::Real, δ::T, ϕ::T) where T <: Real
"""
function jacobian!(ψxx::Real, ψyy::Real, ψxy::Real, θx::Real, θy::Real, δ::T, ϕ::T) where T <: Real
   # Precompute angles
   ϕ = deg2rad(ϕ)
   cos_phi = cos(ϕ)
   sin_phi = sin(ϕ)
   
   # Get radial coordinate
   r = sqrt(θx^2 + θy^2)

   # cos(θ - ϕ) and sin(θ - ϕ)
   cos_value = (θx * cos_phi + θy * sin_phi) / r
   sin_value = (θy * cos_phi - θx * sin_phi) / r

   # cos3(θ - ϕ) and sin3(θ - ϕ)
   cos_3value = 4.0 * cos_value^3 - 3.0 * cos_value
   sin_3value = 3.0 * sin_value - 4.0 * sin_value^3

   # Jacobian
   term1 = 3.0 * (cos_value - cos_3value)
   term2 = (sin_value - 3.0 * sin_3value)
   ψxx_up = ψxx + 0.25 * δ * ((r^2 + θx^2) * term1 + 4.0 * θx * θy * term2 - θy^2 * (cos_value - 9.0 * cos_3value)) / r
   ψyy_up = ψyy + 0.25 * δ * ((r^2 + θy^2) * term1 - 4.0 * θx * θy * term2 - θx^2 * (cos_value - 9.0 * cos_3value)) / r
   ψxy_up = ψxy + 0.25 * δ * (θx * θy * term1 - 2.0 * (θx^2 - θy^2) * term2 + θx * θy * (cos_value - 9.0 * cos_3value)) / r
   return ψxx_up, ψyy_up, ψxy_up
end


"""
    jacobian!(ψxx::ROA, ψyy::ROA, ψxy::ROA, θx::ROA, θy::ROA, δ::T, ϕ::T) where T <: Real
Calculate jacobian at given coordinates for "restricted" third order perturbations corresponding
to SIS lens model and update the jacobian in place.

# Arguments
- `ψxx`: xx-component of jacobian at given coordinates
- `ψyy`: yy-component of jacobian at given coordinates
- `ψxy`: xy-component of jacobian at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `δ` : Amplitude of third order perturbations.
- `ϕ` : External Shear angle (in ``\\rm \\mathbf{degrees}``).
"""
function jacobian!(ψxx::ROA, ψyy::ROA, ψxy::ROA, θx::ROA, θy::ROA, δ::T, ϕ::T) where T <: Real
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

         # cos3(θ - ϕ) and sin3(θ - ϕ)
         cos_3value = 4.0 * cos_value^3 - 3.0 * cos_value
         sin_3value = 3.0 * sin_value - 4.0 * sin_value^3

         # Jacobian
         term1 = 3.0 * (cos_value - cos_3value)
         term2 = (sin_value - 3.0 * sin_3value)
         ψxx[i, j] = ψxx[i, j] + 0.25 * δ * ((r^2 + dx^2) * term1 + 4.0 * dx * dy * term2 - dy^2 * (cos_value - 9.0 * cos_3value)) / r
         ψyy[i, j] = ψyy[i, j] + 0.25 * δ * ((r^2 + dy^2) * term1 - 4.0 * dx * dy * term2 - dx^2 * (cos_value - 9.0 * cos_3value)) / r
         ψxy[i, j] = ψxy[i, j] + 0.25 * δ * (dx * dy * term1 - 2.0 * (dx^2 - dy^2) * term2 + dx * dy * (cos_value - 9.0 * cos_3value)) / r
      end
   end
end

end