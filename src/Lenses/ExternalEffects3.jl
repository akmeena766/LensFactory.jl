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


# --------------------------------------------------------------------------------------------------
# Main functions
# --------------------------------------------------------------------------------------------------
"""
    potential!(ψ::U, θx::S, θy::S, δ::T, ϕ::T) where {U<:Real, S<:Real, T<:Real}
"""
function potential!(ψ::U, θx::S, θy::S, δ::T, ϕ::T) where {U<:Real, S<:Real, T<:Real}
   # Precompute angles
   ϕ = deg2rad(ϕ)
   cos_phi = cos(ϕ)
   sin_phi = sin(ϕ)
   
   # Rotated coordinates
   u = + θx * cos_phi + θy * sin_phi
   v = - θx * sin_phi + θy * cos_phi

   ψ_up = ψ + δ * u * v^2
   return ψ_up
end

"""
    potential!(ψ::U, θx::S, θy::S, δ::T, ϕ::T) where {U<:ROA, S<:ROA, T<:Real}
Calculate potential at given coordinates for "restricted" third order perturbations corresponding
to SIS lens model and update the potential in place.

# Arguments
- `ψ` : Potential at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `δ` : Amplitude of third order perturbations.
- `ϕ` : External Shear angle (in ``\\rm \\mathbf{degrees}``).
"""
function potential!(ψ::U, θx::S, θy::S, δ::T, ϕ::T) where {U<:ROA, S<:ROA, T<:Real}
   # Precompute angles
   ϕ = deg2rad(ϕ)
   cos_phi = cos(ϕ)
   sin_phi = sin(ϕ)
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Rotated coordinates
         u = + θx[i, j] * cos_phi + θy[i, j] * sin_phi
         v = - θx[i, j] * sin_phi + θy[i, j] * cos_phi
 
         ψ[i, j] = ψ[i, j] + δ * u * v^2
      end
   end
end


"""
    deflection!(ψx::U, ψy::U, θx::S, θy::S, δ::T, ϕ::T) where {U<:Real, S<:Real, T<:Real}
"""
function deflection!(ψx::U, ψy::U, θx::S, θy::S, δ::T, ϕ::T) where {U<:Real, S<:Real, T<:Real}
   # Precompute angles
   ϕ_rad = deg2rad(ϕ)
   cos_phi = cos(ϕ_rad)
   sin_phi = sin(ϕ_rad)
 
   # Rotated coordinates
   u = + θx * cos_phi + θy * sin_phi
   v = - θx * sin_phi + θy * cos_phi
 
   # Deflection in the rotated frame
   ψu = δ * v^2
   ψv = 2.0 * δ * u * v
 
   # Rotate back
   ψx_up = ψx + ψu * cos_phi - ψv * sin_phi
   ψy_up = ψy + ψu * sin_phi + ψv * cos_phi
   return ψx_up, ψy_up
end


"""
    deflection!(ψx::U, ψy::U, θx::S, θy::S, δ::T, ϕ::T) where {U<:ROA, S<:ROA, T<:Real}
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
function deflection!(ψx::U, ψy::U, θx::S, θy::S, δ::T, ϕ::T) where {U<:ROA, S<:ROA, T<:Real}
   # Precompute angles
   ϕ_rad = deg2rad(ϕ)
   cos_phi = cos(ϕ_rad)
   sin_phi = sin(ϕ_rad)
 
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Rotated coordinates
         u = + θx[i, j] * cos_phi + θy[i, j] * sin_phi
         v = - θx[i, j] * sin_phi + θy[i, j] * cos_phi
 
         # Deflection in the rotated frame
         ψu = δ * v^2
         ψv = 2.0 * δ * u * v
 
         # Rotate back
         ψx[i, j] = ψx[i, j] + ψu * cos_phi - ψv * sin_phi
         ψy[i, j] = ψy[i, j] + ψu * sin_phi + ψv * cos_phi
      end
   end
   return nothing
end


"""
    jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, δ::T, ϕ::T) where {U<:Real, S<:Real, T<:Real}
"""
function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, δ::T, ϕ::T) where {U<:Real, S<:Real, T<:Real}
   # Precompute angles
   ϕ_rad = deg2rad(ϕ)
   cos_phi = cos(ϕ_rad)
   sin_phi = sin(ϕ_rad)
   sin_2phi = sin(2.0 * ϕ_rad)
   cos_2phi = cos(2.0 * ϕ_rad)
 
   # Rotated coordinates
   u = + θx * cos_phi + θy * sin_phi
   v = - θx * sin_phi + θy * cos_phi
 
   # Jacobian in the rotated frame (ψ_uu = 0)
   ψvv = 2.0 * δ * u
   ψuv = 2.0 * δ * v
 
   # Rotate back
   ψxx_up = ψxx + ψvv * sin_phi^2 - ψuv * sin_2phi
   ψyy_up = ψyy + ψvv * cos_phi^2 + ψuv * sin_2phi
   ψxy_up = ψxy - 0.5 * sin_2phi * ψvv + cos_2phi * ψuv
   return ψxx_up, ψyy_up, ψxy_up
end


"""
    jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, δ::T, ϕ::T) where {U<:ROA, S<:ROA, T<:Real}
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
function jacobian!(ψxx::U, ψyy::U, ψxy::U, θx::S, θy::S, δ::T, ϕ::T) where {U<:ROA, S<:ROA, T<:Real}
   # Precompute angles
   ϕ_rad = deg2rad(ϕ)
   cos_phi = cos(ϕ_rad)
   sin_phi = sin(ϕ_rad)
   sin_2phi = sin(2.0 * ϕ_rad)
   cos_2phi = cos(2.0 * ϕ_rad)
 
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         # Rotated coordinates
         u = + θx[i, j] * cos_phi + θy[i, j] * sin_phi
         v = - θx[i, j] * sin_phi + θy[i, j] * cos_phi
 
         # Jacobian in the rotated frame (ψ_uu = 0)
         ψvv = 2.0 * δ * u
         ψuv = 2.0 * δ * v
 
         # Rotate back
         ψxx[i, j] = ψxx[i, j] + ψvv * sin_phi^2 - ψuv * sin_2phi
         ψyy[i, j] = ψyy[i, j] + ψvv * cos_phi^2 + ψuv * sin_2phi
         ψxy[i, j] = ψxy[i, j] - 0.5 * sin_2phi * ψvv + cos_2phi * ψuv
      end
   end
   return nothing
end

end