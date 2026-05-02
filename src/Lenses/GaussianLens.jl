module GaussianLens

# Inbuilt packages to use
using SpecialFunctions

# LensFactory modules to import
using ..Constants

# Functions to export
export potential!
export deflection!
export jacobian!

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   κs = (2.0 * CONST_G * mass * MASS_SUN / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   κs = κs * θs^2

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr2 = dx^2 + dy^2

   ψ_up = ψ + κs * (log(dr2) - expinti(-0.5 * dr2))
   return ψ_up
end

"""
    potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
Calculate potential at given coordinates for a Gaussian lens and update the potential (ψ) in place.
The lensing potential is given as,
```math
ψ(θ_x, θ_y) = 2 \\, κ_s \\, θ_s^2 \\left[ \\ln\\left( \\frac{|\\pmb{θ} - \\pmb{θ}_c|}{θ_s} \\right) 
      - \\frac{1}{2} \\, \\mathrm{Ei} \\left(- \\frac{|\\pmb{θ} - \\pmb{θ}_c|^2}{2 \\, θ_s^2} \\right) \\right],
```
where ``\\mathrm{Ei}(x)`` is the exponential integral function.

# Arguments
- `ψ`: Potential at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `D_d::RV`: ADD from the observer to the lens (in ``\\rm \\mathbf{meters}``).
- `θxc::RV`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc::RV`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::RV`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `θs::RV`: Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `nothing`: Updates the potential (ψ) in place.
"""
function potential!(ψ::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κs = (2.0 * CONST_G * mass * MASS_SUN / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   κs = κs * θs^2

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr2 = dx^2 + dy^2

         ψ[i, j] = ψ[i, j] + κs * (log(dr2) - expinti(-0.5 * dr2))
      end
   end
end


"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   κs = (2.0 * CONST_G * mass * MASS_SUN / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   κs = 2.0 * κs * θs

   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dr2 = dx^2 + dy^2

   ψx_up = ψx + κs * (1.0 - exp(-0.5 * dr2)) * dx / dr2
   ψy_up = ψy + κs * (1.0 - exp(-0.5 * dr2)) * dy / dr2
   return ψx_up, ψy_up
end

"""
    deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
Calculate deflection at given coordinates for a Gaussian lens and update the deflection components
(ψx, ψy) in place. The deflection angle components are given as,
```math
\\begin{align*} 
ψ_x(θ_x, θ_y) &= 2 \\, κ_s \\, θ_s \\left[ 1 - \\exp\\left(- \\frac{|\\pmb{θ} - 
               \\pmb{θ}_c|^2}{2 \\, θ_s^2} \\right) \\right] \\frac{θ_x - θ_{x,c}}{|\\pmb{θ} - \\pmb{θ}_c|^2}, \\\\
ψ_y(θ_x, θ_y) &= 2 \\, κ_s \\, θ_s \\left[ 1 - \\exp\\left(- \\frac{|\\pmb{θ} - 
               \\pmb{θ}_c|^2}{2 \\, θ_s^2} \\right) \\right] \\frac{θ_y - θ_{y,c}}{|\\pmb{θ} - \\pmb{θ}_c|^2}.
\\end{align*} 
```

# Arguments
- `ψx`: x-component of deflection at given coordinates
- `ψy`: y-component of deflection at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `D_d::RV`: ADD from the observer to the lens (in ``\\rm \\mathbf{meters}``).
- `θxc::RV`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc::RV`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::RV`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `θs::RV`: Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `nothing`: Updates the deflection (ψx, ψy) in place.
"""
function deflection!(ψx::T, ψy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κs = (2.0 * CONST_G * mass * MASS_SUN / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   κs = 2.0 * κs * θs

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dr2 = dx^2 + dy^2

         ψx[i, j] = ψx[i, j] + κs * (1.0 - exp(-0.5 * dr2)) * dx / dr2
         ψy[i, j] = ψy[i, j] + κs * (1.0 - exp(-0.5 * dr2)) * dy / dr2
      end
   end
end


"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: RV
   κs = (2.0 * CONST_G * mass * MASS_SUN / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   
   dx = (θx - θxc) / θs
   dy = (θy - θyc) / θs
   dx2 = dx^2
   dy2 = dy^2
   dr = sqrt(dx2 + dy2)
   dr2 = dr^2
   dr3 = dr * dr2

   exp_term = exp(-0.5 * dr2)
   α_r = κs * 2.0 *  (1.0 - exp_term) / dr
   κ_r = κs * exp_term

   ψxx_up = ψxx + 2.0 * κ_r * dx2 / dr2 - α_r * (dx2 - dy2) / dr3
   ψyy_up = ψyy + 2.0 * κ_r * dy2 / dr2 + α_r * (dx2 - dy2) / dr3
   ψxy_up = ψxy + 2.0 * (κ_r - α_r / dr) * dx * dy / dr2
   return ψxx_up, ψyy_up, ψxy_up
end

"""
    jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
Calculate Jacobian at given coordinates for a Gaussian lens and update the jacobian components 
(ψxx, ψyy, ψxy) in place. The jacobian components are given as.


# Arguments
- `ψxx`: x-component of Jacobian at given coordinates
- `ψyy`: y-component of Jacobian at given coordinates
- `ψxy`: xy-component of Jacobian at given coordinates
- `θx`: x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`: y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `D_d::RV`: ADD from the observer to the lens (in ``\\rm \\mathbf{meters}``).
- `θxc::RV`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `θyc::RV`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::RV`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `θs::RV`: Scale radius i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `nothing`: Updates the jacobian (ψxx, ψyy, ψxy) in place.
"""
function jacobian!(ψxx::T, ψyy::T, ψxy::T, θx::T, θy::T, D_d::RV, θxc::RV, θyc::RV, mass::RV, θs::RV) where T <: ROA
   κs = (2.0 * CONST_G * mass * MASS_SUN / CONST_C^2) / (D_d * θs^2 * ANGLE_ARCSEC^2)
   
   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         dx = (θx[i, j] - θxc) / θs
         dy = (θy[i, j] - θyc) / θs
         dx2 = dx^2
         dy2 = dy^2
         dr = sqrt(dx2 + dy2)
         dr2 = dr^2
         dr3 = dr * dr2

         exp_term = exp(-0.5 * dr2)
         α_r = κs * 2.0 * (1.0 - exp_term) / dr
         κ_r = κs * exp_term

         ψxx[i, j] = ψxx[i, j] + 2.0 * κ_r * dx2 / dr2 - α_r * (dx2 - dy2) / dr3
         ψyy[i, j] = ψyy[i, j] + 2.0 * κ_r * dy2 / dr2 + α_r * (dx2 - dy2) / dr3
         ψxy[i, j] = ψxy[i, j] + 2.0 * (κ_r - α_r / dr) * dx * dy / dr2
      end
   end
end

end