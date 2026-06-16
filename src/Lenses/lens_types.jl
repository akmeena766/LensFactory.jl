export init_PointLens
export init_PlummerLens
export init_SISLens
export init_NSISPLens
export init_NSISMDLens
export init_GaussianLens
export init_SersicLens
export init_ExternalEffects
export init_ExternalEffects3
export init_PIEPLens
export init_SIELens
export init_PJELens
export init_HernquistLens
export init_NFWLens
export init_tNFWLens
export init_gNFWLens
export init_EinastoLens
export init_MultiPlummerLens
export init_MultiGaussianLens
export init_MultiPJELens
export init_CompositeLens
export init_MultiPlaneLens

abstract type AbstractLens end


"""
    init_PointLens(D_d::Real  = NaN, 
                   x_c::Real  = 0.0, 
                   y_c::Real  = 0.0, 
                   mass::Real = NaN)
Initialize a point lens with the given parameters.

# Keyword Arguments
- `D_d::Real = NaN`: ADD from observer to lens (in ``\\rm \\mathbf{meters}``).
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
"""
@kwdef struct init_PointLens <: AbstractLens
   _lens_::Symbol = :PointLens
   D_d::Real  = NaN
   x_c::Real  = 0.0
   y_c::Real  = 0.0
   mass::Real = NaN
end


"""
    init_PlummerLens(D_d::Real  = NaN, 
                     x_c::Real  = 0.0, 
                     y_c::Real  = 0.0, 
                     mass::Real = NaN, 
                     x_s::Real  = NaN)
Initialize a Plummer lens with the given parameters.

# Keyword Arguments
- `D_d::Real = NaN`: ADD from observer to lens (in ``\\rm \\mathbf{meters}``).
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s::Real = NaN`: Core radius (in ``\\rm \\mathbf{arcseconds}``).
"""
@kwdef struct init_PlummerLens <: AbstractLens
   _lens_::Symbol = :PlummerLens
   D_d::Real  = NaN
   x_c::Real  = 0.0
   y_c::Real  = 0.0
   mass::Real = NaN
   x_s::Real  = NaN
end


"""
    init_SISLens(x_c::Real = 0.0, 
                 y_c::Real = 0.0, 
                 v_d::Real = NaN)
Initialize a Singular Isothermal Sphere (SIS) lens with the given parameters.

# Keyword Arguments
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d::Real = NaN`: Velocity dispersion (in ``\\rm \\mathbf{km/s}``).
"""
@kwdef struct init_SISLens <: AbstractLens
   _lens_::Symbol = :SISLens
   x_c::Real = 0.0
   y_c::Real = 0.0
   v_d::Real = NaN
end


"""
    init_NSISPLens(x_c::Real = 0.0, 
                   y_c::Real = 0.0, 
                   v_d::Real = NaN, 
                   x_s::Real = NaN)
Initialize a Non-Singular Isothermal Sphere Potential (NSISP) lens with the given parameters.

# Keyword Arguments
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d::Real = NaN`: Velocity dispersion (in ``\\rm \\mathbf{km/s}``).
- `x_s::Real = NaN`: Core radius (in ``\\rm \\mathbf{arcseconds}``).
"""
@kwdef struct init_NSISPLens <: AbstractLens
   _lens_::Symbol = :NSISPLens
   x_c::Real = 0.0
   y_c::Real = 0.0
   v_d::Real = NaN
   x_s::Real = NaN
end


"""
    init_NSISMDLens(x_c::Real = 0.0, 
                    y_c::Real = 0.0, 
                    v_d::Real = NaN, 
                    x_s::Real = NaN)
Initialize a Non-Singular Isothermal Sphere Mass Distribution (NSISMD) lens with the given parameters.

# Keyword Arguments
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d::Real = NaN`: Velocity dispersion (in ``\\rm \\mathbf{km/s}``).
- `x_s::Real = NaN`: Core radius (in ``\\rm \\mathbf{arcseconds}``).
"""
@kwdef struct init_NSISMDLens <: AbstractLens
   _lens_::Symbol = :NSISMDLens
   x_c::Real = 0.0
   y_c::Real = 0.0
   v_d::Real = NaN
   x_s::Real = NaN
end


"""
    init_GaussianLens(D_d::Real = NaN, 
                      x_c::Real = 0.0, 
                      y_c::Real = 0.0, 
                      mass::Real = NaN, 
                      x_s::Real = NaN)
Initialize a Gaussian lens with the given parameters.

# Keyword Arguments
- `D_d::Real = NaN`: ADD from observer to lens (in ``\\rm \\mathbf{meters}``).
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s::Real = NaN`: Scale radius, i.e., standard deviation of the Gaussian (in ``\\rm \\mathbf{arcseconds}``).
"""
@kwdef struct init_GaussianLens <: AbstractLens
   _lens_::Symbol = :GaussianLens
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   mass::Real= NaN
   x_s::Real = NaN
end


"""
    init_SersicLens(D_d::Real  = NaN, 
                    x_c::Real  = 0.0, 
                    y_c::Real  = 0.0, 
                    mass::Real = NaN, 
                    x_e::Real  = NaN, 
                    n::Real    = 4.0)
Initialize a Sersic lens with the given parameters.

# Keyword Arguments
- `D_d::Real = NaN`: ADD from observer to lens (in ``\\rm \\mathbf{meters}``).
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_e::Real = NaN`: Half-mass radius (in ``\\rm \\mathbf{arcseconds}``).
- `n::Real = 4.0`: Sersic index.
"""
@kwdef struct init_SersicLens <: AbstractLens
   _lens_::Symbol = :SersicLens
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   mass::Real= NaN
   x_e::Real = NaN
   n::Real   = 4.0
end


"""
    init_ExternalEffects(kappa::Real = NaN, 
                         gamma::Real = NaN, 
                         angle::Real = NaN)
Initialize constant external effects with the given parameters.

# Keyword Arguments
- `kappa::Real = NaN`: Convergence (dimensionless).
- `gamma::Real = NaN`: Shear amplitude (dimensionless).
- `angle::Real = NaN`: Shear angle (in ``\\rm \\mathbf{degrees}``).
"""
@kwdef struct init_ExternalEffects <: AbstractLens
   _lens_::Symbol = :ExternalEffects
   kappa::Real = NaN
   gamma::Real = NaN
   angle::Real = NaN
end


"""
    init_ExternalEffects3(delta::Real = NaN, 
                          angle::Real = NaN)
Initialize "restricted" third order perturbations assuming SIS as our perturber.

# Keyword Arguments
- `delta::Real = NaN`: Amplitude of third order perturbations (dimensionless).
- `angle::Real = NaN`: Direction of the perturbation (in ``\\rm \\mathbf{degrees}``).
"""
@kwdef struct init_ExternalEffects3 <: AbstractLens
   _lens_::Symbol = :ExternalEffects3
   delta::Real = NaN
   angle::Real = NaN
end


"""
    init_Multipole(delta::Real = NaN, 
                   angle::Real = NaN, 
                   m::Int64    = 2, 
                   n::Real     = 2.0)
Initialize multipole perturbations of order ``m''.

# Keyword Arguments
- `delta::Real = NaN`: Amplitude of multipole perturbations (dimensionless).
- `angle::Real = NaN`: Direction of the perturbation (in ``\\rm \\mathbf{degrees}``).
- `m::Int64 = 2`: Order of the multipole.
   - `m = 1` → dipole (i.e., constant deflection everywhere)
   - `m = 2` → quadrupole (i.e., external shear)
   - `m = 3` → octupole (i.e., triangular assymmetry)
   - `m = 4` → hexadecapole (i.e., boxiness)

- `n::Real = 2.0`: Controls perturbation scaling with radius.
"""
@kwdef struct init_Multipole <: AbstractLens
   _lens_::Symbol = :Multipole
   delta::Real = NaN
   angle::Real = NaN
   m::Int64 = 2
   n::Real = 2.0
end


"""
    init_PIEPLens(x_c::Real = 0.0, 
                  y_c::Real = 0.0, 
                  v_d::Real = NaN, 
                  x_s::Real = NaN, 
                  eps::Real = NaN, 
                  pa::Real  = NaN)
Initialize pseudo isothermal elliptical potential (PIEP) lens with the given parameters.

# Keyword Arguments
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d::Real = NaN`: Velocity dispersion (in ``\\rm \\mathbf{km/s}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `eps::Real = NaN`: Ellipticity (dimensionless).
- `pa::Real = NaN`: Position angle (in ``\\rm \\mathbf{radians}``).
"""
@kwdef struct init_PIEPLens <: AbstractLens
   _lens_::Symbol = :PIEPLens
   x_c::Real = 0.0
   y_c::Real = 0.0
   v_d::Real = NaN
   x_s::Real = NaN
   eps::Real = NaN
   pa::Real  = NaN
end


"""
    init_SIELens(x_c::Real = 0.0, 
                 y_c::Real = 0.0, 
                 v_d::Real = NaN, 
                 x_s::Real = NaN, 
                 eps::Real = NaN, 
                 pa::Real  = NaN)
Initialize singular isothermal ellipsoid (SIE) lens with the given parameters.

# Keyword Arguments
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d::Real = NaN`: Velocity dispersion (in ``\\rm \\mathbf{km/s}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `eps::Real = NaN`: Ellipticity (dimensionless).
- `pa::Real = NaN`: Position angle (in ``\\rm \\mathbf{radians}``).
"""
@kwdef struct init_SIELens <: AbstractLens
   _lens_::Symbol = :SIELens
   x_c::Real = 0.0
   y_c::Real = 0.0
   v_d::Real = NaN
   x_s::Real = NaN
   eps::Real = NaN
   pa::Real  = NaN
end


"""
    init_PJEMDLens(x_c::Real = 0.0, 
                   y_c::Real = 0.0, 
                   v_d::Real = NaN, 
                   x_s::Real = NaN, 
                   x_t::Real = NaN, 
                   eps::Real = NaN, 
                   pa::Real  = NaN)
Initialize Pseudo-Jaffe Ellipsoid (PJE) lens with the given parameters.

# Keyword Arguments
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `v_d::Real = NaN`: Velocity dispersion (in ``\\rm \\mathbf{km/s}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `x_t::Real = NaN`: Truncation radius (in ``\\rm \\mathbf{arcseconds}``).
- `eps::Real = NaN`: Ellipticity (dimensionless).
- `pa::Real = NaN`: Position angle (in ``\\rm \\mathbf{radians}``).
"""
@kwdef struct init_PJELens <: AbstractLens
   _lens_::Symbol = :PJELens
   x_c::Real = 0.0
   y_c::Real = 0.0
   v_d::Real = NaN
   x_s::Real = NaN
   x_t::Real = NaN
   eps::Real = NaN
   pa::Real  = NaN
end


"""
    init_HernquistLens(D_d::Real  = NaN, 
                       x_c::Real  = 0.0, 
                       y_c::Real  = 0.0, 
                       mass::Real = NaN, 
                       x_s::Real  = NaN)
Initialize a Hernquist lens with the given parameters.

# Keyword Arguments
- `D_d::Real = NaN`: ADD from observer to lens (in ``\\rm \\mathbf{meters}``).
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
"""
@kwdef struct init_HernquistLens <: AbstractLens
   _lens_::Symbol = :HernquistLens
   D_d::Real  = NaN
   x_c::Real  = 0.0
   y_c::Real  = 0.0
   mass::Real = NaN
   x_s::Real  = NaN
end


"""
    init_eHernquistMDLens(D_d::Real  = NaN, 
                          x_c::Real  = 0.0, 
                          y_c::Real  = 0.0, 
                          mass::Real = NaN, 
                          x_s::Real  = NaN, 
                          eps::Real  = NaN, 
                          pa::Real   = NaN)
Initialize an elliptical Hernquist mass distribution lens (eHernquistMDLens) with the given 
parameters.

# Keyword Arguments
- `D_d::Real = NaN`: ADD from observer to lens (in ``\\rm \\mathbf{meters}``).
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `eps::Real = NaN`: Ellipticity.
- `pa::Real = NaN`: Position angle (in ``\\rm \\mathbf{radians}``).
"""
@kwdef struct init_eHernquistMDLens <: AbstractLens
   _lens_::Symbol = :eHernquistMDLens
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   mass::Real= NaN
   x_s::Real = NaN
   eps::Real = NaN
   pa::Real  = NaN
end


"""
    init_aHernquistLens(D_d::Real  = NaN, 
                        x_c::Real  = 0.0, 
                        y_c::Real  = 0.0, 
                        mass::Real = NaN, 
                        x_s::Real  = NaN, 
                        eps::Real  = NaN, 
                        pa::Real   = NaN)
Initialize an approximate Hernquist lens (aHernquistLens) with the given parameters based on 
[Oguri (2021)](https://ui.adsabs.harvard.edu/abs/2021PASP..133g4504O/abstract).

# Keyword Arguments
- `D_d::Real = NaN`: ADD from observer to lens (in ``\\rm \\mathbf{meters}``).
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `eps::Real = NaN`: Ellipticity.
- `pa::Real = NaN`: Position angle (in ``\\rm \\mathbf{radians}``).
"""
@kwdef struct init_aHernquistLens <: AbstractLens
   _lens_::Symbol = :aHernquistLens
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   mass::Real = NaN
   x_s::Real  = NaN
   eps::Real  = NaN
   pa::Real   = NaN
end


"""
    init_MultiPlummerLens(n::Int64  = NaN, 
                          D_d::Real = NaN, 
                          x_c  = Vector{<:Real}, 
                          y_c  = Vector{<:Real}, 
                          mass = Vector{<:Real}, 
                          x_s  = Vector{<:Real})
Initialize a Multi-component Plummer lens with the given parameters.

# Keyword Arguments
- `n::Int64 = NaN`: Number of components.
- `D_d::Real = NaN`: Angular diameter distance to the lens (in ``\\rm \\mathbf{arcseconds}``).
- `x_c = Vector{<:Real}()`: Vector of x-coordinates (in ``\\rm \\mathbf{arcseconds}``).
- `y_c = Vector{<:Real}()`: Vector of y-coordinates (in ``\\rm \\mathbf{arcseconds}``).
- `mass= Vector{<:Real}()`: Vector of masses (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s = Vector{<:Real}()`: Vector of scale radii (in ``\\rm \\mathbf{arcseconds}``).
"""
@kwdef struct init_MultiPlummerLens <: AbstractLens
   _lens_::Symbol = :MultiPlummerLens
   D_d::Real  = NaN
   n::Int64 = NaN
   x_c  = Vector{Real}()
   y_c  = Vector{Real}() 
   mass = Vector{Real}()
   x_s  = Vector{Real}()
end


"""
    init_MultiGaussianLens(n::Int64  = NaN, 
                           D_d::Real = NaN, 
                           x_c  = Vector{<:Real}, 
                           y_c  = Vector{<:Real}, 
                           mass = Vector{<:Real}, 
                           x_s  = Vector{<:Real})
Initialize a Multi-component Gaussian lens with the given parameters.

# Keyword Arguments
- `n::Int64 = NaN`: Number of components.
- `D_d::Real = NaN`: Angular diameter distance to the lens (in ``\\rm \\mathbf{arcseconds}``).
- `x_c = Vector{<:Real}()`: Vector of x-coordinates (in ``\\rm \\mathbf{arcseconds}``).
- `y_c = Vector{<:Real}()`: Vector of y-coordinates (in ``\\rm \\mathbf{arcseconds}``).
- `mass= Vector{<:Real}()`: Vector of masses (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s = Vector{<:Real}()`: Vector of scale radii (in ``\\rm \\mathbf{arcseconds}``).
"""
@kwdef struct init_MultiGaussianLens <: AbstractLens
   _lens_::Symbol = :MultiGaussianLens
   D_d::Real  = NaN
   n::Int64 = NaN
   x_c  = Vector{Real}()
   y_c  = Vector{Real}() 
   mass = Vector{Real}()
   x_s  = Vector{Real}()
end


"""
    init_MultiPJELens(n::Int64 = NaN, 
                      x_c = Vector{<:Real}, 
                      y_c = Vector{<:Real}, 
                      v_d = Vector{<:Real}, 
                      x_s = Vector{<:Real}, 
                      x_t = Vector{<:Real}, 
                      eps = Vector{<:Real}, 
                      pa  = Vector{<:Real})
Initialize a Multi-component PJE lens with the given parameters.

# Keyword Arguments
- `n::Int64 = NaN`: Number of components.
- `x_c = Vector{<:Real}()`: Vector of x-coordinates (in ``\\rm \\mathbf{arcseconds}``).
- `y_c = Vector{<:Real}()`: Vector of y-coordinates (in ``\\rm \\mathbf{arcseconds}``).
- `v_d = Vector{<:Real}()`: Vector of velocity dispersions (in ``\\rm \\mathbf{km/s}``).
- `x_s = Vector{<:Real}()`: Vector of scale radii (in ``\\rm \\mathbf{arcseconds}``).
- `x_t = Vector{<:Real}()`: Vector of tidal radii (in ``\\rm \\mathbf{arcseconds}``).
- `eps = Vector{<:Real}()`: Vector of ellipticities.
- `pa = Vector{<:Real}()`: Vector of position angles (in ``\\rm \\mathbf{radians}``).
"""
@kwdef struct init_MultiPJELens <: AbstractLens
   _lens_::Symbol = :MultiPJELens
   n::Int64 = NaN
   x_c  = Vector{Real}()
   y_c  = Vector{Real}() 
   v_d  = Vector{Real}()
   x_s  = Vector{Real}()
   x_t  = Vector{Real}()
   eps  = Vector{Real}()
   pa   = Vector{Real}()
end


# --------------------------------------------------------------------------------------------------
# Lens models that need constructor for various lenses
# --------------------------------------------------------------------------------------------------
"""
    init_NFWLens(cosmology::AbstractCosmology, z_d::Real; 
                 x_c::Real  = 0.0, 
                 y_c::Real  = 0.0, 
                 mass::Real = NaN, 
                 x_s::Real  = NaN, 
                 c::Real    = NaN)
Initialize a Navarro-Frenk-White (NFW) lens with the given parameters. The lens model can be 
initialized with either the concentration `c` or the scale radius `x_s`. **If both are provided, 
`c` will be used to calculate `x_s` and the input `x_s` will be overwritten.**

# Arguments
- `cosmology::AbstractCosmology`: Cosmology object.
- `z_d::Real`: Redshift of the lens.

# Keyword Arguments
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `c::Real = NaN`: Concentration of the lens.
"""
@kwdef struct init_NFWLens <: AbstractLens
   _lens_::Symbol = :NFWLens
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   k_s::Real = NaN
   x_s::Real = NaN
end

function init_NFWLens(cosmology::Cosmology.AbstractCosmology, z_d::Real; x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN, x_s::Real=NaN, c::Real=NaN)
   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift (in kg/m^3)
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass * MASS_SUN / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide concentration (c) or scale radius (x_s) in **parameter_NFWLens**."))
   end
   # 3D characteristic density
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / (log(1.0 + c) - (c / (1.0 + c)))

   # 2D (normalized) characteristic density
   k_s = ρ_s * D_d * x_s * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)

   return init_NFWLens(x_c=x_c, y_c=y_c, D_d=D_d, k_s=k_s, x_s=x_s)
end


"""
    init_eNFWMDLens(cosmology::AbstractCosmology, z_d::Real; 
                    x_c::Real  = 0.0, 
                    y_c::Real  = 0.0, 
                    mass::Real = NaN, 
                    x_s::Real  = NaN, 
                    c::Real    = NaN, 
                    eps::Real  = NaN, 
                    pa::Real   = NaN)
Initialize an elliptical Navarro-Frenk-White mass distribution lens (eNFWMDLens) with the given 
parameters. The lens model can be initialized with either the concentration `c` or the scale radius 
`x_s`. **If both are provided, `c` will be used to calculate `x_s` and the input `x_s` will be 
overwritten.**

# Arguments
- `cosmology::AbstractCosmology`: Cosmology object.
- `z_d::Real`: Redshift of the lens.

# Keyword Arguments
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `c::Real = NaN`: Concentration of the lens.
- `eps::Real = NaN`: Ellipticity.
- `pa::Real = NaN`: Position angle (in ``\\rm \\mathbf{radians}``).
"""
@kwdef struct init_eNFWMDLens <: AbstractLens
   _lens_::Symbol = :eNFWMDLens
   D_d::Real  = NaN
   x_c::Real  = 0.0
   y_c::Real  = 0.0
   k_s::Real  = NaN
   x_s::Real  = NaN
   eps::Real  = NaN
   pa::Real   = NaN
end

function init_eNFWMDLens(cosmology::Cosmology.AbstractCosmology, z_d::Real; x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN, x_s::Real=NaN, c::Real=NaN, eps::Real=NaN, pa::Real=NaN)
   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift (in kg/m^3)
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass * MASS_SUN / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide concentration (c) or scale radius (x_s) in **parameter_NFWLens**."))
   end
   # 3D characteristic density
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / (log(1.0 + c) - (c / (1.0 + c)))

   # 2D (normalized) characteristic density
   k_s = ρ_s * D_d * x_s * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)

   return init_eNFWMDLens(x_c=x_c, y_c=y_c, D_d=D_d, k_s=k_s, x_s=x_s, eps=eps, pa=pa)
end


"""
    init_aNFWLens(cosmology::AbstractCosmology, z_d::Real; 
                  x_c::Real  = 0.0, 
                  y_c::Real  = 0.0, 
                  mass::Real = NaN, 
                  x_s::Real  = NaN, 
                  c::Real    = NaN, 
                  eps::Real  = NaN, 
                  pa::Real   = NaN)
Initialize an approximate Navarro-Frenk-White lens (aNFWLens) with the given parameters based on 
[Oguri (2021)](https://ui.adsabs.harvard.edu/abs/2021PASP..133g4504O/abstract). The lens model can 
be initialized with either the concentration `c` or the scale radius `x_s`. **If both are provided, 
`c` will be used to calculate `x_s` and the input `x_s` will be overwritten.**

# Arguments
- `cosmology::AbstractCosmology`: Cosmology object.
- `z_d::Real`: Redshift of the lens.

# Keyword Arguments
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `c::Real = NaN`: Concentration of the lens.
- `eps::Real = NaN`: Ellipticity.
- `pa::Real = NaN`: Position angle (in ``\\rm \\mathbf{radians}``).
"""
@kwdef struct init_aNFWLens <: AbstractLens
   _lens_::Symbol = :aNFWLens
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   k_s::Real = NaN
   x_s::Real  = NaN
   eps::Real  = NaN
   pa::Real   = NaN
end

function init_aNFWLens(cosmology::Cosmology.AbstractCosmology, z_d::Real; x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN, x_s::Real=NaN, c::Real=NaN, eps::Real=NaN, pa::Real=NaN)
   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift (in kg/m^3)
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass * MASS_SUN / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide concentration (c) or scale radius (x_s) in **init_aNFWLens**."))
   end
   # 3D characteristic density
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / (log(1.0 + c) - (c / (1.0 + c)))

   # 2D (normalized) characteristic density
   k_s = ρ_s * D_d * x_s * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)

   return init_aNFWLens(x_c=x_c, y_c=y_c, D_d=D_d, k_s=k_s, x_s=x_s, eps=eps, pa=pa)
end


"""
    init_tNFWLens(cosmology::AbstractCosmology, z_d::Real; 
                  x_c::Real  = 0.0, 
                  y_c::Real  = 0.0, 
                  mass::Real = NaN, 
                  x_s::Real  = NaN, 
                  c::Real    = NaN, 
                  x_t::Real  = NaN)
Initialize a truncated Navarro-Frenk-White (tNFW) lens with the given parameters. The lens model can
be initialized with either the concentration `c` or the scale radius `x_s`. **If both are provided, 
`c` will be used to calculate `x_s` and the input `x_s` will be overwritten.**

# Arguments
- `cosmology::AbstractCosmology`: Cosmology object.
- `z_d::Real`: Redshift of the lens.

# Keyword Arguments
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `c::Real = NaN`: Concentration of the lens.
- `x_t::Real = NaN`: Truncation radius (in ``\\rm \\mathbf{arcseconds}``).
"""
@kwdef struct init_tNFWLens <: AbstractLens
   _lens_::Symbol = :tNFWLens
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   k_s::Real = NaN
   x_s::Real = NaN
   x_t::Real = NaN
end

function init_tNFWLens(cosmology::Cosmology.AbstractCosmology, z_d::Real; x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN, x_s::Real=NaN, c::Real=NaN, x_t::Real=NaN)
   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift (in kg/m^3)
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass * MASS_SUN / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide concentration (c) or scale radius (x_s) in **parameter_NFWLens**."))
   end
   # 3D characteristic density
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / (log(1.0 + c) - (c / (1.0 + c)))

   # 2D (normalized) characteristic density
   k_s = ρ_s * D_d * x_s * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)

   return init_tNFWLens(x_c=x_c, y_c=y_c, D_d=D_d, k_s=k_s, x_s=x_s, x_t=x_t)
end


"""
    init_gNFWLens(cosmology::AbstractCosmology, z_d::Real; 
                  x_c::Real  = 0.0, 
                  y_c::Real  = 0.0, 
                  mass::Real = NaN, 
                  x_s::Real  = NaN, 
                  c::Real    = NaN, 
                  n::Real    = 1.0)
Initialize a generalized Navarro-Frenk-White (gNFW) lens with the given parameters. The lens model can
be initialized with either the concentration `c` or the scale radius `x_s`. **If both are provided, 
`c` will be used to calculate `x_s` and the input `x_s` will be overwritten.** The parameter `n` 
defines the slope of the density profile.

# Arguments
- `cosmology::AbstractCosmology`: Cosmology object.
- `z_d::Real`: Redshift of the lens.

# Keyword Arguments
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `c::Real = NaN`: Concentration of the lens.
- `n::Real = 1.0`: Slope parameter of the lens.
"""
@kwdef struct init_gNFWLens <: AbstractLens
   _lens_::Symbol = :gNFWLens
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   k_s::Real = NaN
   x_s::Real = NaN
   n::Real   = NaN
end

function init_gNFWLens(cosmology::Cosmology.AbstractCosmology, z_d::Real; x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN, x_s::Real=NaN, c::Real=NaN, n::Real=1.0)
   # Check for valid slope parameter
   if !(0.0 < n < 2.0)
      throw(ArgumentError("Slope parameter outside allowed range n ∈ (0, 2) in **parameter_gNFWLens**."))
   end

   # Integrand function for mass calculation
   function integrand(x::Real, α::Real)
      return x^(2.0 - α) / (1.0 + x)^(3.0 - α)
   end

   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass * MASS_SUN / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide at least c or x_s in **parameter_gNFWLens**."))
   end
   mass_c, _ = quadgk(x -> integrand(x, n), 0, c)
   
   # 3D characteristic density
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / mass_c

   # 2D (normalized) characteristic density
   k_s = ρ_s * D_d * x_s * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)

   return init_gNFWLens(x_c=x_c, y_c=y_c, D_d=D_d, k_s=k_s, x_s=x_s, n=n)
end


"""
    init_EinastoLens(D_d::Real = NaN, 
                     x_c::Real = 0.0, 
                     y_c::Real = 0.0, 
                     k_s::Real = NaN, 
                     x_s::Real = NaN, 
                     n::Real   = 0.2)
Initialize an Einasto lens with the given parameters. The lens model can be initialized with either
the concentration `c` or the scale radius `x_s`. **If both are provided, `c` will be used to 
calculate `x_s` and the input `x_s` will be overwritten.** The parameter `n` defines the slope of 
the density profile.

# Arguments
- `cosmology::AbstractCosmology`: Cosmology object.
- `z_d::Real`: Redshift of the lens.

# Keyword Arguments
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `c::Real = NaN`: Concentration of the lens.
- `n::Real = 0.2`: Slope parameter of the lens.
"""
@kwdef struct init_EinastoLens <: AbstractLens
   _lens_::Symbol = :EinastoLens
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   k_s::Real = NaN
   x_s::Real = NaN
   n::Real   = 0.2
end

function init_EinastoLens(cosmology::Cosmology.AbstractCosmology, z_d::Real; x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN, x_s::Real=NaN, c::Real=NaN, n::Real=0.2)
   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d  = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)

   # Critical density at the lens redshift
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass * MASS_SUN / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide at least c or x_s in **parameter_EinastoLens**."))
   end
   Pax, _ = gamma_inc(3.0 / n, (2.0 / n) * c^n)
   mass_e = (1.0 / n) * (n / 2.0)^(3.0 / n) * gamma(3.0 / n) * Pax
   
   # 3D characteristic density
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / mass_e

   # 2D (normalized) characteristic density
   k_s = ρ_s * D_d * x_s * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)

   return init_EinastoLens(x_c=x_c, y_c=y_c, D_d=D_d, k_s=k_s, x_s=x_s, n=n)
end


# --------------------------------------------------------------------------------------------------
# Composite and Multi-plane lens constructors
# --------------------------------------------------------------------------------------------------
# Dictionary to map lens types to their initialization functions and arguments
const lens_init_functions = Dict{Symbol, Function}(
   :PointLens         => (comp -> init_PointLens(D_d=comp.D_d, x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass)),
   :PlummerLens       => (comp -> init_PlummerLens(D_d=comp.D_d, x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, x_s=comp.x_s)),
   :SISLens           => (comp -> init_SISLens(x_c=comp.x_c, y_c=comp.y_c, v_d=comp.v_d)),
   :NSISPLens         => (comp -> init_NSISPLens(x_c=comp.x_c, y_c=comp.y_c, v_d=comp.v_d, x_s=comp.x_s)),
   :NSISMDLens        => (comp -> init_NSISMDLens(x_c=comp.x_c, y_c=comp.y_c, v_d=comp.v_d, x_s=comp.x_s)),
   :GaussianLens      => (comp -> init_GaussianLens(D_d=comp.D_d, x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, x_s=comp.x_s)),
   :SersicLens        => (comp -> init_SersicLens(D_d=comp.D_d, x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, x_e=comp.x_e, n=comp.n)),
   :ExternalEffects   => (comp -> init_ExternalEffects(kappa=comp.kappa, gamma=comp.gamma, angle=comp.angle)),
   :ExternalEffects3  => (comp -> init_ExternalEffects3(delta=comp.delta, angle=comp.angle)),
   :Multipole         => (comp -> init_Multipole(delta=comp.delta, angle=comp.angle, m=comp.m, n=comp.n)),
   :PIEPLens          => (comp -> init_PIEPLens(x_c=comp.x_c, y_c=comp.y_c, v_d=comp.v_d, x_s=comp.x_s, eps=comp.eps, pa=comp.pa)),
   :SIELens           => (comp -> init_SIELens(x_c=comp.x_c, y_c=comp.y_c, v_d=comp.v_d, x_s=comp.x_s, eps=comp.eps, pa=comp.pa)),
   :PJELens           => (comp -> init_PJELens(x_c=comp.x_c, y_c=comp.y_c, v_d=comp.v_d, x_s=comp.x_s, x_t=comp.x_t, eps=comp.eps, pa=comp.pa)),
   :HernquistLens     => (comp -> init_HernquistLens(D_d=comp.D_d, x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, x_s=comp.x_s)),
   :aHernquistLens    => (comp -> init_aHernquistLens(D_d=comp.D_d, x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, x_s=comp.x_s, eps=comp.eps, pa=comp.pa)),
   :NFWLens           => (comp -> init_NFWLens(comp.cosmology, comp.z_d; x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, c=comp.c)),
   :aNFWLens          => (comp -> init_aNFWLens(comp.cosmology, comp.z_d; x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, c=comp.c, eps=comp.eps, pa=comp.pa)),
   :eNFWMDLens        => (comp -> init_eNFWMDLens(comp.cosmology, comp.z_d; x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, c=comp.c, eps=comp.eps, pa=comp.pa)),
   :tNFWLens          => (comp -> init_tNFWLens(comp.cosmology, comp.z_d; x_c=comp.x_c, y_c=comp.y_c, x_s=comp.x_s, x_t=comp.x_t)),
   :gNFWLens          => (comp -> init_gNFWLens(comp.cosmology, comp.z_d; x_c=comp.x_c, y_c=comp.y_c, x_s=comp.x_s, n=comp.n)),
   :EinastoLens       => (comp -> init_EinastoLens(comp.cosmology, comp.z_d; x_c=comp.x_c, y_c=comp.y_c, x_s=comp.x_s, n=comp.n)),
   :MultiPlummerLens  => (comp -> init_MultiPlummerLens(D_d=comp.D_d, n=comp.n, x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, x_s=comp.x_s)),
   :MultiGaussianLens => (comp -> init_MultiGaussianLens(D_d=comp.D_d, n=comp.n, x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, x_s=comp.x_s)),
   :MultiPJELens      => (comp -> init_MultiPJELens(n=comp.n, x_c=comp.x_c, y_c=comp.y_c, v_d=comp.v_d, x_s=comp.x_s, x_t=comp.x_t, eps=comp.eps, pa=comp.pa))
   )

"""
    init_CompositeLens(lens::Vector{<:NamedTuple})
Initialize a composite lens from a vector of lens components.

# Arguments
- `lens::Vector{<:NamedTuple}`: Vector of lens components.

# Returns
- `CompositeLens`: Composite lens.
"""
@kwdef struct init_CompositeLens <: AbstractLens
   _lens_::Symbol = :CompositeLens
   _components_ = Vector{AbstractLens}()
end

function init_CompositeLens(lens::Vector{<:NamedTuple})
   # Define a compoent vector of known size
   lens_components = Vector{AbstractLens}(undef, length(lens))

   # Run over the lens components in the composite lens
   for (i, component) in enumerate(lens)
      # Check and get the lens component if available otherwise throw an error
      builder = get(lens_init_functions, component.lens, nothing)
      if builder === nothing
         throw(ArgumentError("Unknown lens type: $(component.lens). Available lens models are $(keys(lens_init_functions))"))
      end

      # Wrapper function call
      lens_components[i] = builder(component)
   end
   return init_CompositeLens(_components_=lens_components)
end


@kwdef struct init_MultiPlaneLens <: AbstractLens
   _lens_::Symbol = :MultiPlaneLens
   n_p::Int64   = NaN
   z_d = Vector{<:Real}()
   _plane_ = Vector{AbstractLens}()
end

function init_MultiPlaneLens(lens::Vector{<:NamedTuple})   
   # Get sorted unique lens redshifts
   zd_unique = unique(component.z_d for component in lens)
   sort!(zd_unique)

   if length(zd_unique) < 2
      throw(ArgumentError("Only $(length(zd_unique)) lens plane found. Need >= 2."))
   end

   # Group lens components by redshift
   lens_by_z = Dict{eltype(getfield.(lens, :z_d)), Vector{NamedTuple}}()
   for component in lens
      push!(get!(lens_by_z, component.z_d, []), component)
   end

   # Construct composite lenses for each unique redshift
   lens_components = AbstractLens[init_CompositeLens(lens_by_z[z]) for z in zd_unique]

   return init_MultiPlaneLens(n_p=length(zd_unique), z_d=zd_unique, _plane_=lens_components)
end


# --------------------------------------------------------------------------------------------------
# Parameter functions for various lenses
# --------------------------------------------------------------------------------------------------
"""
    parameter_NFWLens(; cosmology::Cosmology.AbstractCosmology = nothing, 
                        z_d::Real  = NaN, 
                        mass::Real = NaN, 
                        x_s::Real  = NaN, 
                        c::Real    = NaN)
Calculate parameters for NFW lens. The function would either need the concentration `c` or the 
scale radius `x_s`. **If both are provided, `c` will be used to calculate `x_s` and the input `x_s` 
will be overwritten.**

# Arguments
- `cosmology::AbstractCosmology = nothing`: Cosmology object.
- `z_d::Real = NaN`: Redshift of the lens.
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `c::Real = NaN`: Concentration of the lens.

# Returns
- `NamedTuple`: Tuple of lens parameters.
"""
function parameter_NFWLens(; cosmology::Cosmology.AbstractCosmology = nothing, 
                             z_d::Real  = NaN, 
                             mass::Real = NaN, 
                             x_s::Real  = NaN, 
                             c::Real    = NaN)
   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift (in kg/m^3)
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass * MASS_SUN / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide concentration (c) or scale radius (x_s) in **parameter_NFWLens**."))
   end
   # 3D characteristic density
   mass_c = log(1.0 + c) - (c / (1.0 + c))
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / mass_c

   # 2D (normalized) characteristic density
   k_s = ρ_s * D_d * x_s * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)
   return (mass=mass, rho_s=ρ_s, k_s=k_s, c=c, x_s=x_s)
end


"""
    parameter_gNFWLens(; cosmology::Cosmology.AbstractCosmology = nothing, 
                         z_d::Real  = NaN, 
                         mass::Real = NaN, 
                         x_s::Real  = NaN, 
                         c::Real    = NaN, 
                         n::Real    = 1.0)
Calculate parameters for gNFW lens. The function would either need the concentration `c` or the 
scale radius `x_s`. **If both are provided, `c` will be used to calculate `x_s` and the input `x_s` 
will be overwritten.**

# Arguments
- `cosmology::AbstractCosmology = nothing`: Cosmology object.
- `z_d::Real = NaN`: Redshift of the lens.
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `c::Real = NaN`: Concentration of the lens.
- `n::Real = 1.0`: Slope parameter of the lens.

# Returns
- `NamedTuple`: Tuple of lens parameters.
"""
function parameter_gNFWLens(; cosmology::Cosmology.AbstractCosmology=nothing, 
                              z_d::Real  = NaN, 
                              mass::Real = NaN, 
                              x_s::Real  = NaN, 
                              c::Real    = NaN, 
                              n::Real    = 1.0)
   # Check for valid slope parameter
   if !(0.0 < n < 2.0)
      throw(ArgumentError("Slope parameter outside allowed range n ∈ (0, 2) in **parameter_gNFWLens**."))
   end

   # Integrand function for mass calculation
   function integrand(x::Real, α::Real)
      return x^(2.0 - α) / (1.0 + x)^(3.0 - α)
   end

   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass * MASS_SUN / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide at least c or x_s in **parameter_gNFWLens**."))
   end
   mass_c, _ = quadgk(x -> integrand(x, n), 0, c)
   
   # 3D characteristic density
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / mass_c

   # 2D (normalized) characteristic density
   k_s = ρ_s * D_d * x_s * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)
   return (mass=mass, rho_s=ρ_s, k_s=k_s, c=c, x_s=x_s, n=n)
end


"""
    parameter_EinastoLens(; cosmology::Cosmology.AbstractCosmology=nothing, 
                            z_d::Real  = NaN, 
                            mass::Real = NaN, 
                            x_s::Real  = NaN, 
                            c::Real    = NaN, 
                            n::Real    = 0.2)
Calculate parameters of an Einasto lens model. The function would either need the concentration `c` 
or the scale radius `x_s`. **If both are provided, `c` will be used to calculate `x_s` and the input 
`x_s` will be overwritten.**

# Arguments
- `cosmology::AbstractCosmology = nothing`: Cosmology object.
- `z_d::Real = NaN`: Redshift of the lens.
- `mass::Real= NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s::Real = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `c::Real = NaN`: Concentration of the lens.
- `n::Real = 0.2`: Slope parameter of the lens.

# Returns
- `NamedTuple`: Tuple of lens parameters.
"""
function parameter_EinastoLens(; cosmology::Cosmology.AbstractCosmology=nothing, 
                                 z_d::Real  = NaN, 
                                 mass::Real = NaN, 
                                 x_s::Real  = NaN, 
                                 c::Real    = NaN, 
                                 n::Real    = 0.2)
   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d  = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)

   # Critical density at the lens redshift
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass * MASS_SUN / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide at least c or x_s in **parameter_EinastoLens**."))
   end
   Pax, _ = gamma_inc(3.0 / n, (2.0 / n) * c^n)
   mass_e = (1.0 / n) * (n / 2.0)^(3.0 / n) * gamma(3.0 / n) * Pax
   
   # 3D characteristic density
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / mass_e

   # 2D (normalized) characteristic density
   k_s = ρ_s * D_d * x_s * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)
   return (mass=mass, rho_s=ρ_s, k_s=k_s, c=c, x_s=x_s, n=n)
end
