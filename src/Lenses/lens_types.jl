export AbstractLens
export init_PointLens
export init_PlummerLens
export init_SISLens
export init_NSISPLens
export init_NSISMDLens
export init_GaussianLens
export init_SersicLens
export init_PixelLens
export init_ExternalEffects
export init_ExternalEffects3
export init_Multipole
export init_PIEPLens
export init_SIELens
export init_PJELens
export init_HernquistLens
export init_aHernquistLens
export init_eHernquistMDLens
export init_NFWLens
export init_aNFWLens
export init_eNFWMDLens
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
struct init_PointLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   x_c::T
   y_c::T
   mass::T
end
function init_PointLens(; D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN)
   D_d, x_c, y_c, mass = promote(D_d, x_c, y_c, mass)
   T = typeof(D_d)
   return init_PointLens{T}(:PointLens, D_d, x_c, y_c, mass)
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
struct init_PlummerLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   x_c::T
   y_c::T
   mass::T
   x_s::T
end
function init_PlummerLens(; D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN, x_s::Real=NaN)
   D_d, x_c, y_c, mass, x_s = promote(D_d, x_c, y_c, mass, x_s)
   T = typeof(D_d)
   return init_PlummerLens{T}(:PlummerLens, D_d, x_c, y_c, mass, x_s)
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
struct init_SISLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   x_c::T
   y_c::T
   v_d::T
end
function init_SISLens(; x_c::Real=0.0, y_c::Real=0.0, v_d::Real=NaN)
   x_c, y_c, v_d = promote(x_c, y_c, v_d)
   T = typeof(x_c)
   return init_SISLens{T}(:SISLens, x_c, y_c, v_d)
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
struct init_NSISPLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   x_c::T
   y_c::T
   v_d::T
   x_s::T
end
function init_NSISPLens(; x_c::Real=0.0, y_c::Real=0.0, v_d::Real=NaN, x_s::Real=NaN)
   x_c, y_c, v_d, x_s = promote(x_c, y_c, v_d, x_s)
   T = typeof(x_c)
   return init_NSISPLens{T}(:NSISPLens, x_c, y_c, v_d, x_s)
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
struct init_NSISMDLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   x_c::T
   y_c::T
   v_d::T
   x_s::T
end
function init_NSISMDLens(; x_c::Real=0.0, y_c::Real=0.0, v_d::Real=NaN, x_s::Real=NaN)
   x_c, y_c, v_d, x_s = promote(x_c, y_c, v_d, x_s)
   T = typeof(x_c)
   return init_NSISMDLens{T}(:NSISMDLens, x_c, y_c, v_d, x_s)
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
struct init_GaussianLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   x_c::T
   y_c::T
   mass::T
   x_s::T
end
function init_GaussianLens(; D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN, x_s::Real=NaN)
   D_d, x_c, y_c, mass, x_s = promote(D_d, x_c, y_c, mass, x_s)
   T = typeof(D_d)
   return init_GaussianLens{T}(:GaussianLens, D_d, x_c, y_c, mass, x_s)
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
struct init_SersicLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   x_c::T
   y_c::T
   mass::T
   x_e::T
   n::T
end
function init_SersicLens(; D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN, x_e::Real=NaN, n::Real=4.0)
   D_d, x_c, y_c, mass, x_e, n = promote(D_d, x_c, y_c, mass, x_e, n)
   T = typeof(D_d)
   return init_SersicLens{T}(:SersicLens, D_d, x_c, y_c, mass, x_e, n)
end


"""
    init_PixelLens(x_c::Real = 0.0, 
                   y_c::Real = 0.0, 
                   kappa::Real = NaN)
Initialize a square pixel lens with the given parameters.

# Keyword Arguments
- `x_c::Real = 0.0`: x-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `y_c::Real = 0.0`: y-coordinate of the lens (in ``\\rm \\mathbf{arcseconds}``).
- `kappa::Real = NaN`: Convergence (dimensionless).
"""
struct init_PixelLens{T<:Real} <: AbstractLens  
   _lens_::Symbol
   x_c::T
   y_c::T
   kappa::T
   pixel_size::T
end
function init_PixelLens(; x_c::Real=0.0, y_c::Real=0.0, kappa::Real=NaN, pixel_size::Real=NaN)
   x_c, y_c, kappa = promote(x_c, y_c, kappa, pixel_size)
   T = typeof(x_c)
   return init_PixelLens{T}(:PixelLens, x_c, y_c, kappa, pixel_size)
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
struct init_ExternalEffects{T<:Real} <: AbstractLens
   _lens_::Symbol
   kappa::T
   gamma::T
   angle::T
end
function init_ExternalEffects(; kappa::Real=NaN, gamma::Real=NaN, angle::Real=NaN)
   kappa, gamma, angle = promote(kappa, gamma, angle)
   T = typeof(kappa)
   return init_ExternalEffects{T}(:ExternalEffects, kappa, gamma, angle)
end


"""
    init_ExternalEffects3(delta::Real = NaN, 
                          angle::Real = NaN)
Initialize "restricted" third order perturbations assuming SIS as our perturber.

# Keyword Arguments
- `delta::Real = NaN`: Amplitude of third order perturbations (dimensionless).
- `angle::Real = NaN`: Direction of the perturbation (in ``\\rm \\mathbf{degrees}``).
"""
struct init_ExternalEffects3{T<:Real} <: AbstractLens
   _lens_::Symbol
   delta::T
   angle::T
end
function init_ExternalEffects3(; delta::Real=NaN, angle::Real=NaN)
   delta, angle = promote(delta, angle)
   T = typeof(delta)
   return init_ExternalEffects3{T}(:ExternalEffects3, delta, angle)
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
struct init_Multipole{T<:Real} <: AbstractLens
   _lens_::Symbol
   delta::T
   angle::T
   m::Int64
   n::T
end
function init_Multipole(; delta::Real=NaN, angle::Real=NaN, m::Int64=2, n::Real=2.0)
   delta, angle, n = promote(delta, angle, n)
   T = typeof(delta)
   return init_Multipole{T}(:Multipole, delta, angle, m, n)
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
- `pa::Real = NaN`: Position angle (in ``\\rm \\mathbf{deg}``).
"""
struct init_PIEPLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   x_c::T
   y_c::T
   v_d::T
   x_s::T
   eps::T
   pa::T
end
function init_PIEPLens(; x_c::Real=0.0, y_c::Real=0.0, v_d::Real=NaN, x_s::Real=NaN, eps::Real=NaN, pa::Real=NaN)
   x_c, y_c, v_d, x_s, eps, pa = promote(x_c, y_c, v_d, x_s, eps, pa)
   T = typeof(x_c)
   return init_PIEPLens{T}(:PIEPLens, x_c, y_c, v_d, x_s, eps, pa)
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
- `pa::Real = NaN`: Position angle (in ``\\rm \\mathbf{deg}``).
"""
struct init_SIELens{T<:Real} <: AbstractLens
   _lens_::Symbol
   x_c::T
   y_c::T
   v_d::T
   x_s::T
   eps::T
   pa::T
end
function init_SIELens(; x_c::Real=0.0, y_c::Real=0.0, v_d::Real=NaN, x_s::Real=NaN, eps::Real=NaN, pa::Real=NaN)
   x_c, y_c, v_d, x_s, eps, pa = promote(x_c, y_c, v_d, x_s, eps, pa)
   T = typeof(x_c)
   return init_SIELens{T}(:SIELens, x_c, y_c, v_d, x_s, eps, pa)
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
- `pa::Real = NaN`: Position angle (in ``\\rm \\mathbf{deg}``).
"""
struct init_PJELens{T<:Real} <: AbstractLens
   _lens_::Symbol
   x_c::T
   y_c::T
   v_d::T
   x_s::T
   x_t::T
   eps::T
   pa::T
end
function init_PJELens(; x_c::Real=0.0, y_c::Real=0.0, v_d::Real=NaN, x_s::Real=NaN, x_t::Real=NaN, eps::Real=NaN, pa::Real=NaN)
   x_c, y_c, v_d, x_s, x_t, eps, pa = promote(x_c, y_c, v_d, x_s, x_t, eps, pa)
   T = typeof(x_c)
   return init_PJELens{T}(:PJELens, x_c, y_c, v_d, x_s, x_t, eps, pa)
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
struct init_HernquistLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   x_c::T
   y_c::T
   mass::T
   x_s::T
end
function init_HernquistLens(; D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN, x_s::Real=NaN)
   D_d, x_c, y_c, mass, x_s = promote(D_d, x_c, y_c, mass, x_s)
   T = typeof(D_d)
   return init_HernquistLens{T}(:HernquistLens, D_d, x_c, y_c, mass, x_s)
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
- `pa::Real = NaN`: Position angle (in ``\\rm \\mathbf{deg}``).
"""
struct init_eHernquistMDLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   x_c::T
   y_c::T
   mass::T
   x_s::T
   eps::T
   pa::T
end
function init_eHernquistMDLens(; D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN, x_s::Real=NaN, eps::Real=NaN, pa::Real=NaN)
   D_d, x_c, y_c, mass, x_s, eps, pa = promote(D_d, x_c, y_c, mass, x_s, eps, pa)
   T = typeof(D_d)
   return init_eHernquistMDLens{T}(:eHernquistMDLens, D_d, x_c, y_c, mass, x_s, eps, pa)
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
- `pa::Real = NaN`: Position angle (in ``\\rm \\mathbf{deg}``).
"""
struct init_aHernquistLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   x_c::T
   y_c::T
   mass::T
   x_s::T
   eps::T
   pa::T
end
function init_aHernquistLens(; D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN, x_s::Real=NaN, eps::Real=NaN, pa::Real=NaN)
   D_d, x_c, y_c, mass, x_s, eps, pa = promote(D_d, x_c, y_c, mass, x_s, eps, pa)
   T = typeof(D_d)
   return init_aHernquistLens{T}(:aHernquistLens, D_d, x_c, y_c, mass, x_s, eps, pa)
end


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
struct init_NFWLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   x_c::T
   y_c::T
   k_s::T
   x_s::T
end
function init_NFWLens(; D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, k_s::Real=NaN, x_s::Real=NaN)
   D_d, x_c, y_c, k_s, x_s = promote(D_d, x_c, y_c, k_s, x_s)
   T = typeof(D_d)
   return init_NFWLens{T}(:NFWLens, D_d, x_c, y_c, k_s, x_s)
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

   return init_NFWLens(D_d=D_d, x_c=x_c, y_c=y_c, k_s=k_s, x_s=x_s)
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
- `pa::Real = NaN`: Position angle (in ``\\rm \\mathbf{deg}``).
"""
struct init_eNFWMDLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   x_c::T
   y_c::T
   k_s::T
   x_s::T
   eps::T
   pa::T
end
function init_eNFWMDLens(; D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, k_s::Real=NaN, x_s::Real=NaN, eps::Real=NaN, pa::Real=NaN)
   # promote the input arguments
   D_d, x_c, y_c, k_s, x_s, eps, pa = promote(D_d, x_c, y_c, k_s, x_s, eps, pa)
   T = typeof(D_d)
   return init_eNFWMDLens{T}(:eNFWMDLens, D_d, x_c, y_c, k_s, x_s, eps, pa)
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

   return init_eNFWMDLens(D_d=D_d, x_c=x_c, y_c=y_c, k_s=k_s, x_s=x_s, eps=eps, pa=pa)
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
- `pa::Real = NaN`: Position angle (in ``\\rm \\mathbf{deg}``).
"""
struct init_aNFWLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   x_c::T
   y_c::T
   k_s::T
   x_s::T
   eps::T
   pa::T
end
function init_aNFWLens(; D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, k_s::Real=NaN, x_s::Real=NaN, eps::Real=NaN, pa::Real=NaN)
   D_d, x_c, y_c, k_s, x_s, eps, pa = promote(D_d, x_c, y_c, k_s, x_s, eps, pa)
   T = typeof(D_d)
   return init_aNFWLens{T}(:aNFWLens, D_d, x_c, y_c, k_s, x_s, eps, pa)
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
struct init_tNFWLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   x_c::T
   y_c::T
   k_s::T
   x_s::T
   x_t::T
end
function init_tNFWLens(; D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, k_s::Real=NaN, x_s::Real=NaN, x_t::Real=NaN)
   D_d, x_c, y_c, k_s, x_s, x_t = promote(D_d, x_c, y_c, k_s, x_s, x_t)
   T = typeof(D_d)
   return init_tNFWLens{T}(:tNFWLens, D_d, x_c, y_c, k_s, x_s, x_t)
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
struct init_gNFWLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   x_c::T
   y_c::T
   k_s::T
   x_s::T
   n::T
end
function init_gNFWLens(; D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, k_s::Real=NaN, x_s::Real=NaN, n::Real=1.0)
   if !(0.0 < n < 2.0)
      throw(ArgumentError("Slope parameter outside allowed range n ∈ (0, 2) in **parameter_gNFWLens**."))
   end

   D_d, x_c, y_c, k_s, x_s, n = promote(D_d, x_c, y_c, k_s, x_s, n)
   T = typeof(D_d)
   return init_gNFWLens{T}(:gNFWLens, D_d, x_c, y_c, k_s, x_s, n)
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
struct init_EinastoLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   x_c::T
   y_c::T
   k_s::T
   x_s::T
   n::T
end
function init_EinastoLens(; D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, k_s::Real=NaN, x_s::Real=NaN, n::Real=0.2)
   D_d, x_c, y_c, k_s, x_s, n = promote(D_d, x_c, y_c, k_s, x_s, n)
   T = typeof(D_d)
   return init_EinastoLens{T}(:EinastoLens, D_d, x_c, y_c, k_s, x_s, n)
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


"""
    init_MultiPlummerLens(D_d::Real = NaN, 
                          x_c  = Vector{<:Real}, 
                          y_c  = Vector{<:Real}, 
                          mass = Vector{<:Real}, 
                          x_s  = Vector{<:Real})
Initialize a Multi-component Plummer lens with the given parameters.

# Keyword Arguments
- `D_d::Real = NaN`: Angular diameter distance to the lens (in ``\\rm \\mathbf{arcseconds}``).
- `x_c = Vector{<:Real}()`: Vector of x-coordinates (in ``\\rm \\mathbf{arcseconds}``).
- `y_c = Vector{<:Real}()`: Vector of y-coordinates (in ``\\rm \\mathbf{arcseconds}``).
- `mass= Vector{<:Real}()`: Vector of masses (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s = Vector{<:Real}()`: Vector of scale radii (in ``\\rm \\mathbf{arcseconds}``).
"""
struct init_MultiPlummerLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   n::Int64
   x_c::Vector{T}
   y_c::Vector{T}
   mass::Vector{T}
   x_s::Vector{T}
end
function init_MultiPlummerLens(; D_d::Real=NaN, x_c=Vector{<:Real}(), y_c=Vector{<:Real}(), mass=Vector{<:Real}(), x_s=Vector{<:Real}())
   if !(length(x_c) == length(y_c) == length(mass) == length(x_s))
      throw(ArgumentError("x_c, y_c, mass, and x_s must all have the same length (one entry per component); 
            got $(length(x_c)), $(length(y_c)), $(length(mass)), $(length(x_s))."))
   end
   
   T = promote_type(typeof(D_d), eltype(x_c), eltype(y_c), eltype(mass), eltype(x_s))
   D_d  = convert(T, D_d)
   x_c  = Vector{T}(x_c)
   y_c  = Vector{T}(y_c)
   mass = Vector{T}(mass)
   x_s  = Vector{T}(x_s)
   
   n = length(x_c)
   return init_MultiPlummerLens{T}(:MultiPlummerLens, D_d, n, x_c, y_c, mass, x_s)
end


"""
    init_MultiGaussianLens(D_d::Real = NaN, 
                           x_c  = Vector{<:Real}, 
                           y_c  = Vector{<:Real}, 
                           mass = Vector{<:Real}, 
                           x_s  = Vector{<:Real})
Initialize a Multi-component Gaussian lens with the given parameters.

# Keyword Arguments
- `D_d = NaN`: Angular diameter distance to the lens (in ``\\rm \\mathbf{arcseconds}``).
- `x_c = Vector{<:Real}()`: Vector of x-coordinates (in ``\\rm \\mathbf{arcseconds}``).
- `y_c = Vector{<:Real}()`: Vector of y-coordinates (in ``\\rm \\mathbf{arcseconds}``).
- `mass= Vector{<:Real}()`: Vector of masses (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s = Vector{<:Real}()`: Vector of scale radii (in ``\\rm \\mathbf{arcseconds}``).
"""
struct init_MultiGaussianLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   D_d::T
   n::Int64
   x_c::Vector{T}
   y_c::Vector{T}
   mass::Vector{T}
   x_s::Vector{T}
end
function init_MultiGaussianLens(; D_d::Real=NaN, x_c=Vector{<:Real}(), y_c=Vector{<:Real}(), mass=Vector{<:Real}(), x_s=Vector{<:Real}())
   if !(length(x_c) == length(y_c) == length(mass) == length(x_s))
      throw(ArgumentError("x_c, y_c, mass, and x_s must all have the same length (one entry per component); 
            got $(length(x_c)), $(length(y_c)), $(length(mass)), $(length(x_s))."))
   end

   T = promote_type(typeof(D_d), eltype(x_c), eltype(y_c), eltype(mass), eltype(x_s))
   D_d  = convert(T, D_d)
   x_c  = Vector{T}(x_c)
   y_c  = Vector{T}(y_c)
   mass = Vector{T}(mass)
   x_s  = Vector{T}(x_s)
   
   n = length(x_c)
   return init_MultiGaussianLens{T}(:MultiGaussianLens, D_d, n, x_c, y_c, mass, x_s)
end


"""
    init_MultiPixelLens(x_c = Vector{<:Real}(), y_c = Vector{<:Real}(), kappa = Vector{<:Real}())
Initialize a Multi-component Pixel lens with the given parameters.

# Keyword Arguments
- `x_c = Vector{<:Real}()`: Vector of x-coordinates (in ``\\rm \\mathbf{arcseconds}``).
- `y_c = Vector{<:Real}()`: Vector of y-coordinates (in ``\\rm \\mathbf{arcseconds}``).
- `kappa = Vector{<:Real}()`: Vector of convergences (dimensionless).
"""
struct init_MultiPixelLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   n::Int64
   x_c::Vector{T}
   y_c::Vector{T}
   kappa::Vector{T}
   pixel_size::Vector{T}
end
function init_MultiPixelLens(; x_c=Vector{<:Real}(), y_c=Vector{<:Real}(), kappa=Vector{<:Real}(), pixel_size=Vector{<:Real}())
   if !(length(x_c) == length(y_c) == length(kappa) == length(pixel_size))
      throw(ArgumentError("x_c, y_c, kappa, pixel_size must all have the same length (one entry 
      per component); got $(length(x_c)), $(length(y_c)), $(length(kappa)) and $(length(pixel_size))."))
   end

   T = promote_type(eltype(x_c), eltype(y_c), eltype(kappa))
   x_c   = Vector{T}(x_c)
   y_c   = Vector{T}(y_c)
   kappa = Vector{T}(kappa)
   pixel_size = Vector{T}(pixel_size)
   
   n = length(x_c)
   return init_MultiPixelLens{T}(:MultiPixelLens, n, x_c, y_c, kappa, pixel_size)
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
- `pa = Vector{<:Real}()`: Vector of position angles (in ``\\rm \\mathbf{deg}``).
"""
struct init_MultiPJELens{T<:Real} <: AbstractLens
   _lens_::Symbol
   n::Int64
   x_c::Vector{T}
   y_c::Vector{T}
   v_d::Vector{T}
   x_s::Vector{T}
   x_t::Vector{T}
   eps::Vector{T}
   pa::Vector{T}
end
function init_MultiPJELens(; x_c = Vector{<:Real}(), y_c = Vector{<:Real}(), v_d = Vector{<:Real}(), x_s = Vector{<:Real}(), x_t = Vector{<:Real}(), eps = Vector{<:Real}(), pa = Vector{<:Real}())
   if !(length(x_c) == length(y_c) == length(v_d) == length(x_s) == length(x_t) == length(eps) == length(pa))
      throw(ArgumentError("x_c, y_c, v_d, x_s, x_t, eps, and pa must all have the same length 
                           (one entry per component); got $(length(x_c)), $(length(y_c)), 
                           $(length(v_d)), $(length(x_s)), $(length(x_t)), $(length(eps)), 
                           and $(length(pa))."))
   end

   T = promote_type(eltype(x_c), eltype(y_c), eltype(v_d), eltype(x_s), eltype(x_t), eltype(eps), eltype(pa))
   x_c  = Vector{T}(x_c)
   y_c  = Vector{T}(y_c)
   v_d = Vector{T}(v_d)
   x_s  = Vector{T}(x_s)
   x_t  = Vector{T}(x_t)
   eps  = Vector{T}(eps)
   pa   = Vector{T}(pa)

   n = length(x_c)
   return init_MultiPJELens(:MultiPJELens, n, x_c, y_c, v_d, x_s, x_t, eps, pa)
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
   :PixelLens         => (comp -> init_PixelLens(x_c=comp.x_c, y_c=comp.y_c, kappa=comp.kappa, pixel_size=comp.pixel_size)),
   :ExternalEffects   => (comp -> init_ExternalEffects(kappa=comp.kappa, gamma=comp.gamma, angle=comp.angle)),
   :ExternalEffects3  => (comp -> init_ExternalEffects3(delta=comp.delta, angle=comp.angle)),
   :Multipole         => (comp -> init_Multipole(delta=comp.delta, angle=comp.angle, m=comp.m, n=comp.n)),
   :PIEPLens          => (comp -> init_PIEPLens(x_c=comp.x_c, y_c=comp.y_c, v_d=comp.v_d, x_s=comp.x_s, eps=comp.eps, pa=comp.pa)),
   :SIELens           => (comp -> init_SIELens(x_c=comp.x_c, y_c=comp.y_c, v_d=comp.v_d, x_s=comp.x_s, eps=comp.eps, pa=comp.pa)),
   :PJELens           => (comp -> init_PJELens(x_c=comp.x_c, y_c=comp.y_c, v_d=comp.v_d, x_s=comp.x_s, x_t=comp.x_t, eps=comp.eps, pa=comp.pa)),
   :HernquistLens     => (comp -> init_HernquistLens(D_d=comp.D_d, x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, x_s=comp.x_s)),
   :eHernquistMDLens  => (comp -> init_eHernquistMDLens(D_d=comp.D_d, x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, x_s=comp.x_s, eps=comp.eps, pa=comp.pa)),
   :aHernquistLens    => (comp -> init_aHernquistLens(D_d=comp.D_d, x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, x_s=comp.x_s, eps=comp.eps, pa=comp.pa)),
   :NFWLens           => (comp -> init_NFWLens(comp.cosmology, comp.z_d; x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, c=comp.c)),
   :aNFWLens          => (comp -> init_aNFWLens(comp.cosmology, comp.z_d; x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, c=comp.c, eps=comp.eps, pa=comp.pa)),
   :eNFWMDLens        => (comp -> init_eNFWMDLens(comp.cosmology, comp.z_d; x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, c=comp.c, eps=comp.eps, pa=comp.pa)),
   :tNFWLens          => (comp -> init_tNFWLens(comp.cosmology, comp.z_d; x_c=comp.x_c, y_c=comp.y_c, x_s=comp.x_s, x_t=comp.x_t)),
   :gNFWLens          => (comp -> init_gNFWLens(comp.cosmology, comp.z_d; x_c=comp.x_c, y_c=comp.y_c, x_s=comp.x_s, n=comp.n)),
   :EinastoLens       => (comp -> init_EinastoLens(comp.cosmology, comp.z_d; x_c=comp.x_c, y_c=comp.y_c, x_s=comp.x_s, n=comp.n)),
   :MultiPlummerLens  => (comp -> init_MultiPlummerLens(D_d=comp.D_d, x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, x_s=comp.x_s)),
   :MultiGaussianLens => (comp -> init_MultiGaussianLens(D_d=comp.D_d, x_c=comp.x_c, y_c=comp.y_c, mass=comp.mass, x_s=comp.x_s)),
   :MultiPixelLens    => (comp -> init_PixelLens(x_c=comp.x_c, y_c=comp.y_c, kappa=comp.kappa, pixel_size=comp.pixel_size)),
   :MultiPJELens      => (comp -> init_MultiPJELens(x_c=comp.x_c, y_c=comp.y_c, v_d=comp.v_d, x_s=comp.x_s, x_t=comp.x_t, eps=comp.eps, pa=comp.pa))
   )

"""
    init_CompositeLens(lens::Vector{<:NamedTuple})
Initialize a composite lens from a vector of lens components.

# Arguments
- `lens::Vector{<:NamedTuple}`: Vector of lens components.

# Returns
- `CompositeLens`: Composite lens.
"""
struct init_CompositeLens <: AbstractLens
   _lens_::Symbol
   _components_::Vector{AbstractLens}
end
function init_CompositeLens(; _components_::Vector{<:AbstractLens}=Vector{AbstractLens}())
   return init_CompositeLens(:CompositeLens, _components_)
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


"""
    init_MultiPlaneLens(lens::Vector{<:NamedTuple})
Initialize a multi-plane lens from a vector of lens components. Each component must contain the
lens redshift (`z_d`) along with the parameters of the corresponding lens model. Components
sharing the same redshift are grouped into a single (composite) lens plane, and the lens planes
are sorted in increasing redshift. At least two distinct lens planes are required.

# Arguments
- `lens::Vector{<:NamedTuple}`: Vector of lens components.

# Returns
- `MultiPlaneLens`: Multi-plane lens.
"""
struct init_MultiPlaneLens{T<:Real} <: AbstractLens
   _lens_::Symbol
   n_p::Int64
   z_d::Vector{T}
   _plane_::Vector{AbstractLens}
end
function init_MultiPlaneLens(; n_p::Int64=0, z_d::Vector{<:Real}=Float64[], _plane_::Vector{<:AbstractLens}=AbstractLens[])
   if n_p != length(z_d)
      throw(ArgumentError("n_p ($n_p) must match length(z_d) ($(length(z_d)))."))
   end
   
   T = eltype(z_d)
   return init_MultiPlaneLens{T}(:MultiPlaneLens, n_p, z_d, Vector{AbstractLens}(_plane_))
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