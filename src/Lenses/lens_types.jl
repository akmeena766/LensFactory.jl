export init_PointLens
export init_PlummerLens
export init_SISLens
export init_NSISPLens
export init_NSISMDLens
export init_GaussianLens
export init_SersicLens
export init_ExternalEffects
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
    init_PointLens(D_d::RV=NaN, x_c::RV=0.0, y_c::RV=0.0, mass::RV=NaN)

Initialize a point lens with the given parameters.
"""
@kwdef struct init_PointLens <: AbstractLens
   _lens_::Symbol = :PointLens
   _lid_::Int64 = 1
   D_d::RV  = NaN
   x_c::RV  = 0.0
   y_c::RV  = 0.0
   mass::RV = NaN
end


"""
    init_PlummerLens(D_d::RV=NaN, x_c::RV=0.0, y_c::RV=0.0, mass::RV=NaN, x_s::RV=NaN)

Initialize a Plummer lens with the given parameters.
"""
@kwdef struct init_PlummerLens <: AbstractLens
   _lens_::Symbol = :PlummerLens
   _lid_::Int64 = 2
   D_d::RV  = NaN
   x_c::RV  = 0.0
   y_c::RV  = 0.0
   mass::RV = NaN
   x_s::RV  = NaN
end


"""
    init_SISLens(x_c::RV=0.0, y_c::RV=0.0, v_d::RV=NaN)

Initialize a Singular Isothermal Sphere (SIS) lens with the given parameters.
"""
@kwdef struct init_SISLens <: AbstractLens
   _lens_::Symbol = :SISLens
   _lid_::Int64 = 3
   x_c::RV = 0.0
   y_c::RV = 0.0
   v_d::RV = NaN
end


"""
    init_NSISPLens(x_c::RV=0.0, y_c::RV=0.0, v_d::RV=NaN, x_s::RV=NaN)

Initialize a Non-Singular Isothermal Sphere potential (NSISP) lens with the given parameters.
"""
@kwdef struct init_NSISPLens <: AbstractLens
   _lens_::Symbol = :NSISPLens
   _lid_::Int64 = 4
   x_c::RV = 0.0
   y_c::RV = 0.0
   v_d::RV = NaN
   x_s::RV = NaN
end


"""
    init_NSISMDLens(x_c::RV=0.0, y_c::RV=0.0, v_d::RV=NaN, x_s::RV=NaN)

Initialize a Non-Singular Isothermal Sphere mass distribution (NSISMD) lens with the given parameters.
"""
@kwdef struct init_NSISMDLens <: AbstractLens
   _lens_::Symbol = :NSISMDLens
   _lid_::Int64 = 5
   x_c::RV = 0.0
   y_c::RV = 0.0
   v_d::RV = NaN
   x_s::RV = NaN
end


"""
    init_GaussianLens(D_d::RV=NaN, x_c::RV=0.0, y_c::RV=0.0, mass::RV=NaN, x_s::RV=NaN)
Initialize a Gaussian lens with the given parameters.
"""
@kwdef struct init_GaussianLens <: AbstractLens
   _lens_::Symbol = :GaussianLens
   _lid_::Int = 6
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   mass::RV= NaN
   x_s::RV = NaN
end


"""
    init_SersicLens(D_d::RV=NaN, x_c::RV=0.0, y_c::RV=0.0, mass::RV=NaN, x_e::RV=NaN, n::RV=4)
Initialize a Sersic lens with the given parameters.
"""
@kwdef struct init_SersicLens <: AbstractLens
   _lens_::Symbol = :SersicLens
   _lid_::Int = 7
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   mass::RV= NaN
   x_e::RV = NaN
   n::RV   = 4.0
end


"""
    init_ExternalEffects(kappa::RV=NaN, gamma::RV=NaN, angle::RV=NaN)
Initialize constant external effects with the given parameters.
"""
@kwdef struct init_ExternalEffects <: AbstractLens
   _lens_::Symbol = :ExternalEffects
   _lid_::Int = 8
   kappa::RV = NaN
   gamma1::RV = NaN
   gamma2::RV = NaN
end


"""
    init_PIEPLens(x_c::RV=0.0, y_c::RV=0.0, v_d::RV=NaN, x_s::RV=NaN, eps::RV=NaN, pa::RV=NaN)
Initialize pseudo isothermal elliptical potential (PIEP) lens with the given parameters.
"""
@kwdef struct init_PIEPLens <: AbstractLens
   _lens_::Symbol = :PIEPLens
   _lid_::Int = 9
   x_c::RV = 0.0
   y_c::RV = 0.0
   v_d::RV = NaN
   x_s::RV = NaN
   eps::RV = NaN
   pa::RV  = NaN
end


"""
    init_SIELens(x_c::RV=0.0, y_c::RV=0.0, v_d::RV=NaN, x_s::RV=NaN, eps::RV=NaN, pa::RV=NaN)
Initialize singular isothermal ellipsoid (SIE) lens with the given parameters.
"""
@kwdef struct init_SIELens <: AbstractLens
   _lens_::Symbol = :SIELens
   _lid_::Int = 10
   x_c::RV = 0.0
   y_c::RV = 0.0
   v_d::RV = NaN
   x_s::RV = NaN
   eps::RV = NaN
   pa::RV  = NaN
end


"""
    init_PJEMDLens(x_c::RV=0.0, y_c::RV=0.0, v_d::RV=NaN, x_s::RV=NaN, x_t::RV=NaN, eps::RV=NaN, pa::RV=NaN)
Initialize Pseudo-Jaffe Ellipsoid (PJE) lens with the given parameters.
"""
@kwdef struct init_PJELens <: AbstractLens
   _lens_::Symbol = :PJELens
   _lid_::Int = 11
   x_c::RV = 0.0
   y_c::RV = 0.0
   v_d::RV = NaN
   x_s::RV = NaN
   x_t::RV = NaN
   eps::RV = NaN
   pa::RV  = NaN
end


"""
    init_HernquistLens(D_d::RV=NaN, x_c::RV=0.0, y_c::RV=0.0, mass::RV=NaN, x_s::RV=NaN)
Initialize a Hernquist lens with the given parameters.
"""
@kwdef struct init_HernquistLens <: AbstractLens
   _lens_::Symbol = :HernquistLens
   _lid_::Int64 = 12
   D_d::RV  = NaN
   x_c::RV  = 0.0
   y_c::RV  = 0.0
   mass::RV = NaN
   x_s::RV  = NaN
end

"""
    init_NFWLens(cosmology::Cosmology.AbstractCosmology, z_d::RV; x_c::RV=0.0, y_c::RV=0.0, x_s::RV=NaN, c::RV=NaN, mass::RV=NaN)
Initialize a Navarro-Frenk-White (NFW) lens with the given parameters
"""
@kwdef struct init_NFWLens <: AbstractLens
   _lens_::Symbol = :NFWLens
   _lid_::Int64 = 13
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   rho_s::RV = NaN
   x_s::RV = NaN
end


"""
    init_tNFWLens(cosmology::Cosmology.AbstractCosmology, z_d::RV; x_c::RV=0.0, y_c::RV=0.0, x_s::RV=NaN, x_t::RV=NaN, c::RV=NaN, mass::RV=NaN)
Initialize a truncated Navarro-Frenk-White (tNFW) lens with the given parameters.
"""
@kwdef struct init_tNFWLens <: AbstractLens
   _lens_::Symbol = :tNFWLens
   _lid_::Int64 = 14
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   rho_s::RV = NaN
   x_s::RV = NaN
   x_t::RV = NaN
end


"""
    init_gNFWLens(D_d::RV=NaN, x_c::RV=0.0, y_c::RV=0.0, x_s::RV=NaN, c::RV=NaN, rho_s::RV=NaN, mass::RV=NaN, n::RV=-2.0)
Initialize a generalized Navarro-Frenk-White (gNFW) lens with the given parameters.
"""
@kwdef struct init_gNFWLens <: AbstractLens
   _lens_::Symbol = :gNFWLens
   _lid_::Int64 = 15
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   rho_s::RV = NaN
   x_s::RV = NaN
   n::RV   = NaN
end


"""
    init_EinastoLens(D_d::RV=NaN, x_c::RV=0.0, y_c::RV=0.0, x_s::RV=NaN, c::RV=NaN, rho_s::RV=NaN, mass::RV=NaN, n::RV=-2.0)
Initialize an Einasto lens with the given parameters. The lens model can be initialized with either
the concentration `c` or the scale radius `x_s`. 
**If both are provided, `c` will be used to calculate `x_s` and the input `x_s` will be overwritten.**
The parameter `n` defines the slope of the density profile.
"""
@kwdef struct init_EinastoLens <: AbstractLens
   _lens_::Symbol = :EinastoLens
   _lid_::Int64 = 16
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   rho_s::RV = NaN
   x_s::RV = NaN
   n::RV   =-2.0
end


"""
    init_aHernquistLens(D_d::RV=NaN, x_c::RV=0.0, y_c::RV=0.0, mass::RV=NaN, x_s::RV=NaN, eps::RV=NaN, pa::RV=NaN)
Initialize an approximate Hernquist lens (aHernquistLens) with the given parameters based on 
[Oguri (2021)](https://ui.adsabs.harvard.edu/abs/2021PASP..133g4504O/abstract).
"""
@kwdef struct init_aHernquistLens <: AbstractLens
   _lens_::Symbol = :aHernquistLens
   _lid_::Int64 = 17
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   mass::RV = NaN
   x_s::RV  = NaN
   eps::RV  = NaN
   pa::RV   = NaN
end

"""
    init_aNFWLens(D_d::RV=NaN, x_c::RV=0.0, y_c::RV=0.0, rho_s::RV=NaN, x_s::RV=NaN, eps::RV=NaN, pa::RV=NaN)
Initialize an approximate Navarro-Frenk-White lens (aNFWLens) with the given parameters based on 
[Oguri (2021)](https://ui.adsabs.harvard.edu/abs/2021PASP..133g4504O/abstract).
"""
@kwdef struct init_aNFWLens <: AbstractLens
   _lens_::Symbol = :aNFWLens
   _lid_::Int64 = 17
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   rho_s::RV = NaN
   x_s::RV  = NaN
   eps::RV  = NaN
   pa::RV   = NaN
end

"""
    init_eHernquistMDLens(D_d::RV=NaN, x_c::RV=0.0, y_c::RV=0.0, mass::RV=NaN, x_s::RV=NaN, eps::RV=NaN, pa::RV=NaN)
Initialize an elliptical Hernquist mass distribution lens (eHernquistMDLens) with the given parameters.
"""
@kwdef struct init_eHernquistMDLens <: AbstractLens
   _lens_::Symbol = :eHernquistMDLens
   _lid_::Int64 = 18
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   mass::RV = NaN
   x_s::RV  = NaN
   eps::RV  = NaN
   pa::RV   = NaN
end

"""
    init_eNFWMDLens(D_d::RV=NaN, x_c::RV=0.0, y_c::RV=0.0, rho_s::RV=NaN, x_s::RV=NaN, eps::RV=NaN, pa::RV=NaN)
Initialize an elliptical Navarro-Frenk-White mass distribution lens (eNFWMDLens) with the given parameters.
"""
@kwdef struct init_eNFWMDLens <: AbstractLens
   _lens_::Symbol = :eNFWMDLens
   _lid_::Int64 = 19
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   rho_s::RV = NaN
   x_s::RV  = NaN
   eps::RV  = NaN
   pa::RV   = NaN
end

@kwdef struct init_MultiPlummerLens <: AbstractLens
   _lens_::Symbol = :MultiPlummerLens
   _lid_::Int64 = 101
   D_d::RV  = NaN
   n::Int64 = NaN
   x_c  = Vector{<:RV}()
   y_c  = Vector{<:RV}() 
   mass = Vector{<:RV}()
   x_s  = Vector{<:RV}()
end


@kwdef struct init_MultiGaussianLens <: AbstractLens
   _lens_::Symbol = :MultiGaussianLens
   _lid_::Int64 = 102
   D_d::RV  = NaN
   n::Int64 = NaN
   x_c  = Vector{Float64}()
   y_c  = Vector{Float64}() 
   mass = Vector{Float64}()
   x_s  = Vector{Float64}()
end


@kwdef struct init_MultiPJELens <: AbstractLens
   _lens_::Symbol = :MultiPJELens
   _lid_::Int64 = 103
   n::Int64 = NaN
   x_c  = Vector{Float64}()
   y_c  = Vector{Float64}() 
   v_d  = Vector{Float64}()
   x_s  = Vector{Float64}()
   x_t  = Vector{Float64}()
   eps  = Vector{Float64}()
   pa   = Vector{Float64}()
end


@kwdef struct init_CompositeLens <: AbstractLens
   _lens_::Symbol = :CompositeLens
   _lid_::Int64 = 111
   _components_ = Vector{AbstractLens}()
end


@kwdef struct init_MultiPlaneLens <: AbstractLens
   _lens_::Symbol = :MultiPlaneLens
   _lid_::Int64 = 222
   n_p::Int64   = NaN
   z_d = Vector{<:RV}()
   _plane_ = Vector{AbstractLens}()
end


#---------------------- Composite and Multi-plane lens constructors --------------------------------
# Dictionary to map lens types to their initialization functions and arguments
const lens_init_functions = Dict{Symbol, Function}(
   :PointLens         => (c -> init_PointLens(D_d=c.D_d, x_c=c.x_c, y_c=c.y_c, mass=c.mass)),
   :PlummerLens       => (c -> init_PlummerLens(D_d=c.D_d, x_c=c.x_c, y_c=c.y_c, mass=c.mass, x_s=c.x_s)),
   :SISLens           => (c -> init_SISLens(x_c=c.x_c, y_c=c.y_c, v_d=c.v_d)),
   :NSISPLens         => (c -> init_NSISPLens(x_c=c.x_c, y_c=c.y_c, v_d=c.v_d, x_s=c.x_s)),
   :NSISMDLens        => (c -> init_NSISMDLens(x_c=c.x_c, y_c=c.y_c, v_d=c.v_d, x_s=c.x_s)),
   :GaussianLens      => (c -> init_GaussianLens(D_d=c.D_d, x_c=c.x_c, y_c=c.y_c, mass=c.mass, x_s=c.x_s)),
   :SersicLens        => (c -> init_SersicLens(D_d=c.D_d, x_c=c.x_c, y_c=c.y_c, mass=c.mass, x_e=c.x_e, n=c.n)),
   :ExternalEffects   => (c -> init_ExternalEffects(kappa=c.kappa, gamma1=c.gamma1, gamma2=c.gamma2)),
   :PIEPLens          => (c -> init_PIEPLens(x_c=c.x_c, y_c=c.y_c, v_d=c.v_d, x_s=c.x_s, eps=c.eps, pa=c.pa)),
   :SIELens           => (c -> init_SIELens(x_c=c.x_c, y_c=c.y_c, v_d=c.v_d, x_s=c.x_s, eps=c.eps, pa=c.pa)),
   :PJELens           => (c -> init_PJELens(x_c=c.x_c, y_c=c.y_c, v_d=c.v_d, x_s=c.x_s, x_t=c.x_t, eps=c.eps, pa=c.pa)),
   :HernquistLens     => (c -> init_HernquistLens(D_d=c.D_d, x_c=c.x_c, y_c=c.y_c, mass=c.mass, x_s=c.x_s)),
   :NFWLens           => (c -> init_NFWLens(D_d=c.D_d, x_c=c.x_c, y_c=c.y_c, x_s=c.x_s, rho_s=c.rho_s)),
   :tNFWLens          => (c -> init_tNFWLens(D_d=c.D_d, x_c=c.x_c, y_c=c.y_c, x_s=c.x_s, x_t=c.x_t, rho_s=c.rho_s)),
   :gNFWLens          => (c -> init_gNFWLens(D_d=c.D_d, x_c=c.x_c, y_c=c.y_c, x_s=c.x_s, rho_s=c.rho_s, n=c.n)),
   :EinastoLens       => (c -> init_EinastoLens(D_d=c.D_d, x_c=c.x_c, y_c=c.y_c, x_s=c.x_s, rho_s=c.rho_s, n=c.n)),
   :MultiPlummerLens  => (c -> init_MultiPlummerLens(D_d=c.D_d, n=c.n, x_c=c.x_c, y_c=c.y_c, mass=c.mass, x_s=c.x_s)),
   :MultiGaussianLens => (c -> init_MultiGaussianLens(D_d=c.D_d, n=c.n, x_c=c.x_c, y_c=c.y_c, mass=c.mass, x_s=c.x_s)),
   :MultiPJELens      => (c -> init_MultiPJELens(n=c.n, x_c=c.x_c, y_c=c.y_c, v_d=c.v_d, x_s=c.x_s, x_t=c.x_t, eps=c.eps, pa=c.pa))
   )
# Constructor for composite lens
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


# Constructor for multi-plane lens
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


#---------------------- Parameter functions for various lenses -------------------------------------
# Parameters for NFW lens
function parameter_NFWLens(; cosmology::Cosmology.AbstractCosmology=nothing, z_d::RV=NaN, mass::RV=NaN, x_s::RV=NaN, c::RV=NaN)
   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift (in kg/m^3)
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide concentration (c) or scale radius (x_s) in ** parameter_NFWLens **"))
   end
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / (log(1.0 + c) - (c / (1.0 + c)))

   return (mass=mass, rho_s=ρ_s, c=c, x_s=x_s)
end

# Parameters for tNFW lens
function parameter_tNFWLens(; cosmology::Cosmology.AbstractCosmology=nothing, z_d::RV=NaN, mass::RV=NaN, x_s::RV=NaN, c::RV=NaN, x_t::RV=NaN)
   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide at least c or x_s in ** parameter_tNFWLens **"))
   end
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / (log(1.0 + c) - (c / (1.0 + c)))

   return (mass=mass, rho_s=ρ_s, c=c, x_s=x_s, x_t=x_t)
end

# Parameters for generalized NFW lens
function parameter_gNFWLens(; cosmology::Cosmology.AbstractCosmology=nothing, z_d::RV, mass::RV=NaN, x_s::RV=NaN, c::RV=NaN, n::RV=1.0)
   # Check for valid slope parameter
   if !(0.0 < n < 2.0)
      throw(ArgumentError("Slope parameter outside allowed range n ∈ (0, 2) in ** parameter_gNFWLens **"))
   end

   # Integrand function for mass calculation
   function integrand(x::RV, α::RV)
      return x^(2.0 - α) / (1.0 + x)^(3.0 - α)
   end

   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide at least c or x_s in ** parameter_gNFWLens **"))
   end
   mass, _ = quadgk(x -> integrand(x, n), 0, c)
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / mass

   return (mass=mass, rho_s=ρ_s, c=c, x_s=x_s, n=n)
end

# Parameters for Einasto lens
function parameter_EinastoLens(; cosmology::Cosmology.AbstractCosmology=nothing, z_d::RV, mass::RV=NaN, x_s::RV=NaN, c::RV=NaN, n::RV=0.2)
   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d  = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)

   # Critical density at the lens redshift
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide at least c or x_s in ** parameter_EinastoLens **"))
   end
   Pax, _ = gamma_inc(3.0 / n, (2.0 / n) * c^n)
   m_v = (1.0 / n) * (n / 2.0)^(3.0 / n) * gamma(3.0 / n) * Pax
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / m_v

   return (mass=mass, rho_s=ρ_s, c=c, x_s=x_s, n=n)
end
