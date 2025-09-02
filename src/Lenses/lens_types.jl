export init_PointLens
export init_PlummerLens
export init_SISLens

abstract type AbstractLens end


"""
    init_PointLens(D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN)

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
    init_PlummerLens(D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN, x_s::Real=NaN)

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
    init_SISLens(x_c::Real=0.0, y_c::Real=0.0, v_d::Real=NaN)

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
    init_NSISPLens(x_c::Real=0.0, y_c::Real=0.0, v_d::Real=NaN, x_s::Real=NaN)

Initialize a Non-Singular Isothermal Sphere potential (NSISP) lens with the given parameters.
"""
@kwdef struct init_NSISPLens <: AbstractLens
   _lens_::Symbol = :NSISPLens
   _lid_::Int64 = 4
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   v_d::RV = NaN
   x_s::RV = NaN
end


"""
    init_NSISMDLens(x_c::Real=0.0, y_c::Real=0.0, v_d::Real=NaN, x_s::Real=NaN)

Initialize a Non-Singular Isothermal Sphere mass distribution (NSISMD) lens with the given parameters.
"""
@kwdef struct init_NSISMDLens <: AbstractLens
   _lens_::Symbol = :NSISMDLens
   _lid_::Int64 = 5
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   v_d::RV = NaN
   x_s::RV = NaN
end


@kwdef struct init_EinastoLens <: AbstractLens
   _lens_::Symbol = :EinastoLens
   _lid_::Int64 = 18
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   n_a::RV =-2.0
   x_s::RV = NaN
   c::RV   = NaN
   rho_s::RV = NaN
   mass::RV  = NaN
end


@kwdef struct init_NFWLens <: AbstractLens
   _lens_::Symbol = :NFWLens
   _lid_::Int64 = 11
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   x_s::RV = NaN
   c::RV   = NaN
   rho_s::RV = NaN
   mass::RV  = NaN
end


@kwdef struct init_gNFWLens <: AbstractLens
   _lens_::Symbol = :gNFWLens
   _lid_::Int64 = 14
   D_d::RV = NaN
   x_c::RV = 0.0
   y_c::RV = 0.0
   n_a::RV = NaN
   x_s::RV = NaN
   c::RV   = NaN
   rho_s::RV = NaN
   mass::RV  = NaN
end


@kwdef struct init_CompositeLens <: AbstractLens
   _lens_::Symbol = :CompositeLens
   _lid_::Int64 = 111
   _components_ = Vector{AbstractLens}()
end


@kwdef struct init_MultiPlaneLens <: AbstractLens
   _lens_::Symbol = :MultiPlaneLens
   _lid_::Int64 = 222
   _np_::Int64  = NaN
   _zl_ = Vector{Real}()
   _components_ = Vector{AbstractLens}()
end


# Constructor for NFW lens
function init_NFWLens(cosmology::Cosmology.AbstractCosmology, z_d::RV;
                     x_c::RV=0.0, y_c::RV=0.0, x_s::RV=NaN, c::RV=NaN, mass::RV=NaN)                     
   # ADD to the lens
   D_d = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in meters)
   r_vir = (3.0 * mass / 72.0 / pi^3 / ρ_cz)^(1.0/3.0)

   # Check if concentration is given
   if isfinite(c)
      ρ_s = (18.0 * pi^2 / 3.0) * ρ_cz * c^3 / ( log(1.0 + c) - ( c / (1.0 + c) ) )
      x_s = r_vir / c / D_d
   elseif isfinite(x_s)
      c = r_vir / (x_s * D_d)
      ρ_s = (18.0 * pi^2 / 3.0) * ρ_cz * c^3 / ( log(1.0 + c) - ( c / (1.0 + c) ) )
   else
      throw(ArgumentError("Provide at least c or x_s in ** init_NFWLens **"))
   end
   return init_NFWLens(D_d=D_d, x_c=x_c, y_c=y_c, x_s=x_s, c=c, rho_s=ρ_s, mass=mass)
end


# Constructor for NFW lens
function init_gNFWLens(cosmology::Cosmology.AbstractCosmology, z_d::RV;
                     x_c::RV=0.0, y_c::RV=0.0, x_s::RV=NaN, c::RV=NaN, mass::RV=NaN, n_a::RV=-2.0)                     
   # Calculate the mass
   function intgrand(x::RV)
      return x^(2.0 - n_a) / (1.0 + x)^(3.0 - n_a)
   end

   # ADD to the lens
   D_d = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in meters)
   r_vir= (3.0 * mass / 72.0 / pi^3 / ρ_cz)^(1.0/3.0)

   # Check if concentration is given
   if isfinite(c)
      x_s = r_vir / c / D_d
      mass, _ = quadgk(x -> integrand(x), 0, c)

      ρ_s = (18.0 * pi^2 / 3.0) * ρ_cz * c^3 / mass 
   elseif isfinite(x_s)
      c = r_vir / (x_s * D_d)
      mass, _ = quadgk(x -> integrand(x), 0, c)

      ρ_s = (18.0 * pi^2 / 3.0) * ρ_cz * c^3 / mass
   else
      throw(ArgumentError("Provide at least c or x_s in ** init_gNFWLens **"))
   end
   return init_gNFWLens(D_d=D_d, x_c=x_c, y_c=y_c, x_s=x_s, c=c, rho_s=ρ_s, mass=mass, n_a=n_a)
end


# Constructor for Einasto lens
function init_EinastoLens(cosmology::Cosmology.AbstractCosmology, z_d::RV;
                     x_c::RV=0.0, y_c::RV=0.0, x_s::RV=NaN, c::RV=NaN, mass::RV=NaN, n_a::RV = -2.0)
   # ADD to the lens
   D_d  = Cosmology.angular_diameter_distance(cosmology, 0., z_d)

   # Critical density at the lens redshift
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in meters)
   r_vir= (3.0 * mass / 72.0 / pi^3 / ρ_cz)^(1.0/3.0)

   # Check if concentration is given
   m_v::Float64 = 0.0
   if isfinite(c)
      Pax, _ = gamma_inc( 3.0/n_a, (2.0/n_a)*c^n_a )
      m_v = (1.0/n_a) * (n_a/2.0)^(3.0/n_a) * gamma(3.0/n_a) * Pax

      ρ_s = (18.0 * pi^2 / 3.0) * ρ_cz * c^3 / m_v
   elseif isfinite(x_s)
      c = r_vir / (x_s * D_d)
      Pax, _ = gamma_inc( 3.0/n_a, (2.0/n_a)*c^n_a )
      m_v = (1.0/n_a) * (n_a/2.0)^(3.0/n_a) * gamma(3.0/n_a) * Pax

      ρ_s = (18.0 * pi^2 / 3.0) * ρ_cz * c^3 / m_v
   else
      throw(ArgumentError("Provide at least c or x_s in ** init_EinastoLens **"))
   end
   return init_EinastoLens(D_d=D_d, x_c=x_c, y_c=y_c, n_a=n_a, x_s=x_s, c=c, rho_s=ρ_s, mass=mass)
end


# Dictionary to map lens types to their initialization functions and arguments
const lens_init_functions = Dict(
   :PointLens       => (init_PointLens,        [:D_d, :x_c, :y_c, :mass]),
   :PlummerLens     => (init_PlummerLens,      [:D_d, :x_c, :y_c, :mass, :x_s]),
   :SISLens         => (init_SISLens,          [:x_c, :y_c, :v_d]),
   :NSISPLens       => (init_NSISPLens,        [:D_d, :x_c, :y_c, :v_d, :θ_s]),
   )
# Constructor for composite lens
function init_CompositeLens(lens::Vector{<:NamedTuple})
   # Define an component vector
   lens_components = AbstractLens[]

   # Run over the lens components in the composite lens
   for component in lens
      # Check and get the lens component if available otherwise throw an error
      val = get(lens_init_functions, component.lens, nothing)
      if val === nothing
         throw(ArgumentError("Unknown lens type: $(component.lens). Available lens models are $(keys(lens_init_functions))"))
      end

      # Get the compoent init_* function and arguments
      init_func, init_args = val

      # Create a dict of args to pass to the init_* function
      kwargs = Dict(arg => getproperty(component, arg) for arg in init_args)

      # Push the component into the vector
      push!(lens_components, init_func(; kwargs...))
   end
   return init_CompositeLens(_components_=lens_components)
end


# Constructor for multi-plane lens
function init_MultiPlaneLens(lens::Vector{<:NamedTuple})   
   # Get sorted unique lens redshifts
   zl_unique = unique(component.zl for component in lens)
   sort!(zl_unique)

   # Group lens components by redshift
   lens_by_z = Dict{eltype(getfield.(lens, :zl)), Vector{NamedTuple}}()
   for component in lens
      push!(get!(lens_by_z, component.zl, []), component)
   end

   # Construct composite lenses for each unique redshift
   lens_components = AbstractLens[init_CompositeLens(lens_by_z[z]) for z in zl_unique]

   return init_MultiPlaneLens(_np_=length(zl_unique), _zl_=zl_unique, _components_=lens_components)
end