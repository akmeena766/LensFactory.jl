abstract type AbstractLens end


"""
    init_PointLens(D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, mass::Real=NaN)

Initialize a point lens with the given parameters.
"""
@kwdef struct init_PointLens <: AbstractLens
   _lens_::String = "PointLens"
   _lid_::Int = 1
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   mass::Real= NaN
end


"""
    init_SISLens(D_d::Real=NaN, x_c::Real=0.0, y_c::Real=0.0, v_d::Real=NaN)

Initialize a Singular Isothermal Sphere (SIS) lens with the given parameters.
"""
@kwdef struct init_SISLens <: AbstractLens
   _lens_::String = "SISLens"
   _lid_::Int = 2
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   v_d::Real = NaN
end


@kwdef struct init_ExternalEffects <: AbstractLens
   _lens_::String = "ExternalEffects"
   _lid_::Int = 1
   k_ext::Real  = NaN
   g1_ext::Real = NaN
   g2_ext::Real = NaN
end


@kwdef struct init_PlummerLens <: AbstractLens
   _lens_::String = "PlummerLens"
   _lid_::Int = 3
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   x_s::Real = NaN
   mass::Real= NaN
end





@kwdef struct init_NSISPLens <: AbstractLens
   _lens_::String = "NSISPLens"
   _lid_::Int = 5
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   x_s::Real = NaN
   v_d::Real = NaN
end


@kwdef struct init_NSISMDLens <: AbstractLens
   _lens_::String = "NSISMDLens"
   _lid_::Int = 6
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   x_s::Real = NaN
   v_d::Real = NaN
end


@kwdef struct init_GaussianLens <: AbstractLens
   _lens_::String = "GaussianLens"
   _lid_::Int = 20
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   mass::Real= NaN
   x_s::Real = NaN
end


@kwdef struct init_PowerLawlens <: AbstractLens
   _lens_::String = "PowerLawLens"
   _lid_::Int = 21
   x_c::Real = 0.0
   y_c::Real = 0.0
   x_s::Real = NaN
   n_a::Real = NaN
end


@kwdef struct init_SersicLens <: AbstractLens
   _lens_::String = "SersicLens"
   _lid_::Int = 22
   x_c::Real = 0.0
   y_c::Real = 0.0
   x_e::Real = NaN
   mass::Real= NaN
end


@kwdef struct init_HernquistLens <: AbstractLens
   _lens_::String = "HernquistLens"
   _lid_::Int = 10
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   x_s::Real = NaN
   mass::Real= NaN
end


@kwdef struct init_EinastoLens <: AbstractLens
   _lens_::String = "EinastoLens"
   _lid_::Int = 18
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   n_a::Real =-2.0
   x_s::Real = NaN
   c::Real   = NaN
   rho_s::Real=NaN
   mass::Real= NaN
end


@kwdef struct init_NFWLens <: AbstractLens
   _lens_::String = "NFWLens"
   _lid_::Int = 11
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   x_s::Real = NaN
   c::Real   = NaN
   rho_s::Real=NaN
   mass::Real= NaN
end


@kwdef struct init_gNFWLens <: AbstractLens
   _lens_::String = "gNFWLens"
   _lid_::Int = 14
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   x_s::Real = NaN
   c::Real   = NaN
   rho_s::Real=NaN
   mass::Real= NaN
   n_a::Real = NaN
end


@kwdef struct init_tNFWLens <: AbstractLens
   _lens_::String = "tNFWLens"
   _lid_::Int = 16
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
end


@kwdef struct init_PIEPLens <: AbstractLens
   _lens_::String = "PIEPLens"
   _lid_::Int = 7
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   x_s::Real = 0.0
   v_d::Real = NaN
   eps::Real = NaN
   p_a::Real = 0.0
   n_a::Real = 0.5 
end


@kwdef struct init_SIEMDLens <: AbstractLens
   _lens_::String = "SIEMDLens"
   _lid_::Int = 8
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   x_s::Real = 0.0
   v_d::Real = NaN
   eps::Real = NaN
   p_a::Real = 0.0
end


@kwdef struct init_PJEMDLens <: AbstractLens
   _lens_::String = "PJEMDLens"
   _lid_::Int = 9
   D_d::Real = NaN
   x_c::Real = 0.0
   y_c::Real = 0.0
   x_s::Real = NaN
   x_t::Real = NaN
   v_d::Real = NaN
   eps::Real = NaN
   p_a::Real = NaN
end


@kwdef struct init_CompositeLens <: AbstractLens
   _lens_::String = "CompositeLens"
   _lid_::Int = 111
   _components_ = Vector{AbstractLens}()
end


@kwdef struct init_MultiPlaneLens <: AbstractLens
   _lens_::String = "MultiPlaneLens"
   _lid_::Int = 222
   _np_::Int = NaN
   _zl_ = Vector{Real}()
   _components_ = Vector{AbstractLens}()
end


# Constructor for NFW lens
function init_NFWLens(cosmology::Cosmology.AbstractCosmology, z_d::Real;
                     x_c::Real=0.0, y_c::Real=0.0, x_s::Real=NaN, c::Real=NaN, mass::Real=NaN)                     
   # ADD to the lens
   D_d::Real = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift
   ρ_cz::Real = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in meters)
   r_vir::Real= (3.0 * mass / 72.0 / pi^3 / ρ_cz)^(1.0/3.0)

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
function init_gNFWLens(cosmology::Cosmology.AbstractCosmology, z_d::Real;
                     x_c::Real=0.0, y_c::Real=0.0, 
                     x_s::Real=NaN, c::Real=NaN, mass::Real=NaN, n_a::Real=-2.0)                     
   # Calculate the mass
   function intgrand(x::Real)::Real
      return x^(2.0 - n_a) / (1.0 + x)^(3.0 - n_a)
   end

   # ADD to the lens
   D_d::Real = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift
   ρ_cz::Real = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in meters)
   r_vir::Real= (3.0 * mass / 72.0 / pi^3 / ρ_cz)^(1.0/3.0)

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
function init_EinastoLens(cosmology::Cosmology.AbstractCosmology, z_d::Real;
                           x_c::Real=0.0, y_c::Real=0.0,
                           x_s::Real=NaN, c::Real=NaN, mass::Real=NaN, n_a::Real = -2.0)
   # ADD to the lens
   D_d::Real  = Cosmology.angular_diameter_distance(cosmology, 0., z_d)

   # Critical density at the lens redshift
   ρ_cz::Real = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in meters)
   r_vir::Real= (3.0 * mass / 72.0 / pi^3 / ρ_cz)^(1.0/3.0)

   # Check if concentration is given
   m_v::Real = 0.0
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
lens_init_functions = Dict(
   "ExternalEffects" => (init_ExternalEffects, [:k_ext, :g1_ext, :g2_ext]),
   "PointLens"       => (init_PointLens,       [:x_c, :y_c, :mass]),
   "PlummerLens"     => (init_PlummerLens,     [:x_c, :y_c, :x_s, :mass]),
   "SISLens"         => (init_SISLens,         [:x_c, :y_c, :v_d]),
   "NSISPLens"       => (init_NSISPLens,       [:x_c, :y_c, :x_s, :v_d]),
   "NSISMDLens"      => (init_NSISMDLens,      [:x_c, :y_c, :x_s, :v_d]),
   "PIEPLens"        => (init_PIEPLens,        [:x_c, :y_c, :x_s, :v_d, :eps, :p_a, :n_a]),
   "SIEMDLens"       => (init_SIEMDLens,       [:x_c, :y_c, :x_s, :v_d, :eps, :p_a]),
   "PJEMDLens"       => (init_PJEMDLens,       [:x_c, :y_c, :x_s, :x_t, :v_d, :eps, :p_a]),
   "HernquistLens"   => (init_HernquistLens,   [:x_c, :y_c, :x_s, :mass]),
   "NFWLens"         => (init_NFWLens,         [:x_c, :y_c, :x_s, :c, :mass])
)

# Constructor for composite lens
function init_CompositeLens(D_d::Real, lens::Vector{<:Any})
   # Empty component vector
   lens_components = AbstractLens[]

   # Run over the available lens components
   for component in lens
      # Check if the lens component is available
      if haskey(lens_init_functions, component.lens)            
         # Get the compoent init_* function and arguments
         init_func, init_args = lens_init_functions[component.lens]         
         
         # Create a dict of args to pass to the init_* function
         kwargs = Dict(arg => getproperty(component, arg) for arg in init_args)
         
         # Push the component into the vector
         push!(lens_components, init_func(; D_d, kwargs...))  
      else   
         throw(ArgumentError("Unknown lens type ** $(component.lens) **"))
      end
   end
   return init_CompositeLens(_components_=lens_components)
end


# Constructor for multi-plane lens
function init_MultiPlaneLens(lens::Vector{<:Any})   
   # Get all lens redshifts
   zl_all = Vector{Real}()
   for component in lens
      # Get the redshifts of each plane
      push!(zl_all, component.zl)
   end

   # Get (sorted) unique lens redshifts
   zl_unique::Vector{Real} = unique( sort(zl_all) )

   lens_components = Vector{AbstractLens}()
   indi_components = Vector{AbstractLens}()
   for z in zl_unique
      indi_components = []
      for component in lens
         if component.zl == z            
            push!(indi_components, component)
         end
      end
      push!(lens_components, init_CompositeLens(indi_components[1].D_d, indi_components))
   end
   return init_MultiPlaneLens(_np_=length(zl_unique), _zl_=zl_unique, _components_=lens_components)
end