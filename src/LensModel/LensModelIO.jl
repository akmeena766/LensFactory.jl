module LensModelIO


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using YAML
using DelimitedFiles


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..Constants
using ..Cosmology
using ..Lenses
using ..LFUtils


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export ModelConfig
export Observation
export Parameter
export SourceConfig
export ScalingRelation
export LensConfig
export NMConfig
export GDConfig
export OptimizerConfig
export MHConfig
export AIESConfig
export SamplerConfig


# --------------------------------------------------------------------------------------------------
# Abstract types
# --------------------------------------------------------------------------------------------------
abstract type AbstractLensConfig end
abstract type AbstractOptimizerConfig <: AbstractLensConfig end
abstract type AbstractMCMCConfig <: AbstractLensConfig end


# --------------------------------------------------------------------------------------------------
# Abstract type: Observation
# --------------------------------------------------------------------------------------------------
@kwdef struct Observation <: AbstractLensConfig
   modeler::String
   lens::String
   z_d::Float64
   reference::NTuple{2,Float64}
   pixel_scale::Float64
   FOV::NTuple{2,Float64}
end

# --------------------------------------------------------------------------------------------------
# Abstract type: Parameter
# --------------------------------------------------------------------------------------------------
@kwdef struct Parameter <: AbstractLensConfig
   owner::Symbol
   name::Symbol
   refer::Float64
   lower::Float64
   upper::Float64
   key::Tuple{Symbol,Symbol} = (owner, name)
end

# --------------------------------------------------------------------------------------------------
# Abstract type: Lens component
# --------------------------------------------------------------------------------------------------
@kwdef struct LensComponent <: AbstractLensConfig
   owner::Symbol
   name::Symbol
end

@kwdef struct GalaxyComponent <: AbstractLensConfig
   n::Int64
   x_c::Vector{Float64}
   y_c::Vector{Float64}
   obs_mag::Vector{Float64}
   eps::Vector{Float64}
   pa::Vector{Float64}
end

@kwdef struct ScalingRelation <: AbstractLensConfig
   ref_mag::Float64
   ref_sigma::Float64
   ref_core::Float64
   ref_cut::Float64
   slope_sigma::Float64
   slope_core::Float64
   slope_cut::Float64
end

@kwdef struct LensConfig <: AbstractLensConfig
   components::Vector{LensComponent}
   galaxies::Union{Nothing, GalaxyComponent} = nothing
   scaling::Union{Nothing, ScalingRelation}  = nothing
end

# --------------------------------------------------------------------------------------------------
# Abstract type: Source model
# --------------------------------------------------------------------------------------------------
@kwdef struct Knot <: AbstractLensConfig
   x::Vector{Float64}
   y::Vector{Float64}
   σx::Vector{Float64}
   σy::Vector{Float64}
   σθ::Vector{Float64}
   use_parity::Bool = false
   parity::Vector{Int64} = []
end

@kwdef struct Source <: AbstractLensConfig
   knots::Vector{Knot}
end

@kwdef struct SourceConfig <: AbstractLensConfig
   sources::Vector{Source}
   use_parity::Bool      = false
   parity_force::Float64 = 0.0
end


# --------------------------------------------------------------------------------------------------
# Abstract type: Optimizer
# --------------------------------------------------------------------------------------------------
@kwdef struct NMConfig <: AbstractOptimizerConfig
   max_iter::Int64    = 10000
   tolerance::Float64 = 1e-6
end

@kwdef struct GDConfig <: AbstractOptimizerConfig
   max_iter::Int64        = 10000
   learning_rate::Float64 = 0.1
   tolerance::Float64     = 1e-6
end

@kwdef struct OptimizerConfig <: AbstractOptimizerConfig
   method::Symbol     = :NM
   run_mode::Symbol   = :random
   max_runs::Int64    = 100
   tolerance::Float64 = 1e-3
   config::AbstractOptimizerConfig = NMConfig()
end


# --------------------------------------------------------------------------------------------------
# Abstract type: MCMC
# --------------------------------------------------------------------------------------------------
# Metropolis-Hastings (MH) sampler
@kwdef struct MHConfig <: AbstractMCMCConfig
   n_steps::Int = 10000
   n_adapt::Int = div(n_steps, 10)
end

# Affine-Invariant Ensemble Sampler (AIES)
@kwdef struct AIESConfig <: AbstractMCMCConfig
   n_steps::Int64 = 100000
   a::Float64     = 2.0
end

@kwdef struct MCMCConfig <: AbstractMCMCConfig
   method::Symbol  = :AIES
   n_chains::Int64 = 1
   config::AbstractMCMCConfig = MHConfig()
end

@kwdef struct SamplerConfig <: AbstractLensConfig
   scheme::Symbol = :SourcePlane
   verbose::Bool  = true
   optimizer::Union{Nothing, AbstractOptimizerConfig} = nothing
   mcmc::Union{Nothing, AbstractMCMCConfig}           = nothing
end


# --------------------------------------------------------------------------------------------------
# Abstract type: Final Model Configuration
# --------------------------------------------------------------------------------------------------
@kwdef struct ModelConfig <: AbstractLensConfig
   observation::Observation
   cosmology::Cosmology.AbstractCosmology
   lens_config::LensConfig
   source_config::SourceConfig
   parameters::Vector{Parameter}
   free_param_idxs::Vector{Int64}
   sampler::SamplerConfig
end


# --------------------------------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------------------------------
# Efficient in-place conversion of Dict keys from String to Symbol
# Uses multiple dispatch to avoid type checking overhead
function _symbolize!(x::Dict)
   # Collect string keys first to avoid modification during iteration
   str_keys = String[]

   for k in keys(x)
      k isa String && push!(str_keys, k)
   end

   # Convert string keys to symbols and recursively symbolize values
   for k in str_keys
      v = x[k]
      delete!(x, k)
      x[Symbol(k)] = v
      _symbolize!(v)
   end
   return nothing
end

# Optimized vector handling using foreach
_symbolize!(x::Vector) = (foreach(_symbolize!, x); nothing)

# Base case: do nothing for other types
_symbolize!(x::Any) = nothing


# Ensure that required key exist in the given section
@inline function _require(dict::Dict, key::Symbol)::Bool
   return haskey(dict, key) || error("Missing required key ** $key ** in ** $dict **")
end


# Convert positions to arcsec offsets relative to the observation's reference point
@inline function _to_arcsec(observation::Observation, x, y)
   ref_ra, ref_dec = observation.reference

   if ref_ra == 0.0 && ref_dec == 0.0
      return x, y
   end

   return AstrometricOps.gnomonic_offsets_arcsec(ref_ra, ref_dec, x, y)
end


# Internal function: Extract parameter values
function _extract_param_range(x::Union{Int64,Float64,Dict})::Tuple{Float64,Float64,Float64}
   # Extract parameter values in case of Int64 or Float64
   if x isa Union{Int64,Float64}
      val = Float64(x)
      return val, val, val
   elseif x isa Dict
      # Make sure that value key exists
      haskey(x, :value) || error("Missing required key ** value ** in parameter definition.")

      # Extract reference value
      refer = Float64(x[:value])

      # Extract bounds (avoid repeated haskey calls)
      if haskey(x, :range)
         range = x[:range]
         lower = Float64(range[1])
         upper = Float64(range[2])
      else
         lower = refer
         upper = refer
      end

      return refer, lower, upper
   else
      throw(ArgumentError("Unknown parameter value type. Allowed types are Int64, Float64, Dict."))
   end
end

# --------------------------------------------------------------------------------------------------
# ---------------- Read Observation ----------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
function _observation(dict::Dict)
   # Read observation dict
   obs_dict = dict[:observation]

   # Required keywords
   _require(obs_dict, :z_d)
   _require(obs_dict, :reference)
   _require(obs_dict, :pixel_scale)
   _require(obs_dict, :FOV)

   # Make sure that reference is tuple and convert it to Float64
   ref_point = obs_dict[:reference]
   ref = (Float64(ref_point[1]), Float64(ref_point[2]))

   # Make sure that FOV size is a single value or tuple
   # If it is a single value --> Create a square grid based on that
   FOV = obs_dict[:FOV]
   if FOV isa Union{Int64,Float64}
      FOV = (Float64(FOV), Float64(FOV))
   elseif FOV isa Vector
      FOV = (Float64(FOV[1]), Float64(FOV[2]))
   else
      error("Unknown type of FOV. Allowed types are Int64, Float64, Vector.")
   end

   # Construct struct with the given inputs
   obs = Observation(
      modeler     = get(obs_dict, :modeler, "LensFactory"),
      lens        = get(obs_dict, :lens, "WhoKnows"),
      z_d         = obs_dict[:z_d],
      reference   = ref,
      pixel_scale = obs_dict[:pixel_scale],
      FOV         = FOV
   )
   return obs
end

# --------------------------------------------------------------------------------------------------
# ---------------- Read Cosmology ------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
function _cosmology!(input_dict::Dict, params::Vector{Parameter})
   # Manual tuple of symbols for cosmology part
   cosmo_params = (:H0, :w, :Omega_m0, :Omega_r0, :Omega_w0, :Omega_k0)

   # Get default cosmology
   cosmology = Cosmology.init_cosmology()

   # Check if cosmology section exists otherwise throw a warning and fall back to default cosmology
   if haskey(input_dict, :cosmology)
      cosmo_dict = input_dict[:cosmology]
      for param in cosmo_params
         if haskey(cosmo_dict, param)
            r, l, u = _extract_param_range(cosmo_dict[param])
         else
            # Use default cosmology value
            default_value = getfield(cosmology, param)
            r = l = u = default_value
         end
         push!(params, Parameter(owner=:cosmology, name=param, refer=r, lower=l, upper=u))
      end
   else
      @warn "Cosmology missing in input - using default values. 
      See example input file here: https://github.com/akmeena766/LensFactory_Examples."
      # Add default cosmology parameters
      @inbounds for param in cosmo_params
         default_value = getfield(cosmology, param)
         push!(params, Parameter(owner=:cosmology, name=param, refer=default_value, lower=default_value, upper=default_value))
      end
   end
   return cosmology
end

# --------------------------------------------------------------------------------------------------
# ---------------- Read Lens Model -----------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
const NO_POSITION     = Set([:ExternalEffects, :ExternalEffects3, :Multipole])
const REQUIRE_ADD     = Set([:PointLens, :PlummerLens, :GaussianLens, :SersicLens, :HernquistLens, 
                             :NFWLens, :tNFWLens, :gNFWLens, :EinastoLens, :aHernquistLens, 
                             :aNFWLens, :eHernquistMDLens, :eNFWMDLens, 
                             :MultiPlummerLens, :MultiGaussianLens])
const REQUIRE_SCALING = Set([:MultiPJELens])
function _lensmodel!(dict::Dict, params::Vector{Parameter}, observation::Observation, Dol_ref::Float64)
   lens_dict = dict[:lens_model]

   # Make sure that total number of lenses is greater than zero
   _require(lens_dict, :total_lenses)
   if lens_dict[:total_lenses] <= 0
      error("Total number of lenses must be greater than zero.")
   end

   # Construct a composite lens using initial values
   n_lenses = lens_dict[:total_lenses]

   # Initialize lens name vector
   lens_name = Vector{LensComponent}(undef, n_lenses)

   # Initialize galaxy component and scaling flag
   galaxy_comp = nothing
   needs_scaling = false

   # Single plane vs. multiplane lensing
   multiplane = get(lens_dict, :multiplane, false)
   if multiplane == false
      # Single plane lensing
      for i in 1:n_lenses
         lens_id = Symbol(:lens, i)
         indi_lens_dict = lens_dict[lens_id]
         name = Symbol(indi_lens_dict[:lens])

         # Store lens model name in lens_name vector
         lens_name[i] = LensComponent(owner=lens_id, name=name)

         # Add distance parameters
         if name ∈ REQUIRE_ADD
            push!(params, Parameter(owner=lens_id, name=:D_d, refer=Dol_ref, lower=Dol_ref, upper=Dol_ref))
         end

         if name ∉ REQUIRE_SCALING
            # --- Lens position parameters (always provided) ---
            if name ∉ NO_POSITION
               rx, lx, ux = _extract_param_range(indi_lens_dict[:x_c])
               ry, ly, uy = _extract_param_range(indi_lens_dict[:y_c])

               # Check if a valid (RA, Dec) is provided as reference or (0, 0) is used
               # reference = (0, 0) ⇒ Lens positions are provided in arcseconds
               # reference = (RA, Dec) ⇒ Lens positions are provided in RA and Dec. Conversion needed.
               x_lens, y_lens = _to_arcsec(observation, rx, ry)

               # Add lens position parameters to the parameter vector
               push!(params, Parameter(owner=lens_id, name=:x_c, refer=x_lens, lower=lx, upper=ux))
               push!(params, Parameter(owner=lens_id, name=:y_c, refer=y_lens, lower=ly, upper=uy))
            end

            # --- Remaining lens parameters ---
            for (k, v) in indi_lens_dict
               k ∈ (:lens, :x_c, :y_c) && continue

               # Extract parameter values and bounds
               r, l, u = _extract_param_range(v)

               # Add parameter to the reference, lower, and upper vectors
               push!(params, Parameter(owner=lens_id, name=k, refer=r, lower=l, upper=u))
            end
         else
            # Update scaling flag
            needs_scaling = true

            # Get file name
            file_name = indi_lens_dict[:galaxy_file]

            # Load catalog and create the GalaxyComponent
            catalog_data = readdlm(file_name)
            catalog_data = Float64.(catalog_data)

            # Check if a valid (RA, Dec) is provided as reference or (0, 0) is used
            # reference = (0, 0) ⇒ Lens positions are provided in arcseconds
            # reference = (RA, Dec) ⇒ Lens positions are provided in RA and Dec. Conversion needed.
            x_lens, y_lens = _to_arcsec(observation, catalog_data[:, 2], catalog_data[:, 3])

            # Total number of galaxies
            catalog_n = size(catalog_data, 1)

            # Get fixed galaxy components
            galaxy_comp = GalaxyComponent(
               n       = catalog_n,
               x_c     = x_lens,
               y_c     = y_lens, 
               obs_mag = catalog_data[:, 4],
               eps     = catalog_data[:, 5],
               pa      = catalog_data[:, 6]
            )
         end
      end

      # Initialize scaling object
      scaling_obj = nothing

      # Get scaling relations
      if needs_scaling
         # Get scaling relations dictionary
         scaling_dict = dict[:scaling_relation]
         
         # Owner of scaling parameters
         owner = :scaling
         
         # Get reference magnitude
         rm, lm, um = _extract_param_range(scaling_dict[:ref_mag])
         push!(params, Parameter(owner=:scaling, name=:ref_mag, refer=rm, lower=lm, upper=um))

         # Get reference sigma_star
         rs, ls, us = _extract_param_range(scaling_dict[:ref_sigma])
         push!(params, Parameter(owner=:scaling, name=:ref_sigma, refer=rs, lower=ls, upper=us))

         # Get reference rcore_star
         rc, lc, uc = _extract_param_range(scaling_dict[:ref_core])
         push!(params, Parameter(owner=:scaling, name=:ref_core, refer=rc, lower=lc, upper=uc))

         # Get reference rcut_star
         rt, lt, ut = _extract_param_range(scaling_dict[:ref_cut])
         push!(params, Parameter(owner=:scaling, name=:ref_cut, refer=rt, lower=lt, upper=ut))

         # Get reference slope_sigma
         ss, ls, us = _extract_param_range(scaling_dict[:slope_sigma])
         push!(params, Parameter(owner=:scaling, name=:slope_sigma, refer=ss, lower=ls, upper=us))

         # Get reference slope_core
         sc, lc, uc = _extract_param_range(scaling_dict[:slope_core])
         push!(params, Parameter(owner=:scaling, name=:slope_core, refer=sc, lower=lc, upper=uc))

         # Get reference slope_cut
         st, lt, ut = _extract_param_range(scaling_dict[:slope_cut])
         push!(params, Parameter(owner=:scaling, name=:slope_cut, refer=st, lower=lt, upper=ut))

         scaling_obj = ScalingRelation(
            ref_mag     = rm,
            ref_sigma   = rs,
            ref_core    = rc,
            ref_cut     = rt,
            slope_sigma = ss,
            slope_core  = sc,
            slope_cut   = st
         )
      end
      
      return LensConfig(components=lens_name, galaxies=galaxy_comp, scaling=scaling_obj)
   else
      # Multi-plane lensing not yet implemented
      error("Multi-plane lensing support is not yet implemented.")
   end
end


# --------------------------------------------------------------------------------------------------
# ---------------- Read Source Model ---------------------------------------------------------------
# --------------------------------------------------------------------------------------------------  
function _source_direct!(dict::Dict, cosmo::Cosmology.AbstractCosmology, observation::Observation, params::Vector{Parameter})
   # Get source dictionary from the full dictionary
   source_dict = dict[:source]

   # Make sure that total number of sources is present and greater than zero
   _require(source_dict, :total_sources)
   if source_dict[:total_sources] <= 0
      error("Total number of sources must be greater than zero.")
   end

   # Check if parity is enforced
   use_parity = get!(source_dict, :use_parity, false)
   parity_force = 0.0
   if use_parity
      parity_force = source_dict[:parity_force]
   end

   # Get number of sources
   n_source = source_dict[:total_sources]

   # Get reference lens redshift
   z_d = dict[:observation][:z_d]

   # Run over sources
   sources = Vector{Source}(undef, n_source)
   for i in 1:n_source
      source_id = Symbol(:source, i)
      indi_source_dict = source_dict[source_id]

      # Extract parameter values and bounds
      r, l, u = _extract_param_range(indi_source_dict[:z_s])

      # Convert redshift to adis
      D_ds = Cosmology.angular_diameter_distance(cosmo, z_d, r)
      D_os = Cosmology.angular_diameter_distance(cosmo, 0.0, r)
      adis_r = D_ds / D_os

      if l == u
         adis_l = adis_r
         adis_u = adis_r
      else
         D_ds = Cosmology.angular_diameter_distance(cosmo, z_d, l)
         D_os = Cosmology.angular_diameter_distance(cosmo, 0.0, l)
         adis_l = D_ds / D_os

         D_ds = Cosmology.angular_diameter_distance(cosmo, z_d, u)
         D_os = Cosmology.angular_diameter_distance(cosmo, 0.0, u)
         adis_u = D_ds / D_os
      end

      # Add redshift parameter to the reference, lower, and upper vectors
      push!(params, Parameter(owner=source_id, name=Symbol(:adis, i), refer=adis_r, lower=adis_l, upper=adis_u))

      # Run over knots
      n_knot = indi_source_dict[:total_knots]
      knots = Vector{Knot}(undef, n_knot)
      for k in 1:n_knot
         knot_id = Symbol(:knot, k)
         # Read knot image position(s)
         x = indi_source_dict[knot_id][:x]
         y = indi_source_dict[knot_id][:y]

         # Assert that the number of values is the same for (x, y, σx, σy, θ)
         if length(x) != length(y)
            error("Inconsistent knot position dimensions in source-$i, knot-$k")
         end

         # Read knot image error(s)
         if haskey(indi_source_dict[knot_id], :sigma)
            σx = indi_source_dict[knot_id][:sigma]
            σy = indi_source_dict[knot_id][:sigma]
            σθ = 0.0 .* indi_source_dict[knot_id][:sigma]
         elseif haskey(indi_source_dict[knot_id], :sigma_x) &&
                haskey(indi_source_dict[knot_id], :sigma_y) &&
                haskey(indi_source_dict[knot_id], :sigma_theta)
            σx = indi_source_dict[knot_id][:sigma_x]
            σy = indi_source_dict[knot_id][:sigma_y]
            σθ = indi_source_dict[knot_id][:sigma_theta]
         else
            error("Invalid knot error parameters in source-$i, knot-$k")
         end
         # Assert that the number of values is the same for (x, y, σx, σy, σθ)
         if length(x) != length(σx) || length(x) != length(σy) || length(x) != length(σθ)
            error("Inconsistent knot position and error dimensions in source-$i, knot-$k")
         end

         # Read knot parity if use_parity is true and parity is present
         if use_parity && haskey(indi_source_dict[knot_id], :parity)
            p_temp = indi_source_dict[knot_id][:parity]
            if length(p_temp) != length(x)
               error("Inconsistent knot parity dimensions in source-$i, knot-$k")
            end

            # Check if all values are 1 or -1. otherwise set knot_parity to false
            if all(abs.(p_temp) .== 1)
               knot_parity = true
               p = p_temp
            else
               knot_parity = false
               p = 0 .* p_temp
            end
         else
            knot_parity = false
            p = 0 .* x
         end


         # Convert to arcsec if reference is not (0, 0)
         x_arcsec, y_arcsec = _to_arcsec(observation, x, y)
         knots[k] = Knot(x=x_arcsec, σx=σx, y=y_arcsec, σy=σy, σθ=deg2rad.(σθ), use_parity=knot_parity, parity=p)
      end
      sources[i] = Source(knots=knots)
   end
   return SourceConfig(sources=sources, use_parity=use_parity, parity_force=parity_force)
end


function _source_from_file!(dict::Dict, cosmo::Cosmology.AbstractCosmology, observation::Observation, params::Vector{Parameter})
   # Get source dictionary from the full dictionary
   source_dict = dict[:source]

   # Get file name
   file_name = source_dict[:source_file]

   file_data = readdlm(file_name)
   file_data = Float64.(file_data)

   # --- Parameters directly read from the YAML file -----------------------------------------------
   # Get total number of sources from first column in file
   n_source = Int64(maximum(file_data[:, 1]))

   # Check if parity is enforced
   use_parity = get!(source_dict, :use_parity, false)
   parity_force = 0.0
   if use_parity
      parity_force = source_dict[:parity_force]
   end

   # Get reference lens redshift
   z_d = dict[:observation][:z_d]

   # --- Parameters read from the source file ------------------------------------------------------
   sources = Vector{Source}(undef, n_source)
   for i in 1:n_source
      # Generate source ID
      source_id = Symbol(:source, i)

      # Get individual source data from file
      mask = file_data[:, 1] .== i
      source_data = file_data[mask, :]

      # Get source redshift
      if haskey(source_dict, source_id)
         # Get individual source dict from YAML file
         indi_source_dict = source_dict[source_id]

         # Extract parameter values and bounds
         r, l, u = _extract_param_range(indi_source_dict[:z_s])

         # Convert redshift to adis
         D_ds = Cosmology.angular_diameter_distance(cosmo, z_d, r)
         D_os = Cosmology.angular_diameter_distance(cosmo, 0.0, r)
         adis_r = D_ds / D_os

         if l == u
            adis_l = adis_r
            adis_u = adis_r
         else
            D_ds = Cosmology.angular_diameter_distance(cosmo, z_d, l)
            D_os = Cosmology.angular_diameter_distance(cosmo, 0.0, l)
            adis_l = D_ds / D_os

            D_ds = Cosmology.angular_diameter_distance(cosmo, z_d, u)
            D_os = Cosmology.angular_diameter_distance(cosmo, 0.0, u)
            adis_u = D_ds / D_os
         end

         # Add redshift parameter to the reference, lower, and upper vectors
         push!(params, Parameter(owner=source_id, name=Symbol(:adis, i), refer=adis_r, lower=adis_l, upper=adis_u))
      else
         # Get the source redshift from the file
         z_s = source_data[1, 5]

         # Convert redshift to adis
         D_ds = Cosmology.angular_diameter_distance(cosmo, z_d, z_s)
         D_os = Cosmology.angular_diameter_distance(cosmo, 0.0, z_s)
         adis_r = D_ds / D_os

         # Redshift fixed. Set lower and upper bounds to the same value
         adis_l = adis_r
         adis_u = adis_r

         # Add redshift parameter to the reference, lower, and upper vectors
         push!(params, Parameter(owner=source_id, name=Symbol(:adis, i), refer=adis_r, lower=adis_l, upper=adis_u))
      end

      # Run over knots
      n_knot = Int64(maximum(source_data[:, 2]))
      knots = Vector{Knot}(undef, n_knot)
      for k in 1:n_knot
         # Generate knot id
         knot_id = Symbol(:knot, k)

         # Get individual knot data from file
         mask = source_data[:, 2] .== k
         knot_data = source_data[mask, :]

         # Get knot parameters
         x = knot_data[:, 3]
         y = knot_data[:, 4]

         σx = knot_data[:, 6]
         σy = knot_data[:, 7]
         σθ = knot_data[:, 8]

         # Read knot parity if use_parity is true otherwise assign zero
         if use_parity
            # Get parity values from file
            p_temp = knot_data[:, 9]
            # Check if all values are 1 or -1. otherwise set knot_parity to false
            if all(abs.(p_temp) .== 1)
               knot_parity = true
               p = p_temp
            else
               knot_parity = false
               p = 0 .* p_temp
            end
         else
            knot_parity = false
            p = 0 .* x
         end

         # Convert to arcsec if reference is not (0, 0)
         x_arcsec, y_arcsec = _to_arcsec(observation, x, y)
         knots[k] = Knot(x=x_arcsec, σx=σx, y=y_arcsec, σy=σy, σθ=deg2rad.(σθ), use_parity=knot_parity, parity=p)
      end
      sources[i] = Source(knots=knots)
   end
   return SourceConfig(sources=sources, use_parity=use_parity, parity_force=parity_force)
end

function _source!(dict::Dict, cosmo::Cosmology.AbstractCosmology, observation::Observation, params::Vector{Parameter})
   # Get source dictionary from the full dictionary
   source_dict = dict[:source]

   # Determine if we need to read source details from a file
   from_file = get(source_dict, :from_file, false)

   if from_file
      source_config = _source_from_file!(dict, cosmo, observation, params)
   else
      source_config = _source_direct!(dict, cosmo, observation, params)
   end
   return source_config
end

# --------------------------------------------------------------------------------------------------
# ---------------- Read Optimizer ------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
const OPTIMIZER_META_KEYS = Set([:enabled, :convergence])

function _infer_optimizer_method(dict::Dict)
   found = nothing
   for k in keys(dict)
      k in OPTIMIZER_META_KEYS && continue
      if found !== nothing
         error("Multiple entries specified: $(found), $(k)")
      end
      found = k
   end
   found === nothing && error("No entry specified")
   return found
end

function _optimizer!(sampling_dict::Dict)
   optimizer = get(sampling_dict, :optimizer, nothing)

   # Check if optimizer is missing or explicitly disabled
   if optimizer === nothing || !get(optimizer, :enabled, true)
      return nothing
   end

   # Infer optimizer method from keys
   method = _infer_optimizer_method(optimizer)
   config = optimizer[method]

   # Get convergence parameters (with defaults)
   convergence = get(optimizer,   :convergence, Dict())
   run_mode    = get(convergence, :run_mode,    :random)
   max_runs    = get(convergence, :max_runs,    100)
   tolerance   = get(convergence, :tolerance,   1.0E-3)

   algorithm_config =
      if method == :NM
         NMConfig(; config...)
      elseif method == :GD
         GDConfig(; config...)
      else
         error("Unknown optimizer method: $method")
      end

   return OptimizerConfig(
      method    = method,
      run_mode  = run_mode,
      max_runs  = max_runs,
      tolerance = tolerance,
      config    = algorithm_config
   )
end

# --------------------------------------------------------------------------------------------------
# ---------------- Read MCMC -----------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
const MCMC_META_KEYS = Set([:enabled, :n_chains])

function _infer_mcmc_method(dict::Dict)
   found = nothing
   for k in keys(dict)
      k in MCMC_META_KEYS && continue
      if found !== nothing
         error("Multiple entries specified: $(found), $(k)")
      end
      found = k
   end
   found === nothing && error("No entry specified")
   return found
end

function _mcmc!(sampling_dict::Dict)
   mcmc = get(sampling_dict, :mcmc, nothing)

   # Check if mcmc is missing or explicitly disabled
   if mcmc === nothing || !get(mcmc, :enabled, true)
      return nothing
   end

   # Infer mcmc method from keys
   method = _infer_mcmc_method(mcmc)
   n_chains = get(mcmc, :n_chains, 1)
   config = mcmc[method]

   algorithm_config =
      if method == :MH
         MHConfig(; config...)
      elseif method == :AIES
         AIESConfig(; config...)
      else
         error("Unknown MCMC method: $method")
      end

   return MCMCConfig(
      method   = method,
      n_chains = n_chains,
      config   = algorithm_config
   )
end

# Internal function: Process Lens Model section
const SAMPLING_SCHEMES = Set([:SourcePlane, :ImagePlane_Fast])
function _sampling!(dict::Dict)
   sampling_dict = dict[:sampling]

   # Get sampler details
   _require(sampling_dict, :scheme)
   scheme = Symbol(sampling_dict[:scheme])
   
   # Validate sampler scheme
   if !(scheme in SAMPLING_SCHEMES)
      error("Unsupported sampling scheme: $(scheme). Supported schemes: $(SAMPLING_SCHEMES)")
   end
   
   # Get verbose flag
   verbose = get!(sampling_dict, :verbose, true)

   # Get optimizer details
   optimizer_params = _optimizer!(sampling_dict)

   # Get MCMC details
   mcmc_params = _mcmc!(sampling_dict)

   return SamplerConfig(
      scheme=Symbol(sampling_dict[:scheme]),
      verbose=verbose,
      optimizer=optimizer_params,
      mcmc=mcmc_params
   )
end

# --------------------------------------------------------------------------------------------------
# Read the input YAML file and construct the ModelConfig struct
# --------------------------------------------------------------------------------------------------
function _read_input(filename::AbstractString)
   dict = YAML.load_file(filename)
   _symbolize!(dict)

   # Define reference, lower, and upper vectors
   params = Parameter[]

   # Get basic observation details
   observation = _observation(dict)

   # Get cosmology and its parameters
   cosmology = _cosmology!(dict, params)

   # Calculate reference angular diameter distance
   Dol_ref = Cosmology.angular_diameter_distance(cosmology, 0.0, observation.z_d)

   # Get lens model and its parameters
   lens_config = _lensmodel!(dict, params, observation, Dol_ref)

   # Get source model and its parameters
   source_config = _source!(dict, cosmology, observation, params)

   # Get sampling details
   sampler = _sampling!(dict)

   # Identify free parameters
   free_param_idxs = findall(p -> p.lower != p.upper, params)

   return ModelConfig(
      observation=observation,
      cosmology=cosmology,
      lens_config=lens_config,
      source_config=source_config,
      parameters=params,
      free_param_idxs=free_param_idxs,
      sampler=sampler
   )
end

end