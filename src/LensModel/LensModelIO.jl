module LensModelIO


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using YAML


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
export read_input
export ModelConfig
export Observation
export Parameter
export SourceModel
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
   reference::NTuple{2, Float64}
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
   key::Tuple{Symbol, Symbol} = (owner, name)
end


# --------------------------------------------------------------------------------------------------
# Abstract type: Source model
# --------------------------------------------------------------------------------------------------
@kwdef struct Knot <: AbstractLensConfig
   x::Vector{Float64}
   σx::Vector{Float64}
   y::Vector{Float64}
   σy::Vector{Float64}
end

@kwdef struct Source <: AbstractLensConfig
   knots::Vector{Knot}
end

@kwdef struct SourceModel <: AbstractLensConfig
   sources::Vector{Source}
end


# --------------------------------------------------------------------------------------------------
# Abstract type: Optimizer
# --------------------------------------------------------------------------------------------------
@kwdef struct NMConfig <: AbstractOptimizerConfig
   max_iter::Int64 = 10000
   tolerance::Float64 = 1e-6
end

@kwdef struct GDConfig <: AbstractOptimizerConfig
   max_iter::Int64 = 10000
   learning_rate::Float64 = 0.1
   tolerance::Float64 = 1e-6
end

@kwdef struct OptimizerConfig <: AbstractOptimizerConfig
   method::Symbol = :NM
   run_mode::String = "random"
   max_runs::Int64 = 100
   tolerance::Float64 = 1e-3
   config::AbstractOptimizerConfig = NMConfig()
end


# --------------------------------------------------------------------------------------------------
# Abstract type: MCMC
# --------------------------------------------------------------------------------------------------
# Metropolis-Hastings
@kwdef struct MHConfig <: AbstractMCMCConfig
   n_steps::Int = 10000
   proposal_scale::Float64 = 0.1
   burn_in::Int = 100
end

# Hamiltonian Monte Carlo
@kwdef struct HMCConfig <: AbstractMCMCConfig
   n_steps::Int = 10000
   step_size::Float64 = 0.1
   leapfrog_steps::Int = 10
   burn_in::Int = 100
end

# Affine-Invariant Ensemble Sampler
@kwdef struct AIESConfig <: AbstractMCMCConfig
   n_walkers::Int = 100
   n_steps::Int = 10000
   a::Float64 = 2.0
   burn_in::Int = 100
end

@kwdef struct MCMCConfig <: AbstractMCMCConfig
   method::Symbol = :MH
   n_chains::Int64 = 1
   config::AbstractMCMCConfig = MHConfig()
end

@kwdef struct SamplerConfig <: AbstractLensConfig
   scheme::Symbol = :SourcePlane
   optimizer::Union{Nothing, AbstractOptimizerConfig} = nothing
   mcmc::Union{Nothing, AbstractMCMCConfig} = nothing
end


# --------------------------------------------------------------------------------------------------
# Abstract type: Final Model Configuration
# --------------------------------------------------------------------------------------------------
@kwdef struct ModelConfig <: AbstractLensConfig
   observation::Observation
   cosmology::Cosmology.AbstractCosmology
   lens_model::Lenses.AbstractLens
   source_model::SourceModel
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

# Assign default value to optional key if it does not exist
@inline function _optional!(dict::Dict, key::Symbol, default::Union{Bool, String})
   if !haskey(dict, key)
      dict[key] = default
   end
   return nothing
end

# Internal function: Extract parameter values
function _extract_param_range(x::Union{Int64, Float64, Dict})::Tuple{Float64, Float64, Float64}
   # Extract parameter values in case of Int64 or Float64
   if x isa Union{Int64, Float64}
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
   obs_dict = dict[:observation]
   _optional!(obs_dict, :modeler, "LensFactory")
   _optional!(obs_dict, :lens, "WhoKnows")
   
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
   if FOV isa Union{Int64, Float64}
      FOV = (Float64(FOV), Float64(FOV))
   elseif FOV isa Vector
      FOV = (Float64(FOV[1]), Float64(FOV[2]))
   else
      error("Unknown type of FOV. Allowed types are Int64, Float64, Vector.")
   end

   # Construct struct with the given inputs
   obs = Observation(
      modeler     = obs_dict[:modeler],
      lens        = obs_dict[:lens],
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
      @inbounds for param in cosmo_params
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
function _lensmodel!(dict::Dict, params::Vector{Parameter})
   lens_model_dict = dict[:lens_model]
   
   # Check if we do single or multiplane lensing. Default: single plane
   _optional!(lens_model_dict, :multiplane, false)
   
   # Make sure that total number of lenses is greater than zero
   _require(lens_model_dict, :total_lenses)
   if lens_model_dict[:total_lenses] <= 0
      error("Total number of lenses must be greater than zero.")
   end

   # Construct a composite lens using initial values
   n_lenses = lens_model_dict[:total_lenses]
   lens_vector = NamedTuple[]

   # Reference point
   ref_ra  = dict[:observation][:reference][1]
   ref_dec = dict[:observation][:reference][2]
   use_ref = !(ref_ra == 0.0 && ref_dec == 0.0)

   if lens_model_dict[:multiplane] == false
      # Single plane lensing
      for i in 1:n_lenses
         lens_id = Symbol(:lens, i)
         lens_dict = lens_model_dict[lens_id]

         # Initialize lens parameter Dict
         lens_params = Dict{Symbol, Union{Symbol, Float64}}()
         lens_params[:lens] = Symbol(lens_dict[:lens])

         # --- Lens position parameters (always provided) ---
         rx, lx, ux = _extract_param_range(lens_dict[:x_c])
         ry, ly, uy = _extract_param_range(lens_dict[:y_c])
         
         if use_ref
            x_ref, y_ref = AstrometricOps.gnomonic_offsets_arcsec(ref_ra, ref_dec, rx, ry)
         else
            x_ref, y_ref = rx, ry
         end

         # Add lens position parameters to the lens_params Dict
         lens_params[:x_c] = x_ref
         lens_params[:y_c] = y_ref

         # Add lens position parameters to the parameter vector
         push!(params, Parameter(owner=lens_id, name=:x_c, refer=x_ref, lower=lx, upper=ux))
         push!(params, Parameter(owner=lens_id, name=:y_c, refer=y_ref, lower=ly, upper=uy))

         # --- Remaining lens parameters ---
         for (k, v) in lens_dict
            k ∈ (:lens, :x_c, :y_c) && continue

            # Extract parameter values and bounds
            r, l, u = _extract_param_range(v)

            # Add reference parameter to the lens_params Dict
            lens_params[k] = r

            # Add parameter to the reference, lower, and upper vectors
            push!(params, Parameter(owner=lens_id, name=k, refer=r, lower=l, upper=u))
         end
         # Store lens tuple
         push!(lens_vector, (; lens_params...))
      end
      composite_lens = Lenses.init_CompositeLens(lens_vector)
      return composite_lens
   else
      # Multi-plane lensing not yet implemented
      error("Multi-plane lensing support is not yet implemented.")
   end
end

# --------------------------------------------------------------------------------------------------
# ---------------- Read Source Model ---------------------------------------------------------------
# --------------------------------------------------------------------------------------------------  
function _source!(dict::Dict, cosmo::Cosmology.AbstractCosmology, params::Vector{Parameter})
   source_dict = dict[:source]
   # Make sure that total number of sources is present and greater than zero
   _require(source_dict, :total_sources)
   if source_dict[:total_sources] <= 0
      error("Total number of sources must be greater than zero.")
   end

   # Get number of sources
   n_source = source_dict[:total_sources]

   # Get reference lens redshift
   z_d = dict[:observation][:z_d]

   # Run over sources
   sources = Vector{Source}(undef, n_source)
   adis_ref = Vector{Float64}(undef, n_source)
   for i in 1:n_source
      source_id = Symbol(:source, i)
      individual_source_dict = source_dict[source_id]

      # Extract parameter values and bounds
      r, l, u = _extract_param_range(individual_source_dict[:z_s])

      # Convert redshift to adis
      D_ds = Cosmology.angular_diameter_distance(cosmo, z_d, r)
      D_os = Cosmology.angular_diameter_distance(cosmo, 0.0, r)
      adis_r = D_ds / D_os
      
      # Store reference adis
      adis_ref[i] = adis_r

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
      n_knot = individual_source_dict[:total_knots]
      knots = Vector{Knot}(undef, n_knot)
      for k in 1:n_knot
         x = individual_source_dict[Symbol(:x, k)]
         y = individual_source_dict[Symbol(:y, k)]

         # Assert that the number of knots is the same for x and y
         if length(x[:value]) != length(x[:sigma]) || 
            length(x[:value]) != length(y[:value]) || 
            length(x[:value]) != length(y[:sigma])
            error("Inconsistent knot dimensions in source-$i, knot-$k")
         end

         ref_ra = dict[:observation][:reference][1]
         ref_dec = dict[:observation][:reference][2]
         if ref_ra == 0.0 || ref_dec == 0.0
            knots[k] = Knot(x  = x[:value], σx = x[:sigma], y  = y[:value], σy = y[:sigma])
         else
            x_arcsec, y_arcsec = AstrometricOps.gnomonic_offsets_arcsec(ref_ra, ref_dec, x[:value], y[:value])
            knots[k] = Knot(x  = x_arcsec, σx = x[:sigma], y  = y_arcsec, σy = y[:sigma])
         end
      end

      sources[i] = Source(knots=knots)
   end
   return SourceModel(sources=sources)
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
   convergence = get(optimizer, :convergence, Dict())
   run_mode = get(convergence, :run_mode, :random)
   max_runs = get(convergence, :max_runs, 100)
   tolerance = get(convergence, :tolerance, 1e-3)
   
   algorithm_config = 
   if method == :NM
      NMConfig(
         max_iter    = config[:max_iter],
         tolerance   = config[:tolerance]
      )
   elseif method == :GD
      GDConfig(
         learning_rate = config[:learning_rate],
         max_iter      = config[:max_iter],
         tolerance     = config[:tolerance]
      )
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
   config = mcmc[method]
   
   algorithm_config = 
   if method == :MH
      MHConfig(
         n_steps        = config[:n_steps],
         proposal_scale = config[:proposal_scale],
         burn_in        = get(config, :burn_in, 100)
      )
   elseif method == :HMC
      HMCConfig(
         n_steps        = config[:n_steps],
         step_size      = config[:step_size],
         leapfrog_steps = config[:leapfrog_steps],
         burn_in        = get(config, :burn_in, 100)
      )
   elseif method == :AIES
      AIESConfig(
         n_walkers      = config[:n_walkers],
         n_steps        = config[:n_steps],
         a              = get(config, :a, 2.0),
         burn_in        = get(config, :burn_in, 100)
      )
   else
      error("Unknown MCMC method: $method")
   end

   return MCMCConfig(
      method   = method,
      config   = algorithm_config
   )
end

# Internal function: Process Lens Model section
function _sampling!(dict::Dict)
   sampling_dict = dict[:sampling]

   # Get sampler details
   _require(sampling_dict, :scheme)

   # Get optimizer details
   optimizer_params = _optimizer!(sampling_dict)

   # Get MCMC details
   mcmc_params = _mcmc!(sampling_dict)

   return SamplerConfig(
         scheme = Symbol(sampling_dict[:scheme]),
         optimizer = optimizer_params,
         mcmc = mcmc_params
      )
end

# --------------------------------------------------------------------------------------------------
# ---------------- Read the input file -------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
function read_input(filename::AbstractString)
   dict = YAML.load_file(filename)
   _symbolize!(dict)

   # Define reference, lower, and upper vectors
   params = Parameter[]

   # Get basic observation details
   observation = _observation(dict)

   # Get cosmology and its parameters
   cosmology = _cosmology!(dict, params)

   # Get lens model and its parameters
   lens = _lensmodel!(dict, params)

   # Get source model and its parameters
   source = _source!(dict, cosmology, params)

   # Get sampling details
   sampler = _sampling!(dict)

   # Identify free parameters
   free_param_idxs = findall(p -> p.lower != p.upper, params)

   return ModelConfig(
      observation = observation,
      cosmology = cosmology,
      lens_model = lens,
      source_model = source,
      parameters = params,
      free_param_idxs = free_param_idxs,
      sampler = sampler
   )
end

end