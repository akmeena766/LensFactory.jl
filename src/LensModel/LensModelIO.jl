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
# NOTE: Structs use @kwdef for keyword construction but carry NO input defaults - every field is a
# required keyword. All default values are defined in the IO reader functions below (via get(...)).
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
# Abstract type: Observation
# --------------------------------------------------------------------------------------------------
@kwdef struct Observation <: AbstractLensConfig
   modeler::String
   lens::String
   z_d::Float64
   D_d::Float64
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
   refer::Real
   lower::Real
   upper::Real
   key::Tuple{Symbol, Symbol} = (owner, name)
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
   galaxies::Union{Nothing, GalaxyComponent}
   scaling::Union{Nothing, ScalingRelation}
end

# --------------------------------------------------------------------------------------------------
# Abstract type: Source model
# --------------------------------------------------------------------------------------------------
@kwdef struct Knot <: AbstractLensConfig
   x::Vector{Float64}    # Measured x-position
   y::Vector{Float64}    # Measured y-position
   σx::Vector{Float64}   # Positional error along semi-major axis
   σy::Vector{Float64}   # Positional error along semi-minor axis
   σθ::Vector{Float64}   # Positional error ellipse PA (CCW wrt x-axis)
   parity::Vector{Int64} # Parity vector
   m::Vector{Float64}    # Observed image magnitudes
   σm::Vector{Float64}   # Magnitude errors
   td::Vector{Float64}   # Observed time delays (days; arbitrary zero-point)
   σ_td::Vector{Float64} # Time-delay errors (days)
end

@kwdef struct Source <: AbstractLensConfig
   knots::Vector{Knot}
end

@kwdef struct SourceConfig <: AbstractLensConfig
   sources::Vector{Source}
   use_parity::Bool
   parity_force::Float64
   use_flux::Bool
   use_time_delay::Bool
   pixel_fd::Float64
end


# --------------------------------------------------------------------------------------------------
# Abstract type: Optimizer
# --------------------------------------------------------------------------------------------------
@kwdef struct NMConfig <: AbstractOptimizerConfig
   max_iter::Int64
   tolerance::Float64
end

@kwdef struct GDConfig <: AbstractOptimizerConfig
   max_iter::Int64
   learning_rate::Float64
   tolerance::Float64
end

@kwdef struct OptimizerConfig <: AbstractOptimizerConfig
   method::Symbol
   run_mode::Symbol
   max_runs::Int64
   tolerance::Float64
   config::AbstractOptimizerConfig
end


# --------------------------------------------------------------------------------------------------
# Abstract type: MCMC
# --------------------------------------------------------------------------------------------------
# Metropolis-Hastings (MH) sampler
@kwdef struct MHConfig <: AbstractMCMCConfig
   n_walkers::Int64
   n_steps::Int64
   n_adapt::Int64
end

# Affine-Invariant Ensemble Sampler (AIES)
@kwdef struct AIESConfig <: AbstractMCMCConfig
   n_walkers::Int64
   n_steps::Int64
   a::Float64
end

@kwdef struct MCMCConfig <: AbstractMCMCConfig
   method::Symbol
   config::AbstractMCMCConfig
end

@kwdef struct SamplerConfig <: AbstractLensConfig
   scheme::Symbol
   verbose::Bool
   optimizer::Union{Nothing, AbstractOptimizerConfig}
   mcmc::Union{Nothing, AbstractMCMCConfig}
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

@inline function _adis_ratio(cosmo::Cosmology.AbstractCosmology, z_d::RV, z_s::RV)
   D_ds = Cosmology.angular_diameter_distance(cosmo, z_d, z_s)
   D_os = Cosmology.angular_diameter_distance(cosmo, 0.0, z_s)
   return D_ds / D_os
end

function _zs_to_adis(cosmo::Cosmology.AbstractCosmology, z_d::RV, r::RV, l::RV, u::RV)
   adis_r = _adis_ratio(cosmo, z_d, r)

   if l == u
      adis_l = adis_r
      adis_u = adis_r
   else
      adis_l = _adis_ratio(cosmo, z_d, l)
      adis_u = _adis_ratio(cosmo, z_d, u)
   end

   return adis_r, adis_l, adis_u
end

# Extract parameter values
function _extract_param_range(x::Union{Int64, Float64, Dict})::Tuple{Real, Real, Real}
   # Extract parameter values in case of Int64 or Float64
   if x isa Union{Int64,Float64}
      val = x
      return val, val, val
   elseif x isa Dict
      # Make sure that value key exists
      haskey(x, :value) || error("Missing required key ** value ** in parameter definition.")

      # Extract reference value
      refer = x[:value]

      # Extract bounds (avoid repeated haskey calls)
      if haskey(x, :range)
         range = x[:range]
         lower = range[1]
         upper = range[2]
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
function _observation(dict::Dict,  cosmo::Cosmology.AbstractCosmology)
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

   # Angular diameter distance to the lens - computed once here so it is not recalculated
   D_d = Cosmology.angular_diameter_distance(cosmo, 0.0, Float64(obs_dict[:z_d]))

   # Construct struct with the given inputs
   obs = Observation(
      modeler     = get(obs_dict, :modeler, "LensFactory"),
      lens        = get(obs_dict, :lens, "WhoKnows"),
      z_d         = obs_dict[:z_d],
      D_d         = D_d,
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

            # Free cosmological parameters are not supported (yet): distances cached in
            # Observation and lens D_d / cosmology objects are fixed at the reference cosmology,
            # so freeing them would give an inconsistent model.
            if l != u
               error("Cosmological parameter ** $param ** cannot be free (range given). " *
                     "Free cosmological parameters are not supported yet.")
            end
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
function _lensmodel!(dict::Dict, params::Vector{Parameter}, observation::Observation)
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
   multiplane = get!(lens_dict, :multiplane, false)
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
            push!(params, Parameter(owner = lens_id, 
                                    name  = :D_d, 
                                    refer = observation.D_d, 
                                    lower = observation.D_d, 
                                    upper = observation.D_d))
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
            catalog_data = readdlm(file_name, comments=true, comment_char='#')
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
         
         # Names of all scaling-relation parameters, in the order expected by ScalingRelation
         scaling_params = (:ref_mag, :ref_sigma, :ref_core, :ref_cut, :slope_sigma, :slope_core, :slope_cut)

         # Reference values, keyed by parameter name, used to build the ScalingRelation struct
         refer_values = Dict{Symbol,Float64}()

         for param in scaling_params
            r, l, u = _extract_param_range(scaling_dict[param])
            refer_values[param] = r
            push!(params, Parameter(owner=:scaling, name=param, refer=r, lower=l, upper=u))
         end

         scaling_obj = ScalingRelation(; refer_values...)
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
# Convert a YAML list entry into a Vector{Float64} / Vector{Int64}. Scalars are deliberately
# rejected: measurements and errors must be fully specified, one entry per image, even if identical.
@inline _as_float_vec(x::AbstractVector) = Float64.(x)
@inline _as_float_vec(x::Any) = error("Expected a list with one entry per image, got: ** $x **")

@inline _as_int_vec(x::AbstractVector) = Int64.(x)
@inline _as_int_vec(x::Any) = error("Expected a list with one entry per image, got: ** $x **")

# Read an optional per-knot measurement (value + error), e.g. flux (magnitudes) or time delays.
# Both quantities are RELATIVE measurements (the source magnitude / delay zero-point is profiled
# out in the likelihood), therefore:
#   - One entry per image is required, aligned with the image order of (x, y).
#   - Images without a measurement are marked with error <= 0 (their value entry is ignored).
#   - At least two measured images are required; a single measurement carries no information
#     and the knot is treated as having no measurement (with a warning).
# Returns a pair of empty vectors when the measurement is absent or uninformative.
function _knot_measurement(knot_dict::Dict, val_key::Symbol, err_key::Symbol, n_img::Int64, i::Int64, k::Int64)
   # Measurement not provided --> return empty vectors
   if !haskey(knot_dict, val_key)
      return Float64[], Float64[]
   end

   # Error entry is required whenever the value entry is present
   if !haskey(knot_dict, err_key)
      error("Missing ** $err_key ** for ** $val_key ** in source-$i, knot-$k")
   end

   v = _as_float_vec(knot_dict[val_key])
   σ = _as_float_vec(knot_dict[err_key])

   # One entry per image, always
   if length(v) != n_img || length(σ) != n_img
      error("** $val_key / $err_key ** must have one entry per image ($n_img) in source-$i, knot-$k. " *
            "Mark images without a measurement using $err_key <= 0.")
   end

   # Count measured images (error > 0)
   n_meas = count(>(0.0), σ)
   if n_meas == 0
      return Float64[], Float64[]
   elseif n_meas == 1
      @warn "Only one measured image for ** $val_key ** in source-$i, knot-$k. Relative " *
            "measurements need at least two measured images - this knot will be ignored."
      return Float64[], Float64[]
   end
   return v, σ
end

function _source_direct!(dict::Dict, cosmo::Cosmology.AbstractCosmology, observation::Observation, params::Vector{Parameter})
   # Get source dictionary from the full dictionary
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

   # Track whether any knot carries flux / time-delay measurements
   any_parity = false
   any_flux   = false
   any_td     = false

   # Run over sources
   sources = Vector{Source}(undef, n_source)
   for i in 1:n_source
      source_id = Symbol(:source, i)
      indi_source_dict = source_dict[source_id]

      # Extract parameter values and bounds
      r, l, u = _extract_param_range(indi_source_dict[:z_s])

      # Convert redshift to adis
      adis_r, adis_l, adis_u = _zs_to_adis(cosmo, z_d, r, l, u)

      # Add adis parameter to the reference, lower, and upper vectors
      push!(params, Parameter(owner = source_id, 
                              name  = Symbol(:adis, i), 
                              refer = adis_r, 
                              lower = adis_l, 
                              upper = adis_u))

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

         # Read optional parity values: one entry per image, +1 / -1, with 0 marking images
         # without parity information. Parity is an absolute per-image quantity, so a single
         # measured image is already informative (no minimum count).
         p = zeros(Int64, length(x))
         if haskey(indi_source_dict[knot_id], :parity)
            p_temp = indi_source_dict[knot_id][:parity]
            if length(p_temp) != length(x)
               error("** parity ** must have one entry per image ($(length(x))) in source-$i, knot-$k. " *
                     "Mark images without parity information with 0.")
            end
            if !all(v -> v in (-1, 0, 1), p_temp)
               error("** parity ** values must be +1, -1, or 0 (no information) in source-$i, knot-$k")
            end
            p = p_temp
         end
         if any(!=(0), p)
            any_parity = true
         end

         # Read optional flux measurements (magnitudes) and their errors
         m, σm = _knot_measurement(indi_source_dict[knot_id], :mag, :sigma_mag, length(x), i, k)
         if !isempty(m)
            any_flux = true
         end

         # Read optional time-delay measurements (days, arbitrary zero-point) and their errors
         td, σ_td = _knot_measurement(indi_source_dict[knot_id], :td, :sigma_td, length(x), i, k)
         if !isempty(td)
            any_td = true
         end

         # Convert to arcsec if reference is not (0, 0)
         x_arcsec, y_arcsec = _to_arcsec(observation, x, y)
         knots[k] = Knot(x      = x_arcsec, 
                         σx     = σx, 
                         y      = y_arcsec, 
                         σy     = σy, 
                         σθ     = deg2rad.(σθ), 
                         parity = p,
                         m      = m, 
                         σm     = σm, 
                         td     = td, 
                         σ_td   = σ_td)
      end
      sources[i] = Source(knots=knots)
   end
   # Enable parity / flux / time-delay terms if data is present, unless explicitly disabled in the input.
   # Parity
   if get(source_dict, :use_parity, false)
      use_parity = any_parity
   else
      use_parity = false
   end

   # Flux
   if get(source_dict, :use_flux, false)
      use_flux = any_flux
   else
      use_flux = false
   end

   # Time Delay
   if get(source_dict, :use_time_delay, false)
      use_td = any_td
   else
      use_td = false
   end

   if get(source_dict, :use_parity, false) === true && !any_parity
      @warn "use_parity is set but no knot provides ** parity ** values - parity term disabled."
   end

   if get(source_dict, :use_flux, false) === true && !any_flux
      @warn "use_flux is set but no knot provides ** mag / sigma_mag ** - flux chi2 disabled."
   end

   if get(source_dict, :use_time_delay, false) === true && !any_td
      @warn "use_time_delay is set but no knot provides ** td / sigma_td ** - time-delay chi2 disabled."
   end

   # Parity penalty strength (default defined here, override with parity_force in the source section)
   parity_force = Float64(get(source_dict, :parity_force, 1000.0))

   # Finite-difference step (arcsec) used in the flux chi2 magnification correction
   pixel_fd = Float64(get(source_dict, :pixel_fd, observation.pixel_scale))

   return SourceConfig(sources        = sources, 
                       use_parity     = use_parity, 
                       parity_force   = parity_force,
                       use_flux       = use_flux, 
                       use_time_delay = use_td, 
                       pixel_fd       = pixel_fd)
end


function _source_from_file!(dict::Dict, cosmo::Cosmology.AbstractCosmology, observation::Observation, params::Vector{Parameter})
   # Get source dictionary from the full dictionary
   source_dict = dict[:source]

   # Get file name
   file_name = source_dict[:source_file]
   file_data = readdlm(file_name, comments=true, comment_char='#')
   file_data = Float64.(file_data)

   # --- Parameters directly read from the YAML file -----------------------------------------------
   # Get total number of sources from first column in file
   n_source = Int64(maximum(file_data[:, 1]))

   # Get reference lens redshift
   z_d = dict[:observation][:z_d]

   # Number of columns in the source file
   # Column 1 : source id, 
   # Column 2 : knot id
   # Column 3 : x
   # Column 4 : y
   # Column 5 : z_s
   # Column 6 : sigma_x
   # Column 7 : sigma_y
   # Column 8 : sigma_theta
   # Column 9 : (optional) parity (parity: +1 / -1, with 0 marking images without parity information)
   # Column 10: (optional) magnitude
   # Column 11: (optional) sigma_mag (sigma_mag <= 0 --> no flux data for knot)
   # Column 12: (optional) td
   # Column 13: (optional) sigma_td (sigma_td  <= 0 --> no time-delay data for knot)
   n_col = size(file_data, 2)

   # Track whether any knot carries parity / flux / time-delay measurements
   any_parity = false
   any_flux   = false
   any_td     = false

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
         adis_r, adis_l, adis_u = _zs_to_adis(cosmo, z_d, r, l, u)

         # Add redshift parameter to the reference, lower, and upper vectors
         push!(params, Parameter(owner = source_id, 
                                 name  = Symbol(:adis, i), 
                                 refer = adis_r, 
                                 lower = adis_l, 
                                 upper = adis_u))
      else
         # Get the source redshift from the file; treat it as fixed (lower == upper == refer)
         z_s = source_data[1, 5]
         adis_r, adis_l, adis_u = _zs_to_adis(cosmo, z_d, z_s, z_s, z_s)

         # Add redshift parameter to the reference, lower, and upper vectors
         push!(params, Parameter(owner = source_id, 
                                 name  = Symbol(:adis, i), 
                                 refer = adis_r, 
                                 lower = adis_l, 
                                 upper = adis_u))
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

         # Read optional parity values from column 9: +1 / -1, with 0 marking images without
         # parity information. A single measured image is already informative (absolute quantity).
         p = zeros(Int64, length(x))
         if n_col >= 9
            p_temp = knot_data[:, 9]
            if !all(v -> v in (-1.0, 0.0, 1.0), p_temp)
               error("** parity ** values (column 9) must be +1, -1, or 0 (no information) in source-$i, knot-$k")
            end
            p = Int64.(p_temp)
         end
         if any(!=(0), p)
            any_parity = true
         end
         
         # Read optional flux measurements (magnitudes). sigma <= 0 marks unmeasured images;
         # at least two measured images are required (relative measurement).
         m  = Float64[]
         σm = Float64[]
         if n_col >= 11
            σm_temp = knot_data[:, 11]
            n_meas  = count(>(0.0), σm_temp)
            if n_meas >= 2
               m  = knot_data[:, 10]
               σm = σm_temp
               any_flux = true
            elseif n_meas == 1
               @warn "Only one measured image for ** mag ** in source-$i, knot-$k - ignored (needs >= 2)."
            end
         end

         # Read optional time-delay measurements (days). Same rules as flux above.
         td   = Float64[]
         σ_td = Float64[]
         if n_col >= 13
            σtd_temp = knot_data[:, 13]
            n_meas   = count(>(0.0), σtd_temp)
            if n_meas >= 2
               td   = knot_data[:, 12]
               σ_td = σtd_temp
               any_td = true
            elseif n_meas == 1
               @warn "Only one measured image for ** td ** in source-$i, knot-$k - ignored (needs >= 2)."
            end
         end

         # Convert to arcsec if reference is not (0, 0)
         x_arcsec, y_arcsec = _to_arcsec(observation, x, y)
         knots[k] = Knot(x      = x_arcsec, 
                         σx     = σx, 
                         y      = y_arcsec, 
                         σy     = σy, 
                         σθ     = deg2rad.(σθ), 
                         parity = p,
                         m      = m, 
                         σm     = σm, 
                         td     = td, 
                         σ_td   = σ_td)
      end
      sources[i] = Source(knots=knots)
   end
   # Enable parity / flux / time-delay terms if data is present, unless explicitly disabled in the input.
   # Parity
   if get(source_dict, :use_parity, false)
      use_parity = any_parity
   else
      use_parity = false
   end

   # Flux
   if get(source_dict, :use_flux, false)
      use_flux = any_flux
   else
      use_flux = false
   end

   # Time Delay
   if get(source_dict, :use_time_delay, false)
      use_time_delay = any_td
   else
      use_time_delay = false
   end

   if get(source_dict, :use_parity, false) === true && !any_parity
      @warn "use_parity is set but no knot provides ** parity ** values - parity term disabled."
   end

   if get(source_dict, :use_flux, false) === true && !any_flux
      @warn "use_flux is set but no knot provides ** mag / sigma_mag ** - flux chi2 disabled."
   end

   if get(source_dict, :use_time_delay, false) === true && !any_td
      @warn "use_time_delay is set but no knot provides ** td / sigma_td ** - time-delay chi2 disabled."
   end

   # Parity penalty strength (default defined here, override with parity_force in the source section)
   parity_force = Float64(get(source_dict, :parity_force, 1000.0))

   # Finite-difference step (arcsec) used in the flux chi2 magnification correction
   pixel_fd = Float64(get(source_dict, :pixel_fd, observation.pixel_scale))

   return SourceConfig(sources        = sources, 
                       use_parity     = use_parity, 
                       parity_force   = parity_force,
                       use_flux       = use_flux, 
                       use_time_delay = use_time_delay, 
                       pixel_fd       = pixel_fd)
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

# Infer method
function _infer_method(dict::Dict, meta_keys::Set{Symbol})
   found = nothing
   for k in keys(dict)
      k in meta_keys && continue
      if found !== nothing 
         error("Multiple modeling methods found. Please clarify.")
      end
      found = k
   end
   if found === nothing
      error("No modeling method found.")
   end
   return found
end

function _optimizer!(sampling_dict::Dict)
   optimizer = get(sampling_dict, :optimizer, nothing)

   # Check if optimizer is missing or explicitly disabled
   if optimizer === nothing || !get(optimizer, :enabled, true)
      return nothing
   end

   # Infer optimizer method from keys
   method = _infer_method(optimizer, OPTIMIZER_META_KEYS)

   # Method section may be empty (e.g. "NM:" with nothing below it) --> use all defaults
   config = optimizer[method]
   if config === nothing 
      config = Dict{Symbol,Any}()
   end

   algorithm_config =
      if method == :NM
         NMConfig(
            max_iter  = Int64(get(config, :max_iter, 10000)),
            tolerance = Float64(get(config, :tolerance, 1e-6))
         )
      elseif method == :GD
         GDConfig(
            max_iter      = Int64(get(config, :max_iter, 10000)),
            learning_rate = Float64(get(config, :learning_rate, 0.1)),
            tolerance     = Float64(get(config, :tolerance, 1e-6))
         )
      else
         error("Unknown optimizer method: $method")
      end
   
   # Convergence parameters (defaults defined here)
   convergence = get(optimizer, :convergence, Dict{Symbol,Any}())
   if convergence === nothing
      convergence = Dict{Symbol,Any}()
   end

   if haskey(convergence, :run_mode)
      convergence[:run_mode] = Symbol(convergence[:run_mode])
   end

   run_mode  = Symbol(get(convergence, :run_mode, :random))
   max_runs  = Int64(get(convergence, :max_runs, 100))
   tolerance = Float64(get(convergence, :tolerance, 1e-3))

   return OptimizerConfig(method    = method, 
                          run_mode  = run_mode,
                          max_runs  = max_runs,
                          tolerance = tolerance,
                          config    = algorithm_config)
end

# --------------------------------------------------------------------------------------------------
# ---------------- Read MCMC -----------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
const MCMC_META_KEYS = Set([:enabled])
function _mcmc!(sampling_dict::Dict)
   mcmc = get(sampling_dict, :mcmc, nothing)

   # Check if mcmc is missing or explicitly disabled
   if mcmc === nothing || !get(mcmc, :enabled, true)
      return nothing
   end

   # Infer mcmc method from keys
   method = _infer_method(mcmc, MCMC_META_KEYS)
   
   # Method section may be empty (e.g. "AIES:" with nothing below it) --> use all defaults
   config = mcmc[method]
   if config === nothing
      config = Dict{Symbol,Any}()
   end

   algorithm_config =
      if method == :MH
         n_steps = Int64(get(config, :n_steps, 10000))
         MHConfig(
            n_walkers = Int64(get(config, :n_walkers, 1)),
            n_steps   = n_steps,
            n_adapt   = Int64(get(config, :n_adapt, div(n_steps, 10)))
         )
      elseif method == :AIES
         AIESConfig(
            n_walkers = Int64(get(config, :n_walkers, 100)),
            n_steps   = Int64(get(config, :n_steps, 100000)),
            a         = Float64(get(config, :a, 2.0))
         )
      else
         error("Unknown MCMC method: $method")
      end
   
   return MCMCConfig(method = method, config = algorithm_config)
end

# Internal function: Process Lens Model section
const SAMPLING_SCHEMES = Set([:SourcePlane, :ImagePlane])
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
      scheme    = scheme,
      verbose   = verbose,
      optimizer = optimizer_params,
      mcmc      = mcmc_params
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

   # Get cosmology and its parameters
   cosmology = _cosmology!(dict, params)

   # Get basic observation details
   observation = _observation(dict, cosmology)

   # Get lens model and its parameters
   lens_config = _lensmodel!(dict, params, observation)

   # Get source model and its parameters
   source_config = _source!(dict, cosmology, observation, params)

   # Get sampling details
   sampler = _sampling!(dict)

   # Identify free parameters
   free_param_idxs = findall(p -> p.lower != p.upper, params)

   return ModelConfig(
      observation     = observation,
      cosmology       = cosmology,
      lens_config     = lens_config,
      source_config   = source_config,
      parameters      = params,
      free_param_idxs = free_param_idxs,
      sampler         = sampler
   )
end

end