# --------------------------------------------------------------------------------------------------
# Sampling - Optimizer/MCMC/sampler structs and readers for the sampling section.
# This file is include-d by LensModelIO.jl and lives inside the LensModelIO module.
# --------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export NMConfig
export OptimizerConfig
export MHConfig
export AIESConfig
export SamplerConfig


# --------------------------------------------------------------------------------------------------
# Abstract type: Optimizer
# --------------------------------------------------------------------------------------------------
@kwdef struct NMConfig <: AbstractOptimizerConfig
   max_iter::Int64
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
# ---------------- Read Optimizer ------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------
const OPTIMIZER_META_KEYS = Set([:enabled, :convergence])

# Build all multi-plane distance (ratio) matrices - one-time cost at input reading, since the
# cosmology, the lens plane redshifts, and the source redshifts are all fixed
function _build_multiplane_distances(cosmo::Cosmology.AbstractCosmology, z_lenses::Vector{Float64}, zs_all::Vector{Float64})
   # Sorted unique plane redshifts (must match the grouping in Lenses.init_MultiPlaneLens)
   z_planes = sort(unique(z_lenses))
   n_p = length(z_planes)
   n_s = length(zs_all)

   # Plane-plane distance ratios (source independent)
   adis_ij = zeros(n_p, n_p)
   for ni in 1:n_p
      for nj in 1:ni-1
         adis_ij[nj, ni] = Cosmology.angular_diameter_distance(cosmo, z_planes[nj], z_planes[ni]) /
                           Cosmology.angular_diameter_distance(cosmo, 0.0, z_planes[ni])
      end
   end

   # Per-source quantities
   adis_is = Vector{Vector{Float64}}(undef, n_s)
   D_ij_td = Vector{Matrix{Float64}}(undef, n_s)
   for si in 1:n_s
      zs = zs_all[si]
      D_os = Cosmology.angular_diameter_distance(cosmo, 0.0, zs)

      # Plane-source distance ratios (planes at or behind the source do not deflect it)
      v = zeros(n_p)
      for ni in 1:n_p
         if z_planes[ni] < zs
            v[ni] = Cosmology.angular_diameter_distance(cosmo, z_planes[ni], zs) / D_os
         end
      end
      adis_is[si] = v

      # Distance matrix truncated to the planes in front of this source (time-delay core
      # infers the number of planes from the matrix size)
      m = count(z -> z < zs, z_planes)
      z_all = [0.0; z_planes[1:m]; zs]
      D = zeros(m + 2, m + 2)
      for i in 1:m+2
         for j in i+1:m+2
            D[i, j] = Cosmology.angular_diameter_distance(cosmo, z_all[i], z_all[j])
         end
      end
      D_ij_td[si] = D
   end

   return MultiPlaneDistances(z_planes = z_planes, 
                              adis_ij  = adis_ij, 
                              adis_is  = adis_is, 
                              D_ij_td  = D_ij_td)
end

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
