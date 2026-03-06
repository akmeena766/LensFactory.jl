module LensModelUtils

# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using Random
using StatsBase


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..Constants
using ..Cosmology
using ..Lenses
using ..LFUtils
using ..LensModelIO
using ..Likelihood

# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export free_parameters
export θ_initializer
export param_dict
export adis_current
export build_lens
export lens_quantities
export free_parameter_names
export get_best_parameters
export get_best_parameters_with_errors
export check_parity
export get_best_fit_rms
export save_best_fit


# --------------------------------------------------------------------------------------------------
# Parameter-space → Physical-space transformation
# --------------------------------------------------------------------------------------------------
const PARAM_TRANSFORM = Dict{Symbol,Function}(
   # Velocity dispersion: km/s -> m/s
   :v_d => x -> x * 1.0E3,
   
   # Mass:Log10(M/M☉) -> kg
   :mass => x -> 10^x * MASS_SUN,
)

function transform_params!(pvals::Dict{Tuple{Symbol,Symbol}, Float64})
   for key in collect(keys(pvals))
      name = key[2]
      if haskey(PARAM_TRANSFORM, name)
         pvals[key] = PARAM_TRANSFORM[name](pvals[key])
      end
   end
end


# --------------------------------------------------------------------------------------------------
# Get free parameters
# --------------------------------------------------------------------------------------------------
function free_parameters(model::ModelConfig)
   return @view model.parameters[model.free_param_idxs]
end


# --------------------------------------------------------------------------------------------------
# Get parameter values
# --------------------------------------------------------------------------------------------------
function θ_reference(model::ModelConfig)
   getfield.(free_parameters(model), :refer)
end

function θ_random(model::ModelConfig)
   θ = Vector{Float64}(undef, length(model.free_param_idxs))
   for (i, p) in enumerate(free_parameters(model))
      θ[i] = rand() * (p.upper - p.lower) + p.lower
   end
   return θ
end

function θ_grid(model::ModelConfig, max_runs::Int64)
   fps = free_parameters(model)
   nd = length(fps)
   
   # Automatically determine points per dimension
   n_per_dim = ceil(Int, max_runs^(1/nd))
   total = n_per_dim^nd

   # Generate grid
   θ_grid = Vector{Vector{Float64}}(undef, total)
   for idx in 1:total
      θ = Vector{Float64}(undef, nd)
      i = idx - 1
      for (k, p) in enumerate(fps)
         ik = i % n_per_dim
         i ÷= n_per_dim
         θ[k] = p.lower + (ik / (n_per_dim - 1)) * (p.upper - p.lower)
      end
      θ_grid[idx] = θ
   end
   return θ_grid
end

function θ_lhs(model::ModelConfig, max_runs::Int64)
   fps = free_parameters(model)
   nd = length(fps)

   lhs_matrix = zeros(max_runs, nd)
   for k in 1:nd
      # Divide [0, 1] into n intervals and pick a random point in each
      intervals = ((0:max_runs-1) .+ rand(max_runs)) ./ max_runs
      
      # shuffle for randomness
      shuffle!(intervals)
      
      # Scale to parameter bounds
      lhs_matrix[:, k] = fps[k].lower .+ intervals .* (fps[k].upper - fps[k].lower)
   end
   # Convert each row to a vector
   return [vec(lhs_matrix[i, :]) for i in 1:max_runs]
end

function θ_initializer(model::ModelConfig)
   opt = model.sampler.optimizer

   run_mode = Symbol(opt.run_mode)
   max_runs = opt.max_runs

   if run_mode === :reference
      return [θ_reference(model) for _ in 1:max_runs]
   elseif run_mode === :jitter
      return [θ_reference(model) .+ 1e-3 .* randn(length(model.free_param_idxs)) for _ in 1:max_runs]
   elseif run_mode === :random
      return [θ_random(model) for _ in 1:max_runs]
   elseif run_mode === :grid
      return θ_grid(model, max_runs)
   elseif run_mode === :lhs
      return θ_lhs(model, max_runs)
   else
      error("Unknown run_mode: $run_mode")
   end
end


# --------------------------------------------------------------------------------------------------
# Parameter dictionary (θ -> pvals)
# --------------------------------------------------------------------------------------------------
function param_dict(model::ModelConfig, θ::Vector{Float64}, param_ref::Dict{Tuple{Symbol,Symbol}, Float64})
   # Make a copy of the reference dictionary
   pvals = copy(param_ref)

   # Update free parameters
   @inbounds for (i, idx) in enumerate(model.free_param_idxs)
      p = model.parameters[idx]
      pvals[p.key] = θ[i]
   end
   return pvals
end


# --------------------------------------------------------------------------------------------------
# Angular-diameter distance (pvals -> adis)
# --------------------------------------------------------------------------------------------------
# @inline function _get_adis(pvals, adis_ref, key)
#    return get(pvals, key, adis_ref)
# end

# @inline function adis_current(model::ModelConfig, pvals)
#    return (_get_adis(pvals, model.adis_ref[i], model.parameters[i].key) for i in eachindex(model.adis_ref))
# end


function adis_current(model::ModelConfig, pvals::Dict{Tuple{Symbol,Symbol},Float64})
   nsrc = length(model.source_config.sources)
   adis = Vector{Float64}(undef, nsrc)

   @inbounds for i in 1:nsrc
      key = (Symbol(:source, i), Symbol(:adis, i))
      adis[i] = pvals[key]
   end
   return adis
end


# --------------------------------------------------------------------------------------------------
# Build lens model from physical paramerters
# --------------------------------------------------------------------------------------------------
REQUIRE_COSMO = Set([:NFWLens, :eNFWMDLens, :aNFWLens, :tNFWLens, :gNFWLens, :EinastoLens])
REQUIRE_SCALING = Set([:MultiPJELens])
function build_lens(model::ModelConfig, pvals::Dict{Tuple{Symbol,Symbol}, Float64})
   # Update scaling if it has free parameters
   scaling_params = Dict{Symbol, Any}() 
   for (k, v) in pvals
      if k[1] == :scaling
         scaling_params[k[2]] = v
      end
   end
   
   if !isempty(scaling_params)
      updated_scaling = ScalingRelation(; scaling_params...)
   else
      updated_scaling = model.lens_config.scaling
   end
   
   # Determine the number of components from the lens model container
   n_lens = length(model.lens_config.components)

   # Transform parameters
   transform_params!(pvals)

   # Initialize an empty vector to store lens parameters
   lens_vector = NamedTuple[]

   # Iterate over each lens component
   components = model.lens_config.components

   for i in 1:n_lens
      lens_id = Symbol(:lens, i)
      lens_params = Dict{Symbol, Union{Symbol, Int64, Float64, Vector{Float64}, Cosmology.AbstractCosmology}}()
      
      # Get lens name and add it to the lens_params dictionary
      for (k, v) in enumerate(components)
         if v.owner == lens_id
            lens_params[:lens] = v.name
         end
      end

      # Get other lens parameters and add it to the lens_params dictionary
      for (k, v) in pvals
         if k[1] == lens_id
            lens_params[k[2]] = v
         end
      end

      # Check if this lens model requires cosmology and lens redshift
      if lens_params[:lens] ∈ REQUIRE_COSMO
         lens_params[:cosmology] = model.cosmology
         lens_params[:z_d] = model.observation.z_d
      end

      # Check if this lens model requires galaxy component
      if lens_params[:lens] ∈ REQUIRE_SCALING
         # Total number of galaxies
         lens_params[:n] = model.lens_config.galaxies.n
         
         # Position of the galaxies (vector)
         lens_params[:x_c] = model.lens_config.galaxies.x_c
         lens_params[:y_c] = model.lens_config.galaxies.y_c

         # Ellipticity and PA of galaxies (vector)
         lens_params[:eps] = model.lens_config.galaxies.eps
         lens_params[:pa] = model.lens_config.galaxies.pa

         # Observed magnitudes (vector)
         obs_mag = model.lens_config.galaxies.obs_mag

         # Reference magnitude
         ref_mag = updated_scaling.ref_mag

         # L/L⋆ (vector)
         l_lstar = @. 10.0^(-0.4 * (obs_mag - ref_mag))

         # Velocity dispersion (vector)
         lens_params[:v_d] = @. updated_scaling.ref_sigma * 1.0E3 * l_lstar^updated_scaling.slope_sigma
         
         # Core radius (vector)
         lens_params[:x_s] = @. updated_scaling.ref_core * l_lstar^updated_scaling.slope_core

         # Truncation radius (vector)
         lens_params[:x_t] = @. updated_scaling.ref_cut * l_lstar^updated_scaling.slope_cut
      end
      push!(lens_vector, (; lens_params...))
   end
   return Lenses.init_CompositeLens(lens_vector)
end


# --------------------------------------------------------------------------------------------------
# All lensing quantities needed for log-likelihood 
# --------------------------------------------------------------------------------------------------
function lens_quantities(model::ModelConfig, lens::Lenses.AbstractLens)
   # Count the total number of knots in the lens model
   n_knots = sum(length(s.knots) for s in model.source_config.sources)

   # Allocate outputs (vector of vectors)
   ψ_all  = Vector{Vector{Float64}}(undef, n_knots)
   αx_all = Vector{Vector{Float64}}(undef, n_knots)
   αy_all = Vector{Vector{Float64}}(undef, n_knots)
   A_all  = Vector{NTuple{4, Vector{Float64}}}(undef, n_knots)

   kid = 1
   for src in model.source_config.sources
      for knot in src.knots
         # One image system knot positions
         x = knot.x
         y = knot.y

         # Potential
         ψ_all[kid] = Lenses.get_potential(lens, x, y)

         # Deflection
         αx_all[kid], αy_all[kid] = Lenses.get_deflection(lens, x, y)

         # Deformation tensor
         ψxx, ψyy, ψxy = Lenses.get_jacobian(lens, x, y)
         A_all[kid] = (ψxx, ψxy, copy(ψxy), ψyy)

         # Increment
         kid = kid + 1
      end
   end
   return ψ_all, αx_all, αy_all, A_all
end


# --------------------------------------------------------------------------------------------------
# Free parameter names
# --------------------------------------------------------------------------------------------------
function free_parameter_names(model::ModelConfig)
   return [(p.owner, p.name) for p in model.parameters[model.free_param_idxs]]
end


# --------------------------------------------------------------------------------------------------
# Get best-fit parameters from the full optimization/MCMC results
# --------------------------------------------------------------------------------------------------
function get_best_parameters(results, chains=nothing)
   # Initialize outputs
   best_θ = nothing
   best_chi2 = -Inf

   # Case A: Input is from the Parallel Optimizer (Vector of NamedTuples)
   if results isa Vector && eltype(results) <: NamedTuple
      # We already sorted the results in descending order of logL, so results[1] is the best fit
      best_run = results[1]
      best_θ = best_run.θ
      best_chi2 = best_run.f

   # Case B: Input is from MCMC (ll_history matrix)
   elseif results isa Matrix{Float64} && chains !== nothing
      # Get (step, chain) of the best logL by finding the index of the maximum value in the logL matrix
      best_idx = argmax(results) 
      step, chain_num = best_idx[1], best_idx[2]
      
      best_θ = chains[step, chain_num, :]
      best_chi2 = results[step, chain_num]
   else
      error("Please provide either Optimization results or both MCMC ll_history and chains.")
   end
   return best_θ, best_chi2
end


function get_best_parameters_with_errors(chains::Array{Float64, 3}, chi2::Matrix{Float64}; burn_in=0.2, thinning=100)
   # Get chain details 
   n_steps, n_chains, n_params = size(chains)
    
   # Extract Best-Fit (Maximum Likelihood Estimate)
   best_idx = argmax(chi2) 
   step_bf, chain_bf = best_idx[1], best_idx[2]
   
   # Best parameter and LogL values
   best_θ = chains[step_bf, chain_bf, :]
   best_chi2 = chi2[step_bf, chain_bf]

   # Remove Burn-in
   start_idx = Int(floor(n_steps * burn_in)) + 1
   thinned_chain = chains[start_idx:thinning:end, :, :]

   # Reshape the thinned chain into a flat array
   flat_chain = reshape(thinned_chain, :, n_params)
    
   # Calculate (asymmetric) errors as we are defining errors relative to the Best-Fit value
   lower_err = zeros(n_params)
   upper_err = zeros(n_params)

   for i in 1:n_params
      # Get 16th and 84th percentiles of the posterior
      q16, q84 = StatsBase.quantile(flat_chain[:, i], [0.16, 0.84])
        
      # Asymmetric error: distance from best-fit to the quantiles
      lower_err[i] = best_θ[i] - q16
      upper_err[i] = q84 - best_θ[i]
   end

   # Return as a NamedTuple for easy access
   return best_θ, lower_err, upper_err, best_chi2
end

end