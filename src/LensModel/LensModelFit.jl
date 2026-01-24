module LensModelFit


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..Constants
using ..Lenses

using ..LensModelIO
using ..LensModelUtils
using ..Likelihood

using ..NelderMead

using ..MH


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export fit_model


# --------------------------------------------------------------------------------------------------
# Parameter-space → Physical-space transformation
# --------------------------------------------------------------------------------------------------
const PARAM_TRANSFORM = Dict{Symbol,Function}(
   # Velocity dispersion: km/s -> m/s
   :v_d => x -> x * 1.0E3,
   
   # Mass:Log10(M/M☉) -> kg
   :m => x -> 10^x * MASS_SUN,
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
# Log-likelihood
# --------------------------------------------------------------------------------------------------
function log_likelihood(model::ModelConfig, θ::Vector{Float64}, param_ref::Dict{Tuple{Symbol,Symbol}, Float64})
   # Merge θ (free parameters) with param_ref (fixed parameters)
   pvals = LensModelUtils.param_dict(model, θ, param_ref)

   # Transform parameters (from sample space to physical space)
   transform_params!(pvals)

   # Build lens model
   lens_model = LensModelUtils.build_lens(model, pvals)

   # Get angular-diameter distance ratios
   adis = LensModelUtils.adis_current(model, pvals)

   # Calculate deflection at image positions
   ψ_all, αx_all, αy_all, A_all, P_all = LensModelUtils.lens_quantities(model, lens_model)

   # Calculate position likelihood
   pos_ll = Likelihood.LogL_position(model, adis, αx_all, αy_all, A_all)

   
   return pos_ll
end


# --------------------------------------------------------------------------------------------------
# Log-prior
# --------------------------------------------------------------------------------------------------
function log_prior(model::ModelConfig, θ::Vector{Float64})
   @inbounds for (x, p) in zip(θ, free_parameters(model))
      (x < p.lower || x > p.upper) && return -Inf
   end
   return 0.0
end


# --------------------------------------------------------------------------------------------------
# Posterior (may return -Inf)
# --------------------------------------------------------------------------------------------------
function log_posterior(model::ModelConfig, θ::Vector{Float64}, param_ref::Dict{Tuple{Symbol,Symbol}, Float64})
   # Calculate log-prior (returns -Inf if any parameter is out of bounds)
   lp = log_prior(model, θ)
   if lp == -Inf
      return -Inf
   end

   # Calculate log-likelihood
   ll = log_likelihood(model, θ, param_ref)
   
   # Return log-posterior
   return lp + ll
end


# --------------------------------------------------------------------------------------------------
# Optimizer-safe objective (finite value)
# --------------------------------------------------------------------------------------------------
function objective(model::ModelConfig, θ::Vector{Float64}, param_ref::Dict{Tuple{Symbol,Symbol}, Float64})
   # Calculate log-prior (returns a large negative value if any parameter is out of bounds)
   lp = log_prior(model, θ)
   if lp == -Inf
      return -1e300
   end

   # Calculate log-likelihood
   ll = log_likelihood(model, θ, param_ref)
   
   # Return negative log-posterior
   return lp + ll
end


# --------------------------------------------------------------------------------------------------
# Run Optimizer
# --------------------------------------------------------------------------------------------------
function run_optimizer(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol}, Float64}, opt::OptimizerConfig, cfg::NMConfig, verbose::Bool)
   θ_initial = θ_initializer(model)

   best_θ   = nothing
   best_val = -Inf

   # Store results to check for convergence
   converged_results = Vector{NamedTuple{(:θ, :f), Tuple{Vector{Float64}, Float64}}}(undef, 0)

   if verbose
      println("\nRunning Nelder-Mead optimizer...")
   end
   for θ0 in θ_initial
      θ = copy(θ0)

      # Call optimizer
      θ_opt, fmax, _, _, converged = nmsmax(x -> objective(model, x, param_ref), θ; tol = cfg.tolerance, max_its = cfg.max_iter)

      # Process only if the simplex size reached tolerance
      if converged
         push!(converged_results, (θ=copy(θ_opt), f=fmax))

         # Update best parameters 
         if fmax > best_val
            best_val = fmax
            best_θ   = copy(θ_opt)
         end
      end
   end

   # Statistics and consistency checks
   total_runs = opt.max_runs
   total_converged = length(converged_results)
   same_best_count = 0
   
   for result in converged_results
      # Calculate relative difference for parameters (avoiding division by zero)
      # Using 1e-8 as a floor to handle parameters that are exactly 0.0
      rel_diff_θ = abs.(result.θ .- best_θ) ./ (max.(abs.(best_θ), abs.(result.θ)) .+ 1e-8)
    
      # Calculate relative difference for the LogL (objective value)
      rel_diff_f = abs(result.f - best_val) / (abs(best_val) + 1e-8)
      # println(rel_diff_θ, " ", rel_diff_f)
      # Check if all parameters are within percentage tolerance (e.g., 0.1% = 0.001)
      if all(rel_diff_θ .< opt.tolerance) && (rel_diff_f < opt.tolerance)
         same_best_count += 1
      end
   end

   # Calculate the percentage
   conv_rate = (total_converged / total_runs) * 100
   stability_rate = (same_best_count / total_converged) * 100

   # Print statistics
   if verbose
      w = 12
      # Expanded table to include the stability percentage
      col_h = 
         "| " * rpad("Total", 8) * 
         "| " * rpad("Convergence (%)", 12) * 
         "| " * rpad("Stability (%)", 12) * 
         "| " * rpad("Best LogL", w) * 
         "|"
      
      val_r = 
         "| " * rpad(total_runs, 8) * 
         "| " * rpad("$(round(conv_rate, digits=1))%", 15) * 
         "| " * rpad("$(round(stability_rate, digits=1))%", 13) * 
         "| " * rpad(round(best_val, digits=4), w) * 
         "|"
   
      line = "-" ^ length(col_h)
      println(line); println(col_h); println(line); println(val_r); println(line)
   end
   return converged_results
end

# Main optimizer function (wrapper to use multiple dispatch)
function run_optimizer(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol},Float64}, verbose::Bool)
   opt = model.sampler.optimizer

   opt === nothing && return nothing

   # Initial vectors in free-parameter space
   θ_initial = LensModelUtils.θ_initializer(model)

   return run_optimizer(model, param_ref, opt, opt.config, verbose)
end


# --------------------------------------------------------------------------------------------------
# Run MCMC
# --------------------------------------------------------------------------------------------------
function get_unique_seeds(results::Vector{@NamedTuple{θ::Vector{Float64}, f::Float64}}, n_chains::Int64, tol::Float64=1E-2)
   # Sort by LogL (highest/best first)
   # results should be a vector of structs/objects with .θ and .f (LogL)
   sorted_res = sort(results, by = x -> x.f, rev = true)

   # Start with the absolute best
   seeds = [sorted_res[1].θ]

   # Add unique seeds
   for i in eachindex(sorted_res)[2:end]
      # Get current result
      current_θ = sorted_res[i].θ

      # Check if this result is far enough away from existing seeds
      is_distinct = all(seeds) do s
         rel_diff = abs.(current_θ .- s) ./ (max.(abs.(s), abs.(current_θ)) .+ 1e-8)
         return any(rel_diff .> tol) 
      end
        
      if is_distinct
         push!(seeds, current_θ)
      end

      # Exit if we have reached the required number of chains
      length(seeds) >= n_chains && break
   end

   # Jitter fill if we are short on unique peaks
   while length(seeds) < n_chains
      push!(seeds, seeds[1] .+ (abs.(seeds[1]) .* 1e-3 .+ 1e-5) .* randn(length(seeds[1])))
   end

   # Return unique seeds
   return seeds
end

function get_best_seeds(results::Vector{@NamedTuple{θ::Vector{Float64}, f::Float64}}, n_chains::Int64; jitter::Float64=0.001)
   # Sort by LogL (highest/best first)
   # results should be a vector of structs/objects with .θ and .f (LogL)
   sorted_res = sort(results, by = x -> x.f, rev = true)
 
   n_params = length(sorted_res[1].θ)
    
    # Create seeds by adding a tiny amount of Gaussian noise to the best solution.
    # We use relative jitter (jitter * best_theta) so that parameters with 
    # large absolute values (like v_d ≈ 250) get scaled appropriately.
    seeds = [
        sorted_res[1].θ .+ (jitter .* abs.(sorted_res[1].θ) .* randn(n_params)) 
        for _ in 1:n_chains
    ]
    
    return seeds
end

function run_mcmc(model::ModelConfig, mcmc_config::MHConfig, param_ref::Dict{Tuple{Symbol,Symbol},Float64}, θ_start::Vector{Vector{Float64}}, verbose::Bool)
   return MH.mh_runner(x -> log_posterior(model, x, param_ref), θ_start, mcmc_config.n_steps, mcmc_config.n_adapt)
end

function run_mcmc(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol},Float64}, θ_start::Union{Nothing, Vector{@NamedTuple{θ::Vector{Float64}, f::Float64}}}, verbose::Bool)   
   if verbose
      println("\nRunning MCMC...")
   end

   # Get MCMC config for multiple dispatch
   mcmc = model.sampler.mcmc

   # Initialize seeds (No optimizer → random initial parameters)
   seeds = if θ_start !== nothing
      get_best_seeds(θ_start, mcmc.n_chains)
   else
      [LensModelUtils.θ_random(model) for _ in 1:mcmc.n_chains]
   end

   return run_mcmc(model, mcmc.config, param_ref, seeds, verbose)
end


# --------------------------------------------------------------------------------------------------
# Fit model
# --------------------------------------------------------------------------------------------------
function fit_model(model::ModelConfig)
   # Extract sampler
   sampler = model.sampler
   verbose = sampler.verbose
   
   # Freeze reference once
   param_ref = Dict(p.key => p.refer for p in model.parameters)

   # Initial guess
   θ_start = nothing

   # Optimization
   if sampler.optimizer !== nothing
      θ_start = run_optimizer(model, param_ref, verbose)
   end

   # MCMC
   if sampler.mcmc !== nothing
      return run_mcmc(model, param_ref, θ_start, verbose)
   end

   # Nothing to do
   error("Neither optimizer nor MCMC enabled.")
end

end