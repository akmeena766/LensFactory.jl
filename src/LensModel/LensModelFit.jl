module LensModelFit


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using Base.Threads
using ProgressMeter
using JLD2
using Dates


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
include("./NM.jl")
using .NM

include("./MH.jl")
using .MH

include("./AIES.jl")
using .AIES

using ..Constants
using ..Lenses
using ..LensModelIO
using ..LensModelUtils
using ..Likelihood


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export _fit_model


# --------------------------------------------------------------------------------------------------
# Log-likelihood
# --------------------------------------------------------------------------------------------------
function log_likelihood(model::ModelConfig, θ::Vector{Float64}, param_ref::Dict{Tuple{Symbol,Symbol}, <:Real})
   # Merge θ (free parameters) with param_ref (fixed parameters)
   pvals = LensModelUtils.param_dict(model, θ, param_ref)

   # Build lens model
   lens_model = LensModelUtils.build_lens(model, pvals)
   
   # Get angular-diameter distance ratios
   adis = LensModelUtils.adis_current(model, pvals)

   # Calculate position likelihood
   logL = 0.0
   if model.sampler.scheme == :SourcePlane
      # Calculate deflection at image positions
      ψ_all, αx_all, αy_all, A_all = LensModelUtils.lens_quantities(model, lens_model)

      # Calculate position likelihood
      logL_position = Likelihood.logL_sourceplane(model, adis, αx_all, αy_all, A_all)

      # Calculate parity likelihood
      if model.source_config.use_parity
         logL_parity = Likelihood.logL_sourceplane_parity(model, adis, A_all)
         logL = logL_position + logL_parity
      end      
   elseif model.sampler.scheme == :ImagePlane
      error("Image plane sampling is not implemented yet.")
      
      # Calculate deflection at image positions
      _, αx_all, αy_all, A_all = LensModelUtils.lens_quantities(model, lens_model)
      
   else
      error("Unsupported sampling scheme: $(model.sampler.scheme)")
   end
   
   return logL
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
function log_posterior(model::ModelConfig, θ::Vector{Float64}, param_ref::Dict{Tuple{Symbol,Symbol}, <:Real})
   # Calculate log-prior (returns -Inf if any parameter is out of bounds)
   lp = log_prior(model, θ)
   if lp == -Inf
      return -Inf
   end

   # Calculate log-likelihood
   logL = log_likelihood(model, θ, param_ref)
   
   # Return log-posterior
   return lp + logL
end


# --------------------------------------------------------------------------------------------------
# Optimizer-safe objective (finite value)
# --------------------------------------------------------------------------------------------------
function objective(model::ModelConfig, θ::Vector{Float64}, param_ref::Dict{Tuple{Symbol,Symbol}, <:Real})
   # Calculate log-prior (returns a large negative value if any parameter is out of bounds)
   lp = log_prior(model, θ)
   if lp == -Inf
      return -1e300
   end

   # Calculate log-likelihood
   logL = log_likelihood(model, θ, param_ref)
   
   # Return negative log-posterior
   return lp + logL
end


# --------------------------------------------------------------------------------------------------
# Run Optimizer
# --------------------------------------------------------------------------------------------------
function run_optimizer(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol}, <:Real}, opt::OptimizerConfig, cfg::NMConfig, verbose::Bool)
   # Initialize free parameter vector
   θ_initial = θ_initializer(model; run_mode=opt.run_mode, max_runs=opt.max_runs)


   # Number of runs
   n_runs = length(θ_initial)

   # Store results to check for convergence
   results = Vector{Union{Nothing, NamedTuple{(:θ, :f), Tuple{Vector{Float64}, Float64}}}}(nothing, n_runs)

   # Progress bar setup
   p = Progress(n_runs; dt=0.1, desc="Running Nelder-Mead Optimizer... ")
   
   @threads for i in 1:n_runs
      # Copying input ensures thread isolation
      θ0 = copy(θ_initial[i])

      # Call optimizer
      θ_opt, fmax, _, _, converged = nmsmax(x -> objective(model, x, param_ref), θ0; tol = cfg.tolerance, max_its = cfg.max_iter)

      # Write to the specific memory slot reserved for individual run
      if converged
         results[i] = (θ=copy(θ_opt), f=fmax)
      end
      
      # Update progress bar
      next!(p)
   end

   # Filter out the 'nothing' entries (failed convergences) and collect
   converged_results = collect(skipmissing([r === nothing ? missing : r for r in results]))
   if isempty(converged_results)
      error("No optimization runs converged.")
   end

   # Sort results (Best log-posterior first)
   sort!(converged_results, by = x -> x.f, rev = true)
   best_θ = converged_results[1].θ
   best_val = converged_results[1].f

   # Statistics and consistency checks
   total_runs = opt.max_runs
   total_converged = length(converged_results)
   same_best_count = 0
   
   for result in converged_results
      # Calculate relative difference for parameters (avoiding division by zero)
      # Using 1e-8 as a floor to handle parameters that are exactly 0.0
      rel_diff_θ = abs.(result.θ .- best_θ) ./ (max.(abs.(best_θ), abs.(result.θ)) .+ 1e-8)
    
      # Calculate relative difference for the log-posterior (objective value)
      rel_diff_f = abs(result.f - best_val) / (abs(best_val) + 1e-8)

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
         "| " * rpad("Best χ²", w) * 
         "|"
      
      val_r = 
         "| " * rpad(total_runs, 8) * 
         "| " * rpad("$(round(conv_rate, digits=1))%", 15) * 
         "| " * rpad("$(round(stability_rate, digits=1))%", 13) * 
         "| " * rpad(round(-best_val, digits=4), w) * 
         "|"
   
      line = "-" ^ length(col_h)
      println(line); println(col_h); println(line); println(val_r); println(line)
   end
   return converged_results
end

# Main optimizer function (wrapper to use multiple dispatch)
function run_optimizer(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol}, <:Real}, verbose::Bool)
   opt = model.sampler.optimizer
   return run_optimizer(model, param_ref, opt, opt.config, verbose)
end


# --------------------------------------------------------------------------------------------------
# Run MCMC
# --------------------------------------------------------------------------------------------------
function get_unique_seeds(results::Vector{@NamedTuple{θ::Vector{Float64}, f::Float64}}, n_chains::Int64, tol::Float64=1E-2)
   # Sort by log-posterior (highest/best first)
   # results should be a vector of structs/objects with .θ and .f (log-posterior)
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
   # Sort by log-posterior (highest/best first)
   # results should be a vector of structs/objects with .θ and .f (log-posterior)
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

function run_mcmc(model::ModelConfig, mcmc_config::MHConfig, param_ref::Dict{Tuple{Symbol,Symbol}, <:Real}, θ_start::Vector{Vector{Float64}}, verbose::Bool)
   return MH.mh_runner(x -> log_posterior(model, x, param_ref), θ_start, mcmc_config.n_steps, mcmc_config.n_adapt, verbose)
end

function run_mcmc(model::ModelConfig, mcmc_config::AIESConfig, param_ref::Dict{Tuple{Symbol,Symbol}, <:Real}, θ_start::Vector{Vector{Float64}}, verbose::Bool)
   return AIES.aies_runner(x -> log_posterior(model, x, param_ref), θ_start, mcmc_config.n_steps; a=mcmc_config.a)
end

function run_mcmc(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol}, <:Real}, θ_start::Union{Nothing, Vector{@NamedTuple{θ::Vector{Float64}, f::Float64}}}, verbose::Bool)   
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
function _fit_model(model::ModelConfig; save::Bool=true, file_name::Union{String, Nothing}=nothing)
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
   chains = nothing
   logL = nothing
   if sampler.mcmc !== nothing
      chains, logL = run_mcmc(model, param_ref, θ_start, verbose)
   end

   # Nothing to do
   if sampler.optimizer === nothing && sampler.mcmc === nothing
      error("Neither optimizer nor MCMC enabled.")
   end

   # Save best fit
   if save
      if file_name === nothing
         file_name = "$(model.observation.lens)_$(Dates.today()).jld2"
      end
      date = Dates.now(UTC)
      jldsave(file_name; Date      = Dates.Date(date),
                         Time      = Dates.Time(date),
                         model     = model,
                         optimizer = θ_start,
                         chains    = chains,
                         logL      = logL)
   end
   return nothing
end

end