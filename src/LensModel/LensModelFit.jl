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
# Build lens model from physical paramerters
# --------------------------------------------------------------------------------------------------
function build_lens(model::ModelConfig, pvals::Dict{Tuple{Symbol,Symbol}, Float64})
   # Determine the number of components from the lens model container
   n_lens = length(model.lens_config.components)

   # Initialize an empty vector to store lens parameters
   lens_vector = NamedTuple[]

   # Iterate over each lens component
   components = model.lens_config.components

   for i in 1:n_lens
      lens_id = Symbol(:lens, i)
      lens_params = Dict{Symbol, Union{Symbol, Float64}}()
      
      for (k, v) in enumerate(components)
         if v.owner == lens_id
            lens_params[:lens] = v.name
         end
      end

      for (k, v) in pvals
         if k[1] == lens_id
            lens_params[k[2]] = v
         end
      end
      push!(lens_vector, (; lens_params...))
   end
   return Lenses.init_CompositeLens(lens_vector)
end


# --------------------------------------------------------------------------------------------------
# Log-likelihood
# --------------------------------------------------------------------------------------------------
function log_likelihood(model::ModelConfig, θ::Vector{Float64}, param_ref::Dict{Tuple{Symbol,Symbol}, Float64})
   # 1. Merge θ (free parameters) with param_ref (fixed parameters)
   # param_dict is a utility that creates a full parameter mapping 
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
   pos_ll = Likelihood.LogL_position(model, adis, αx_all, αy_all)

   
   return 0.0
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
   return - (lp + ll)
end


# --------------------------------------------------------------------------------------------------
# Run Optimizer
# --------------------------------------------------------------------------------------------------
function run_optimizer(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol}, Float64}, opt::OptimizerConfig, cfg::NMConfig, verbose::Bool)
   θ_initial = θ_initializer(model)

   best_θ   = nothing
   best_val = -Inf

   # Store results to check for convergence
   converged_results = []

   if verbose
      println("Running Nelder-Mead optimizer...")
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
      if all(abs.(result.θ .- best_θ) .< opt.tolerance) && isapprox(result.f, best_val, atol=opt.tolerance)
         same_best_count += 1
      end
   end
   
   # Print statistics
   if verbose
      w = 12
      println("Optimization Summary: ", rpad("Total Runs", w), "| ", rpad("Converged", w), "| ", rpad("Same Best", w), "| ", rpad("Best LogL", w))
      println("Values:               ", rpad(total_runs, w), "| ", rpad(total_converged, w), "| ", rpad(same_best_count, w), "| ", rpad(total_converged > 0 ? round(best_val, digits=4) : "N/A", w))
   end

   return best_θ, best_val
end

function run_optimizer(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol},Float64}, verbose::Bool)
   opt = model.sampler.optimizer

   opt === nothing && return nothing

   # Initial vectors in free-parameter space
   θ_initial = LensModelUtils.θ_initializer(model)

   return run_optimizer(model, param_ref, opt, opt.config, verbose)
end


# # --------------------------------------------------------------------------------------------------
# # Run MCMC
# # --------------------------------------------------------------------------------------------------
# function run_mcmc(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol},Float64}, )
#    if θ_start === nothing
#       # No optimizer → initialize internally
#       θ_start = θ_reference(model)
#       # or: θ_random(model)
#    end
# end


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
      θ_start, _ = run_optimizer(model, param_ref, verbose)
   else
      θ_start = θ_initializer(model)[1]
   end
   
   # # MCMC
   # if sampler.mcmc !== nothing
   #    return run_mcmc(model, θ_start)
   # end

   # # Nothing to do
   # error("Neither optimizer nor MCMC enabled.")
end

end