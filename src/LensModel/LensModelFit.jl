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

include("./NelderMead.jl")
using .NelderMead


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


# # --------------------------------------------------------------------------------------------------
# # Build lens model from physical paramerters
# # --------------------------------------------------------------------------------------------------
# function build_lens(model::ModelConfig, pvals::Dict{Tuple{Symbol,Symbol}, Float64})
#    n_lens = length(model.lens_model.lenses)
#    lens_vector = NamedTuple[]

#    for i in 1:n_lens
#       lens_id = Symbol(:lens, i)
#       lens_params = Dict{Symbol, Union{Symbol, Float64}}()
      
#       for (k, v) in pvals
#          if k[1] == lens_id
#             lens_params[k[2]] = v
#          end
#       end
      
#    end
# end


# --------------------------------------------------------------------------------------------------
# Log-likelihood
# --------------------------------------------------------------------------------------------------
function log_likelihood(model::ModelConfig, θ::Vector{Float64})

end


# --------------------------------------------------------------------------------------------------
# Log-prior
# --------------------------------------------------------------------------------------------------
function log_prior(model::ModelConfig, θ::Vector{Float64})
   for (x, p) in zip(θ, free_parameters(model))
      if x < p.lower || x > p.upper
         return -Inf
      end
   end
   return 0.0
end


# --------------------------------------------------------------------------------------------------
# Posterior (may return -Inf)
# --------------------------------------------------------------------------------------------------
function log_posterior(model::ModelConfig, θ::Vector{Float64})
    log_prior(model, θ) + log_likelihood(model, θ)
end


# --------------------------------------------------------------------------------------------------
# Optimizer-safe objective (finite value)
# --------------------------------------------------------------------------------------------------
function objective(model::ModelConfig, θ::Vector{Float64})
    lp = log_prior(model, θ)
    lp == -Inf && return -1e300
    return lp + log_likelihood(model, θ)
end


# --------------------------------------------------------------------------------------------------
# Run Optimizer
# --------------------------------------------------------------------------------------------------
function run_optimizer(model::ModelConfig, opt::OptimizerConfig, cfg::NMConfig)
   θ_initial = θ_initializer(model)

   best_θ   = nothing
   best_val = -Inf

   # Store results to check for convergence
   converged_results = []

   println("Running optimizer...")
   for θ0 in θ_initial
      θ = copy(θ0)

      θ_opt, fmax, _, _, converged = nmsmax(x -> objective(model, x), θ; tol = cfg.tolerance, max_its = cfg.max_iter)

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
   println("Optimization Summary: Total Runs | Converged | Same Best | Best LogL")
   println("Values:               $total_runs | 
                                 $total_converged | 
                                 $same_best_count | 
                                 $(total_converged > 0 ? round(best_val, digits=4) : "N/A")")

   return best_θ, best_val
end

function run_optimizer(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol},Float64})
   opt = model.sampler.optimizer

   opt === nothing && return nothing

   # Initial vectors in free-parameter space
   θ_initial = LensModelUtils.θ_initializer(model)

   return run_optimizer(model, opt, opt.config)
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
   
   # Freeze reference once
   param_ref = Dict(p.key => p.refer for p in model.parameters)

   # Initial guess
   θ_start = nothing

   # Optimization
   if sampler.optimizer !== nothing
      θ_start, _ = run_optimizer(model, param_ref)
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