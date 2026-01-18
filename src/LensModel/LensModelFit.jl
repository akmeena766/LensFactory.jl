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


# # --------------------------------------------------------------------------------------------------
# # Log-likelihood
# # --------------------------------------------------------------------------------------------------
# function log_likelihood(model::ModelConfig, θ::Vector{Float64})

# end


# # --------------------------------------------------------------------------------------------------
# # Log-prior
# # --------------------------------------------------------------------------------------------------
# function log_prior(model::ModelConfig, θ::Vector{Float64})
#    for (x, p) in zip(θ, free_parameters(model))
#       if x < p.lower || x > p.upper
#          return -Inf
#       end
#    end
#    return 0.0
# end


# --------------------------------------------------------------------------------------------------
# Run Optimizer
# --------------------------------------------------------------------------------------------------
function run_optimizer(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol},Float64})
   opt = model.sampler.optimizer
   opt === nothing && return nothing

   # Initial vectors in free-parameter space
   θ0 = LensModelUtils.θ_initializer(model)
   println(θ0)
   θ_best = nothing
   ll_best = -Inf
   return θ_best, ll_best   
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
   end
   
   # # MCMC
   # if sampler.mcmc !== nothing
   #    return run_mcmc(model, θ_start)
   # end

   # # Nothing to do
   # error("Neither optimizer nor MCMC enabled.")
end

end