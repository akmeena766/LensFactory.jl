module LensModelFit

using ..LensModel

function log_likelihood(model::ModelConfig, θ::Vector{Float64})

end


function log_prior(model::ModelConfig, θ::Vector{Float64})
   for (x, p) in zip(θ, free_parameters(model))
      if x < p.lower || x > p.upper
         return -Inf
      end
   end
   return 0.0
end

function run_optimizer(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol},Float64})
   opt = model.sampler.optimizer
   opt === nothing && return nothing

   # --------------------------------------------------
   # Initial point in free-parameter space
   # --------------------------------------------------
   θ0 = θ_initializer(model)

      
end

function run_mcmc(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol},Float64}, )
   
end

function fit_model(model::ModelConfig)
   # Extract sampler
   sampler = model.sampler
   
   # --------------------------------------------------
   # Freeze reference once
   # --------------------------------------------------
   param_ref = Dict(p.key => p.refer for p in model.parameters)


   # --------------------------------------------------
   # Optional optimization
   # --------------------------------------------------
   if sampler.optimizer !== nothing
      θ_start, _ = run_optimizer(model, param_ref)
   end

   # --------------------------------------------------
   # MCMC
   # --------------------------------------------------
   if sampler.mcmc !== nothing
      return run_mcmc(model, θ_start)
   end

   # --------------------------------------------------
   # Nothing to do
   # --------------------------------------------------
   error("Neither optimizer nor MCMC enabled.")
end

end