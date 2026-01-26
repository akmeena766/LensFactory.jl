module AIES


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using Random
using ProgressMeter
using Base.Threads


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export aies_runner


function aies_runner(log_posterior::Function, seeds::Vector{Vector{Float64}}, n_steps::Integer, a::Number=2.0, rng::Random.AbstractRNG=Random.default_rng())
   # Determine dimensions from seeds
   n_chains = length(seeds)
   n_params = length(seeds[1])

   # --- Strict Validation Checks ---
   if isodd(n_chains)
      throw(ArgumentError("n_chains must be even for the AIES batch update (currently $n_chains)."))
   end

   if n_chains < 2 * n_params
      throw(ArgumentError("n_chains ($n_chains) < 2 * n_params ($n_params)."))
   end
   # --------------------------------
   
   # Pre-allocate output: [Steps, Chains, Params]
   chain = Array{Float64}(undef, n_steps, n_chains, n_params)
   llhood_history = Array{Float64}(undef, n_steps, n_chains)

   # Progress bar setup
   p = Progress(n_steps * n_chains; dt=0.1, desc="Sampling Posterior... ", barlen=50)

   # Matrix size: (n_params x n_chains)
   x = hcat(seeds...)

   # Initial Likelihood Evaluations
   last_llhoods = Array{Float64}(undef, n_chains)
   @threads for w in 1:n_chains
      last_llhoods[w] = log_posterior(x[:, w])
      next!(p)
   end    

   # Store initial state
   for w in 1:n_chains
      chain[1, w, :] .= x[:, w]
      llhood_history[1, w] = last_llhoods[w]
   end

   # Batch setup for AIES Parallel Update
   half = div(n_chains, 2)
   batch1 = 1:half
   batch2 = (half + 1):n_chains
   divisions = [(batch1, batch2), (batch2, batch1)]

   for i in 2:n_steps
      for (active, inactive) in divisions
         n_active = length(active)
         new_positions = Array{Float64}(undef, n_params, n_active)
         new_ll = Array{Float64}(undef, n_active)
         zs = Array{Float64}(undef, n_active)

         @threads for k in eachindex(active)
            chain_idx = active[k]
            
            # Z ~ g(z) ∝ 1/√z
            z = ((a - 1) * rand(rng) + 1)^2 / a
            zs[k] = z

            # Pick random partner from the stationary set
            j = rand(rng, inactive)

            # Y = X_j + Z(X_i - X_j)
            @views proposal = x[:, j] + z * (x[:, chain_idx] - x[:, j])

            new_positions[:, k] .= proposal
            new_ll[k] = log_posterior(proposal)
         end

         for (k, chain_idx) in enumerate(active)
            z = zs[k]
            log_ratio = (n_params - 1) * log(z) + new_ll[k] - last_llhoods[chain_idx]

            if log(rand(rng)) < log_ratio
               x[:, chain_idx] .= new_positions[:, k]
               last_llhoods[chain_idx] = new_ll[k]
            end

            # Record results using the correct chain index
            chain[i, chain_idx, :] .= x[:, chain_idx]
            llhood_history[i, chain_idx] = last_llhoods[chain_idx]
            
            next!(p)
         end
      end
   end 
   return chain, llhood_history
end

end