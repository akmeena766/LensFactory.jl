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
export autocorrelation
export acceptance_rate

# """
#     aies_runner(log_posterior, seeds, n_steps, a, rng) -> (chain, llhood)
# Affine-Invariant Ensemble Sampler (Goodman & Weare 2010).
 
# # Arguments
# - `log_posterior::Function`: accepts a `Vector{Float64}` and returns a scalar log-posterior.
# - `seeds::Vector{Vector{Float64}}`: initial walker positions. Must satisfy:
#     - `length(seeds)` is even
#     - `length(seeds) >= 2 × n_params`
# - `n_steps::Integer`: number of MCMC steps *per walker* to return.

# # Keyword Arguments
# - `a::Number = 2.0`: stretch move scale parameter. The default rarely needs changing.
# - `rng::AbstractRNG`: master RNG used for the sequential acceptance draws and for seeding
#     the per-thread proposal RNGs. Pass a seeded `Xoshiro(seed)` for reproducibility.
 
# # Returns
# - `chain::Array{Float64, 3}`: sampled positions, shape `(n_steps, n_walkers, n_params)`.
#    Row 1 contains the seed positions; rows 2:n_steps are MCMC samples.
#    The seed row is discarded naturally when burn-in is applied downstream.
# - `llhood::Matrix{Float64}`: log-likelihood at each sample, shape `(n_steps, n_walkers)`.
 
# # Notes
# - Stretch move draws `z ~ g(z) ∝ 1/√z` on `[1/a, a]` via exact CDF inversion:
#     `z = ((a − 1)·u + 1)² / a`, `u ~ Uniform(0, 1)`.
# - Proposals for each batch are generated in parallel using `@threads`.  Each thread
#   holds its own `Xoshiro` RNG (independently seeded from `rng`) to avoid data races.
# - The MH acceptance step is kept sequential to preserve detailed balance.
# """
function aies_runner(log_posterior::Function, 
                     seeds::Vector{Vector{Float64}}, 
                     n_steps::Integer;
                     a::Number               = 2.0, 
                     rng::Random.AbstractRNG = Random.default_rng())
   # Determine dimensions from seeds
   n_walkers = length(seeds)
   n_params  = length(seeds[1])


   # -----------------------------------------------------------------------------------------------
   # Strict Validation Checks
   # -----------------------------------------------------------------------------------------------
   if isodd(n_walkers)
      throw(ArgumentError("n_walkers must be even for the AIES batch update (currently $n_walkers)."))
   end

   if n_walkers < 2 * n_params
      throw(ArgumentError("n_walkers ($n_walkers) < 2 x n_params ($(2*n_params))."))
   end
   

   # -----------------------------------------------------------------------------------------------
   # Per-walker RNGs
   # Each walker gets its own RNG independently seeded from the master rng.
   # Indexed by k (1:n_active) inside @threads — avoids threadid() which is
   # unstable under Julia >=1.9's task scheduler and can exceed nthreads(),
   # causing a BoundsError on a per-thread RNG vector.
   # -----------------------------------------------------------------------------------------------
   walker_rngs = [Xoshiro(rand(rng, UInt64)) for _ in 1:n_walkers]

   
   # -----------------------------------------------------------------------------------------------
   # Pre-allocate output: [Steps, Chains, Params]
   # Row 1 holds the seed positions; rows 2:n_steps hold the MCMC samples.
   # The seed row is intentionally kept — it will be discarded as part of burn-in
   # by the caller (get_best_fit_parameters applies a burn-in fraction).
   # -----------------------------------------------------------------------------------------------
   _chain  = Array{Float64}(undef, n_steps, n_walkers, n_params)
   _llhood = Array{Float64}(undef, n_steps, n_walkers)


   # -----------------------------------------------------------------------------------------------
   # Initialize current state from seeds
   # Column-per-walker: (n_params x n_walkers) for cheap column slicing
   # -----------------------------------------------------------------------------------------------
   x = hcat(seeds...)

   # Initial Likelihood Evaluations
   last_llhoods = Array{Float64}(undef, n_walkers)
   @threads for w in 1:n_walkers
      last_llhoods[w] = log_posterior(x[:, w])
   end

   # Store initial state in first row
   for w in 1:n_walkers
      _chain[1, w, :] .= x[:, w]
      _llhood[1, w]    = last_llhoods[w]
   end


   # -----------------------------------------------------------------------------------------------
   # Batch setup: walkers are split into two complementary sets.
   # Each step, one set generates proposals stretched against the other.
   # -----------------------------------------------------------------------------------------------
   half      = div(n_walkers, 2)
   batch1    = 1:half
   batch2    = (half + 1):n_walkers
   divisions = [(batch1, batch2), (batch2, batch1)]

   # Pre-allocate per-batch workspace once (avoids repeated allocation in the hot loop).
   # Both batches have the same size (n_walkers is even), so n_half covers both.
   n_half        = half
   new_positions = Matrix{Float64}(undef, n_params, n_half)
   new_ll        = Vector{Float64}(undef, n_half)
   zs            = Vector{Float64}(undef, n_half)


   # -----------------------------------------------------------------------------------------------
   # Main sampling loop
   # Progress bar tracks (n_steps - 1) × n_walkers walker-steps.
   # Row 1 is seed-filling only — no ticks for that.
   # -----------------------------------------------------------------------------------------------
   p = Progress((n_steps - 1) * n_walkers; dt=0.1, desc="Sampling Posterior (AIES)... ")

   for i in 2:n_steps
      for (active, inactive) in divisions
         n_active = length(active)

         # -- Parallel proposal generation ---------------------------------------------------------
         # Using thread local RNGs         
         @threads for k in 1:n_active
            walker = active[k]
            
            # Stretch move: sample Z ~ g(z) ∝ 1/√z on [1/a, a]
            U     = rand(walker_rngs[k])
            z     = ((a - 1) * U + 1)^2 / a
            zs[k] = z

            # Pick random partner from the complemetary (inactive) set
            j = rand(walker_rngs[k], inactive)

            # Proposal: Y = X_j + Z(X_i - X_j)
            @views new_positions[:, k] .= x[:, j] + z * (x[:, walker] - x[:, j])
            new_ll[k] = log_posterior(new_positions[:, k])
         end

         # -- Sequential acceptence ----------------------------------------------------------------
         # The master rng is used here; no race condition since this loop is sequential.
         for k in 1:n_active
            walker = active[k]

            # Goodman & Weare acceptance ratio
            log_ratio = (n_params - 1) * log(zs[k]) + new_ll[k] - last_llhoods[walker]

            if log(rand(rng)) < log_ratio
               @views x[:, walker] .= new_positions[:, k]
               last_llhoods[walker] = new_ll[k]
            end

            # Record results using the correct chain index
            @views _chain[i, walker, :] .= x[:, walker]
            _llhood[i, walker]           = last_llhoods[walker]
            
            # Update the progress bar
            next!(p)
         end
      end
   end 
   return _chain, _llhood
end

end