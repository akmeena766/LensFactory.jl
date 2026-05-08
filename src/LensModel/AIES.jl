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
   # Per-thread RNGs
   # Each thread gets its own RNG independently seeded from the master rng.
   # This avoids data races in @threads blocks while keeping the run reproducible
   # when a seeded master rng is supplied.
   # -----------------------------------------------------------------------------------------------
   nt   = nthreads()
   rngs = [Xoshiro(rand(rng, UInt64)) for _ in 1:nt]

   
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
      next!(p)
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
   half      = div(n_chains, 2)
   batch1    = 1:half
   batch2    = (half + 1):n_chains
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
            tid    = threadid()
            walker = active[k]
            
            # Stretch move: sample Z ~ g(z) ∝ 1/√z on [1/a, a]
            U     = rand(rngs[tid])
            z     = ((a - 1) * U + 1)^2 / a
            zs[k] = z

            # Pick random partner from the complemetary (inactive) set
            j = rand(rngs[tid], inactive)

            # Proposal: Y = X_j + Z(X_i - X_j)
            @views new_positions[:, k] .= x[:, j] + z * (x[:, chain_idx] - x[:, j])
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


function autocorrelation(chains::Array{Float64, 3}; 
                         param_names :: Union{Vector{Tuple{Symbol, Symbol}}, Nothing} = nothing, 
                         burn_in     :: Float64                                       = 0.2)
   # Dimension
   n_steps, n_walkers, n_params = size(chains)

   # Calculate burn-in offset
   start_idx     = max(1, Int(floor(n_steps * burn_in)) + 1)
   n_steps_post  = n_steps - start_idx + 1
   total_samples = n_steps_post * n_walkers

   # -- Printing the UI ----------------------------------------------------------------------------
   println("\n" * "-"^77)
   @printf("| %-14s | %-16s | %-12s | %-12s | %-10s  |\n", "Owner", "Parameter", "Tau (τ)", "ESS", "ESS %")
   println("-"^77 * "\n")

   for i in 1:n_params
      tau_total = 0.0
      
      # Calculate Tau per chain to avoid artificial "jumps" from flattening
      for w in 1:n_walkers
         @views chain_data = chains[start_idx:end, w, i]
         
         # Limit lags to N/5 — the standard recommendation beyond which
         # autocorrelation estimates become unreliable due to sample noise.
         max_lag = min(length(chain_data) ÷ 5, 2000)
         ac      = autocor(chain_data, 0:max_lag)
         
         # Integrated Autocorrelation Time (IAT or τ)
         # Truncate at the first negative lag to avoid noise inflation.
         # Note: for very noisy chains this may truncate early; a full
         # automated-windowing estimator (Sokal 1989) can be added later.
         idx     = findfirst(val -> val <= 0.0, ac)
         stop_at = isnothing(idx) ? length(ac) : idx
         
         # Tau formula: 1 + 2 * sum(autocorrelations)
         tau_total += 1.0 + 2.0 * sum(@view ac[2:stop_at])
      end
      
      avg_tau = tau_total / n_walkers
      ess     = total_samples / avg_tau
      ess_per = (ess / total_samples) * 100
      
      # Identify parameter labels
      owner  = (param_names !== nothing && i <= length(param_names)) ? string(param_names[i][1]) : "Lens"
      p_name = (param_names !== nothing && i <= length(param_names)) ? string(param_names[i][2]) : "theta_$i"
      
      # 4. Print Row
      @printf("| %-14s | %-16s | %-12.1f | %-12d | %-9.2f%% | \n", owner, p_name, avg_tau, round(Int, ess), ess_per)
   end
   println("-"^77 * "\n")
end


function acceptance_rate(chains::Array{Float64, 3}; burn_in::Float64=0.2)
   # Dimension
   n_steps, n_walkers, _ = size(chains)
   
   # Calculate burn-in offset
   start_idx = max(1, Int(floor(n_steps * burn_in)) + 1)

   # Calculate per-walker acceptance
   # A rejection is an exact bitwise copy of the previous row, so checking
   # whether any parameter changed is unambiguous and handles all edge cases.
   walker_rates = zeros(n_walkers)
   for w in 1:n_walkers
      accepted = 0
      for s in (start_idx + 1):n_steps
         @views if any(chains[s, w, :] .!= chains[s-1, w, :])
               accepted += 1
         end
      end
      walker_rates[w] = (accepted / (n_steps - start_idx)) * 100
   end

   # Summary Statistics
   avg_acc = StatsBase.mean(walker_rates)
   min_acc = minimum(walker_rates)
   max_acc = maximum(walker_rates)

   # --Printing the UI -----------------------------------------------------------------------------
   println("\n" * "─"^60)
   println(" ACCEPTANCE RATE DIAGNOSTICS")
   println("─"^60)

   @printf("%-22s %6.2f%%\n", " Overall Average Rate:", avg_acc)
   @printf("%-22s %6.2f%%\n", " Lowest Chain Rate:",    min_acc)
   @printf("%-22s %6.2f%%\n", " Highest Chain Rate:",   max_acc)
   println("─"^60)

   # Visual Sparkline
   print(" Ensemble Spread: [")
   blocks = [" ", "▂", "▃", "▄", "▅", "▆", "▇", "█"]
   for rate in chain_rates
      # Map rate 0-100% to block index 1-8
      b_idx = clamp(Int(ceil(rate / 12.5)), 1, 8)
      print(blocks[b_idx])
   end
   println("]")

   # Guidance Logic
   println("─"^60)
   @printf(" Status: ")
   if avg_acc < 15.0
      println("⚠️  LOW (Stiff). Consider decreasing stretch parameter 'a'.")
   elseif avg_acc > 50.0
      println("⚠️  HIGH (Baby steps). Consider increasing stretch parameter 'a'.")
   else
      println("✅ HEALTHY. The sampler is mixing well.")
   end
   println("─"^60 * "\n")
end

end