module MH


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using Random
using StatsBase
using Base.Threads
using ProgressMeter


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export mh_runner
export report_rates


# --------------------------------------------------------------------------------------------------
# Functions
# --------------------------------------------------------------------------------------------------
function _mh_runner_adaptive(log_posterior, start_θ::Vector{Float64}, n_steps::Int64, n_adapt::Int64, initial_σ::Vector{Float64}, p::Progress)
   # Number of free parameters
   n_params = length(start_θ)
   
   # Pre-allocate memory for the full chain
   full_chain = zeros(Float64, n_steps, n_params)
   
   # Pre-allocate memory for the log-likelihood history
   ll_history = zeros(Float64, n_steps)
    
   # Current state
   current_θ = copy(start_θ)
   proposal_θ = similar(current_θ)
   current_logp = log_posterior(current_θ)
   σ = copy(initial_σ)
    
   # Target acceptance rate
   target_rate   = 0.234
   target_window = (0.23, 0.24)
   block_size    = n_adapt
   num_blocks    = cld(n_steps, block_size)
   rate_history  = zeros(Float64, num_blocks)

   # Adaptation lock variables
   adaptation_locked = false
   stable_blocks = 0
   required_stable_blocks = 2

   for b in 1:num_blocks
      block_accepted = 0

      # Number of steps in this block (last block may be shorter)
      i0 = (b - 1) * block_size
      n_block = min(block_size, n_steps - i0)

      for i in 1:n_block
         # Calculate absolute index in the full chain
         idx = i0 + i
            
         # Proposal Step: In-place to avoid allocations
         for j in 1:n_params
            proposal_θ[j] = current_θ[j] + randn() * σ[j]
         end
            
         # Likelihood Evaluation
         proposal_logp = log_posterior(proposal_θ)

         # Metropolis-Hastings Acceptance Criterion
         # Uses log-space to prevent underflow
         if log(rand()) < (proposal_logp - current_logp)
            current_θ .= proposal_θ
            current_logp = proposal_logp
            block_accepted += 1
         end
            
         # Record State
         @views full_chain[idx, :] .= current_θ
         ll_history[idx] = current_logp

         next!(p)
      end
        
      # --- Adaptation Block ---
      # Update σ based acceptance rate over the last block
      current_rate = block_accepted / n_block
      rate_history[b] = current_rate
        
      if !adaptation_locked
         # Check if we are in the "Goldilocks" zone (23% - 24%)
         if target_window[1] <= current_rate <= target_window[2]
            stable_blocks += 1
         else
            # Outside window: Update σ and reset stability counter
            γ = 0.5 / sqrt(b)
            σ .*= exp.(γ .* (current_rate .- target_rate))
            stable_blocks = 0 
         end
        
         # Lock σ if we've been stable long enough
         if stable_blocks >= required_stable_blocks
            adaptation_locked = true
         end
      end
   end
   return full_chain, ll_history, rate_history
end


function report_rates(all_rates::Matrix{Float64}, target::Float64=0.234)
   num_blocks, n_chains = size(all_rates)
    
   # Focus on the last 20% of the adaptation
   tail_idx = Int(floor(0.8 * num_blocks))
    
   header = "| " * rpad("Chain", 6) * "| " * rpad("Final Mean", 12) * "| " * rpad("Status", 10) * " |"
   line = "-"^length(header)
   
   println(line)
   println(header)
   println(line)

   for c in 1:n_chains
      final_rates = @view all_rates[tail_idx:end, c]
      mu = round(StatsBase.mean(final_rates), digits=3)
        
      status = abs(mu - target) < 0.05 ? "Tuned" : "Drifting"
        
      row = "| " * rpad(c, 6) * "| " * rpad(mu, 12) * "| " * rpad(status, 10) * " |"
      println(row)
   end
   println(line)
end


function mh_runner(log_posterior::Function, seeds::Vector{Vector{Float64}}, n_steps::Int64, n_adapt::Int64, verbose::Bool)
   # Number of chains
   n_chains = length(seeds)
   
   # Number of parameters
   n_params = length(seeds[1])
   
   # Number of blocks
   num_blocks = cld(n_steps, n_adapt)
    
   # Initial heuristic for σ: 2% of starting value or a small floor
   # This acts as the starting point for the adaptive blocks.
   σ_initial = max.(abs.(seeds[1]) .* 0.02, 1e-4)

   # Pre-allocate 3D Tensor: [Iteration, Parameter, Chain]
   all_chains = zeros(Float64, n_steps, n_chains, n_params)
   all_rates = zeros(Float64, num_blocks, n_chains)
   all_llhoods = zeros(Float64, n_steps, n_chains) # NEW: Likelihood tensor

   # Run Metropolis-Hastings sampler
   p = Progress(n_steps * n_chains; dt=0.1, desc="Sampling Posterior... ", barlen=50)
   @threads for c in 1:n_chains
      # Each thread runs its own independent adaptive sampler
      chains, ll_history, rate_history = _mh_runner_adaptive(log_posterior, seeds[c], n_steps, n_adapt, σ_initial, p)
      @views all_chains[:, c, :] .= chains
      all_rates[:, c]   .= rate_history
      all_llhoods[:, c] .= ll_history
   end

   # Check rate stability
   if verbose
      report_rates(all_rates)
   end

   return all_chains, all_llhoods
end

end