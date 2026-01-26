module MH


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using Random
using StatsBase
using .Threads
using ProgressMeter


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export mh_runner
export report_rates
export get_diagnostics


function _mh_runner_adaptive(log_posterior, start_θ::Vector{Float64}, n_steps::Int64, n_adapt::Int64, initial_σ::Vector{Float64}, p::Progress)
   # Number of free parameters
   n_params = length(start_θ)
   
   # Pre-allocate memory for the full chain
   full_chain = zeros(Float64, n_steps, n_params)
    
   # Current state
   current_θ = copy(start_θ)
   proposal_θ = similar(current_θ)
   current_logp = log_posterior(current_θ)
   σ = copy(initial_σ)
    
   # Target acceptance rate
   target_rate = 0.234
   target_window = (0.23, 0.24)
   block_size = n_adapt
   num_blocks = div(n_steps, block_size)
   rate_history = zeros(Float64, num_blocks)

   # Adaptation lock variables
   adaptation_locked = false
   stable_blocks = 0
   required_stable_blocks = 2

   for b in 1:num_blocks
      block_accepted = 0

      for i in 1:block_size
         # Calculate absolute index in the full chain
         idx = (b-1)*block_size + i
            
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

         next!(p)
      end
        
      # --- Adaptation Block ---
      # Update σ based on the performance of the last 10,000 steps
      current_rate = block_accepted / block_size
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
   return full_chain, rate_history
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


function get_diagnostics(chains::Array{Float64, 3}; param_names=nothing, burn_in::Float64=0.3)
   # 1. Map dimensions to the [Param, Chain, Step] convention
   n_steps, n_chains, n_params = size(chains)
   
   # Calculate burn-in offset for the 3rd dimension
   start_idx = max(1, Int(floor(n_steps * burn_in)) + 1)
   n_steps_post = n_steps - start_idx + 1
   total_samples = n_steps_post * n_chains
   
   # Table UI setup
   println("\n" * "-"^77)
   header = "| " * rpad("Owner", 14) * 
            "| " * rpad("Parameter", 16) * 
            "| " * rpad("Tau (τ)", 12) * 
            "| " * rpad("ESS", 12) * 
            "| " * rpad("ESS %", 10) * "  |"
   println(header)
   println("-" * "─"^75 * "-")
    
   for i in 1:n_params
      # 2. Extract and Flatten: [Chain, Step] -> 1D Vector
      # We loop through chains and then steps to preserve temporal order per chain
      @views flat_data = vec(chains[start_idx:end, :, i])
      total_samples = length(flat_data)
      
      # Identify parameter labels
      owner = param_names !== nothing ? string(param_names[i][1]) : "Unknown"
      p_name = param_names !== nothing ? string(param_names[i][2]) : "theta_$i"
      
      # 3. Statistical Calculations
      # Limit lags to prevent excessive computation on 5M steps
      max_lag = min(total_samples ÷ 5, 10000)
      ac = StatsBase.autocor(flat_data, 0:max_lag)
      
      # Integrated Autocorrelation Time (Tau)
      # Find where autocorrelation drops below noise level (0.05)
      idx = findfirst(val -> val < 0.05, ac)
      stop_at = isnothing(idx) ? length(ac) : idx
      tau = 1 + 2 * sum(@view ac[2:stop_at])
      
      # Effective Sample Size
      ess = total_samples / tau
      ess_per = (ess / total_samples) * 100
      
      # 4. Print Row
      row = "| " * rpad(owner, 14) * 
            "| " * rpad(p_name, 16) * 
            "| " * rpad(string(round(tau, digits=1)), 12) * 
            "| " * rpad(string(round(Int, ess)), 12) * 
            "| " * rpad(string(round(ess_per, digits=3)) * "%", 10) * 
            "  |"
      println(row)
   end
   println("-"^77 * "\n")
end


function mh_runner(log_posterior::Function, seeds::Vector{Vector{Float64}}, n_steps::Int64, n_adapt::Int64)
   # Number of chains
   n_chains = length(seeds)
   
   # Number of parameters
   n_params = length(seeds[1])
   
   # Number of blocks
   num_blocks = div(n_steps, n_adapt)
    
   # Initial heuristic for σ: 2% of starting value or a small floor
   # This acts as the starting point for the adaptive blocks.
   σ_initial = max.(abs.(seeds[1]) .* 0.02, 1e-4)

   # Pre-allocate 3D Tensor: [Iteration, Parameter, Chain]
   all_chains = zeros(Float64, n_steps, n_chains, n_params)
   all_rates = zeros(Float64, num_blocks, n_chains)

   # Run Metropolis-Hastings sampler
   p = Progress(n_steps * n_chains; dt=0.1, desc="Sampling Posterior... ", barlen=50)
   @threads for c in 1:n_chains
      # Each thread runs its own independent adaptive sampler
      chains, rate_history = _mh_runner_adaptive(log_posterior, seeds[c], n_steps, n_adapt, σ_initial, p)
      @views all_chains[:, c, :] .= chains
      all_rates[:, c] .= rate_history
   end

   # Check rate stability
   report_rates(all_rates)

   return all_chains, all_rates
end

end