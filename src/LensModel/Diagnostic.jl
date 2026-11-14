module Diagnostic


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to use
# --------------------------------------------------------------------------------------------------
using StatsBase
using Statistics
using Printf


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..LensModelUtils: _burn_in_start


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export calculate_gr
export print_gr_report
export autocorrelation
export acceptance_rate


# --------------------------------------------------------------------------------------------------
# Gelman-Rubin diagnostic
# --------------------------------------------------------------------------------------------------
function calculate_gr(chains::Array{Float64, 3}; burn_in::Float64=0.2)
   # Dimensions
   n_steps, n_chains, n_params = size(chains)
   
   # Check number of chains
   if n_chains < 2
      error("At least 2 chains are required for Gelman-Rubin diagnostic.")
   end

   # Remove burn-in
   start_idx = _burn_in_start(n_steps, burn_in)
   samples   = chains[start_idx:end, :, :]

   # Post burn in size
   n     = size(samples, 1)
   m     = n_chains
   R_hat = zeros(Float64, n_params)

   for p in 1:n_params
      @views param_samples = samples[:, :, p]
        
      # W: Within-chain variance
      chain_vars = Statistics.var(param_samples, dims=1)
      W = Statistics.mean(chain_vars)
        
      # B: Between-chain variance
      chain_means = Statistics.mean(param_samples, dims=1)
      grand_mean = Statistics.mean(chain_means)
      B = (n / (m - 1)) * sum((chain_means .- grand_mean).^2)
        
      # V_hat: Pooled variance estimate
      V_hat = ((n - 1) / n) * W + (1 / n) * B
      R_hat[p] = sqrt(V_hat / W)
   end
   return R_hat
end

function print_gr_report(chains::Array{Float64, 3}; free_param_names=nothing, burn_in=0.2)
   # Calculate R̂
   r_hats = calculate_gr(chains, burn_in=burn_in)
   n_params = length(r_hats)

   # Table Formatting
   w_owner = 15
   w_param = 15
   w_rhat  = 10
   w_stat  = 12

   # Build Header
   header = "| " * rpad("Owner", w_owner) * " | " * rpad("Parameter", w_param) * " | " * rpad("R-hat", w_rhat) * " | " * rpad("Status", w_stat) * " |"
   border = "-"^length(header)

   total_inner_width = w_owner + w_param + w_rhat + w_stat + (3 * 3)
   title = "GELMAN-RUBIN CONVERGENCE DIAGNOSTIC"
   padding = total_inner_width - length(title)
   left_pad = div(padding, 2)
   right_pad = padding - left_pad + 2
   centered_title = " "^left_pad * title * " "^right_pad

   # -- Printing the UI ----------------------------------------------------------------------------
   println("\n" * border)
   println("|" * centered_title * "|")
   println(border)
   println(header)
   println(border)

   for i in 1:n_params
        # Extract metadata from the parameter object
        if free_param_names !== nothing
            owner = string(free_param_names[i][1])
            param = string(free_param_names[i][2])
        else
            owner = "Unknown"
            param = "theta_$i"
        end

        val   = round(r_hats[i], digits=4)
        
        status = val < 1.1 ? "Converged" : "FAILED"
        
        # Build Row
        row = "| " * rpad(owner, w_owner) * " | " * rpad(param, w_param) * " | " * rpad(val, w_rhat)   * " | " * rpad(status, w_stat) * " |"
        println(row)
    end
    println(border)
   
   # Global Summary
   inner_width = length(header) - 4
   if all(filter(!isnan, r_hats) .< 1.1)
        println("| " * rpad("✅ SUCCESS: Global convergence reached.", inner_width) * " |")
    else
        println("| " * rpad("⚠️  WARNING: High variance detected between chains.", inner_width) * " |")
    end
   println(border * "\n")
end


# --------------------------------------------------------------------------------------------------
# Time series diagnostics
# --------------------------------------------------------------------------------------------------
# τ = 1 + 2 Σ_{k=1}^{K} ρ_k, where ρ_k = Σ_t z_t z_{t+k} / Σ_t z_t² (z = x - mean(x)) and K is the
# last lag before ρ first becomes ≤ 0 (or `max_lag` if it never does). Lags are computed one at a
# time and the loop stops at the first non-positive ρ_k, so a well-mixed chain only needs a few lags.
function _integrated_autocorr_time(x::AbstractVector{Float64}, max_lag::Int)
   n = length(x)
   μ = sum(x) / n
   z = x .- μ
   c0 = 0.0
   @inbounds @simd for t in 1:n
      c0 += z[t] * z[t]
   end
   c0 == 0.0 && return 1.0

   acc = 0.0
   for k in 1:min(max_lag, n - 1)
      ck = 0.0
      @inbounds @simd for t in 1:(n - k)
         ck += z[t] * z[t + k]
      end
      ρ = ck / c0
      ρ <= 0.0 && break
      acc += ρ
   end
   return 1.0 + 2.0 * acc
end


# Compute and print the Integrated Autocorrelation Time (IAT / τ) and Effective Sample Size (ESS) 
# for each free parameter across all walkers.
# - IAT is estimated per walker and averaged, avoiding artificial decorrelation from flattening 
#   walkers into a single chain.
# - Autocorrelation is truncated at the first negative lag (Sokal windowing, simplified). This is 
#   conservative — a full automated-windowing estimator can be substituted later for very long chains.
# - Lags are capped at `min(N ÷ 5, 2000)` where `N` is the post-burn-in length per walker — the 
#   standard recommendation beyond which lag estimates become noise-dominated.
# - A healthy ESS% is typically > 1–2%. Very low values (< 0.5%) suggest the chain needs more steps 
#   or the sampler is poorly mixed.
function autocorrelation(chains::Array{Float64, 3}; 
                         free_param_names::Union{Vector{Tuple{Symbol, Symbol}}, Nothing} = nothing, 
                         burn_in::Float64                                                = 0.2)
   # Dimension
   n_steps, n_walkers, n_params = size(chains)

   # Calculate burn-in offset
   start_idx         = _burn_in_start(n_steps, burn_in)
   post_burnin_steps = n_steps - start_idx + 1
   total_samples     = post_burnin_steps * n_walkers

   # Store IAT and ESS
   avg_taus = zeros(Float64, n_params)
   ess_vals = zeros(Float64, n_params)

   # -- Per-walker IAT, computed in parallel -------------------------------------------------------
   # Every (walker, parameter) pair is independent; each iteration writes only to its own slot.
   taus = Matrix{Float64}(undef, n_walkers, n_params)
   @threads for k in 1:(n_walkers * n_params)
      w, i = fldmod1(k, n_walkers)
      @views chain_data = chains[start_idx:end, w, i]

      # Limit lags to N/5 — beyond this, estimates are noise-dominated.
      max_lag = min(length(chain_data) ÷ 5, 2000)
      taus[w, i] = _integrated_autocorr_time(chain_data, max_lag)
   end



   # -- Printing the UI ----------------------------------------------------------------------------
   println("\n" * "-"^78)
   @printf("| %-14s | %-16s | %-12s | %-12s | %-8s |\n", "Owner", "Parameter", "Tau (τ)", "ESS", "ESS [%]")
   println("-"^78)

   for i in 1:n_params
      avg_tau = sum(@view taus[:, i]) / n_walkers
      ess     = total_samples / avg_tau
      ess_per = (ess / total_samples) * 100
      
      # Identify parameter labels
      owner  = (free_param_names !== nothing && i <= length(free_param_names)) ? string(free_param_names[i][1]) : "Lens"
      p_name = (free_param_names !== nothing && i <= length(free_param_names)) ? string(free_param_names[i][2]) : "theta_$i"
      
      # Print Row
      ess_str = @sprintf("%.2f %%", ess_per)
      @printf("| %-14s | %-16s | %-12.1f | %-12d | %-8s | \n", owner, p_name, avg_tau, round(Int, ess), ess_str)
   end
   println("-"^78 * "\n")
   return nothing
end


# Compute and print the acceptance rate for each walker, along with an ensemble sparkline and 
# sampler health guidance. Acceptance is detected by checking whether any parameter changed between 
# consecutive steps. A rejected proposal is an exact bitwise copy of the previous row, so 
# `any(chains[s, w, :] .!= chains[s-1, w, :])` is unambiguous. A healthy AIES acceptance rates 
# typically fall in the range `[20%, 40%]` for well-conditioned posteriors.
function acceptance_rate(chains::Array{Float64, 3}; burn_in::Float64=0.2)
   # Dimension
   n_steps, n_walkers, _ = size(chains)
   
   # Calculate burn-in offset
   start_idx = _burn_in_start(n_steps, burn_in)

   # Calculate per-walker acceptance
   # A rejection is an exact bitwise copy of the previous row, so checking
   # whether any parameter changed is unambiguous and handles all edge cases.
   walker_rates = zeros(n_walkers)
   for w in 1:n_walkers
      accepted = 0
      for s in (start_idx + 1):n_steps
         if @views any(chains[s, w, :] .!= chains[s-1, w, :])
               accepted += 1
         end
      end
      walker_rates[w] = (accepted / (n_steps - start_idx)) * 100
   end

   # Summary Statistics
   avg_acc = Statistics.mean(walker_rates)
   min_acc = minimum(walker_rates)
   max_acc = maximum(walker_rates)

   # --Printing the UI -----------------------------------------------------------------------------
   title = "ACCEPTANCE RATE DIAGNOSTICS"
   padding = 78 - length(title)
   left_pad = div(padding, 2)
   right_pad = padding - left_pad + 2
   centered_title = " "^left_pad * title * " "^right_pad
   
   println("\n" * "-"^78)
   println(centered_title)
   println("-"^78)

   @printf("%-22s %6.2f%%\n", " Overall Average Rate:", avg_acc)
   @printf("%-22s %6.2f%%\n", " Lowest Chain Rate:",    min_acc)
   @printf("%-22s %6.2f%%\n", " Highest Chain Rate:",   max_acc)
   println("-"^78)

   # Visual Sparkline
   print(" Ensemble Spread: [")
   blocks = [" ", "▂", "▃", "▄", "▅", "▆", "▇", "█"]
   for rate in walker_rates
      # Map rate 0-100% to block index 1-8
      b_idx = clamp(Int(ceil(rate / 12.5)), 1, 8)
      print(blocks[b_idx])
   end
   println("]")

   # Guidance Logic
   println("-"^78)
   @printf(" Status: ")
   if avg_acc < 20.0
      println("⚠️  LOW (Stiff). Consider decreasing stretch parameter 'a'.")
   elseif avg_acc > 50.0
      println("⚠️  HIGH (Baby steps). Consider increasing stretch parameter 'a'.")
   else
      println("✅ HEALTHY. The sampler is mixing well.")
   end
   println("-"^78 * "\n")
   return nothing
end

end