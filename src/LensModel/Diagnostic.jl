module Diagnostic


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to use
# --------------------------------------------------------------------------------------------------
using StatsBase
using Printf

# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export calculate_gr
export print_gr_report
export time_series_diagnostics
export acceptance_diagnostics


# --------------------------------------------------------------------------------------------------
# Optimizer convergence
# --------------------------------------------------------------------------------------------------
function optimizer_convergence(optimizer::Vector{@NamedTuple{θ::Vector{Float64}, f::Float64}}; tolerance::Float64=0.01, verbose::Bool=true)
   # Best parameter vector and chi2
   best_θ = converged_results[1].θ
   best_val = converged_results[1].f

   for result in converged_results
      # Calculate relative difference for parameters (avoiding division by zero)
      # Using 1e-8 as a floor to handle parameters that are exactly 0.0
      rel_diff_θ = abs.(result.θ .- best_θ) ./ (max.(abs.(best_θ), abs.(result.θ)) .+ 1e-8)
    
      # Calculate relative difference for the chi2 (objective value)
      rel_diff_f = abs(result.f - best_val) / (abs(best_val) + 1e-8)

      # Check if all parameters are within percentage tolerance (e.g., 0.1% = 0.001)
      if all(rel_diff_θ .< opt.tolerance) && (rel_diff_f < opt.tolerance)
         same_best_count += 1
      end
   end

   # Calculate the percentage
   conv_rate = (total_converged / total_runs) * 100
   stability_rate = (same_best_count / total_converged) * 100

   # Print statistics
   if verbose
      w = 12
      # Expanded table to include the stability percentage
      col_h = 
         "| " * rpad("Total", 8) * 
         "| " * rpad("Convergence (%)", 12) * 
         "| " * rpad("Stability (%)", 12) * 
         "| " * rpad("Best χ²", w) * 
         "|"
      
      val_r = 
         "| " * rpad(total_runs, 8) * 
         "| " * rpad("$(round(conv_rate, digits=1))%", 15) * 
         "| " * rpad("$(round(stability_rate, digits=1))%", 13) * 
         "| " * rpad(round(-best_val, digits=4), w) * 
         "|"
   
      line = "-" ^ length(col_h)
      println(line); println(col_h); println(line); println(val_r); println(line)
   end
end


# --------------------------------------------------------------------------------------------------
# Gelman-Rubin diagnostic
# --------------------------------------------------------------------------------------------------
function calculate_gr(chains::Array{Float64, 3}; burn_in::Float64=0.2)
   n_steps, n_chains, n_params = size(chains)
   
   if n_chains < 2
      error("At least 2 chains are required for Gelman-Rubin diagnostic.")
   end

   # Remove burn-in (the adaptation blocks)
   start_idx = Int(floor(burn_in * n_steps))
   samples = chains[start_idx:end, :, :]

   n = n_steps - start_idx + 1
   m = n_chains
   r_hats = zeros(Float64, n_params)

   for p in 1:n_params
      @views param_samples = samples[start_idx:end, :, p]
        
      # W: Within-chain variance
      chain_vars = StatsBase.var(param_samples, dims=1)
      W = StatsBase.mean(chain_vars)
        
      # B: Between-chain variance
      chain_means = StatsBase.mean(param_samples, dims=1)
      grand_mean = StatsBase.mean(chain_means)
      B = (n / (m - 1)) * sum((chain_means .- grand_mean).^2)
        
      # V_hat: Pooled variance estimate
      V_hat = ((n - 1) / n) * W + (1 / n) * B
      r_hats[p] = sqrt(V_hat / W)
   end
   return r_hats
end

function print_gr_report(chains::Array{Float64, 3}; param_names=nothing, burn_in=0.2)
   # Calculate R̂
   r_hats = calculate_gr(chains, burn_in=burn_in)
   n_params = length(r_hats)

   # Handle missing names
   if param_names === nothing
      # Generate names: ["θ₁", "θ₂", ...]
      display_names = ["theta_" * string(i) for i in 1:n_params]
   else
      display_names = param_names
   end

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

    println("\n" * border)
    println("|" * centered_title * "|")
    println(border)
    println(header)
    println(border)

   for i in 1:n_params
        # Extract metadata from the parameter object
        if param_names !== nothing
            owner = string(param_names[i][1])
            param = string(param_names[i][2])
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
# """
#     autocorrelation(chains; param_names, burn_in)
# Compute and print the Integrated Autocorrelation Time (IAT / τ) and Effective
# Sample Size (ESS) for each free parameter across all walkers.
 
# # Arguments
# - `chains::Array{Float64, 3}`: MCMC chain array of shape `(n_steps, n_walkers, n_params)`,
#    as returned by `aies_runner`.
 
# # Keyword Arguments
# - `param_names::Union{Vector{Tuple{Symbol, Symbol}}, Nothing} = nothing`: optional vector
#   of `(owner, param)` tuples identifying each parameter column. If `nothing`, labels default to
#   `"Lens"` / `"theta_i"`.
# - `burn_in::Float64 = 0.2`: fraction of steps to discard as burn-in before computing statistics. 
#    `start_idx = floor(n_steps * burn_in) + 1`.
 
# # Returns
# - `nothing`: Prints a formatted table with one row per parameter and the following columns:

# ```
# | Owner | Parameter | Tau (τ) | ESS | ESS % |
# ```
 
# - `Owner`     : First element of the `param_names` tuple (e.g. `:SIE`, `:src`).
# - `Parameter` : Second element of the `param_names` tuple (e.g. `:b`, `:x`).
# - `Tau (τ)`   : Average IAT across walkers. Values >> 1 indicate strong autocorrelation;
#                 τ ≈ 1 means nearly independent samples.
# - `ESS`       : Effective Sample Size = total post-burn-in samples / τ.
# - `ESS %`     : ESS as a percentage of total post-burn-in samples.
 
# # Notes
# - IAT is estimated per walker and averaged, avoiding artificial decorrelation
#   from flattening walkers into a single chain.
# - Autocorrelation is truncated at the first negative lag (Sokal windowing,
#   simplified). This is conservative — a full automated-windowing estimator
#   can be substituted later for very long chains.
# - Lags are capped at `min(N ÷ 5, 2000)` where `N` is the post-burn-in length
#   per walker — the standard recommendation beyond which lag estimates become
#   noise-dominated.
# - A healthy ESS% is typically > 1–2%. Very low values (< 0.5%) suggest the
#   chain needs more steps or the sampler is poorly mixed.
# """
function autocorrelation(chains::Array{Float64, 3}; 
                         param_names :: Union{Vector{Tuple{Symbol, Symbol}}, Nothing} = nothing, 
                         burn_in     :: Float64                                       = 0.2)
   # Dimension
   n_steps, n_walkers, n_params = size(chains)

   # Calculate burn-in offset
   start_idx     = max(1, Int(floor(n_steps * burn_in)) + 1)
   total_samples = (n_steps - start_idx + 1) * n_walkers

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
         stop_at = isnothing(idx) ? length(ac) : idx - 1
         
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


# """
#     acceptance_rate(chains; burn_in)
# Compute and print the acceptance rate for each walker, along with an ensemble
# sparkline and sampler health guidance.
 
# # Arguments
# - `chains::Array{Float64, 3}` : MCMC chain array of shape `(n_steps, n_walkers, n_params)`,
#   as returned by `aies_runner`.
 
# # Keyword Arguments
# - `burn_in::Float64 = 0.2` : fraction of steps to discard before computing rates.
#   `start_idx = floor(n_steps * burn_in) + 1`.
 
# # Returns
# `nothing`: Prints a diagnostic panel with the following sections:
 
# - `Overall Average Rate` : mean acceptance rate across all walkers (%).
# - `Lowest Walker Rate`   : minimum acceptance rate across walkers (%).
# - `Highest Walker Rate`  : maximum acceptance rate across walkers (%).
# - `Ensemble Spread`      : sparkline (`▁` to `█`) showing per-walker acceptance
#                            visually. Each block maps 12.5% of acceptance rate.
#                            `▁` indicates near-zero acceptance; `█` indicates ~100%.
# - `Status`               : health guidance based on the average rate:
#     - `< 15%`  → LOW  — chain is stiff, consider decreasing stretch parameter `a`.
#     - `> 50%`  → HIGH — steps are too small, consider increasing `a`.
#     - `15–50%` → HEALTHY — sampler is mixing well.
 
# # Notes
# - Acceptance is detected by checking whether any parameter changed between
#   consecutive steps. A rejected proposal is an exact bitwise copy of the
#   previous row, so `any(chains[s, w, :] .!= chains[s-1, w, :])` is unambiguous.
# - Healthy AIES acceptance rates typically fall in the range `[20%, 40%]` for
#   well-conditioned posteriors.
# """
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
         if @views any(chains[s, w, :] .!= chains[s-1, w, :])
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
   for rate in walker_rates
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