module Diagnostic


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to use
# --------------------------------------------------------------------------------------------------
using StatsBase


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
function time_series_diagnostics(chains::Array{Float64, 3}; param_names=nothing, burn_in::Float64=0.2)
   # 1. Map dimensions: [Step, Chain, Param]
   n_steps, n_chains, n_params = size(chains)
   
   # Calculate burn-in offset
   start_idx = max(1, Int(floor(n_steps * burn_in)) + 1)
   n_steps_post = n_steps - start_idx + 1
   
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
      tau_total = 0.0
      
      # Calculate Tau per chain to avoid artificial "jumps" from flattening
      for c in 1:n_chains
         @views chain_data = chains[start_idx:end, c, i]
         
         # 2. Autocorrelation using StatsBase
         # Limit lags to 2000; if it hasn't decayed by then, the chain is stuck.
         max_lag = min(length(chain_data) ÷ 5, 2000)
         ac = autocor(chain_data, 0:max_lag)
         
         # 3. Integrated Autocorrelation Time (Tau)
         # Sum until the autocorrelation becomes negative or noise-dominated
         idx = findfirst(val -> val <= 0.0, ac)
         stop_at = isnothing(idx) ? length(ac) : idx
         
         # Tau formula: 1 + 2 * sum(autocorrelations)
         tau_total += 1.0 + 2.0 * sum(@view ac[2:stop_at])
      end
      
      avg_tau = tau_total / n_chains
      total_samples = n_steps_post * n_chains
      ess = total_samples / avg_tau
      ess_per = (ess / total_samples) * 100
      
      # Identify parameter labels
      owner = (param_names !== nothing && i <= length(param_names)) ? string(param_names[i][1]) : "Lens"
      p_name = (param_names !== nothing && i <= length(param_names)) ? string(param_names[i][2]) : "theta_$i"
      
      # 4. Print Row
      row = "| " * rpad(owner, 14) * 
            "| " * rpad(p_name, 16) * 
            "| " * rpad(string(round(avg_tau, digits=1)), 12) * 
            "| " * rpad(string(round(Int, ess)), 12) * 
            "| " * rpad(string(round(ess_per, digits=2)) * "%", 10) * "  |"
      println(row)
   end
   println("-"^77 * "\n")
end


# --------------------------------------------------------------------------------------------------
# Acceptance rate diagnostics
# --------------------------------------------------------------------------------------------------
function acceptance_diagnostics(chains::Array{Float64, 3}; burn_in::Float64=0.2)
   n_steps, n_chains, _ = size(chains)
   start_idx = max(1, Int(floor(n_steps * burn_in)) + 1)
   
   # Calculate per-chain acceptance
   chain_rates = zeros(n_chains)
   for c in 1:n_chains
      accepted = 0
      for s in (start_idx + 1):n_steps
         @views if chains[s, c, 1] != chains[s-1, c, 1]
               accepted += 1
         end
      end
      chain_rates[c] = (accepted / (n_steps - start_idx)) * 100
   end

   # Summary Statistics
   avg_acc = StatsBase.mean(chain_rates)
   min_acc = minimum(chain_rates)
   max_acc = maximum(chain_rates)
    
   # Printing the UI
   println("\n" * "─"^60)
   println(" ACCEPTANCE RATE DIAGNOSTICS")
   println("─"^60)
   
   println(" Overall Average Rate:  ", lpad(round(avg_acc, digits=2), 6), "%")
   println(" Lowest Chain Rate:     ", lpad(round(min_acc, digits=2), 6), "%")
   println(" Highest Chain Rate:    ", lpad(round(max_acc, digits=2), 6), "%")
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
   print(" Status: ")
   if avg_acc < 15.0
      println("⚠️  LOW (Stiff). Consider decreasing stretch parameter 'a'.")
   elseif avg_acc > 60.0
      println("⚠️  HIGH (Baby steps). Consider increasing stretch parameter 'a'.")
   else
      println("✅ HEALTHY. The sampler is mixing well.")
   end
   println("─"^60 * "\n")
end

end