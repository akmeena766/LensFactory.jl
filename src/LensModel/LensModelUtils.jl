module LensModelUtils

# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using Random
using StatsBase

# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..LensModelIO
using ..Lenses

# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export free_parameters
export θ_initializer
export param_dict
export adis_current
export build_lens
export lens_quantities
export free_parameter_names
export calculate_gr, print_gr_report
export time_series_diagnostics
export acceptance_diagnostics
export get_best_fit, get_best_fit_with_errors


function free_parameters(model::ModelConfig)
   return @view model.parameters[model.free_param_idxs]
end

function θ_reference(model::ModelConfig)
   getfield.(free_parameters(model), :refer)
end

function θ_random(model::ModelConfig)
   θ = Vector{Float64}(undef, length(model.free_param_idxs))
   for (i, p) in enumerate(free_parameters(model))
      θ[i] = rand() * (p.upper - p.lower) + p.lower
   end
   return θ
end

function θ_grid(model::ModelConfig, max_runs::Int64)
   fps = free_parameters(model)
   nd = length(fps)
   
   # Automatically determine points per dimension
   n_per_dim = ceil(Int, max_runs^(1/nd))
   total = n_per_dim^nd

   # Generate grid
   θ_grid = Vector{Vector{Float64}}(undef, total)
   for idx in 1:total
      θ = Vector{Float64}(undef, nd)
      i = idx - 1
      for (k, p) in enumerate(fps)
         ik = i % n_per_dim
         i ÷= n_per_dim
         θ[k] = p.lower + (ik / (n_per_dim - 1)) * (p.upper - p.lower)
      end
      θ_grid[idx] = θ
   end
   return θ_grid
end

function θ_lhs(model::ModelConfig, max_runs::Int64)
   fps = free_parameters(model)
   nd = length(fps)

   lhs_matrix = zeros(max_runs, nd)
   for k in 1:nd
      # Divide [0, 1] into n intervals and pick a random point in each
      intervals = ((0:max_runs-1) .+ rand(max_runs)) ./ max_runs
      
      # shuffle for randomness
      shuffle!(intervals)
      
      # Scale to parameter bounds
      lhs_matrix[:, k] = fps[k].lower .+ intervals .* (fps[k].upper - fps[k].lower)
   end
   # Convert each row to a vector
   return [vec(lhs_matrix[i, :]) for i in 1:max_runs]
end

function θ_initializer(model::ModelConfig)
   opt = model.sampler.optimizer

   run_mode = Symbol(opt.run_mode)
   max_runs = opt.max_runs

   if run_mode === :reference
      return [θ_reference(model) for _ in 1:max_runs]
   elseif run_mode === :jitter
      return [θ_reference(model) .+ 1e-3 .* randn(length(model.free_param_idxs)) for _ in 1:max_runs]
   elseif run_mode === :random
      return [θ_random(model) for _ in 1:max_runs]
   elseif run_mode === :grid
      return θ_grid(model, max_runs)
   elseif run_mode === :lhs
      return θ_lhs(model, max_runs)
   else
      error("Unknown run_mode: $run_mode")
   end
end


# --------------------------------------------------------------------------------------------------
# Parameter dictionary (θ -> pvals)
# --------------------------------------------------------------------------------------------------
function param_dict(model::ModelConfig, θ::Vector{Float64}, param_ref::Dict{Tuple{Symbol,Symbol}, Float64})
   # Make a copy of the reference dictionary
   pvals = copy(param_ref)

   # Update free parameters
   @inbounds for (i, idx) in enumerate(model.free_param_idxs)
      p = model.parameters[idx]
      pvals[p.key] = θ[i]
   end
   return pvals
end


# --------------------------------------------------------------------------------------------------
# Angular-diameter distance (pvals -> adis)
# --------------------------------------------------------------------------------------------------
@inline function _get_adis(pvals, adis_ref, key)
   return get(pvals, key, adis_ref)
end

@inline function adis_current(model::ModelConfig, pvals)
   return (_get_adis(pvals, model.adis_ref[i], model.parameters[i].key) for i in eachindex(model.adis_ref))
end


function adis_current(model::ModelConfig, pvals::Dict{Tuple{Symbol,Symbol},Float64})
   nsrc = length(model.source_config.sources)
   adis = Vector{Float64}(undef, nsrc)

   @inbounds for i in 1:nsrc
      key = (Symbol(:source, i), Symbol(:adis, i))
      adis[i] = pvals[key]
   end
   return adis
end


# --------------------------------------------------------------------------------------------------
# Build lens model from physical paramerters
# --------------------------------------------------------------------------------------------------
function build_lens(model::ModelConfig, pvals::Dict{Tuple{Symbol,Symbol}, Float64})
   # Determine the number of components from the lens model container
   n_lens = length(model.lens_config.components)

   # Initialize an empty vector to store lens parameters
   lens_vector = NamedTuple[]

   # Iterate over each lens component
   components = model.lens_config.components

   for i in 1:n_lens
      lens_id = Symbol(:lens, i)
      lens_params = Dict{Symbol, Union{Symbol, Float64}}()
      
      for (k, v) in enumerate(components)
         if v.owner == lens_id
            lens_params[:lens] = v.name
         end
      end

      for (k, v) in pvals
         if k[1] == lens_id
            lens_params[k[2]] = v
         end
      end
      push!(lens_vector, (; lens_params...))
   end
   return Lenses.init_CompositeLens(lens_vector)
end


# --------------------------------------------------------------------------------------------------
# All lensing quantities needed for log-likelihood 
# --------------------------------------------------------------------------------------------------
function lens_quantities(model::ModelConfig, lens::Lenses.AbstractLens)
   # Count the total number of knots in the lens model
   n_knots = sum(length(s.knots) for s in model.source_config.sources)

   # Allocate outputs (vector of vectors)
   ψ_all  = Vector{Vector{Float64}}(undef, n_knots)
   αx_all = Vector{Vector{Float64}}(undef, n_knots)
   αy_all = Vector{Vector{Float64}}(undef, n_knots)
   A_all  = Vector{NTuple{4, Vector{Float64}}}(undef, n_knots)

   kid = 1
   for src in model.source_config.sources
      for knot in src.knots
         # One image system knot positions
         x = knot.x
         y = knot.y

         # Potential
         ψ_all[kid] = Lenses.get_potential(lens, x, y)

         # Deflection
         αx_all[kid], αy_all[kid] = Lenses.get_deflection(lens, x, y)

         # Deformation tensor
         ψxx, ψyy, ψxy = Lenses.get_jacobian(lens, x, y)
         A_all[kid] = (ψxx, ψxy, copy(ψxy), ψyy)
                  
         # Increment
         kid = kid + 1
      end
   end
   return ψ_all, αx_all, αy_all, A_all
end


# --------------------------------------------------------------------------------------------------
# Free parameter names
# --------------------------------------------------------------------------------------------------
function free_parameter_names(model::ModelConfig)
   return [(p.owner, p.name) for p in model.parameters[model.free_param_idxs]]
end


# --------------------------------------------------------------------------------------------------
# Gelman-Rubin diagnostic
# --------------------------------------------------------------------------------------------------
function calculate_gr(chains::Array{Float64, 3}; burn_in::Float64=0.3)
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

function print_gr_report(chains::Array{Float64, 3}; param_names=nothing, burn_in=0.3)
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
function time_series_diagnostics(chains::Array{Float64, 3}; param_names=nothing, burn_in::Float64=0.3)
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


function acceptance_diagnostics(chains::Array{Float64, 3}; burn_in::Float64=0.3)
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


# --------------------------------------------------------------------------------------------------
# Get best-fit from the full optimization/MCMC results
# --------------------------------------------------------------------------------------------------
function get_best_fit(results, chains=nothing)
   best_θ = nothing
   best_logL = -Inf

   # Case A: Input is from the Parallel Optimizer (Vector of NamedTuples)
   if results isa Vector && eltype(results) <: NamedTuple
      # We already sorted the results in descending order of logL, so results[1] is the best fit
      best_run = results[1]
      best_θ = best_run.θ
      best_logL = best_run.f

   # Case B: Input is from MCMC (ll_history matrix)
   elseif results isa Matrix{Float64} && chains !== nothing
      # Get (step, chain) of the best logL by finding the index of the maximum value in the logL matrix
      best_idx = argmax(results) 
      step, chain_num = best_idx[1], best_idx[2]
      
      best_θ = chains[step, chain_num, :]
      best_logL = results[step, chain_num]
   else
      error("Please provide either Optimization results or both MCMC ll_history and chains.")
    end

   # Calculate Chi2: χ² = -2 * logL
   # Note: If your objective includes priors (log-posterior), 
   # this technically gives you the MAP (Maximum A Posteriori) chi2.
   chi2 = -2.0 * best_logL

   return best_θ, best_logL, chi2
end


function get_best_fit_with_errors(chains::Array{Float64, 3}, lls::Matrix{Float64}; burn_in=0.2, thinning=100)
   # Get chain details 
   n_steps, n_chains, n_params = size(chains)
    
   # Extract Best-Fit (Maximum Likelihood Estimate)
   best_idx = argmax(lls) 
   step_bf, chain_bf = best_idx[1], best_idx[2]
   
   # Best parameter and LogL values
   best_θ = chains[step_bf, chain_bf, :]
   best_logL = lls[step_bf, chain_bf]

   # Remove Burn-in
   start_idx = Int(floor(n_steps * burn_in)) + 1
   thinned_chain = chains[start_idx:thinning:end, :, :]

   # Reshape the thinned chain into a flat array
   flat_chain = reshape(thinned_chain, :, n_params)
    
   # Calculate (asymmetric) errors as we are defining errors relative to the Best-Fit value
   lower_err = zeros(n_params)
   upper_err = zeros(n_params)

   for i in 1:n_params
      # Get 16th and 84th percentiles of the posterior
      q16, q84 = StatsBase.quantile(flat_chain[:, i], [0.16, 0.84])
        
      # Asymmetric error: distance from best-fit to the quantiles
      lower_err[i] = best_θ[i] - q16
      upper_err[i] = q84 - best_θ[i]
   end

   # Return as a NamedTuple for easy access
   return best_θ, lower_err, upper_err, best_logL, -2.0 * best_logL
end


# --------------------------------------------------------------------------------------------------
# Get best-fit model RMS and parity check
# --------------------------------------------------------------------------------------------------
function check_parity(chains::Array{Float64, 3}, lls::Matrix{Float64})
   
end

function get_best_fit_rms(chains::Array{Float64, 3}, lls::Matrix{Float64}; check_parity::Bool=false)
   
end

end