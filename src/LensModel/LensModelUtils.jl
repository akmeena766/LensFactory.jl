module LensModelUtils

# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using Random
using Statistics

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
   P_all  = Vector{Vector{Int64}}(undef, n_knots)

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
         
         # Parity
         P_all[kid] = sign.( @. ψxx * ψyy - ψxy * ψxy )
         
         # Increment
         kid = kid + 1
      end
   end
   return ψ_all, αx_all, αy_all, A_all, P_all
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
function calculate_gr(chains::Array{Float64, 3}; burn_in::Float64=0.5)
   n_iter, n_params, n_chains = size(chains)
   
   if n_chains < 2
      return fill(NaN, n_params)
   end

   # Remove burn-in (the adaptation blocks)
   start_idx = Int(floor(burn_in * n_iter))
   samples = chains[start_idx:end, :, :]

   n = size(samples, 1) # Samples per chain
   m = n_chains        # Number of chains
   r_hats = zeros(Float64, n_params)

   for p in 1:n_params
      param_samples = samples[:, p, :] # [n, m]
        
      # W: Within-chain variance
      chain_vars = Statistics.var(param_samples, dims=1)
      W = Statistics.mean(chain_vars)
        
      # B: Between-chain variance
      chain_means = Statistics.mean(param_samples, dims=1)
      grand_mean = Statistics.mean(chain_means)
      B = (n / (m - 1)) * sum((chain_means .- grand_mean).^2)
        
      # V_hat: Pooled variance estimate
      V_hat = ((n - 1) / n) * W + (1 / n) * B
      r_hats[p] = sqrt(V_hat / W)
   end
   return r_hats
end

function print_gr_report(chains::Array{Float64, 3}; param_names=nothing, burn_in=0.5)
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
        owner = string(param_names[i][1])
        param = string(param_names[i][2])
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
   
end