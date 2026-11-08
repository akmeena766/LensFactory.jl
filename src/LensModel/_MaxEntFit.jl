# --------------------------------------------------------------------------------------------------
# Helper functions
# --------------------------------------------------------------------------------------------------
_table_row(cells, widths) = join("| " * rpad(cell, w) for (cell, w) in zip(cells, widths)) * "|"

function _smooth_map(map::Matrix{Float64}, σ::Float64)
   if σ <= 0
      return copy(map)
   end
 
   # Normalised 1D Gaussian weights, truncated at 3σ
   h = ceil(Int64, 3 * σ)
   w = [exp(-0.5 * (k / σ)^2) for k in -h:h]
   w = w ./ sum(w)
 
   nx, ny = size(map)
   tmp = similar(map)
   out = similar(map)
   @inbounds for j in 1:ny, i in 1:nx
      s = 0.0
      for k in -h:h
         s = s + w[k + h + 1] * map[clamp(i + k, 1, nx), j]
      end
      tmp[i, j] = s
   end
   @inbounds for j in 1:ny, i in 1:nx
      s = 0.0
      for k in -h:h
         s = s + w[k + h + 1] * tmp[i, clamp(j + k, 1, ny)]
      end
      out[i, j] = s
   end
   return out
end


# --------------------------------------------------------------------------------------------------
# Log-likelihood
# --------------------------------------------------------------------------------------------------
function log_likelihood(model::ModelConfig, κ::Matrix{Float64}, K::NamedTuple, adis::Vector{Float64})
   if model.sampler.scheme != :SourcePlane || model.source_config.use_flux || model.source_config.use_time_delay
      error("MaxEnt supports only image positions in the source plane (no flux or time delays).")
   end
 
   # Deflection and deformation tensor at the observed image positions
   αx_all, αy_all = LensModelUtils.lens_quantities_def(model, K, κ)
   A_all          = LensModelUtils.lens_quantities_jac(model, K, κ)
 
   # Calculate position likelihood
   logL, _, _ = Likelihood.logL_sourceplane(model, adis, αx_all, αy_all, A_all)
   return logL
end


# MaxEnt: gradient of log_likelihood with respect to κ (an nx × ny matrix like κ). The derivatives
# with respect to the lens quantities at each image are spread over the pixels with the kernel.
function log_likelihood_grad(model::ModelConfig, κ::Matrix{Float64}, K::NamedTuple, adis::Vector{Float64})
   # Deflection and deformation tensor at the observed image positions
   αx_all, αy_all = LensModelUtils.lens_quantities_def(model, K, κ)
   A_all          = LensModelUtils.lens_quantities_jac(model, K, κ)
 
   # Derivatives with respect to the lens quantities, one vector per knot
   gαx, gαy, gψxx, gψyy, gψxy = Likelihood.logL_sourceplane_grad(model, adis, αx_all, αy_all, A_all)
 
   # Kernel columns follow the knot order, so concatenating the knots gives one value per column
   g = zeros(length(κ))
   mul!(g, K.αx,  reduce(vcat, gαx),  1.0, 1.0)
   mul!(g, K.αy,  reduce(vcat, gαy),  1.0, 1.0)
   mul!(g, K.ψxx, reduce(vcat, gψxx), 1.0, 1.0)
   mul!(g, K.ψyy, reduce(vcat, gψyy), 1.0, 1.0)
   mul!(g, K.ψxy, reduce(vcat, gψxy), 1.0, 1.0)
   return reshape(g, size(κ))
end



# --------------------------------------------------------------------------------------------------
# Log-prior
# --------------------------------------------------------------------------------------------------
function log_prior(model::ModelConfig, κ::Matrix{Float64}, κp::Matrix{Float64})
   α = model.lens_config.alpha
   return -α * sum(@. κp - κ + κ * log(κ / κp))
end


function log_prior_grad(model::ModelConfig, κ::Matrix{Float64}, κp::Matrix{Float64})
   α = model.lens_config.alpha
   return @. -α * log(κ / κp)
end



# --------------------------------------------------------------------------------------------------
# Posterior (may return -Inf)
# --------------------------------------------------------------------------------------------------
# MaxEnt: log-posterior of the convergence map κ and its gradient with respect to κ
function log_posterior(model::ModelConfig, κ::Matrix{Float64}, κp::Matrix{Float64}, K::NamedTuple, adis::Vector{Float64})
   return log_prior(model, κ, κp) + log_likelihood(model, κ, K, adis)
end


function log_posterior_grad(model::ModelConfig, κ::Matrix{Float64}, κp::Matrix{Float64}, K::NamedTuple, adis::Vector{Float64})
   return log_prior_grad(model, κ, κp) + log_likelihood_grad(model, κ, K, adis)
end


function _neg_log_posterior(model::ModelConfig, θ::Vector{Float64}, κp::Matrix{Float64}, K::NamedTuple, adis::Vector{Float64})
   κ = exp.(reshape(θ, size(κp)))
   return -log_posterior(model, κ, κp, K, adis)
end


function _neg_log_posterior_grad!(G::Vector{Float64}, model::ModelConfig, θ::Vector{Float64}, κp::Matrix{Float64}, K::NamedTuple, adis::Vector{Float64})
   κ = exp.(reshape(θ, size(κp)))
   G .= vec(-log_posterior_grad(model, κ, κp, K, adis) .* κ)
   return G
end


# --------------------------------------------------------------------------------------------------
# Optimizer run
# --------------------------------------------------------------------------------------------------
function run_optimizer(model::ModelConfig, param_ref::Dict{Tuple{Symbol,Symbol}, <:Real}, method_cfg::MaxEntOptimizer, opt_cfg::Dict{Symbol, Any}, verbose::Bool)
   # Angular-diameter distance ratios at the reference parameters (fixed during the fit)
   cosmo = LensModelUtils.current_cosmology(model, param_ref)
   adis  = LensModelUtils.adis_current(model, param_ref, cosmo)

   # Kernel at the observed image positions (built once)
   K = LensModelUtils.build_kernel(model)

   # Default map: prior and starting map of the first run
   κp = copy(model.lens_config.pixels.κ)
   if any(<=(0.0), κp)
      error("The default map must be positive in every pixel.")
   end
   κ = copy(κp)

   # L-BFGS settings (sampling[:optimizer][:LBFGS])
   method  = Optim.LBFGS(linesearch   = Optim.LineSearches.HagerZhang())
   options = Optim.Options(iterations = Int64(get(opt_cfg, :iterations, 1000)),
                           f_reltol   = Float64(get(opt_cfg, :tolerance, 1e-6)))

   # Results table: header now, one row per run
   names  = ("Run", "Iterations", "Converged", "χ²", "Log-prior", "Log-posterior", "Change", "Time (s)")
   widths = (5, 12, 11, 14, 14, 15, 10, 10)
   line   = "-" ^ length(_table_row(names, widths))
   if verbose
      println(line)
      println(_table_row(names, widths))
      println(line)
   end
   
   runs      = NamedTuple[]
   converged = false
   κ_last    = κ
   for run in 1:method_cfg.n_runs
      # This run
      t_start = time()
      res = Optim.optimize(θ      -> _neg_log_posterior(model, θ, κp, K, adis),
                           (G, θ) -> _neg_log_posterior_grad!(G, model, θ, κp, K, adis),
                           vec(log.(κ)), 
                           method, 
                           options)
      t_run = time() - t_start

      # Result of this run
      κ      = exp.(reshape(Optim.minimizer(res), size(κp)))
      lprior = log_prior(model, κ, κp)
      logL   = log_likelihood(model, κ, K, adis)
      lpost  = lprior + logL

      # Relative change of the map from the previous run
      if run > 1
         change = sqrt(sum(abs2, κ .- κ_last) / sum(abs2, κ))
      else
         change = NaN
      end
      κ_last = κ

      push!(runs, (log_posterior  = lpost, 
                   log_likelihood = logL, 
                   log_prior      = lprior,
                   iterations     = Optim.iterations(res), 
                   converged      = Optim.converged(res),
                   change         = change,
                   time           = t_run))
      
      # Save this run
      if method_cfg.save_runs
         mkpath(method_cfg.save_path)
         jldsave(joinpath(method_cfg.save_path, "$(model.observation.lens)_MaxEnt_run$(run).jld2");
                 kappa = κ, prior = κp, run = runs[end])
      end

      if verbose
         if Optim.converged(res)
            run_converged = "yes"
         else
            run_converged = "no"
         end

         if run > 1
            run_change = @sprintf("%.2e", change)
         else
            run_change = "-"
         end
         row = (run, 
                Optim.iterations(res), 
                run_converged, 
                @sprintf("%.4f", -2.0 * logL),
                @sprintf("%.4f", lprior), 
                @sprintf("%.4f", lpost), 
                run_change,
                @sprintf("%.2f", t_run))
         println(_table_row(row, widths))
      end

      # Stop when the map no longer changes between runs
      if run > 1 && change <= method_cfg.tolerance
         converged = true
         break
      end

      # Smoothed result: prior and starting map of the next run
      if run < method_cfg.n_runs
         κp = _smooth_map(κ, method_cfg.smoothing)
         κ  = copy(κp)
      end
   end

   if verbose
      t_total = @sprintf("%.2f", sum(r.time for r in runs))
      if converged
         summary = "✅ Converged after $(length(runs)) runs in $(t_total) s (change ≤ $(method_cfg.tolerance))."
      else
         summary = "❌ Not converged after $(length(runs)) runs in $(t_total) s (n_runs = $(method_cfg.n_runs))."
      end
      println(line)
      println(_table_row((summary,), (length(line) - 3,)))
      println(line)
   end

   # Map and prior of the last run
   return (kappa = κ, prior = κp, runs = runs, converged = converged)
end
