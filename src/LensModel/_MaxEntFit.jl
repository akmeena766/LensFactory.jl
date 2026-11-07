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
   κ     = copy(κp)
   nx_ny = size(κp)

   # L-BFGS settings (sampling[:optimizer][:LBFGS])
   method  = Optim.LBFGS(linesearch   = Optim.LineSearches.HagerZhang())
   options = Optim.Options(iterations = Int64(get(opt_cfg, :iterations, 1000)),
                           f_reltol   = Float64(get(opt_cfg, :tolerance, 1e-6)))

   # Results table: header now, one row per run
   cols = (("Run", 5), ("Iterations", 12), ("Converged", 11), ("χ²", 14), ("Log-prior", 14), 
           ("Log-posterior", 15), ("Change", 10), ("Time (s)", 10))
   line = "-" ^ sum(w + 2 for (_, w) in cols)
   if verbose
      println(line)
      println(join("| " * rpad(name, w) for (name, w) in cols) * "|")
      println(line)
   end

   runs      = NamedTuple[]
   converged = false
   for run in 1:method_cfg.n_runs
      t_start = time()
      
   end
end
