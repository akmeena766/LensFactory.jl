function _get_levels(dens; quantiles=[0.393, 0.865, 0.989])
   sorted_dens = sort(vec(dens), rev=true)
   cumulative_dens = cumsum(sorted_dens) ./ sum(sorted_dens)
   return [sorted_dens[findfirst(x -> x >= q, cumulative_dens)] for q in quantiles]
end

# --------------------------------------------------------------------------------------------------
# Corner Plot
# --------------------------------------------------------------------------------------------------
PARAM_NAME = Dict(:lens1 => "L1", :lens2 => "L2", :lens3 => "L3", :lens4 => "L4", :lens5 => "L5", 
                  :lens6 => "L6", :lens7 => "L7", :lens8 => "L8", :lens9 => "L9", :lens10 => "L10",
                  :scaling => "S", 
                  :x_c => "x_c", :y_c => "y_c", :eps => "\\epsilon", :pa => "\\phi", 
                  :v_d => "v_d", :mass => "\\log[M]", :c => "c",
                  :gamma1 => "\\gamma_1", :gamma2 => "\\gamma_2",
                  :ref_sigma => "\\sigma_{\\star}", :ref_core => "\\theta_{c,\\star}", :ref_cut => "\\theta_{t,\\star}")


"""
    LensFactory.LensModel.plot_corner(chains, lls; param_names=nothing, burn_in=0.3, thinning=100)
Generates a corner plot for the given MCMC chains and log-likelihoods.
# Arguments
- `chains`: MCMC chains of shape (n_steps, n_chains, n_params).
- `lls`: Log-likelihoods corresponding to the chains, of shape (n_steps, n_chains).
- `param_names=nothing`: Optional list of parameter names for labeling the axes.
- `burn_in=0.3`: Fraction of the initial samples to discard as burn-in.
- `thinning=100`: Interval for thinning the chains to reduce autocorrelation.
# Returns
- A Makie figure object containing the corner plot.
"""
function LensFactory.LensModel.plot_corner(chains, lls; param_names=nothing, burn_in=0.3, thinning=100)
   # Get chain details 
   n_steps, n_chains, n_params = size(chains)
   
   # Get best-fit parameter and errors
   best_θ, lower_err, upper_err, _, _ = LensFactory.LensModel.get_best_fit_with_errors(chains, lls; burn_in=burn_in, thinning=thinning)

   # Remove Burn-in
   start_idx = Int(floor(n_steps * burn_in)) + 1
   thinned_chain = chains[start_idx:thinning:end, :, :]

   # Reshape the thinned chain into a flat array
   flat_chain = reshape(thinned_chain, :, n_params)
    
   # Initialize figure
   fig = Figure(size=(1200, 1200), figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))

   for i in 1:n_params
      # Get parameter labels
      p_name = PARAM_NAME[param_names[i][2]]

      # Marginal Stats for Title
      title_str = L"%$(p_name) = %$(round(best_θ[i], digits=2))^{+ %$(round(upper_err[i], digits=2))}_{-%$(round(lower_err[i], digits=2))}"
      
      for j in 1:i
         p_name_i = L"%$(PARAM_NAME[param_names[i][1]]): %$(PARAM_NAME[param_names[i][2]])"
         p_name_j = L"%$(PARAM_NAME[param_names[j][1]]): %$(PARAM_NAME[param_names[j][2]])"
         ax = Axis(fig[i, j]; 
                  title = (i == j) ? title_str : "",
                  xlabel = (i == n_params) ? p_name_j : "",
                  ylabel = (j == 1) ? p_name_i : "",
                  titlesize = 18,
                  xticklabelrotation = π/4)
         
         set_plotKws!(ax)
         
         if i == j
               # Diagonal: 1D Density
               density!(ax, flat_chain[:, i], color=(:dodgerblue, 0.5), strokewidth=2)
               
               # Sigma lines
               vlines!(ax, [best_θ[i] - lower_err[i], best_θ[i], best_θ[i] + upper_err[i]], color=:black, linestyle=:dash)

               ylims!(ax, 0, nothing)
         else
               # 2D Corner: KDE Contours
               # Note: j is x-axis (column), i is y-axis (row)
               k = kde((flat_chain[:, j], flat_chain[:, i]))
               levels = _get_levels(k.density)
               
               contour_levels = sort(unique([levels..., maximum(k.density) + 1e-10]))
               contourf!(ax, k.x, k.y, k.density, levels = contour_levels, 
                        colormap = [:white, (:dodgerblue, 0.3), (:dodgerblue, 0.6), (:dodgerblue, 0.9)], extendlow = :auto)
               contour!(ax, k.x, k.y, k.density, levels = levels, color = :black, linewidth = 1.2)
         end

         # Clean up decorations for inner plots
         if i < n_params hidexdecorations!(ax, grid = false, ticks = false) end
         if j > 1 hideydecorations!(ax, grid = false, ticks = false) end
      end
   end
   
   save("./corner.png", fig)
   return fig
end

"""
    LensFactory.LensModel.plot_trace(chains; param_names=nothing, burn_in=0.0, thinning=1)
Generates trace plots for the given MCMC chains.
# Arguments
- `chains`: MCMC chains of shape (n_steps, n_chains, n_params).
- `param_names=nothing`: Optional list of parameter names for labeling the axes.
- `burn_in=0.0`: Fraction of the initial samples to discard as burn-in.
- `thinning=1`: Interval for thinning the chains to reduce autocorrelation.
# Returns
- A Makie figure object containing the trace plots.
"""
function LensFactory.LensModel.plot_trace(chains; param_names=nothing, burn_in=0.0, thinning=1)
   # Adapt to [n_params, n_chains, n_steps]
   n_steps, n_chains, n_params = size(chains)
    
   # Calculate start index based on burn-in (slicing the 3rd dimension)
   start_idx = max(1, Int(floor(n_steps * burn_in)) + 1)
   steps_range = start_idx:thinning:n_steps
    
   # Create the figure
   fig = Figure(size=(1000, 200 * n_params))

   for i in 1:n_params
      # Extract the parameter name
      p_name = param_names !== nothing ? string(param_names[i][2]) : "theta_$i"

      # Create Axis
      ax = Axis(fig[i, 1], 
               xscale = log10,
               ylabel = p_name,
               xlabel = i == n_params ? "Iteration (log10)" : "",
               xminorticksvisible = true, 
               xminorticks = IntervalsBetween(9),
               rightspinevisible = false,
               topspinevisible = false)
        
      # Plot each chain correctly
      for c in 1:n_chains
         # Extract: Parameter i, Chain c, and the range of Steps
         @views trace_data = chains[steps_range, c, i]
         lines!(ax, collect(steps_range), trace_data, alpha = 0.4, linewidth = 1)
      end
        
      # Formatting
      if i < n_params
         hidexdecorations!(ax, grid = false, ticks = false)
      end
   end

   save("./trace.png", fig)
   return fig
end

"""
    LensFactory.LensModel.plot_best_model(θx, θy, model, chains, lls)
Generates a plot comparing the observed image positions with the predicted image positions from the best-fit lens model.
# Arguments
- `θx`: Matrix of x-coordinates of the image plane grid.
- `θy`: Matrix of y-coordinates of the image plane grid.
- `model`: The lens model configuration.
- `chains`: MCMC chains of shape (n_steps, n_chains, n_params).
- `lls`: Log-likelihoods corresponding to the chains, of shape (n_steps, n_chains).
# Returns
- A Makie figure object containing the comparison plot.
"""
function LensFactory.LensModel.plot_best_model(θx::Matrix{<:RV}, θy::Matrix{<:RV}, model::LensModel.ModelConfig, chains::Array{Float64, 3}, lls::Matrix{Float64})
   # Get the best parameters based on likelihood (lls)
   best_θ, _, _ = LensModel.get_best_fit(lls, chains)

   # Get list of parameters for the lens model
   param_ref = Dict(p.key => p.refer for p in model.parameters)
   
   # Replace free parameter values by best-fit values
   pvals = LensModel.param_dict(model, best_θ, param_ref)

   # Get best-fit model
   best_model = LensModel.build_lens(model, pvals)

   # Get angular-diameter distance ratios
   adis = LensModel.adis_current(model, pvals)

   # Calculate deformation at all image positions
   ψ_all, αx_all, αy_all, A_all = LensModel.lens_quantities(model, best_model)

   # Identity tuple
   I4 = (1.0, 0.0, 0.0, 1.0)

   # Initialize empty figure
   fig = Figure(size=(600, 600), figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))
   ax = Axis(fig[1, 1])
   
   # Calculate RMS for each image
   sid = 1
   kid = 1
   for src in model.source_config.sources
      # Get angular-diameter distance ratio for this source
      adis_value = adis[sid]
      
      # Generate source id
      src_id = Symbol(:src, sid)
      for knot in src.knots
         # Generate knot id
         knot_id = Symbol(:knot, kid)

         # Knot positions and measurement errors
         x  = knot.x
         y  = knot.y
         σx = knot.σx
         σy = knot.σy
         σθ = knot.σθ

         # Plot observed positions of knots
         scatter!(ax, x, y, markersize=20, marker=:circle, color=:transparent, strokecolor=:black, strokewidth=2)

         # Number of images for this knot
         n = length(x)

         # Deflection vector at the knot positions
         αx = @. adis_value * αx_all[kid]
         αy = @. adis_value * αy_all[kid]

         # Deformation tensor at the knot positions
         A = @. adis_value * A_all[kid]
         for i in eachindex(A)
            @. A[i] = I4[i] - A[i]
         end

         # Individual source positions using broadcasting
         βx_ind = @. x - αx
         βy_ind = @. y - αy

         # Get weighted source position (Section 4.1 in https://arxiv.org/pdf/astro-ph/0102340)
         βx_model, βy_model, _ = LensModel.Likelihood._weighted_position(βx_ind, βy_ind, A, σx, σy, σθ, n)

         # Get image positions
         predicted_image = Lenses.get_image(best_model, θx, θy, adis_value, (βx_model, βy_model))

         # Plot predicted image positions
         scatter!(ax, first.(predicted_image), last.(predicted_image), markersize=20, marker=:diamond, color=:transparent, strokecolor=:red, strokewidth=2)
         kid = kid + 1
      end
      sid = sid + 1
   end

   xlims!(ax, minimum(θx), maximum(θx))
   ylims!(ax, minimum(θy), maximum(θy))
   return fig
end