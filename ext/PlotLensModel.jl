function _get_levels(dens; quantiles=[0.393, 0.865, 0.989])
   sorted_dens = sort(vec(dens), rev=true)
   cumulative_dens = cumsum(sorted_dens) ./ sum(sorted_dens)
   return [sorted_dens[findfirst(x -> x >= q, cumulative_dens)] for q in quantiles]
end

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
      p_name = string(param_names[i][2])

      # Marginal Stats for Title
      title_str = L"%$(p_name) = %$(round(best_θ[i], digits=2))^{+ %$(round(upper_err[i], digits=2))}_{- %$(round(lower_err[i], digits=2))}"
      
      for j in 1:i
         p_name_j = string(param_names[j][2])
         
         ax = Axis(fig[i, j]; 
                  title = (i == j) ? title_str : "",
                  xlabel = (i == n_params) ? p_name_j : "",
                  ylabel = (j == 1) ? p_name : "",
                  titlesize = 18)
         
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