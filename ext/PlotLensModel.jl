# --------------------------------------------------------------------------------------------------
# Calculate the density threshold values corresponding to specific cumulative 
# probability contours (quantiles).
#
# Arguments
# - `dens::` A multi-dimensional array or vector containing density values (e.g., from a KDE).
#
# Keyword Arguments
# - `quantiles::Vector{Float64} = [0.393, 0.865, 0.989]`: A sorted list of target cumulative 
# probabilities. 
#
# Returns
# - `Vector{Float64}`: The threshold density values for each requested quantile, ordered the same 
# as `quantiles`
# --------------------------------------------------------------------------------------------------
function _get_levels(dens; quantiles=[0.393, 0.865, 0.989])
   sorted_dens = sort(vec(dens), rev=true)
   cumulative_dens = cumsum(sorted_dens) ./ sum(sorted_dens)
   return [sorted_dens[findfirst(x -> x >= q, cumulative_dens)] for q in quantiles]
end

# --------------------------------------------------------------------------------------------------
# Turns the user-facing `param_indices` argument into a concrete vector of
# integer positions into `free_parameter_names` / the third dimension of `chains`.
#
# Arguments
# - `free_parameter_names::Vector{Tuple{Symbol, Symbol}}`: list of free parameter names.
#
# Arguments
# - `param_ids`: can be:
#   - `Nothing`:                       every parameter
#   - `Symbol`:                        filter by component (1st tuple element), e.g. :lens1
#   - `Vector{Symbol}`:                filter by component (1st tuple element), e.g. :lens1 
#                                      or `[:lens1, :lens3]`
#   - `Vector{Tuple{Symbol, Symbol}}`: exact matches, e.g. [(:lens1, :eps), (:lens2, :mass)]
#   - `Vector{Int64}`:                 direct positions, e.g. [1, 3, 5]
#
# Returns
# - `Vector{Int64}`: indices corresponding to `free_parameter_names` based on `param_ids`
# --------------------------------------------------------------------------------------------------
function _resolve_param_indices(free_parameter_names::AbstractVector{Tuple{Symbol, Symbol}}, 
                                param_ids::Union{Nothing, Symbol, AbstractVector{Symbol}, AbstractVector{Int64}, AbstractVector{Tuple{Symbol, Symbol}}})
   # Get total number of free parameters
   n_total = length(free_parameter_names)
   
   # If no indices are provided, return all indices
   if isnothing(param_ids)
      return collect(1:n_total)
   end

   # Allow a single bare symbol, e.g. param_indices = :lens1
   selection = param_ids isa Symbol ? [param_ids] : param_ids

   if all(x -> x isa Integer, selection)
      # Direct positions, e.g. [1, 3, 5]
      idx = collect(Int64, selection)
      @assert all(1 .<= idx .<= n_total) "param_ids out of range 1:$(n_total)"

   elseif all(x -> x isa Symbol, selection)
      # Filter by component name only, e.g. :lens1 -> every (:lens1, *) entry
      idx = findall(p -> p[1] in selection, free_parameter_names)
      isempty(idx) && error("No entries in free_parameter_names match component(s) $(selection)")
   else
      # Full (component, param) pairs to look up directly
      idx = Int[]
      for entry in selection
         pos = findfirst(==(entry), free_parameter_names)
         isnothing(pos) && error("Could not find $(entry) in free_parameter_names")
         push!(idx, pos)
      end
   end
   return idx
end

# --------------------------------------------------------------------------------------------------
# Corner Plot
# --------------------------------------------------------------------------------------------------
PARAM_NAME = Dict{Symbol, String}(:lens1 => "L1", :lens2 => "L2", :lens3 => "L3", :lens4 => "L4", 
                                  :lens5 => "L5", :lens6 => "L6", :lens7 => "L7", :lens8 => "L8", 
                                  :lens9 => "L9", :lens10 => "L10",
                                  :x_c => "x_c", :y_c => "y_c", :eps => "\\epsilon", :pa => "\\phi", 
                                  :x_s => "x_s", :x_t => "x_t",
                                  :v_d => "v_d", :mass => "\\log[M]", :c => "c",
                                  :gamma => "\\gamma", :angle => "\\theta", :delta => "\\delta",
                                  :ref_sigma => "\\sigma_{\\star}", 
                                  :ref_core => "\\theta_{c,\\star}", 
                                  :ref_cut => "\\theta_{t,\\star}",
                                  :ref_mag     => "m_{\\star}",
                                  :slope_sigma => "\\lambda",
                                  :slope_core  => "\\beta",
                                  :slope_cut   => "\\alpha")

function get_param_name(key::Symbol)
   # Check static dictionary first
   haskey(PARAM_NAME, key) && return PARAM_NAME[key]
   
   # Match :lensN pattern
   m = match(r"^lens(\d+)$", string(key))
   if m !== nothing
      n = m.captures[1]
      return "L$n"
   end

   # Match :sourceN pattern
   m = match(r"^source(\d+)$", string(key))
   if m !== nothing
      n = m.captures[1]
      return "S$n"
   end

   # Match :adisN pattern
   m = match(r"^adis(\d+)$", string(key))
   if m !== nothing
      n = m.captures[1]
      return "a_{dis}"
   end

   # Match :scalingN pattern
   m = match(r"^scaling(\d+)$", string(key))
   if m !== nothing
      n = m.captures[1]
      return "L_{$n}^{\\star}"
   end
   
   error("Unknown parameter name: $key")
end

"""
    LensFactory.LensModel.plot_corner(chains, logL; free_parameter_names=nothing,
                                                    burn_in::Float64  = 0.3,
                                                    thin::Int64       = 100,
                                                    save_plot::Bool   = true,
                                                    plot_name::String = "./corner.png",
                                                    resolution::Int64 = 2)
Generates a corner plot for the given MCMC chains and log-likelihood values.

# Arguments
- `chains`: MCMC chains of shape (n_steps, n_chains, n_params).
- `logL`: Log-likelihood values corresponding to the chains, of shape (n_steps, n_chains).

# Keyword arguments
- `free_parameter_names=nothing`: Optional list of parameter names for labeling the axes.
- `burn_in::Float64 = 0.3`: Fraction of the initial samples to discard as burn-in.
- `thinn::Int64     = 100`: Interval for thinning the chains to reduce autocorrelation.
- `save_plot::Bool  = true`: Whether to save the plot as "corner.png".
   - `plot_name::String = "./corner.png"`: Filename for saving the plot.
   - `resolution::Int64 = 2`: Resolution for saving the plot.

# Returns
- A Makie figure object containing the corner plot.
"""
function LensFactory.LensModel.plot_corner(chains, logL; 
                                           free_parameter_names = nothing,
                                           plot_parameters      = nothing,
                                           burn_in::Float64     = 0.3,
                                           thin::Int64          = 100,
                                           save_plot::Bool      = true,
                                           plot_name::String    = "./corner.png",
                                           resolution::Int64    = 2)
   # Get chain details 
   n_steps, n_chains, n_params = size(chains)

   # Get parameter indices to plot
   plot_idx = _resolve_param_indices(free_parameter_names, plot_parameters)
   n_params = length(plot_idx)
   
   # Get best-fit parameter and errors
   best_θ_all, _, lower_err_all, upper_err_all = LensFactory.LensModel.get_best_fit_parameters(logL; 
                                                   chains       = chains, 
                                                    with_errors = true, 
                                                    burn_in     = burn_in, 
                                                    thin        = thin)

   # Get best-fit and errors for selected parameters
   best_θ    = best_θ_all[plot_idx]
   lower_err = lower_err_all[plot_idx]
   upper_err = upper_err_all[plot_idx]

   # Keep only the requested labels, in the requested order
   free_parameter_names = free_parameter_names[plot_idx]
   
   # Remove Burn-in
   start_idx = Int(floor(n_steps * burn_in)) + 1
   thinned_chain = chains[start_idx:thin:end, :, plot_idx]

   # Reshape the thinned chain into a flat array
   flat_chain = reshape(thinned_chain, :, n_params)
    
   # Dynamic figure sampling
   N = n_params
   base_size = max(1000, N * 90) 
   font_size = max(6, 20 - (N / 3))
   tick_size = max(5, 15 - (N / 3))
   
   # Create figure and pad in all directions
   fig = Figure(size = (base_size, base_size), 
               rowgap = 0, 
               colgap = 0, 
               figure_padding = 25,
               fontsize = font_size, 
               fonts=(; regular="Times New Roman"))
   
   for i in 1:n_params
      # Marginal Stats for Title
      title_str = L"%$(round(best_θ[i], digits=2))^{+ %$(round(upper_err[i], digits=2))}_{%$(round(lower_err[i], digits=2))}"
      
      for j in 1:i
         owner_i = get_param_name(free_parameter_names[i][1])
         param_i = get_param_name(free_parameter_names[i][2])

         owner_j = get_param_name(free_parameter_names[j][1])
         param_j = get_param_name(free_parameter_names[j][2])

         p_name_i = LaTeXString("\$$owner_i: $param_i\$")
         p_name_j = LaTeXString("\$$owner_j: $param_j\$")

         ax = Axis(fig[i, j];
                  xtickalign = 1, ytickalign = 1,
                  xticksize = 2, yticksize = 2,
                  xticklabelsize = tick_size, 
                  yticklabelsize = tick_size,
                  xgridvisible = false, ygridvisible = false,
                  titlevisible = false,
                  xticks = LinearTicks(3),
                  yticks = LinearTicks(3),
                  xlabel = (i == n_params) ? p_name_j : "",
                  ylabel = (j == 1) ? p_name_i : "",
                  xticklabelrotation = π/4)
         
         tightlimits!(ax)
         
         if i == j
            # Diagonal: 1D Density
            density!(ax, flat_chain[:, i], color=(:dodgerblue, 0.5), strokewidth=2)
               
            # Sigma lines
            vlines!(ax, [best_θ[i] + lower_err[i], best_θ[i], best_θ[i] + upper_err[i]], color=:black, linestyle=:dash)
            
            # Manual label creation and placing it above the diagonal plots
            Label(fig[i, j, Top()], title_str;
                    tellheight = false,
                    padding = (0, 0, 23, 0),
                    fontsize = font_size,
                    font = :bold,
                    halign = :center)
         else
            # Off-diagonal
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
         if i < n_params 
            hidexdecorations!(ax, grid = false, ticks = false)
         end
         
         if j > 1 
            hideydecorations!(ax, grid = false, ticks = false)
         end
      end
   end

   # Reset gaps to zero
   rowgap!(fig.layout, 0)
   colgap!(fig.layout, 0)

   if save_plot
      save(plot_name, fig, px_per_unit=resolution)
   end
   return fig
end


"""
    LensFactory.LensModel.plot_trace(chains)
Generates trace plots for the given MCMC chains.

# Arguments
- `chains`: MCMC chains of shape (n_steps, n_chains, n_params).

# Keyword arguments
- `free_parameter_names = nothing`: Optional list of parameter names for labeling the axes.
- `burn_in::Float64 = 0.0`: Fraction of the initial samples to discard as burn-in.
- `thin::Int64 = 1`: Interval for thinning the chains to reduce autocorrelation.
- `save_plot::Bool = true`: Whether to save the plot as "corner_optimizer.png".
   - `plot_name::String = "./trace.png"`: Filename for saving the plot.
   - `resolution::Int64 = 2`: Resolution for saving the plot.

# Returns
- A Makie figure object containing the trace plots.
"""
function LensFactory.LensModel.plot_trace(chains; 
                              free_parameter_names=nothing, 
                              burn_in::Float64  = 0.0,
                              thin::Int64       = 1,
                              save_plot::Bool   = true, 
                              plot_name::String = "./trace.png", 
                              resolution::Int64 = 2)
   # Adapt to [n_params, n_chains, n_steps]
   n_steps, n_chains, n_params = size(chains)
    
   # Calculate start index based on burn-in (slicing the 3rd dimension)
   start_idx = max(1, Int(floor(n_steps * burn_in)) + 1)
   steps_range = start_idx:thin:n_steps
    
   # Create the figure
   fig = Figure(size=(1000, 200 * n_params))

   for i in 1:n_params
      # Extract the parameter name
      p_name = free_parameter_names !== nothing ? string(free_parameter_names[i][2]) : "theta_$i"

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
   
   if save_plot
      save(plot_name, fig, px_per_unit=resolution)
   end
   return fig
end


function _ellipse(ex::RV, ey::RV, eθ::RV; x0::RV=0.0, y0::RV=0.0, n::Int64=100)
   # Numer of points
   t = range(0, 2π, length=n)

   # Generate points
   x = ex .* cos.(t)
   y = ey .* sin.(t)

   # Rotate
   θ_rad = deg2rad(eθ)
   x_rot = x .* cos(θ_rad) - y .* sin(θ_rad)
   y_rot = x .* sin(θ_rad) + y .* cos(θ_rad)

   # Translate
   x_rot .+= x0
   y_rot .+= y0

   return x_rot, y_rot
end

"""
    LensFactory.LensModel.plot_best_model(model, chains, logL)
Generates a plot comparing the observed image positions with the predicted image positions from the best-fit lens model.

# Arguments
- `model`: The lens model configuration.
- `chains`: MCMC chains of shape (n_steps, n_chains, n_params).
- `logL`: Log-likelihood values corresponding to the chains, of shape (n_steps, n_chains).

# Keyword arguments
- `source::Union{Nothing, Integer} = nothing`: Which system to plot. `nothing` plots every system and
   draws the critical curves/caustics at the redshift `z_s`. An integer plots only that system and
   draws the critical curves/caustics at that system's own angular-diameter distance ratio, in which
   case `z_s` is ignored.
- `z_s::RV = 1.5`: Source redshift.
- `plot_error::Bool = true`: Whether to plot the error ellipses for the observed image positions.
- `plot_critical::Bool = true`: Whether to plot the critical curves.
   - `critical_tan_kws::NamedTuple = (color=:red, linewidth=2, linestyle=:solid)`: Keyword arguments for plotting the tangential critical curve.
   - `critical_rad_kws::NamedTuple = (color=:red, linewidth=2, linestyle=:dash)`: Keyword arguments for plotting the radial critical curve.
- `plot_caustics::Bool = true`: Whether to plot the caustics.
   - `caustic_tan_kws::NamedTuple = (color=:green, linewidth=2, linestyle=:solid)`: Keyword arguments for plotting the tangential caustic.
   - `caustic_rad_kws::NamedTuple = (color=:green, linewidth=2, linestyle=:dash)`: Keyword arguments for plotting the radial caustic.
- `save_plot::Bool = true`: Whether to save the plot as "best_model.png".
   - `plot_name::String = "./best_model.png"`: Filename for saving the plot.
   - `resolution::Int64 = 2`: Resolution for saving the plot.

# Returns
- A Makie figure object containing the comparison plot.
"""
function LensFactory.LensModel.plot_best_model(model::LensModel.ModelConfig, 
                                              chains::Array{Float64, 3}, 
                                              logL::Matrix{Float64};
                                              source::Union{Nothing, Int64} = nothing,
                                              z_s::RV = 1.5,
                                              plot_errors::Bool = true,
                                              plot_critical::Bool = true,
                                              critical_rad_kws::NamedTuple = (color=:red, linewidth=2, linestyle=:dash),
                                              critical_tan_kws::NamedTuple = (color=:red, linewidth=2, linestyle=:solid),
                                              plot_caustics::Bool = true,
                                              caustic_tan_kws::NamedTuple = (color=:green, linewidth=2, linestyle=:solid),
                                              caustic_rad_kws::NamedTuple = (color=:green, linewidth=2, linestyle=:dash),
                                              save_plot::Bool = true,
                                              plot_name::String = "./best_model.png",
                                              resolution::Int64 = 2)
   # Get the best parameters based on minimum log-likelihood
   best_θ, _ = LensFactory.LensModel.get_best_fit_parameters(logL; chains=chains)

   # Get list of parameters for the lens model
   param_ref = Dict(p.key => p.refer for p in model.parameters)
   
   # Replace free parameter values by best-fit values
   pvals = LensFactory.LensModel.LensModelUtils.param_dict(model, best_θ, param_ref)
   
   # Update cosmology
   cosmo = LensFactory.LensModel.LensModelUtils.current_cosmology(model, pvals)
   
   # Get best-fit model
   best_model = LensFactory.LensModel.LensModelUtils.build_lens(model, pvals, cosmo)

   # Generate grid
   FOV = model.observation.FOV
   pixel_scale = model.observation.pixel_scale
   x_grid, y_grid = Lenses.get_meshgrid(0.5 * FOV[1], 0.5 * FOV[2], pixel_scale)

   # Get angular-diameter distance ratios
   adis = LensFactory.LensModel.adis_current(model, pvals, cosmo)

   # Systems to plot: a single system, or all of them
   n_sources = length(model.source_config.sources)
   if isnothing(source)
      plot_sids = collect(1:n_sources)
   else
      plot_sids = [Int64(source)]
   end

   # Number of knots preceding each source, so `kid` stays aligned with αx_all / A_all
   knot_offsets = cumsum([0; [length(s.knots) for s in model.source_config.sources]])

   # Calculate deflection and deformation tensor at all image positions
   αx_all, αy_all = LensFactory.LensModel.LensModelUtils.lens_quantities_def(model, best_model)
   A_all          = LensFactory.LensModel.LensModelUtils.lens_quantities_jac(model, best_model)

   # Identity tuple
   I4 = (1.0, 0.0, 0.0, 1.0)

   # Initialize empty figure
   fig = Figure(size=(600, 600), figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))
   ax = Axis(fig[1, 1])

   # Set plot keywords
   LensFactory.Lenses.set_plotKws!(ax)
   
   # Calculate RMS for each image
   first_plot = true
   for sid in plot_sids
      # Source configuration for this system
      src = model.source_config.sources[sid]

      # Get angular-diameter distance ratio for this source
      adis_value = adis[sid]

      # Knot counter, offset into the global knot list
      kid = knot_offsets[sid]
      
      for knot in src.knots
         # Generate knot id
         kid = kid + 1

         # Knot positions and measurement errors
         x  = knot.x
         y  = knot.y
         σx = knot.σx
         σy = knot.σy
         σθ = knot.σθ

         # Plot observed positions of knots
         scatter!(ax, x, y; markersize  = 20, 
                            marker      = :circle, 
                            color       = :transparent, 
                            strokecolor = :black, 
                            strokewidth = 2, 
                            label       = first_plot ? "Observed" : nothing)
         if plot_errors
            for i in 1:size(x, 1)
               e_x, e_y = _ellipse(σx[i], σy[i], σθ[i]; x0=x[i], y0=y[i])
               lines!(ax, e_x, e_y, color=:black, linewidth=2)
            end
         end

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
         β_ind = @. tuple(x - αx, y - αy)

         # Get weighted source position
         β_model, _, _ = LensModel.Likelihood._weighted_position(β_ind, A, σx, σy, σθ, n)
         βx_model, βy_model = β_model

         # Get image positions
         predicted_image = Lenses.get_image(best_model, x_grid, y_grid, adis_value, (βx_model, βy_model))

         # Plot predicted image positions
         scatter!(ax, first.(predicted_image), last.(predicted_image); 
                     markersize  = 20,
                     marker      = :diamond, 
                     color       = :transparent, 
                     strokecolor = :dodgerblue, 
                     strokewidth = 2, 
                     label       = first_plot ? "Predicted" : nothing)
         
         # Update the label to false
         first_plot = false
      end
   end

   if plot_critical || plot_caustics
      if isnothing(source)
         # No single system to inherit from: use the requested source redshift
         Dls = Cosmology.angular_diameter_distance(cosmo, model.observation.z_d, z_s)
         Dos = Cosmology.angular_diameter_distance(cosmo, 0.0, z_s)
         adis_crit = Dls / Dos
      else
         # Curves at the redshift of the system being plotted
         adis_crit = adis[source]
      end      

      if plot_critical
         # Get critical curves
         critical_tan, critical_rad = Lenses.get_critical_curve(best_model, x_grid, y_grid, adis_crit)

         # Plot tangential critical curve
         for curve in critical_tan
            lines!(ax, first.(curve), last.(curve); critical_tan_kws...)
         end

         # Plot radial critical curve
         for curve in critical_rad
            lines!(ax, first.(curve), last.(curve); critical_rad_kws...)
         end
      end

      if plot_caustics
         # Get caustics
         caustic_tan, caustic_rad = Lenses.get_caustic(best_model, x_grid, y_grid, adis_crit)

         # Plot tangential caustic
         for curve in caustic_tan
            lines!(ax, first.(curve), last.(curve); caustic_tan_kws...)
         end

         # Plot radial caustic
         for curve in caustic_rad
            lines!(ax, first.(curve), last.(curve); caustic_rad_kws...)
         end
      end
   end

   # Axis limits
   xlims!(ax, minimum(x_grid), maximum(x_grid))
   ylims!(ax, minimum(y_grid), maximum(y_grid))

   # Set axis labels and limits
   ax.xlabel = L"\theta_1~\text{(in arcseconds)}"
   ax.ylabel = L"\theta_2~\text{(in arcseconds)}"
   
   # Legend
   axislegend(ax)

   if save_plot
      save(plot_name, fig, px_per_unit=resolution)
   end
   return fig
end


function LensFactory.LensModel.plot_image_scatter(model::LensModel.ModelConfig,
                                                  chains::Array{Float64, 3},
                                                  logL::Matrix{Float64};
                                                  save_plot::Bool        = true,
                                                  plot_name::String      = "./predicted_scatter.png",
                                                  resolution::Int64      = 2,
                                                  point_kws::NamedTuple  = (markersize  = 18, 
                                                                            marker      = :circle, 
                                                                            color       = :transparent, 
                                                                            strokecolor = :black, 
                                                                            strokewidth = 2))
   # Get the best parameters based on minimum log-likelihood
   best_θ, _ = LensFactory.LensModel.get_best_fit_parameters(logL; chains=chains)
 
   # Get list of parameters for the lens model
   param_ref = Dict(p.key => p.refer for p in model.parameters)
 
   # Replace free parameter values by best-fit values
   pvals = LensFactory.LensModel.LensModelUtils.param_dict(model, best_θ, param_ref)
   
   # Update cosmology
   cosmo = LensFactory.LensModel.LensModelUtils.current_cosmology(model, pvals)
   
   # Get best-fit model
   best_model = LensFactory.LensModel.LensModelUtils.build_lens(model, pvals, cosmo)

   # Generate grid
   FOV = model.observation.FOV
   pixel_scale = model.observation.pixel_scale
   x_grid, y_grid = Lenses.get_meshgrid(0.5 * FOV[1], 0.5 * FOV[2], pixel_scale)
 
   # Get angular-diameter distance ratios
   adis = LensFactory.LensModel.adis_current(model, pvals, cosmo)
 
   # Calculate deflection and deformation tensor at all image positions
   αx_all, αy_all = LensFactory.LensModel.LensModelUtils.lens_quantities_def(model, best_model)
   A_all          = LensFactory.LensModel.LensModelUtils.lens_quantities_jac(model, best_model)

   # Identity tuple
   I4 = (1.0, 0.0, 0.0, 1.0)

   # Collect residuals across every knot/image
   dx_all = Float64[]
   dy_all = Float64[]

   sid = 1
   kid = 1
   for src in model.source_config.sources
      # Angular-diameter distance ratio for this source
      adis_value = adis[sid]
 
      for knot in src.knots
         # Knot positions and measurement errors
         x  = knot.x
         y  = knot.y
         σx = knot.σx
         σy = knot.σy
         σθ = knot.σθ

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
         β_ind = @. tuple(x - αx, y - αy)

         # Get weighted source position
         β_model, _, _ = LensModel.Likelihood._weighted_position(β_ind, A, σx, σy, σθ, n)
         βx_model, βy_model = β_model
 
         # Predicted image positions
         predicted_image = Lenses.get_image(best_model, x_grid, y_grid, adis_value, (βx_model, βy_model))

         # Convert predicted to mutable arrays for iterative removal
         pred_x = Float64[p[1] for p in predicted_image]
         pred_y = Float64[p[2] for p in predicted_image]

         # Matching observed images to predicted images
         for i in 1:n
            if isempty(pred_x)
               push!(results, "MISSING")
               continue
            end

            # Calculate distances to all remaining candidates
            dx = @. pred_x .- x[i]
            dy = @. pred_y .- y[i]
            dist_sq = @. dx^2 + dy^2

            # Find the closest predicted image index
            best_idx = argmin(dist_sq)

            d2 = dist_sq[best_idx]
            dist = sqrt(d2)

            # Residuals (observed - predicted), assumed same ordering/length as knot.x, knot.y
            push!(dx_all, x[i] - pred_x[best_idx])
            push!(dy_all, y[i] - pred_y[best_idx])

            # Remove this candidate so it can't be matched twice
            deleteat!(pred_x, best_idx)
            deleteat!(pred_y, best_idx)
         end
         kid = kid + 1
      end
      sid = sid + 1
   end

   # Initialize figure
   fig, ax = LensFactory.Lenses.plot_sky(2, 2; xlabel=L"\Delta\theta_1~\text{(arcsec)}", 
                                               ylabel=L"\Delta\theta_2~\text{(arcsec)}")

   # Residuals
   scatter!(ax, dx_all, dy_all; point_kws..., label="Residuals (obs - pred)")
  
   # Legend
   axislegend(ax)

   if save_plot
      save(plot_name, fig, px_per_unit=resolution)
   end
 
   return fig, ax
end