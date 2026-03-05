"""
    LensFactory.MultiPlane.plot_image_plane(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::ROA, θy::ROA, zs::RV)

# Arguments
- `cosmology::Cosmology.AbstractCosmology` -- Cosmology
- `lens::Lenses.AbstractLens` -- The lens model to plot
- `θx::ROA` -- x-coordinates grid
- `θy::ROA` -- y-coordinates grid
- `zs::RV` -- Source redshift

# Keyword arguments
- `two_panel::Bool = false` -- Whether to create a two-panel plot with source plane on the left and image plane on the right
- `plot_caustic::Bool = true`
   - `caustic_kws::NamedTuple = (color = :green, linewidth = 2)`
- `plot_critical::Bool =true`
   - `critical_kws::NamedTuple = (color = :red, linewidth = 2)`
- `source::Union{Nothing, NTuple{2, RV}, Matrix{<:RV}} = nothing`
   - `source_kws::NamedTuple = (color=:red, markersize=10, marker=:star5, heatmap=cgrad([:white, :blue]))`
   - `image_kws::NamedTuple = (color=:blue, markersize=10, marker=:star5, heatmap=cgrad([:white, :red]))`
- `save_plot::Bool = false`
   - `plot_name::String = "image_plane.png"`
   - `resolution::Int = 2`

# Returns
- `fig`: A Makie figure object containing the image plane plot.
- `ax`: The axis object of the plot for further customization.
"""
function LensFactory.MultiPlane.plot_image_plane(cosmology::Cosmology.AbstractCosmology, 
                           lens::Lenses.AbstractLens, θx::Matrix{<:RV}, θy::Matrix{<:RV}, zs::RV;
                           two_panel::Bool = false,
                           plot_caustic::Bool = true,
                           caustic_kws::NamedTuple = (color = :green, linewidth = 2),
                           plot_critical::Bool =true,
                           critical_kws::NamedTuple = (color = :red, linewidth = 2),
                           source::Union{Nothing, NTuple{2, RV}, Matrix{<:RV}} = nothing,
                           source_kws::NamedTuple = (color=:red, markersize=10, marker=:star5, heatmap=cgrad([:white, :blue])),
                           image_kws::NamedTuple = (color=:blue, markersize=10, marker=:star5, heatmap=cgrad([:white, :red])),
                           save_plot::Bool = false,
                           plot_name::String = "image_plane.png",
                           resolution::Int = 2)
   
   if two_panel
      # Initialize empty figure
      fig = Figure(size=(800, 400), figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))

      # Axis for source plane
      ax1 = Axis(fig[1, 1])

      # Plot source and its images
      if source !== nothing
         if isa(source, NTuple{2, RV})
            scatter!(ax1, source[1], source[2], color=source_kws.color, markersize=source_kws.markersize, marker=source_kws.marker)
         elseif isa(source, Matrix{<:RV})
            heatmap!(ax1, θx[:,1], θy[1,:], source, colormap=source_kws.heatmap)
         else
            error("Invalid source type: $(typeof(source)). Must be NTuple{2, RV} or Matrix{<:RV}.")
         end
      end

      # Get caustics and plot
      if plot_caustic
         # Get caustics
         caustic_curve = MultiPlane.get_caustic(cosmology, lens, θx, θy, zs)

         # Plot tangential caustic
         for curve in caustic_curve
            lines!(ax1, first.(curve), last.(curve); caustic_kws...)
         end
      end

      # Set plot keywords
      set_plotKws!(ax1)

      # Set axis labels and limits
      ax1.xlabel = L"\theta_1~\text{(in arcseconds)}"
      ax1.ylabel = L"\theta_2~\text{(in arcseconds)}"
      xlims!(minimum(θx), maximum(θx))
      ylims!(minimum(θy), maximum(θy))

      # Axis for image plane
      ax2 = Axis(fig[1, 2])

      # Plot source and its images
      if source !== nothing
         # Get the image positions
         image = MultiPlane.get_image(cosmology, lens, θx, θy, zs, source)
         if isa(source, NTuple{2, RV})
            scatter!(ax2, first.(image), last.(image); color=image_kws.color, markersize=image_kws.markersize, marker=image_kws.marker)
         elseif isa(source, Matrix{<:RV})
            heatmap!(ax2, θx[:,1], θy[1,:], image, colormap=image_kws.heatmap)
         else
            ArgumentError("Invalid source type: $(typeof(source)). Must be NTuple{2, RV} or Matrix{<:RV}.")
         end
      end

      # Get critical curves
      if plot_critical
         # Get critical curves
         criical_curve = MultiPlane.get_critical_curve(cosmology, lens, θx, θy, zs)

         # Plot tangential critical curve
         for curve in criical_curve
            lines!(ax2, first.(curve), last.(curve); critical_kws...)
         end
      end

      # Set plot keywords
      set_plotKws!(ax2)

      # Set axis labels and limits
      ax2.xlabel = L"\theta_1~\text{(in arcseconds)}"
      ax2.ylabel = L"\theta_2~\text{(in arcseconds)}"
      xlims!(minimum(θx), maximum(θx))
      ylims!(minimum(θy), maximum(θy))

      # Save plot
      if save_plot
         save(plot_name, fig, px_per_unit=resolution)
      end
      return fig, [ax1, ax2]
   else
      # Initialize empty figure
      fig = Figure(size=(400, 400), figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))
      
      # Plot source + image plane
      ax = Axis(fig[1, 1])

      if source !== nothing
         # Get the image positions
         image = MultiPlane.get_image(cosmology, lens, θx, θy, zs, source)

         if isa(source, NTuple{2, RV})
            scatter!(ax, source[1], source[2], color=source_kws.color, markersize=source_kws.markersize, marker=source_kws.marker)
            scatter!(ax, first.(image), last.(image), color=image_kws.color, markersize=image_kws.markersize, marker=image_kws.marker)
         elseif isa(source, Matrix{<:RV})
            heatmap!(ax, θx[:,1], θy[1,:], source, colormap=source_kws.heatmap, alpha=1.0)                                    
            heatmap!(ax, θx[:,1], θy[1,:], image, colormap=image_kws.heatmap, alpha=0.8)
         else
            error("Invalid source type: $(typeof(source)). Must be NTuple{2, RV} or Matrix{<:RV}.")
         end
      end

      # Get caustics and plot
      if plot_caustic
         # Get caustics
         caustic_curve = MultiPlane.get_caustic(cosmology, lens, θx, θy, zs)

         # Plot tangential caustic
         for curve in caustic_curve
            lines!(ax, first.(curve), last.(curve), color=caustic_kws.color, linewidth=caustic_kws.linewidth, linestyle=:solid)
         end
      end

      if plot_critical
         # Get critical curves
         critical_curve = MultiPlane.get_critical_curve(cosmology, lens, θx, θy, zs)

         # Plot tangential critical curve
         for curve in critical_curve
            lines!(ax, first.(curve), last.(curve), color=critical_kws.color, linewidth=critical_kws.linewidth, linestyle=:solid)
         end
      end
      # Set plot keywords
      set_plotKws!(ax)

      # Set axis labels and limits
      ax.xlabel = L"\theta_1~\text{(in arcseconds)}"
      ax.ylabel = L"\theta_2~\text{(in arcseconds)}"
      xlims!(minimum(θx), maximum(θx))
      ylims!(minimum(θy), maximum(θy))

      if save_plot
         save(plot_name, fig, px_per_unit=resolution)
      end
      return fig, ax
   end
end
