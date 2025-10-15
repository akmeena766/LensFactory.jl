module PlotExt

using LensFactory
using LensFactory.Constants
using Makie


"""
    LensFactory.Lenses.plot_image_plane(lens::Lenses.AbstractLens, θx::ROA, θy::ROA, adis::RV)

# Additional keyword arguments and their default values
- `two_panel::Bool = false` -- Whether to create a two-panel plot with source plane on the left and image plane on the right
- `plot_caustic::Bool = true`
   - `caustic_kws::NamedTuple = (color_tan = :green, color_rad = :green, linewidth = 2)`
- `plot_critical::Bool =true`
   - `critical_kws::NamedTuple = (color_tan = :red, color_rad = :red, linewidth = 2)`
- `source::Union{Nothing, NTuple{2, RV}, Matrix{<:RV}} = nothing`
   - `source_kws::NamedTuple = (color=:red, markersize=10, marker=:star5, heatmap=cgrad([:white, :blue]))`
   - `image_kws::NamedTuple = (color=:blue, markersize=10, marker=:star5, heatmap=cgrad([:white, :red]))`
- `save_plot::Bool = false`
   - `plot_name::String = "image_plane.png"`
   - `resolution::Int = 2`
"""
function LensFactory.Lenses.plot_image_plane(lens::Lenses.AbstractLens, θx::ROA, θy::ROA, adis::Float64;
                           two_panel::Bool = false,
                           plot_caustic::Bool = true,
                           caustic_kws::NamedTuple = (color_tan = :green, color_rad = :green, linewidth = 2),
                           plot_critical::Bool =true,
                           critical_kws::NamedTuple = (color_tan = :red, color_rad = :red, linewidth = 2),
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
         caustic_tan, caustic_rad = Lenses.get_caustic(lens, θx, θy, adis)

         # Plot tangential caustic
         for curve in caustic_tan
            lines!(ax1, first.(curve), last.(curve), color=caustic_kws.color_tan, linewidth=caustic_kws.linewidth,linestyle=:solid)
         end

         # Plot radial caustic
         for curve in caustic_rad
            lines!(ax1, first.(curve), last.(curve), color=caustic_kws.color_rad, linewidth=caustic_kws.linewidth, linestyle=:dash)
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
         image = Lenses.get_image(lens, θx, θy, adis, source)
         if isa(source, NTuple{2, RV})
            scatter!(ax2, first.(image), last.(image), color=image_kws.color, markersize=image_kws.markersize, marker=image_kws.marker)
         elseif isa(source, Matrix{<:RV})
            heatmap!(ax2, θx[:,1], θy[1,:], image, colormap=image_kws.heatmap)
         else
            ArgumentError("Invalid source type: $(typeof(source)). Must be NTuple{2, RV} or Matrix{<:RV}.")
         end
      end

      # Get critical curves
      if plot_critical
         # Get critical curves
         crit_tan, crit_rad = Lenses.get_critical_curve(lens, θx, θy, adis)

         # Plot tangential critical curve
         for curve in crit_tan
            lines!(ax2, first.(curve), last.(curve), color=critical_kws.color_tan, linewidth=critical_kws.linewidth, linestyle=:solid)
         end

         # Plot radial critical curve
         for curve in crit_rad
            lines!(ax2, first.(curve), last.(curve), color=critical_kws.color_rad, linewidth=critical_kws.linewidth, linestyle=:dash)
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
         image = Lenses.get_image(lens, θx, θy, adis, source)

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

      if plot_caustic
         # Get caustics
         caustic_tan, caustic_rad = Lenses.get_caustic(lens, θx, θy, adis)

         # Plot tangential caustic
         for curve in caustic_tan
            lines!(ax, first.(curve), last.(curve), color=caustic_kws.color_tan, linewidth=caustic_kws.linewidth, linestyle=:solid)
         end

         # Plot radial caustic
         for curve in caustic_rad
            lines!(ax, first.(curve), last.(curve), color=caustic_kws.color_rad, linewidth=caustic_kws.linewidth, linestyle=:dash)
         end
      end

      if plot_critical
         # Get critical curves
         crit_tan, crit_rad = Lenses.get_critical_curve(lens, θx, θy, adis)

         # Plot tangential critical curve
         for curve in crit_tan
            lines!(ax, first.(curve), last.(curve), color=critical_kws.color_tan, linewidth=critical_kws.linewidth, linestyle=:solid)
         end

         # Plot radial critical curve
         for curve in crit_rad
            lines!(ax, first.(curve), last.(curve), color=critical_kws.color_rad, linewidth=critical_kws.linewidth, linestyle=:dash)
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

"""
    LensFactory.Lenses.plot_surface_density(lens::Lenses.AbstractLens, θx::ROA, θy::ROA, adis::Float64; D_d::Float6=NaN)

# Additional keyword arguments and their default values
- `unit::Symbol = :kg_m2` -- Unit for surface density. Options are :convergence, :kg_m2, :msun_pc2, :msun_arcsec2
- `figure_size::NTuple{2, RV} = (500, 400)`
- `heatmap_kws::NamedTuple = (colormap=:cubehelix, colorrange=(0, 6))`
- `plot_contour::Bool = false`
   - `contour_kws::NamedTuple = (levels=0.5:0.2:1.5, labels=false)`
- `save_plot::Bool = false`
   - `plot_name::String = "surface_density.png"`
   - `resolution::Int = 2`
"""
function LensFactory.Lenses.plot_surface_density(lens::Lenses.AbstractLens, θx::ROA, θy::ROA, adis::Float64; 
                              D_d::Float64 = NaN,
                              unit::Symbol = :kg_m2,
                              figure_size::NTuple{2, RV} = (500, 400),
                              heatmap_kws::NamedTuple = (colormap=:cubehelix, colorrange=(0, 6)),
                              plot_contour::Bool = false,
                              contour_kws::NamedTuple = (levels=0.5:0.2:1.5, labels=false),
                              save_plot::Bool = false,
                              plot_name::String = "surface_density.png",
                              resolution::Int = 2)
   # Get jacobian and rescale according to adis
   ψxx, ψyy, _ = Lenses.get_jacobian(lens, θx, θy)
   ψxx .*= adis
   ψyy .*= adis

   if unit == :convergence
      cb_label = L"\kappa"
   elseif unit == :kg_m2
      cb_label = L"\text{Log_{10} Σ (in kg/m}^2\text{)}"
   elseif unit == :msun_pc2
      cb_label = L"\text{Log_{10} Σ (in M}_{\odot}\text{/pc}^2\text{)}"
   elseif unit == :msun_arcsec2
      cb_label = L"\text{Log_{10} Σ (in M}_{\odot}\text{/arcsec}^2\text{)}"
   else
      ArgumentError("Invalid unit: $unit. Must be :convergence or :kg_m2 or :msun_pc2 or :msun_arcsec2.")
   end

   # Get critical density
   if unit == :convergence
      Σ_cr = 1.0
   else
      Σ_cr = Lenses.get_critical_density(D_d=D_d, adis=adis, unit=unit)
   end

   # Get surface density
   Σ = 0.5 .* (ψxx .+ ψyy) .* Σ_cr

   # Initialize empty figure
   fig = Figure(size=figure_size, figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))

   # Plot source + image plane
   ax = Axis(fig[1, 1])

   # Plot surface density
   if unit == :convergence
      hm = heatmap!(ax, θx[:,1], θy[1,:], Σ; heatmap_kws...)
   else
      hm = heatmap!(ax, θx[:,1], θy[1,:], log10.(Σ); heatmap_kws...)
   end

   # Colorbar specification
   cb = Colorbar(fig[1, 2], hm; label=cb_label, labelpadding=5, width=20, tickalign=1, ticksize=10, tickwidth=1.5, labelrotation=3*pi/2)
   colgap!(fig.layout, 5)
   
   # Plot contours
   if plot_contour
      contour!(ax, θx[:,1], θy[1,:], Σ; contour_kws...)
   end

   # Set plot keywords
   set_plotKws!(ax)

   # Set axis labels and limits
   ax.xlabel = L"\theta_1 \text{(in arcseconds)}"
   ax.ylabel = L"\theta_2 \text{(in arcseconds)}"
   xlims!(minimum(θx), maximum(θx))
   ylims!(minimum(θy), maximum(θy))

   if save_plot
      save(plot_name, fig, px_per_unit=resolution)
   end
   return fig, ax
end


"""
    LensFactory.Lenses.plot_magnification_map(lens::Lenses.AbstractLens, θx::ROA, θy::ROA, adis::Float64)

# Additional keyword arguments and their default values
- `plane::Symbol = :image` -- Whether to plot |μ| map in `:image` or `:source` plane
- `rays_per_pixel::Int64 = 1` -- Only relevant if `plane = :source`
- `figure_size::NTuple{2, RV} = (500, 400)`
- `heatmap_kws::NamedTuple = (colormap=:binary, colorrange=(1, 100))`
- `save_plot::Bool = false`
   - `plot_name::String = "magnification_map.png"`
   - `resolution::Int = 2`
"""
function LensFactory.Lenses.plot_magnification_map(lens::Lenses.AbstractLens, θx::ROA, θy::ROA, adis::Float64;
                              plane::Symbol = :image,
                              rays_per_pixel::Int64 = 1,
                              figure_size::NTuple{2, RV} = (500, 400),
                              heatmap_kws::NamedTuple = (colormap=:binary, colorrange=(1, 100)),
                              save_plot::Bool = false,
                              plot_name::String = "magnification_map.png",
                              resolution::Int = 2)
   # Get the magnification map
   if plane == :image
      μ = Lenses.get_magnification_image(lens, θx, θy, adis)
   elseif plane == :source
      μ = Lenses.get_magnification_source(lens, θx, θy, adis, rays_per_pixel=rays_per_pixel)
   else
      ArgumentError("Invalid plane: $plane. Must be :image or :source.")
   end

   # Initialize empty figure
   fig = Figure(size=figure_size, figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))

   # Axis for the plot
   ax = Axis(fig[1, 1])

   # Plot the magnification map
   hm = heatmap!(ax, θx[:,1], θy[1,:], abs.(μ); heatmap_kws...)

   # Colorbar specification
   cb = Colorbar(fig[1, 2], hm; label=L"|μ|", labelpadding=5, width=20, tickalign=1, ticksize=10, tickwidth=1.5, labelrotation=3*pi/2)
   colgap!(fig.layout, 5)

   # Set plot keywords
   set_plotKws!(ax)

   # Set axis labels and limits
   ax.xlabel = L"\theta_1 \text{(in arcseconds)}"
   ax.ylabel = L"\theta_2 \text{(in arcseconds)}"
   xlims!(minimum(θx), maximum(θx))
   ylims!(minimum(θy), maximum(θy))

   if save_plot
      save(plot_name, fig, px_per_unit=resolution)
   end
   return fig, ax
end


"""
    LensFactory.Lenses.plot_magnification_profile(lens::Lenses.AbstractLens, θx::ROA, θy::ROA, adis::Float64)

# Additional keyword arguments and their default values
- `plane::Symbol = :image` -- Whether to plot |μ| map in `:image` or `:source` plane
- `rays_per_pixel::Int64 = 1` -- Only relevant if `plane = :source`
- `mu_range::StepRange{<:RV, <:RV} = 1:5:500`
- `unit::Float64 = ANGLE_ARCMIN^2`
- `figure_size::NTuple{2, RV} = (500, 400)`
- `plot_kws::NamedTuple = (color=:black, linewidth=2, linestyle=:solid)`
- `save_plot::Bool = false`
   - `plot_name::String = "magnification_profile.png"`
   - `resolution::Int = 2`
"""
function LensFactory.Lenses.plot_magnification_profile(lens::Lenses.AbstractLens, θx::ROA, θy::ROA, adis::Float64;
                              plane::Symbol = :image,
                              rays_per_pixel::Int64 = 1,
                              mu_range::StepRange{<:RV, <:RV} = 1:5:500,
                              unit::Float64 = ANGLE_ARCMIN^2,
                              figure_size::NTuple{2, RV} = (500, 400),
                              plot_kws::NamedTuple = (color=:black, linewidth=2, linestyle=:solid),
                              save_plot::Bool = false,
                              plot_name::String = "magnification_profile.png",
                              resolution::Int = 2)
   # Get the magnification map
   if plane == :image
      μ = Lenses.get_magnification_image(lens, θx, θy, adis)
   elseif plane == :source
      μ = Lenses.get_magnification_source(lens, θx, θy, adis, rays_per_pixel=rays_per_pixel)
   else
      ArgumentError("Invalid plane: $plane. Must be :image or :source.")
   end

   # Pixel size
   pixel_h::Float64 = abs(θx[2, 1] - θx[1, 1])

   # Get the magnification bins
   μ_bins = collect(mu_range)

   # Flatten the magnification map and get the area for each bin
   μ_flatten = abs.(vec(μ))
   μ_area = [sum(μ_flatten .>= t) * pixel_h^2 for t in μ_bins]

   # Convert to proper units
   μ_area .*= (ANGLE_ARCSEC^2 / unit)

   # Initialize empty figure
   fig = Figure(size=figure_size, figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))

   # Axis for the plot
   ax = Axis(fig[1, 1])

   lines!(ax, μ_bins, μ_area; plot_kws...)

   # Set plot keywords
   set_plotKws!(ax)

   # Set axis labels and limits
   ax.xlabel = L"|\mu|"
   ax.ylabel = L"\text{Area ≥ |μ| [arcmin^2]}"
   xlims!(minimum(μ_bins), maximum(μ_bins))
   ylims!(1E-5, maximum(μ_area)*2)

   ax.xscale = log10
   ax.yscale = log10

   if save_plot
      save(plot_name, fig, px_per_unit=resolution)
   end
   return fig, ax
end

"""
    LensFactory.MultiPlane.plot_image_plane(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::ROA, θy::ROA, zs::RV)

# Additional keyword arguments and their default values
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
"""
function LensFactory.MultiPlane.plot_image_plane(cosmology::Cosmology.AbstractCosmology, lens::Lenses.AbstractLens, θx::ROA, θy::ROA, zs::RV;
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


function set_plotKws!(ax)
   ax.xtickalign = 1
   ax.xticksmirrored = true
   ax.ytickalign = 1
   ax.yticksmirrored = true

   ax.xminorticksvisible = true
   ax.xminortickalign = 1
   ax.xminorticksize = 6
   ax.xminorgridwidth = 2
   ax.yminorticksvisible = true
   ax.yminortickalign = 1
   ax.yminorticksize = 6
   ax.yminorgridwidth = 2

   ax.xticksize = 10
   ax.xtickwidth = 2
   ax.yticksize = 10
   ax.ytickwidth = 2
end


end