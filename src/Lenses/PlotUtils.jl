module PlotUtils


# Julia inbuilt packages to import
using GLMakie


# LensFactory modules to use
using ..Constants
using ..Lenses


# Functions to export
export plot_image_plane


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


function plot_image_plane(lens::Lenses.AbstractLens, θx::ROA, θy::ROA, adis::RV; 
                           two_panel::Bool = false,
                           plot_caustic::Bool = true,
                           caustic_kws::NamedTuple = (color_tan = :red, color_rad = :green, linewidth = 2),
                           plot_critical::Bool =true,
                           critical_kws::NamedTuple = (color_tan = :red, color_rad = :green, linewidth = 2),
                           plot_source::Bool = true,
                           plot_image::Bool = true,
                           savePlot::Bool = false,
                           plotName::String = "image_plane.png",
                           resolution::Int = 2)    
   # Initialize empty figure
   if two_panel
      fig = Figure(size=(800, 400), figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))

      # Plot source plane
      ax1 = Axis(fig[1, 1])

      # Get caustics and plot
      if plot_caustic
         # Get caustics
         caustic_tan, caustic_rad = Lenses.get_caustic(lens, θx, θy, adis)

         # Plot tangential caustic
         for curve in caustic_tan
            lines!(ax1, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
                                       color=caustic_kws.color_tan, 
                                       linewidth=caustic_kws.linewidth)
         end

         # Plot radial caustic
         for curve in caustic_rad
            lines!(ax1, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
                                       color=caustic_kws.color_rad, 
                                       linewidth=caustic_kws.linewidth)
         end
      end
   
      
      # Set plot keywords
      set_plotKws!(ax1)

      # Set axis labels and limits
      ax1.xlabel = L"\theta_1 \text{(in arcseconds)}"
      ax1.ylabel = L"\theta_2 \text{(in arcseconds)}"
      xlims!(minimum(θx) / ANGLE_ARCSEC, maximum(θx) / ANGLE_ARCSEC)
      ylims!(minimum(θy) / ANGLE_ARCSEC, maximum(θy) / ANGLE_ARCSEC)

      # Plot image plane
      ax2 = Axis(fig[1, 2])

      # Get critical curves
      if plot_critical
         # Get critical curves
         crit_tan, crit_rad = Lenses.get_critical_curve(lens, θx, θy, adis)

         # Plot tangential critical curve
         for curve in crit_tan
            lines!(ax2, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
                                       color=critical_kws.color_tan, 
                                       linewidth=critical_kws.linewidth)
         end

         # Plot radial critical curve
         for curve in crit_rad
            lines!(ax2, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
                                       color=critical_kws.color_rad, 
                                       linewidth=critical_kws.linewidth)
         end
      end

      # Set plot keywords
      set_plotKws!(ax2)

      # Set axis labels and limits
      ax2.xlabel = L"\theta_1 \text{(in arcseconds)}"
      ax2.ylabel = L"\theta_2 \text{(in arcseconds)}"
      xlims!(minimum(θx) / ANGLE_ARCSEC, maximum(θx) / ANGLE_ARCSEC)
      ylims!(minimum(θy) / ANGLE_ARCSEC, maximum(θy) / ANGLE_ARCSEC)

      # Save plot
      if savePlot
         save(plotName, fig, px_per_unit=resolution)
      end
      return fig, [ax1, ax2]
   else
      fig = Figure(size=(400, 400), figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))
      
      # Plot source + image plane
      ax = Axis(fig[1, 1])

      if plot_caustic
         # Get caustics
         caustic_tan, caustic_rad = Lenses.get_caustic(lens, θx, θy, adis)

         # Plot tangential caustic
         for curve in caustic_tan
            lines!(ax, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
                                       color=caustic_kws.color_tan, 
                                       linewidth=caustic_kws.linewidth)
         end

         # Plot radial caustic
         for curve in caustic_rad
            lines!(ax, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
                                       color=caustic_kws.color_rad, 
                                       linewidth=caustic_kws.linewidth)
         end
      end

      if plot_critical
         # Get critical curves
         crit_tan, crit_rad = Lenses.get_critical_curve(lens, θx, θy, adis)

         # Plot tangential critical curve
         for curve in crit_tan
            lines!(ax, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
                                       color=critical_kws.color_tan, 
                                       linewidth=critical_kws.linewidth)
         end

         # Plot radial critical curve
         for curve in crit_rad
            lines!(ax, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
                                       color=critical_kws.color_rad, 
                                       linewidth=critical_kws.linewidth)
         end
      end

      # Set plot keywords
      set_plotKws!(ax)

      # Set axis labels and limits
      ax.xlabel = L"\theta_1 \text{(in arcseconds)}"
      ax.ylabel = L"\theta_2 \text{(in arcseconds)}"
      xlims!(minimum(θx) / ANGLE_ARCSEC, maximum(θx) / ANGLE_ARCSEC)
      ylims!(minimum(θy) / ANGLE_ARCSEC, maximum(θy) / ANGLE_ARCSEC)

      if savePlot
         save(plotName, fig, px_per_unit=resolution)
      end
      return fig, ax
   end

end


end