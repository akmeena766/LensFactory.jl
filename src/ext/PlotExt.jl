module PlotExt

using LensFactory
using Makie


# Functions to export
# export plot_image_plane
# export plot_surface_density


# function plot_image_plane(lens::Lenses.AbstractLens, θx::ROA, θy::ROA, adis::RV; 
#                            two_panel::Bool = false,
#                            plot_caustic::Bool = true,
#                            caustic_kws::NamedTuple = (color_tan = :red, color_rad = :green, linewidth = 2),
#                            plot_critical::Bool =true,
#                            critical_kws::NamedTuple = (color_tan = :red, color_rad = :green, linewidth = 2),
#                            plot_source::Bool = true,
#                            plot_image::Bool = true,
#                            save_plot::Bool = false,
#                            plot_name::String = "image_plane.png",
#                            resolution::Int = 2)

#    if two_panel
#       # Initialize empty figure
#       fig = Figure(size=(800, 400), figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))

#       # Plot source plane
#       ax1 = Axis(fig[1, 1])

#       # Get caustics and plot
#       if plot_caustic
#          # Get caustics
#          caustic_tan, caustic_rad = Lenses.get_caustic(lens, θx, θy, adis)

#          # Plot tangential caustic
#          for curve in caustic_tan
#             lines!(ax1, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
#                                        color=caustic_kws.color_tan, 
#                                        linewidth=caustic_kws.linewidth)
#          end

#          # Plot radial caustic
#          for curve in caustic_rad
#             lines!(ax1, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
#                                        color=caustic_kws.color_rad, 
#                                        linewidth=caustic_kws.linewidth)
#          end
#       end
   
      
#       # Set plot keywords
#       set_plotKws!(ax1)

#       # Set axis labels and limits
#       ax1.xlabel = L"\theta_1 \text{(in arcseconds)}"
#       ax1.ylabel = L"\theta_2 \text{(in arcseconds)}"
#       xlims!(minimum(θx) / ANGLE_ARCSEC, maximum(θx) / ANGLE_ARCSEC)
#       ylims!(minimum(θy) / ANGLE_ARCSEC, maximum(θy) / ANGLE_ARCSEC)

#       # Plot image plane
#       ax2 = Axis(fig[1, 2])

#       # Get critical curves
#       if plot_critical
#          # Get critical curves
#          crit_tan, crit_rad = Lenses.get_critical_curve(lens, θx, θy, adis)

#          # Plot tangential critical curve
#          for curve in crit_tan
#             lines!(ax2, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
#                                        color=critical_kws.color_tan, 
#                                        linewidth=critical_kws.linewidth)
#          end

#          # Plot radial critical curve
#          for curve in crit_rad
#             lines!(ax2, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
#                                        color=critical_kws.color_rad, 
#                                        linewidth=critical_kws.linewidth)
#          end
#       end

#       # Set plot keywords
#       set_plotKws!(ax2)

#       # Set axis labels and limits
#       ax2.xlabel = L"\theta_1 \text{(in arcseconds)}"
#       ax2.ylabel = L"\theta_2 \text{(in arcseconds)}"
#       xlims!(minimum(θx) / ANGLE_ARCSEC, maximum(θx) / ANGLE_ARCSEC)
#       ylims!(minimum(θy) / ANGLE_ARCSEC, maximum(θy) / ANGLE_ARCSEC)

#       # Save plot
#       if save_plot
#          save(plot_name, fig, px_per_unit=resolution)
#       end
#       return fig, [ax1, ax2]
#    else
#       # Initialize empty figure
#       fig = Figure(size=(400, 400), figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))
      
#       # Plot source + image plane
#       ax = Axis(fig[1, 1])

#       if plot_caustic
#          # Get caustics
#          caustic_tan, caustic_rad = Lenses.get_caustic(lens, θx, θy, adis)

#          # Plot tangential caustic
#          for curve in caustic_tan
#             lines!(ax, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
#                                        color=caustic_kws.color_tan, 
#                                        linewidth=caustic_kws.linewidth)
#          end

#          # Plot radial caustic
#          for curve in caustic_rad
#             lines!(ax, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
#                                        color=caustic_kws.color_rad, 
#                                        linewidth=caustic_kws.linewidth)
#          end
#       end

#       if plot_critical
#          # Get critical curves
#          crit_tan, crit_rad = Lenses.get_critical_curve(lens, θx, θy, adis)

#          # Plot tangential critical curve
#          for curve in crit_tan
#             lines!(ax, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
#                                        color=critical_kws.color_tan, 
#                                        linewidth=critical_kws.linewidth)
#          end

#          # Plot radial critical curve
#          for curve in crit_rad
#             lines!(ax, first.(curve) ./ ANGLE_ARCSEC, last.(curve) ./ ANGLE_ARCSEC, 
#                                        color=critical_kws.color_rad, 
#                                        linewidth=critical_kws.linewidth)
#          end
#       end

#       # Set plot keywords
#       set_plotKws!(ax)

#       # Set axis labels and limits
#       ax.xlabel = L"\theta_1 \text{(in arcseconds)}"
#       ax.ylabel = L"\theta_2 \text{(in arcseconds)}"
#       xlims!(minimum(θx) / ANGLE_ARCSEC, maximum(θx) / ANGLE_ARCSEC)
#       ylims!(minimum(θy) / ANGLE_ARCSEC, maximum(θy) / ANGLE_ARCSEC)

#       if save_plot
#          save(plot_name, fig, px_per_unit=resolution)
#       end
#       return fig, ax
#    end
# end


# function plot_surface_density(lens::Lenses.AbstractLens, θx::ROA, θy::ROA; 
#                               D_d::RV=NaN, adis::RV=NaN,
#                               unit::Symbol = :kg_m2,
#                               figure_size::Tuple{Int, Int} = (500, 400),
#                               plot_contour::Bool = false,
#                               save_plot::Bool = false,
#                               plot_name::String = "surface_density.png",
#                               resolution::Int = 2)
#    # Get jacobian and rescale according to adis
#    ψxx, ψyy, _ = Lenses.get_jacobian(lens, θx, θy)
#    ψxx .*= adis
#    ψyy .*= adis

#    # Get critical density
#    Σ_cr = Lenses.get_critical_density(D_d=D_d, adis=adis, unit=unit)

#    # Get surface density
#    Σ = 0.5 .* (ψxx .+ ψyy) .* Σ_cr

#    # Initialize empty figure
#    fig = Figure(size=figure_size, figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))

#    # Plot source + image plane
#    ax = Axis(fig[1, 1])

#    # Plot surface density
#    hm = heatmap!(ax, θx[:,1] ./ ANGLE_ARCSEC, θy[1,:] ./ ANGLE_ARCSEC, log10.(Σ);
#                                     colormap=:Greys, 
#                                     colorrange=(log10.(minimum(Σ)), 6))

#    # Colorbar specification
#    cb = Colorbar(fig[1, 2], hm; label=L"\text{Log_{10} Σ}", 
#                                     labelpadding=5, 
#                                     width=20, 
#                                     tickalign=1, 
#                                     ticksize=10, 
#                                     tickwidth=1.5, 
#                                     labelrotation=3*pi/2)
#    if unit == :kg_m2
#       cb.label = L"\text{Log_{10} Σ (in kg/m}^2\text{)}"
#    elseif unit == :msun_pc2
#       cb.label = L"\text{Log_{10} Σ (in M}_{\odot}\text{/pc}^2\text{)}"
#    elseif unit == :msun_arcsec2
#       cb.label = L"\text{Log_{10} Σ (in M}_{\odot}\text{/arcsec}^2\text{)}"
#    else
#       error("Invalid unit: $unit. Must be :kg_m2 or :msun_pc2 or :msun_arcsec2.")
#    end

#    # Plot contours
#    if plot_contour
#       levels = collect(0.5:0.5:2)
#       contour!(ax, θx[:,1] ./ ANGLE_ARCSEC, θy[1,:] ./ ANGLE_ARCSEC, Σ ./ Σ_cr; 
#                                     levels=levels, 
#                                     labels=true, 
#                                     color=:black,
#                                     labelsize = 15)
#    end

#    # Set plot keywords
#    set_plotKws!(ax)

#    # Set axis labels and limits
#    ax.xlabel = L"\theta_1 \text{(in arcseconds)}"
#    ax.ylabel = L"\theta_2 \text{(in arcseconds)}"
#    xlims!(minimum(θx) / ANGLE_ARCSEC, maximum(θx) / ANGLE_ARCSEC)
#    ylims!(minimum(θy) / ANGLE_ARCSEC, maximum(θy) / ANGLE_ARCSEC)

#    if save_plot
#       save(plot_name, fig, px_per_unit=resolution)
#    end
#    return fig, ax
# end

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