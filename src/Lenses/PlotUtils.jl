# module PlotUtils


# # Julia inbuilt packages to import
# using GLMakie


# # LensFactory modules to use
# using ..Constants
# using ..Lenses


# # Functions to export
# plot_image_plane


# function set_plotKws!(ax)
#    ax.xtickalign = 1
#    ax.xticksmirrored = true
#    ax.ytickalign = 1
#    ax.yticksmirrored = true

#    ax.xminorticksvisible = true
#    ax.xminortickalign = 1
#    ax.xminorticksize = 6
#    ax.xminorgridwidth = 2
#    ax.yminorticksvisible = true
#    ax.yminortickalign = 1
#    ax.yminorticksize = 6
#    ax.yminorgridwidth = 2

#    ax.xticksize = 10
#    ax.xtickwidth = 2
#    ax.yticksize = 10
#    ax.ytickwidth = 2
# end


# function plot_image_plane(lens::AbstractLens, θx::ROA, θy::ROA, adis::RV; 
#                            two_panel::Bool = false,
#                            plot_caustic::Bool = true,
#                            caustic_kws::NamedTuple = (color_tan = :red, color_rad = :green, linewidth = 2),
#                            plot_critical::Bool =true,
#                            critical_kws::NamedTuple = (color_tan = :red, color_rad = :green, linewidth = 2),
#                            plot_source::Bool = true,
#                            plot_image::Bool = true,
#                            sourceColor = (0, 0, 0),
#                            imageColor = (0, 0, 0),
#                            savePlot::Bool = false,
#                            plotName::String = "image_plane.png",
#                            resolution::Int = 2)    
#    # Initialize empty figure
#    if two_panel
#       fig = Figure(size=(800, 400), figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))

#       # Plot source plane
#       ax1 = Axis(fig[1, 1])

#       set_plotKws!(ax1)
#       ax1.xlabel = L"\theta_1 \text{(in arcseconds)}"
#       ax1.ylabel = L"\theta_2 \text{(in arcseconds)}"
#       xlims!(minimum(θx) / ANGLE_ARCSEC, maximum(θx) / ANGLE_ARCSEC)
#       ylims!(minimum(θy) / ANGLE_ARCSEC, maximum(θy) / ANGLE_ARCSEC)

#       if savePlot
#          save(plotName, fig, px_per_unit=resolution)
#       end
#       return fig, [ax1, ax2]
#    else
#       fig = Figure(size=(400, 400), figure_padding=15, fontsize=20, fonts=(; regular="Times New Roman"))
      
#       # Plot source + image plane
#       ax = Axis(fig[1, 1])

#       if savePlot
#          save(plotName, fig, px_per_unit=resolution)
#       end
#       return fig, ax
#    end

# end


# end