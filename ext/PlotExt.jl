module PlotExt


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using Makie
using StatsBase
using KernelDensity

# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using LensFactory
using LensFactory.Constants
using LensFactory.LensModel


# --------------------------------------------------------------------------------------------------
# Function to set the plot keywords for all the plots in this module
# --------------------------------------------------------------------------------------------------
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


# --------------------------------------------------------------------------------------------------
# Plot functions for basic single-plane lensing
# --------------------------------------------------------------------------------------------------
include("./PlotLenses.jl")


# --------------------------------------------------------------------------------------------------
# Plot functions for multi-plane lensing
# --------------------------------------------------------------------------------------------------
include("./PlotMultiPlane.jl")


# --------------------------------------------------------------------------------------------------
# Plot functions for lens modelling
# --------------------------------------------------------------------------------------------------
include("./PlotLensModel.jl")


end