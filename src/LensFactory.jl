module LensFactory

# Global constants
include("./LensFactoryUtils/Constants.jl")
export Constants

# Module for cosmology
include("./Cosmology/Cosmology.jl")
export Cosmology

# Module for source profiles
include("./Lenses/Sources.jl")
export Sources

# Module for lens models
include("./Lenses/Lenses.jl")
export Lenses

# Plotting functions (see extension for more information)
export plot_image_plane
export plot_surface_density
export plot_magnification_map
export plot_magnification_profile

function plot_image_plane end
function plot_surface_density end
function plot_magnification_map end
function plot_magnification_profile end

end
