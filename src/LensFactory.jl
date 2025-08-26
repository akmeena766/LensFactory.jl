"""
    LensFactory

Testing :-)

"""
module LensFactory


# export cpu_vs_gpu
# """
#     cpu_vs_gpu::Symbol = :cpu

# A symbol that controls whether computations run on CPU or GPU.
# """
# cpu_vs_gpu::Symbol = :cpu


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


export plot_image_plane
export plot_surface_density
export plot_magnification_map

function plot_image_plane end
function plot_surface_density end
function plot_magnification_map end


end
