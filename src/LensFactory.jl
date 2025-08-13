"""
    LensFactory

Testing :-)

"""
module LensFactory


export cpu_vs_gpu
"""
    cpu_vs_gpu::Symbol = :cpu

A symbol that controls whether computations run on CPU or GPU.
"""
cpu_vs_gpu::Symbol = :cpu


# Global constants
include("./LensFactoryUtils/Constants.jl")
export Constants

# Modules for contour finding on a 2D grid
include("./LensFactoryUtils/ContourFinder.jl")
export ContourFinder

# Module to get intersection points of two contours
include("./LensFactoryUtils/IntersectionFinder.jl")
export IntersectionFinder

# Module for cosmology
include("./Cosmology/Cosmology.jl")
export Cosmology

# Module for source profiles
include("./Lenses/Sources.jl")
export Sources

# Module for lens models
include("./Lenses/Lenses.jl")
export Lenses

end
