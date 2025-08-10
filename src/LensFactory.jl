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


# Module for cosmology
include("./Cosmology/Cosmology.jl")
export Cosmology


# Module for lens models
include("./Lenses/Lenses.jl")
export Lenses

end
