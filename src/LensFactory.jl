"""
    LensFactory

Testing :-)

"""
module LensFactory


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
