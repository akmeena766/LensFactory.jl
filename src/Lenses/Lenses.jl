"""
    Lenses

Testing :-)

"""
module Lenses


# Julia inbuilt functions to import


# Using cosmology from one level up (i.e., LensFactory.Main)
using ..Constants
using ....Cosmology


# Include the lens types files
include("./lens_types.jl")
include("./PointLens.jl")

# Various lensing function to export
export get_meshgrid


"""
    get_meshgrid(θx::RV, θy::RV, dθ::RV) --> Tuple{Matrix{<:Float64}, Matrix{<:Float64}}

    Generate a meshgrid of coordinates on which various quantities can be evaluated.
    (-θx, -θy) +--- dθ --- dθ ---+ (+θx, -θy)
               |        |        |
               dθ       |        dθ
               |        |        |
               +--- dθ --- dθ ---+
               |        |        |
               dθ       |        dθ
               |        |        |
    (-θx, +θy) +--- dθ --- dθ ---+ (+θx, +θy)
"""
function get_meshgrid(θx::RV, θy::RV, dθ::RV)#::Tuple{Matrix{<:Float64}, Matrix{<:Float64}}
   # Making sure that grid and pixel size are positive
   if θx <= 0 || θy <= 0 || dθ <= 0
      throw(ArgumentError("All arguments must be positive."))
   end

   # Number of pixels along x- and y-directions
   nx::Int64 = floor(Int64, 2.0*θx/dθ + 1.5)
   ny::Int64 = floor(Int64, 2.0*θy/dθ + 1.5)

   # Initialize an empty nx x ny grid
   grid_x = Matrix{Float64}(undef, ny, nx)
   grid_y = Matrix{Float64}(undef, ny, nx)

   # Filling the empty grid with positions
   @inbounds for i = 1:nx     # Loop over x-dimension (i.e., rows)
      x_val = - θx + (i - 1.0) * dθ
      @inbounds for j = 1:ny  # Loop over y-dimension (i.e., columns)
         grid_x[j, i] = x_val
         grid_y[j, i] = - θy + (j - 1.0) * dθ
      end
   end
   return grid_x, grid_y
end

end