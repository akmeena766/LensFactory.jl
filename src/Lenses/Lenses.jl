module Lenses


# Julia inbuilt functions to import


# Using cosmology from one level up (i.e., LensFactory.Main)
using ..Cosmology
using ..Constants


# Include the lens types files
include("./lens_types.jl")


# Various lensing function to export
export get_meshgrid


"""
     get_meshgrid(θx::RV, θy::RV, pixel_size::RV)

Generate a meshgrid for the x and y dimensions and pixel size.
"""
function get_meshgrid(θx::RV, θy::RV, pixel_size::RV)::Tuple{Matrix{<:Float64}, Matrix{<:Float64}}
   if θx <= 0 || θy <= 0 || pixel_size <= 0
      throw(ArgumentError("All arguments must be positive."))
   end

   nx::Int64 = floor(Int64, 2.0*θx/pixel_size + 1 + 0.5)
   ny::Int64 = floor(Int64, 2.0*θy/pixel_size + 1 + 0.5)

   grid_x = Matrix{Float64}(undef, ny, nx)
   grid_y = Matrix{Float64}(undef, ny, nx)

   @inbounds for i = 1:nx     # Loop over x-dimension (i.e., rows)
      @inbounds for j = 1:ny  # Loop over y-dimension (i.e., columns)
         grid_x[j, i] = - θx + float(i - 1) * pixel_size
         grid_y[j, i] = - θy + float(j - 1) * pixel_size
      end
   end
   return grid_x, grid_y
end

end