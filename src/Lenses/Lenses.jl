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
export get_critical_density
export get_potential

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


"""
    get_critical_density(Dd::RV, Dds::RV, Ds::RV; unit="kg_m2")::RV

Calculate the critical surface density ``(\\Sigma_{\\rm cr})`` given the angular diameter distances.
The result can be returned in different units: "kg\\_m2" ``{\\rm (i.e., kg/m^2)}``, 
"msun\\_pc2" ``{\\rm (i.e., M_{\\odot}/pc^2)}``, 
or "msun\\_arcsec2" ``{\\rm (i.e., M_{\\odot}/arcsec^2)}``.
"""   
function get_critical_density(Dd::RV, Dds::RV, Ds::RV; unit::String="kg_m2")::RV
   # Calculate Σ_cr in kg/m^2
   Σ_cr::Float64 = ( CONST_C^2 / 4.0 / π / CONST_G ) * ( Ds / Dd / Dds )
   
   # Convert to the requested unit
   if unit == "kg_m2"
      return Σ_cr
   elseif unit == "msun_pc2"
      return Σ_cr * ( DIST_PC^2 / MASS_SUN )
   elseif unit == "msun_arcsec2"
      return Σ_cr * ( Dol^2 * ANGLE_ARCSEC^2 / MASS_SUN )
   else
      error( "Invalid unit: $unit. Must be 'kg_m2' or 'msun_pc2' or 'msun_arcsec2'." )
   end
end


function get_potential(lens::AbstractLens, θ_x::ROA, θ_y::ROA)
   # Check if the input coordinates are of the same type and size
   if typeof(θ_x) != typeof(θ_y) || size(θ_x) != size(θ_y)
      throw(ArgumentError("Input coordinates must be of the same type and size."))
   end

   # Initialize zero-valued potential array
   ψ::ROA = zero(θ_x)

   if lens._lens_ == :CompositeLens
      for component in lens._components_
         potential_helper!(ψ, component, θ_x, θ_y)
      end
      return ψ
   else
      potential_helper!(ψ, lens, θ_x, θ_y)
   end
end


# Define lens_map globally (module-level or script-level)
const potential_map = Dict(
   :PointLens => (PointLens, [:D_d, :x_c, :y_c, :mass])
)
function potential_helper!(ψ::ROA, lens::AbstractLens, θ_x::ROA, θ_y::ROA)
   # Check if the lens type is in the potential_map otherwise throw an error
   entry = get(potential_map, lens._lens_, nothing)
   if entry === nothing
      throw(ArgumentError("Unknown lens type ** $(lens._lens_type) **"))
   end

   # Get the function and arguments from the map
   module_name, properties = entry

   # Extract fields from lens
   args = [getfield(lens, p) for p in properties]

   # Call the potential function for specific model
   return getfield(module_name, :potential!)(ψ, θ_x, θ_y, args...)
end


end