"""
    Lenses

Testing :-)

"""
module Lenses


# Julia inbuilt functions to import


# LensFactory modules to use
import ..cpu_vs_gpu

using ..Constants
using ..Cosmology
using ..ContourFinder
using ..IntersectionFinder


# Include the lens types files
include("./lens_types.jl")
include("./PointLens.jl")


# Various lensing function to export
export get_meshgrid
export get_critical_density
export get_potential
export get_deflection
export get_jacobian
export get_time_delay
export get_magnification
export get_image

"""
    get_meshgrid(θx::RV, θy::RV, dθ::RV) --> Tuple{Matrix{<:Float64}, Matrix{<:Float64}}

Generate a meshgrid of coordinates on which various quantities can be evaluated. At present,
this function only generates square pixels. In future if the need arises, it can be extended 
to generate rectangular pixels as well.
"""
function get_meshgrid(θx::RV, θy::RV, dθ::RV)::Tuple{Matrix{<:Float64}, Matrix{<:Float64}}
   # Making sure that grid and pixel size are positive
   if θx <= 0 || θy <= 0 || dθ <= 0
      throw(ArgumentError("All arguments must be positive."))
   end

   # Number of pixels along x- and y-directions
   nx::Int64 = round(Int64, 2.0*θx/dθ + 1.0)
   ny::Int64 = round(Int64, 2.0*θy/dθ + 1.0)

   # Initialize an empty nx x ny grid
   grid_x = Matrix{Float64}(undef, nx, ny)
   grid_y = Matrix{Float64}(undef, nx, ny)

   # Filling the empty grid with positions
   @inbounds for j = 1:ny     # Loop over y-dimension (i.e., columns)
      y_val = - θy + (j - 1.0) * dθ
      @inbounds for i = 1:nx  # Loop over y-dimension (i.e., rows)
         grid_x[i, j] = - θx + (i - 1.0) * dθ
         grid_y[i, j] = y_val
      end
   end
   return grid_x, grid_y
end


"""
    get_critical_density(Dd::RV, Dds::RV, Ds::RV; unit="kg_m2")::RV

Calculate the critical surface density ``(\\Sigma_{\\rm cr})`` given the angular diameter distances.
The result can be returned in different units,
- "kg\\_m2" ``\\Rightarrow{\\rm kg/m^2}``, 
- "msun\\_pc2" ``\\Rightarrow{\\rm M_{\\odot}/pc^2}``, 
- "msun\\_arcsec2" ``\\Rightarrow{\\rm M_{\\odot}/arcsec^2}``.
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


# Define lens_map globally (module-level or script-level)
const potential_map = Dict(
   :PointLens => (PointLens, [:D_d, :x_c, :y_c, :mass])
)
function potential_helper!(ψ::ROA, lens::AbstractLens, θ_x::ROA, θ_y::ROA)
   # Check if the lens type is in the potential_map otherwise throw an error
   entry = get(potential_map, lens._lens_, nothing)
   if entry === nothing
      throw(ArgumentError("Unknown lens type ** $(lens._lens_) **"))
   end

   # Get the function and arguments from the map
   module_name, properties = entry

   # Extract fields from lens
   args = [getfield(lens, p) for p in properties]

   # Call the potential function for specific model
   return getfield(module_name, :potential!)(ψ, θ_x, θ_y, args...)
end

"""
    get_potential(lens::AbstractLens, θ_x::ROA, θ_y::ROA) --> ROA
"""
function get_potential(lens::AbstractLens, θ_x::ROA, θ_y::ROA)::ROA
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
      return ψ
   end
end


# Define lens_map globally (module-level or script-level)
const deflection_map = Dict(
   :PointLens => (PointLens, [:D_d, :x_c, :y_c, :mass])
)
function deflection_helper!(ψx::ROA, ψy::ROA, lens::AbstractLens, θ_x::ROA, θ_y::ROA)
   # Check if the lens type is in the deflection_map otherwise throw an error
   entry = get(deflection_map, lens._lens_, nothing)
   if entry === nothing
      throw(ArgumentError("Unknown lens type ** $(lens._lens_) **"))
   end

   # Get the function and arguments from the map
   module_name, properties = entry

   # Extract fields from lens
   args = [getfield(lens, p) for p in properties]

   # Call the potential function for specific model
   return getfield(module_name, :deflection!)(ψx, ψy, θ_x, θ_y, args...)
end

"""
    get_deflection(lens::AbstractLens, θ_x::ROA, θ_y::ROA) --> Tuple{ROA, ROA}

Calculates the deflection angles (i.e., the gradient of the potential) for a given lens model. 
Returns a tuple of deflection components, i.e., ``(\\psi_x, \\psi_y)``.
"""
function get_deflection(lens::AbstractLens, θ_x::ROA, θ_y::ROA)::Tuple{ROA, ROA}
   # Check if the input coordinates are of the same type and size
   if typeof(θ_x) != typeof(θ_y) || size(θ_x) != size(θ_y)
      throw(ArgumentError("Input coordinates must be of the same type and size."))
   end

   # Initialize zero-valued potential array
   ψx::ROA = zero(θ_x)
   ψy::ROA = zero(θ_x)

   if lens._lens_ == :CompositeLens
      for component in lens._components_
         deflection_helper!(ψx, ψy, component, θ_x, θ_y)
      end
      return ψx, ψy
   else
      deflection_helper!(ψx, ψy, lens, θ_x, θ_y)
      return ψx, ψy
   end
end



# Define lens_map globally (module-level or script-level)
const jacobian_map = Dict(
   :PointLens => (PointLens, [:D_d, :x_c, :y_c, :mass])
)
function jacobian_helper!(ψxx::ROA, ψyy::ROA, ψxy::ROA, lens::AbstractLens, θ_x::ROA, θ_y::ROA)
   # Check if the lens type is in the jacobian_map otherwise throw an error
   entry = get(jacobian_map, lens._lens_, nothing)
   if entry === nothing
      throw(ArgumentError("Unknown lens type ** $(lens._lens_) **"))
   end

   # Get the function and arguments from the map
   module_name, properties = entry

   # Extract fields from lens
   args = [getfield(lens, p) for p in properties]

   # Call the jacobian function for specific model
   return getfield(module_name, :jacobian!)(ψxx, ψyy, ψxy, θ_x, θ_y, args...)
end

"""
    get_jacobian(lens::AbstractLens, θ_x::ROA, θ_y::ROA) --> Tuple{ROA, ROA, ROA}

Calculates the jacobian (i.e., deformation tensor) of the lens mapping for a given lens model.
The jacobian is a ``2\\times2`` matrix composed of the second derivatives of the potential, 
which is given as,
```math
\\mathcal{A} =
\\begin{pmatrix}
\\psi_{xx} & \\psi_{xy} \\\\
\\psi_{xy} & \\psi_{yy}
\\end{pmatrix}.
```
Since the jacobian is symmetric (for single lens plane), only three components are returned, 
i.e., ``(\\psi_{xx}, \\psi_{yy}, \\psi_{xy})``.
"""
function get_jacobian(lens::AbstractLens, θ_x::ROA, θ_y::ROA)::Tuple{ROA, ROA, ROA}
   # Check if the input coordinates are of the same type and size
   if typeof(θ_x) != typeof(θ_y) || size(θ_x) != size(θ_y)
      throw(ArgumentError("Input coordinates must be of the same type and size."))
   end

   # Initialize zero-valued potential array
   ψxx::ROA = zero(θ_x)
   ψyy::ROA = zero(θ_x)
   ψxy::ROA = zero(θ_x)

   if lens._lens_ == :CompositeLens
      for component in lens._components_
         jacobian_helper!(ψxx, ψyy, ψxy, component, θ_x, θ_y)
      end
      return ψxx, ψyy, ψxy
   else
      jacobian_helper!(ψxx, ψyy, ψxy, lens, θ_x, θ_y)
      return ψxx, ψyy, ψxy
   end
end


"""
    get_time_delay(lens::AbstractLens, θ_x::ROA, θ_y::ROA) --> ROA

Calculates the time delay for a given lens model. The corresponding expression is given as,
```math
t_d(\\pmb{\\theta}; \\pmb{\\beta}) = \\frac{1+z_l}{\\rm c} \\frac{D_d D_s}{D_{ds}}
   \\left[ \\frac{(\\pmb{\\theta} - \\pmb{\\beta})^2}{2} - \\frac{D_{ds}}{D_s} \\psi(\\pmb{\\theta}) \\right]
```
"""
function get_time_delay(lens::AbstractLens, θ_x::ROA, θ_y::ROA, zl::RV, adis::Float64, β::NTuple{2, RV})::ROA
   # Constant multiplicative factor
   constant_factor::Float64 =  ( (1.0 + zl) / CONST_C ) * lens.D_d / adis

   # Initialize zero-valued arrays to store time delay
   ϕ::ROA = zero(θ_x)

   # Get time delay components
   ϕ_potential::ROA = get_potential(lens, θ_x, θ_y)

   ax2, ax1 = axes(θ_x, 1), axes(θ_x, 2)
   @inbounds for i in ax1
      @inbounds for j in ax2
         ϕ[j, i] = constant_factor * ( 0.5 * ((θ_x[j, i] - β[1])^2 + (θ_y[j, i] - β[2])^2) - adis * ϕ_potential[j, i] )
      end
   end
   return ϕ
end


"""
    get_magnification(lens::AbstractLens, θ_x::ROA, θ_y::ROA) --> ROA

Calculates the magnification for a given lens model. The corresponding expression is given as,
```math
\\mu = \\frac{1}{\\det \\mathcal{A}} = \\frac{1}{(1 - \\kappa)^2 - \\gamma^2}
```
"""
function get_magnification(lens::AbstractLens, θ_x::ROA, θ_y::ROA, adis::Float64)::ROA
   # Get the jacobian components
   ψxx, ψyy, ψxy = get_jacobian(lens, θ_x, θ_y)

   # Scale the deformation tensor
   ψxx .*= adis
   ψyy .*= adis
   ψxy .*= adis

   # Magnification is the inverse of the determinant
   return 1.0 ./ (1.0 .+ ψxx .* ψyy .- ψxx .- ψyy .- ψxy.^2)
end


"""
    get_image(lens::AbstractLens, θ_x::ROA, θ_y::ROA, adis::Float64, β::NTuple{2, RV})::Vector{NTuple{2, RV}}
    get_image(lens::AbstractLens, θ_x::ROA, θ_y::ROA, adis::Float64, β::Matrix{<:RV})::Matrix{<:RV}

Calculates the image position for a given lens model by solving the lens equation,
```math
\\pmb{\\beta} = \\pmb{\\theta} - \\frac{D_{ds}}{D_s} \\nabla \\psi(\\pmb{\\theta})
```
"""
function get_image(lens::AbstractLens, θ_x::ROA, θ_y::ROA, adis::Float64, β::NTuple{2, RV})::Vector{NTuple{2, RV}}
   # Get the potential gradient
   ψx, ψy = get_deflection(lens, θ_x, θ_y)

   # Get grid for contour
   RXC = ContourFinder.get_contour(θ_x, θ_y, β[1] .- θ_x .+ adis .* ψx, 0.0)
   RYC = ContourFinder.get_contour(θ_x, θ_y, β[2] .- θ_y .+ adis .* ψy, 0.0)

   # Initialize empty Vector of tuples to store image positions
   image_position::Vector{NTuple{2, RV}} = []
   for contour_1 in RXC
      for contour_2 in RYC
         # Find the intersection points
         intersect_points = ContourFinder.get_intersection( first.(contour_1), last.(contour_1), first.(contour_2), last.(contour_2) )
         
         # Store the intersection points in the image_position vector
         for point in intersect_points
            push!(image_position, point)
         end
      end
   end
   return image_position
end

function get_image(lens::AbstractLens, θ_x::ROA, θ_y::ROA, adis::Float64, β::Matrix{<:RV})::Matrix{<:RV}
   # Get the potential gradient
   ψx, ψy = get_deflection(lens, θ_x, θ_y)

   # Create an empty image map
   image_map::Matrix{<:RV} = zero(θ_x)

   # Grid size
   ny, nx = size(θ_x)
   pixel_h::Float64 = abs(θ_x[2, 1] - θ_x[1, 1])

   beta_x::Float64 = 0
   beta_y::Float64 = 0
   pixel_x::Int64 = 0
   pixel_y::Int64 = 0

   # Loop over the image plane and assign values from source
   ax1, ax2 = axes(θ_x, 1), axes(θ_x, 2)

   @inbounds for j in ax2
      @inbounds for i in ax1
         # Get source plane position in radians
         beta_x = θ_x[i, j] - adis * ψx[i, j]
         beta_y = θ_y[i, j] - adis * ψy[i, j]

         # Get pixel position from radians
         pixel_x = round(Int64, beta_x / pixel_h + 0.5 * nx + 1.0)
         pixel_y = round(Int64, beta_y / pixel_h + 0.5 * ny + 1.0)

         # make sure pixel position is within bounds
         if (1 <= pixel_x <= nx) && (1 <= pixel_y <= ny)
            image_map[i, j] = β[pixel_x, pixel_y]
         end
      end
   end
   return image_map
end


end