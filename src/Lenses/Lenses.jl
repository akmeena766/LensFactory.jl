"""
    Lenses
"""
module Lenses

# Julia inbuilt functions to import

# LensFactory modules to use
using ..Constants
using ..Cosmology

# Modules for contour finding on a 2D grid
include("../LensFactoryUtils/ContourFinder.jl")
using .ContourFinder

# Module to get intersection points of two contours
include("../LensFactoryUtils/IntersectionFinder.jl")
using .IntersectionFinder

include("../LensFactoryUtils/PolygonOps.jl")
using .PolygonOps

# Include the lens types files
include("./lens_types.jl")
include("./PointLens.jl")
include("./SISLens.jl")
include("./PlummerLens.jl")
include("./NSISPLens.jl")
include("./NSISMDLens.jl")

# Various lensing function to export
export get_meshgrid
export get_critical_density
export get_potential
export get_deflection
export get_jacobian
export get_time_delay
export get_magnification_image
export get_magnification_source
export get_image
export get_critical_curve
export get_caustic
export get_critical_area
export get_einstein_angle

"""
    get_meshgrid(θx::RV, θy::RV, dθ::RV)

Generate a meshgrid of coordinates on which various quantities can be evaluated. At present,
this function only generates square pixels. In future if the need arises, it can be extended 
to generate rectangular pixels as well.
"""
function get_meshgrid(θx::RV, θy::RV, dθ::RV)
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
    get_critical_density(Dd::RV, Dds::RV, Ds::RV; unit::Symbol=:kg_m2)

Calculate the critical surface density,
```math
Σ_{\\rm cr} = \\frac{c^2}{4 π {\\rm G}} \\frac{D_s}{D_d D_{ds}},
```
given the angular diameter distances. The result can be returned in different units,
- `:kg_m2` ``\\Rightarrow{\\rm kg/m^2}``, 
- `:msun_pc2` ``\\Rightarrow{\\rm M_⊙/pc^2}``, 
- `:msun_arcsec2` ``\\Rightarrow{\\rm M_⊙/arcsec^2}``.
"""   
function get_critical_density(; D_d::RV=NaN, adis::RV=NaN, unit::Symbol=:kg_m2)
   # Calculate Σ_cr in kg/m^2
   Σ_cr::Float64 = ( CONST_C^2 / 4.0 / π / CONST_G ) * ( 1.0 / D_d / adis )
   
   # Convert to the requested unit
   if unit == :kg_m2
      return Σ_cr
   elseif unit == :msun_pc2
      return Σ_cr * ( DIST_PC^2 / MASS_SUN )
   elseif unit == :msun_arcsec2
      return Σ_cr * ( D_d^2 * ANGLE_ARCSEC^2 / MASS_SUN )
   else
      error( "Invalid unit: $unit. Must be 'kg_m2' or 'msun_pc2' or 'msun_arcsec2'." )
   end
end


# Define lens_map globally (module-level or script-level)
const lens_map = Dict(
   :PointLens   => (PointLens,     [:D_d, :x_c, :y_c, :mass]),
   :PlummerLens => (PlummerLens,   [:D_d, :x_c, :y_c, :mass, :x_s]),
   :SISLens     => (SISLens,       [:x_c, :y_c, :v_d]),
   :NSISPLens   => (NSISPLens,     [:x_c, :y_c, :v_d, :θ_s]),
   :NSISMDLens   => (NSISMDLens,     [:x_c, :y_c, :v_d, :θ_s])
)
function potential_helper!(ψ::T, lens::AbstractLens, θx::T, θy::T) where T <: ROA
   # Check if the lens type is in the lens_map otherwise throw an error
   entry = get(lens_map, lens._lens_, nothing)
   if entry === nothing
      throw(ArgumentError("Unknown lens type ** $(lens._lens_) **"))
   end

   # Get the function and arguments from the map
   module_name, properties = entry

   # Extract fields from lens
   args = [getfield(lens, p) for p in properties]

   # Call the potential function for specific model
   return getfield(module_name, :potential!)(ψ, θx, θy, args...)
end

"""
    get_potential(lens::AbstractLens, θx::ROA, θy::ROA) --> ROA
"""
function get_potential(lens::AbstractLens, θx::T, θy::T) where T <: Union{RV, ROA}
   # If RV is passed, covert to vector
   θx = isa(θx, RV) ? [θx] : θx
   θy = isa(θy, RV) ? [θy] : θy

   # Check if the input coordinates are of the same type and size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same type and size."))
   end

   # Initialize zero-valued potential array
   ψ::ROA = zero(θx)

   if lens._lens_ == :CompositeLens
      for component in lens._components_
         potential_helper!(ψ, component, θx, θy)
      end
      return ψ
   else
      potential_helper!(ψ, lens, θx, θy)
      return ψ
   end
end

function deflection_helper!(ψx::T, ψy::T, lens::AbstractLens, θx::T, θy::T) where T <: ROA
   # Check if the lens type is in the lens_map otherwise throw an error
   entry = get(lens_map, lens._lens_, nothing)
   if entry === nothing
      throw(ArgumentError("Unknown lens type ** $(lens._lens_) **"))
   end

   # Get the function and arguments from the map
   module_name, properties = entry

   # Extract fields from lens
   args = [getfield(lens, p) for p in properties]

   # Call the potential function for specific model
   return getfield(module_name, :deflection!)(ψx, ψy, θx, θy, args...)
end

"""
    get_deflection(lens::AbstractLens, θx::ROA, θy::ROA) --> Tuple{ROA, ROA}

Calculates the deflection angles (i.e., the gradient of the potential) for a given lens model. 
Returns a tuple of deflection components, i.e., ``(ψ_x, ψ_y)``.
"""
function get_deflection(lens::AbstractLens, θx::T, θy::T) where T <: Union{RV, ROA}
   # If RV is passed, covert to vector
   θx = isa(θx, RV) ? [θx] : θx
   θy = isa(θy, RV) ? [θy] : θy

   # Check if the input coordinates are of the same type and size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same type and size."))
   end

   # Initialize zero-valued potential array
   ψx::ROA = zero(θx)
   ψy::ROA = zero(θx)

   if lens._lens_ == :CompositeLens
      for component in lens._components_
         deflection_helper!(ψx, ψy, component, θx, θy)
      end
      return ψx, ψy
   else
      deflection_helper!(ψx, ψy, lens, θx, θy)
      return ψx, ψy
   end
end


function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::AbstractLens, θx::T, θy::T) where T <: ROA
   # Check if the lens type is in the lens_map otherwise throw an error
   entry = get(lens_map, lens._lens_, nothing)
   if entry === nothing
      throw(ArgumentError("Unknown lens type ** $(lens._lens_) **"))
   end

   # Get the function and arguments from the map
   module_name, properties = entry

   # Extract fields from lens
   args = [getfield(lens, p) for p in properties]

   # Call the jacobian function for specific model
   return getfield(module_name, :jacobian!)(ψxx, ψyy, ψxy, θx, θy, args...)
end

"""
    get_jacobian(lens::AbstractLens, θx::ROA, θy::ROA) --> Tuple{ROA, ROA, ROA}

Calculates the jacobian (i.e., deformation tensor) of the lens mapping for a given lens model.
The jacobian is a ``2\\times2`` matrix composed of the second derivatives of the potential, 
which is given as,
```math
\\mathcal{A} =
\\begin{pmatrix}
ψ_{xx} & ψ_{xy} \\\\
ψ_{xy} & ψ_{yy}
\\end{pmatrix}.
```
Since the jacobian is symmetric (for single lens plane), only three components are returned,
i.e., ``(ψ_{xx}, ψ_{yy}, ψ_{xy})``.
"""
function get_jacobian(lens::AbstractLens, θx::T, θy::T) where T <: Union{RV, ROA}
   # If RV is passed, covert to vector
   θx = isa(θx, RV) ? [θx] : θx
   θy = isa(θy, RV) ? [θy] : θy

   # Check if the input coordinates are of the same type and size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same type and size."))
   end
   
   # Initialize zero-valued potential array
   ψxx::ROA = zero(θx)
   ψyy::ROA = zero(θx)
   ψxy::ROA = zero(θx)

   if lens._lens_ == :CompositeLens
      for component in lens._components_
         jacobian_helper!(ψxx, ψyy, ψxy, component, θx, θy)
      end
      return ψxx, ψyy, ψxy
   else
      jacobian_helper!(ψxx, ψyy, ψxy, lens, θx, θy)
      return ψxx, ψyy, ψxy
   end
end


"""
    get_time_delay(lens::AbstractLens, θ_x::ROA, θ_y::ROA) --> ROA

Calculates the time delay for a given lens model. The corresponding expression is given as,
```math
t_d(\\pmb{θ}; \\pmb{β}) = \\frac{1+z_l}{\\rm c} \\frac{D_d D_s}{D_{ds}}
   \\left[ \\frac{(\\pmb{θ} - \\pmb{β})^2}{2} - \\frac{D_{ds}}{D_s} \\psi(\\pmb{θ}) \\right]
```
"""
function get_time_delay(lens::AbstractLens, θx::T, θy::T, zl::RV, adis::Float64, β::NTuple{2, RV}) where T <: ROA
   # Constant multiplicative factor
   constant_factor::Float64 =  ( (1.0 + zl) / CONST_C ) * lens.D_d / adis

   # Initialize zero-valued arrays to store time delay
   ϕ::ROA = zero(θ_x)

   # Get time delay components
   ϕ_potential::ROA = get_potential(lens, θ_x, θ_y)

   ax1, ax2 = axes(θ_x, 1), axes(θ_x, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ϕ[i, j] = constant_factor * (0.5 * ((θx[i, j] - β[1])^2 + (θy[i, j] - β[2])^2) - adis * ϕ_potential[i, j]) 
      end
   end
   return ϕ
end


"""
    get_magnification(lens::AbstractLens, θx::ROA, θy::ROA) --> ROA

Calculates the magnification for a given lens model. The corresponding expression is given as,
```math
\\mu = \\frac{1}{\\det \\mathcal{A}} = \\frac{1}{(1 - κ)^2 - γ^2}
```
"""
function get_magnification_image(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: ROA
   # Get the jacobian components
   ψxx, ψyy, ψxy = get_jacobian(lens, θx, θy)

   # Scale the deformation tensor
   ψxx .*= adis
   ψyy .*= adis
   ψxy .*= adis

   # Magnification is the inverse of the determinant of jacobian
   return 1.0 ./ (1.0 .+ ψxx .* ψyy .- ψxx .- ψyy .- ψxy.^2)
end


function get_magnification_source(lens::AbstractLens, θx::T, θy::T, adis::Float64; rays_per_pixel::Int64=1) where T <: Matrix{<:RV}
   # Deflection field
   αx, αy = get_deflection(lens, θx, θy)

   # Pixel size
   nx, ny = size(θx)
   pixel_h::Float64 = abs(θx[2, 1] - θx[1, 1])

   # Area per pixel in the image plane
   area_per_ray::Float64 = 1.0 / rays_per_pixel

   # Initialize an empty source plane magnification map
   μ_source::Matrix{<:RV} = zero(θx)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for _ in ax2
      @inbounds  for _ in ax1
         # Draw random pixel numbers equal to rays_per_pixel
         rand_x = 1.0 .+ rand(Float64, rays_per_pixel) .* (nx - 1.0)
         rand_y = 1.0 .+ rand(Float64, rays_per_pixel) .* (ny - 1.0)

         for k in 1:rays_per_pixel
            # Get the source plane position
            βx = interpolation(rand_x[k], rand_y[k], θx) - adis * interpolation(rand_x[k], rand_y[k], αx)
            βy = interpolation(rand_x[k], rand_y[k], θy) - adis * interpolation(rand_x[k], rand_y[k], αy)

            # Get the corresponding pixel values
            βx_p = round(Int64, βx/pixel_h + nx/2.0)
            βy_p = round(Int64, βy/pixel_h + ny/2.0)

            # make sure pixel position is within bounds
            if (1 <= βx_p <= nx) && (1 <= βy_p <= ny)
               μ_source[βx_p, βy_p] += area_per_ray
            end
         end
      end
   end
   return μ_source
end


"""
    get_image(lens::AbstractLens, θx::ROA, θy::ROA, adis::Float64, β::NTuple{2, RV})::Vector{NTuple{2, RV}}
    get_image(lens::AbstractLens, θx::ROA, θy::ROA, adis::Float64, β::Matrix{<:RV})::Matrix{<:RV}

Calculates the image position for a given lens model by solving the lens equation,
```math
\\pmb{β} = \\pmb{θ} - \\frac{D_{ds}}{D_s} ∇\\psi(\\pmb{θ})
```
"""
function get_image(lens::AbstractLens, θx::T, θy::T, adis::Float64, β::NTuple{2, RV}) where T <: Matrix{<:RV}
   # Get the potential gradient
   ψx, ψy = get_deflection(lens, θx, θy)

   # Get grid for contour
   RXC = ContourFinder.get_contour(θx, θy, β[1] .- θx .+ adis .* ψx, 0.0)
   RYC = ContourFinder.get_contour(θx, θy, β[2] .- θy .+ adis .* ψy, 0.0)

   # Initialize empty Vector of tuples to store image positions
   image_position::Vector{NTuple{2, RV}} = []
   for contour_1 in RXC
      for contour_2 in RYC
         # Find the intersection points
         intersect_points = IntersectionFinder.get_intersection( first.(contour_1), last.(contour_1), first.(contour_2), last.(contour_2) )
         
         # Store the intersection points in the image_position vector
         for point in intersect_points
            push!(image_position, point)
         end
      end
   end
   return image_position
end

function get_image(lens::AbstractLens, θx::T, θy::T, adis::Float64, β::T) where T <: Matrix{<:RV}
   # Get the potential gradient
   ψx, ψy = get_deflection(lens, θx, θy)

   # Create an empty image map
   image_map::Matrix{<:RV} = zero(θx)

   # Grid size
   nx, ny = size(θx)
   pixel_h::Float64 = abs(θx[2, 1] - θx[1, 1])

   beta_x::Float64 = 0
   beta_y::Float64 = 0
   pixel_x::Int64 = 0
   pixel_y::Int64 = 0

   # Loop over the image plane and assign values from source
   ax1, ax2 = axes(θx, 1), axes(θx, 2)

   @inbounds for j in ax2
      @inbounds for i in ax1
         # Get source plane position in radians
         βx = θx[i, j] - adis * ψx[i, j]
         βy = θy[i, j] - adis * ψy[i, j]

         # Get pixel position from radians
         pixel_x = round(Int64, βx / pixel_h + 0.5 * nx + 1.0)
         pixel_y = round(Int64, βy / pixel_h + 0.5 * ny + 1.0)

         # make sure pixel position is within bounds
         if (1 <= pixel_x <= nx) && (1 <= pixel_y <= ny)
            image_map[i, j] = β[pixel_x, pixel_y]
         end
      end
   end
   return image_map
end


function get_critical_curve(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:RV}
   # Get the jacobian components
   ψxx, ψyy, ψxy = get_jacobian(lens, θx, θy)

   # Scale the deformation tensor
   ψxx .*= adis
   ψyy .*= adis
   ψxy .*= adis

   # Convergence and shear components
   κ::ROA = 0.5 * (ψxx + ψyy)
   γ1::ROA = 0.5 * (ψxx - ψyy)
   γ2::ROA = ψxy

   # Get the zero eigenvalue contours
   critical_tan = ContourFinder.get_contour(θx, θy, 1.0 .- κ .- sqrt.(γ1.^2 .+ γ2.^2), 0)
   critical_rad = ContourFinder.get_contour(θx, θy, 1.0 .- κ .+ sqrt.(γ1.^2 .+ γ2.^2), 0)

   return critical_tan, critical_rad   
end


function get_caustic(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:RV}
   # Generate critical curves
   critical_tan, critical_rad = get_critical_curve(lens, θx, θy, adis)

   # Get tangential caustics
   caustics_tan::Vector{Vector{Vector{Float64}}} = []
   for curve in critical_tan
      ψ_x, ψ_y = get_deflection(lens, first.(curve), last.(curve))
      src_x = first.(curve) .- adis .* ψ_x
      src_y =  last.(curve) .- adis .* ψ_y
      push!(caustics_tan, [[x, y] for (x, y) in zip(src_x, src_y)])
   end
 
   # Get radial caustics
   caustics_rad::Vector{Vector{Vector{Float64}}} = []
   for curve in critical_rad
      ψ_x, ψ_y = get_deflection(lens, first.(curve), last.(curve))
      src_x = first.(curve) .- adis .* ψ_x
      src_y =  last.(curve) .- adis .* ψ_y
      push!(caustics_rad, [[x, y] for (x, y) in zip(src_x, src_y)])
   end
   return caustics_tan, caustics_rad
end


function get_critical_area(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:RV}
   area::Float64 = 0.0

   # Get tangential critical curves
   critical_tan, _ = get_critical_curve(lens, θx, θy, adis)

   # Run a loop over all tangential critical curves
   for curve in critical_tan
      area += shoelace(curve)
   end
   return area
end


function get_einstein_angle(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:RV}
   return √(get_critical_area(lens, θx, θy, adis) / π)
end


function interpolation(x::Float64, y::Float64, df::Matrix{<:RV})::Float64
   px::Int64 = floor(Int64, x)
   py::Int64 = floor(Int64, y)

   df_00::Float64 = df[px + 0, py + 0]
   df_01::Float64 = df[px + 0, py + 1]
   df_10::Float64 = df[px + 1, py + 0]
   df_11::Float64 = df[px + 1, py + 1]

   df_interpolated::Float64 = df_00 * (px + 1 - x) * (py + 1 - y) + 
                              df_01 * (px + 1 - x) * (y - py - 0) +
                              df_10 * (x - px - 0) * (py + 1 - y) + 
                              df_11 * (x - px - 0) * (y - py - 0)
   return df_interpolated
end

end