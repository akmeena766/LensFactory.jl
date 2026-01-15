module Lenses

# Julia inbuilt functions to import
using QuadGK
using SpecialFunctions

# LensFactory modules to use
using ..Constants
using ..Cosmology
using ..LFUtils

# Include the lens types files
include("./lens_types.jl")
include("./PointLens.jl")
include("./SISLens.jl")
include("./PlummerLens.jl")
include("./NSISPLens.jl")
include("./NSISMDLens.jl")
include("./GaussianLens.jl")
include("./SersicLens.jl")
include("./ExternalEffects.jl")
include("./PIEPLens.jl")
include("./SIELens.jl")
include("./PJELens.jl")
include("./HernquistLens.jl")
include("./NFWLens.jl")
include("./tNFWLens.jl")
include("./gNFWLens.jl")
include("./EinastoLens.jl")
include("./aHernquistLens.jl")
include("./aNFWLens.jl")
include("./eHernquistMDLens.jl")
include("./eNFWMDLens.jl")
include("./MultiPlummerLens.jl")
include("./MultiGaussianLens.jl")
include("./MultiPJELens.jl")

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

# Plotting functions (see ../../ext folder for functions)
export plot_sky
export plot_image_plane
export plot_surface_density
export plot_magnification_map
export plot_magnification_profile

function plot_sky end
function plot_image_plane end
function plot_surface_density end
function plot_magnification_map end
function plot_magnification_profile end


"""
    get_meshgrid(θx::RV, θy::RV, dθ::RV) --> Tuple{Matrix{Float64}, Matrix{Float64}}

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
   nx = floor(Int64, 2.0*θx/dθ + 1)
   ny = floor(Int64, 2.0*θy/dθ + 1)

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
    get_critical_density(D_d::Float64, adis::Float64; unit::Symbol=:kg_m2) --> Float64
"""
function get_critical_density(D_d::Float64, adis::Float64; unit::Symbol=:kg_m2)
   # Calculate Σ_cr in kg/m^2 based type of input parameters
   Σ_cr = (CONST_C^2 / (4.0 * π * CONST_G)) * (1.0 / (D_d * adis))
   
   # Convert to the requested unit
   if unit == :kg_m2
      return Σ_cr
   elseif unit == :msun_pc2
      return Σ_cr * ( DIST_PC^2 / MASS_SUN )
   elseif unit == :msun_arcsec2
      return Σ_cr * ( D_d^2 * ANGLE_ARCSEC^2 / MASS_SUN )
   else
      throw(ArgumentError("Invalid unit: $unit. Must be 'kg_m2' or 'msun_pc2' or 'msun_arcsec2'."))
   end
end

"""
    get_critical_density(D_d::Float64, D_ds::Float64, D_s::Float64s; unit::Symbol=:kg_m2) --> Float64

Calculate the critical surface density,
```math
Σ_{\\rm cr} = \\frac{c^2}{4 π {\\rm G}} \\frac{D_s}{D_d D_{ds}},
```
given the angular diameter distances. The result can be returned in different units,
- `:kg_m2` ``\\Rightarrow{\\rm kg/m^2}``, 
- `:msun_pc2` ``\\Rightarrow{\\rm M_⊙/pc^2}``, 
- `:msun_arcsec2` ``\\Rightarrow{\\rm M_⊙/arcsec^2}``.

"""
function get_critical_density(D_d::Float64, D_ds::Float64, D_s::Float64; unit::Symbol=:kg_m2)
   # Calculate Σ_cr in kg/m^2 based type of input parameters
   Σ_cr = (CONST_C^2 / (4.0 * π * CONST_G)) * (D_s / (D_d * D_ds))
   
   # Convert to the requested unit
   if unit == :kg_m2
      return Σ_cr
   elseif unit == :msun_pc2
      return Σ_cr * ( DIST_PC^2 / MASS_SUN )
   elseif unit == :msun_arcsec2
      return Σ_cr * ( D_d^2 * ANGLE_ARCSEC^2 / MASS_SUN )
   else
      throw(ArgumentError("Invalid unit: $unit. Must be 'kg_m2' or 'msun_pc2' or 'msun_arcsec2'."))
   end
end


"""
    get_potential(lens::AbstractLens, θx::T, θy::T) where T <: RV --> RV
"""
function get_potential(lens::AbstractLens, θx::T, θy::T) where T <: RV
   # Promote the input coordinates from Int64 to Float64
   if typeof(θx) === Int64 || typeof(θy) === Int64
      θx = Float64(θx)
      θy = Float64(θy)
   end

   # Initialize zero-valued potential scalar
   ψ = zero(θx)
   
   if lens._lens_ == :CompositeLens
      for component in lens._components_
         ψ = potential_helper!(ψ, component, θx, θy)
      end
      return ψ
   else
      ψ = potential_helper!(ψ, lens, θx, θy)
      return ψ
   end
end

"""
    get_potential(lens::AbstractLens, θx::T, θy::T) where T <: ROA --> ROA
"""
function get_potential(lens::AbstractLens, θx::T, θy::T) where T <: Union{ROA, Vector{Int64}}
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
   end

   # Promote both only if either is Int64
   if eltype(θx) === Int64 || eltype(θy) === Int64
      θx = Float64.(θx)
      θy = Float64.(θy)
   end

   # Initialize zero-valued potential array
   ψ = zero(θx)

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


"""
    get_deflection(lens::AbstractLens, θx::T, θy::T) where T <: RV --> Tuple{RV, RV}
"""
function get_deflection(lens::AbstractLens, θx::T, θy::T) where T <: RV
   # Promote the input coordinates from Int64 to Float64
   if typeof(θx) === Int64 || typeof(θy) === Int64
      θx = Float64(θx)
      θy = Float64(θy)
   end

   # Initialize zero-valued deflection scalars
   ψx = zero(θx)
   ψy = zero(θy)

   # Calculate deflection based on lens type
   if lens._lens_ == :CompositeLens
      for component in lens._components_
         ψx, ψy = deflection_helper!(ψx, ψy, component, θx, θy)
      end
      return ψx, ψy
   else
      ψx, ψy = deflection_helper!(ψx, ψy, lens, θx, θy)
      return ψx, ψy
   end
end

"""
    get_deflection(lens::AbstractLens, θx::T, θy::T) where T <: ROA --> Tuple{ROA, ROA}
Calculates the deflection angles (i.e., the gradient of the potential) for a given lens model. 
Returns a tuple of deflection components, i.e., ``(ψ_x, ψ_y)``.
"""
function get_deflection(lens::AbstractLens, θx::T, θy::T) where T <: ROA
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
   end
   
   # Promote both only if either is Int64
   if eltype(θx) === Int64 || eltype(θy) === Int64
      θx = Float64.(θx)
      θy = Float64.(θy)
   end

   # Initialize zero-valued potential array
   ψx = zero(θx)
   ψy = zero(θx)

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


"""
    get_jacobian(lens::AbstractLens, θx::T, θy::T) where T <: RV --> Tuple{RV, RV, RV}
"""
function get_jacobian(lens::AbstractLens, θx::T, θy::T) where T <: RV
   # Promote the input coordinates from Int64 to Float64
   if typeof(θx) === Int64 || typeof(θy) === Int64
      θx = Float64(θx)
      θy = Float64(θy)
   end

   # Initialize zero-valued deflection scalars
   ψxx = zero(θx)
   ψyy = zero(θy)
   ψxy = zero(θx)

   # Promote the input coordinates from Int64 to Float64
   ψxx, ψyy, ψxy, θx, θy = promote(ψxx, ψyy, ψxy, θx, θy)

   if lens._lens_ == :CompositeLens
      for component in lens._components_
         ψxx, ψyy, ψxy = jacobian_helper!(ψxx, ψyy, ψxy, component, θx, θy)
      end
      return ψxx, ψyy, ψxy
   else
      ψxx, ψyy, ψxy = jacobian_helper!(ψxx, ψyy, ψxy, lens, θx, θy)
      return ψxx, ψyy, ψxy
   end
end

"""
    get_jacobian(lens::AbstractLens, θx::T, θy::T) where T <: ROA --> Tuple{ROA, ROA, ROA}
Calculates the jacobian (i.e., deformation tensor) of the lens mapping for a given lens model.
Since the jacobian is symmetric (for single lens plane), only three components are returned,
i.e., ``(ψ_{xx}, ψ_{yy}, ψ_{xy})``.
"""
function get_jacobian(lens::AbstractLens, θx::T, θy::T) where T <: ROA
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
   end
   
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
   end
   
   # Promote both only if either is Int64
   if eltype(θx) === Int64 || eltype(θy) === Int64
      θx = Float64.(θx)
      θy = Float64.(θy)
   end

   # Initialize zero-valued potential array
   ψxx = zero(θx)
   ψyy = zero(θy)
   ψxy = zero(θx)

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
    get_time_delay(lens::AbstractLens, θx::T, θy::T, adis::Float64, z_d::RV, D_d::RV, β::NTuple{2, RV}) where T <: RV --> RV
"""
function get_time_delay(lens::AbstractLens, θx::T, θy::T, adis::Float64, z_d::RV, D_d::RV, β::NTuple{2, RV}) where T <: RV
   # Constant multiplicative factor
   constant_factor =  (1.0 + z_d) / CONST_C * (D_d / adis) * ANGLE_ARCSEC^2

   # Get potential at each point
   ϕ_potential = get_potential(lens, θx, θy)

   # Get time delay
   ϕ = constant_factor * (0.5 * ((θx - β[1])^2 + (θy - β[2])^2) - adis * ϕ_potential) 
   return ϕ
end

"""
    get_time_delay(lens::AbstractLens, θx::T, θy::T, adis::Float64, z_d::RV, D_d::RV, β::NTuple{2, RV}) where T <: ROA --> ROA
Calculates the time delay for a given lens model. The corresponding expression is given as,
```math
t_d(\\pmb{θ}; \\pmb{β}) = \\frac{1+z_l}{\\rm c} \\frac{D_d D_s}{D_{ds}} \\theta_0^2
   \\left[ \\frac{(\\pmb{θ} - \\pmb{β})^2}{2} - \\frac{D_{ds}}{D_s} \\psi(\\pmb{θ}) \\right],
```
where ``\\theta_0`` is normalizing angular unit.
"""
function get_time_delay(lens::AbstractLens, θx::T, θy::T, adis::Float64, z_d::RV, D_d::RV, β::NTuple{2, RV}) where T <: ROA
   # Constant multiplicative factor
   constant_factor =  (1.0 + z_d) / CONST_C * (D_d / adis) * ANGLE_ARCSEC^2

   # Get potential at each point
   ϕ_potential = get_potential(lens, θx, θy)

   return @. constant_factor * (0.5 * ((θx - β[1])^2 + (θy - β[2])^2) - adis * ϕ_potential)
end


"""
    get_magnification_image(lens::AbstractLens, θx::T, θy::T) where T <: RV
"""
function get_magnification_image(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: RV
   # Get the jacobian components
   ψxx, ψyy, ψxy = get_jacobian(lens, θx, θy)

   # Scale the deformation tensor
   ψxx = adis * ψxx
   ψyy = adis * ψyy
   ψxy = adis * ψxy

   # μ = 1 / det(A)
   return 1.0 / (1.0 + ψxx * ψyy - ψxx - ψyy - ψxy^2)
end

"""
    get_magnification_image(lens::AbstractLens, θx::T, θy::T) where T <: ROA 
"""
function get_magnification_image(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: ROA
   # Get the jacobian components
   ψxx, ψyy, ψxy = get_jacobian(lens, θx, θy)

   # Scale the deformation tensor
   @. ψxx = adis * ψxx
   @. ψyy = adis * ψyy
   @. ψxy = adis * ψxy

   # μ = 1 / det(A)
   return @. 1.0 / (1.0 + ψxx * ψyy - ψxx - ψyy - ψxy^2)
end


function get_magnification_source(lens::AbstractLens, θx::T, θy::T, adis::Float64; rays_per_pixel::Int64=1) where T <: Matrix{<:RV}
   # Deflection field
   ψx, ψy = get_deflection(lens, θx, θy)

   # Pixel size
   nx, ny = size(θx)
   pixel_h = abs(θx[2, 1] - θx[1, 1])

   # Area per pixel in the image plane
   area_per_ray = 1.0 / rays_per_pixel

   # Initialize an empty source plane magnification map
   μ_source = zero(θx)

   # Allocate arrays ONCE outside all loops
   rand_x = Vector{Float64}(undef, rays_per_pixel)
   rand_y = Vector{Float64}(undef, rays_per_pixel)

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for _ in ax2
      @inbounds  for _ in ax1
         # Generate random values and scale with explicit SIMD vectorization
         @simd for k in 1:rays_per_pixel
            rand_x[k] = 1.0 + rand(Float64) * (nx - 1.0)
            rand_y[k] = 1.0 + rand(Float64) * (ny - 1.0)
         end

         for k in 1:rays_per_pixel
            # Get the source plane position
            βx = PolygonOps.interpolation(rand_x[k], rand_y[k], θx) - adis * PolygonOps.interpolation(rand_x[k], rand_y[k], ψx)
            βy = PolygonOps.interpolation(rand_x[k], rand_y[k], θy) - adis * PolygonOps.interpolation(rand_x[k], rand_y[k], ψy)

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
    get_image(lens::AbstractLens, θx::ROA, θy::ROA, adis::Float64, β::NTuple{2, RV}) --> Vector{NTuple{2, RV}}
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
         intersect_points = IntersectionFinder.get_intersection(first.(contour_1), last.(contour_1), first.(contour_2), last.(contour_2))
         
         # Store the intersection points in the image_position vector
         for point in intersect_points
            push!(image_position, point)
         end
      end
   end
   return image_position
end

"""
    get_image(lens::AbstractLens, θx::ROA, θy::ROA, adis::Float64, β::Matrix{<:RV}) --> Matrix{<:RV}

Calculates the image position for a given lens model by solving the lens equation,
```math
\\pmb{β} = \\pmb{θ} - \\frac{D_{ds}}{D_s} ∇\\psi(\\pmb{θ})
```   
In case of a point source, given by ``\\pmb{β} = (β_1,\\:β_2)``, the current implementation finds 
the intersection points of contours corresponding to ``\\pmb{β} - \\pmb{θ} + a_{\\rm dis} ∇\\psi(\\pmb{θ}) = 0``.
On the other hand, for an extended source, the function does **inverse ray shooting** to construct
the image plane.
"""

function get_image(lens::AbstractLens, θx::T, θy::T, adis::Float64, β::T) where T <: Matrix{<:RV}
   # Get the potential gradient
   ψx, ψy = get_deflection(lens, θx, θy)

   # Create an empty image map
   image_map = zero(θx)

   # Grid size
   nx, ny = size(θx)
   pixel_h = abs(θx[2, 1] - θx[1, 1])

   βx = 0.0
   βy = 0.0
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

"""
    get_critical_curve(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:RV}
"""
function get_critical_curve(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:RV}
   # Get the jacobian components
   ψxx, ψyy, ψxy = get_jacobian(lens, θx, θy)

   # Scale the deformation tensor
   @. ψxx = adis * ψxx
   @. ψyy = adis * ψyy
   @. ψxy = adis * ψxy

   # Convergence and shear components
   κ  = 0.5 .* (ψxx .+ ψyy)
   γ1 = 0.5 .* (ψxx .- ψyy)
   γ2 = ψxy

   # Get the zero eigenvalue contours
   critical_tan = ContourFinder.get_contour(θx, θy, 1.0 .- κ .- sqrt.(γ1.^2 .+ γ2.^2), 0)
   critical_rad = ContourFinder.get_contour(θx, θy, 1.0 .- κ .+ sqrt.(γ1.^2 .+ γ2.^2), 0)

   return critical_tan, critical_rad   
end


"""
    get_caustic(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:RV}
"""
function get_caustic(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:RV}
   # Generate critical curves
   critical_tan, critical_rad = get_critical_curve(lens, θx, θy, adis)

   # Get tangential caustics
   caustics_tan = Vector{Vector{Vector{Float64}}}(undef, length(critical_tan))
   for (idx, curve) in enumerate(critical_tan)
      ψ_x, ψ_y = get_deflection(lens, first.(curve), last.(curve))
      src_x = first.(curve) .- adis .* ψ_x
      src_y =  last.(curve) .- adis .* ψ_y
      caustics_tan[idx] = [[x, y] for (x, y) in zip(src_x, src_y)]
   end
 
   # Get radial caustics
   caustics_rad = Vector{Vector{Vector{Float64}}}(undef, length(critical_rad))
   for (idx, curve) in enumerate(critical_rad)
      ψ_x, ψ_y = get_deflection(lens, first.(curve), last.(curve))
      src_x = first.(curve) .- adis .* ψ_x
      src_y =  last.(curve) .- adis .* ψ_y
      caustics_rad[idx] = [[x, y] for (x, y) in zip(src_x, src_y)]
   end
   return caustics_tan, caustics_rad
end


"""
    get_critical_area(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:RV} --> Float64
"""
function get_critical_area(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:RV}
   # Get tangential critical curves
   critical_tan, _ = get_critical_curve(lens, θx, θy, adis)

   # Run a loop over all tangential critical curves
   area = 0.0
   for curve in critical_tan
      area += shoelace(curve)
   end
   return area
end


"""
    get_einstein_angle(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:RV} --> Float64
"""
function get_einstein_angle(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:RV}
   return sqrt(get_critical_area(lens, θx, θy, adis) / π)
end


#--------------------- Potential functions for specific lens models -------------------------------#
@inline function potential_helper!(ψ::T, lens::init_PointLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return PointLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass)
end

@inline function potential_helper!(ψ::T, lens::init_PlummerLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return PlummerLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function potential_helper!(ψ::T, lens::init_SISLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return SISLens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.v_d)
end

@inline function potential_helper!(ψ::T, lens::init_NSISPLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return NSISPLens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s)
end

@inline function potential_helper!(ψ::T, lens::init_NSISMDLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return NSISMDLens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s)
end

@inline function potential_helper!(ψ::T, lens::init_GaussianLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return GaussianLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function potential_helper!(ψ::T, lens::init_SersicLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return SersicLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_e, lens.n)
end

@inline function potential_helper!(ψ::T, lens::init_ExternalEffects, θx::T, θy::T) where T <: Union{RV, ROA}
   return ExternalEffects.potential!(ψ, θx, θy, lens.kappa, lens.gamma, lens.angle)
end

@inline function potential_helper!(ψ::T, lens::init_PIEPLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return PIEPLens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.eps, lens.pa)
end

@inline function potential_helper!(ψ::T, lens::init_SIELens, θx::T, θy::T) where T <: Union{RV, ROA}
   return SIELens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.eps, lens.pa)
end

@inline function potential_helper!(ψ::T, lens::init_PJELens, θx::T, θy::T) where T <: Union{RV, ROA}
   return PJELens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.x_t, lens.eps, lens.pa)
end

@inline function potential_helper!(ψ::T, lens::init_HernquistLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return HernquistLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function potential_helper!(ψ::T, lens::init_NFWLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return NFWLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s)
end

@inline function potential_helper!(ψ::T, lens::init_tNFWLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return tNFWLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.x_t)
end

@inline function potential_helper!(ψ::T, lens::init_gNFWLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return gNFWLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.n)
end

@inline function potential_helper!(ψ::T, lens::init_EinastoLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return EinastoLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.n)
end

@inline function potential_helper!(ψ::T, lens::init_aHernquistLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return aHernquistLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.eps, lens.pa)
end

@inline function potential_helper!(ψ::T, lens::init_aNFWLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return aNFWLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.eps, lens.pa)
end

@inline function potential_helper!(ψ::T, lens::init_eHernquistMDLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return eHernquistMDLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.eps, lens.pa)
end

@inline function potential_helper!(ψ::T, lens::init_eNFWMDLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return eNFWMDLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.eps, lens.pa)
end

@inline function potential_helper!(ψ::T, lens::init_MultiPlummerLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return MultiPlummerLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.n)
end

@inline function potential_helper!(ψ::T, lens::init_MultiGaussianLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return MultiGaussianLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.n)
end

@inline function potential_helper!(ψ::T, lens::init_MultiPJELens, θx::T, θy::T) where T <: Union{RV, ROA}
   return MultiPJELens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.x_t, lens.eps, lens.pa, lens.n)
end


#--------------------- Deflection functions for specific lens models -------------------------------#
@inline function deflection_helper!(ψx::T, ψy::T, lens::init_PointLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return PointLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_PlummerLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return PlummerLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_SISLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return SISLens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.v_d)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_NSISPLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return NSISPLens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_NSISMDLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return NSISMDLens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_GaussianLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return GaussianLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_SersicLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return SersicLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_e, lens.n)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_ExternalEffects, θx::T, θy::T) where T <: Union{RV, ROA}
   return ExternalEffects.deflection!(ψx, ψy, θx, θy, lens.kappa, lens.gamma, lens.angle)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_PIEPLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return PIEPLens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.eps, lens.pa)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_SIELens, θx::T, θy::T) where T <: Union{RV, ROA}
   return SIELens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.eps, lens.pa)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_PJELens, θx::T, θy::T) where T <: Union{RV, ROA}
   return PJELens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.x_t, lens.eps, lens.pa)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_HernquistLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return HernquistLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_NFWLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return NFWLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_tNFWLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return tNFWLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.x_t)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_gNFWLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return gNFWLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.n)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_EinastoLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return EinastoLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.n)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_aHernquistLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return aHernquistLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.eps, lens.pa)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_aNFWLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return aNFWLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.eps, lens.pa)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_eHernquistMDLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return eHernquistMDLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.eps, lens.pa)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_eNFWMDLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return eNFWMDLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.eps, lens.pa)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_MultiPlummerLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return MultiPlummerLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.n)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_MultiGaussianLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return MultiGaussianLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.n)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_MultiPJELens, θx::T, θy::T) where T <: Union{RV, ROA}
   return MultiPJELens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.x_t, lens.eps, lens.pa, lens.n)
end


#--------------------- Deformation tensor for various lens models ---------------------------------#
@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_PointLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return PointLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_PlummerLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return PlummerLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_SISLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return SISLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.v_d)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_NSISPLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return NSISPLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_NSISMDLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return NSISMDLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_GaussianLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return GaussianLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_SersicLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return SersicLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_e, lens.n)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_ExternalEffects, θx::T, θy::T) where T <: Union{RV, ROA}
   return ExternalEffects.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.kappa, lens.gamma, lens.angle)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_PIEPLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return PIEPLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.eps, lens.pa)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_SIELens, θx::T, θy::T) where T <: Union{RV, ROA}
   return SIELens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.eps, lens.pa)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_PJELens, θx::T, θy::T) where T <: Union{RV, ROA}
   return PJELens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.x_t, lens.eps, lens.pa)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_HernquistLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return HernquistLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_NFWLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return NFWLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_tNFWLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return tNFWLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.x_t)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_gNFWLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return gNFWLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.n)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_EinastoLens, θx::T, θy::T) where T <: Union{RV, ROA}
   EinastoLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.n)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_aHernquistLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return aHernquistLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.eps, lens.pa)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_aNFWLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return aNFWLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.eps, lens.pa)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_eHernquistMDLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return eHernquistMDLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.eps, lens.pa)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_eNFWMDLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return eNFWMDLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.rho_s, lens.x_s, lens.eps, lens.pa)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_MultiPlummerLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return MultiPlummerLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.n)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_MultiGaussianLens, θx::T, θy::T) where T <: Union{RV, ROA}
   return MultiGaussianLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.n)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_MultiPJELens, θx::T, θy::T) where T <: Union{RV, ROA}
   return MultiPJELens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.x_t, lens.eps, lens.pa, lens.n)
end


end