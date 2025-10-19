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

# Module for various polygon operations
include("../LensFactoryUtils/PolygonOps.jl")
using .PolygonOps

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
export plot_image_plane
export plot_surface_density
export plot_magnification_map
export plot_magnification_profile

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


# Define lens_map globally (module-level)
const lens_map = Dict(
   :PointLens         => (PointLens,         [:D_d, :x_c, :y_c, :mass]),
   :PlummerLens       => (PlummerLens,       [:D_d, :x_c, :y_c, :mass, :x_s]),
   :SISLens           => (SISLens,           [:x_c, :y_c, :v_d]),
   :NSISPLens         => (NSISPLens,         [:x_c, :y_c, :v_d, :x_s]),
   :NSISMDLens        => (NSISMDLens,        [:x_c, :y_c, :v_d, :x_s]),
   :GaussianLens      => (GaussianLens,      [:D_d, :x_c, :y_c, :mass, :x_s]),
   :SersicLens        => (SersicLens,        [:D_d, :x_c, :y_c, :mass, :x_e, :n]),
   :ExternalEffects   => (ExternalEffects,   [:kappa, :gamma1, :gamma2]),
   :PIEPLens          => (PIEPLens,          [:x_c, :y_c, :v_d, :x_s, :eps, :pa]),
   :SIELens           => (SIELens,           [:x_c, :y_c, :v_d, :x_s, :eps, :pa]),
   :PJELens           => (PJELens,           [:x_c, :y_c, :v_d, :x_s, :x_t, :eps, :pa]),
   :HernquistLens     => (HernquistLens,     [:D_d, :x_c, :y_c, :mass, :x_s]),
   :NFWLens           => (NFWLens,           [:D_d, :x_c, :y_c, :rho_s, :x_s]),
   :MultiPlummerLens  => (MultiPlummerLens,  [:D_d, :x_c, :y_c, :mass, :x_s, :n]),
   :MultiGaussianLens => (MultiGaussianLens, [:D_d, :x_c, :y_c, :mass, :x_s, :n]),
   :MultiPJELens      => (MultiPJELens,      [:x_c, :y_c, :v_d, :x_s, :x_t, :eps, :pa, :n])
)
function potential_helper!(ψ::T, lens::AbstractLens, θx::T, θy::T) where T <: Union{RV, ROA}
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
    get_potential(lens::AbstractLens, θx::T, θy::T) where T <: RV --> RV
"""
function get_potential(lens::AbstractLens, θx::T, θy::T) where T <: RV
   # Initialize zero-valued potential array
   ψ = 0.0

   # Promote the input coordinates from Int64 to Float64
   ψ, θx, θy = promote(ψ, θx, θy)
   
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
function get_potential(lens::AbstractLens, θx::T, θy::T) where T <: ROA
   # Check if the input coordinates are of the same size
   if size(θx) != size(θy)
      throw(ArgumentError("Input coordinates must be of the same size."))
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

function deflection_helper!(ψx::T, ψy::T, lens::AbstractLens, θx::T, θy::T) where T <: Union{RV, ROA}
   # Check if the lens type is in the lens_map otherwise throw an error
   entry = get(lens_map, lens._lens_, nothing)
   if entry === nothing
      throw(ArgumentError("Unknown lens type ** $(lens._lens_) **"))
   end

   # Get the function and arguments from the map
   module_name, properties = entry

   # Extract fields from lens
   args = [getfield(lens, p) for p in properties]

   # Call the potential function for specific model)
   return getfield(module_name, :deflection!)(ψx, ψy, θx, θy, args...)
end

"""
    get_deflection(lens::AbstractLens, θx::T, θy::T) where T <: RV --> Tuple{RV, RV}
"""
function get_deflection(lens::AbstractLens, θx::T, θy::T) where T <: RV
   # Initialize zero-valued potential array
   ψx = 0.0
   ψy = 0.0

   # Promote the input coordinates from Int64 to Float64
   ψx, ψy, θx, θy = promote(ψx, ψy, θx, θy)

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


function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::AbstractLens, θx::T, θy::T) where T <: Union{RV, ROA}
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
    get_jacobian(lens::AbstractLens, θx::T, θy::T) where T <: RV --> Tuple{RV, RV, RV}
"""
function get_jacobian(lens::AbstractLens, θx::T, θy::T) where T <: RV
   # Initialize zero-valued potential array
   ψxx = 0.0
   ψyy = 0.0
   ψxy = 0.0

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
   
   # Initialize zero-valued potential array
   ψxx = zero(θx)
   ψyy = zero(θx)
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
    get_time_delay(lens::AbstractLens, θx::T, θy::T, z_d::RV, D_d::RV, adis::Float64, β::NTuple{2, RV}) where T <: RV --> RV
"""
function get_time_delay(lens::AbstractLens, θx::T, θy::T, z_d::RV, D_d, adis::Float64, β::NTuple{2, RV}) where T <: RV
   # Constant multiplicative factor
   constant_factor =  (1.0 + z_d) / CONST_C * (D_d / adis) * ANGLE_ARCSEC^2

   # Get potential at each point
   ϕ_potential = get_potential(lens, θx, θy)

   # Get time delay
   ϕ = constant_factor * (0.5 * ((θx - β[1])^2 + (θy - β[2])^2) - adis * ϕ_potential) 
   return ϕ
end

"""
    get_time_delay(lens::AbstractLens, θx::T, θy::T, z_d::RV, D_d::RV, adis::Float64, β::NTuple{2, RV}) where T <: ROA --> ROA
Calculates the time delay for a given lens model. The corresponding expression is given as,
```math
t_d(\\pmb{θ}; \\pmb{β}) = \\frac{1+z_l}{\\rm c} \\frac{D_d D_s}{D_{ds}} \\theta_0^2
   \\left[ \\frac{(\\pmb{θ} - \\pmb{β})^2}{2} - \\frac{D_{ds}}{D_s} \\psi(\\pmb{θ}) \\right],
```
where ``\\theta_0`` is normalizing angular unit.
"""
function get_time_delay(lens::AbstractLens, θx::T, θy::T, z_d::RV, D_d::RV, adis::Float64, β::NTuple{2, RV}) where T <: ROA
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

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   @inbounds for _ in ax2
      @inbounds  for _ in ax1
         # Draw random pixel numbers equal to rays_per_pixel
         rand_x = 1.0 .+ rand(Float64, rays_per_pixel) .* (nx - 1.0)
         rand_y = 1.0 .+ rand(Float64, rays_per_pixel) .* (ny - 1.0)

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
         intersect_points = IntersectionFinder.get_intersection( first.(contour_1), last.(contour_1), first.(contour_2), last.(contour_2) )
         
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
   caustics_tan = Vector{Vector{Vector{Float64}}}()
   for curve in critical_tan
      ψ_x, ψ_y = get_deflection(lens, first.(curve), last.(curve))
      src_x = first.(curve) .- adis .* ψ_x
      src_y =  last.(curve) .- adis .* ψ_y
      push!(caustics_tan, [[x, y] for (x, y) in zip(src_x, src_y)])
   end
 
   # Get radial caustics
   caustics_rad = Vector{Vector{Vector{Float64}}}()
   for curve in critical_rad
      ψ_x, ψ_y = get_deflection(lens, first.(curve), last.(curve))
      src_x = first.(curve) .- adis .* ψ_x
      src_y =  last.(curve) .- adis .* ψ_y
      push!(caustics_rad, [[x, y] for (x, y) in zip(src_x, src_y)])
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

end