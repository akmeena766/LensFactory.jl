module Lenses


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using QuadGK
using StatsBase
using SpecialFunctions


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..Constants
using ..Cosmology
using ..LFUtils

# Include the lens model files
include("./lens_types.jl")
include("./PointLens.jl")
include("./SISLens.jl")
include("./PlummerLens.jl")
include("./NSISPLens.jl")
include("./NSISMDLens.jl")
include("./GaussianLens.jl")
include("./SersicLens.jl")
include("./PixelLens.jl")
include("./ExternalEffects.jl")
include("./ExternalEffects3.jl")
include("./Multipole.jl")
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
include("./MultiPixelLens.jl")
include("./MultiPJELens.jl")


# --------------------------------------------------------------------------------------------------
# Functions to export
# --------------------------------------------------------------------------------------------------
export get_meshgrid
export get_critical_density
export get_potential
export get_deflection
export get_jacobian
export get_time_delay
export get_magnification_image
export get_magnification_source
export get_image
export get_kappa_gamma
export get_critical_curve
export get_caustic
export get_critical_area
export get_einstein_angle
export get_image_multiplicity
export get_radial_profile
export get_mass_profile
export shear_cartesian2polar, shear_polar2cartesian
export ellipticity_cartesian2polar, ellipticity_polar2cartesian
export parameter_NFWLens, parameter_gNFWLens, parameter_EinastoLens


# --------------------------------------------------------------------------------------------------
# Plotting functions (see ../../ext folder for functions)
# --------------------------------------------------------------------------------------------------
export plot_sky
export set_plotKws!
export plot_image_plane
export plot_surface_density
export plot_magnification_map
export plot_magnification_profile

function plot_sky end
function set_plotKws! end
function plot_image_plane end
function plot_surface_density end
function plot_magnification_map end
function plot_magnification_profile end


"""
    get_meshgrid(θx::Real, θy::Real, dθ::Real) --> Tuple(Matrix{<:Real}, Matrix{<:Real})
Generate a meshgrid of coordinates on which various quantities can be evaluated. At present, this 
function only generates square pixels.

# Arguments
- `θx`: Half-size of the grid in x-direction (in ``\\rm \\mathbf{arcseconds}``)
- `θy`: Half-size of the grid in y-direction (in ``\\rm \\mathbf{arcseconds}``)
- `dθ`: Pixel size (in ``\\rm \\mathbf{arcseconds}``)

# Returns
- `grid_x`: x-coordinates of the grid (in ``\\rm \\mathbf{arcseconds}``)
- `grid_y`: y-coordinates of the grid (in ``\\rm \\mathbf{arcseconds}``)
"""
function get_meshgrid(θx::Real, θy::Real, dθ::Real)
   # Promote to common type
   θx, θy, dθ = promote(θx, θy, dθ)
   T = typeof(θx)

   # Making sure that grid and pixel size are positive
   if θx <= 0 || θy <= 0 || dθ <= 0
      throw(ArgumentError("All arguments must be positive."))
   end

   # Number of pixels along x- and y-directions
   nx = floor(Int, 2.0 * θx / dθ + 1)
   ny = floor(Int, 2.0 * θy / dθ + 1)

   # Initialize an empty nx x ny grid
   grid_x = Matrix{T}(undef, nx, ny)
   grid_y = Matrix{T}(undef, nx, ny)

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
    get_critical_density(D_d::Real, D_ds::Real, D_s::Real; 
                         unit::Symbol=:kg_m2) --> Real
Calculate the critical surface density,
```math
Σ_{\\rm cr} = \\frac{c^2}{4 π {\\rm G}} \\frac{D_s}{D_d D_{ds}}.
``` 

# Arguments
- `D_d` : Angular diameter distance to the lens (in ``\\rm \\mathbf{meters}``)
- `D_ds`: Angular diameter distance to the source (in ``\\rm \\mathbf{meters}``)
- `D_s` : Angular diameter distance to the lens (in ``\\rm \\mathbf{meters}``)

# Keyword Arguments
- `unit = :kg_m2`: Output units of the critical surface density.
   - `:kg_m2` ``\\Rightarrow{\\rm \\mathbf{kg/m^2}}``
   - `:msun_pc2` ``\\Rightarrow{\\rm \\mathbf{M_⊙/pc^2}}``
   - `:msun_arcsec2` ``\\Rightarrow{\\rm \\mathbf{M_⊙/arcsec^2}}``

# Returns
- `Σ_cr::Real`: Critical surface density in the requested unit
"""
function get_critical_density(D_d::Real, D_ds::Real, D_s::Real; unit::Symbol=:kg_m2)
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
      throw(ArgumentError("Invalid unit: $unit. Must be ':kg_m2' or ':msun_pc2' or ':msun_arcsec2'."))
   end
end


# --------------------------------------------------------------------------------------------------
# 
# --------------------------------------------------------------------------------------------------
lens_eltype(lens::AbstractLens) = typeof(lens).parameters[1]
lens_eltype(lens::init_CompositeLens)  = mapreduce(lens_eltype, promote_type, lens._components_)
lens_eltype(lens::init_MultiPlaneLens) = mapreduce(lens_eltype, promote_type, lens._plane_)



# --------------------------------------------------------------------------------------------------
# 
# --------------------------------------------------------------------------------------------------
"""
    get_potential(lens::AbstractLens, θx::T, θy::T) where T <: Real --> Real
"""
function get_potential(lens::AbstractLens, θx::Real, θy::Real)
   # Promote to common type
   θx, θy = promote(θx, θy)
   T = typeof(θx)

   # Initialize zero-valued potential scalar
   ψ = zero(T)
   
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
Calculate lensing potential at the given angular coordinates for the given lens model,
```math
ψ(\\pmb{θ}) = \\frac{4{\\rm G}}{\\rm c^2} \\frac{1}{D_d} \\int_{\\mathbb{R}^2} d^2 \\pmb{θ}' \\, Σ(\\pmb{θ}') \\, \\ln|\\pmb{θ} - \\pmb{θ'}|,
```
where `` Σ(\\pmb{θ})`` is in units of ``\\rm{\\mathbf{kg/arcsec^2}}``.

# Arguments
- `lens`: Lens model.
- `θx`  : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `ψ`  : Lensing potential at the given angular coordinate(s).
"""
function get_potential(lens::AbstractLens, θx::T, θy::T) where T <: ROA
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
   OutT = promote_type(eltype(θx), eltype(θy), lens_eltype(lens))
   ψ = zeros(OutT, size(θx))

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
    get_deflection(lens::AbstractLens, θx::T, θy::T) where T <: Real --> Real
"""
function get_deflection(lens::AbstractLens, θx::Real, θy::Real)
   # Promote to common type
   θx, θy = promote(θx, θy)
   T = typeof(θx)

   # Initialize zero-valued deflection scalars
   ψx = zero(T)
   ψy = zero(T)

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
Calculate (vector) deflection angle at the given angular coordinate(s) for a given lens model,
```math
\\pmb{α}(\\pmb{θ}) = \\pmb{∇} ψ(\\pmb{θ}) 
= \\frac{4{\\rm G}}{\\rm c^2} \\frac{1}{D_d} \\int_{\\mathbb{R}^2} d^2 \\pmb{θ}' \\, Σ(\\pmb{θ}') \\frac{\\pmb{θ} - \\pmb{θ}'}{|\\pmb{θ} - \\pmb{θ}'|^2}.
```

# Arguments
- `lens`: Lens model.
- `θx`  : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `αx`  : x-component of the deflection angle (in ``\\rm \\mathbf{arcseconds}``).
- `αy`  : y-component of the deflection angle (in ``\\rm \\mathbf{arcseconds}``).
"""
function get_deflection(lens::AbstractLens, θx::T, θy::T) where T <: Union{ROA, Vector{Int64}}
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
   OutT = promote_type(eltype(θx), eltype(θy), lens_eltype(lens))
   ψx = zeros(OutT, size(θx))
   ψy = zeros(OutT, size(θy))

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
    get_jacobian(lens::AbstractLens, θx::T, θy::T) where T <: Real --> Tuple{Real, Real, Real}
"""
function get_jacobian(lens::AbstractLens, θx::Real, θy::Real)
   # Promote to common type
   θx, θy = promote(θx, θy)
   T = typeof(θx)

   # Initialize zero-valued deflection scalars
   ψxx = zero(T)
   ψyy = zero(T)
   ψxy = zero(T)

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
    get_jacobian(lens::AbstractLens, θx::T, θy::T) where T <: ROA --> Tuple(ROA, ROA, ROA)
Calculate jacobian (i.e., deformation tensor) of the lens mapping for a given lens model,
```math
\\mathbb{A}(\\pmb{θ}) = 
\\begin{pmatrix}
   ψ_{xx} & ψ_{xy} \\\\
   ψ_{xy} & ψ_{yy}
\\end{pmatrix}.
```

Since the jacobian is symmetric (for single lens plane), only three components are returned,
i.e., ``(ψ_{xx}, ψ_{yy}, ψ_{xy})``.

# Arguments
- `lens`: Lens model.
- `θx`  : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `ψxx` : xx-component of the jacobian.
- `ψyy` : yy-component of the jacobian.
- `ψxy` : xy-component of the jacobian.
"""
function get_jacobian(lens::AbstractLens, θx::T, θy::T) where T <: Union{ROA, Vector{Int64}}
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
   OutT = promote_type(eltype(θx), eltype(θy), lens_eltype(lens))
   ψxx = zeros(OutT, size(θx))
   ψyy = zeros(OutT, size(θy))
   ψxy = zeros(OutT, size(θx))

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


# --------------------------------------------------------------------------------------------------
# 
# --------------------------------------------------------------------------------------------------
"""
    get_time_delay(lens::AbstractLens, θx::T, θy::T, 
                   adis::Real, 
                   z_d::Real, 
                   D_d::Real, 
                   β::NTuple{2, Real}) where T <: Real --> Real
"""
function get_time_delay(lens::AbstractLens, θx::T, θy::T, adis::Real, z_d::Real, D_d::Real, β::NTuple{2, Real}) where T <: Real
   # Constant multiplicative factor
   constant_factor =  (1.0 + z_d) / CONST_C * (D_d / adis) * ANGLE_ARCSEC^2

   # Get potential at each point
   ϕ_potential = get_potential(lens, θx, θy)

   # Get time delay
   ϕ = constant_factor * (0.5 * ((θx - β[1])^2 + (θy - β[2])^2) - adis * ϕ_potential) 
   return ϕ
end

"""
    get_time_delay(lens::AbstractLens, θx::T, θy::T, 
                   adis::Real, 
                   z_d::Real, 
                   D_d::Real, 
                   β::NTuple{2, Real}) where T <: ROA --> ROA
Calculate time delay for a given lens model and source position. The corresponding expression is given as,
```math
t_d(\\pmb{θ}; \\pmb{β}) = \\frac{1+z_d}{\\rm c} \\frac{D_d D_s}{D_{ds}} \\theta_0^2
   \\left[ \\frac{(\\pmb{θ} - \\pmb{β})^2}{2} - \\frac{D_{ds}}{D_s} \\psi(\\pmb{θ}) \\right],
```
where ``\\theta_0`` is normalizing angular unit. Since all the angular coordinates are in arcseconds,
for our case, ``\\mathbf{\\theta_0 = 1~\\rm \\mathbf{arcsecond}}``.

# Arguments
- `lens`: Lens model.
- `θx`  : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `adis`: Distance ratio (i.e., ``D_{ds}/D_s``).
- `z_d` : Lens redshift.
- `D_d` : Angular diameter distance to the lens (in ``\\rm \\mathbf{meters}``).
- `β`   : Source angular position (in ``\\rm \\mathbf{arcseconds}``).

# Returns
- `t_d`: Time delay at the given angular coordinate(s) (in ``\\rm \\mathbf{seconds}``).
"""
function get_time_delay(lens::AbstractLens, θx::T, θy::T, adis::Real, z_d::Real, D_d::Real, β::NTuple{2, Real}) where T <: ROA
   # Constant multiplicative factor
   constant_factor =  (1.0 + z_d) / CONST_C * (D_d / adis) * ANGLE_ARCSEC^2

   # Get potential at each point
   ϕ_potential = get_potential(lens, θx, θy)

   return @. constant_factor * (0.5 * ((θx - β[1])^2 + (θy - β[2])^2) - adis * ϕ_potential)
end

"""
    get_kappa_gamma(lens::AbstractLens, θx::T, θy::T, 
                    adis::Real) where T <: Real --> Real
"""
function get_kappa_gamma(lens::AbstractLens, θx::T, θy::T, adis::Real) where T <: Real
   # Get the jacobian components
   ψxx, ψyy, ψxy = get_jacobian(lens, θx, θy)

   # Scale the deformation tensor
   ψxx = adis * ψxx
   ψyy = adis * ψyy
   ψxy = adis * ψxy

   # Convergence and shear components
   κ  = 0.5 * (ψxx + ψyy)
   γ1 = 0.5 * (ψxx - ψyy)
   γ2 = ψxy

   return κ, γ1, γ2
end

"""
    get_kappa_gamma(lens::AbstractLens, θx::T, θy::T, 
                    adis::Real) where T <: ROA --> ROA
Calculate convergence (``\\kappa``) and shear components (``\\gamma_1``, ``\\gamma_2``) at the 
given angular coordinate(s) for a given lens model. The corresponding formulae are,
```math
\\begin{equation}
\\begin{split}
   κ   &= (ψ_{xx} + ψ_{yy})/2, \\\\
   γ_1 &= (ψ_{xx} - ψ_{yy})/2, \\\\
   γ_2 &= ψ_{xy}.
\\end{split}
\\end{equation}
```

# Arguments
- `lens`: Lens model.
- `θx`  : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `adis`: Distance ratio (i.e., ``D_{ds}/D_s``).

# Returns
- `κ`   : Convergence at the given angular coordinate(s).
- `γ1`  : First component of shear at the given angular coordinate(s).
- `γ2`  : Second component of shear at the given angular coordinate(s).
"""
function get_kappa_gamma(lens::AbstractLens, θx::T, θy::T, adis::Real) where T <: ROA
   # Get the jacobian components
   ψxx, ψyy, ψxy = get_jacobian(lens, θx, θy)

   # Convergence and shear components
   κ  = 0.5 .* adis .* (ψxx .+ ψyy)
   γ1 = 0.5 .* adis .* (ψxx .- ψyy)
   γ2 =        adis .* ψxy

   return κ, γ1, γ2
end


"""
    get_magnification_image(lens::AbstractLens, θx::T, θy::T, 
                            adis::Real) where T <: Real --> Real
"""
function get_magnification_image(lens::AbstractLens, θx::T, θy::T, adis::Real) where T <: Real
   # Get the jacobian components
   ψxx, ψyy, ψxy = get_jacobian(lens, θx, θy)

   # Convergence and shear components
   κ  = 0.5 * adis * (ψxx + ψyy)
   γ1 = 0.5 * adis * (ψxx - ψyy)
   γ2 =        adis * ψxy

   # μ = 1 / det(A)
   return 1.0 / ( (1.0 - κ)^2 - (γ1^2 + γ2^2))
end


"""
    get_magnification_image(lens::AbstractLens, θx::T, θy::T, 
                            adis::Real) where T <: ROA --> ROA
Calculate signed magnification at the given angular coordinate(s) for a given lens model,
```math
\\mu(\\pmb{θ}) = \\frac{1}{det\\left[ \\mathbb{I} - a_{\\rm dis} \\, \\mathbb{A} \\right]},
```
where ``\\mathbb{I}`` is the identity matrix.

# Arguments
- `lens`: Lens model.
- `θx`  : x-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-coordinate(s) (in ``\\rm \\mathbf{arcseconds}``).
- `adis`: Distance ratio (i.e., ``D_{ds}/D_s``).

# Returns
- `μ`   : Magnification at the given angular coordinate(s).
"""
function get_magnification_image(lens::AbstractLens, θx::T, θy::T, adis::Real) where T <: ROA
   # Get the jacobian components
   ψxx, ψyy, ψxy = get_jacobian(lens, θx, θy)

   # Convergence and shear components
   κ  = 0.5 .* adis .* (ψxx .+ ψyy)
   γ1 = 0.5 .* adis .* (ψxx .- ψyy)
   γ2 =        adis .* ψxy

   # μ = 1 / det(1 - A)
   return @. 1.0 / ( (1.0 - κ)^2 - (γ1^2 + γ2^2))
end


"""
    get_magnification_source(lens::AbstractLens, θx::T, θy::T, 
                             adis::Real; 
                             rays_per_pixel::Int64 = 1) where T <: Matrix{<:Real} --> Matrix{<:Real}
Calculates the magnification map in source plane using inverse ray shooting (IRS) for a given 
lens model. The number of average rays per pixel can be specified using the `rays_per_pixel` 
keyword argument. **This function is not optimized for speed and is only intended to visualize 
the magnification map.**

# Arguments
- `lens`: Lens model.
- `θx`  : x-grid (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-grid (in ``\\rm \\mathbf{arcseconds}``).
- `adis`: Distance ratio (i.e., ``D_{ds}/D_s``).

# Keyword Arguments
- `rays_per_pixel = 1`: Average number of rays per pixel.

# Returns
- `μ`   : Magnification map in source plane.
"""
function get_magnification_source(lens::AbstractLens, θx::T, θy::T, adis::Real; 
                                  rays_per_pixel::Int64 = 1) where T <: Matrix{<:Real}
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
            # Interpolated deflection valies
            ψx_interp = PolygonOps.bilinear_interpolation(rand_x[k], rand_y[k], ψx)
            ψy_interp = PolygonOps.bilinear_interpolation(rand_x[k], rand_y[k], ψy)

            # Get the source plane position
            βx = PolygonOps.bilinear_interpolation(rand_x[k], rand_y[k], θx) - adis * ψx_interp
            βy = PolygonOps.bilinear_interpolation(rand_x[k], rand_y[k], θy) - adis * ψy_interp

            # Get the corresponding pixel values
            βx_p = round(Int64, βx/pixel_h + 0.5*nx + 0.5)
            βy_p = round(Int64, βy/pixel_h + 0.5*ny + 0.5)

            # make sure pixel position is within bounds
            if (1 <= βx_p <= nx) && (1 <= βy_p <= ny)
               μ_source[βx_p, βy_p] = μ_source[βx_p, βy_p] + area_per_ray
            end
         end
      end
   end
   return μ_source
end


"""
    get_image(lens::AbstractLens, θx::T, θy::T, 
              adis::Real, 
              β::NTuple{2, Real}) where T <: AbstractMatrix{<:Real} --> Vector{NTuple{2, Real}}
Calculate image positions for a given lens model and point source position. To get the image 
positions, this implementation finds the intersection points of contours corresponding to,
```math
\\pmb{β} - \\pmb{θ} + a_{\\rm dis} \\, \\pmb{α}(\\pmb{θ}) = 0,
```
where ``\\pmb{β}`` is the source position, ``\\pmb{θ}`` is the image plane grid, ``a_{\\rm dis}`` 
is the distance ratio (i.e., ``D_{ds}/D_s``), and ``\\pmb{α}(\\pmb{θ})`` is the deflection angle. 
"""
function get_image(lens::AbstractLens, θx::T, θy::T, adis::Real, β::NTuple{2, Real}) where T <: AbstractMatrix{<:Real}
   # Get the potential gradient
   ψx, ψy = get_deflection(lens, θx, θy)

   # Get grid for contour
   RXC = ContourFinder.get_contour(θx, θy, β[1] .- θx .+ adis .* ψx, 0.0)
   RYC = ContourFinder.get_contour(θx, θy, β[2] .- θy .+ adis .* ψy, 0.0)

   # Initialize empty Vector of tuples to store image positions
   image_position::Vector{NTuple{2, Real}} = []
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
    get_image(lens::AbstractLens, θx::T, θy::T, 
              adis::Real, 
              β::T) where T <: Matrix{<:Real} --> Matrix{<:Real}
Calculate lensed  image map for an extended source given a lens model. This function employs 
inverse ray shooting (IRS) to construct the image plane map.

# Arguments
- `lens`: Lens model.
- `θx`  : x-grid (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-grid (in ``\\rm \\mathbf{arcseconds}``).
- `adis`: Distance ratio (i.e., ``D_{ds}/D_s``).
- `β`   : Either point source position (in ``\\rm \\mathbf{arcseconds}``) or source intensity map.

# Returns
- `img` : Image positions (in ``\\rm \\mathbf{arcseconds}``).
"""
function get_image(lens::AbstractLens, θx::T, θy::T, adis::Float64, β::T) where T <: Matrix{<:Real}
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
         pixel_x = round(Int64, βx / pixel_h + 0.5*nx + 0.5)
         pixel_y = round(Int64, βy / pixel_h + 0.5*ny + 0.5)

         # make sure pixel position is within bounds
         if (1 <= pixel_x <= nx) && (1 <= pixel_y <= ny)
            image_map[i, j] = β[pixel_x, pixel_y]
         end
      end
   end
   return image_map
end


"""
    get_critical_curve(θx::T, θy::T, adis::Float64,
                       ψxx::T, ψyy::T, ψxy::T) where T <: Matrix{<:Real} --> Tuple{Vector{Vector{Vector{Float64}}}, Vector{Vector{Vector{Float64}}}}
Calculate critical curves from a deformation tensor that has already been computed. For this 
function, no lens is needed.
 
# Arguments
- `θx`  : x-grid (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-grid (in ``\\rm \\mathbf{arcseconds}``).
- `adis`: Distance ratio (i.e., ``D_{ds}/D_s``).
- `ψxx` : xx-component of the deformation tensor.
- `ψyy` : yy-component of the deformation tensor.
- `ψxy` : xy-component of the deformation tensor.
 
# Returns
- `critical_tan`: Tangential critical curve(s).
- `critical_rad`: Radial critical curve(s).
"""
function get_critical_curve(θx::T, θy::T, adis::Float64, ψxx::T, ψyy::T, ψxy::T) where T <: Matrix{<:Real}
   # Convergence and shear components, scaled out of place so the caller keeps its tensor
   κ  = 0.5 .* adis .* (ψxx .+ ψyy)
   γ1 = 0.5 .* adis .* (ψxx .- ψyy)
   γ2 = adis .* ψxy
   γ  = sqrt.(γ1.^2 .+ γ2.^2)

   # Get the zero eigenvalue contours
   critical_tan = ContourFinder.get_contour(θx, θy, 1.0 .- κ .- γ, 0)
   critical_rad = ContourFinder.get_contour(θx, θy, 1.0 .- κ .+ γ, 0)
 
   return critical_tan, critical_rad
end

"""
    get_critical_curve(lens::AbstractLens, θx::T, θy::T, 
                       adis::Real) where T <: Matrix{<:Real} --> Tuple{Vector{Vector{Vector{Float64}}}, Vector{Vector{Vector{Float64}}}}
Calculate critical curves for a given lens model. This function essentially runs marching squares
algorithm to find the zero eigenvalue contours.

# Arguments
- `lens`: Lens model.
- `θx`  : x-grid (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-grid (in ``\\rm \\mathbf{arcseconds}``).
- `adis`: Distance ratio (i.e., ``D_{ds}/D_s``).

# Returns
- `crit_tan`: Tangential critical curves (in ``\\rm \\mathbf{arcseconds}``).
- `crit_rad`: Radial critical curves (in ``\\rm \\mathbf{arcseconds}``).
"""
function get_critical_curve(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:Real}
   # Get the jacobian components
   ψxx, ψyy, ψxy = get_jacobian(lens, θx, θy)

   return get_critical_curve(θx, θy, adis, ψxx, ψyy, ψxy)
end


# Map a set of critical curves into the source plane
function _map_to_source(lens::AbstractLens, curves::Vector{<:Vector{<:Vector{Float64}}}, adis::Float64)
   caustics = Vector{Vector{Vector{Float64}}}(undef, length(curves))
   for (idx, curve) in enumerate(curves)
      ψ_x, ψ_y = get_deflection(lens, first.(curve), last.(curve))
      src_x = first.(curve) .- adis .* ψ_x
      src_y =  last.(curve) .- adis .* ψ_y
      caustics[idx] = [[x, y] for (x, y) in zip(src_x, src_y)]
   end
   return caustics
end


"""
    get_caustic(lens::AbstractLens, θx::T, θy::T, adis::Float64, 
                ψxx::T, ψyy::T, ψxy::T) where T <: Matrix{<:Real} --> Tuple{Vector{Vector{Vector{Float64}}}, Vector{Vector{Vector{Float64}}}}
Calculate caustics for a given lens model. The function first gets the critical curves and then 
maps them to the source plane using the lens equation. **The functon runs over all tangential 
and radial critical curves.**

# Arguments
- `lens`: Lens model.
- `θx`  : x-grid (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-grid (in ``\\rm \\mathbf{arcseconds}``).
- `adis`: Distance ratio (i.e., ``D_{ds}/D_s``).
- `ψxx` : xx-component of the deformation tensor.
- `ψyy` : yy-component of the deformation tensor.
- `ψxy` : xy-component of the deformation tensor.

# Returns
- `caustics_tan`: Tangential caustics (in ``\\rm \\mathbf{arcseconds}``).
- `caustics_rad`: Radial caustics (in ``\\rm \\mathbf{arcseconds}``).
"""
function get_caustic(lens::AbstractLens, θx::T, θy::T, adis::Float64, 
                     ψxx::T, ψyy::T, ψxy::T) where T <: Matrix{<:Real}
   # Calculate critical curves
   critical_tan, critical_rad = get_critical_curve(θx, θy, adis, ψxx, ψyy, ψxy)

   # Map the critical curves to the source plane
   caustics_tan = _map_to_source(lens, critical_tan, adis)
   caustics_rad = _map_to_source(lens, critical_rad, adis)

   return caustics_tan, caustics_rad
end

"""
    get_caustic(lens::AbstractLens, θx::T, θy::T, 
                adis::Real) where T <: Matrix{<:Real} --> Tuple{Vector{Vector{Vector{Float64}}}, Vector{Vector{Vector{Float64}}}}
Calculate caustics for a given lens model. The function first gets the critical curves and then 
maps them to the source plane using the lens equation. **The functon runs over all tangential 
critical curves.**

# Arguments
- `lens`: Lens model.
- `θx`  : x-grid (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-grid (in ``\\rm \\mathbf{arcseconds}``).
- `adis`: Distance ratio (i.e., ``D_{ds}/D_s``).

# Returns
- `caus_tan`: Tangential caustic curves (in ``\\rm \\mathbf{arcseconds}``).
- `caus_rad`: Radial caustic curves (in ``\\rm \\mathbf{arcseconds}``).
"""
function get_caustic(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:Real}
   # Get deformation tensor components
   ψxx, ψyy, ψxy = get_jacobian(lens, θx, θy)
   
   return get_caustic(lens, θx, θy, adis, ψxx, ψyy, ψxy)
end


"""
    get_critical_area(lens::AbstractLens, θx::T, θy::T, 
                      adis::Real) where T <: Matrix{<:Real} --> Real
Calculate the total angular area enclosed by tangential critical curve(s). The function runs shoelace
algorithm to calculate the area.

# Arguments
- `lens`: Lens model.
- `θx`  : x-grid (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-grid (in ``\\rm \\mathbf{arcseconds}``).
- `adis`: Distance ratio (i.e., ``D_{ds}/D_s``).

# Returns
- `area`: Total angular area enclosed by tangential critical curves (in ``\\rm \\mathbf{arcseconds^2}``).
"""
function get_critical_area(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:Real}
   # Get tangential critical curves
   critical_tan, _ = get_critical_curve(lens, θx, θy, adis)

   # Run a loop over all tangential critical curves
   area = 0.0
   for curve in critical_tan
      area = area + PolygonOps.shoelace(curve)
   end
   return area
end


"""
    get_einstein_angle(lens::AbstractLens, θx::T, θy::T, 
                       adis::Real) where T <: Matrix{<:Real} --> Real
Calculate the Einstein radius (i.e., ``θ_E``) for an arbitrary lens model, which is defined as,
```math
θ_E = \\sqrt{\\frac{A_{\\rm critical}}{π}},
```
where ``A_{\\rm critical}`` is the total angular area enclosed by the tangential critical curve(s).
Again, it is important to note that we are running over all tangential critical curves.

# Arguments
- `lens`: Lens model.
- `θx`  : x-grid (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-grid (in ``\\rm \\mathbf{arcseconds}``).
- `adis`: Distance ratio (i.e., ``D_{ds}/D_s``).

# Returns
- `θ_E` : Einstein angle (in ``\\rm \\mathbf{arcseconds}``).
"""
function get_einstein_angle(lens::AbstractLens, θx::T, θy::T, adis::Float64) where T <: Matrix{<:Real}
   return sqrt(get_critical_area(lens, θx, θy, adis) / π)
end


function _is_closed_curve(curve::Vector{<:Vector{<:Real}})
   if length(curve) < 4
      return false
   end 
   return curve[1] == curve[end]
end

function _caustic_winding(caustics::Vector{<:Vector{<:Vector{<:Real}}}, βx::Real, βy::Real; 
                          verbose::Bool = true)
   point = [Float64(βx), Float64(βy)]
 
   w = Int64[]
   for curve in caustics
      if !_is_closed_curve(curve)
         if verbose
            @warn "Skipping an open caustic; enlarge the grid (θx, θy) so that all " *
                   "critical curves are enclosed by it."
         end
         continue
      end
      push!(w, PolygonOps.winding_number(point, curve))
   end
   return w
end

"""
    get_image_multiplicity(caustics::Vector{<:Vector{<:Vector{<:Real}}}, βx::Real, βy::Real;
                           n_far::Int    = 1,
                           verbose::Bool = true) --> Int
Calculate the number of images of a point source located at ``(β_x, β_y)`` from the caustics alone.

A source lying far outside every caustic has a single image. Moving in towards the lens, the number
of images changes by ``\\pm 2`` at every caustic crossing, so the multiplicity anywhere in the source
plane follows from counting the crossings along the way in,
 
```math
N(\\pmb{β}) = N_\\infty + 2 \\sum_i \\left| w_i(\\pmb{β}) \\right|,
```
 
where ``w_i`` is the winding number of the ``i``-th caustic about ``\\pmb{β}`` and
``N_\\infty = 1``. Using winding numbers rather than a raw crossing count takes care of the sign of
each crossing automatically, so the answer does not depend on the direction one comes in from and
stays correct for caustics that self-intersect (swallowtails) or that overlap each other (e.g. a
naked cusp configuration, where the tangential caustic sticks out of the radial one). Note that the 
caustics are the only input: the image plane is not needed, and the direction in which
the caustics happen to be traced is irrelevant because only ``|w_i|`` enters.

!!! warning
    The lens is assumed to be non-singular, so that ``N`` is odd everywhere 
    [1981ApJ...244L...1B](@cite). A singular lens (e.g.[`PointLens`](PointLens.md), 
    [`SISLens`](SISLens.md)) destroys the central image and produces 
    cuts [1987ApJ...312...22K](@cite) that also change the multiplicity but are not images of 
    critical curves, and hence are invisible here. For such models the count returned is that of 
    the corresponding non-singular lens.
 
# Arguments
- `caustics`: Caustic curves in the source plane, e.g. from [`get_caustic`](@ref).
- `βx`      : x-coordinate of the source (in ``\\rm \\mathbf{arcseconds}``).
- `βy`      : y-coordinate of the source (in ``\\rm \\mathbf{arcseconds}``).
 
# Keyword Arguments
- `n_far     = 1`: Number of images of a source lying far outside all caustics.
- `verbose   = true`: Warn when an open (i.e., unclosed) caustic is discarded.
 
# Returns
- `N`: Number of images of the point source.
"""
function get_image_multiplicity(caustics::Vector{<:Vector{<:Vector{<:Real}}}, βx::Real, βy::Real;
                         n_far::Int    = 1,
                         verbose::Bool = true)
   # Get the winding numbers for all tangential caustics
   w = _caustic_winding(caustics, βx, βy; verbose = verbose)
   return n_far + 2 * sum(abs, w; init = 0)
end


"""
    get_image_multiplicity(lens::AbstractLens, θx::T, θy::T, adis::Float64;
                           source::Union{Nothing, NTuple{2, Real}} = nothing,
                           n_far::Int    = 1,
                           verbose::Bool = true) where T <: Matrix{<:Real}
Calculate the number of images for a given lens model. The function first gets the tangential and
radial caustics with [`get_caustic`](@ref) and then counts the caustic crossings on the way in from
infinity; see the caustic based method above for the details.
 
The shape of the output follows the `source` keyword argument,
 
- `source = (βx, βy)` ``\\Rightarrow`` number of images of that point source,
- `source = nothing`  ``\\Rightarrow`` multiplicity map over the whole grid.
 
In the latter case the grid `(θx, θy)` is reused as the source plane grid, so that `N[i, j]` is the
number of images of a point source sitting at `(θx[i, j], θy[i, j])`. This needs no extra input: a
grid that contains every critical curve, as it must, also contains every caustic, since the
deflection pulls the critical curves inwards. For a multiplicity map on a different (e.g. zoomed in)
grid, get the caustics with [`get_caustic`](@ref) and pass them to the caustic based method above.
 
The image plane grid must be large enough to contain **all** critical curves, otherwise the
corresponding caustics are open and are dropped from the count.
 
# Arguments
- `lens`: Lens model.
- `θx`  : x-grid (in ``\\rm \\mathbf{arcseconds}``).
- `θy`  : y-grid (in ``\\rm \\mathbf{arcseconds}``).
- `adis`: Distance ratio (i.e., ``D_{ds}/D_s``).
 
# Keyword Arguments
- `source = nothing`: Source angular position (in ``\\rm \\mathbf{arcseconds}``).
   * If ``\\mathbf{nothing}``, then the multiplicity is calculated at every pixel of the grid.
- `n_far  = 1`: Number of images of a source lying far outside all caustics.
 
# Returns
- `N`: Number of images, either for the requested source position or over the whole grid.
"""
function get_image_multiplicity(lens::AbstractLens, θx::T, θy::T, adis::Float64;
                                source::Union{Nothing, NTuple{2, Real}} = nothing,
                                n_far::Int    = 1,
                                verbose::Bool = true) where T <: Matrix{<:Real}
   # Get tangential and rafial caustics
   caustics_tan, caustics_rad = get_caustic(lens, θx, θy, adis)
   caustics = vcat(caustics_tan, caustics_rad)

   # Drop the open caustics so that we do not get warnings later
   n_open = count(curve -> !_is_closed_curve(curve), caustics)
   if n_open > 0
      @warn "Skipping $(n_open) open caustic(s); enlarge the grid (θx, θy) so that all critical " *
            "curves are enclosed by it."
      caustics = [curve for curve in caustics if _is_closed_curve(curve)]
   end

   # If no source position is provided, run over every pixel of the grid
   if isnothing(source)
      N = Matrix{Int}(undef, size(θx))
      Threads.@threads for k in eachindex(θx)
         N[k] = get_image_multiplicity(caustics, θx[k], θy[k]; n_far = n_far, verbose = false)
      end
      return N
   elseif source isa NTuple{2, Real}
      return get_image_multiplicity(caustics, source[1], source[2]; n_far = n_far, verbose = false)
   else
      throw(ArgumentError("`source` must either be `nothing` or a tuple of coordinates (βx, βy), " *
                          "got a $(typeof(source))."))
   end
end


"""
    get_radial_profile(kappa::T, θx::T, θy::T; 
                       origin::Union{Tuple{Float64, Float64}, Nothing} = nothing, 
                       n_bin::Int64     = 50, 
                       bin_type::Symbol = :log) --> Tuple(Vector{Real}, Vector{Real}, Vector{Real})
Calculate the radial convergence profile. To convert it into physicsal units, multiply it by the 
critical density ``Σ_{\\rm cr}``.


# Arguments
- `kappa`: Convergence map.
- `θx`   : x-grid (in ``\\rm \\mathbf{arcseconds}``).
- `θy`   : y-grid (in ``\\rm \\mathbf{arcseconds}``).

# Keyword Arguments
- `origin   = nothing`: Center (in ``\\rm \\mathbf{arcseconds}``).
- `n_bin    = 50`: Number of radial bins.
- `bin_type = :log`: Type of radial binning (i.e., ``:log`` or ``:linear``).

# Returns
- `centers`: Radial bin centers (in ``\\rm \\mathbf{arcseconds}``).
- `profile`: Radial profile (in ``\\rm \\mathbf{arcseconds}``).
- `edges`  : Radial bin edges (in ``\\rm \\mathbf{arcseconds}``).
"""
function get_radial_profile(kappa::T, θx::T, θy::T; 
                            origin::Union{Tuple{Float64, Float64}, Nothing} = nothing, 
                            n_bin::Int64     = 50, 
                            bin_type::Symbol = :log) where T <: Matrix{<:Real}
   # Get the radial 1D grid based the input 2D grid and origin. If origin is not provided, then the 
   # largest value pixel will be chosen as the center.
   if origin === nothing
      px, py = argmax(kappa).I
      x0, y0 = θx[px, py], θy[px, py]
   else
      x0, y0 = origin
   end

   # Get the maximum possible radial length
   θ = @. sqrt((θx - x0)^2 + (θy - y0)^2)
   θ_min = abs(θx[2, 1] - θx[1, 1])
   θ_max = maximum(θ)

   # Create radial grid
   if bin_type == :log
      # Bin edges
      edges = 10.0 .^ range(start=log10(θ_min), stop=log10(θ_max), length=n_bin)
      
      # Geometric mean as bin center for log binning
      centers = sqrt.(edges[1:end-1] .* edges[2:end])
   else
      # Bin edges
      edges = collect(range(start=θ_min, stop=θ_max, length=n_bin))
      
      # Arithmetic mean of bin edges
      centers = (edges[1:end-1] .+ edges[2:end]) ./ 2.0
   end

   # Efficient binning with StatsBase
   kappa_clean = replace(kappa, NaN => 0.0)
   h_sum = fit(Histogram, vec(θ), weights(vec(kappa_clean)), edges)
   h_count = fit(Histogram, vec(θ), edges)
   profile = [count > 0 ? s / count : 0.0 for (s, count) in zip(h_sum.weights, h_count.weights)]
   return centers, profile, edges
end


"""
    get_mass_profile(kappa::T, θx::T, θy::T, 
                     D_d::Real; 
                     origin::Union{Tuple{Float64, Float64}, Nothing} = nothing, 
                     n_bin::Int64     = 50, 
                     bin_type::Symbol = :log) where T <: Matrix{<:Real} --> Tuple{Vector{Float64}, Vector{Float64}}
Calculate the cumulative mass enclosed within a given angular radius ``θ`` for a given lens model. 
While converting the input convergence map into the physical units, it is assumed that source is 
at infinity (i.e., ``a_{\\rm dis} = 1``). Hence, if the input convergence is for any finite source 
redshift, then divide the output mass values by ``a_{\\rm dis}``.

# Arguments
- `kappa`: Convergence map.
- `θx`   : x-grid (in ``\\rm \\mathbf{arcseconds}``).
- `θy`   : y-grid (in ``\\rm \\mathbf{arcseconds}``).
- `D_d`  : Angular diameter distance to the lens (in ``\\rm \\mathbf{meters}``).

# Keyword Arguments
- `origin   = nothing`: Center (in ``\\rm \\mathbf{arcseconds}``).
   - If ``\\rm \\mathbf{nothing}``, then the center is set to be the location of the maximum value 
      of the convergence map.
- `n_bin    = 50`: Number of radial bins.
- `bin_type = :log`: Type of radial binning (i.e., `:log` or `:linear`).

# Returns
- `centers`: Radial centers (in ``\\rm \\mathbf{arcseconds}``).
- `mass`   : Cumulative mass (in ``\\rm \\mathbf{M_\\odot}``).
"""
function get_mass_profile(kappa::T, θx::T, θy::T, D_d::Float64; 
                          origin::Union{Tuple{Float64, Float64}, Nothing} = nothing, 
                          n_bin::Int64     = 50, 
                          bin_type::Symbol = :log) where T <: Matrix{<:Real}
   # Calculate convergence radial profile (in units of Σ_cr)
   centers, profile, edges = get_radial_profile(kappa, θx, θy; origin=origin, n_bin=n_bin, bin_type=bin_type)

   # Get critical density for a source at z = ∞
   Σ_cr = get_critical_density(D_d, 1.0, 1.0; unit=:msun_arcsec2)

   # Mass profile unit conversion (Σ_cr --> M⊙/arcsec^2)
   profile .= profile .* Σ_cr

   # Calcualte bin areas
   θ_in  = edges[1:end-1]
   θ_out = edges[2:end]
   bin_areas = @. π * (θ_out^2 - θ_in^2)

   # Cumulative mass sum
   mass = cumsum(profile .* bin_areas)
   return edges[2:end], mass
end


"""
    shear_cartesian2polar(γ1::Real, γ2::Real) --> Tuple{Real, Real}
Converts the Cartesian components of the shear (i.e., ``γ_1`` and ``γ_2``) to polar components
(i.e., ``γ`` and ``φ``) using the relations,
```math
\\begin{align*}
γ &= \\sqrt{γ_1^2 + γ_2^2}, \\\\
φ &= \\frac{1}{2} \\tan^{-1}\\left(\\frac{γ_2}{γ_1}\\right).
\\end{align*}
```

# Arguments
- `γ1`: Cartesian component of the shear (i.e., ``γ_1``).
- `γ2`: Cartesian component of the shear (i.e., ``γ_2``).

# Returns
- `γ` : Polar component of the shear (i.e., ``γ``).
- `φ` : Polar component of the shear (i.e., ``φ`` in ``\\rm \\mathbf{degrees}``).
"""
function shear_cartesian2polar(γ1::Real, γ2::Real)
   return hypot(γ1, γ2), 0.5 * rad2deg(atan(γ2, γ1))
end


"""
    shear_polar2cartesian(γ::Real, phi::Real) --> Tuple{Real, Real}
Converts the polar components of the shear (i.e., ``γ`` and ``φ``) to Cartesian components
(i.e., ``γ_1`` and ``γ_2``) using the relations,
```math
\\begin{align*}
γ_1 &= γ \\cos(2φ), \\\\
γ_2 &= γ \\sin(2φ).
\\end{align*}
```

# Arguments
- `γ` : Polar component of the shear (i.e., ``γ``).
- `φ` : Polar component of the shear (i.e., ``φ`` in ``\\rm \\mathbf{degrees}``).

# Returns
- `γ1`: Cartesian component of the shear (i.e., ``γ_1``).
- `γ2`: Cartesian component of the shear (i.e., ``γ_2``).
"""
function shear_polar2cartesian(γ::Real, phi::Real)
   return γ * cos(2.0 * deg2rad(phi)), γ * sin(2.0 * deg2rad(phi))
end


"""
    ellipticity_cartesian2polar(e1::Real, e2::Real) --> Tuple{Real, Real}
Converts the Cartesian components of the ellipticity (i.e., ``e_1`` and ``e_2``) to polar components
(i.e., ``e`` and ``φ``) using the relations,
```math
\\begin{align*}
e &= \\sqrt{e_1^2 + e_2^2}, \\\\
φ &= \\frac{1}{2} \\tan^{-1}\\left(\\frac{e_2}{e_1}\\right).
\\end{align*}
```

# Arguments
- `e1`: Cartesian component of the ellipticity (i.e., ``e_1``).
- `e2`: Cartesian component of the ellipticity (i.e., ``e_2``).

# Returns
- `e` : Polar component of the ellipticity (i.e., ``e``).
- `φ` : Polar component of the ellipticity (i.e., ``φ`` in ``\\rm \\mathbf{degrees}``).
"""
function ellipticity_cartesian2polar(e1::Real, e2::Real)
   return hypot(e1, e2), 0.5 * rad2deg(atan(e2, e1))
end


"""
    ellipticity_polar2cartesian(e::Real, phi::Real) --> Tuple{Real, Real}
Converts the polar components of the ellipticity (i.e., ``e`` and ``φ``) to Cartesian components
(i.e., ``e_1`` and ``e_2``) using the relations,
```math
\\begin{align*}
e_1 &= e \\cos(2φ), \\\\
e_2 &= e \\sin(2φ).
\\end{align*}
```

# Arguments
- `e` : Polar component of the ellipticity (i.e., ``e``).
- `φ` : Polar component of the ellipticity (i.e., ``φ`` in ``\\rm \\mathbf{degrees}``).

# Returns
- `e1`: Cartesian component of the ellipticity (i.e., ``e_1``).
- `e2`: Cartesian component of the ellipticity (i.e., ``e_2``).
"""
function ellipticity_polar2cartesian(e::Real, phi::Real)
   return e * cos(2.0 * deg2rad(phi)), e * sin(2.0 * deg2rad(phi))
end


# --------------------------------------------------------------------------------------------------
# Parameter functions for various lenses
# --------------------------------------------------------------------------------------------------
"""
    parameter_NFWLens(; cosmology::Cosmology.AbstractCosmology = nothing, 
                        z_d::Real  = NaN, 
                        mass::Real = NaN, 
                        x_s::Real  = NaN, 
                        c::Real    = NaN) --> NamedTuple
Calculate parameters for NFW lens. The function would either need the concentration `c` or the 
scale radius `x_s`. **If both are provided, `c` will be used to calculate `x_s` and the input `x_s` 
will be overwritten.**

# Keyword Arguments
- `cosmology = nothing`: Cosmology object.
- `z_d       = NaN`: Redshift of the lens.
- `mass      = NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s       = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `c         = NaN`: Concentration of the lens.

# Returns
- Named tuple of NFW lens parameters.
"""
function parameter_NFWLens(; cosmology::Cosmology.AbstractCosmology = nothing, 
                             z_d::Real  = NaN, 
                             mass::Real = NaN, 
                             x_s::Real  = NaN, 
                             c::Real    = NaN)
   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift (in kg/m^3)
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass * MASS_SUN / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide concentration (c) or scale radius (x_s) in **parameter_NFWLens**."))
   end
   # 3D characteristic density
   mass_c = log(1.0 + c) - (c / (1.0 + c))
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / mass_c

   # 2D (normalized) characteristic density
   k_s = ρ_s * D_d * x_s * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)
   return (mass=mass, rho_s=ρ_s, k_s=k_s, c=c, x_s=x_s)
end


"""
    parameter_gNFWLens(; cosmology::Cosmology.AbstractCosmology = nothing, 
                         z_d::Real  = NaN, 
                         mass::Real = NaN, 
                         x_s::Real  = NaN, 
                         c::Real    = NaN, 
                         n::Real    = 1.0) --> NamedTuple
Calculate parameters for gNFW lens. The function would either need the concentration `c` or the 
scale radius `x_s`. **If both are provided, `c` will be used to calculate `x_s` and the input `x_s` 
will be overwritten.**

# Keyword Arguments
- `cosmology = nothing`: Cosmology object.
- `z_d       = NaN`: Redshift of the lens.
- `mass      = NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s       = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `c         = NaN`: Concentration of the lens.
- `n         = 1.0`: Slope parameter of the lens.

# Returns
- Named tuple of gNFW lens parameters.
"""
function parameter_gNFWLens(; cosmology::Cosmology.AbstractCosmology=nothing, 
                              z_d::Real  = NaN, 
                              mass::Real = NaN, 
                              x_s::Real  = NaN, 
                              c::Real    = NaN, 
                              n::Real    = 1.0)
   # Check for valid slope parameter
   if !(0.0 < n < 2.0)
      throw(ArgumentError("Slope parameter outside allowed range n ∈ (0, 2) in **parameter_gNFWLens**."))
   end

   # Integrand function for mass calculation
   function integrand(x::Real, α::Real)
      return x^(2.0 - α) / (1.0 + x)^(3.0 - α)
   end

   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)
   
   # Critical density at the lens redshift
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass * MASS_SUN / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide at least c or x_s in **parameter_gNFWLens**."))
   end
   mass_c, _ = quadgk(x -> integrand(x, n), 0, c)
   
   # 3D characteristic density
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / mass_c

   # 2D (normalized) characteristic density
   k_s = ρ_s * D_d * x_s * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)
   return (mass=mass, rho_s=ρ_s, k_s=k_s, c=c, x_s=x_s, n=n)
end


"""
    parameter_EinastoLens(; cosmology::Cosmology.AbstractCosmology=nothing, 
                            z_d::Real  = NaN, 
                            mass::Real = NaN, 
                            x_s::Real  = NaN, 
                            c::Real    = NaN, 
                            n::Real    = 0.2) --> NamedTuple
Calculate parameters of an Einasto lens model. The function would either need the concentration `c` 
or the scale radius `x_s`. **If both are provided, `c` will be used to calculate `x_s` and the input 
`x_s` will be overwritten.**

# Keyword Arguments
- `cosmology = nothing`: Cosmology object.
- `z_d       = NaN`: Redshift of the lens.
- `mass      = NaN`: Mass of the lens (in ``\\rm \\mathbf{M_\\odot}``).
- `x_s       = NaN`: Scale radius (in ``\\rm \\mathbf{arcseconds}``).
- `c         = NaN`: Concentration of the lens.
- `n         = 0.2`: Slope parameter of the lens.

# Returns
- Named tuple of Einasto lens parameters.
"""
function parameter_EinastoLens(; cosmology::Cosmology.AbstractCosmology=nothing, 
                                 z_d::Real  = NaN, 
                                 mass::Real = NaN, 
                                 x_s::Real  = NaN, 
                                 c::Real    = NaN, 
                                 n::Real    = 0.2)
   # Overdensity value
   Δ_z = 200.0

   # ADD to the lens
   D_d  = Cosmology.angular_diameter_distance(cosmology, 0.0, z_d)

   # Critical density at the lens redshift
   ρ_cz = Cosmology.rho_cz(cosmology, z_d)

   # Virial radius of the lens (in ANGLE_ARCSEC)
   θ_vir = (3.0 * mass * MASS_SUN / 4.0 / pi / Δ_z / ρ_cz)^(1.0/3.0) / D_d / ANGLE_ARCSEC

   # Check if concentration is given
   if isfinite(c)
      x_s = θ_vir / c
   elseif isfinite(x_s)
      c = θ_vir / x_s
   else
      throw(ArgumentError("Provide at least c or x_s in **parameter_EinastoLens**."))
   end
   Pax, _ = gamma_inc(3.0 / n, (2.0 / n) * c^n)
   mass_e = (1.0 / n) * (n / 2.0)^(3.0 / n) * gamma(3.0 / n) * Pax
   
   # 3D characteristic density
   ρ_s = (Δ_z / 3.0) * ρ_cz * c^3 / mass_e

   # 2D (normalized) characteristic density
   k_s = ρ_s * D_d * x_s * ANGLE_ARCSEC / (CONST_C^2 / 4.0 / pi / CONST_G / D_d)
   return (mass=mass, rho_s=ρ_s, k_s=k_s, c=c, x_s=x_s, n=n)
end


# --------------------------------------------------------------------------------------------------
# -------------------- Potential functions for specific lens models --------------------------------
# --------------------------------------------------------------------------------------------------
@inline function potential_helper!(ψ::T, lens::init_PointLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PointLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass)
end

@inline function potential_helper!(ψ::T, lens::init_PlummerLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PlummerLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function potential_helper!(ψ::T, lens::init_SISLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return SISLens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.v_d)
end

@inline function potential_helper!(ψ::T, lens::init_NSISPLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return NSISPLens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s)
end

@inline function potential_helper!(ψ::T, lens::init_NSISMDLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return NSISMDLens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s)
end

@inline function potential_helper!(ψ::T, lens::init_GaussianLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return GaussianLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function potential_helper!(ψ::T, lens::init_SersicLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return SersicLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_e, lens.n)
end

@inline function potential_helper!(ψ::T, lens::init_PixelLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PixelLens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.kappa, lens.pixel_size)
end

@inline function potential_helper!(ψ::T, lens::init_ExternalEffects, θx::T, θy::T) where T <: Union{Real, ROA}
   return ExternalEffects.potential!(ψ, θx, θy, lens.kappa, lens.gamma, lens.angle)
end

@inline function potential_helper!(ψ::T, lens::init_ExternalEffects3, θx::T, θy::T) where T <: Union{Real, ROA}
   return ExternalEffects3.potential!(ψ, θx, θy, lens.delta, lens.angle)
end

@inline function potential_helper!(ψ::T, lens::init_Multipole, θx::T, θy::T) where T <: Union{Real, ROA}
   return Multipole.potential!(ψ, θx, θy, lens.delta, lens.angle, lens.m, lens.n)
end

@inline function potential_helper!(ψ::T, lens::init_PIEPLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PIEPLens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.eps, lens.pa)
end

@inline function potential_helper!(ψ::T, lens::init_SIELens, θx::T, θy::T) where T <: Union{Real, ROA}
   return SIELens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.eps, lens.pa)
end

@inline function potential_helper!(ψ::T, lens::init_PJELens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PJELens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.x_t, lens.eps, lens.pa)
end

@inline function potential_helper!(ψ::T, lens::init_HernquistLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return HernquistLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function potential_helper!(ψ::T, lens::init_NFWLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return NFWLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s)
end

@inline function potential_helper!(ψ::T, lens::init_tNFWLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return tNFWLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.x_t)
end

@inline function potential_helper!(ψ::T, lens::init_gNFWLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return gNFWLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.n)
end

@inline function potential_helper!(ψ::T, lens::init_EinastoLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return EinastoLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.n)
end

@inline function potential_helper!(ψ::T, lens::init_aHernquistLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return aHernquistLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.eps, lens.pa)
end

@inline function potential_helper!(ψ::T, lens::init_aNFWLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return aNFWLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.eps, lens.pa)
end

@inline function potential_helper!(ψ::T, lens::init_eHernquistMDLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return eHernquistMDLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.eps, lens.pa)
end

@inline function potential_helper!(ψ::T, lens::init_eNFWMDLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return eNFWMDLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.eps, lens.pa)
end

@inline function potential_helper!(ψ::T, lens::init_MultiPlummerLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return MultiPlummerLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.n)
end

@inline function potential_helper!(ψ::T, lens::init_MultiGaussianLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return MultiGaussianLens.potential!(ψ, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.n)
end

@inline function potential_helper!(ψ::T, lens::init_MultiPixelLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return MultiPixelLens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.kappa, lens.pixel_size, lens.n)
end

@inline function potential_helper!(ψ::T, lens::init_MultiPJELens, θx::T, θy::T) where T <: Union{Real, ROA}
   return MultiPJELens.potential!(ψ, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.x_t, lens.eps, lens.pa, lens.n)
end


# --------------------------------------------------------------------------------------------------
# -------------------- Deflection functions for specific lens models -------------------------------
# --------------------------------------------------------------------------------------------------
@inline function deflection_helper!(ψx::T, ψy::T, lens::init_PointLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PointLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_PlummerLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PlummerLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_SISLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return SISLens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.v_d)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_NSISPLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return NSISPLens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_NSISMDLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return NSISMDLens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_GaussianLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return GaussianLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_SersicLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return SersicLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_e, lens.n)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_PixelLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PixelLens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.kappa, lens.pixel_size)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_ExternalEffects, θx::T, θy::T) where T <: Union{Real, ROA}
   return ExternalEffects.deflection!(ψx, ψy, θx, θy, lens.kappa, lens.gamma, lens.angle)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_ExternalEffects3, θx::T, θy::T) where T <: Union{Real, ROA}
   return ExternalEffects3.deflection!(ψx, ψy, θx, θy, lens.delta, lens.angle)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_Multipole, θx::T, θy::T) where T <: Union{Real, ROA}
   return Multipole.deflection!(ψx, ψy, θx, θy, lens.delta, lens.angle, lens.m, lens.n)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_PIEPLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PIEPLens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.eps, lens.pa)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_SIELens, θx::T, θy::T) where T <: Union{Real, ROA}
   return SIELens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.eps, lens.pa)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_PJELens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PJELens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.x_t, lens.eps, lens.pa)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_HernquistLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return HernquistLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_NFWLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return NFWLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_tNFWLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return tNFWLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.x_t)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_gNFWLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return gNFWLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.n)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_EinastoLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return EinastoLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.n)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_aHernquistLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return aHernquistLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.eps, lens.pa)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_aNFWLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return aNFWLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.eps, lens.pa)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_eHernquistMDLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return eHernquistMDLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.eps, lens.pa)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_eNFWMDLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return eNFWMDLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.eps, lens.pa)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_MultiPlummerLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return MultiPlummerLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.n)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_MultiGaussianLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return MultiGaussianLens.deflection!(ψx, ψy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.n)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_MultiPixelLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return MultiPixelLens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.kappa, lens.pixel_size, lens.n)
end

@inline function deflection_helper!(ψx::T, ψy::T, lens::init_MultiPJELens, θx::T, θy::T) where T <: Union{Real, ROA}
   return MultiPJELens.deflection!(ψx, ψy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.x_t, lens.eps, lens.pa, lens.n)
end


# --------------------------------------------------------------------------------------------------
# -------------------- Deformation tensor for various lens models ----------------------------------
# --------------------------------------------------------------------------------------------------
@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_PointLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PointLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_PlummerLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PlummerLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_SISLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return SISLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.v_d)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_NSISPLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return NSISPLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_NSISMDLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return NSISMDLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_GaussianLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return GaussianLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_SersicLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return SersicLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_e, lens.n)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_PixelLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PixelLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.kappa, lens.pixel_size)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_ExternalEffects, θx::T, θy::T) where T <: Union{Real, ROA}
   return ExternalEffects.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.kappa, lens.gamma, lens.angle)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_ExternalEffects3, θx::T, θy::T) where T <: Union{Real, ROA}
   return ExternalEffects3.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.delta, lens.angle)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_Multipole, θx::T, θy::T) where T <: Union{Real, ROA}
   return Multipole.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.delta, lens.angle, lens.m, lens.n)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_PIEPLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PIEPLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.eps, lens.pa)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_SIELens, θx::T, θy::T) where T <: Union{Real, ROA}
   return SIELens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.eps, lens.pa)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_PJELens, θx::T, θy::T) where T <: Union{Real, ROA}
   return PJELens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.x_t, lens.eps, lens.pa)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_HernquistLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return HernquistLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_NFWLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return NFWLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_tNFWLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return tNFWLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.x_t)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_gNFWLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return gNFWLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.n)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_EinastoLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return EinastoLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.n)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_aHernquistLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return aHernquistLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.eps, lens.pa)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_aNFWLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return aNFWLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.eps, lens.pa)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_eHernquistMDLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return eHernquistMDLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.eps, lens.pa)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_eNFWMDLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return eNFWMDLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.k_s, lens.x_s, lens.eps, lens.pa)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_MultiPlummerLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return MultiPlummerLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.n)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_MultiGaussianLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return MultiGaussianLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.D_d, lens.x_c, lens.y_c, lens.mass, lens.x_s, lens.n)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_MultiPixelLens, θx::T, θy::T) where T <: Union{Real, ROA}
   return MultiPixelLens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.kappa, lens.pixel_size, lens.n)
end

@inline function jacobian_helper!(ψxx::T, ψyy::T, ψxy::T, lens::init_MultiPJELens, θx::T, θy::T) where T <: Union{Real, ROA}
   return MultiPJELens.jacobian!(ψxx, ψyy, ψxy, θx, θy, lens.x_c, lens.y_c, lens.v_d, lens.x_s, lens.x_t, lens.eps, lens.pa, lens.n)
end


end