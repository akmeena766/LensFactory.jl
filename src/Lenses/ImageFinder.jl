module ImageFinder


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ...LFUtils: ContourFinder
using ...LFUtils: IntersectionFinder


# --------------------------------------------------------------------------------------------------
# Various function to export
# --------------------------------------------------------------------------------------------------
export find_images_contour
export find_images_adaptive


# --------------------------------------------------------------------------------------------------
# Find images using contour method
# --------------------------------------------------------------------------------------------------
"""
    find_images_contour(θx::T, θy::T, βx::T, βy::T, 
                        β::NTuple{2, Real}) where T <: AbstractMatrix{<:Real} --> Vector{NTuple{2, Float64}}
Find the solutions ``\\pmb{θ}`` of the lens equation ``\\pmb{β}(\\pmb{θ}) = \\pmb{β}`` on a grid
as the intersection points of the zero contours of ``β_x(\\pmb{θ}) - β_x`` and
``β_y(\\pmb{θ}) - β_y``. The contours are traced with marching squares and intersected as
piecewise-linear curves, so the accuracy is set by the grid resolution.

# Arguments
- `θx`, `θy`: Image-plane grid (same layout as `Lenses.get_meshgrid`).
- `βx`, `βy`: Grid mapped to the source plane.
- `β`: Source position.

# Returns
- `Vector{NTuple{2, Float64}}`: Image positions as `(x, y)` tuples.
"""
function find_images_contour(θx::T, θy::T, βx::T, βy::T, β::NTuple{2, Real}) where T <: AbstractMatrix{<:Real}
   # Check if the input matrices have the same axes
   if !(axes(θx) == axes(θy) == axes(βx) == axes(βy))
      throw(ArgumentError("θx, θy, βx, and βy must have the same axes."))
   end

   # Zero contours of both components of the lens equation
   RXC = ContourFinder.get_contour(θx, θy, βx .- β[1], 0.0)
   RYC = ContourFinder.get_contour(θx, θy, βy .- β[2], 0.0)

   # Their intersections are the images
   images = NTuple{2, Float64}[]
   for contour_1 in RXC, contour_2 in RYC
      # A contour needs at least one segment to intersect
      (length(contour_1) < 2 || length(contour_2) < 2) && continue
      for point in IntersectionFinder.get_intersection(first.(contour_1), last.(contour_1), first.(contour_2), last.(contour_2))
         push!(images, (Float64(point[1]), Float64(point[2])))
      end
   end
   return images
end


# --------------------------------------------------------------------------------------------------
# Find images using adaptive method
# --------------------------------------------------------------------------------------------------
# A point in the image or source plane
const _Pt = NTuple{2, Float64}

# A triangle struct in the image plane (`θ`), its vertices mapped to the source plane (`β`), the 
# estimated error in the source plane (`δ`) of approximating the lens mapping inside the triangle 
# by the linear map through its vertices, and its refinement `level`.
struct _Tri
   θ::NTuple{3, _Pt}
   β::NTuple{3, _Pt}
   δ::Float64
   level::Int
end

@inline function _isfinite(p::_Pt)
   return isfinite(p[1]) && isfinite(p[2])
end

@inline function _isfinite(t::NTuple{3, _Pt})
   return _isfinite(t[1]) && _isfinite(t[2]) && _isfinite(t[3])
end

@inline function _finite_or_zero(x::Float64)
   return isfinite(x) ? x : 0.0
end

@inline function _mid(a::_Pt, b::_Pt)
   return (0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2]))
end

@inline function _dist(a::_Pt, b::_Pt)
   return hypot(a[1] - b[1], a[2] - b[2])
end

@inline function _centroid(t::NTuple{3, _Pt})
   return ((t[1][1] + t[2][1] + t[3][1]) / 3.0, (t[1][2] + t[2][2] + t[3][2]) / 3.0)
end

@inline function _edge(t::NTuple{3, _Pt})
   return max(_dist(t[1], t[2]), _dist(t[2], t[3]), _dist(t[3], t[1]))
end

# Floor of the source-plane residual β(θ) - β = θ - α(θ) - β at θ. Its rounding error scales with 
# the largest terms, i.e. |θ| (≈ |α| for multiple images) and |β|, which matters for cluster-scale 
# lenses where |θ| ~ 100″ while |β| may be small.
@inline function _residual_floor(θ::_Pt, β::_Pt)
   return eps(Float64) * (1.0 + abs(θ[1]) + abs(θ[2]) + abs(β[1]) + abs(β[2]))
end

@inline function _map(lens_map::F, p::_Pt) where F
   bx, by = lens_map(p[1], p[2])
   return (Float64(bx), Float64(by))
end

@inline function _det(lens_jac::G, p::_Pt) where G
   a11, a12, a21, a22 = lens_jac(p[1], p[2])
   return Float64(a11 * a22 - a12 * a21)
end



"""
    find_images_adaptive(lens_map::F, lens_jac::G,
                         θx::T, θy::T, βx::T, βy::T, β::NTuple{2, Real};
                         tol_images::Real = 1.0E-8,
                         tol_refine::Real = 1.0E-3,
                         tol_source::Real = 1.0E-6,
                         max_levels::Int  = 20,
                         max_newton::Int  = 20,
                         pad::Real = 0.1,
                         cfac::Real = 1.0,
                         max_candidates::Int = 100_000) --> Vector{NTuple{2, Float64}}
Find the solutions ``\\pmb{θ}`` of the lens equation ``\\pmb{β}(\\pmb{θ}) = \\pmb{β}`` on a grid
as the intersection points of the zero contours of ``β_x(\\pmb{θ}) - β_x`` and
``β_y(\\pmb{θ}) - β_y``. The contours are traced with marching squares and intersected as
piecewise-linear curves, so the accuracy is set by the grid resolution.

# Arguments
- `θx`, `θy`: Image-plane grid (same layout as `Lenses.get_meshgrid`).
- `βx`, `βy`: Grid mapped to the source plane.
- `β`: Source position.

# Returns
- `Vector{NTuple{2, Float64}}`: Image positions as `(x, y)` tuples.
"""
function find_images_adaptive(lens_map::F, lens_jac::G,
                              θx::T, θy::T, βx::T, βy::T, β::NTuple{2, Real};
                              tol_images::Real = 1.0E-8,
                              tol_refine::Real = 1.0E-3,
                              tol_source::Real = 1.0E-6,
                              max_levels::Int  = 20,
                              max_newton::Int  = 20) where {T <: AbstractMatrix{<:Real}, F, G}
   # Check inputs
   if !(axes(θx) == axes(θy) == axes(βx) == axes(βy))
      throw(ArgumentError("θx, θy, βx, and βy must have the same axes."))
   end

   if size(θx, 1) < 2 || size(θx, 2) < 2
      throw(ArgumentError("The grid must have at least 2 × 2 points."))
   end

   if !(tol_images > 0) || !(tol_refine > 0) || !(tol_source > 0) || max_levels < 0 || max_newton < 0
      throw(ArgumentError("tol_images, tol_refine, and tol_source must be positive; max_levels and max_newton must be non-negative."))
   end

   # Create a point-type
   βs::_Pt = (Float64(β[1]), Float64(β[2]))
   tol_θ_f = Float64(tol_images)
   tol_β_f = Float64(tol_source)
   pad_f   = Float64(pad)
   cfac_f  = Float64(cfac)
   tol_r_f = max(Float64(tol_refine), tol_f)

   # -----------------------------------------------------------------------------------------------
   #  Level 0: two triangles per grid cell
   # -----------------------------------------------------------------------------------------------


end


end