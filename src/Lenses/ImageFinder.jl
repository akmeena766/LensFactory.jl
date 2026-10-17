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

   
end


end