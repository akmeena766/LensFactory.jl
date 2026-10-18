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

# Does a critical curve (det A = 0) possibly pass through triangle `t`? If so, the triangle may
# contain two images, and it is not handed over to Newton-Raphson before reaching `tol`.
@inline function _crosses_critical(lens_jac::G, t::NTuple{3, _Pt}) where G
   d1 = _det(lens_jac, t[1])
   d2 = _det(lens_jac, t[2])
   d3 = _det(lens_jac, t[3])
   if !(isfinite(d1) && isfinite(d2) && isfinite(d3))
      return true
   end
   return min(d1, d2, d3) <= 0.0 <= max(d1, d2, d3)
end

# Deviation of the mapped midpoint `M` of an edge from the midpoint of the mapped edge `AB`, i.e., a
# measure of the non-linearity of the lens mapping along that edge (zero if it cannot be computed)
@inline function _deviation(A::_Pt, B::_Pt, M::_Pt)
   return _finite_or_zero(hypot(M[1] - 0.5 * (A[1] + B[1]), M[2] - 0.5 * (A[2] + B[2])))
end

# Distance between point `p` and segment `ab`
@inline function _segment_distance(p::_Pt, a::_Pt, b::_Pt)
   dx = b[1] - a[1]
   dy = b[2] - a[2]
   L2 = dx^2 + dy^2

   # Degenerate segment: a single point
   if L2 == 0.0
      return _dist(p, a)
   end
   
   # Position of the projection of p along the segment (0 at a, 1 at b), clamped to the segment
   s = ((p[1] - a[1]) * dx + (p[2] - a[2]) * dy) / L2
   s = clamp(s, 0.0, 1.0)
   return hypot(p[1] - (a[1] + s * dx), p[2] - (a[2] + s * dy))
end

# Estimate the non-linearity of the lens mapping at every grid node from second differences
# (`c`, along x and y) and in every cell from the mixed difference (`m`); non-finite values → 0
function _grid_curvature(βx::AbstractMatrix{<:Real}, βy::AbstractMatrix{<:Real})
   ax1, ax2 = axes(βx, 1), axes(βx, 2)
   
   c = zeros(Float64, size(βx))
   m = zeros(Float64, size(βx))
   
   i0, i1 = first(ax1), last(ax1)
   j0, j1 = first(ax2), last(ax2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ci = i - i0 + 1
         cj = j - j0 + 1
         if i0 < i < i1
            v = hypot(βx[i+1, j] - 2βx[i, j] + βx[i-1, j], βy[i+1, j] - 2βy[i, j] + βy[i-1, j])
            c[ci, cj] = max(c[ci, cj], _finite_or_zero(Float64(v)))
         end
         if j0 < j < j1
            v = hypot(βx[i, j+1] - 2βx[i, j] + βx[i, j-1], βy[i, j+1] - 2βy[i, j] + βy[i, j-1])
            c[ci, cj] = max(c[ci, cj], _finite_or_zero(Float64(v)))
         end
         if i < i1 && j < j1
            v = hypot(βx[i+1, j+1] - βx[i+1, j] - βx[i, j+1] + βx[i, j],
                      βy[i+1, j+1] - βy[i+1, j] - βy[i, j+1] + βy[i, j])
            m[ci, cj] = _finite_or_zero(Float64(v))
         end
      end
   end
   return c, m
end

# Test whether the source position `β` may have a pre-image inside triangle `t`. This is the case if
# `β` lies inside the source-plane triangle `t.β`, where the barycentric coordinates are allowed to 
# be as small as `-pad`, and `β` may also lie up to `t.δ` outside of it (to account for the 
# non-linearity of the mapping). Triangles that collapse to a line (e.g., when straddling a critical 
# curve) are tested against their padded bounding box. Triangles with a non-finite vertex (i.e., 
# touching a singular lens centre) are always kept, so that images arbitrarily close to the 
# singularity can still be found, unless all their vertices are non-finite (e.g., inside a region 
# where the mapping is undefined).
@inline function _keep(t::_Tri, β::_Pt, pad::Float64)
   # If any vertex is non-finite, keep the triangle unless all are non-finite.
   # This allows finding images near singularities.
   if !_isfinite(t.β)
      return _isfinite(t.β[1]) || _isfinite(t.β[2]) || _isfinite(t.β[3])
   end

   # Extract triangle details
   (x1, y1), (x2, y2), (x3, y3) = t.β
   βx, βy = β
   δ = t.δ

   # Padded bounding box test (cheap rejection)
   xmin, xmax = min(x1, x2, x3), max(x1, x2, x3)
   ymin, ymax = min(y1, y2, y3), max(y1, y2, y3)
   wx, wy = pad * (xmax - xmin) + δ, pad * (ymax - ymin) + δ
   if βx < xmin - wx || βx > xmax + wx || βy < ymin - wy || βy > ymax + wy
      return false
   end

   # Degenerate (zero area) triangle: the bounding box test is all we can do
   d = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
   if abs(d) <= 1e-14 * max((xmax - xmin)^2 + (ymax - ymin)^2, floatmin(Float64))
      return true
   end

   # Barycentric coordinates
   l1 = ((y2 - y3) * (βx - x3) + (x3 - x2) * (βy - y3)) / d
   l2 = ((y3 - y1) * (βx - x3) + (x1 - x3) * (βy - y3)) / d
   l3 = 1.0 - l1 - l2
   if (l1 >= -pad) && (l2 >= -pad) && (l3 >= -pad)
      return true
   end

   # Outside: accept if within the non-linearity error δ of the triangle
   if δ > 0.0
      A, B, C = t.β
      dist_edge = min(_segment_distance(β, A, B), _segment_distance(β, B, C), _segment_distance(β, C, A))
      if dist_edge <= δ
         return true
      end
   end
   return false
end

# Split a triangle into four congruent children using its edge midpoints. The error estimate of the
# children is updated from the measured non-linearity along the parent edges.
@inline function _subdivide(t::_Tri, lens_map::F, cfac::Float64) where F
   a, b, c = t.θ
   A, B, C = t.β
   ab, bc, ca = _mid(a, b), _mid(b, c), _mid(c, a)
   AB, BC, CA = _map(lens_map, ab), _map(lens_map, bc), _map(lens_map, ca)

   # Midpoint deviation e ≈ |∂²β| L² / 8 for the parent ⇒ |∂²β| (L/2)² ≈ 2 e for the children
   e = max(_deviation(A, B, AB), _deviation(B, C, BC), _deviation(C, A, CA))
   δ = max(0.25 * t.δ, 2.0 * cfac * e)
   l = t.level + 1

   return (_Tri((a, ab, ca), (A, AB, CA), δ, l),
           _Tri((ab, b, bc), (AB, B, BC), δ, l),
           _Tri((ca, bc, c), (CA, BC, C), δ, l),
           _Tri((ab, bc, ca), (AB, BC, CA), δ, l))
end

# Newton-Raphson for β(θ) = β. The step is damped (halved) until the source-plane residual
# decreases, and the iterate is never allowed to move further than `R` from the starting point.
# Returns the best point and its residual.
function _newton(lens_map::F, lens_jac::G, θ0::_Pt, β::_Pt, R::Float64, max_iter::Int) where {F, G}
   # Initial guess
   x, y = θ0

   # Map the inital guess to source plane
   bx, by = _map(lens_map, (x, y))

   # Resisual with respect to actual source position (β)
   fx, fy = bx - β[1], by - β[2]
   r = hypot(fx, fy)
   r_floor = 4.0 * _residual_floor(θ0, β)

   for _ in 1:max_iter
      # Get out if the residual is small enough
      (isfinite(r) && r > r_floor) || break

      # Calculate jacobian (∂β/∂θ)
      a11, a12, a21, a22 = lens_jac(x, y)
      det = a11 * a22 - a12 * a21
      
      # Get out if jacobian is (not finite) or zero
      (isfinite(det) && det != 0.0) || break

      # Full Newton step
      dx = ( a22 * fx - a12 * fy) / det
      dy = (-a21 * fx + a11 * fy) / det

      # Damped line search
      s = 1.0
      improved = false
      for _ in 1:30
         xn = x - s * dx
         yn = y - s * dy
         if hypot(xn - θ0[1], yn - θ0[2]) <= R
            bxn, byn = _map(lens_map, (xn, yn))
            fxn, fyn = bxn - β[1], byn - β[2]
            rn = hypot(fxn, fyn)
            if isfinite(rn) && rn < r
               x, y, fx, fy, r = xn, yn, fxn, fyn, rn
               improved = true
               break
            end
         end
         s *= 0.5
      end
      improved || break
   end
   return (x, y), r
end


"""
    find_images_adaptive(lens_map::F, lens_jac::G,
                         θx::T, θy::T, βx::T, βy::T, β::NTuple{2, Real};
                         tol_images::Real    = 1.0E-8,
                         tol_refine::Real    = 1.0E-3,
                         tol_source::Real    = 1.0E-6,
                         max_levels::Int     = 20,
                         max_newton::Int     = 20,
                         pad::Real           = 0.1,
                         cfac::Real          = 1.0,
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
                              tol_images::Real    = 1.0E-8,
                              tol_refine::Real    = 1.0E-3,
                              tol_source::Real    = 1.0E-6,
                              max_levels::Int     = 20,
                              max_newton::Int     = 20,
                              pad::Real           = 0.1,
                              cfac::Real          = 1.0,
                              max_candidates::Int = 100_000) where {T <: AbstractMatrix{<:Real}, F, G}
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
   c, m = _grid_curvature(βx, βy)
   active = _Tri[]

   ax1, ax2 = axes(θx, 1), axes(θx, 2)
   i0 = first(ax1) 
   j0 = first(ax2)

   @inbounds for j in j0:last(ax2)-1
      @inbounds for i in i0:last(ax1)-1
         ci = i - i0 + 1
         cj = j - j0 + 1
         δ  = cfac_f * max(c[ci, cj], c[ci+1, cj], c[ci, cj+1], c[ci+1, cj+1], m[ci, cj])

         v00 = (Float64(θx[i,   j  ]), Float64(θy[i,   j  ]))
         v10 = (Float64(θx[i+1, j  ]), Float64(θy[i+1, j  ]))
         v01 = (Float64(θx[i,   j+1]), Float64(θy[i,   j+1]))
         v11 = (Float64(θx[i+1, j+1]), Float64(θy[i+1, j+1]))
         b00 = (Float64(βx[i,   j  ]), Float64(βy[i,   j  ]))
         b10 = (Float64(βx[i+1, j  ]), Float64(βy[i+1, j  ]))
         b01 = (Float64(βx[i,   j+1]), Float64(βy[i,   j+1]))
         b11 = (Float64(βx[i+1, j+1]), Float64(βy[i+1, j+1]))
         
         t1 = _Tri((v00, v10, v11), (b00, b10, b11), δ, 0)
         t2 = _Tri((v00, v11, v01), (b00, b11, b01), δ, 0)
         _keep(t1, βs, pad_f) && push!(active, t1)
         _keep(t2, βs, pad_f) && push!(active, t2)
      end
   end

   # -----------------------------------------------------------------------------------------------
   # Level 1: Refinement + Newton-Raphson
   # -----------------------------------------------------------------------------------------------
   sols = Tuple{_Pt, Float64}[]
   while !isempty(active)
      next = _Tri[]
      for t in active
         # Tringle size
         edge = _edge(t.θ)

         # Check if it is smaller than tolerence in image plane or we exceeded the maximum 
         # refinement iterations
         is_last = (edge <= tol_θ_f) || (t.level >= max_levels)

         # Newton-Raphson once the triangle is small enough, and always on the last level
         if _isfinite(t.β) && (is_last || (edge <= tol_r_f && max_newton > 0))
            if is_last || !_crosses_critical(lens_jac, t.θ)
               # Newton-Raphson from the centroid, accepted only if it stays in this triangle
               θn, r = _newton(lens_map, lens_jac, _centroid(t.θ), βs, 2.0 * edge, max_newton)
               in_t = isfinite(r) && _inside(θn, t.θ, pad_f)

               # Residual below which Newton-Raphson is considered converged (never below rounding)
               r_conv = min(tol_β_f, max(1e-3 * tol_θ_f, 40.0 * _residual_floor(θn, βs)))
               if in_t && r <= r_conv
                  push!(sols, (θn, r))
                  continue
               end

               # Smallest triangle: fall back to the best point available
               if is_last
                  if in_t && r <= β_tol_f
                     push!(sols, (θn, r))
                  else
                     θc = _centroid(t.θ)
                     bc = _map(lens_map, θc)
                     rc = hypot(bc[1] - βs[1], bc[2] - βs[2])
                     (isfinite(rc) && rc <= β_tol_f) && push!(sols, (θc, rc))
                  end
                  continue
               end
            end
         end
         is_last && continue

         # Refine
         for child in _subdivide(t, lens_map, cfac_f)
            _keep(child, βs, pad_f) && push!(next, child)
         end
      end
      if length(next) > max_candidates
         error("Adaptive image finder: more than $max_candidates candidate triangles. This " *
               "usually means that the solutions are not isolated points (e.g., an Einstein " *
               "ring for a source exactly behind an axisymmetric lens). Otherwise, try a larger " *
               "`tol_images` or a larger `max_candidates`.")
      end
      active = next
   end

   # -----------------------------------------------------------------------------------------------
   # Merge duplicates (keep the one with the smallest residual)
   # -----------------------------------------------------------------------------------------------
   tol_merge = 10.0 * tol_f
   sort!(sols, by = last)
   images = _Pt[]
   for (p, _) in sols
      if all(q -> _dist(p, q) > tol_merge, images)
         push!(images, p)
      end
   end
   return images
end


end