module PolygonOps


# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using StatsBase
using LinearAlgebra


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# Various function to export
# --------------------------------------------------------------------------------------------------
export bilinear_interpolation
export shoelace
export hao_sun
export winding_number
export fit_ellipse


"""
    bilinear_interpolation(x::Real, y::Real, df::AbstractMatrix{<:Real})
Bilinear interpolation at (x, y) given in pixel coordinates.

# Arguments
- `x::Float64`: x pixel coordinate
- `y::Float64`: y pixel coordinate
- `df::Matrix{<:Float64}`: Data matrix

# Returns
- `Float64`: Interpolated value at (x, y) position
"""
function bilinear_interpolation(x::Real, y::Real, df::AbstractMatrix{<:Real})
   # Data matrix dimesions
   nx, ny = size(df)
   
   # Check if the pixel is inside the grid range
   if x < 1 || x ≥ nx || y < 1 || y ≥ ny
      throw(ArgumentError("Point ($x, $y) is outside the valid interpolation range [1, $(nx)), [1, $(ny))"))
   end
   
   # Lower-left pixel position
   px = floor(Int, x)
   py = floor(Int, y)

   # Fractional offsets
   dx = x - px
   dy = y - py

   # Funtion values at the vertices
   f00 = df[px + 0, py + 0]
   f01 = df[px + 0, py + 1]
   f10 = df[px + 1, py + 0]
   f11 = df[px + 1, py + 1]

   return f00 * (1 - dx) * (1 - dy) + f01 * (1 - dx) * dy + f10 * dx * (1 - dy) + f11 * dx * dy
end


"""
    shoelace(polygon::Vector{<:Vector{<:Real}})
Calculate the area of a 2D polygon using the shoelace formula.

# Arguments
- `polygon::Vector{Vector{<:Real}}`: A vector of `[x, y]` coordinates representing the polygon vertices

# Returns
- The area of the polygon
"""
function shoelace(polygon::Vector{<:Vector{<:Real}})
   # Close the polygon if needed (copy to avoid mutating the caller's input)
   if polygon[1] != polygon[end]
      polygon = vcat(polygon, [polygon[1]])
   end

   # Calculate the area
   area = 0.0
   for i in 1:(length(polygon)-1)
      x1, y1 = polygon[i]
      x2, y2 = polygon[i+1]
      area += x1*y2 - x2*y1
   end
   return 0.5 * abs(area)
end


"""
    hao_sun(point, polygon)
Algorithm to determine if a point is in the given polygon. Taken from [sym10100477](@citet)

# Arguments
- `point::Vector{Float64}`: The (x, y) coordinate of the point to check
- `polygon::Vector{Vector{Float64}}`: A vector of [x, y] coordinates representing the polygon vertices

# Returns
- `Int`: 
   - +1: inside
   - 0: outside
   - -1: on the edge
"""
function hao_sun(point::Vector{Float64}, polygon::Vector{Vector{Float64}})
   k = 0

   xp = point[1]
   yp = point[2]

   @inbounds for i in 1:length(polygon)-1
      v1 = polygon[i][2] - yp
      v2 = polygon[i+1][2] - yp

      if v1 < 0 && v2 < 0 || v1 > 0 && v2 > 0
         continue
      end

      u1 = polygon[i][1]   - xp
      u2 = polygon[i+1][1] - xp

      f = (u1 * v2) - (u2 * v1)

      if v2 > 0 && v1 <= 0
         if f > 0
            k += 1
         elseif f == 0
            return -1
         end
      elseif v1 > 0 && v2 <= 0
         if f < 0
            k += 1
         elseif f == 0
            return -1
         end
      elseif v2 == 0 && v1 < 0
         if f == 0
            return -1
         end   
      elseif v1 == 0 && v2 < 0
         if f == 0
            return -1
         end 
      elseif v1 == 0 && v2 == 0
         if u2 <= 0 && u1 >= 0
            return -1
         elseif u1 <= 0 && u2 >= 0
            return -1
         end
      end
   end
   
   if k % 2 == 0
      return 0
   end
   return 1
end


"""
    winding_number(point::Vector{<:Real}, polygon::Vector{<:Vector{<:Real}})
Calculate the winding number of a closed `polygon` about the given `point`, i.e. the (signed)
number of times the polygon wraps around the point. A ray is shot from the point towards ``+x``
and the crossings of the polygon are counted `+1` (polygon crossing the ray upwards) or `-1`
(crossing it downwards).

The winding number is undefined for a point lying exactly on an edge, in which case the returned
value depends on which side of the edge the floating point arithmetic lands.

# Arguments
- `point`  : The (x, y) coordinate of the point to check
- `polygon`: A vector of [x, y] coordinates representing the polygon vertices
 
# Returns
- `Int`: Winding number about the given point (`0` if the point lies outside the polygon)
"""
function winding_number(point::Vector{<:Real}, polygon::Vector{<:Vector{<:Real}})
   xp = point[1]
   yp = point[2]

   # The polygon is closed by the loop below. A repeated last vertex needs no special treatment,
   # since the resulting zero-length edge cannot cross the ray and contributes nothing.
   n = length(polygon)
   if n < 3
      return 0
   end

   wn = 0
   @inbounds for i in 1:n
      # Index of the next vertex, wrapping around at the last one
      if i == n
         j = 1
      else
         j = i + 1
      end

      x1, y1 = polygon[i][1], polygon[i][2]
      x2, y2 = polygon[j][1], polygon[j][2]

      # Positive if the point lies to the left of the directed edge (i --> j)
      side = (x2 - x1) * (yp - y1) - (xp - x1) * (y2 - y1)

      if y1 <= yp
         if y2 > yp && side > 0
            wn = wn + 1
         end
      else
         if y2 <= yp && side < 0
            wn = wn - 1
         end
      end
   end
   return wn
end


"""
    fit_ellipse(x::Vector{Float64}, y::Vector{Float64})
Fits an ellipse to a set of points (x, y) using the least squares method with a quadratic constraint.

# Arguments
- `x::Vector{Float64}`: The x-coordinates of the points
- `y::Vector{Float64}`: The y-coordinates of the points

# Returns
- `cx`: The x-coordinate of the center of the ellipse
- `cy`: The y-coordinate of the center of the ellipse
- `rx`: The semi-major axis length of the ellipse
- `ry`: The semi-minor axis length of the ellipse
- `θ`: The rotation angle of the ellipse (in degrees)
"""
function fit_ellipse(x::Vector{Float64}, y::Vector{Float64})
   # --- Normalize ---
   mx = mean(x)
   my = mean(y)

   sx = (maximum(x) - minimum(x)) / 2
   sy = (maximum(y) - minimum(y)) / 2

   xn = (x .- mx) ./ sx
   yn = (y .- my) ./ sy

   # --- Fit on normalized data ---
   D = hcat(xn.^2, xn.*yn, yn.^2, xn, yn, ones(length(xn)))
   S = D' * D

   # Constraint matrix for b² - 4ac = -1
   C = zeros(6, 6)
   C[1,3] = 2.0
	C[3,1] = 2.0
	C[2,2] = -1.0

   # Eigen Decomposition
   eigenvalues, eigenvectors = eigen(S, C)

   # We want the only positive eigenvalue for this specific C matrix
   valid_idx = findall(λ -> isreal(λ) && isfinite(real(λ)) && real(λ) > 0, eigenvalues)
   isempty(valid_idx) && error("No valid ellipse solution found")

   # Get the coefficients for the best fit
   idx = valid_idx[argmin(real(eigenvalues[valid_idx]))]
   a, b, c, d, e, f = real(eigenvectors[:, idx])

   # --- Denormalize ---
   a2 = a / sx^2
   b2 = b / (sx * sy)
   c2 = c / sy^2
   d2 = d / sx - 2a * mx / sx^2 - b * my / (sx * sy)
   e2 = e / sy - 2c * my / sy^2 - b * mx / (sx * sy)
   f2 = (a * mx^2 / sx^2 + b * mx * my / (sx * sy) + c * my^2 / sy^2 - d * mx / sx - e * my / sy + f)

   # --- Geometric parameters ---

   # Find Center via the intersection of partial derivatives
   cx, cy = [2.0*a2  b2; 
             b2      2.0*c2] \ [-d2, -e2]

   # Calculate the constant term relative to the center
   f_center = a2 * cx^2 + b2 * cx * cy + c2 * cy^2 + d2 * cx + e2 * cy + f2

   # Extract axes and rotation from the quadratic part
   eval2, evec2 = eigen([a2       b2/2.0; 
                         b2/2.0   c2])

   # Eigenvalues λ are related to axes lengths by r = sqrt(|-f_center / λ|)
   λ1, λ2 = eval2[1], eval2[2]
   r1 = sqrt(abs(-f_center / λ1))
   r2 = sqrt(abs(-f_center / λ2))

   if r1 >= r2
      rx, ry = r1, r2
      θ = atan(evec2[2,1], evec2[1,1])
   else
      rx, ry = r2, r1
      θ = atan(evec2[2,2], evec2[1,2])
   end
   θ = rad2deg(θ)

   return cx, cy, rx, ry, θ
end


end