module AstrometricOps

# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using LinearAlgebra


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------
using ..Constants

# --------------------------------------------------------------------------------------------------
# Various function to export
# --------------------------------------------------------------------------------------------------
export gnomonic_offsets_arcsec
export gnomonic_offsets_radec


# --------------------------------------------------------------------------------------------------
# Internal functions
# --------------------------------------------------------------------------------------------------
# Element type of the array outputs. The conversion constant is a Float64, so the result is at least
# Float64; abstract element types (e.g., a Vector{Real} read from YAML) also fall back to Float64.
@inline function _output_eltype(ref1::Real, ref2::Real, a::ROA, b::ROA)
   T = float(promote_type(typeof(ref1), typeof(ref2), eltype(a), eltype(b), Float64))
   return isconcretetype(T) ? T : Float64
end


# --------------------------------------------------------------------------------------------------
# Gnomonic projection offsets (north up, east left)
# --------------------------------------------------------------------------------------------------
"""
    gnomonic_offsets_arcsec(ra_ref::Float64, dec_ref::Float64, ra_cat::Float64, dec_cat::Float64)
"""
function gnomonic_offsets_arcsec(ra_ref::Real, dec_ref::Real, ra_cat::Real, dec_cat::Real)
   # Degree --> Radian
   ra0  = deg2rad(float(ra_ref))
   dec0 = deg2rad(float(dec_ref))
   ra   = deg2rad(float(ra_cat))
   dec  = deg2rad(float(dec_cat))

   # Denominator
   Δra  = ra - ra0
   cosc = sin(dec0) * sin(dec) + cos(dec0) * cos(dec) * cos(Δra)

   # In radians
   x = -(cos(dec) * sin(Δra)) / cosc
   y = +(cos(dec0) * sin(dec) - sin(dec0) * cos(dec) * cos(Δra)) / cosc

   return x / ANGLE_ARCSEC, y / ANGLE_ARCSEC
end


"""
    gnomonic_offsets_arcsec(ra_ref::Real, dec_ref::Real, ra_cat::ROA, dec_cat::ROA)
Project catalog positions onto the tangent plane at a reference position using the gnomonic
projection and return the angular offsets in arcseconds. The convention is north up, east left
(i.e., `x` increases towards west).

Scalar inputs `(ra_cat::Real, dec_cat::Real)` are also accepted, returning scalar offsets. Arrays
of any shape (vectors, matrices) are accepted and the output has the same shape as the input.

# Arguments
- `ra_ref`: Reference right ascension (in degrees)
- `dec_ref`: Reference declination (in degrees)
- `ra_cat`: Catalog right ascensions (in degrees)
- `dec_cat`: Catalog declinations (in degrees)

# Returns
- `x`: Offsets along the x-axis (in arcseconds)
- `y`: Offsets along the y-axis (in arcseconds)
"""
function gnomonic_offsets_arcsec(ra_ref::Real, dec_ref::Real, ra_cat::ROA, dec_cat::ROA)
   if size(ra_cat) != size(dec_cat)
      throw(ArgumentError("(RA, Dec) catalog must have the same size; got $(size(ra_cat)) and $(size(dec_cat))."))
   end
   
   # Output element type follows the inputs (Float64 unless they are of higher precision)
   T = _output_eltype(ra_ref, dec_ref, ra_cat, dec_cat)
   x = Array{T}(undef, size(ra_cat))
   y = Array{T}(undef, size(ra_cat))
   
   ax1, ax2 = axes(ra_cat, 1), axes(ra_cat, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         x[i, j], y[i, j] = gnomonic_offsets_arcsec(ra_ref, dec_ref, ra_cat[i, j], dec_cat[i, j])
      end
   end
   return x, y
end


"""
    gnomonic_offsets_radec(ra_ref::Real, dec_ref::Real, x_as::Real, y_as::Real)
"""
function gnomonic_offsets_radec(ra_ref::Real, dec_ref::Real, x_as::Real, y_as::Real)
   # Reference in Degrees --> Radians
   ra0  = deg2rad(float(ra_ref))
   dec0 = deg2rad(float(dec_ref))
 
   # Offsets in Arcsec --> Radians
   x = float(x_as) * ANGLE_ARCSEC
   y = float(y_as) * ANGLE_ARCSEC
 
   # Tangent-plane radius. In the gnomonic projection this equals tan(c), where c is the true
   # angular distance from the reference, so c = atan(rho) (NOT c = rho).
   rho = sqrt(x^2 + y^2)
 
   # At the exact centre the expressions below are 0/0; the answer is the reference itself
   if rho == 0.0
      return float(ra_ref), float(dec_ref)
   end
 
   c    = atan(rho)
   srho = sin(c)
   crho = cos(c)
 
   # Calculate Dec (clamp avoids domain errors in asin from float precision)
   dec = asin(clamp(crho * sin(dec0) + (y * srho * cos(dec0) / rho), -1.0, 1.0))
 
   # Calculate RA. The projection uses a north-up, east-left convention (x increases towards
   # west), so the east-positive standard coordinate is -x. atan2 handles the quadrant.
   ra = ra0 + atan(-x * srho, rho * cos(dec0) * crho - y * sin(dec0) * srho)
 
   return rad2deg(ra), rad2deg(dec)
end

"""
    gnomonic_offsets_radec(ra_ref::Real, dec_ref::Real, x_as::ROA, y_as::ROA)
Inverse gnomonic projection: convert tangent-plane offsets (in arcseconds) around a reference
position back to (RA, Dec) sky coordinates. Follows the same north up, east left convention as
[`gnomonic_offsets_arcsec`](@ref).

Scalar inputs `(x_as::Real, y_as::Real)` are also accepted, returning scalar coordinates. Arrays
(vectors or matrices) are accepted and the output has the same shape as the input.

# Arguments
- `ra_ref` : Reference right ascension (in degrees)
- `dec_ref`: Reference declination (in degrees)
- `x_as`   : Offsets along the x-axis (in arcseconds)
- `y_as`   : Offsets along the y-axis (in arcseconds)

# Returns
- `ra` : Right ascensions (in degrees)
- `dec`: Declinations (in degrees)
"""
function gnomonic_offsets_radec(ra_ref::Real, dec_ref::Real, x_as::ROA, y_as::ROA)
   if size(x_as) != size(y_as)
      throw(ArgumentError("Input catalog must have the same size; got $(size(x_as)) and $(size(y_as))."))
   end

   # Output element type follows the inputs (Float64 unless they are of higher precision)
   T   = _output_eltype(ra_ref, dec_ref, x_as, y_as)
   ra  = Array{T}(undef, size(x_as))
   dec = Array{T}(undef, size(x_as))

   ax1, ax2 = axes(x_as, 1), axes(x_as, 2)
   @inbounds for j in ax2
      @inbounds for i in ax1
         ra[i, j], dec[i, j] = gnomonic_offsets_radec(ra_ref, dec_ref, x_as[i, j], y_as[i, j])
      end
   end
   return ra, dec
end

end