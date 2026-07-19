module AstrometricOps

# --------------------------------------------------------------------------------------------------
# Julia inbuilt functions to import
# --------------------------------------------------------------------------------------------------
using LinearAlgebra


# --------------------------------------------------------------------------------------------------
# LensFactory modules to use
# --------------------------------------------------------------------------------------------------


# --------------------------------------------------------------------------------------------------
# Various function to export
# --------------------------------------------------------------------------------------------------
export gnomonic_offsets_arcsec
export gnomonic_offsets_radec


# --------------------------------------------------------------------------------------------------
# Gnomonic projection offsets (north up, east left)
# --------------------------------------------------------------------------------------------------
"""
    gnomonic_offsets_arcsec(ra_ref::Float64, dec_ref::Float64, ra_cat::Float64, dec_cat::Float64)
"""
function gnomonic_offsets_arcsec(ra_ref::Float64, dec_ref::Float64, ra_cat::Float64, dec_cat::Float64)
   x, y = gnomonic_offsets_arcsec(ra_ref, dec_ref, [ra_cat], [dec_cat])
   return x[1], y[1]
end


"""
    gnomonic_offsets_arcsec(ra_ref::Float64, dec_ref::Float64, ra_cat::Vector{Float64}, dec_cat::Vector{Float64})
Project catalog positions onto the tangent plane at a reference position using the gnomonic
projection and return the angular offsets in arcseconds. The convention is north up, east left
(i.e., `x` increases towards west).

Scalar inputs `(ra_cat::Float64, dec_cat::Float64)` are also accepted, returning scalar offsets.

# Arguments
- `ra_ref::Float64`: Reference right ascension (in degrees)
- `dec_ref::Float64`: Reference declination (in degrees)
- `ra_cat::Vector{Float64}`: Catalog right ascensions (in degrees)
- `dec_cat::Vector{Float64}`: Catalog declinations (in degrees)

# Returns
- `x::Vector{Float64}`: Offsets along the x-axis (in arcseconds)
- `y::Vector{Float64}`: Offsets along the y-axis (in arcseconds)
"""
function gnomonic_offsets_arcsec(ra_ref::Float64, dec_ref::Float64, ra_cat::Vector{Float64}, dec_cat::Vector{Float64})
   deg2rad = π / 180.0
   rad2as  = 180.0 * 3600.0 / π

   # Reference in Degree --> Radian
   ra0  = ra_ref  * deg2rad
   dec0 = dec_ref * deg2rad
   
   # Catalog in Degree --> Radian
   ra   = ra_cat  .* deg2rad
   dec  = dec_cat .* deg2rad

   # Denomintor
   Δra  = ra .- ra0
   cosc = @. sin(dec0)*sin(dec) + cos(dec0)*cos(dec)*cos(Δra)

   # In radians
   x = @. -(cos(dec) * sin(Δra)) / cosc      # radians
   y = @. +(cos(dec0)*sin(dec) - sin(dec0)*cos(dec)*cos(Δra)) / cosc

   return x .* rad2as, y .* rad2as
end


"""
    gnomonic_offsets_radec(ra_ref::Float64, dec_ref::Float64, x_as::Float64, y_as::Float64)
"""
function gnomonic_offsets_radec(ra_ref::Float64, dec_ref::Float64, x_as::Float64, y_as::Float64)
   x, y = gnomonic_offsets_radec(ra_ref, dec_ref, [x_as], [y_as])
   return x[1], y[1]
end

"""
    gnomonic_offsets_radec(ra_ref::Float64, dec_ref::Float64, x_as::Vector{Float64}, y_as::Vector{Float64})
Inverse gnomonic projection: convert tangent-plane offsets (in arcseconds) around a reference
position back to (RA, Dec) sky coordinates. Follows the same north up, east left convention as
[`gnomonic_offsets_arcsec`](@ref).

Scalar inputs `(x_as::Float64, y_as::Float64)` are also accepted, returning scalar coordinates.

# Arguments
- `ra_ref::Float64`: Reference right ascension (in degrees)
- `dec_ref::Float64`: Reference declination (in degrees)
- `x_as::Vector{Float64}`: Offsets along the x-axis (in arcseconds)
- `y_as::Vector{Float64}`: Offsets along the y-axis (in arcseconds)

# Returns
- `ra::Vector{Float64}`: Right ascensions (in degrees)
- `dec::Vector{Float64}`: Declinations (in degrees)
"""
function gnomonic_offsets_radec(ra_ref::Float64, dec_ref::Float64, x_as::Vector{Float64}, y_as::Vector{Float64})
   deg2rad = π / 180.0
   as2rad  = π / (180.0 * 3600.0)

   # Reference in Degrees --> Radians
   ra0  = ra_ref  * deg2rad
   dec0 = dec_ref * deg2rad
   
   # Offsets in Arcsec --> Radians
   x = x_as .* as2rad
   y = y_as .* as2rad

   # Angular distance rho
   rho = @. sqrt(x^2 + y^2)
   
   # Avoid division by zero at the exact center
   # For very small rho, sin(rho)/rho is approximately 1
   srho = @. sin(rho)
   crho = @. cos(rho)

   # Calculate Dec
   # Note: we use clamp to avoid domain errors in asin from float precision
   term_dec = @. crho * sin(dec0) + (y * srho * cos(dec0) / rho)
   replace!(term_dec, NaN => sin(dec0)) # Handle rho=0 case
   dec = @. asin(clamp(term_dec, -1.0, 1.0))

   # Calculate RA
   num_ra = @. x * srho
   den_ra = @. rho * cos(dec0) * crho - y * sin(dec0) * srho
   
   # atan2 handles the quadrant correctly
   # If rho is 0, Δra should be 0
   ra = @. ra0 + atan(num_ra, den_ra)

   return ra ./ deg2rad, dec ./ deg2rad
end

end