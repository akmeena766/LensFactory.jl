module AstrometricOps

# --------------------------------------------------------------------------------------------------
# Gnomonic projection offsets (north up, east left)
# --------------------------------------------------------------------------------------------------
function gnomonic_offsets_arcsec(ra_ref::Float64, dec_ref::Float64, ra_cat::Vector{Float64}, dec_cat::Vector{Float64})
   deg2rad = π / 180.0
   rad2as  = 180.0 * 3600.0 / π

   # Degree --> Radian
   ra0  = ra_ref  * deg2rad
   dec0 = dec_ref * deg2rad
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

function gnomonic_offsets_arcsec(ra_ref::Float64, dec_ref::Float64, ra_cat::Float64, dec_cat::Float64)
   x, y = gnomonic_offsets_arcsec(ra_ref, dec_ref, [ra_cat], [dec_cat])
   return x[1], y[1]
end

end