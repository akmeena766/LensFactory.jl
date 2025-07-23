using Test
using LensFactory

# Tests for default flat cosmology module (AGAINST ASTROPY)
cosmo = Cosmology.init_cosmology()

com_dist = Cosmology.comoving_distance_radial(cosmo, 0.0, 1.0)
@test com_dist/Constants.DIST_MPC ≈ 3303.8288058874678

lum_dist = Cosmology.luminosity_distance(cosmo, 1)
@test lum_dist/Constants.DIST_MPC ≈ 6607.6576117749355

ang_dist = Cosmology.angular_diameter_distance(cosmo, 0.0, 1.0)
@test ang_dist/Constants.DIST_MPC ≈ 1651.9144029437339