using Test
using LensFactory

# Tests for cosmology module (AGAINST ASTROPY)
cosmo = Cosmology.init_cosmology()

# Default flat cosmology
ang_dist = Cosmology.angular_diameter_distance(cosmo, 0.0, 1.0)
@test ang_dist/Constants.DIST_MPC ≈ 1651.9144029437339