using Test
using LensFactory

# Tests for cosmology module
cosmo = Cosmology.init_cosmology()

# Default flat cosmology
ang_dist = Cosmology.angular_diameter_distance(cosmo, 0.0, 1.0)
round(ang_dist/Constants.DIST_MPC, sigdigits=10) == round(1651.9144029437339, sigdigits=10)
