using Test
using LensFactory

# Tests for default flat cosmology module (AGAINST ASTROPY)
cosmo = Cosmology.init_cosmology()

a = Cosmology.scale_factor(1.0)
@test a == 0.5

ez = Cosmology.Ez(cosmo, 1.0)
@test ez ≈ 1.7606816861659007

Hz = Cosmology.hubble_parameter(cosmo, 1.0)
@test Hz ≈ 123.24771803161305

Ht = Cosmology.hubble_time(cosmo.H0)
@test Ht ≈ 13.968460309725561

age = Cosmology.age(cosmo, 1.0)
@test age ≈ 5.751646943448283

lb_time = Cosmology.lookback_time(cosmo, 1.0)
@test lb_time ≈ 7.715337003613594

com_dist = Cosmology.comoving_distance_radial(cosmo, 0.0, 1.0)
@test com_dist/Constants.DIST_MPC ≈ 3303.8288058874678

lum_dist = Cosmology.luminosity_distance(cosmo, 1)
@test lum_dist/Constants.DIST_MPC ≈ 6607.6576117749355

ang_dist = Cosmology.angular_diameter_distance(cosmo, 0.0, 1.0)
@test ang_dist/Constants.DIST_MPC ≈ 1651.9144029437339

com_vol_element = Cosmology.comoving_volume_element(cosmo, 1.0)
@test com_vol_element ≈ 26.550755712253192

com_vol = Cosmology.comoving_volume(cosmo, 1.0)
@test com_vol ≈ 151.05712532061932