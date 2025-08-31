@testset "flat LCDM" begin
   # Tests for default flat LCDM cosmology module (AGAINST ASTROPY)
   cosmo = Cosmology.init_cosmology()
   z = 1.0

   @test Cosmology.scale_factor(z) ≈ 0.5
   @test Cosmology.Ez(cosmo, z) ≈ 1.7606816861659007
   @test Cosmology.hubble_parameter(cosmo, z) ≈ 123.24771803161305
   @test Cosmology.hubble_time(cosmo.H0) ≈ 13.968460309725561

   @test Cosmology.age(cosmo, z) ≈ 5.751646943448283
   @test Cosmology.lookback_time(cosmo, z) ≈ 7.715337003613594

   @test Cosmology.rho_cz(cosmo, z) ≈ 2.85320091612148e-26
   @test Cosmology.Omega_mz(cosmo, z) ≈ 0.7741935483870968
   @test Cosmology.Omega_rz(cosmo, z) ≈ 0.0
   @test Cosmology.Omega_kz(cosmo, z) ≈ 0.0
   @test Cosmology.Omega_wz(cosmo, z) ≈ 0.2258064516129032

   @test Cosmology.comoving_distance_radial(cosmo, 0.0, z)/Constants.DIST_MPC ≈ 3303.8288058874678
   @test Cosmology.luminosity_distance(cosmo, z)/Constants.DIST_MPC ≈ 6607.6576117749355
   @test Cosmology.angular_diameter_distance(cosmo, 0.0, z)/Constants.DIST_MPC ≈ 1651.9144029437339

   @test Cosmology.comoving_volume_element(cosmo, z) ≈ 26.550755712253192
   @test Cosmology.comoving_volume(cosmo, z) ≈ 151.05712532061932
end

@testset "open LCDM" begin
   # Tests for default flat LCDM cosmology module (AGAINST ASTROPY)
   cosmo = Cosmology.init_cosmology(Omega_m0=0.26)
   z = 1.0

   @test Cosmology.Ez(cosmo, z) ≈ 1.7146428199482249
   @test Cosmology.hubble_parameter(cosmo, z) ≈ 120.02499739637574

   @test Cosmology.age(cosmo, z) ≈ 6.022600655462329
   @test Cosmology.lookback_time(cosmo, z) ≈ 7.796137793333031

   @test Cosmology.rho_cz(cosmo, z) ≈ 2.7059389333539216e-26

   @test Cosmology.comoving_distance_radial(cosmo, 0.0, z)/Constants.DIST_MPC ≈ 3343.5622067797044
   @test Cosmology.comoving_distance_transverse(cosmo, 0.0, 1)/Constants.DIST_MPC ≈ 3357.164778591662
   @test Cosmology.luminosity_distance(cosmo, z)/Constants.DIST_MPC ≈ 6714.329557183324
   @test Cosmology.angular_diameter_distance(cosmo, 0.0, z)/Constants.DIST_MPC ≈ 1678.582389295831

   @test Cosmology.comoving_volume_element(cosmo, z) ≈ 28.15103157576382
   @test Cosmology.comoving_volume(cosmo, z) ≈ 157.3382098743564
end

@testset "closed CDM" begin
   # Tests for default flat LCDM cosmology module (AGAINST ASTROPY)
   cosmo = Cosmology.init_cosmology(Omega_m0=0.34)
   z = 1.0

   @test Cosmology.Ez(cosmo, z) ≈ 1.8055470085267789
   @test Cosmology.hubble_parameter(cosmo, z) ≈ 126.38829059687453

   @test Cosmology.age(cosmo, z) ≈ 5.516141610883366
   @test Cosmology.lookback_time(cosmo, z) ≈ 7.638261424137102

   @test Cosmology.rho_cz(cosmo, z) ≈ 3.0004628988890413e-26

   @test Cosmology.comoving_distance_radial(cosmo, 0.0, z)/Constants.DIST_MPC ≈ 3266.04304016069
   @test Cosmology.comoving_distance_transverse(cosmo, 0.0, 1)/Constants.DIST_MPC ≈ 3253.394978599389
   @test Cosmology.luminosity_distance(cosmo, z)/Constants.DIST_MPC ≈ 6506.789957198778
   @test Cosmology.angular_diameter_distance(cosmo, 0.0, z)/Constants.DIST_MPC ≈ 1626.6974892996946

   @test Cosmology.comoving_volume_element(cosmo, z) ≈ 25.10657361038687
   @test Cosmology.comoving_volume(cosmo, z) ≈ 145.2558165125011
end