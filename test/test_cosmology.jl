@testset "Cosmology" begin
   @testset "flat LCDM" begin
      # Tests for default flat (Ωk = 0) LCDM cosmology module (AGAINST ASTROPY)
      cosmo = Cosmology.init_cosmology()
      z = 1.0

      @test Cosmology.scale_factor(z) ≈ 0.5 atol=1e-15 rtol=1e-15
      @test Cosmology.Ez(cosmo, z) ≈ 1.7606816861659007 atol=1e-15 rtol=1e-15
      @test Cosmology.hubble_parameter(cosmo, z) ≈ 123.24771803161305 atol=1e-15 rtol=1e-15
      @test Cosmology.hubble_time(cosmo.H0) ≈ 13.968460309725561 atol=1e-15 rtol=1e-15

      @test Cosmology.age(cosmo, z) ≈ 5.751646943448283 atol=1e-8 rtol=1e-8
      @test Cosmology.lookback_time(cosmo, z) ≈ 7.715337003613594 atol=1e-8 rtol=1e-8

      @test Cosmology.rho_cz(cosmo, z) ≈ 2.85320091612148e-26 atol=1e-15 rtol=1e-15
      @test Cosmology.Omega_mz(cosmo, z) ≈ 0.7741935483870968 atol=1e-15 rtol=1e-15
      @test Cosmology.Omega_rz(cosmo, z) ≈ 0.0 atol=1e-15 rtol=1e-15
      @test Cosmology.Omega_kz(cosmo, z) ≈ 0.0 atol=1e-15 rtol=1e-15
      @test Cosmology.Omega_wz(cosmo, z) ≈ 0.2258064516129032 atol=1e-15 rtol=1e-15

      @test Cosmology.comoving_distance_radial(cosmo, 0.0, z)/Constants.DIST_MPC ≈ 3303.8288058874678 atol=1e-8 rtol=1e-8
      @test Cosmology.luminosity_distance(cosmo, z)/Constants.DIST_MPC ≈ 6607.6576117749355 atol=1e-8 rtol=1e-8
      @test Cosmology.angular_diameter_distance(cosmo, 0.0, z)/Constants.DIST_MPC ≈ 1651.9144029437339 atol=1e-8 rtol=1e-8
      @test Cosmology.distance_modulus(cosmo, z) ≈ 44.10023765554372 atol=1e-8 rtol=1e-8
      @test Cosmology.angular_scale(cosmo, z) ≈ 8.008707025690128 atol=1e-8 rtol=1e-8

      @test Cosmology.comoving_volume_element(cosmo, z) ≈ 26.550755712253192 atol=1e-8 rtol=1e-8
      @test Cosmology.comoving_volume(cosmo, z) ≈ 151.05712532061932 atol=1e-8 rtol=1e-8
   end

   @testset "open LCDM" begin
      # Tests for open (Ωk > 0) LCDM cosmology module (AGAINST ASTROPY)
      cosmo = Cosmology.init_cosmology(Omega_m0=0.26)
      z = 1.0

      @test Cosmology.Ez(cosmo, z) ≈ 1.7146428199482249 atol=1e-15 rtol=1e-15
      @test Cosmology.hubble_parameter(cosmo, z) ≈ 120.02499739637574 atol=1e-15 rtol=1e-15

      @test Cosmology.age(cosmo, z) ≈ 6.022600655462329 atol=1e-8 rtol=1e-8
      @test Cosmology.lookback_time(cosmo, z) ≈ 7.796137793333031 atol=1e-8 rtol=1e-8

      @test Cosmology.rho_cz(cosmo, z) ≈ 2.7059389333539216e-26 atol=1e-15 rtol=1e-15

      @test Cosmology.comoving_distance_radial(cosmo, 0.0, z)/Constants.DIST_MPC ≈ 3343.5622067797044 atol=1e-8 rtol=1e-8
      @test Cosmology.comoving_distance_transverse(cosmo, 0.0, 1)/Constants.DIST_MPC ≈ 3357.164778591662 atol=1e-8 rtol=1e-8
      @test Cosmology.luminosity_distance(cosmo, z)/Constants.DIST_MPC ≈ 6714.329557183324 atol=1e-8 rtol=1e-8
      @test Cosmology.angular_diameter_distance(cosmo, 0.0, z)/Constants.DIST_MPC ≈ 1678.582389295831 atol=1e-8 rtol=1e-8

      @test Cosmology.comoving_volume_element(cosmo, z) ≈ 28.15103157576382 atol=1e-8 rtol=1e-8
      @test Cosmology.comoving_volume(cosmo, z) ≈ 157.3382098743564 atol=1e-8 rtol=1e-8
   end

   @testset "closed CDM" begin
      # Tests for closed (Ωk < 0) LCDM cosmology module (AGAINST ASTROPY)
      cosmo = Cosmology.init_cosmology(Omega_m0=0.34)
      z = 1.0

      @test Cosmology.Ez(cosmo, z) ≈ 1.8055470085267789 atol=1e-15 rtol=1e-15
      @test Cosmology.hubble_parameter(cosmo, z) ≈ 126.38829059687453 atol=1e-15 rtol=1e-15

      @test Cosmology.age(cosmo, z) ≈ 5.516141610883366 atol=1e-8 rtol=1e-8
      @test Cosmology.lookback_time(cosmo, z) ≈ 7.638261424137102 atol=1e-8 rtol=1e-8

      @test Cosmology.rho_cz(cosmo, z) ≈ 3.0004628988890413e-26 atol=1e-15 rtol=1e-15

      @test Cosmology.comoving_distance_radial(cosmo, 0.0, z)/Constants.DIST_MPC ≈ 3266.04304016069 atol=1e-8 rtol=1e-8
      @test Cosmology.comoving_distance_transverse(cosmo, 0.0, 1)/Constants.DIST_MPC ≈ 3253.394978599389 atol=1e-8 rtol=1e-8
      @test Cosmology.luminosity_distance(cosmo, z)/Constants.DIST_MPC ≈ 6506.789957198778 atol=1e-8 rtol=1e-8
      @test Cosmology.angular_diameter_distance(cosmo, 0.0, z)/Constants.DIST_MPC ≈ 1626.6974892996946 atol=1e-8 rtol=1e-8

      @test Cosmology.comoving_volume_element(cosmo, z) ≈ 25.10657361038687 atol=1e-8 rtol=1e-8
      @test Cosmology.comoving_volume(cosmo, z) ≈ 145.2558165125011 atol=1e-8 rtol=1e-8
   end
end