# --------------------------------------------------------------------------------------------------
# Test function exports
# --------------------------------------------------------------------------------------------------
using LensFactory.Cosmology

RTOL_TIGHT = 1.0e-6
RTOL_PHYS  = 1.0e-4
RTOL_RHOC  = 1.0e-6 


@testset "Cosmology" begin
   # -----------------------------------------------------------------------------------------------
   # Test function exports
   # -----------------------------------------------------------------------------------------------
   @testset "module exports" begin
      exported_names = (:init_cosmology, 
                        :scale_factor, :hubble_time, :age, :lookback_time, 
                        :hubble_distance, :hubble_parameter, 
                        :Omega_mz, :Omega_rz, :Omega_wz, :Omega_kz, :rho_cz, 
                        :comoving_distance_radial, 
                        :comoving_distance_transverse,
                        :angular_diameter_distance, 
                        :luminosity_distance, :distance_modulus,
                        :angular_scale, :comoving_volume_element, :comoving_volume, 
                        :adis2zs, :zs2adis)
      for nm in exported_names
         @test isdefined(Cosmology, nm)
         @test nm in names(Cosmology)
      end
   end


   # -----------------------------------------------------------------------------------------------
   # flat_lcdm:   standard flat Lambda-CDM                  (Om0=0.3, Ow0=0.7, Ok0=0 auto)
   # eds:         Einstein-de Sitter, matter-only, flat     (Om0=1.0)
   # milne:       empty universe, pure curvature            (Ok0=1.0, auto since all else=0)
   # rad_only:    pure radiation, flat                      (Or0=1.0)
   # open_univ:   open, self-consistent (sum=1)             (Ok0=+0.2, explicit)
   # closed_univ: closed, self-consistent (sum=1)           (Ok0=-0.2, explicit)
   # wcdm:        flat, non-(-1) equation of state          (w=-0.8)
   # -----------------------------------------------------------------------------------------------
   flat_lcdm   = init_cosmology(H0=70.0, w=-1.0, Omega_m0=0.3, Omega_r0=0.0, Omega_w0=0.7)
   eds         = init_cosmology(H0=70.0, w=-1.0, Omega_m0=1.0, Omega_r0=0.0, Omega_w0=0.0)
   milne       = init_cosmology(H0=70.0, w=-1.0, Omega_m0=0.0, Omega_r0=0.0, Omega_w0=0.0)
   rad_only    = init_cosmology(H0=70.0, w=-1.0, Omega_m0=0.0, Omega_r0=1.0, Omega_w0=0.0)
   open_univ   = init_cosmology(H0=70.0, w=-1.0, Omega_m0=0.3, Omega_r0=0.0, Omega_w0=0.5, Omega_k0=0.2)
   closed_univ = init_cosmology(H0=70.0, w=-1.0, Omega_m0=0.3, Omega_r0=0.0, Omega_w0=0.9, Omega_k0=-0.2)
   wcdm        = init_cosmology(H0=70.0, w=-0.8, Omega_m0=0.3, Omega_r0=0.0, Omega_w0=0.7)
   self_consistent_fixtures = (flat_lcdm, eds, milne, rad_only, open_univ, closed_univ, wcdm)


   # -----------------------------------------------------------------------------------------------
   @testset "init_cosmology default" begin
      # Default constructor: flat LCDM, H0=70, Om0=0.3, Ow0=0.7, Ok0 auto=0
      c = init_cosmology()
      @test c isa AbstractCosmology
      @test c.H0 == 70.0
      @test c.w == -1.0
      @test c.Omega_m0 == 0.3
      @test c.Omega_r0 == 0.0
      @test c.Omega_w0 == 0.7
      @test c.Omega_k0 == 0.0   # round(1.0 - 0.3 - 0.0 - 0.7, digits=6)

      # The `round(..., digits=6)` in the default should clean up floating-point noise
      c2 = init_cosmology(Omega_m0=0.25, Omega_w0=0.65)
      @test c2.Omega_k0 == 0.1

      # Explicit Omega_k0 overrides the auto-flat default
      c3 = init_cosmology(Omega_m0=0.3, Omega_w0=0.7, Omega_k0=0.2)
      @test c3.Omega_k0 == 0.2
      @test !isapprox(Cosmology.Ez(c3, 0.0), 1.0; rtol=1.0e-6)
      @test isapprox(Cosmology.Ez(c3, 0.0), 1.0954451150103321; rtol=RTOL_TIGHT)

      # An explicit Omega_k0 that happens to still sum to 1 is self-consistent.
      c4 = init_cosmology(Omega_m0=0.3, Omega_w0=0.6, Omega_k0=0.1)
      @test isapprox(Cosmology.Ez(c4, 0.0), 1.0; rtol=RTOL_TIGHT)

      # Spot-check the 7 fixtures used throughout this file.
      @test flat_lcdm.Omega_k0 == 0.0
      @test milne.Omega_k0     == 1.0
      @test rad_only.Omega_k0  == 0.0
      @test eds.Omega_k0       == 0.0
      @test wcdm.w             == -0.8
   end

   # -----------------------------------------------------------------------------------------------
   @testset "scale_factor" begin
      @test scale_factor(0.0) == 1.0
      @test scale_factor(1.0) == 0.5
      @test scale_factor(3.0) == 0.25
      @test scale_factor(9.0) == 0.1

      # a(z) * (1+z) == 1 for any z > -1
      for z in (0.0, 0.2, 1.0, 2.5, 10.0, 100.0)
         @test isapprox(scale_factor(z) * (1.0 + z), 1.0; rtol=RTOL_TIGHT)
      end
   end


   # -----------------------------------------------------------------------------------------------
   @testset "Ez, H(z)" begin
      # E(0) == 1 and H(0) == H0 for every self-consistent cosmology.
      for c in self_consistent_fixtures
         @test isapprox(Cosmology.Ez(c, 0.0), 1.0; rtol=RTOL_TIGHT)
         @test isapprox(hubble_parameter(c, 0.0), c.H0; rtol=RTOL_TIGHT)
      end

      # flat_lcdm: independently-derived E(z), H(z)
      @test isapprox(Cosmology.Ez(flat_lcdm, 0.5),     1.3086252328302401; rtol=RTOL_TIGHT)
      @test isapprox(Cosmology.Ez(flat_lcdm, 1.0),     1.7606816861659007; rtol=RTOL_TIGHT)
      @test isapprox(Cosmology.Ez(flat_lcdm, 2.0),     2.966479394838265;  rtol=RTOL_TIGHT)
      @test isapprox(hubble_parameter(flat_lcdm, 1.0), 123.24771803161305; rtol=RTOL_TIGHT)

      # eds (Om0=1, flat): E(z) = (1+z)^1.5
      @test isapprox(Cosmology.Ez(eds, 1.0), 2.8284271247461903; rtol=RTOL_TIGHT)
      @test isapprox(Cosmology.Ez(eds, 3.0), 8.0;                rtol=RTOL_TIGHT)

      # milne (Ok0=1 only): E(z) = 1+z, exactly
      for z in (0.0, 1.0, 4.0, 20.0)
         @test isapprox(Cosmology.Ez(milne, z), 1.0 + z; rtol=RTOL_TIGHT)
      end

      # rad_only (Or0=1, flat): E(z) = (1+z)^2, exactly
      for z in (0.0, 1.0, 3.0, 10.0)
         @test isapprox(Cosmology.Ez(rad_only, z), (1.0 + z)^2; rtol=RTOL_TIGHT)
      end

      # wcdm (w=-0.8, flat): independently-derived E(z), H(z)
      @test isapprox(Cosmology.Ez(wcdm, 1.0),     1.8603767351150353; rtol=RTOL_TIGHT)
      @test isapprox(Cosmology.Ez(wcdm, 2.0),     3.074610126740012;  rtol=RTOL_TIGHT)
      @test isapprox(hubble_parameter(wcdm, 1.0), 130.22637145805248; rtol=RTOL_TIGHT)

      # open_univ / closed_univ: independently-derived E(z)
      @test isapprox(Cosmology.Ez(open_univ, 1.0),   1.9235384061671346; rtol=RTOL_TIGHT)
      @test isapprox(Cosmology.Ez(closed_univ, 1.0), 1.5811388300841898; rtol=RTOL_TIGHT)

      # White-box check of the non-exported integrands (accessible via the
      # module-qualified name since they're defined inside `Cosmology`).
      for c in self_consistent_fixtures
         for z in (0.0, 0.5, 2.0)
            ez = Cosmology.Ez(c, z)
            @test isapprox(Cosmology._time_integrand(c, z), 1.0 / ((1.0 + z) * ez); rtol=RTOL_TIGHT)
            @test isapprox(Cosmology._distance_integrand(c, z), 1.0 / ez; rtol=RTOL_TIGHT)
         end
      end
   end


   # -----------------------------------------------------------------------------------------------
   @testset "hubble_time and hubble_distance" begin
      @test isapprox(hubble_time(70.0), 13.96846030972556;        rtol=RTOL_PHYS)
      @test isapprox(hubble_distance(70.0) / DIST_MPC, 4282.7494; rtol=RTOL_PHYS)

      # Both scale as 1/H0.
      @test isapprox(hubble_time(140.0), hubble_time(70.0) / 2.0;         rtol=RTOL_TIGHT)
      @test isapprox(hubble_distance(140.0), hubble_distance(70.0) / 2.0; rtol=RTOL_TIGHT)
   end


   # -----------------------------------------------------------------------------------------------
   @testset "Age and Lookback_time" begin
      # flat_lcdm: independently-derived reference ages.
      @test isapprox(age(flat_lcdm, 0.0), 13.466983947061875;         rtol=RTOL_PHYS)
      @test isapprox(age(flat_lcdm, 1.0), 5.751646943448284;          rtol=RTOL_PHYS)
      @test isapprox(lookback_time(flat_lcdm, 1.0), 7.71533700361359; rtol=RTOL_PHYS)

      # age(z) + lookback_time(z) == age(0) for any z and cosmology (the two
      # integrals split the same integrand at the point z).
      for c in self_consistent_fixtures
         a0 = age(c, 0.0)
         for z in (0.0, 0.3, 1.0, 2.0, 5.0)
            @test isapprox(age(c, z) + lookback_time(c, z), a0; rtol=1.0e-5)
         end
      end

      # Analytic closed-form ages for special-case cosmologies. age(0)/tH is
      # dimensionless and exact regardless of which constants define tH.
      @test isapprox(age(eds, 0.0) / hubble_time(eds.H0), 2.0 / 3.0;     rtol=RTOL_TIGHT)
      @test isapprox(age(milne, 0.0) / hubble_time(milne.H0), 1.0;       rtol=RTOL_TIGHT)
      @test isapprox(age(rad_only, 0.0) / hubble_time(rad_only.H0), 0.5; rtol=RTOL_TIGHT)

      # wcdm: independently-derived reference ages.
      @test isapprox(age(wcdm, 0.0), 12.999633440210864; rtol=RTOL_PHYS)
      @test isapprox(age(wcdm, 1.0), 5.583353720501408; rtol=RTOL_PHYS)

      # age decreases and lookback_time increases monotonically with z.
      zs   = (0.0, 0.5, 1.0, 2.0, 5.0)
      ages = [age(flat_lcdm, z) for z in zs]
      lbs  = [lookback_time(flat_lcdm, z) for z in zs]
      @test issorted(ages, rev=true)
      @test issorted(lbs)
   end

   # -----------------------------------------------------------------------------------------------
   @testset "rho_cz" begin
      @test isapprox(rho_cz(flat_lcdm, 0.0), 9.20387392297252e-27;                  rtol=RTOL_RHOC)
      @test isapprox(rho_cz(flat_lcdm, 0.0; unit=:kg_m3), 9.20387392297252e-27;     rtol=RTOL_RHOC)
      @test isapprox(rho_cz(flat_lcdm, 0.0; unit=:msun_pc3), 1.3599294735045703e-7; rtol=RTOL_RHOC)
      @test isapprox(rho_cz(flat_lcdm, 1.0), 2.8532009161214806e-26;                rtol=RTOL_RHOC)

      # rho_c scales as H(z)^2 regardless of which unit is requested.
      ratio_expected = (hubble_parameter(flat_lcdm, 1.0) / hubble_parameter(flat_lcdm, 0.0))^2
      @test isapprox(rho_cz(flat_lcdm, 1.0) / rho_cz(flat_lcdm, 0.0), ratio_expected; rtol=RTOL_TIGHT)
      ratio_expected_msun = rho_cz(flat_lcdm, 1.0; unit=:msun_pc3) / rho_cz(flat_lcdm, 0.0; unit=:msun_pc3)
      @test isapprox(ratio_expected_msun, ratio_expected; rtol=RTOL_TIGHT)

      # Invalid unit symbol must raise ArgumentError.
      @test_throws ArgumentError rho_cz(flat_lcdm, 0.0; unit=:bad_unit)
   end


    @testset "Omega_mz, Omega_rz, Omega_wz, Omega_kz" begin
      # Sum of the four density parameters is exactly 1 at every redshift, for
      # ANY cosmology -- a mathematical identity built into the definitions,
      # holding even when the cosmology is not self-consistent at z=0.
      inconsistent = init_cosmology(Omega_m0=0.3, Omega_w0=0.7, Omega_k0=0.2)  # sum != 1 at z=0
      for c in (self_consistent_fixtures..., inconsistent)
         for z in (0.0, 0.5, 1.0, 3.0, 10.0)
            total = Omega_mz(c, z) + Omega_rz(c, z) + Omega_wz(c, z) + Omega_kz(c, z)
            @test isapprox(total, 1.0; rtol=RTOL_TIGHT)
         end
      end

      # At z=0, each Omega_xz reduces to its corresponding Omega_x0 -- true for
      # any cosmology with E(0)=1, i.e. all 7 of the self-consistent fixtures.
      for c in self_consistent_fixtures
         @test isapprox(Omega_mz(c, 0.0), c.Omega_m0; rtol=RTOL_TIGHT)
         @test isapprox(Omega_rz(c, 0.0), c.Omega_r0; rtol=RTOL_TIGHT)
         @test isapprox(Omega_wz(c, 0.0), c.Omega_w0; rtol=RTOL_TIGHT)
         @test isapprox(Omega_kz(c, 0.0), c.Omega_k0; rtol=RTOL_TIGHT)
      end

      # Pure single-component universes: that component's Omega_xz stays
      # exactly 1 at every redshift (nothing else to dilute it).
      for z in (0.0, 1.0, 5.0, 50.0)
         @test isapprox(Omega_mz(eds, z), 1.0; rtol=RTOL_TIGHT)
         @test isapprox(Omega_kz(milne, z), 1.0; rtol=RTOL_TIGHT)
         @test isapprox(Omega_rz(rad_only, z), 1.0; rtol=RTOL_TIGHT)
      end

      # flat_lcdm: independently-derived component values at z=1.
      @test isapprox(Omega_mz(flat_lcdm, 1.0), 0.3 * 8.0 / 3.1; rtol=RTOL_TIGHT)
      @test isapprox(Omega_wz(flat_lcdm, 1.0), 0.7 / 3.1; rtol=RTOL_TIGHT)

      # open_univ: independently-derived component values at z=1 (Ez^2 = 3.7).
      @test isapprox(Omega_mz(open_univ, 1.0), 0.3 * 8.0 / 3.7; rtol=RTOL_TIGHT)
      @test isapprox(Omega_kz(open_univ, 1.0), 0.2 * 4.0 / 3.7; rtol=RTOL_TIGHT)
      @test isapprox(Omega_wz(open_univ, 1.0), 0.5 / 3.7; rtol=RTOL_TIGHT)
   end


   @testset "comoving_distance_radial" begin
      @test isapprox(comoving_distance_radial(flat_lcdm, 0.0, 1.0) / DIST_MPC, 3303.8288058874678; rtol=RTOL_PHYS)
      @test isapprox(comoving_distance_radial(flat_lcdm, 0.0, 3.0) / DIST_MPC, 6355.685436297465;  rtol=RTOL_PHYS)
      @test isapprox(comoving_distance_radial(wcdm, 0.0, 1.0) / DIST_MPC, 3165.897835253956;       rtol=RTOL_PHYS)

      # Additivity: D_C(z1,z2) + D_C(z2,z3) == D_C(z1,z3) (the integral over
      # [z1,z3] splits exactly at z2), for any cosmology.
      for c in self_consistent_fixtures
         d12 = comoving_distance_radial(c, 0.3, 0.8)
         d23 = comoving_distance_radial(c, 0.8, 1.5)
         d13 = comoving_distance_radial(c, 0.3, 1.5)
         @test isapprox(d12 + d23, d13; rtol=1.0e-5)
      end

      # Distance to itself is zero; swapping the redshift order negates it.
      @test isapprox(comoving_distance_radial(flat_lcdm, 1.0, 1.0), 0.0; atol=1.0)
      @test isapprox(comoving_distance_radial(flat_lcdm, 0.0, 1.0),
                     -comoving_distance_radial(flat_lcdm, 1.0, 0.0); rtol=1.0e-5)

      # Monotonically increasing in the upper limit.
      @test comoving_distance_radial(flat_lcdm, 0.0, 2.0) > comoving_distance_radial(flat_lcdm, 0.0, 1.0)
    end


    @testset "comoving_distance_transverse" begin
      # Flat universes: D_M == D_C exactly, for any z1, z2.
      for c in (flat_lcdm, eds, rad_only, wcdm)
         @test isapprox(comoving_distance_transverse(c, 0.0, 1.0), comoving_distance_radial(c, 0.0, 1.0); rtol=RTOL_TIGHT)
         @test isapprox(comoving_distance_transverse(c, 0.5, 2.0), comoving_distance_radial(c, 0.5, 2.0); rtol=RTOL_TIGHT)
      end

      # Open universe (Ωk0 > 0): D_M > D_C.
      @test comoving_distance_transverse(open_univ, 0.0, 1.0) > comoving_distance_radial(open_univ, 0.0, 1.0)
      @test isapprox(comoving_distance_radial(open_univ, 0.0, 1.0) / DIST_MPC,     3124.8660057639263; rtol=RTOL_PHYS)
      @test isapprox(comoving_distance_transverse(open_univ, 0.0, 1.0) / DIST_MPC, 3180.6153031869667; rtol=RTOL_PHYS)
      @test isapprox(comoving_distance_transverse(open_univ, 0.0, 2.0) / DIST_MPC, 5051.33366192017;   rtol=RTOL_PHYS)
      @test isapprox(comoving_distance_transverse(open_univ, 0.5, 2.0) / DIST_MPC, 3073.801723605768;  rtol=RTOL_PHYS)

      # Closed universe (Ωk0 < 0): D_M < D_C.
      @test comoving_distance_transverse(closed_univ, 0.0, 1.0) < comoving_distance_radial(closed_univ, 0.0, 1.0)
      @test isapprox(comoving_distance_radial(closed_univ, 0.0, 1.0) / DIST_MPC, 3527.7088629922177;     rtol=RTOL_PHYS)
      @test isapprox(comoving_distance_transverse(closed_univ, 0.0, 1.0) / DIST_MPC, 3448.4652020155568; rtol=RTOL_PHYS)
      @test isapprox(comoving_distance_transverse(closed_univ, 0.0, 2.0) / DIST_MPC, 5299.867073374233;  rtol=RTOL_PHYS)
      @test isapprox(comoving_distance_transverse(closed_univ, 0.5, 2.0) / DIST_MPC, 3559.008556395156;  rtol=RTOL_PHYS)

      # Continuity across the |Ωk0| > 1e-6 branch threshold: a tiny curvature
      # on either side of the threshold should produce a result extremely
      # close to the flat-universe value (no jump discontinuity).
      near_flat_pos = init_cosmology(Omega_m0=0.3, Omega_w0=0.7, Omega_k0=1.0e-5)
      near_flat_neg = init_cosmology(Omega_m0=0.3, Omega_w0=0.7, Omega_k0=-1.0e-5)
      flat_val = comoving_distance_transverse(flat_lcdm, 0.0, 1.0)
      @test isapprox(comoving_distance_transverse(near_flat_pos, 0.0, 1.0), flat_val; rtol=1.0e-4)
      @test isapprox(comoving_distance_transverse(near_flat_neg, 0.0, 1.0), flat_val; rtol=1.0e-4)
   
      # Milne (Ωk0 = 1): a fully analytic exercise of the open-curvature branch.
      # E(z) = 1+z  =>  D_C = D_H ln(1+z), and with Ωk = 1,
      # D_M(0,z) = D_H sinh(D_C/D_H) = D_H sinh(ln(1+z)) = D_H [(1+z) - 1/(1+z)] / 2.
      dH_milne = hubble_distance(milne.H0)
      for z in (0.5, 1.0, 3.0, 9.0)
         dM_analytic = dH_milne * ((1.0 + z) - 1.0 / (1.0 + z)) / 2.0
         @test isapprox(comoving_distance_transverse(milne, 0.0, z), dM_analytic; rtol=RTOL_PHYS)
         @test isapprox(comoving_distance_radial(milne, 0.0, z), dH_milne * log(1.0 + z); rtol=RTOL_PHYS)
      end

      # z1 > z2 gives a negative radial distance, and D_M carries the sign through
      # sinh (odd function), so the transverse distance is negated as well.
      @test isapprox(comoving_distance_transverse(open_univ, 1.0, 0.0),
                     -comoving_distance_transverse(open_univ, 0.0, 1.0); rtol=1.0e-5)

      # White-box: the internal _transverse_from_radial applied to the radial
      # comoving distance must reproduce comoving_distance_transverse(0, z) for
      # every curvature sign (flat, open, closed, Milne).
      for c in (flat_lcdm, open_univ, closed_univ, milne)
         for z in (0.5, 1.0, 2.0)
            χ = comoving_distance_radial(c, 0.0, z)
            @test isapprox(Cosmology._transverse_from_radial(c, χ),
                           comoving_distance_transverse(c, 0.0, z); rtol=RTOL_TIGHT)
         end
      end
   end


   @testset "luminosity_distance, angular_diameter_distance, distance_modulus" begin
      # flat_lcdm: independently-derived reference values.
      @test isapprox(luminosity_distance(flat_lcdm, 0.5) / DIST_MPC, 2832.938093900109;  rtol=RTOL_PHYS)
      @test isapprox(luminosity_distance(flat_lcdm, 1.0) / DIST_MPC, 6607.657611774936;  rtol=RTOL_PHYS)
      @test isapprox(luminosity_distance(flat_lcdm, 2.0) / DIST_MPC, 15539.586223228122; rtol=RTOL_PHYS)
      @test isapprox(luminosity_distance(flat_lcdm, 3.0) / DIST_MPC, 25422.741745189862; rtol=RTOL_PHYS)

      @test isapprox(angular_diameter_distance(flat_lcdm, 0.0, 0.5) / DIST_MPC, 1259.083597288937;  rtol=RTOL_PHYS)
      @test isapprox(angular_diameter_distance(flat_lcdm, 0.0, 1.0) / DIST_MPC, 1651.914402943734;  rtol=RTOL_PHYS)
      @test isapprox(angular_diameter_distance(flat_lcdm, 0.0, 2.0) / DIST_MPC, 1726.6206914697912; rtol=RTOL_PHYS)
      @test isapprox(angular_diameter_distance(flat_lcdm, 0.0, 3.0) / DIST_MPC, 1588.9213590743664; rtol=RTOL_PHYS)
      @test isapprox(angular_diameter_distance(flat_lcdm, 1.0, 2.0) / DIST_MPC, 625.344422840635;   rtol=RTOL_PHYS)

      @test isapprox(distance_modulus(flat_lcdm, 0.5), 42.261185421540894; rtol=RTOL_PHYS)
      @test isapprox(distance_modulus(flat_lcdm, 1.0), 44.10023765554372;  rtol=RTOL_PHYS)
      @test isapprox(distance_modulus(flat_lcdm, 2.0), 45.95719725271018;  rtol=RTOL_PHYS)
      @test isapprox(distance_modulus(flat_lcdm, 3.0), 47.026111928689645; rtol=RTOL_PHYS)

      @test isapprox(luminosity_distance(wcdm, 1.0) / DIST_MPC, 6331.795670507912; rtol=RTOL_PHYS)

      # Etherington distance-duality: D_L(z) == (1+z)^2 * D_A(0,z). This is an
      # algebraic consequence of how the two functions are both built from
      # comoving_distance_transverse(0,z), so it should hold for ANY cosmology.
      for c in self_consistent_fixtures
         for z in (0.5, 1.0, 2.0, 5.0)
            dl = luminosity_distance(c, z)
            da = angular_diameter_distance(c, 0.0, z)
            @test isapprox(dl, (1.0 + z)^2 * da; rtol=RTOL_TIGHT)
         end
      end

      # angular_diameter_distance is NOT additive in general (D_A(0,2) !=
      # D_A(0,1) + D_A(1,2)), unlike comoving_distance_radial.
      da01 = angular_diameter_distance(flat_lcdm, 0.0, 1.0)
      da12 = angular_diameter_distance(flat_lcdm, 1.0, 2.0)
      da02 = angular_diameter_distance(flat_lcdm, 0.0, 2.0)
      @test !isapprox(da01 + da12, da02; rtol=1.0e-3)
      @test isapprox(da01 + da12, 2277.258825784369 * DIST_MPC; rtol=RTOL_PHYS)

      # distance_modulus(z) == 5*log10(D_L/pc) - 5, derived directly from D_L.
      for z in (0.5, 1.0, 2.0)
         dl_pc = luminosity_distance(flat_lcdm, z) / DIST_PC
         @test isapprox(distance_modulus(flat_lcdm, z), 5.0 * log10(dl_pc) - 5.0; rtol=RTOL_TIGHT)
      end
   end


   @testset "angular_scale" begin
      @test isapprox(angular_scale(flat_lcdm, 0.5), 6.104209536262862; rtol=RTOL_PHYS)
      @test isapprox(angular_scale(flat_lcdm, 1.0), 8.00870702569013;  rtol=RTOL_PHYS)
      @test isapprox(angular_scale(flat_lcdm, 2.0), 8.370893333113619; rtol=RTOL_PHYS)

      # Direct definition check: d = D_A(0,z)/kpc * 1 arcsec [rad]
      for z in (0.5, 1.0, 2.0)
         da_kpc = angular_diameter_distance(flat_lcdm, 0.0, z) / DIST_KPC
         @test isapprox(angular_scale(flat_lcdm, z), da_kpc * ANGLE_ARCSEC; rtol=RTOL_TIGHT)
      end
   end


   @testset "comoving_volume_element and comoving_volume" begin
      # flat_lcdm: independently-derived reference values.
      @test isapprox(comoving_volume(flat_lcdm, 0.5), 28.217990639160497; rtol=RTOL_PHYS)
      @test isapprox(comoving_volume(flat_lcdm, 1.0), 151.05712532061938; rtol=RTOL_PHYS)
      @test isapprox(comoving_volume(flat_lcdm, 2.0), 582.1611191031681;  rtol=RTOL_PHYS)
      @test isapprox(comoving_volume(flat_lcdm, 3.0), 1075.414263975425;  rtol=RTOL_PHYS)

      @test isapprox(comoving_volume_element(flat_lcdm, 0.5), 11.673444513055363; rtol=RTOL_PHYS)
      @test isapprox(comoving_volume_element(flat_lcdm, 1.0), 26.55075571225319;  rtol=RTOL_PHYS)
      @test isapprox(comoving_volume_element(flat_lcdm, 2.0), 38.736262797681135; rtol=RTOL_PHYS)

      # Open / closed universe: independently-derived reference values, plus
      # confirmation that the closed-form comoving_volume formula stays
      # numerically well-behaved (no asin/asinh domain errors) for both signs
      # of curvature.
      @test isapprox(comoving_volume(open_univ, 1.0),           130.56487941799563; rtol=RTOL_PHYS)
      @test isapprox(comoving_volume_element(open_univ, 1.0),   22.523925864640308; rtol=RTOL_PHYS)
      @test isapprox(comoving_volume(closed_univ, 1.0),         178.96693726248202; rtol=RTOL_PHYS)
      @test isapprox(comoving_volume_element(closed_univ, 1.0), 32.211010875458136; rtol=RTOL_PHYS)

      # Flat universe: comoving_volume reduces to (4*pi/3)*D_M^3 exactly.
      for z in (0.5, 1.0, 2.0)
         dm = comoving_distance_transverse(flat_lcdm, 0.0, z) / DIST_GPC
         @test isapprox(comoving_volume(flat_lcdm, z), (4.0 * pi / 3.0) * dm^3; rtol=RTOL_TIGHT)
      end

      # Continuity across the |Ωk0| > 1e-6 branch threshold.
      near_flat_pos = init_cosmology(Omega_m0=0.3, Omega_w0=0.7, Omega_k0=1.0e-5)
      flat_val = comoving_volume(flat_lcdm, 1.0)
      @test isapprox(comoving_volume(near_flat_pos, 1.0), flat_val; rtol=1.0e-3)

      # comoving_volume should be monotonically increasing with z.
      @test comoving_volume(flat_lcdm, 2.0) > comoving_volume(flat_lcdm, 1.0) > comoving_volume(flat_lcdm, 0.5)
   end


   @testset "zs2adis" begin
      # zs2adis returns a_dis = D_ds/D_s = D_A(z_d, z_s)/D_A(0, z_s). The (1+z_s)
      # factors cancel, so it must equal the angular-diameter-distance ratio
      # exactly, for ANY cosmology and any z_d < z_s.
      for c in self_consistent_fixtures
         for (z_d, z_s) in ((0.3, 1.0), (0.5, 2.0), (1.0, 3.0), (0.2, 5.0))
            direct = angular_diameter_distance(c, z_d, z_s) / angular_diameter_distance(c, 0.0, z_s)
            @test isapprox(zs2adis(c, z_d, z_s), direct; rtol=RTOL_TIGHT)
         end
      end

      # Independently-derived reference ratios (these are the same anchors the
      # adis2zs round-trip uses below).
      @test isapprox(zs2adis(flat_lcdm, 0.5, 2.0), 0.6353907944259781; rtol=RTOL_PHYS)
      @test isapprox(zs2adis(flat_lcdm, 1.0, 3.0), 0.4801774192569025; rtol=RTOL_PHYS)

      # Passing the precomputed radial comoving distance χ(z_d) via `dC` must give
      # the identical result (it only skips the recomputation of one integral).
      χ_d = comoving_distance_radial(flat_lcdm, 0.0, 0.5)
      @test isapprox(zs2adis(flat_lcdm, 0.5, 2.0; dC=χ_d), zs2adis(flat_lcdm, 0.5, 2.0); rtol=RTOL_TIGHT)
      for c in (open_univ, closed_univ, wcdm)
         χ = comoving_distance_radial(c, 0.0, 0.4)
         @test isapprox(zs2adis(c, 0.4, 2.5; dC=χ), zs2adis(c, 0.4, 2.5); rtol=RTOL_TIGHT)
      end

      # a_dis lies in (0, 1) and, for a fixed lens, increases monotonically with
      # the source redshift toward a supremum below 1 (D_ds/D_s grows as the
      # source recedes, since χ_s saturates at the finite comoving horizon).
      ratios = [zs2adis(flat_lcdm, 0.5, z_s) for z_s in (0.6, 1.0, 2.0, 5.0, 20.0)]
      @test all(0.0 .< ratios .< 1.0)
      @test issorted(ratios)

      # Requires z_s > z_d: equal or inverted redshifts must error.
      @test_throws ErrorException zs2adis(flat_lcdm, 1.0, 1.0)
      @test_throws ErrorException zs2adis(flat_lcdm, 2.0, 1.0)
   end


   @testset "adis2zs" begin
      # Forward/backward round-trip: pick a target z_s, compute its true distance
      # ratio directly, then confirm adis2zs recovers z_s from that ratio.
      z_d, z_s_true = 0.5, 2.0
      adis_true = angular_diameter_distance(flat_lcdm, z_d, z_s_true) / angular_diameter_distance(flat_lcdm, 0.0, z_s_true)
      @test isapprox(adis_true, 0.6353907944259781; rtol=RTOL_PHYS)
      z_s_recovered = adis2zs(flat_lcdm, z_d, adis_true)
      @test isapprox(z_s_recovered, z_s_true; atol=1.0e-3)

      z_d2, z_s_true2 = 1.0, 3.0
      adis_true2 = angular_diameter_distance(flat_lcdm, z_d2, z_s_true2) / angular_diameter_distance(flat_lcdm, 0.0, z_s_true2)
      @test isapprox(adis_true2, 0.4801774192569025; rtol=RTOL_PHYS)
      z_s_recovered2 = adis2zs(flat_lcdm, z_d2, adis_true2)
      @test isapprox(z_s_recovered2, z_s_true2; atol=1.0e-3)

      # Forced non-convergence: for z_d=0.5 in flat_lcdm, the distance ratio
      # adis(z_s) -> ~0.8501 as z_s -> 100 (its supremum on the search domain),
      # so a target of 1.0 can never be bracketed and the bisection must run out
      # of iterations and raise an error.
      @test_throws ErrorException adis2zs(flat_lcdm, 0.5, 1.0; max_iter=50)

      # zs2adis / adis2zs are mutual inverses: for any cosmology and any lens,
      # feeding zs2adis's output back through adis2zs must recover the source
      # redshift. This exercises the bisection on open, closed and w-CDM models,
      # not just flat_lcdm.
      for c in self_consistent_fixtures
         for (z_d, z_s) in ((0.3, 1.5), (0.5, 2.0), (1.0, 4.0))
            adis = zs2adis(c, z_d, z_s)
            @test isapprox(adis2zs(c, z_d, adis), z_s; atol=1.0e-3)
         end
      end

      # A source just behind the lens (a_dis -> 1) should be recovered as z_s
      # slightly above z_d.
      z_close = 0.55
      adis_close = zs2adis(flat_lcdm, 0.5, z_close)
      @test isapprox(adis2zs(flat_lcdm, 0.5, adis_close), z_close; atol=1.0e-3)

      # Forced non-convergence: for z_d=0.5 in flat_lcdm, the distance ratio
      # adis(z_s) -> ~0.8501 as z_s -> 100 (its supremum on the search domain),
      # so a target of 1.0 can never be bracketed and the bisection must run out
      # of iterations and raise an error.
      @test_throws ErrorException adis2zs(flat_lcdm, 0.5, 1.0; max_iter=50)
   end
end