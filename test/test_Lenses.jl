# --------------------------------------------------------------------------------------------------
# Testing individual lens models
# --------------------------------------------------------------------------------------------------
# Cosmology   
cosmo = Cosmology.init_cosmology()

# Lens and source redshifts
zl = 0.5
zs = 1.5

# ADDs and distance ratio
Dol = Cosmology.angular_diameter_distance(cosmo, 0., zl)
Dls = Cosmology.angular_diameter_distance(cosmo, zl, zs)
Dos = Cosmology.angular_diameter_distance(cosmo, 0., zs)
adis = Dls/Dos

# Points (in image plane) to evaluate lensing quantities
xt1, yt1 = 1.0, 1.0
xt2, yt2 = 1.0, 0.0

@testset "Lenses (Individual)" begin
   include("./test_Lenses/test_Lenses_PointLens.jl")
   include("./test_Lenses/test_Lenses_PlummerLens.jl")
   include("./test_Lenses/test_Lenses_SISLens.jl")
   include("./test_Lenses/test_Lenses_NSISPLens.jl")
   include("./test_Lenses/test_Lenses_NSISMDLens.jl")
   include("./test_Lenses/test_Lenses_GaussianLens.jl")
   include("./test_Lenses/test_Lenses_SersicLens.jl")
   include("./test_Lenses/test_Lenses_PixelLens.jl")
   include("./test_Lenses/test_Lenses_ExternalEffects.jl")
   include("./test_Lenses/test_Lenses_ExternalEffects3.jl")
   include("./test_Lenses/test_Lenses_Multipole.jl")
   include("./test_Lenses/test_Lenses_PIEPLens.jl")
   include("./test_Lenses/test_Lenses_SIELens.jl")
   include("./test_Lenses/test_Lenses_PJELens.jl")
   include("./test_Lenses/test_Lenses_HernquistLens.jl")
   include("./test_Lenses/test_Lenses_eHernquistMDLens.jl")
   include("./test_Lenses/test_Lenses_aHernquistLens.jl")
   include("./test_Lenses/test_Lenses_NFWLens.jl")
   include("./test_Lenses/test_Lenses_eNFWMDLens.jl")
   include("./test_Lenses/test_Lenses_aNFWLens.jl")
   include("./test_Lenses/test_Lenses_tNFWLens.jl")
   # include("./test_Lenses/test_Lenses_gNFWLens.jl")
   include("./test_Lenses/test_Lenses_EinastoLens.jl")
   include("./test_Lenses/test_Lenses_MultiPlummerLens.jl")
   include("./test_Lenses/test_Lenses_MultiGaussianLens.jl")
   include("./test_Lenses/test_Lenses_MultiPixelLens.jl")
   include("./test_Lenses/test_Lenses_MultiPJELens.jl")
   include("./test_Lenses/test_Lenses_CompositeLens.jl")
end

# ------- Testing module-level API functions in src/Lenses/Lenses.jl -------------------------------
# These tests rely on the SIS lens, for which everything is known analytically:
#   - deflection amplitude    : α = θE (constant)
#   - convergence / shear     : κ = γ = θE / (2 r)
#   - tangential critical line: circle of radius θE
#   - radial critical line    : none (λ_rad = 1 everywhere)
#   - tangential caustic      : degenerate point at the origin
#   - images of source (β, 0) : (β + θE, 0) and (β - θE, 0) for β < θE
# A zero-amplitude ExternalEffects lens is used as an exact "no lensing" reference.

# Effective (adis-scaled) Einstein angle of the v_d = 200 km/s SIS lens
θE_raw = Lenses.SISLens.einstein_angle(D_ds=1.0, D_s=1.0, v_d=200.0)
θE_eff = adis * θE_raw

# Lenses used throughout
sis_lens  = Lenses.init_SISLens(v_d=200.0)
zero_lens = Lenses.init_ExternalEffects(kappa=0.0, gamma=0.0, angle=0.0)

# Image-plane grid; pixel size chosen so that no grid point hits (0, 0) exactly,
# where the SIS deflection is singular
grid_x, grid_y = Lenses.get_meshgrid(2.0, 2.0, 0.03)
ngx, ngy = size(grid_x)


@testset "Lenses (API)" begin
   @testset "Meshgrid" begin
      # 2D grid constructor: grid/pixel size <=0 error handeling
      @test_throws ArgumentError Lenses.get_meshgrid(1, 1, -1)
      @test_throws ArgumentError Lenses.get_meshgrid(1, 0, +1)

      # 2D grid constructor: size 
      xp, yp, pp = 1, 2, 0.3
      xg, yg = Lenses.get_meshgrid(xp, yp, pp)
      nx, ny = size(xg)
      @test size(xg) == size(yg)
      @test size(xg, 1) == size(collect(-xp:pp:xp), 1)
      @test size(xg, 2) == size(collect(-yp:pp:yp), 1)
   end

   @testset "Σ_cr" begin
      D_d, D_ds, D_s = 1.0, 1.0, 2.0
      Σcr = Lenses.get_critical_density(D_d, D_ds, D_s)
      
      # Σ_cr > 0 
      @test Σcr > 0.0

      # Unit conversion: kg/m^2 --> msun/pc^2
      @test Σcr * (DIST_PC^2 / MASS_SUN) ≈ Lenses.get_critical_density(D_d, D_ds, D_s; unit=:msun_pc2)

      # Unit conversion: kg/m^2 --> msun/arcsec^2
      @test Σcr * (D_d^2 * ANGLE_ARCSEC^2 / MASS_SUN) ≈ Lenses.get_critical_density(D_d, D_ds, D_s; unit=:msun_arcsec2)

      # Unknown unit error handeling
      @test_throws ArgumentError Lenses.get_critical_density(D_d, D_ds, D_s; unit=:invalid)
   end

   @testset "Input validation" begin
      # Mismatched coordinate sizes
      @test_throws ArgumentError Lenses.get_potential(sis_lens, [1.0, 2.0], [1.0])
      @test_throws ArgumentError Lenses.get_deflection(sis_lens, [1.0, 2.0], [1.0])
      @test_throws ArgumentError Lenses.get_jacobian(sis_lens, [1.0, 2.0], [1.0])

      # Int inputs are promoted to Float64 and give identical results
      @test Lenses.get_potential(sis_lens, 1, 1) ≈ Lenses.get_potential(sis_lens, 1.0, 1.0) atol=1e-15
      dint = Lenses.get_deflection(sis_lens, [1, 2], [1, 0])
      dflt = Lenses.get_deflection(sis_lens, [1.0, 2.0], [1.0, 0.0])
      @test dint[1] ≈ dflt[1] atol=1e-15
      @test dint[2] ≈ dflt[2] atol=1e-15
   end

   @testset "get_kappa_gamma" begin
      # Analytic SIS values at (1, 0): κ = γ = θE_eff/2, γ oriented along x ⇒ γ1 = -κ, γ2 = 0
      κ, γ1, γ2 = Lenses.get_kappa_gamma(sis_lens, 1.0, 0.0, adis)
      @test κ  ≈  0.5 * θE_eff atol=1e-14
      @test γ1 ≈ -0.5 * θE_eff atol=1e-14
      @test γ2 ≈  0.0          atol=1e-14

      # Consistency with the jacobian (scalar)
      ψxx, ψyy, ψxy = adis .* Lenses.get_jacobian(sis_lens, 1.0, 1.0)
      κb, γ1b, γ2b  = Lenses.get_kappa_gamma(sis_lens, 1.0, 1.0, adis)
      @test κb  ≈ 0.5 * (ψxx + ψyy) atol=1e-15
      @test γ1b ≈ 0.5 * (ψxx - ψyy) atol=1e-15
      @test γ2b ≈ ψxy               atol=1e-15

      # Array version matches scalar version
      κv, γ1v, γ2v = Lenses.get_kappa_gamma(sis_lens, [1.0, 1.0], [0.0, 1.0], adis)
      @test κv[1]  ≈ κ   atol=1e-15
      @test γ1v[1] ≈ γ1  atol=1e-15
      @test γ2v[1] ≈ γ2  atol=1e-15
      @test κv[2]  ≈ κb  atol=1e-15
      @test γ1v[2] ≈ γ1b atol=1e-15
      @test γ2v[2] ≈ γ2b atol=1e-15
   end

   @testset "get_time_delay" begin
      β = (0.3, -0.2)

      # No lens and θ = β  ⇒  zero time delay
      @test Lenses.get_time_delay(zero_lens, β[1], β[2], adis, zl, Dol, β) ≈ 0.0 atol=1e-25

      # No lens ⇒ pure geometric term
      cf = (1.0 + zl) / CONST_C * (Dol / adis) * ANGLE_ARCSEC^2
      @test Lenses.get_time_delay(zero_lens, 1.0, 1.0, adis, zl, Dol, β) ≈
            cf * 0.5 * ((1.0 - β[1])^2 + (1.0 - β[2])^2) rtol=1e-12

      # SIS: geometric + potential terms
      ψ = Lenses.get_potential(sis_lens, 1.0, 1.0)
      @test Lenses.get_time_delay(sis_lens, 1.0, 1.0, adis, zl, Dol, β) ≈
            cf * (0.5 * ((1.0 - β[1])^2 + (1.0 - β[2])^2) - adis * ψ) rtol=1e-12

      # Array version matches scalar version
      tdv = Lenses.get_time_delay(sis_lens, [1.0, -0.8], [1.0, 0.6], adis, zl, Dol, β)
      @test tdv[1] ≈ Lenses.get_time_delay(sis_lens,  1.0, 1.0, adis, zl, Dol, β) rtol=1e-14
      @test tdv[2] ≈ Lenses.get_time_delay(sis_lens, -0.8, 0.6, adis, zl, Dol, β) rtol=1e-14
   end

   @testset "get_magnification_image" begin
      # Analytic SIS magnification: μ = 1 / (1 - θE_eff/r)
      for (x, y) in ((1.0, 0.0), (1.0, 1.0), (-0.9, 0.4))
         r = hypot(x, y)
         @test Lenses.get_magnification_image(sis_lens, x, y, adis) ≈ 1.0 / (1.0 - θE_eff / r) atol=1e-12
      end

      # Array version matches scalar version
      μv = Lenses.get_magnification_image(sis_lens, [1.0, 1.0], [0.0, 1.0], adis)
      @test μv[1] ≈ Lenses.get_magnification_image(sis_lens, 1.0, 0.0, adis) atol=1e-15
      @test μv[2] ≈ Lenses.get_magnification_image(sis_lens, 1.0, 1.0, adis) atol=1e-15

      # Inside the Einstein radius the parity is negative
      @test Lenses.get_magnification_image(sis_lens, 0.5 * θE_eff, 0.0, adis) < 0.0
   end

   @testset "get_magnification_source" begin
      # Without a lens, every ray lands back on the grid: total "flux" is conserved
      μ_src = Lenses.get_magnification_source(zero_lens, grid_x, grid_y, adis)
      @test size(μ_src) == size(grid_x)
      @test all(μ_src .>= 0.0)
      @test sum(μ_src) ≈ ngx * ngy rtol=0.01

      # Same with more rays per pixel
      μ_src4 = Lenses.get_magnification_source(zero_lens, grid_x, grid_y, adis; rays_per_pixel=4)
      @test sum(μ_src4) ≈ ngx * ngy rtol=0.01

      # With a lens, some rays scatter outside the grid: flux can only decrease
      μ_sis = Lenses.get_magnification_source(sis_lens, grid_x, grid_y, adis)
      @test all(μ_sis .>= 0.0)
      @test sum(μ_sis) <= ngx * ngy + 1e-8
   end

   @testset "get_image (point)" begin
      β = (0.2, 0.0)
      images = Lenses.get_image(sis_lens, grid_x, grid_y, adis, β)

      # SIS with 0 < β < θE ⇒ two images at (β ± θE_eff, 0)
      @test length(images) >= 2
      @test any(p -> isapprox(p[1], β[1] + θE_eff; atol=0.05) && isapprox(p[2], 0.0; atol=0.05), images)
      @test any(p -> isapprox(p[1], β[1] - θE_eff; atol=0.05) && isapprox(p[2], 0.0; atol=0.05), images)
   end

   @testset "get_image (IRS)" begin
      # Gaussian source intensity map
      source_map = exp.(-0.5 .* (grid_x.^2 .+ grid_y.^2) ./ 0.3^2)

      # Without a lens, IRS must return the source map identically
      image_map = Lenses.get_image(zero_lens, grid_x, grid_y, adis, source_map)
      @test size(image_map) == size(source_map)
      @test image_map ≈ source_map atol=1e-14

      # With the SIS lens, the map changes but stays bounded by the source maximum
      image_sis = Lenses.get_image(sis_lens, grid_x, grid_y, adis, source_map)
      @test maximum(image_sis) <= maximum(source_map)
      @test minimum(image_sis) >= 0.0
   end

   @testset "CC, Caustics, and θE" begin
      crit_tan, crit_rad = Lenses.get_critical_curve(sis_lens, grid_x, grid_y, adis)

      # SIS: one tangential critical curve (circle of radius θE_eff), no radial one
      @test length(crit_tan) >= 1
      @test length(crit_rad) == 0
      for curve in crit_tan
         radii = hypot.(first.(curve), last.(curve))
         @test all(r -> isapprox(r, θE_eff; atol=0.02), radii)
      end

      # Tangential caustic degenerates to a point at the origin
      caus_tan, caus_rad = Lenses.get_caustic(sis_lens, grid_x, grid_y, adis)
      @test length(caus_tan) >= 1
      for curve in caus_tan
         @test all(p -> hypot(p[1], p[2]) < 0.05, curve)
      end

      # Area enclosed by the tangential critical curve and derived Einstein angle
      @test Lenses.get_critical_area(sis_lens, grid_x, grid_y, adis) ≈ π * θE_eff^2 rtol=0.01
      @test Lenses.get_einstein_angle(sis_lens, grid_x, grid_y, adis) ≈ θE_eff rtol=0.01
   end

   @testset "get_*_profile" begin
      # Constant convergence map ⇒ radial profile is flat (empty bins return 0)
      κ_const = fill(0.5, size(grid_x))
      for bin_type in (:linear, :log)
         centers, profile, edges = Lenses.get_radial_profile(κ_const, grid_x, grid_y;
                                                            origin=(0.0, 0.0), n_bin=40, bin_type=bin_type)
         @test length(edges)   == 40
         @test length(centers) == 39
         @test length(profile) == 39
         @test all(p -> p == 0.0 || isapprox(p, 0.5; atol=1e-12), profile)
         @test all(centers .> 0.0)
         @test issorted(edges)
      end

      # κ = 1 everywhere ⇒ M(<θ) = Σ_cr(z_s = ∞) ⋅ π θ²
      κ_one = fill(1.0, size(grid_x))
      radii, mass = Lenses.get_mass_profile(κ_one, grid_x, grid_y, Dol;
                                          origin=(0.0, 0.0), n_bin=40, bin_type=:linear)
      Σ_cr = Lenses.get_critical_density(Dol, 1.0, 1.0; unit=:msun_arcsec2)

      @test length(radii) == length(mass) == 39
      @test issorted(mass)  # cumulative mass must be non-decreasing

      # Compare at a radius fully contained inside the (square) grid
      idx = findlast(r -> r <= 1.5, radii)
      @test mass[idx] ≈ Σ_cr * π * radii[idx]^2 rtol=0.05
   end

   @testset "cartesian2polar" begin
      # Known values
      @test all(Lenses.shear_cartesian2polar(0.1, 0.0) .≈ (0.1,  0.0))
      @test all(Lenses.shear_cartesian2polar(0.0, 0.1) .≈ (0.1, 45.0))

      g1, g2 = Lenses.shear_polar2cartesian(0.1, 90.0)
      @test g1 ≈ -0.1 atol=1e-15
      @test g2 ≈  0.0 atol=1e-15

      # Roundtrips
      γ, φ = Lenses.shear_cartesian2polar(0.03, -0.04)
      γ1, γ2 = Lenses.shear_polar2cartesian(γ, φ)
      @test γ1 ≈  0.03 atol=1e-15
      @test γ2 ≈ -0.04 atol=1e-15

      e, φe = Lenses.ellipticity_cartesian2polar(-0.2, 0.1)
      e1, e2 = Lenses.ellipticity_polar2cartesian(e, φe)
      @test e1 ≈ -0.2 atol=1e-15
      @test e2 ≈  0.1 atol=1e-15

      # Shear and ellipticity conversions use the same spin-2 convention
      @test all(Lenses.ellipticity_cartesian2polar(0.05, 0.12) .≈ Lenses.shear_cartesian2polar(0.05, 0.12))
   end

   @testset "parameter_*" begin
      # NFW: c ↔ x_s roundtrip
      pN = Lenses.parameter_NFWLens(cosmology=cosmo, z_d=zl, mass=1.0E12, c=5.0)
      @test pN.x_s > 0.0
      @test pN.k_s > 0.0
      @test pN.rho_s > 0.0
      pN2 = Lenses.parameter_NFWLens(cosmology=cosmo, z_d=zl, mass=1.0E12, x_s=pN.x_s)
      @test pN2.c   ≈ 5.0     rtol=1e-12
      @test pN2.k_s ≈ pN.k_s rtol=1e-12

      # NFW: must provide either c or x_s
      @test_throws ArgumentError Lenses.parameter_NFWLens(cosmology=cosmo, z_d=zl, mass=1.0E12)

      # gNFW with n = 1 reduces to NFW
      pG = Lenses.parameter_gNFWLens(cosmology=cosmo, z_d=zl, mass=1.0E12, c=5.0, n=1.0)
      @test pG.x_s   ≈ pN.x_s   rtol=1e-8
      @test pG.rho_s ≈ pN.rho_s rtol=1e-8
      @test pG.k_s   ≈ pN.k_s   rtol=1e-8

      # gNFW: slope must lie in (0, 2)
      @test_throws ArgumentError Lenses.parameter_gNFWLens(cosmology=cosmo, z_d=zl, mass=1.0E12, c=5.0, n=2.5)
      @test_throws ArgumentError Lenses.parameter_gNFWLens(cosmology=cosmo, z_d=zl, mass=1.0E12, c=5.0, n=0.0)

      # Einasto: c ↔ x_s roundtrip
      pE = Lenses.parameter_EinastoLens(cosmology=cosmo, z_d=zl, mass=1.0E12, c=5.0, n=0.2)
      @test pE.x_s > 0.0
      @test pE.k_s > 0.0
      pE2 = Lenses.parameter_EinastoLens(cosmology=cosmo, z_d=zl, mass=1.0E12, x_s=pE.x_s, n=0.2)
      @test pE2.c ≈ 5.0 rtol=1e-12
      @test_throws ArgumentError Lenses.parameter_EinastoLens(cosmology=cosmo, z_d=zl, mass=1.0E12)
   end
end