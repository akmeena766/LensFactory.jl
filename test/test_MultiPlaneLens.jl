@testset "MultiPlane" begin
   # ---- Fixtures ---------------------------------------------------------------------------------
   # flat ΛCDM (H0=70, Ωm=0.3, ΩΛ=0.7)
   cosmo = Cosmology.init_cosmology()

   # two lens planes + a source plane
   z1, z2, zs = 0.3, 0.6, 1.2
   Dd1 = Cosmology.angular_diameter_distance(cosmo, 0.0, z1)
   Dd2 = Cosmology.angular_diameter_distance(cosmo, 0.0, z2)

   # One point lens per plane, offset so the test point sits away from both singularities.
   multi_lens = [(lens = :PointLens, z_d = z1, D_d = Dd1, x_c = -1.0, y_c = 0.0, mass = 1.0e12),
                 (lens = :PointLens, z_d = z2, D_d = Dd2, x_c = +1.0, y_c = 0.0, mass = 1.0e12)]
   lens = Lenses.init_MultiPlaneLens(multi_lens)

   θx, θy = 0.3, 0.5                                   # evaluation point (arcsec)

   # Precomputed distance ratios (internal helper, shared by the core methods)
   D_ij, adis_ij, adis_is = MultiPlane._distances(cosmo, lens, zs)


   # ---- Potential --------------------------------------------------------------------------------
   @testset "potential: core vs cosmology, scalar vs array" begin
      ψ_c = MultiPlane.get_potential(lens, θx, θy, adis_ij, adis_is)
      ψ_w = MultiPlane.get_potential(cosmo, lens, θx, θy, zs)
      @test ψ_c ≈ ψ_w
      @test isfinite(ψ_w)

      X = fill(θx, 2, 2); Y = fill(θy, 2, 2)
      ψA = MultiPlane.get_potential(cosmo, lens, X, Y, zs)
      @test ψA[1, 1] ≈ ψ_w
   end


   # ---- Deflection -------------------------------------------------------------------------------
   @testset "deflection: core vs cosmology wrapper" begin
      αx_c, αy_c = MultiPlane.get_deflection(lens, θx, θy, adis_ij, adis_is)
      αx_w, αy_w = MultiPlane.get_deflection(cosmo, lens, θx, θy, zs)
      @test αx_c ≈ αx_w
      @test αy_c ≈ αy_w
      @test isfinite(αx_w) && isfinite(αy_w)
   end

   @testset "deflection: scalar vs array" begin
      αx_s, αy_s = MultiPlane.get_deflection(cosmo, lens, θx, θy, zs)
      X = fill(θx, 2, 2); Y = fill(θy, 2, 2)
      αxA, αyA = MultiPlane.get_deflection(cosmo, lens, X, Y, zs)
      @test size(αxA) == (2, 2)
      @test αxA[1, 1] ≈ αx_s
      @test αyA[1, 1] ≈ αy_s
   end


   # ---- Jacobian (deformation tensor) ------------------------------------------------------------
   @testset "jacobian: core vs cosmology, scalar vs array" begin
      j_c = MultiPlane.get_jacobian(lens, θx, θy, adis_ij, adis_is)     # (ψxx, ψyy, ψxy, ψyx)
      j_w = MultiPlane.get_jacobian(cosmo, lens, θx, θy, zs)
      @test length(j_w) == 4
      @test all(j_c .≈ j_w)

      X = fill(θx, 2, 2); Y = fill(θy, 2, 2)
      ψxxA, ψyyA, ψxyA, ψyxA = MultiPlane.get_jacobian(cosmo, lens, X, Y, zs)
      @test ψxxA[1, 1] ≈ j_w[1]
      @test ψyyA[1, 1] ≈ j_w[2]
      @test ψxyA[1, 1] ≈ j_w[3]
      @test ψyxA[1, 1] ≈ j_w[4]
   end

   @testset "jacobian == d(deflection)/dθ  (central differences)" begin
      h = 1.0e-4
      αxpx, αypx = MultiPlane.get_deflection(cosmo, lens, θx + h, θy, zs)
      αxmx, αymx = MultiPlane.get_deflection(cosmo, lens, θx - h, θy, zs)
      αxpy, αypy = MultiPlane.get_deflection(cosmo, lens, θx, θy + h, zs)
      αxmy, αymy = MultiPlane.get_deflection(cosmo, lens, θx, θy - h, zs)

      dαx_dx = (αxpx - αxmx) / (2h)
      dαy_dx = (αypx - αymx) / (2h)
      dαx_dy = (αxpy - αxmy) / (2h)
      dαy_dy = (αypy - αymy) / (2h)

      ψxx, ψyy, ψxy, ψyx = MultiPlane.get_jacobian(cosmo, lens, θx, θy, zs)
      @test ψxx ≈ dαx_dx atol = 1e-6
      @test ψyy ≈ dαy_dy atol = 1e-6
      @test ψxy ≈ dαx_dy atol = 1e-6   # ψ_xy = ∂α_x/∂θ_y
      @test ψyx ≈ dαy_dx atol = 1e-6   # ψ_yx = ∂α_y/∂θ_x (≠ ψ_xy in general for multi-plane)
   end


   # ---- Magnification (image plane) --------------------------------------------------------------
   @testset "magnification_image == 1/det(A)" begin
      ψxx, ψyy, ψxy, ψyx = MultiPlane.get_jacobian(cosmo, lens, θx, θy, zs)
      detA = 1.0 + ψxx * ψyy - ψxx - ψyy - ψxy * ψyx
      μ = MultiPlane.get_magnification_image(cosmo, lens, θx, θy, zs)
      @test isfinite(μ)
      @test μ ≈ 1.0 / detA
   end


   # ---- Fixtures for the additional API tests ----------------------------------------------------
   Dos = Cosmology.angular_diameter_distance(cosmo, 0.0, zs)
   D1s = Cosmology.angular_diameter_distance(cosmo, z1, zs)
   Dd2s = Cosmology.angular_diameter_distance(cosmo, z2, zs)
   adis1 = D1s / Dos

   # Multi-plane system whose 2nd plane is massless ⇒ must reduce to a single-plane point lens
   sp_limit_mp = Lenses.init_MultiPlaneLens(
      [(lens = :PointLens, z_d = z1, D_d = Dd1, x_c =  0.0, y_c =  0.0, mass = 1.0e11),
       (lens = :PointLens, z_d = z2, D_d = Dd2, x_c = 10.0, y_c = 10.0, mass = 0.0)])
   sp_limit_single = Lenses.init_PointLens(D_d = Dd1, x_c = 0.0, y_c = 0.0, mass = 1.0e11)

   # Fully massless two-plane system (exact "no lensing" reference)
   massless_mp = Lenses.init_MultiPlaneLens(
      [(lens = :PointLens, z_d = z1, D_d = Dd1, x_c = 10.0, y_c =  10.0, mass = 0.0),
       (lens = :PointLens, z_d = z2, D_d = Dd2, x_c = 10.0, y_c = -10.0, mass = 0.0)])

   # Two-plane system with moderate masses for the map-based tests (θE ≈ 0.6" each)
   map_mp = Lenses.init_MultiPlaneLens(
      [(lens = :PointLens, z_d = z1, D_d = Dd1, x_c = -1.0, y_c = 0.0, mass = 1.0e11),
       (lens = :PointLens, z_d = z2, D_d = Dd2, x_c = +1.0, y_c = 0.0, mass = 1.0e11)])

   # Image-plane grid; pixel size chosen so that no grid point hits the lens centers (±1, 0)
   gX, gY = Lenses.get_meshgrid(3.0, 3.0, 0.03)
   ngx, ngy = size(gX)


   # ---- Single-plane limit ------------------------------------------------------------------------
   @testset "single-plane limit (massless 2nd plane)" begin
      for (x, y) in ((0.3, 0.5), (-0.7, 0.2))
         # Potential
         @test MultiPlane.get_potential(cosmo, sp_limit_mp, x, y, zs) ≈
               adis1 * Lenses.get_potential(sp_limit_single, x, y) rtol=1e-10

         # Deflection
         αx_m, αy_m = MultiPlane.get_deflection(cosmo, sp_limit_mp, x, y, zs)
         αx_s, αy_s = Lenses.get_deflection(sp_limit_single, x, y)
         @test αx_m ≈ adis1 * αx_s rtol=1e-10
         @test αy_m ≈ adis1 * αy_s rtol=1e-10

         # Jacobian (must also be symmetric in this limit)
         ψxx_m, ψyy_m, ψxy_m, ψyx_m = MultiPlane.get_jacobian(cosmo, sp_limit_mp, x, y, zs)
         ψxx_s, ψyy_s, ψxy_s = Lenses.get_jacobian(sp_limit_single, x, y)
         @test ψxx_m ≈ adis1 * ψxx_s rtol=1e-10
         @test ψyy_m ≈ adis1 * ψyy_s rtol=1e-10
         @test ψxy_m ≈ adis1 * ψxy_s rtol=1e-10
         @test ψyx_m ≈ ψxy_m rtol=1e-10

         # Magnification
         @test MultiPlane.get_magnification_image(cosmo, sp_limit_mp, x, y, zs) ≈
               Lenses.get_magnification_image(sp_limit_single, x, y, adis1) rtol=1e-10
      end

      # Time delay: the multi-plane sum telescopes to the single-plane value in flat ΛCDM, but
      # only ON-SHELL, i.e., when θ and β are related by the lens equation. So derive β from θ.
      θpx, θpy = 0.3, 0.5
      αx_p, αy_p = Lenses.get_deflection(sp_limit_single, θpx, θpy)
      β_on = (θpx - adis1 * αx_p, θpy - adis1 * αy_p)
      @test MultiPlane.get_time_delay(cosmo, sp_limit_mp, θpx, θpy, zs, β_on) ≈
            Lenses.get_time_delay(sp_limit_single, θpx, θpy, adis1, z1, Dd1, β_on) rtol=1e-6
   end


   # ---- Time delay --------------------------------------------------------------------------------
   @testset "time delay: massless limit, scalar vs array, symmetry" begin
      β = (0.25, -0.15)

      # No deflectors and θ = β ⇒ zero time delay
      @test MultiPlane.get_time_delay(cosmo, massless_mp, β[1], β[2], zs, β) ≈ 0.0 atol=1e-20

      # No deflectors ⇒ pure geometric term, contributed by the last plane
      cf2 = (1.0 + z2) / CONST_C * (Dd2 * Dos / Dd2s) * ANGLE_ARCSEC^2
      @test MultiPlane.get_time_delay(cosmo, massless_mp, 1.0, 0.4, zs, β) ≈
            cf2 * 0.5 * ((1.0 - β[1])^2 + (0.4 - β[2])^2) rtol=1e-10

      # Scalar vs array versions agree (two-plane fixture with real masses)
      tds = MultiPlane.get_time_delay(cosmo, lens, 0.3, 0.5, zs, β)
      tda = MultiPlane.get_time_delay(cosmo, lens, fill(0.3, 2, 2), fill(0.5, 2, 2), zs, β)
      @test size(tda) == (2, 2)
      @test tda[1, 1] ≈ tds rtol=1e-12

      # Fixture lenses both sit on the x-axis ⇒ y → -y mirror symmetry
      @test MultiPlane.get_time_delay(cosmo, lens, 0.3, +0.5, zs, (β[1], +β[2])) ≈
            MultiPlane.get_time_delay(cosmo, lens, 0.3, -0.5, zs, (β[1], -β[2])) rtol=1e-12
   end


   # ---- Extended source imaging and IRS -----------------------------------------------------------
   @testset "extended source imaging and source magnification (IRS)" begin
      source_map = exp.(-0.5 .* (gX.^2 .+ gY.^2) ./ 0.5^2)

      # No deflectors: image map is identical to the source map, IRS conserves the flux
      @test MultiPlane.get_image(cosmo, massless_mp, gX, gY, zs, source_map) ≈ source_map atol=1e-14
      μ0 = MultiPlane.get_magnification_source(cosmo, massless_mp, gX, gY, zs)
      @test size(μ0) == size(gX)
      @test all(μ0 .>= 0.0)
      @test sum(μ0) ≈ ngx * ngy rtol=0.01

      # Real deflectors: maps stay bounded, IRS flux can only leave the grid
      img = MultiPlane.get_image(cosmo, map_mp, gX, gY, zs, source_map)
      @test maximum(img) <= maximum(source_map)
      @test minimum(img) >= 0.0
      μ1 = MultiPlane.get_magnification_source(cosmo, map_mp, gX, gY, zs; rays_per_pixel=2)
      @test all(μ1 .>= 0.0)
      @test sum(μ1) <= ngx * ngy + 1e-8
   end


   # ---- Point-source image positions --------------------------------------------------------------
   @testset "point-source image positions satisfy the lens equation" begin
      β = (0.1, 0.05)
      images = MultiPlane.get_image(cosmo, map_mp, gX, gY, zs, β)
      @test length(images) >= 1

      # Every image away from the lens centers must map back onto the source position
      n_checked = 0
      for p in images
         if hypot(p[1] + 1.0, p[2]) > 0.2 && hypot(p[1] - 1.0, p[2]) > 0.2
            αx, αy = MultiPlane.get_deflection(cosmo, map_mp, p[1], p[2], zs)
            @test hypot(β[1] - (p[1] - αx), β[2] - (p[2] - αy)) < 0.05
            n_checked += 1
         end
      end
      @test n_checked >= 1
   end


   # ---- Critical curves and caustics --------------------------------------------------------------
   @testset "critical curves and caustics" begin
      crit = MultiPlane.get_critical_curve(cosmo, map_mp, gX, gY, zs)
      @test length(crit) >= 1

      # det(A) ≈ 0 along the critical curves (checked on a subsample of points). Points within
      # 0.3" of a point-mass center are skipped: there the gridded det(A) field is far too steep
      # for the marching-squares interpolation to place the contour accurately.
      n_checked = 0
      for curve in crit
         for p in curve[1:5:end]
            if hypot(p[1] + 1.0, p[2]) > 0.3 && hypot(p[1] - 1.0, p[2]) > 0.3
               ψxx, ψyy, ψxy, ψyx = MultiPlane.get_jacobian(cosmo, map_mp, p[1], p[2], zs)
               @test abs(1.0 + ψxx * ψyy - ψxx - ψyy - ψxy * ψyx) < 0.1
               n_checked += 1
            end
         end
      end
      @test n_checked >= 1

      # One caustic per critical curve, with a matching number of points
      caus = MultiPlane.get_caustic(cosmo, map_mp, gX, gY, zs)
      @test length(caus) == length(crit)
      @test all(length(caus[i]) == length(crit[i]) for i in eachindex(crit))
   end
end
