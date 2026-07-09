@testset "MultiPlane" begin
   # ---- Fixtures ---------------------------------------------------------------------------------
   cosmo = Cosmology.init_cosmology()                 # flat ΛCDM (H0=70, Ωm=0.3, ΩΛ=0.7)

   z1, z2, zs = 0.3, 0.6, 1.2                          # two lens planes + a source plane
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
end
