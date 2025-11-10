#!!!!!!!!!!!!!! Testing/cross-checked AGAINST Glafic !!!!!!!!!!!!!!
@testset "aHernquist lens" begin
   mass = 1E11 * MASS_SUN
   x_s = 0.3
   eps = 0.0
   pa = 0.0
   lens = Lenses.init_aHernquistLens(D_d=Dol, mass=mass, x_s=x_s, eps=eps, pa=pa)

   pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   kappa  = 0.5 * (jac1[1] + jac1[2])
   gamma1 = 0.5 * (jac1[1] - jac1[2])
   gamma2 = jac1[3]
   @test dex1[1] ≈ 0.135504 atol=1e-4 rtol=1e-4
   @test dex1[2] ≈ 0.135504 atol=1e-4 rtol=1e-4
   @test kappa  ≈ +0.018715 atol=1e-4 rtol=1e-4
   @test gamma1 ≈ +0.000000 atol=1e-4 rtol=1e-4
   @test gamma2 ≈ -0.116789 atol=1e-4 rtol=1e-4

   pot2 = adis  * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   kappa  = 0.5 * (jac2[1] + jac2[2])
   gamma1 = 0.5 * (jac2[1] - jac2[2])
   gamma2 = jac2[3]
   @test dex2[1] ≈ 0.242642 atol=1e-4 rtol=1e-4
   @test dex2[2] ≈ 0.000000 atol=1e-4 rtol=1e-4
   @test kappa  ≈ +0.044271 atol=1e-4 rtol=1e-4
   @test gamma1 ≈ -0.198371 atol=1e-4 rtol=1e-4
   @test gamma2 ≈ +0.000000 atol=1e-4 rtol=1e-4

   potc = adis .* Lenses.get_potential(lens, [xt1, xt2], [yt1, yt2])
   dexc = adis .* Lenses.get_deflection(lens, [xt1, xt2], [yt1, yt2])
   jacc = adis .* Lenses.get_jacobian(lens, [xt1, xt2], [yt1, yt2])
   @test potc[1] ≈ pot1 atol=1e-15 rtol=1e-15
   @test potc[2] ≈ pot2 atol=1e-15 rtol=1e-15

   @test dexc[1][1] ≈ dex1[1] atol=1e-15 rtol=1e-15
   @test dexc[2][1] ≈ dex1[2] atol=1e-15 rtol=1e-15
   @test dexc[1][2] ≈ dex2[1] atol=1e-15 rtol=1e-15
   @test dexc[2][2] ≈ dex2[2] atol=1e-15 rtol=1e-15

   @test jacc[1][1] ≈ jac1[1] atol=1e-15 rtol=1e-15
   @test jacc[1][2] ≈ jac2[1] atol=1e-15 rtol=1e-15
   @test jacc[2][1] ≈ jac1[2] atol=1e-15 rtol=1e-15
   @test jacc[2][2] ≈ jac2[2] atol=1e-15 rtol=1e-15
   @test jacc[3][1] ≈ jac1[3] atol=1e-15 rtol=1e-15
   @test jacc[3][2] ≈ jac2[3] atol=1e-15 rtol=1e-15


   # Get aHernquist lens with ellipticity
   eps = 0.3
   pa = 45.0
   lens = Lenses.init_aHernquistLens(D_d=Dol, mass=mass, x_s=x_s, eps=eps, pa=pa)
   pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   kappa  = 0.5 * (jac1[1] + jac1[2])
   gamma1 = 0.5 * (jac1[1] - jac1[2])
   gamma2 = jac1[3]
   @test dex1[1] ≈ 0.135978 atol=1e-4 rtol=1e-4
   @test dex1[2] ≈ 0.135978 atol=1e-4 rtol=1e-4
   @test kappa  ≈ +0.029328 atol=1e-4 rtol=1e-4
   @test gamma1 ≈ +0.000000 atol=1e-4 rtol=1e-4
   @test gamma2 ≈ -0.126183 atol=1e-4 rtol=1e-4

   pot2 = adis  * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   kappa  = 0.5 * (jac2[1] + jac2[2])
   gamma1 = 0.5 * (jac2[1] - jac2[2])
   gamma2 = jac2[3]
   @test dex2[1] ≈ +0.240930 atol=1e-4 rtol=1e-4
   @test dex2[2] ≈ -0.016681 atol=1e-4 rtol=1e-4
   @test kappa  ≈ +0.041059 atol=1e-4 rtol=1e-4
   @test gamma1 ≈ -0.193343 atol=1e-4 rtol=1e-4
   @test gamma2 ≈ +0.019136 atol=1e-4 rtol=1e-4

   potc = adis .* Lenses.get_potential(lens, [xt1, xt2], [yt1, yt2])
   dexc = adis .* Lenses.get_deflection(lens, [xt1, xt2], [yt1, yt2])
   jacc = adis .* Lenses.get_jacobian(lens, [xt1, xt2], [yt1, yt2])
   @test potc[1] ≈ pot1 atol=1e-15 rtol=1e-15
   @test potc[2] ≈ pot2 atol=1e-15 rtol=1e-15

   @test dexc[1][1] ≈ dex1[1] atol=1e-15 rtol=1e-15
   @test dexc[2][1] ≈ dex1[2] atol=1e-15 rtol=1e-15
   @test dexc[1][2] ≈ dex2[1] atol=1e-15 rtol=1e-15
   @test dexc[2][2] ≈ dex2[2] atol=1e-15 rtol=1e-15

   @test jacc[1][1] ≈ jac1[1] atol=1e-15 rtol=1e-15
   @test jacc[1][2] ≈ jac2[1] atol=1e-15 rtol=1e-15
   @test jacc[2][1] ≈ jac1[2] atol=1e-15 rtol=1e-15
   @test jacc[2][2] ≈ jac2[2] atol=1e-15 rtol=1e-15
   @test jacc[3][1] ≈ jac1[3] atol=1e-15 rtol=1e-15
   @test jacc[3][2] ≈ jac2[3] atol=1e-15 rtol=1e-15
end