#!!!!!!!!!!!!!! Testing/cross-checked AGAINST Glafic !!!!!!!!!!!!!!
@testset "Einasto lens" begin
   mass = 1E11 * MASS_SUN
   lens = Lenses.init_EinastoLens(cosmo, zl; mass=mass, c=6, n=0.2)

   pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   kappa = 0.5 * (jac1[1] + jac1[2])
   gamma1 = 0.5 * (jac1[1] - jac1[2])
   gamma2 = jac1[3]
   @test dex1[1] ≈ 0.031899 atol=1e-4 rtol=1e-4
   @test dex1[2] ≈ 0.031899 atol=1e-4 rtol=1e-4
   @test kappa  ≈ +0.018809 atol=1e-4 rtol=1e-4
   @test gamma1 ≈ +0.000000 atol=1e-4 rtol=1e-4
   @test gamma2 ≈ -0.013090 atol=1e-4 rtol=1e-4


   pot2 = adis  * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   kappa = 0.5 * (jac2[1] + jac2[2])
   gamma1 = 0.5 * (jac2[1] - jac2[2])
   gamma2 = jac2[3]
   @test dex2[1] ≈ 0.041787 atol=1e-4 rtol=1e-4
   @test dex2[2] ≈ 0.000000 atol=1e-4 rtol=1e-4
   @test kappa  ≈ +0.026430 atol=1e-4 rtol=1e-4
   @test gamma1 ≈ -0.015356 atol=1e-4 rtol=1e-4
   @test gamma2 ≈ -0.000000 atol=1e-4 rtol=1e-4
   
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

   @test_throws ArgumentError Lenses.init_EinastoLens(cosmo, zl, mass=mass)
end