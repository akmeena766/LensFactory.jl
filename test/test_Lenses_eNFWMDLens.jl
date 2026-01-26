#!!!!!!!!!!!!!! Testing/cross-checked AGAINST Glafic !!!!!!!!!!!!!!
@testset "eNFWMD lens" begin
   mass = 1E11 * MASS_SUN
   eps_old = 0.0
   eps_new = eps_old / (2.0 - eps_old)
   # Get NFW lens parameters
   param = Lenses.parameter_NFWLens(cosmology=cosmo, z_d=zl, mass=mass, c=6)
   lens = Lenses.init_eNFWMDLens(D_d=Dol, rho_s=param.rho_s, x_s=param.x_s, eps=eps_new, pa=0.0)

   lens1 = Lenses.init_NFWLens(D_d=Dol, rho_s=param.rho_s, x_s=param.x_s)
   pot1_1 = adis  * Lenses.get_potential(lens1, xt1, yt1)
   dex1_1 = adis .* Lenses.get_deflection(lens1, xt1, yt1)
   jac1_1 = adis .* Lenses.get_jacobian(lens1, xt1, yt1)
   pot2_1 = adis  * Lenses.get_potential(lens1, xt2, yt2)
   dex2_1 = adis .* Lenses.get_deflection(lens1, xt2, yt2)
   jac2_1 = adis .* Lenses.get_jacobian(lens1, xt2, yt2)

   pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   @test pot1 ≈ pot1_1 atol=1e-6 rtol=1e-6
   @test dex1[1] ≈ dex1_1[1] atol=1e-6 rtol=1e-6
   @test dex1[2] ≈ dex1_1[2] atol=1e-6 rtol=1e-6
   @test jac1[1] ≈ jac1_1[1] atol=1e-6 rtol=1e-6
   @test jac1[2] ≈ jac1_1[2] atol=1e-6 rtol=1e-6
   @test jac1[3] ≈ jac1_1[3] atol=1e-6 rtol=1e-6

   pot2 = adis  * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   @test pot2 ≈ pot2_1 atol=1e-6 rtol=1e-6
   @test dex2[1] ≈ dex2_1[1] atol=1e-6 rtol=1e-6
   @test dex2[2] ≈ dex2_1[2] atol=1e-6 rtol=1e-6
   @test jac2[1] ≈ jac2_1[1] atol=1e-6 rtol=1e-6
   @test jac2[2] ≈ jac2_1[2] atol=1e-6 rtol=1e-6
   @test jac2[3] ≈ jac2_1[3] atol=1e-6 rtol=1e-6

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

   # Get eNFWMD lens with ellipticity
   eps_old = 0.3
   eps_new = eps_old / (2.0 - eps_old)
   lens = Lenses.init_eNFWMDLens(D_d=Dol, rho_s=param.rho_s, x_s=param.x_s, eps=eps_new, pa=45.0)

   # pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   kappa  = 0.5 * (jac1[1] + jac1[2])
   gamma1 = 0.5 * (jac1[1] - jac1[2])
   gamma2 = jac1[3]

   @test dex1[1] ≈ 0.028554 atol=1e-4 rtol=1e-4
   @test dex1[2] ≈ 0.028554 atol=1e-4 rtol=1e-4
   @test kappa  ≈ +0.022748 atol=1e-4 rtol=1e-4
   @test gamma1 ≈ +0.000000 atol=1e-4 rtol=1e-4
   @test gamma2 ≈ -0.015482 atol=1e-4 rtol=1e-4

   # pot2 = adis  * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   kappa  = 0.5 * (jac2[1] + jac2[2])
   gamma1 = 0.5 * (jac2[1] - jac2[2])
   gamma2 = jac2[3]
   @test dex2[1] ≈ +0.040098 atol=1e-4 rtol=1e-4
   @test dex2[2] ≈ -0.005782 atol=1e-4 rtol=1e-4
   @test kappa  ≈ +0.025821 atol=1e-4 rtol=1e-4
   @test gamma1 ≈ -0.013248 atol=1e-4 rtol=1e-4
   @test gamma2 ≈ -0.001227 atol=1e-4 rtol=1e-4

   # potc = adis .* Lenses.get_potential(lens, [xt1, xt2], [yt1, yt2])
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