#!!!!!!!!!!!!!! Testing AGAINST Glafic !!!!!!!!!!!!!!
#! Keep in mind that coordinates are rotated by 90 degree in Glafic compared to LensFactory. Hence,
#! 0° in LensFactory is equal to 90° in Glafic.
@testset "SIE lens" begin
   v_d = 200
   x_s = 0.5
   eps_old = 0.2
   eps_new = eps_old / (2.0 - eps_old)

   lens = Lenses.init_SIELens(v_d=v_d, x_s=x_s, eps=eps_new, pa=0)

   pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   kappa =  0.5 * (jac1[1] + jac1[2])
   gamma1 = 0.5 * (jac1[1] - jac1[2])
   gamma2 = jac1[3]
   @test dex1[1] ≈ 0.296095 atol=1e-4 rtol=1e-4
   @test dex1[2] ≈ 0.356473 atol=1e-4 rtol=1e-4
   @test kappa  ≈ 0.215727  atol=1e-4 rtol=1e-4
   @test gamma1 ≈ -0.006024 atol=1e-4 rtol=1e-4
   @test gamma2 ≈ -0.107085 atol=1e-4 rtol=1e-4

   pot2 = adis  * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   kappa =  0.5 * (jac2[1] + jac2[2])
   gamma1 = 0.5 * (jac2[1] - jac2[2])
   gamma2 = jac2[3]
   @test dex2[1] ≈ 0.376719 atol=1e-4 rtol=1e-4
   @test dex2[2] ≈ 0.000000 atol=1e-4 rtol=1e-4
   @test kappa  ≈ 0.319282  atol=1e-4 rtol=1e-4
   @test gamma1 ≈ -0.139998 atol=1e-4 rtol=1e-4
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


   lens = Lenses.init_SIELens(v_d=v_d, x_s=x_s, eps=eps_new, pa=30)
   pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   kappa =  0.5 * (jac1[1] + jac1[2])
   gamma1 = 0.5 * (jac1[1] - jac1[2])
   gamma2 = jac1[3]
   @test dex1[1] ≈ 0.294167 atol=1e-4 rtol=1e-4
   @test dex1[2] ≈ 0.326012 atol=1e-4 rtol=1e-4
   @test kappa  ≈ 0.236711  atol=1e-4 rtol=1e-4
   @test gamma1 ≈ -0.003419 atol=1e-4 rtol=1e-4
   @test gamma2 ≈ -0.127479 atol=1e-4 rtol=1e-4

   pot2 = adis  * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   kappa =  0.5 * (jac2[1] + jac2[2])
   gamma1 = 0.5 * (jac2[1] - jac2[2])
   gamma2 = jac2[3]
   @test dex2[1] ≈ 0.390306 atol=1e-4 rtol=1e-4
   @test dex2[2] ≈ -0.034752 atol=1e-4 rtol=1e-4
   @test kappa  ≈ 0.303440  atol=1e-4 rtol=1e-4
   @test gamma1 ≈ -0.123715 atol=1e-4 rtol=1e-4
   @test gamma2 ≈ -0.011411 atol=1e-4 rtol=1e-4
end