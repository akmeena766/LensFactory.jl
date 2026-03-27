#!!!!!!!!!!!!!! Testing AGAINST Glafic !!!!!!!!!!!!!!
#! Keep in mind that coordinates are rotated by 90 degree in Glafic compared to LensFactory. Hence,
#! 0° in LensFactory is equal to 90° in Glafic.
@testset "ExternalEffects3 lens" begin
   delta = 0.2
   angle = 90.0
   
   lens = Lenses.init_ExternalEffects3(delta=delta, angle=angle)
   pot1 = Lenses.get_potential(lens, xt1, yt1)
   dex1 = Lenses.get_deflection(lens, xt1, yt1)
   jac1 = Lenses.get_jacobian(lens, xt1, yt1)
   kappa =  0.5 * (jac1[1] + jac1[2])
   gamma1 = 0.5 * (jac1[1] - jac1[2])
   gamma2 = jac1[3]
   @test dex1[1] ≈ 0.4 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ 0.2 atol=1e-15 rtol=1e-15
   @test kappa  ≈ 0.2  atol=1e-15 rtol=1e-15
   @test gamma1 ≈ 0.2 atol=1e-15 rtol=1e-15
   @test gamma2 ≈ 0.4 atol=1e-15 rtol=1e-15

   pot2 = Lenses.get_potential(lens, xt2, yt2)
   dex2 = Lenses.get_deflection(lens, xt2, yt2)
   jac2 = Lenses.get_jacobian(lens, xt2, yt2)
   @test dex2[1] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test dex2[2] ≈ 0.2 atol=1e-15 rtol=1e-15
   @test jac2[1] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac2[2] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac2[3] ≈ 0.4 atol=1e-15 rtol=1e-15

   potc = Lenses.get_potential(lens, [xt1, xt2], [yt1, yt2])
   dexc = Lenses.get_deflection(lens, [xt1, xt2], [yt1, yt2])
   jacc = Lenses.get_jacobian(lens, [xt1, xt2], [yt1, yt2])
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