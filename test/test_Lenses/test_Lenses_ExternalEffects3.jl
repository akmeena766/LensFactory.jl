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
   @test kappa  ≈ 0.2  atol=1e-4 rtol=1e-4
   @test gamma1 ≈ 0.2 atol=1e-4 rtol=1e-4
   @test gamma2 ≈ 0.4 atol=1e-4 rtol=1e-4


end