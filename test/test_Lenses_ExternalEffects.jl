#!!!!!!!!!!!!!! Testing AGAINST formulae !!!!!!!!!!!!!!
@testset "ExternalEffects lens" begin
   kappa, gamma1, gamma2 = 0, 0.1, 0.1
   lens = Lenses.init_ExternalEffects(kappa=kappa, gamma1=gamma1, gamma2=gamma2)
   pot1 = Lenses.get_potential(lens, xt1, yt1)
   dex1 = Lenses.get_deflection(lens, xt1, yt1)
   jac1 = Lenses.get_jacobian(lens, xt1, yt1)
   @test pot1 ≈ 0.5 * (kappa + gamma1) * xt1^2 + 0.5 * (kappa - gamma1) * yt1^2 + gamma2 * xt1 * yt1 atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ (kappa + gamma1) * xt1 + gamma2 * yt1 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ (kappa - gamma1) * yt1 + gamma2 * xt1 atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ (kappa + gamma1) atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ (kappa - gamma1) atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ gamma2 atol=1e-15 rtol=1e-15

   kappa, gamma1, gamma2 = 0.5, 0.0, 0.1
   lens = Lenses.init_ExternalEffects(kappa=kappa, gamma1=gamma1, gamma2=gamma2)
   pot1 = Lenses.get_potential(lens, xt1, yt1)
   dex1 = Lenses.get_deflection(lens, xt1, yt1)
   jac1 = Lenses.get_jacobian(lens, xt1, yt1)
   @test pot1 ≈ 0.5 * (kappa + gamma1) * xt1^2 + 0.5 * (kappa - gamma1) * yt1^2 + gamma2 * xt1 * yt1 atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ (kappa + gamma1) * xt1 + gamma2 * yt1 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ (kappa - gamma1) * yt1 + gamma2 * xt1 atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ (kappa + gamma1) atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ (kappa - gamma1) atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ gamma2 atol=1e-15 rtol=1e-15

   kappa, gamma1, gamma2 = 0.5, 0.2, 0.0
   lens = Lenses.init_ExternalEffects(kappa=kappa, gamma1=gamma1, gamma2=gamma2)
   pot1 = Lenses.get_potential(lens, xt1, yt1)
   dex1 = Lenses.get_deflection(lens, xt1, yt1)
   jac1 = Lenses.get_jacobian(lens, xt1, yt1)
   @test pot1 ≈ 0.5 * (kappa + gamma1) * xt1^2 + 0.5 * (kappa - gamma1) * yt1^2 + gamma2 * xt1 * yt1 atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ (kappa + gamma1) * xt1 + gamma2 * yt1 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ (kappa - gamma1) * yt1 + gamma2 * xt1 atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ (kappa + gamma1) atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ (kappa - gamma1) atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ gamma2 atol=1e-15 rtol=1e-15

   kappa, gamma1, gamma2 = 0, 0.1, 0.1
   lens = Lenses.init_ExternalEffects(kappa=kappa, gamma1=gamma1, gamma2=gamma2)
   pot1 = Lenses.get_potential(lens, xt2, yt2)
   dex1 = Lenses.get_deflection(lens, xt2, yt2)
   jac1 = Lenses.get_jacobian(lens, xt2, yt2)
   @test pot1 ≈ 0.5 * (kappa + gamma1) * xt2^2 + 0.5 * (kappa - gamma1) * yt2^2 + gamma2 * xt2 * yt2 atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ (kappa + gamma1) * xt2 + gamma2 * yt2 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ (kappa - gamma1) * yt2 + gamma2 * xt2 atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ (kappa + gamma1) atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ (kappa - gamma1) atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ gamma2 atol=1e-15 rtol=1e-15

   kappa, gamma1, gamma2 = 0.5, 0.0, 0.1
   lens = Lenses.init_ExternalEffects(kappa=kappa, gamma1=gamma1, gamma2=gamma2)
   pot1 = Lenses.get_potential(lens, xt2, yt2)
   dex1 = Lenses.get_deflection(lens, xt2, yt2)
   jac1 = Lenses.get_jacobian(lens, xt2, yt2)
   @test pot1 ≈ 0.5 * (kappa + gamma1) * xt2^2 + 0.5 * (kappa - gamma1) * yt2^2 + gamma2 * xt2 * yt2 atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ (kappa + gamma1) * xt2 + gamma2 * yt2 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ (kappa - gamma1) * yt2 + gamma2 * xt2 atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ (kappa + gamma1) atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ (kappa - gamma1) atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ gamma2 atol=1e-15 rtol=1e-15

   kappa, gamma1, gamma2 = 0.5, 0.2, 0.0
   lens = Lenses.init_ExternalEffects(kappa=kappa, gamma1=gamma1, gamma2=gamma2)
   pot1 = Lenses.get_potential(lens, xt2, yt2)
   dex1 = Lenses.get_deflection(lens, xt2, yt2)
   jac1 = Lenses.get_jacobian(lens, xt2, yt2)
   @test pot1 ≈ 0.5 * (kappa + gamma1) * xt2^2 + 0.5 * (kappa - gamma1) * yt2^2 + gamma2 * xt2 * yt2 atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ (kappa + gamma1) * xt2 + gamma2 * yt2 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ (kappa - gamma1) * yt2 + gamma2 * xt2 atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ (kappa + gamma1) atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ (kappa - gamma1) atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ gamma2 atol=1e-15 rtol=1e-15
end