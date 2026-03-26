#!!!!!!!!!!!!!! Testing AGAINST formulae !!!!!!!!!!!!!!
@testset "ExternalEffects lens" begin
   kappa = 0.0
   gamma = 0.1
   angle = 45.0
   gamma1 = gamma * cos(2.0 * deg2rad(angle))
   gamma2 = gamma * sin(2.0 * deg2rad(angle))


   lens = Lenses.init_ExternalEffects(kappa=kappa, gamma=gamma, angle=angle)
   pot1 = Lenses.get_potential(lens, xt1, yt1)
   dex1 = Lenses.get_deflection(lens, xt1, yt1)
   jac1 = Lenses.get_jacobian(lens, xt1, yt1)
   @test pot1 ≈ 0.5 * (kappa + gamma1) * xt1^2 + 0.5 * (kappa - gamma1) * yt1^2 + gamma2 * xt1 * yt1 atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ (kappa + gamma1) * xt1 + gamma2 * yt1 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ (kappa - gamma1) * yt1 + gamma2 * xt1 atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ (kappa + gamma1) atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ (kappa - gamma1) atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ gamma2 atol=1e-15 rtol=1e-15

   pot2 = Lenses.get_potential(lens, xt2, yt2)
   dex2 = Lenses.get_deflection(lens, xt2, yt2)
   jac2 = Lenses.get_jacobian(lens, xt2, yt2)
   @test pot2 ≈ 0.5 * (kappa + gamma1) * xt2^2 + 0.5 * (kappa - gamma1) * yt2^2 + gamma2 * xt2 * yt2 atol=1e-15 rtol=1e-15
   @test dex2[1] ≈ (kappa + gamma1) * xt2 + gamma2 * yt2 atol=1e-15 rtol=1e-15
   @test dex2[2] ≈ (kappa - gamma1) * yt2 + gamma2 * xt2 atol=1e-15 rtol=1e-15
   @test jac2[1] ≈ (kappa + gamma1) atol=1e-15 rtol=1e-15
   @test jac2[2] ≈ (kappa - gamma1) atol=1e-15 rtol=1e-15
   @test jac2[3] ≈ gamma2 atol=1e-15 rtol=1e-15

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

   
   kappa = 0.0
   gamma = 0.1
   angle = 30.0
   gamma1 = gamma * cos(2.0 * deg2rad(angle))
   gamma2 = gamma * sin(2.0 * deg2rad(angle))

   lens = Lenses.init_ExternalEffects(kappa=kappa, gamma=gamma, angle=angle)
   pot1 = Lenses.get_potential(lens, xt1, yt1)
   dex1 = Lenses.get_deflection(lens, xt1, yt1)
   jac1 = Lenses.get_jacobian(lens, xt1, yt1)
   @test pot1 ≈ 0.5 * (kappa + gamma1) * xt1^2 + 0.5 * (kappa - gamma1) * yt1^2 + gamma2 * xt1 * yt1 atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ (kappa + gamma1) * xt1 + gamma2 * yt1 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ (kappa - gamma1) * yt1 + gamma2 * xt1 atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ (kappa + gamma1) atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ (kappa - gamma1) atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ gamma2 atol=1e-15 rtol=1e-15

   pot2 = Lenses.get_potential(lens, xt2, yt2)
   dex2 = Lenses.get_deflection(lens, xt2, yt2)
   jac2 = Lenses.get_jacobian(lens, xt2, yt2)
   @test pot2 ≈ 0.5 * (kappa + gamma1) * xt2^2 + 0.5 * (kappa - gamma1) * yt2^2 + gamma2 * xt2 * yt2 atol=1e-15 rtol=1e-15
   @test dex2[1] ≈ (kappa + gamma1) * xt2 + gamma2 * yt2 atol=1e-15 rtol=1e-15
   @test dex2[2] ≈ (kappa - gamma1) * yt2 + gamma2 * xt2 atol=1e-15 rtol=1e-15
   @test jac2[1] ≈ (kappa + gamma1) atol=1e-15 rtol=1e-15
   @test jac2[2] ≈ (kappa - gamma1) atol=1e-15 rtol=1e-15
   @test jac2[3] ≈ gamma2 atol=1e-15 rtol=1e-15

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


   kappa = 0.0
   gamma = 0.3
   angle = 45
   gamma1 = gamma * cos(2.0 * deg2rad(angle))
   gamma2 = gamma * sin(2.0 * deg2rad(angle))

   lens = Lenses.init_ExternalEffects(kappa=kappa, gamma=gamma, angle=angle)
   pot1 = Lenses.get_potential(lens, xt1, yt1)
   dex1 = Lenses.get_deflection(lens, xt1, yt1)
   jac1 = Lenses.get_jacobian(lens, xt1, yt1)
   @test pot1 ≈ 0.5 * (kappa + gamma1) * xt1^2 + 0.5 * (kappa - gamma1) * yt1^2 + gamma2 * xt1 * yt1 atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ (kappa + gamma1) * xt1 + gamma2 * yt1 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ (kappa - gamma1) * yt1 + gamma2 * xt1 atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ (kappa + gamma1) atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ (kappa - gamma1) atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ gamma2 atol=1e-15 rtol=1e-15

   pot2 = Lenses.get_potential(lens, xt2, yt2)
   dex2 = Lenses.get_deflection(lens, xt2, yt2)
   jac2 = Lenses.get_jacobian(lens, xt2, yt2)
   @test pot2 ≈ 0.5 * (kappa + gamma1) * xt2^2 + 0.5 * (kappa - gamma1) * yt2^2 + gamma2 * xt2 * yt2 atol=1e-15 rtol=1e-15
   @test dex2[1] ≈ (kappa + gamma1) * xt2 + gamma2 * yt2 atol=1e-15 rtol=1e-15
   @test dex2[2] ≈ (kappa - gamma1) * yt2 + gamma2 * xt2 atol=1e-15 rtol=1e-15
   @test jac2[1] ≈ (kappa + gamma1) atol=1e-15 rtol=1e-15
   @test jac2[2] ≈ (kappa - gamma1) atol=1e-15 rtol=1e-15
   @test jac2[3] ≈ gamma2 atol=1e-15 rtol=1e-15

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
   
   
   kappa = 0.0
   gamma = 0.5
   angle = 45
   
   gamma1 = gamma * cos(2.0 * deg2rad(angle))
   gamma2 = gamma * sin(2.0 * deg2rad(angle))

   lens = Lenses.init_ExternalEffects(kappa=kappa, gamma=gamma, angle=angle)
   pot1 = Lenses.get_potential(lens, xt1, yt1)
   dex1 = Lenses.get_deflection(lens, xt1, yt1)
   jac1 = Lenses.get_jacobian(lens, xt1, yt1)
   @test pot1 ≈ 0.5 * (kappa + gamma1) * xt1^2 + 0.5 * (kappa - gamma1) * yt1^2 + gamma2 * xt1 * yt1 atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ (kappa + gamma1) * xt1 + gamma2 * yt1 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ (kappa - gamma1) * yt1 + gamma2 * xt1 atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ (kappa + gamma1) atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ (kappa - gamma1) atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ gamma2 atol=1e-15 rtol=1e-15

   pot2 = Lenses.get_potential(lens, xt2, yt2)
   dex2 = Lenses.get_deflection(lens, xt2, yt2)
   jac2 = Lenses.get_jacobian(lens, xt2, yt2)
   @test pot2 ≈ 0.5 * (kappa + gamma1) * xt2^2 + 0.5 * (kappa - gamma1) * yt2^2 + gamma2 * xt2 * yt2 atol=1e-15 rtol=1e-15
   @test dex2[1] ≈ (kappa + gamma1) * xt2 + gamma2 * yt2 atol=1e-15 rtol=1e-15
   @test dex2[2] ≈ (kappa - gamma1) * yt2 + gamma2 * xt2 atol=1e-15 rtol=1e-15
   @test jac2[1] ≈ (kappa + gamma1) atol=1e-15 rtol=1e-15
   @test jac2[2] ≈ (kappa - gamma1) atol=1e-15 rtol=1e-15
   @test jac2[3] ≈ gamma2 atol=1e-15 rtol=1e-15

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