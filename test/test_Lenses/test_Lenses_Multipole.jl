@testset "Multipole" begin
   delta = 0.1
   angle = 90.0
   m = 2
   n = 2.0

   # First point
   lens = Lenses.init_Multipole(delta=delta, angle=angle, m=m, n=n)
   pot1 = Lenses.get_potential(lens,  xt1, yt1)
   dex1 = Lenses.get_deflection(lens, xt1, yt1)
   jac1 = Lenses.get_jacobian(lens,   xt1, yt1)

   lens2 = Lenses.init_ExternalEffects(kappa=0.0, gamma=delta, angle=angle)
   pot21 = Lenses.get_potential(lens2,  xt1, yt1)
   dex21 = Lenses.get_deflection(lens2, xt1, yt1)
   jac21 = Lenses.get_jacobian(lens2,   xt1, yt1)
   
   @test pot1 ≈ pot21 atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ dex21[1] atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ dex21[2] atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ jac21[1] atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ jac21[2] atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ jac21[3] atol=1e-15 rtol=1e-15

   # Test on origin
   pot1_ = Lenses.get_potential(lens,  0.0, 0.0)
   dex1_ = Lenses.get_deflection(lens, 0.0, 0.0)
   jac1_ = Lenses.get_jacobian(lens,   0.0, 0.0)
   @test pot1_ == 0.0
   @test dex1_ == (0.0, 0.0)
   @test jac1_ == (0.0, 0.0, 0.0)

   # Test on origin
   pot1_ = Lenses.get_potential(lens,  [0.0, 0.0], [0.0, 0.0])
   dex1_ = Lenses.get_deflection(lens, [0.0, 0.0], [0.0, 0.0])
   jac1_ = Lenses.get_jacobian(lens,   [0.0, 0.0], [0.0, 0.0])
   @test pot1_ == [0.0, 0.0]
   @test dex1_ == ([0.0, 0.0], [0.0, 0.0])
   @test jac1_ == ([0.0, 0.0], [0.0, 0.0], [0.0, 0.0])

   # Second point
   pot2 = Lenses.get_potential(lens,  xt2, yt2)
   dex2 = Lenses.get_deflection(lens, xt2, yt2)
   jac2 = Lenses.get_jacobian(lens,   xt2, yt2)

   pot22 = Lenses.get_potential(lens2,  xt2, yt2)
   dex22 = Lenses.get_deflection(lens2, xt2, yt2)
   jac22 = Lenses.get_jacobian(lens2,   xt2, yt2)
   
   @test pot2 ≈ pot22 atol=1e-15 rtol=1e-15
   @test dex2[1] ≈ dex22[1] atol=1e-15 rtol=1e-15
   @test dex2[2] ≈ dex22[2] atol=1e-15 rtol=1e-15
   @test jac2[1] ≈ jac22[1] atol=1e-15 rtol=1e-15
   @test jac2[2] ≈ jac22[2] atol=1e-15 rtol=1e-15
   @test jac2[3] ≈ jac22[3] atol=1e-15 rtol=1e-15
   
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

   # -----------------------------------------------------------------------------------------------
   #!!!!!!!!!!!!!! Testing AGAINST Glafic !!!!!!!!!!!!!!
   #! Keep in mind that coordinates are rotated by 90 degree in Glafic compared to LensFactory. 
   #! Hence, 0° in LensFactory is equal to 90° in Glafic.
   delta = 0.1
   angle = 90.0
   m = 3
   n = 2
   lens = Lenses.init_Multipole(delta=delta, angle=angle, m=m, n=n)   
   
   # First point
   pot1 = Lenses.get_potential(lens,  xt1, yt1)
   dex1 = Lenses.get_deflection(lens, xt1, yt1)
   jac1 = Lenses.get_jacobian(lens,   xt1, yt1)
   kappa =  0.5 * (jac1[1] + jac1[2])
   gamma1 = 0.5 * (jac1[1] - jac1[2])
   gamma2 = jac1[3]
   @test dex1[1] ≈ -0.117851 atol=1e-6 rtol=1e-6
   @test dex1[2] ≈ +0.023570 atol=1e-6 rtol=1e-6
   @test kappa  ≈ +0.058926  atol=1e-6 rtol=1e-6
   @test gamma1 ≈ -0.070711  atol=1e-6 rtol=1e-6
   @test gamma2 ≈ -0.106066  atol=1e-6 rtol=1e-6

   # Second point
   pot2 = Lenses.get_potential(lens,  xt2, yt2)
   dex2 = Lenses.get_deflection(lens, xt2, yt2)
   jac2 = Lenses.get_jacobian(lens,   xt2, yt2)
   kappa =  0.5 * (jac2[1] + jac2[2])
   gamma1 = 0.5 * (jac2[1] - jac2[2])
   gamma2 = jac2[3]
   @test dex2[1] ≈ -0.000000 atol=1e-6 rtol=1e-6
   @test dex2[2] ≈ -0.100000 atol=1e-6 rtol=1e-6
   @test kappa  ≈ +0.000000  atol=1e-6 rtol=1e-6
   @test gamma1 ≈ -0.000000  atol=1e-6 rtol=1e-6
   @test gamma2 ≈ -0.100000  atol=1e-6 rtol=1e-6

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
