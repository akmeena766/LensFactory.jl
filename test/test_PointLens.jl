@testset "Point lens" begin
   # Einstein angle
   D_d, D_ds, D_s, mass = 1.0, 1.0, 1.0, 1.0
   θE = Lenses.PointLens.einstein_angle(D_d=D_d, D_ds=D_ds, D_s=D_s, mass=mass)
   @test θE > 0.0
   @test θE ≈ sqrt(4.0 *CONST_G * mass / CONST_C^2 * (D_ds / D_d / D_s) ) / ANGLE_ARCSEC atol=1e-15 rtol=1e-15

   # Create a point lens
   lens = Lenses.init_PointLens(D_d=Dol, mass=1E11*MASS_SUN)

   pot = adis * Lenses.get_potential(lens, xt1, yt1)
   dex = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   mag = Lenses.get_magnification_image(lens, xt1, yt1, adis)

   @test pot ≈ 0.12714991581219895 atol=1e-15 rtol=1e-15
   @test dex[1] ≈ 0.18343855299170855 atol=1e-15 rtol=1e-15
   @test dex[2] ≈ 0.18343855299170855 atol=1e-15 rtol=1e-15
   @test jac[1] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac[2] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac[3] ≈ -0.18343855299170858 atol=1e-15 rtol=1e-15 
   @test mag ≈ 1.0348214336131885 atol=1e-15 rtol=1e-15

   pot = adis * Lenses.get_potential(lens, xt2, yt2)
   dex = adis .* Lenses.get_deflection(lens, xt1, yt2)
   jac = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   mag = Lenses.get_magnification_image(lens, xt2, yt2, adis)

   @test pot ≈ 0.0 atol=1e-15 rtol=1e-15
   @test dex[1] ≈ 0.36687710598341716 atol=1e-15 rtol=1e-15
   @test dex[2] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac[1] ≈ -0.36687710598341716 atol=1e-15 rtol=1e-15
   @test jac[2] ≈ +0.36687710598341716 atol=1e-15 rtol=1e-15
   @test jac[3] ≈ 0.0 atol=1e-15 rtol=1e-15 
   @test mag ≈ 1.155533424947028 atol=1e-15 rtol=1e-15
end