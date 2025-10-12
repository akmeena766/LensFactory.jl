@testset "Gaussian lens" begin
   # lens properties
   mass = 1.0E11 * MASS_SUN
   x_s = 0.5
   lens = Lenses.init_GaussianLens(D_d=Dol, mass=mass, x_s=x_s)

   pot = adis * Lenses.get_potential(lens, xt1, yt1)
   dex = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   mag = Lenses.get_magnification_image(lens, xt1, yt1, adis)
   @test dex[1] ≈ 0.1800787586968404 atol=1e-15 rtol=1e-15
   @test dex[2] ≈ 0.1800787586968404 atol=1e-15 rtol=1e-15
   @test jac[1] ≈ +0.01343917717947224 atol=1e-15 rtol=1e-15
   @test jac[2] ≈ +0.01343917717947224 atol=1e-15 rtol=1e-15
   @test jac[3] ≈ -0.16663958151736816 atol=1e-15 rtol=1e-15
   @test mag ≈ 1.0576039797648285 atol=1e-15 rtol=1e-15


   pot = adis * Lenses.get_potential(lens, xt2, yt2)
   dex = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   mag = Lenses.get_magnification_image(lens, xt2, yt2, adis)
   @test dex[1] ≈ 0.3172256889321225 atol=1e-15 rtol=1e-15
   @test dex[2] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac[1] ≈ -0.11862002072694444 atol=1e-15 rtol=1e-15
   @test jac[2] ≈ +0.3172256889321225 atol=1e-15 rtol=1e-15
   @test jac[3] ≈ -0.0 atol=1e-15 rtol=1e-15
   @test mag ≈ 1.3093032302758327 atol=1e-15 rtol=1e-15

   # Lenstronomy gives different potential values. However, the relevant quantities is potential
   # difference, which is what I am going to test here.
   pot1 = adis * Lenses.get_potential(lens, xt1, yt1)
   pot2 = adis * Lenses.get_potential(lens, xt2, yt2)
   pot3 = adis * Lenses.get_potential(lens, 1.5, 2.5)
   @test pot1 - pot2 ≈ +0.118872955824665 atol=1e-15 rtol=1e-15
   @test pot1 - pot3 ≈ -0.26472744601185416 atol=1e-15 rtol=1e-15
   @test pot2 - pot3 ≈ -0.38360040183651917 atol=1e-15 rtol=1e-15
end