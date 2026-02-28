#!!!!!!!!!!!!!! Testing AGAINST Composite lens !!!!!!!!!!!!!!
@testset "Multi-Gaussian lens" begin
   xc1 = 0.0; yc1 = 0.0; mass1 = 1.0E11 * MASS_SUN; xs1 = 0.5
   xc2 = 0.5; yc2 = 0.5; mass2 = 1.0E11 * MASS_SUN; xs2 = 0.5
   
   lens1 = Lenses.init_CompositeLens([(lens=:GaussianLens, D_d=Dol, x_c=xc1, y_c=yc1, mass=mass1, x_s=xs1),
                                      (lens=:GaussianLens, D_d=Dol, x_c=xc2, y_c=yc2, mass=mass2, x_s=xs2)])
   pot11 = adis  * Lenses.get_potential(lens1, xt1, yt1)
   dex11 = adis .* Lenses.get_deflection(lens1, xt1, yt1)
   jac11 = adis .* Lenses.get_jacobian(lens1, xt1, yt1)
   mag11 = Lenses.get_magnification_image(lens1, xt1, yt1, adis)

   lens2 = Lenses.init_MultiGaussianLens(D_d=Dol, x_c=[xc1, xc2], y_c=[yc1, yc2], mass=[mass1, mass2], x_s=[xs1, xs2], n=2)
   pot21 = adis  * Lenses.get_potential(lens2, xt1, yt1)
   dex21 = adis .* Lenses.get_deflection(lens2, xt1, yt1)
   jac21 = adis .* Lenses.get_jacobian(lens2, xt1, yt1)
   mag21 = Lenses.get_magnification_image(lens2, xt1, yt1, adis)

   @test pot11 ≈ pot21  atol=1e-15 rtol=1e-15
   @test dex11[1] ≈ dex21[1]  atol=1e-15 rtol=1e-15
   @test dex11[2] ≈ dex21[2]  atol=1e-15 rtol=1e-15
   @test jac11[1] ≈ jac21[1]  atol=1e-15 rtol=1e-15
   @test jac11[2] ≈ jac21[2]  atol=1e-15 rtol=1e-15
   @test jac11[3] ≈ jac21[3]  atol=1e-15 rtol=1e-15


   pot12 = adis  * Lenses.get_potential(lens1, xt2, yt2)
   dex12 = adis .* Lenses.get_deflection(lens1, xt2, yt2)
   jac12 = adis .* Lenses.get_jacobian(lens1, xt2, yt2)
   mag12 = Lenses.get_magnification_image(lens1, xt2, yt2, adis)

   pot22 = adis  * Lenses.get_potential(lens2, xt2, yt2)
   dex22 = adis .* Lenses.get_deflection(lens2, xt2, yt2)
   jac22 = adis .* Lenses.get_jacobian(lens2, xt2, yt2)
   mag22 = Lenses.get_magnification_image(lens2, xt2, yt2, adis)

   @test pot12 ≈ pot22  atol=1e-15 rtol=1e-15
   @test dex12[1] ≈ dex22[1]  atol=1e-15 rtol=1e-15
   @test dex12[2] ≈ dex22[2]  atol=1e-15 rtol=1e-15
   @test jac12[1] ≈ jac22[1]  atol=1e-15 rtol=1e-15
   @test jac12[2] ≈ jac22[2]  atol=1e-15 rtol=1e-15
   @test jac12[3] ≈ jac22[3]  atol=1e-15 rtol=1e-15


   potc = adis .* Lenses.get_potential(lens2, [xt1, xt2], [yt1, yt2])
   dexc = adis .* Lenses.get_deflection(lens2, [xt1, xt2], [yt1, yt2])
   jacc = adis .* Lenses.get_jacobian(lens2, [xt1, xt2], [yt1, yt2])
   @test potc[1] ≈ pot21 atol=1e-15 rtol=1e-15
   @test potc[2] ≈ pot22 atol=1e-15 rtol=1e-15
   
   @test dexc[1][1] ≈ dex21[1] atol=1e-15 rtol=1e-15
   @test dexc[2][1] ≈ dex21[2] atol=1e-15 rtol=1e-15
   @test dexc[1][2] ≈ dex22[1] atol=1e-15 rtol=1e-15
   @test dexc[2][2] ≈ dex22[2] atol=1e-15 rtol=1e-15

   @test jacc[1][1] ≈ jac21[1] atol=1e-15 rtol=1e-15
   @test jacc[2][1] ≈ jac21[2] atol=1e-15 rtol=1e-15
   @test jacc[3][1] ≈ jac21[3] atol=1e-15 rtol=1e-15
   @test jacc[1][2] ≈ jac22[1] atol=1e-15 rtol=1e-15
   @test jacc[2][2] ≈ jac22[2] atol=1e-15 rtol=1e-15
   @test jacc[3][2] ≈ jac22[3] atol=1e-15 rtol=1e-15
end