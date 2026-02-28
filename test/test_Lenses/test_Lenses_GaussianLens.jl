#!!!!!!!!!!!!!! Testing AGAINST Lenstronomy + Glafic !!!!!!!!!!!!!!
# lensModel = LensModel(lens_model_list=['GAUSSIAN'])
# gaussian = {'center_x': 0, 'center_y': 0, 'sigma':0.5, 'amp': 0.18343855299170853 * 2 * np.pi}
# tx, ty = [1.0, 1.0], [1.0, 0.0]
# pot = lensModel.potential(tx, ty, [gaussian])
# dex = lensModel.alpha(tx, ty, [gaussian])
# hes = lensModel.hessian(tx, ty, [gaussian])
# mag = lensModel.magnification(tx, ty, [gaussian])

@testset "Gaussian lens" begin
   # lens properties
   mass = 1.0E11 * MASS_SUN
   x_s = 0.5
   lens = Lenses.init_GaussianLens(D_d=Dol, mass=mass, x_s=x_s)

   pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   mag1 = Lenses.get_magnification_image(lens, xt1, yt1, adis)
   @test dex1[1] ≈ 0.1800787586968404 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ 0.1800787586968404 atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ +0.01343917717947224 atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ +0.01343917717947224 atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ -0.16663958151736816 atol=1e-15 rtol=1e-15
   @test mag1 ≈ 1.0576039797648285 atol=1e-15 rtol=1e-15

   pot2 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   mag2 = Lenses.get_magnification_image(lens, xt2, yt2, adis)
   @test dex2[1] ≈ 0.3172256889321225 atol=1e-15 rtol=1e-15
   @test dex2[2] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac2[1] ≈ -0.11862002072694444 atol=1e-15 rtol=1e-15
   @test jac2[2] ≈ +0.3172256889321225 atol=1e-15 rtol=1e-15
   @test jac2[3] ≈ -0.0 atol=1e-15 rtol=1e-15
   @test mag2 ≈ 1.3093032302758327 atol=1e-15 rtol=1e-15

   # Lenstronomy gives different potential values. However, the relevant quantities is potential
   # difference, which is what I am going to test here.
   pot1 = adis * Lenses.get_potential(lens, xt1, yt1)
   pot2 = adis * Lenses.get_potential(lens, xt2, yt2)
   pot3 = adis * Lenses.get_potential(lens, 1.5, 2.5)
   @test pot1 - pot2 ≈ +0.118872955824665 atol=1e-15 rtol=1e-15
   @test pot1 - pot3 ≈ -0.26472744601185416 atol=1e-15 rtol=1e-15
   @test pot2 - pot3 ≈ -0.38360040183651917 atol=1e-15 rtol=1e-15

   potc = adis  * Lenses.get_potential(lens, [xt1, xt2], [yt1, yt2])
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