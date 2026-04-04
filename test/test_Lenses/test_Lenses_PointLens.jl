#!!!!!!!!!!!!!! Testing AGAINST LENSTRONOMY !!!!!!!!!!!!!!
# lensModel = LensModel(lens_model_list=['POINT_MASS'])
# point = {'center_x': 0, 'center_y': 0, 'theta_E': 0.6057038104415533}
# tx, ty = [1.0, 1.0], [1.0, 0.0]
# pot = lensModel.potential(tx, ty, [point])
# dex = lensModel.alpha(tx, ty, [point])
# hes = lensModel.hessian(tx, ty, [point])
# mag = lensModel.magnification(tx, ty, [point])

@testset "Point lens" begin
   # Einstein angle
   D_d, D_ds, D_s, mass = 1.0, 1.0, 1.0, 1.0
   θE = Lenses.PointLens.einstein_angle(D_d=D_d, D_ds=D_ds, D_s=D_s, mass=mass)
   @test θE > 0.0
   @test θE ≈ sqrt(4.0 *CONST_G * mass * MASS_SUN / CONST_C^2 * (D_ds / D_d / D_s) ) / ANGLE_ARCSEC atol=1e-15 rtol=1e-15

   # Create a point lens
   lens = Lenses.init_PointLens(D_d=Dol, mass=1E11)

   pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   mag1 = Lenses.get_magnification_image(lens, xt1, yt1, adis)
   @test pot1 ≈ 0.12714991581219895 atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ 0.18343855299170855 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ 0.18343855299170855 atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ -0.18343855299170858 atol=1e-15 rtol=1e-15 
   @test mag1 ≈ 1.0348214336131885 atol=1e-15 rtol=1e-15

   pot2 = adis * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt1, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   mag2 = Lenses.get_magnification_image(lens, xt2, yt2, adis)
   @test pot2 ≈ 0.0 atol=1e-15 rtol=1e-15
   @test dex2[1] ≈ 0.36687710598341716 atol=1e-15 rtol=1e-15
   @test dex2[2] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac2[1] ≈ -0.36687710598341716 atol=1e-15 rtol=1e-15
   @test jac2[2] ≈ +0.36687710598341716 atol=1e-15 rtol=1e-15
   @test jac2[3] ≈ 0.0 atol=1e-15 rtol=1e-15 
   @test mag2 ≈ 1.155533424947028 atol=1e-15 rtol=1e-15

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
end