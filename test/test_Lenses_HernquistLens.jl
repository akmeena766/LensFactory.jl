#!!!!!!!!!!!!!! Testing AGAINST Lenstronomy + Glafic !!!!!!!!!!!!!!
# from lenstronomy.Cosmo.lens_cosmo import LensCosmo
# from astropy.cosmology import FlatLambdaCDM
# import lenstronomy.Util.constants as const

# cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Ob0=0.05)
# lens_cosmo = LensCosmo(z_lens=0.5, z_source=1.5, cosmo=cosmo)
# print("rs = 0.3 in arcseconds:", 0.3*const.arcsec*lens_cosmo.dd)
# sigma0, rs_angle = lens_cosmo.hernquist_phys2angular(mass=10**11, rs=0.0018312628608788694)
# m_tot, rs = lens_cosmo.hernquist_angular2phys(sigma0=sigma0, rs_angle=rs_angle)

# lensModel = LensModel(lens_model_list=['HERNQUIST'])
# hernquist = {'center_x': 0, 'center_y': 0, "Rs": rs_angle, "sigma0": sigma0}
# tx, ty = [1.0, 1.0], [1.0, 0.0]
# pot = lensModel.potential(tx, ty, [hernquist])
# dex = lensModel.alpha(tx, ty, [hernquist])
# hes = lensModel.hessian(tx, ty, [hernquist])
# mag = lensModel.magnification(tx, ty, [hernquist])

@testset "Hernquist lens" begin
   x_s = 0.3
   mass = 1E11 * MASS_SUN
   lens = Lenses.init_HernquistLens(D_d=Dol, mass=mass, x_s=x_s)

   pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   mag1 = Lenses.get_magnification_image(lens, xt1, yt1, adis)
   @test pot1 ≈ 0.4227508608283786 atol=1e-3 rtol=1e-3
   @test dex1[1] ≈ 0.13553713 atol=1e-3 rtol=1e-3
   @test dex1[2] ≈ 0.13553713 atol=1e-3 rtol=1e-3
   @test jac1[1] ≈ +0.01871943 atol=1e-3 rtol=1e-3
   @test jac1[2] ≈ +0.01871943 atol=1e-3 rtol=1e-3
   @test jac1[3] ≈ -0.1168177 atol=1e-3 rtol=1e-3
   @test mag1 ≈ 1.0534464026800108 atol=1e-3 rtol=1e-3

   pot2 = adis  * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   mag2 = Lenses.get_magnification_image(lens, xt2, yt2, adis)
   @test pot2 ≈ 0.33358287248183577 atol=1e-3 rtol=1e-3
   @test dex2[1] ≈ 0.24270179 atol=1e-3 rtol=1e-3
   @test dex2[2] ≈ 0.0 atol=1e-3 rtol=1e-3
   @test jac2[1] ≈ -0.1541402 atol=1e-3 rtol=1e-3
   @test jac2[2] ≈ +0.24270179 atol=1e-3 rtol=1e-3
   @test jac2[3] ≈ +0.0 atol=1e-3 rtol=1e-3
   @test mag2 ≈ 1.144127692401055 atol=1e-3 rtol=1e-3

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