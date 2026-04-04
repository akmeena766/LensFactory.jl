#!!!!!!!!!!!!!! Testing AGAINST Lenstronomy + Glafic !!!!!!!!!!!!!!
# from lenstronomy.Cosmo.lens_cosmo import LensCosmo
# from astropy.cosmology import FlatLambdaCDM

# cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Ob0=0.0)
# lens_cosmo = LensCosmo(z_lens=0.5, z_source=1.5, cosmo=cosmo)

# R_sersic, n_sersic = 1.0, 4.0
# k_eff = lens_cosmo.sersic_m_star2k_eff(m_star=10**11, R_sersic=R_sersic, n_sersic=n_sersic)

# lensModel = LensModel(lens_model_list=['SERSIC'])
# sersic = {'center_x': 0, 'center_y': 0, 'k_eff': k_eff, 'R_sersic': R_sersic, "n_sersic": n_sersic}
# tx, ty = [1.0], [1.0]
# pot = lensModel.potential(tx, ty, [sersic])
# dex = lensModel.alpha(tx, ty, [sersic])
# hes = lensModel.hessian(tx, ty, [sersic])
# mag = lensModel.magnification(tx, ty, [sersic])
# mag = lensModel.magnification(tx, ty, [sersic])

@testset "Sersic lens" begin
   n = 4
   x_e = 1.0
   mass = 1E11
   lens = Lenses.init_SersicLens(D_d=Dol, mass=mass, x_e=x_e)

   pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   mag1 = Lenses.get_magnification_image(lens, xt1, yt1, adis)
   @test pot1 ≈ 0.2977881520594202 atol=1e-3 rtol=1e-3
   @test dex1[1] ≈ 0.10946562 atol=1e-3 rtol=1e-3
   @test dex1[2] ≈ 0.10946562 atol=1e-3 rtol=1e-3
   @test jac1[1] ≈ +0.02540788 atol=1e-3 rtol=1e-3
   @test jac1[2] ≈ +0.02540788 atol=1e-3 rtol=1e-3
   @test jac1[3] ≈ -0.08405774 atol=1e-3 rtol=1e-3
   @test mag1 ≈ 1.0607107322260219 atol=1e-3 rtol=1e-3

   pot2 = adis  * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   mag2 = Lenses.get_magnification_image(lens, xt2, yt2, adis)
   @test pot2 ≈ 0.22804884511775525 atol=1e-3 rtol=1e-3
   @test dex2[1] ≈ 0.18351348 atol=1e-3 rtol=1e-3
   @test dex2[2] ≈ 0.0 atol=1e-3 rtol=1e-3
   @test jac2[1] ≈ -0.08177834 atol=1e-3 rtol=1e-3
   @test jac2[2] ≈ +0.18351348 atol=1e-3 rtol=1e-3
   @test jac2[3] ≈ +0.0 atol=1e-3 rtol=1e-3
   @test mag2 ≈ 1.1321727534761221 atol=1e-3 rtol=1e-3

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

   #! For n < 0.36, I compared the values with Glafic and got same value upto third decimal point.
   #! Hence, I think we are good to go!!
end