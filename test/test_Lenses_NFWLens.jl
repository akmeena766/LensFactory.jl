#!!!!!!!!!!!!!! Testing/cross-checked AGAINST PyGRALE/Glafic !!!!!!!!!!!!!!
# cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
# z_d, z_s = 0.5, 1.5
# D_d = cosm.getAngularDiameterDistance(z_d)
# D_s = cosm.getAngularDiameterDistance(z_s)
# D_ds = cosm.getAngularDiameterDistance(z_s, z_d)
# adis = D_ds/D_s

# lens = lenses.NFWLens(D_d, { "rho_s": 2.0846284842174704e-22, "theta_s": 2.185063820164341 * cv.ANGLE_ARCSEC})
# pot = lens.getProjectedPotential(D_s, D_ds, [1.0 * cv.ANGLE_ARCSEC, 1.0 * cv.ANGLE_ARCSEC]) / cv.ANGLE_ARCSEC**2
# dex = lens.getAlphaVector(np.array([1.0 * cv.ANGLE_ARCSEC, 1.0 * cv.ANGLE_ARCSEC])) / cv.ANGLE_ARCSEC
# jac = lens.getAlphaVectorDerivatives(np.array([1.0 * cv.ANGLE_ARCSEC, 1.0 * cv.ANGLE_ARCSEC]))
# mag = 1 / lens.getInverseMagnification(D_s, D_ds, np.array([1.0 * cv.ANGLE_ARCSEC, 1.0 * cv.ANGLE_ARCSEC]))
# print("Potential at (x, y): ", list(adis * pot * 6.6743e-11 / cv.CONST_G))
# print("Deflection at (x, y): ", list(adis * dex * 6.6743e-11 / cv.CONST_G))
# print("Deformation at (x, y): ", list(adis * jac * 6.6743e-11 / cv.CONST_G))
# print("Magnfication at (x, y):", list(mag))

@testset "NFW lens" begin
   mass = 1E11 * MASS_SUN
   
   # Get NFW lens parameters
   param = Lenses.parameter_NFWLens(cosmology=cosmo, z_d=zl, mass=mass, c=6)
   @test_throws ArgumentError Lenses.parameter_NFWLens(cosmology=cosmo, z_d=zl, mass=mass)

   lens = Lenses.init_NFWLens(cosmo, zl; mass=mass, c=6)

   pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   @test dex1[1] ≈ 0.031306618316540795 atol=1e-4 rtol=1e-4
   @test dex1[2] ≈ 0.031306618316540795 atol=1e-4 rtol=1e-4
   @test jac1[1] ≈ +0.019069805741985472 atol=1e-4 rtol=1e-4
   @test jac1[2] ≈ +0.019069805741985472 atol=1e-4 rtol=1e-4
   @test jac1[3] ≈ -0.012236812574555332 atol=1e-4 rtol=1e-4

   pot2 = adis  * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   @test dex2[1] ≈ 0.04035349424899914 atol=1e-4 rtol=1e-4
   @test dex2[2] ≈ 0.0 atol=1e-4 rtol=1e-4
   @test jac2[1] ≈ +0.012724243945134012 atol=1e-4 rtol=1e-4
   @test jac2[2] ≈ +0.04035349424899913 atol=1e-4 rtol=1e-4
   @test jac2[3] ≈ +0.0 atol=1e-4 rtol=1e-4

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

   pot1 = adis  * Lenses.get_potential(lens, 1.0, -1.0)
   pot2 = adis  * Lenses.get_potential(lens, 2.0, +0.0)
   pot3 = adis  * Lenses.get_potential(lens, 3.0, +1.0)
   @test pot1 - pot2 ≈ -0.026823145088490218 atol=1e-4 rtol=1e-4
   @test pot2 - pot3 ≈ -0.05538925777883399 atol=1e-4 rtol=1e-4
   @test pot3 - pot1 ≈ +0.0822124028673242 atol=1e-4 rtol=1e-4
end