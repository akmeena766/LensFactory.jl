#!!!!!!!!!!!!!! Testing AGAINST formulae and cross-checked against PyGRALE !!!!!!!!!!!!!!
# cosm = cosmology.Cosmology(0.7, 0.3, 0, 0.7)
# z_d, z_s = 0.5, 1.5
# D_d = cosm.getAngularDiameterDistance(z_d)
# D_s = cosm.getAngularDiameterDistance(z_s)
# D_ds = cosm.getAngularDiameterDistance(z_s, z_d)
# adis = D_ds/D_s

# lens = lenses.NSISLens(D_d, { "velocityDispersion": 200E3, "coreRadius": 0.3 * cv.ANGLE_ARCSEC })
# dex = lens.getAlphaVector(np.array([1.0 * cv.ANGLE_ARCSEC, 0.0 * cv.ANGLE_ARCSEC])) / cv.ANGLE_ARCSEC
# jac = lens.getAlphaVectorDerivatives(np.array([1.0 * cv.ANGLE_ARCSEC, 0.0 * cv.ANGLE_ARCSEC]))


@testset "NSISMD lens" begin
   # lens properties
   v_d = 200E3
   x_s = 0.3

   # Create an NSISMD lens
   lens = Lenses.init_NSISMDLens(v_d=v_d, x_s=x_s)

   constant = 4π * (v_d/CONST_C)^2 * adis / ANGLE_ARCSEC

   # Test Einstein angle
   θE = Lenses.NSISMDLens.einstein_angle(D_ds=Dls, D_s=Dos, v_d=v_d, x_s=x_s)
   @test θE > 0
   @test θE ≈ sqrt(constant^2 - 2.0 * x_s * constant) atol=1e-15 rtol=1e-15

   pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   @test pot1 ≈ constant * (sqrt(x_s^2 + xt1^2 + yt1^2) - x_s * log(sqrt(x_s^2 + xt1^2 + yt1^2) + x_s)) atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ 0.3748291690042468 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ 0.3748291690042468 atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ +0.22630579815399432 atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ +0.22630579815399432 atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ -0.14852337085025255 atol=1e-15 rtol=1e-15

   pot2 = adis  * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   @test pot2 ≈ constant * (sqrt(x_s^2 + xt2^2 + yt2^2) - x_s * log(sqrt(x_s^2 + xt2^2 + yt2^2) + x_s)) atol=1e-15 rtol=1e-15
   @test dex2[1] ≈ 0.48684380361182233 atol=1e-15 rtol=1e-15
   @test dex2[2] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac2[1] ≈ +0.1398935375689341 atol=1e-15 rtol=1e-15
   @test jac2[2] ≈ +0.48684380361182233 atol=1e-15 rtol=1e-15
   @test jac2[3] ≈ +0.0 atol=1e-15 rtol=1e-15

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
