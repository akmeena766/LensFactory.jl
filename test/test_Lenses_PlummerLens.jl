#!!!!!!!!!!!!!! Testing AGAINST formulae and cross-checked against PyGRALE !!!!!!!!!!!!!!
@testset "Plummer lens" begin
   mass = 1.0E11 * MASS_SUN
   x_s = 0.1

   # Create a Plummer lens
   lens = Lenses.init_PlummerLens(D_d=Dol, mass=1E11*MASS_SUN, x_s=x_s)

   θE2 = (4.0 * CONST_G * mass / CONST_C^2) / Dol / ANGLE_ARCSEC^2

   # Test Einstein angle calculation
   θE = Lenses.PlummerLens.einstein_angle(D_d=Dol, D_ds=Dls, D_s=Dos, mass=mass, x_s=x_s)
   @test θE > 0.0
   @test θE ≈ sqrt(θE2 * adis - x_s^2) atol=1e-15 rtol=1e-15

   # Test lensing quantities
   pot1 = adis * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   @test pot1 ≈ adis * 0.5 * θE2 * log(0.1^2 + xt1^2 + yt1^2) atol=1e-15 rtol=1e-15
   @test dex1[1] ≈  adis * θE2 * (xt1) / (0.1^2 + xt1^2 + yt1^2) atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ adis * θE2 * (yt1) / (0.1^2 + xt1^2 + yt1^2) atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ +adis * θE2 * (0.1^2 + yt1^2 - xt1^2) / (0.1^2 + xt1^2 + yt1^2)^2 atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ +adis * θE2 * (0.1^2 + xt1^2 - yt1^2) / (0.1^2 + xt1^2 + yt1^2)^2 atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ -adis * θE2 * 2.0 * xt1 * yt1 / (0.1^2 + xt1^2 + yt1^2)^2 atol=1e-15 rtol=1e-15

   pot2 = adis  * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   @test pot2 ≈ adis * 0.5 * θE2 * log(0.1^2 + xt2^2 + yt2^2) atol=1e-15 rtol=1e-15
   @test dex2[1] ≈ adis * θE2 * (xt2) / (0.1^2 + xt2^2 + yt2^2) atol=1e-15 rtol=1e-15
   @test dex2[2] ≈ adis * θE2 * (yt2) / (0.1^2 + xt2^2 + yt2^2) atol=1e-15 rtol=1e-15
   @test jac2[1] ≈ +adis * θE2 * (0.1^2 + yt2^2 - xt2^2) / (0.1^2 + xt2^2 + yt2^2)^2 atol=1e-15 rtol=1e-15
   @test jac2[2] ≈ +adis * θE2 * (0.1^2 + xt2^2 - yt2^2) / (0.1^2 + xt2^2 + yt2^2)^2 atol=1e-15 rtol=1e-15
   @test jac2[3] ≈ -adis * θE2 * 2.0 * xt2 * yt2 / (0.1^2 + xt2^2 + yt2^2)^2 atol=1e-15 rtol=1e-15

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