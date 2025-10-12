#!!!!!!!!!!!!!! Testing AGAINST formulae !!!!!!!!!!!!!!
@testset "Plummer lens" begin
   # Einstein angle
   mass = 1.0E11 * MASS_SUN
   θE = Lenses.PlummerLens.einstein_angle(D_d=Dol, D_ds=Dls, D_s=Dos, mass=mass, x_s=0.1)
   @test θE > 0.0
   @test θE ≈ sqrt((4.0 * CONST_G * mass / CONST_C^2) * (Dls / Dol / Dos) / ANGLE_ARCSEC^2 - 0.1^2) atol=1e-15 rtol=1e-15

   # Create a Plummer lens
   lens = Lenses.init_PlummerLens(D_d=Dol, mass=1E11*MASS_SUN, x_s=0.1)

   θE2 = (4.0 * CONST_G * mass / CONST_C^2) / Dol / ANGLE_ARCSEC^2

   pot = adis * Lenses.get_potential(lens, xt1, yt1)
   dex = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   @test pot ≈ adis * 0.5 * θE2 * log(0.1^2 + xt1^2 + yt1^2) atol=1e-15 rtol=1e-15
   @test dex[1] ≈ adis * θE2 * (xt1) / (0.1^2 + xt1^2 + yt1^2) atol=1e-15 rtol=1e-15
   @test dex[2] ≈ adis * θE2 * (yt1) / (0.1^2 + xt1^2 + yt1^2) atol=1e-15 rtol=1e-15
   @test jac[1] ≈ +adis * θE2 * (0.1^2 + yt1^2 - xt1^2) / (0.1^2 + xt1^2 + yt1^2)^2 atol=1e-15 rtol=1e-15
   @test jac[2] ≈ +adis * θE2 * (0.1^2 + xt1^2 - yt1^2) / (0.1^2 + xt1^2 + yt1^2)^2 atol=1e-15 rtol=1e-15
   @test jac[3] ≈ -adis * θE2 * 2.0 * xt1 * yt1 / (0.1^2 + xt1^2 + yt1^2)^2 atol=1e-15 rtol=1e-15

   pot = adis * Lenses.get_potential(lens, xt2, yt2)
   dex = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   @test pot ≈ adis * 0.5 * θE2 * log(0.1^2 + xt2^2 + yt2^2) atol=1e-15 rtol=1e-15
   @test dex[1] ≈ adis * θE2 * (xt2) / (0.1^2 + xt2^2 + yt2^2) atol=1e-15 rtol=1e-15
   @test dex[2] ≈ adis * θE2 * (yt2) / (0.1^2 + xt2^2 + yt2^2) atol=1e-15 rtol=1e-15
   @test jac[1] ≈ +adis * θE2 * (0.1^2 + yt2^2 - xt2^2) / (0.1^2 + xt2^2 + yt2^2)^2 atol=1e-15 rtol=1e-15
   @test jac[2] ≈ +adis * θE2 * (0.1^2 + xt2^2 - yt2^2) / (0.1^2 + xt2^2 + yt2^2)^2 atol=1e-15 rtol=1e-15
   @test jac[3] ≈ -adis * θE2 * 2.0 * xt2 * yt2 / (0.1^2 + xt2^2 + yt2^2)^2 atol=1e-15 rtol=1e-15
end