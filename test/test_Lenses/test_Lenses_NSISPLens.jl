#!!!!!!!!!!!!!! Testing AGAINST formulae !!!!!!!!!!!!!!
@testset "NSISP lens" begin
   # Einstein angle
   v_d = 200.0
   x_s = 0.1

   # Create an NSISP lens
   lens = Lenses.init_NSISPLens(v_d=v_d, x_s=x_s)

   constant = 4π * (v_d * 1.0E3/CONST_C)^2 * adis / ANGLE_ARCSEC

   # Test Einstein angle
   θE = Lenses.NSISPLens.einstein_angle(D_ds=Dls, D_s=Dos, v_d=v_d, x_s=x_s)
   @test θE > 0.0
   @test θE ≈ sqrt(constant^2 - x_s^2) atol=1e-15 rtol=1e-15

   # Test lensing quantities
   pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   @test pot1    ≈ constant * sqrt(x_s^2 + xt1^2 + yt1^2) atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ constant * xt1 / sqrt(x_s^2 + xt1^2 + yt1^2) atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ constant * yt1 / sqrt(x_s^2 + xt1^2 + yt1^2) atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ +constant * (x_s^2 + yt1^2) / (x_s^2 + xt1^2 + yt1^2)^(3/2) atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ +constant * (x_s^2 + xt1^2) / (x_s^2 + xt1^2 + yt1^2)^(3/2) atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ -constant * xt1 * yt1 / (x_s^2 + xt1^2 + yt1^2)^(3/2) atol=1e-15 rtol=1e-15

   pot2 = adis  * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   @test pot2    ≈ constant * sqrt(x_s^2 + xt2^2 + yt2^2) atol=1e-15 rtol=1e-15
   @test dex2[1] ≈ constant * xt2 / sqrt(x_s^2 + xt2^2 + yt2^2) atol=1e-15 rtol=1e-15
   @test dex2[2] ≈ constant * yt2 / sqrt(x_s^2 + xt2^2 + yt2^2) atol=1e-15 rtol=1e-15
   @test jac2[1] ≈ +constant * (x_s^2 + yt2^2) / (x_s^2 + xt2^2 + yt2^2)^(3/2) atol=1e-15 rtol=1e-15
   @test jac2[2] ≈ +constant * (x_s^2 + xt2^2) / (x_s^2 + xt2^2 + yt2^2)^(3/2) atol=1e-15 rtol=1e-15
   @test jac2[3] ≈ -constant * xt2 * yt2 / (x_s^2 + xt2^2 + yt2^2)^(3/2) atol=1e-15 rtol=1e-15

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