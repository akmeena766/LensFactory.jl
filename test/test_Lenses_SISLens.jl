#!!!!!!!!!!!!!! Testing AGAINST LENSTRONOMY !!!!!!!!!!!!!!
@testset "SIS lens" begin
   # Einstein angle
   D_d, D_ds, D_s, v_d = 1.0, 1.0, 1.0, 200E3
   θE = Lenses.SISLens.einstein_angle(D_ds=D_ds, D_s=D_s, v_d=v_d)
   @test θE > 0.0
   @test θE ≈ 4.0 * π * (v_d^2 / CONST_C^2) / ANGLE_ARCSEC atol=1e-15 rtol=1e-15

   # Create a SIS lens
   lens = Lenses.init_SISLens(v_d=200E3)

   pot1 = adis  * Lenses.get_potential(lens, xt1, yt1)
   dex1 = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac1 = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   mag1 = Lenses.get_magnification_image(lens, xt1, yt1, adis)
   @test pot1 ≈ 0.9253665947775005 atol=1e-15 rtol=1e-15
   @test dex1[1] ≈ 0.4626832973887502 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ 0.4626832973887502 atol=1e-15 rtol=1e-15
   @test jac1[1] ≈ 0.2313416486943751 atol=1e-15 rtol=1e-15
   @test jac1[2] ≈ 0.2313416486943751 atol=1e-15 rtol=1e-15
   @test jac1[3] ≈ -0.2313416486943751 atol=1e-15 rtol=1e-15 
   @test mag1 ≈ 1.8610997855458498 atol=1e-15 rtol=1e-15

   pot2 = adis  * Lenses.get_potential(lens, xt2, yt2)
   dex2 = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac2 = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   mag2 = Lenses.get_magnification_image(lens, xt2, yt2, adis)
   @test pot2 ≈ 0.6543329942506746 atol=1e-15 rtol=1e-15
   @test dex2[1] ≈ 0.6543329942506746 atol=1e-15 rtol=1e-15
   @test dex2[2] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac2[1] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac2[2] ≈ 0.6543329942506746 atol=1e-15 rtol=1e-15
   @test jac2[3] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test mag2 ≈ 2.892957625019007 atol=1e-15 rtol=1e-15

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
