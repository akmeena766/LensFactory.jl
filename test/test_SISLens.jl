#!!!!!!!!!!!!!! Testing AGAINST LENSTRONOMY !!!!!!!!!!!!!!
@testset "SIS lens" begin
   # Einstein angle
   D_d, D_ds, D_s, v_d = 1.0, 1.0, 1.0, 200E3
   θE = Lenses.SISLens.einstein_angle(D_ds=D_ds, D_s=D_s, v_d=v_d)
   @test θE > 0.0
   @test θE ≈ 4.0 * π * (v_d^2 / CONST_C^2) / ANGLE_ARCSEC atol=1e-15 rtol=1e-15

   # Create a SIS lens
   lens = Lenses.init_SISLens(v_d=200E3)

   pot = adis * Lenses.get_potential(lens, xt1, yt1)
   dex = adis .* Lenses.get_deflection(lens, xt1, yt1)
   jac = adis .* Lenses.get_jacobian(lens, xt1, yt1)
   mag = Lenses.get_magnification_image(lens, xt1, yt1, adis)

   @test pot ≈ 0.9253665947775005 atol=1e-15 rtol=1e-15
   @test dex[1] ≈ 0.4626832973887502 atol=1e-15 rtol=1e-15
   @test dex[2] ≈ 0.4626832973887502 atol=1e-15 rtol=1e-15
   @test jac[1] ≈ 0.2313416486943751 atol=1e-15 rtol=1e-15
   @test jac[2] ≈ 0.2313416486943751 atol=1e-15 rtol=1e-15
   @test jac[3] ≈ -0.2313416486943751 atol=1e-15 rtol=1e-15 
   @test mag ≈ 1.8610997855458498 atol=1e-15 rtol=1e-15

   pot = adis * Lenses.get_potential(lens, xt2, yt2)
   dex = adis .* Lenses.get_deflection(lens, xt2, yt2)
   jac = adis .* Lenses.get_jacobian(lens, xt2, yt2)
   mag = Lenses.get_magnification_image(lens, xt2, yt2, adis)

   @test pot ≈ 0.6543329942506746 atol=1e-15 rtol=1e-15
   @test dex[1] ≈ 0.6543329942506746 atol=1e-15 rtol=1e-15
   @test dex[2] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac[1] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test jac[2] ≈ 0.6543329942506746 atol=1e-15 rtol=1e-15
   @test jac[3] ≈ 0.0 atol=1e-15 rtol=1e-15
   @test mag ≈ 2.892957625019007 atol=1e-15 rtol=1e-15
end
