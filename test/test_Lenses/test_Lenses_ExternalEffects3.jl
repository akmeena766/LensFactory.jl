#!!!!!!!!!!!!!! Testing AGAINST formulae !!!!!!!!!!!!!!
@testset "ExternalEffects3 lens" begin
   delta = 0.2
   angle = 90.0
   lens = Lenses.init_ExternalEffects3(delta=delta, angle=angle)
   pot1 = Lenses.get_potential(lens, xt1, yt1)
   dex1 = Lenses.get_deflection(lens, xt1, yt1)

   @test dex1[1] ≈ 0.4 atol=1e-15 rtol=1e-15
   @test dex1[2] ≈ 0.2 atol=1e-15 rtol=1e-15


end