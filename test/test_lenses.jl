@testset "Lenses" begin
   @testset "Lens type" begin
      lens = Lenses.init_CompositeLens([(lens=:PointLens, D_d=1.0, x_c=0.0, y_c=0.0, mass=1.0)])
      
      # Test on the type and values
      @test typeof(lens) <: Lenses.AbstractLens
      @test lens._lens_ == :CompositeLens
      @test length(lens._components_) == 1
      @test lens._components_[1].mass == 1

      # Bad lens throws error
      @test_throws ArgumentError Lenses.init_CompositeLens([(lens=:BadLens, D_d=1.0, x_c=0.0, y_c=0.0, mass=1.0)])
   end
end