@testset "Composite Lens" begin
   lens = Lenses.init_CompositeLens([(lens=:PointLens, D_d=1.0, x_c=0.0, y_c=0.0, mass=1.0)])
   
   # Composite lens constructor: size and parameters
   @test typeof(lens) <: Lenses.AbstractLens
   @test lens._lens_ == :CompositeLens
   @test length(lens._components_) == 1
   @test lens._components_[1].mass == 1

   # Composite lens constructor: Unknown lens error handeling
   @test_throws ArgumentError Lenses.init_CompositeLens([(lens=:BadLens, D_d=1.0, x_c=0.0, y_c=0.0, mass=1.0)])
end
