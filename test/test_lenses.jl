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

      # Multiplane lens test size
      multi_lens = [(lens=:PointLens, z_d=0.4, D_d=1, x_c= -2, y_c= 0, mass=1),
                    (lens=:PointLens, z_d=0.5, D_d=1, x_c= +2, y_c= 0, mass=1),
                    (lens=:PointLens, z_d=0.5, D_d=1, x_c= +2, y_c= 0, mass=1),
                    (lens=:PointLens, z_d=0.2, D_d=1, x_c= +2, y_c= 0, mass=1)]

      lens = Lenses.init_MultiPlaneLens(multi_lens)
      @test lens._lens_ == :MultiPlaneLens
      @test lens.z_d == [0.2, 0.4, 0.5]
      @test length(lens._plane_) == lens.n_p
      @test lens._plane_[1]._lens_ == :CompositeLens
      @test lens._plane_[2]._lens_ == :CompositeLens
      @test length(lens._plane_[1]._components_) == 1
      @test length(lens._plane_[3]._components_) == 2
   end
end