@testset "Lenses" begin
   @testset "Meshgrid" begin
      # 2D grid constructor: grid/pixel size <=0 error handeling
      @test_throws ArgumentError Lenses.get_meshgrid(1, 1, -1)
      @test_throws ArgumentError Lenses.get_meshgrid(1, 0, +1)

      # 2D grid constructor: size 
      xp, yp, pp = 1, 2, 0.3
      xg, yg = Lenses.get_meshgrid(xp, yp, pp)
      nx, ny = size(xg)
      @test size(xg) == size(yg)
      @test size(xg, 1) == size(collect(-xp:pp:xp), 1)
      @test size(xg, 2) == size(collect(-yp:pp:yp), 1)
   end

   @testset "Lens type" begin
      lens = Lenses.init_CompositeLens([(lens=:PointLens, D_d=1.0, x_c=0.0, y_c=0.0, mass=1.0)])
      
      # Composite lens constructor: size and parameters
      @test typeof(lens) <: Lenses.AbstractLens
      @test lens._lens_ == :CompositeLens
      @test length(lens._components_) == 1
      @test lens._components_[1].mass == 1

      # Composite lens constructor: Unknown lens error handeling
      @test_throws ArgumentError Lenses.init_CompositeLens([(lens=:BadLens, D_d=1.0, x_c=0.0, y_c=0.0, mass=1.0)])

      # Multi-plane lens constructor
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

   @testset "Σ_cr" begin
      D_d, D_ds, D_s = 1.0, 1.0, 2.0
      Σcr = Lenses.get_critical_density(D_d=D_d, D_ds=D_ds, D_s=D_s)
      
      # Σ_cr > 0 
      @test Σcr > 0.0

      # Σ_cr using distance ratio
      @test Σcr ≈ Lenses.get_critical_density(D_d=D_d, adis=0.5)

      # Unit conversion: kg/m^2 --> msun/pc^2
      @test Σcr * (DIST_PC^2 / MASS_SUN) ≈ Lenses.get_critical_density(D_d=D_d, D_ds=D_ds, D_s=D_s; unit=:msun_pc2)

      # Unit conversion: kg/m^2 --> msun/arcsec^2
      @test Σcr * (D_d^2 * ANGLE_ARCSEC^2 / MASS_SUN) ≈ Lenses.get_critical_density(D_d=D_d, D_ds=D_ds, D_s=D_s; unit=:msun_arcsec2)

      # Unknown unit error handeling
      @test_throws ArgumentError Lenses.get_critical_density(D_d=D_d, D_ds=D_ds, D_s=D_s, unit=:invalid)
   end

   
end