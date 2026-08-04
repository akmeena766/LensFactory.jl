@testset "Composite Lens" begin
   lens = Lenses.init_CompositeLens([(lens=:PointLens, D_d=1.0, x_c=0.0, y_c=0.0, mass=1.0)])
   
   # Composite lens constructor: size and parameters
   @test typeof(lens) <: Lenses.AbstractLens
   @test lens._lens_ == :CompositeLens
   @test length(lens._components_) == 1
   @test lens._components_[1].mass == 1

   # Composite lens constructor: Unknown lens error handeling
   @test_throws ArgumentError Lenses.init_CompositeLens([(lens=:BadLens, D_d=1.0, x_c=0.0, y_c=0.0, mass=1.0)])

   # ---- A Pixel component inside a composite equals the sum of the parts ---------------
   pt    = Lenses.init_PointLens(D_d=1.0, x_c=0.3, y_c=0.4, mass=1.0)
   px    = Lenses.init_PixelLens(x_c=0.1, y_c=-0.2, kappa=0.8, pixel_size=0.7)
   mixed = Lenses.init_CompositeLens([(lens=:PointLens, D_d=1.0, x_c=0.3, y_c=0.4, mass=1.0),
                                      (lens=:PixelLens, x_c=0.1, y_c=-0.2, kappa=0.8, pixel_size=0.7)])
   @test length(mixed._components_) == 2
   @test mixed._components_[2]._lens_ == :PixelLens

   @test Lenses.get_potential(mixed, xt1, yt1) ≈
         Lenses.get_potential(pt, xt1, yt1) + Lenses.get_potential(px, xt1, yt1) atol=1e-12
   let dm = Lenses.get_deflection(mixed, xt1, yt1),
       dp = Lenses.get_deflection(pt, xt1, yt1),
       dx = Lenses.get_deflection(px, xt1, yt1)
      @test dm[1] ≈ dp[1] + dx[1] atol=1e-12
      @test dm[2] ≈ dp[2] + dx[2] atol=1e-12
   end
   let jm = Lenses.get_jacobian(mixed, xt1, yt1),
       jp = Lenses.get_jacobian(pt, xt1, yt1),
       jx = Lenses.get_jacobian(px, xt1, yt1)
      @test jm[1] ≈ jp[1] + jx[1] atol=1e-12
      @test jm[2] ≈ jp[2] + jx[2] atol=1e-12
      @test jm[3] ≈ jp[3] + jx[3] atol=1e-12
   end

   # ---- A MultiPixelLens component is built correctly and matches a standalone one -----
   mpx  = Lenses.init_MultiPixelLens(x_c=[0.0, 0.5], y_c=[0.0, 0.5], kappa=[1.0, 0.7], pixel_size=[1.0, 0.6])
   cmpx = Lenses.init_CompositeLens([(lens=:MultiPixelLens, x_c=[0.0, 0.5], y_c=[0.0, 0.5],
                                      kappa=[1.0, 0.7], pixel_size=[1.0, 0.6])])
   @test cmpx._components_[1]._lens_ == :MultiPixelLens
   @test Lenses.get_potential(cmpx, xt1, yt1)   ≈ Lenses.get_potential(mpx, xt1, yt1)   atol=1e-13
   @test Lenses.get_deflection(cmpx, xt1, yt1)[1] ≈ Lenses.get_deflection(mpx, xt1, yt1)[1] atol=1e-13
   @test Lenses.get_jacobian(cmpx, xt1, yt1)[3] ≈ Lenses.get_jacobian(mpx, xt1, yt1)[3] atol=1e-13
end
