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

   #!!!!!!!!!!!!!! Testing individual lens models (AGAINST LENSTRONOMY/PyGRALE) !!!!!!!!!!!!!!
   
   # Cosmology   
   cosmo = Cosmology.init_cosmology()
   
   # Lens and source redshifts
   zl = 0.5
   zs = 1.5

   # ADDs and distance ratio
   Dol = Cosmology.angular_diameter_distance(cosmo, 0., zl)
   Dls = Cosmology.angular_diameter_distance(cosmo, zl, zs)
   Dos = Cosmology.angular_diameter_distance(cosmo, 0., zs)
   adis = Dls/Dos

   # Points (in image plane) to evaluate lensing quantities
   xt1, yt1 = 1.0, 1.0
   xt2, yt2 = 1.0, 0.0
   
   @testset "Point lens" begin
      # Einstein angle
      D_d, D_ds, D_s, mass = 1.0, 1.0, 1.0, 1.0
      θE = Lenses.PointLens.einstein_angle(D_d=D_d, D_ds=D_ds, D_s=D_s, mass=mass)
      @test θE > 0.0
      @test θE ≈ sqrt(4.0 *CONST_G * mass / CONST_C^2 * (D_ds / D_d / D_s) ) / ANGLE_ARCSEC atol=1e-15 rtol=1e-15

      # Create a point lens
      lens = Lenses.init_PointLens(D_d=Dol, mass=1E11*MASS_SUN)

      pot = adis * Lenses.get_potential(lens, xt1, yt1)
      dex = adis .* Lenses.get_deflection(lens, xt1, yt1)
      jac = adis .* Lenses.get_jacobian(lens, xt1, yt1)
      mag = Lenses.get_magnification_image(lens, xt1, yt1, adis)

      @test pot ≈ 0.12714991581219895 atol=1e-15 rtol=1e-15
      @test dex[1] ≈ 0.18343855299170855 atol=1e-15 rtol=1e-15
      @test dex[2] ≈ 0.18343855299170855 atol=1e-15 rtol=1e-15
      @test jac[1] ≈ 0.0 atol=1e-15 rtol=1e-15
      @test jac[2] ≈ 0.0 atol=1e-15 rtol=1e-15
      @test jac[3] ≈ -0.18343855299170858 atol=1e-15 rtol=1e-15 
      @test mag ≈ 1.0348214336131885 atol=1e-15 rtol=1e-15

      pot = adis * Lenses.get_potential(lens, xt2, yt2)
      dex = adis .* Lenses.get_deflection(lens, xt1, yt2)
      jac = adis .* Lenses.get_jacobian(lens, xt2, yt2)
      mag = Lenses.get_magnification_image(lens, xt2, yt2, adis)

      @test pot ≈ 0.0 atol=1e-15 rtol=1e-15
      @test dex[1] ≈ 0.36687710598341716 atol=1e-15 rtol=1e-15
      @test dex[2] ≈ 0.0 atol=1e-15 rtol=1e-15
      @test jac[1] ≈ -0.36687710598341716 atol=1e-15 rtol=1e-15
      @test jac[2] ≈ +0.36687710598341716 atol=1e-15 rtol=1e-15
      @test jac[3] ≈ 0.0 atol=1e-15 rtol=1e-15 
      @test mag ≈ 1.155533424947028 atol=1e-15 rtol=1e-15
   end

   @testset "Plummer lens" begin
      
   end

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
end