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

#!!!!!!!!!!!!!! Testing individual lens models !!!!!!!!!!!!!!
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

include("./test_Lenses_PointLens.jl")
include("./test_Lenses_PlummerLens.jl")
include("./test_Lenses_SISLens.jl")
include("./test_Lenses_NSISPLens.jl")
include("./test_Lenses_NSISMDLens.jl")
include("./test_Lenses_GaussianLens.jl")
include("./test_Lenses_ExternalEffects.jl")
include("./test_Lenses_HernquistLens.jl")
include("./test_Lenses_CompositeLens.jl")
include("./test_Lenses_MultiPlaneLens.jl")