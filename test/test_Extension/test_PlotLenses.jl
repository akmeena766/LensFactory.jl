# --------------------------------------------------------------------------------------------------
# Tests for the single-plane plotting API implemented in `ext/PlotLenses.jl`.
#
# The `Lenses.plot_*` functions are empty stubs in `src/` that only get implementations once the
# `PlotExt` extension loads (triggered by Makie + KernelDensity + LaTeXStrings, loaded in
# `runtests.jl`). These tests check that the extension initialized, that every entry point gained a
# method, and that each one runs and returns valid Makie `Figure`/`Axis` objects.
# --------------------------------------------------------------------------------------------------

@testset "PlotLenses (single-plane plotting)" begin

   # ---- Extension initialization -----------------------------------------------------------------
   # None of the single-plane plots can work unless the extension actually loaded.
   ext = Base.get_extension(LensFactory, :PlotExt)
   @test ext !== nothing
   @test ext isa Module

   # ---- Wiring: every plot_* entry point gained a method from the extension -----------------------
   for f in (Lenses.plot_sky, Lenses.plot_image_plane, Lenses.plot_surface_density,
             Lenses.plot_magnification_map, Lenses.plot_magnification_profile)
      @test !isempty(methods(f))
   end

   # ---- Fixtures ---------------------------------------------------------------------------------
   cosmo    = Cosmology.init_cosmology()
   z_l, z_s = 0.3, 1.5
   D_d      = Cosmology.angular_diameter_distance(cosmo, 0.0, z_l)
   D_s      = Cosmology.angular_diameter_distance(cosmo, 0.0, z_s)
   D_ds     = Cosmology.angular_diameter_distance(cosmo, z_l, z_s)
   adis     = D_ds / D_s

   # Lens centre is offset by half a pixel so that no grid point ever hits the singularity.
   lens   = Lenses.init_PointLens(D_d = D_d, x_c = 0.05, y_c = 0.05, mass = 1.0e12)
   θx, θy = Lenses.get_meshgrid(3.0, 3.0, 0.1)

   # ---- Smoke tests: entry points build valid Figure/Axis objects --------------------------------
   @testset "plot_sky initializes a figure" begin
      fig, ax = Lenses.plot_sky(3.0, 3.0)
      @test fig isa Makie.Figure
      @test ax  isa Makie.Axis

      # Two-panel layout returns a pair of axes.
      fig2, axes = Lenses.plot_sky(3.0, 3.0; two_panel = true)
      @test fig2 isa Makie.Figure
      @test length(axes) == 2
      @test all(a -> a isa Makie.Axis, axes)
   end

   @testset "plot_image_plane initializes a figure" begin
      # Default keywords also exercise the critical-curve / caustic code paths.
      fig, ax = Lenses.plot_image_plane(lens, θx, θy, adis)
      @test fig isa Makie.Figure
      @test ax  isa Makie.Axis

      # Two-panel layout with a point source.
      fig2, axes = Lenses.plot_image_plane(lens, θx, θy, adis;
                                           two_panel = true, source = (0.1, 0.05))
      @test fig2 isa Makie.Figure
      @test length(axes) == 2
      @test all(a -> a isa Makie.Axis, axes)
   end

   @testset "plot_surface_density initializes a figure" begin
      # `:convergence` avoids needing a critical-surface-density calibration.
      fig, ax = Lenses.plot_surface_density(lens, θx, θy, adis; unit = :convergence)
      @test fig isa Makie.Figure
      @test ax  isa Makie.Axis

      # Same call with contours overlaid.
      fig2, ax2 = Lenses.plot_surface_density(lens, θx, θy, adis;
                                              unit = :convergence, plot_contour = true)
      @test fig2 isa Makie.Figure
      @test ax2  isa Makie.Axis
   end

   @testset "plot_magnification_map initializes a figure" begin
      fig, ax = Lenses.plot_magnification_map(lens, θx, θy, adis)
      @test fig isa Makie.Figure
      @test ax  isa Makie.Axis
   end

   @testset "plot_magnification_profile initializes a figure" begin
      fig, ax = Lenses.plot_magnification_profile(lens, θx, θy, adis)
      @test fig isa Makie.Figure
      @test ax  isa Makie.Axis
   end

   # ---- File-saving code path --------------------------------------------------------------------
   # Confirms the `save_plot = true` branch initializes and actually writes a file. The output goes
   # into a temporary directory that is removed automatically on block exit (even if a test fails),
   # so nothing is left behind in the repo.
   @testset "save_plot path writes a file" begin
      mktempdir() do dir
         out = joinpath(dir, "image_plane.png")
         fig, _ = Lenses.plot_image_plane(lens, θx, θy, adis;
                                          plot_critical = false, plot_caustic = false,
                                          save_plot = true, plot_name = out)
         @test fig isa Makie.Figure
         @test isfile(out)
         @test filesize(out) > 0
      end
   end
end