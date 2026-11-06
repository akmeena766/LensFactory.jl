#!!!!!!!!!!!!!! Testing a Multi-component Pixel lens !!!!!!!!!!!!!!
# A MultiPixelLens is a complete rectangular grid of equal square pixels, so it must equal a
# CompositeLens built from the same pixels (primary reference, mirroring
# test_Lenses_MultiPlummerLens.jl). One-pixel grids reproduce golden numbers from 2-D quadrature.

@testset "Multi-Pixel lens" begin
   # 3 × 2 grid of 0.5" pixels; edges at x = -0.65, -0.15, 0.35, 0.85 and y = -0.55, -0.05, 0.45
   p   = 0.5
   x_c = [-0.4 -0.4;  0.1  0.1;  0.6  0.6]
   y_c = [-0.3  0.2; -0.3  0.2; -0.3  0.2]
   κ   = [ 1.0  0.7;  0.4  0.9;  0.2  0.5]
   lens = Lenses.init_MultiPixelLens(x_c=x_c, y_c=y_c, kappa=κ, pixel_size=p)
   comp = Lenses.init_CompositeLens([(lens=:PixelLens, x_c=x_c[i, j], y_c=y_c[i, j], kappa=κ[i, j], pixel_size=p)
                                     for j in 1:2 for i in 1:3])

   # ---- MultiPixelLens == CompositeLens at (1, 1), (1, 0) and inside the grid ----------
   for (x, y) in ((xt1, yt1), (xt2, yt2), (0.2, 0.1))
      @test Lenses.get_potential(lens, x, y) ≈ Lenses.get_potential(comp, x, y) atol=1e-13
      dl = Lenses.get_deflection(lens, x, y); dc = Lenses.get_deflection(comp, x, y)
      @test dl[1] ≈ dc[1] atol=1e-13
      @test dl[2] ≈ dc[2] atol=1e-13
      jl = Lenses.get_jacobian(lens, x, y); jc = Lenses.get_jacobian(comp, x, y)
      @test jl[1] ≈ jc[1] atol=1e-13
      @test jl[2] ≈ jc[2] atol=1e-13
      @test jl[3] ≈ jc[3] atol=1e-13
   end

   # ---- Golden values at (1, 1) (independent 2-D quadrature, summed over two pixels) -----
   # Pixel 1: centre (0, 0), κ = 1, size 1; pixel 2: centre (0.5, 0.5), κ = 0.7, size 0.6.
   # They overlap, so each is a one-pixel grid of its own.
   g1 = Lenses.init_MultiPixelLens(x_c=fill(0.0, 1, 1), y_c=fill(0.0, 1, 1), kappa=fill(1.0, 1, 1), pixel_size=1.0)
   g2 = Lenses.init_MultiPixelLens(x_c=fill(0.5, 1, 1), y_c=fill(0.5, 1, 1), kappa=fill(0.7, 1, 1), pixel_size=0.6)
   d1 = Lenses.get_deflection(g1, xt1, yt1); d2 = Lenses.get_deflection(g2, xt1, yt1)
   j1 = Lenses.get_jacobian(g1, xt1, yt1);   j2 = Lenses.get_jacobian(g2, xt1, yt1)
   @test Lenses.get_potential(g1, xt1, yt1) + Lenses.get_potential(g2, xt1, yt1) ≈ 0.08200543982233413 atol=1e-12
   @test d1[1] + d2[1] ≈  0.24077150926032356  atol=1e-12
   @test d1[2] + d2[2] ≈  0.2407715092603236   atol=1e-12
   @test j1[1] + j2[1] ≈  0.0                  atol=1e-12
   @test j1[2] + j2[2] ≈  0.0                  atol=1e-12
   @test j1[3] + j2[3] ≈ -0.33055395779679986  atol=1e-12

   # ---- A one-pixel grid reduces to a plain PixelLens ----------------------------------
   single = Lenses.init_PixelLens(x_c=0.3, y_c=-0.2, kappa=0.9, pixel_size=0.8)
   multi1 = Lenses.init_MultiPixelLens(x_c=fill(0.3, 1, 1), y_c=fill(-0.2, 1, 1), kappa=fill(0.9, 1, 1), pixel_size=0.8)
   @test Lenses.get_potential(multi1, xt1, yt1)     ≈ Lenses.get_potential(single, xt1, yt1)     atol=1e-13
   @test Lenses.get_deflection(multi1, xt1, yt1)[1] ≈ Lenses.get_deflection(single, xt1, yt1)[1] atol=1e-13
   @test Lenses.get_jacobian(multi1, xt1, yt1)[3]   ≈ Lenses.get_jacobian(single, xt1, yt1)[3]   atol=1e-13

   # ---- Array call reproduces the scalar calls -----------------------------------------
   xs = [xt1, xt2, 0.2]; ys = [yt1, yt2, 0.1]
   potc = Lenses.get_potential(lens, xs, ys)
   dexc = Lenses.get_deflection(lens, xs, ys)
   jacc = Lenses.get_jacobian(lens, xs, ys)
   for k in 1:3
      @test potc[k] ≈ Lenses.get_potential(lens, xs[k], ys[k]) atol=1e-13
      ds = Lenses.get_deflection(lens, xs[k], ys[k]); js = Lenses.get_jacobian(lens, xs[k], ys[k])
      @test dexc[1][k] ≈ ds[1] atol=1e-13
      @test dexc[2][k] ≈ ds[2] atol=1e-13
      @test jacc[1][k] ≈ js[1] atol=1e-13
      @test jacc[2][k] ≈ js[2] atol=1e-13
      @test jacc[3][k] ≈ js[3] atol=1e-13
   end

   # ---- Poisson equation and get_kappa_gamma inside pixel (2, 2) -----------------------
   jin = Lenses.get_jacobian(lens, 0.2, 0.1)
   @test jin[1] + jin[2] ≈ 2.0 * κ[2, 2] atol=1e-12
   κo, γ1o, γ2o = Lenses.get_kappa_gamma(lens, 0.2, 0.1, 1.0)
   @test κo  ≈ κ[2, 2]                atol=1e-12
   @test κo  ≈ 0.5 * (jin[1] + jin[2]) atol=1e-13
   @test γ1o ≈ 0.5 * (jin[1] - jin[2]) atol=1e-13
   @test γ2o ≈ jin[3]                 atol=1e-13

   # ---- Points on shared pixel edges and corners are counted in exactly one pixel ------
   uniform = Lenses.init_MultiPixelLens(x_c=x_c, y_c=y_c, kappa=fill(0.8, 3, 2), pixel_size=p)
   for (x, y) in ((-0.15, -0.3), (0.1, -0.05), (-0.15, -0.05), (0.35, -0.05))   # edges and corners
      ju = Lenses.get_jacobian(uniform, x, y)
      @test 0.5 * (ju[1] + ju[2]) ≈ 0.8 atol=1e-12
   end
   je = Lenses.get_jacobian(lens, -0.15, -0.3)       # edge between κ = 1.0 and κ = 0.4
   @test any(abs(0.5 * (je[1] + je[2]) - k) < 1e-12 for k in (κ[1, 1], κ[2, 1]))

   # ---- A grid from get_meshgrid can be used directly -----------------------------------
   gx, gy = Lenses.get_meshgrid(0.5, 0.25, 0.25)
   mg = Lenses.init_MultiPixelLens(x_c=gx, y_c=gy, kappa=fill(0.5, size(gx)), pixel_size=0.25)
   cg = Lenses.init_CompositeLens([(lens=:PixelLens, x_c=gx[k], y_c=gy[k], kappa=0.5, pixel_size=0.25) for k in eachindex(gx)])
   @test Lenses.get_deflection(mg, xt1, yt1)[1] ≈ Lenses.get_deflection(cg, xt1, yt1)[1] atol=1e-13

   # ---- Constructor checks -------------------------------------------------------------
   @test_throws ArgumentError Lenses.init_MultiPixelLens(x_c=x_c, y_c=y_c, kappa=κ[1:2, :], pixel_size=p)
   @test_throws ArgumentError Lenses.init_MultiPixelLens(x_c=permutedims(x_c), y_c=permutedims(y_c), kappa=permutedims(κ), pixel_size=p)   # x along 2nd dim
   @test_throws ArgumentError Lenses.init_MultiPixelLens(x_c=x_c, y_c=y_c, kappa=κ, pixel_size=0.4)
   @test_throws ArgumentError Lenses.init_MultiPixelLens(x_c=x_c, y_c=y_c, kappa=κ, pixel_size=0.0)
   @test_throws UndefKeywordError Lenses.init_MultiPixelLens(x_c=x_c, y_c=y_c, kappa=κ)
end
