@testset "Sources" begin
   # Common grid: symmetric about the origin, includes (0, 0) exactly
   grid_x, grid_y = Lenses.get_meshgrid(3.0, 3.0, 0.02)
   h = 0.02                                       # pixel size
   r_grid = @. sqrt(grid_x^2 + grid_y^2)          # radial distance of every pixel


   @testset "disk" begin
      # NOTE: θ_r is chosen such that no grid pixel lies exactly on the disk boundary — with,
      # e.g., θ_r = 0.8 many pixels sit at exactly r = θ_r (scaled 3-4-5 triples like
      # 0.48² + 0.64² = 0.8²), and their membership would be decided by the last floating-point
      # bit, which is not translation-invariant.
      θ_r = 0.79
      src = Sources.disk(grid_x, grid_y, θ_r, (0.0, 0.0))
   
      # Binary profile: A inside, 0 outside
      @test all(v -> v == 0.0 || v == 1.0, src)
      @test all(src[r_grid .<= θ_r] .== 1.0)
      @test all(src[r_grid .>  θ_r] .== 0.0)
   
      # Pixel-counted area matches π θ_r²
      @test sum(src) * h^2 ≈ π * θ_r^2 rtol=0.03
   
      # Amplitude scaling
      src_A = Sources.disk(grid_x, grid_y, θ_r, (0.0, 0.0); A=2.5)
      @test src_A ≈ 2.5 .* src
   
      # Translation by an integer number of pixels (0.5" = 25 pixels) shifts the array
      src_shift = Sources.disk(grid_x, grid_y, θ_r, (0.5, 0.0))
      @test sum(src_shift) == sum(src)                       # fully inside the grid
      @test src_shift[26:end, :] == src[1:end-25, :]
    end


   @testset "gaussian" begin
      σx, σy = 0.3, 0.2
      src = Sources.gaussian(grid_x, grid_y, σx, σy, (0.0, 0.0))

      # Total flux integrates to A
      @test sum(src) * h^2 ≈ 1.0 rtol=1e-6

      # Peak value and position
      @test maximum(src) ≈ 1.0 / (2.0 * π * σx * σy) rtol=1e-12
      pk = argmax(src)
      @test abs(grid_x[pk]) < 1e-12
      @test abs(grid_y[pk]) < 1e-12

      # 180° rotation symmetry for a centered profile
      @test src ≈ reverse(reverse(src, dims=1), dims=2) rtol=1e-10

      # Anisotropy: exp(-1/2) drop at 1σ along each axis (nearest grid point to ±1σ)
      ix = argmin(abs.(grid_x[:, 1] .- σx))
      iy = argmin(abs.(grid_y[1, :] .- σy))
      i0 = argmin(abs.(grid_x[:, 1]))
      j0 = argmin(abs.(grid_y[1, :]))
      @test src[ix, j0] / src[i0, j0] ≈ exp(-0.5) rtol=1e-10
      @test src[i0, iy] / src[i0, j0] ≈ exp(-0.5) rtol=1e-10

      # Amplitude scaling
      @test Sources.gaussian(grid_x, grid_y, σx, σy, (0.0, 0.0); A=2.5) ≈ 2.5 .* src

      # Off-center profile peaks at the requested position
      src_off = Sources.gaussian(grid_x, grid_y, σx, σy, (0.5, -1.0))
      pk_off = argmax(src_off)
      @test grid_x[pk_off] ≈ 0.5  atol=h
      @test grid_y[pk_off] ≈ -1.0 atol=h
   end


   @testset "sersic" begin
      # n = 1/2 Sersic is exactly a Gaussian with σ = θ_e / √(2 b_n)
      θ_e = 0.4
      bn = Sources.b_n(0.5)
      σ_eq = θ_e / sqrt(2.0 * bn)
      @test Sources.sersic(grid_x, grid_y, 0.5, θ_e, (0.0, 0.0)) ≈
            Sources.gaussian(grid_x, grid_y, σ_eq, σ_eq, (0.0, 0.0)) rtol=1e-12

      # θ_e is the half-light radius and the total flux is A — this validates b_n(n).
      # (Both b_n branches: n = 0.2 uses the polynomial fit, n = 1 the asymptotic series.)
      for n in (0.2, 1.0)
         src = Sources.sersic(grid_x, grid_y, n, θ_e, (0.0, 0.0))
         flux_total  = sum(src) * h^2
         flux_within = sum(src[r_grid .<= θ_e]) * h^2
         @test flux_total ≈ 1.0 rtol=0.01
         @test flux_within / flux_total ≈ 0.5 atol=0.02
      end

      # 180° rotation symmetry for a centered profile
      src1 = Sources.sersic(grid_x, grid_y, 1.0, θ_e, (0.0, 0.0))
      @test src1 ≈ reverse(reverse(src1, dims=1), dims=2) rtol=1e-10

      # Amplitude scaling
      @test Sources.sersic(grid_x, grid_y, 1.0, θ_e, (0.0, 0.0); A=2.5) ≈ 2.5 .* src1

      # Steeper index ⇒ more concentrated: higher peak, same total flux
      src4 = Sources.sersic(grid_x, grid_y, 4.0, θ_e, (0.0, 0.0))
      @test maximum(src4) > maximum(src1)
   end
end
