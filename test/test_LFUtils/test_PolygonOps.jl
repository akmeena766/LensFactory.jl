@testset "PolygonOps" begin
   # Short alias for the (nested) module under test
   PO = LFUtils.PolygonOps


   # ----------------------------------------------------------------------------------------------
   # shoelace: polygon area
   # ----------------------------------------------------------------------------------------------
   @testset "shoelace" begin

      @testset "unit square" begin
         # Explicitly closed and implicitly open inputs give the same area
         closed = [[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0],[0.0,0.0]]
         open   = [[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0]]
         @test PO.shoelace(closed) ≈ 1.0
         @test PO.shoelace(open)   ≈ 1.0
      end


      @testset "right triangle" begin
         @test PO.shoelace([[0.0,0.0],[4.0,0.0],[0.0,3.0]]) ≈ 6.0
      end


      @testset "orientation independence" begin
         # abs() in the formula makes the result independent of winding direction
         ccw = [[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0]]
         cw  = [[0.0,0.0],[0.0,1.0],[1.0,1.0],[1.0,0.0]]
         @test PO.shoelace(ccw) ≈ PO.shoelace(cw) ≈ 1.0
      end


      @testset "translation invariance" begin
         base    = [[0.0,0.0],[2.0,0.0],[2.0,3.0],[0.0,3.0]]
         shifted = [[p[1]+100.0, p[2]-50.0] for p in base]
         @test PO.shoelace(base) ≈ 6.0
         @test PO.shoelace(shifted) ≈ PO.shoelace(base)
      end


      @testset "degenerate (collinear) polygon has zero area" begin
         @test PO.shoelace([[0.0,0.0],[1.0,0.0],[2.0,0.0]]) ≈ 0.0 atol=1e-15
      end


      @testset "integer coordinates are accepted" begin
         @test PO.shoelace([[0,0],[2,0],[2,2],[0,2]]) ≈ 4.0
      end


      @testset "input is not mutated" begin
         # An open polygon is closed on a copy, so the caller's vector is untouched
         open = [[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0]]
         n0   = length(open)
         PO.shoelace(open)
         @test length(open) == n0
         @test open[end] == [0.0,1.0]
      end
   end


   # ----------------------------------------------------------------------------------------------
   # hao_sun: point-in-polygon (+1 inside, 0 outside, -1 on the boundary)
   # ----------------------------------------------------------------------------------------------
   @testset "hao_sun" begin

      # Closed unit square (first vertex repeated at the end, as the algorithm expects)
      square = [[0.0,0.0],[1.0,0.0],[1.0,1.0],[0.0,1.0],[0.0,0.0]]

      @testset "interior points" begin
         @test PO.hao_sun([0.5,0.5], square) == 1
         @test PO.hao_sun([0.1,0.9], square) == 1
      end


      @testset "exterior points" begin
         @test PO.hao_sun([1.5,0.5], square) == 0
         @test PO.hao_sun([-1.0,-1.0], square) == 0
         @test PO.hao_sun([0.5,2.0], square) == 0
         @test PO.hao_sun([2.0,2.0], square) == 0
      end


      @testset "boundary points (edges)" begin
         @test PO.hao_sun([0.5,0.0], square) == -1   # bottom
         @test PO.hao_sun([1.0,0.5], square) == -1   # right
         @test PO.hao_sun([0.5,1.0], square) == -1   # top
         @test PO.hao_sun([0.0,0.5], square) == -1   # left
      end


      @testset "boundary points (vertices)" begin
         @test PO.hao_sun([0.0,0.0], square) == -1
         @test PO.hao_sun([1.0,1.0], square) == -1
      end


      @testset "concave (L-shaped) polygon" begin
         # An L / notch: the concave region must be classified as outside
         Lshape = [[0.0,0.0],[2.0,0.0],[2.0,2.0],[1.0,2.0],[1.0,1.0],[0.0,1.0],[0.0,0.0]]
         @test PO.hao_sun([0.5,0.5], Lshape) == 1   # in the lower arm
         @test PO.hao_sun([1.5,1.5], Lshape) == 1   # in the upper arm
         @test PO.hao_sun([1.5,0.5], Lshape) == 1   # in the right arm
         @test PO.hao_sun([0.5,1.5], Lshape) == 0   # in the notch -> outside
      end


      @testset "consistency with shoelace membership" begin
         # Every sampled grid point is classified inside/outside consistently
         # with lying within the axis-aligned square (0,1) x (0,1).
         for xi in 0.05:0.1:0.95, yi in 0.05:0.1:0.95
            @test PO.hao_sun([xi,yi], square) == 1
         end
      end
   end


   # ----------------------------------------------------------------------------------------------
   # bilinear_interpolation
   # ----------------------------------------------------------------------------------------------
   @testset "bilinear_interpolation" begin

      @testset "reproduces values at grid nodes" begin
         df = Float64[1 2 3; 4 5 6; 7 8 9]
         for i in 1:2, j in 1:2               # px+1, py+1 must stay in bounds
            @test PO.bilinear_interpolation(float(i), float(j), df) ≈ df[i,j]
         end
      end


      @testset "midpoint is the average of the four corners" begin
         df = Float64[0 0; 0 4]               # df[1,1]=0, df[1,2]=0, df[2,1]=0, df[2,2]=4
         # Centre of the single cell -> mean of the four corner values
         @test PO.bilinear_interpolation(1.5, 1.5, df) ≈ 1.0
      end


      @testset "reproduces a linear field exactly" begin
         # f(i,j) = 2i + 3j - 1 is bilinear, so interpolation is exact everywhere
         nx, ny = 5, 4
         df = Float64[2i + 3j - 1 for i in 1:nx, j in 1:ny]
         for (x,y) in [(1.0,1.0),(2.3,1.7),(4.9,3.5),(3.0,2.0),(1.25,3.9)]
            @test PO.bilinear_interpolation(x, y, df) ≈ 2x + 3y - 1
         end
      end


      @testset "edge along a single axis is linear" begin
         df = Float64[0 0; 10 0]              # varies only along the first index
         @test PO.bilinear_interpolation(1.25, 1.0, df) ≈ 2.5
         @test PO.bilinear_interpolation(1.75, 1.0, df) ≈ 7.5
      end


      @testset "out-of-range queries throw" begin
         df = Float64[1 2 3; 4 5 6; 7 8 9]    # size (3,3): valid x,y in [1,3)
         @test_throws ArgumentError PO.bilinear_interpolation(0.5, 2.0, df)  # x < 1
         @test_throws ArgumentError PO.bilinear_interpolation(2.0, 0.5, df)  # y < 1
         @test_throws ArgumentError PO.bilinear_interpolation(3.0, 2.0, df)  # x == nx
         @test_throws ArgumentError PO.bilinear_interpolation(2.0, 3.0, df)  # y == ny
      end
   end


   # ----------------------------------------------------------------------------------------------
   # fit_ellipse
   # ----------------------------------------------------------------------------------------------
   @testset "fit_ellipse" begin

      # Sample n points on an ellipse of centre (cx,cy), semi-axes (a,b), rotation φ (radians)
      function sample_ellipse(cx, cy, a, b, φ; n=60)
         t  = [2π*(k-1)/n for k in 1:n]
         x  = [cx + a*cos(ti)*cos(φ) - b*sin(ti)*sin(φ) for ti in t]
         y  = [cy + a*cos(ti)*sin(φ) + b*sin(ti)*cos(φ) for ti in t]
         return x, y
      end

      # Ellipse orientation is only defined modulo 180 degrees
      θmod(a) = mod(a, 180.0)


      @testset "axis-aligned ellipse" begin
         x, y = sample_ellipse(2.0, 3.0, 4.0, 2.0, 0.0)
         cx, cy, rx, ry, θ = PO.fit_ellipse(x, y)
         @test cx ≈ 2.0 atol=1e-8
         @test cy ≈ 3.0 atol=1e-8
         @test rx ≈ 4.0 atol=1e-8
         @test ry ≈ 2.0 atol=1e-8
         @test rx ≥ ry                        # semi-major first, by construction
         @test θmod(θ) ≈ 0.0 atol=1e-6
      end


      @testset "rotated ellipse" begin
         x, y = sample_ellipse(2.0, 3.0, 4.0, 2.0, deg2rad(30.0))
         cx, cy, rx, ry, θ = PO.fit_ellipse(x, y)
         @test cx ≈ 2.0 atol=1e-8
         @test cy ≈ 3.0 atol=1e-8
         @test rx ≈ 4.0 atol=1e-8
         @test ry ≈ 2.0 atol=1e-8
         @test θmod(θ) ≈ 30.0 atol=1e-6
      end


      @testset "circle recovers equal axes" begin
         x, y = sample_ellipse(0.0, 0.0, 5.0, 5.0, 0.0)
         cx, cy, rx, ry, _ = PO.fit_ellipse(x, y)
         @test cx ≈ 0.0 atol=1e-8
         @test cy ≈ 0.0 atol=1e-8
         @test rx ≈ 5.0 atol=1e-8
         @test ry ≈ 5.0 atol=1e-8
      end


      @testset "all fitted points satisfy the conic" begin
         # Points on the fitted ellipse must lie (near) on it: the normalized
         # radial coordinate in the ellipse frame should be ~1 for every sample.
         cx0, cy0, a0, b0, φ0 = 1.5, -2.0, 6.0, 3.0, deg2rad(-20.0)
         x, y = sample_ellipse(cx0, cy0, a0, b0, φ0)
         cx, cy, rx, ry, θ = PO.fit_ellipse(x, y)
         φ = deg2rad(θ)
         for i in eachindex(x)
            dx = x[i] - cx
            dy = y[i] - cy
            u  =  dx*cos(φ) + dy*sin(φ)      # rotate into the ellipse frame
            v  = -dx*sin(φ) + dy*cos(φ)
            @test (u/rx)^2 + (v/ry)^2 ≈ 1.0 atol=1e-6
         end
      end
   end
end
