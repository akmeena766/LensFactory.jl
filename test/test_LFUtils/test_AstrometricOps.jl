@testset "AstrometricOps" begin
   # Short alias for the (nested) module under test
   AO = LFUtils.AstrometricOps

   # A reference position well away from the poles and RA = 0 wrap-around
   ra_ref, dec_ref = 150.0, 30.0


   @testset "gnomonic_offsets_arcsec" begin

      @testset "reference position maps to the origin" begin
         x, y = AO.gnomonic_offsets_arcsec(ra_ref, dec_ref, ra_ref, dec_ref)
         @test x ≈ 0.0 atol=1e-12
         @test y ≈ 0.0 atol=1e-12
      end


      @testset "scalar method returns scalars" begin
         x, y = AO.gnomonic_offsets_arcsec(ra_ref, dec_ref, 150.2, 30.15)
         @test x isa Float64
         @test y isa Float64
         # Regression values (computed independently from the projection formula)
         @test x ≈ -622.5973869267 rtol=1e-10
         @test y ≈  540.5470166574 rtol=1e-10
      end


      @testset "scalar and vector methods agree" begin
         ra_cat  = [150.05, 149.98, 150.2, 151.3, 148.7]
         dec_cat = [ 30.03,  29.90,  30.15, 31.20,  28.50]
         xv, yv = AO.gnomonic_offsets_arcsec(ra_ref, dec_ref, ra_cat, dec_cat)

         @test xv isa Vector{Float64}
         @test yv isa Vector{Float64}
         @test length(xv) == length(ra_cat)
         @test length(yv) == length(ra_cat)

         for i in eachindex(ra_cat)
            xs, ys = AO.gnomonic_offsets_arcsec(ra_ref, dec_ref, ra_cat[i], dec_cat[i])
            @test xs ≈ xv[i] rtol=1e-12
            @test ys ≈ yv[i] rtol=1e-12
         end
      end


      @testset "north-up, east-left convention" begin
         # A point to the east (larger RA) has NEGATIVE x (x increases towards west)
         xe, ye = AO.gnomonic_offsets_arcsec(ra_ref, dec_ref, ra_ref + 0.001, dec_ref)
         @test xe < 0.0
         @test abs(ye) < 1e-4                     # essentially no north/south motion

         # A point to the west (smaller RA) has POSITIVE x
         xw, yw = AO.gnomonic_offsets_arcsec(ra_ref, dec_ref, ra_ref - 0.001, dec_ref)
         @test xw > 0.0
         @test xw ≈ -xe rtol=1e-6                 # symmetric about the reference

         # A point to the north (larger Dec) has POSITIVE y
         xn, yn = AO.gnomonic_offsets_arcsec(ra_ref, dec_ref, ra_ref, dec_ref + 0.001)
         @test yn > 0.0
         @test abs(xn) < 1e-9                     # no east/west motion at ΔRA = 0
      end


      @testset "small-angle limits" begin
         # Near the reference the tangent plane is locally flat:
         #   x ≈ -ΔRA·cos(dec_ref),  y ≈ ΔDec   (all in arcsec)
         Δdeg = 1e-3
         as_per_deg = 3600.0

         xe, _  = AO.gnomonic_offsets_arcsec(ra_ref, dec_ref, ra_ref + Δdeg, dec_ref)
         @test xe ≈ -Δdeg * as_per_deg * cosd(dec_ref) rtol=1e-6

         _, yn  = AO.gnomonic_offsets_arcsec(ra_ref, dec_ref, ra_ref, dec_ref + Δdeg)
         @test yn ≈ Δdeg * as_per_deg rtol=1e-6
      end
   end


   @testset "gnomonic_offsets_radec" begin

      @testset "origin maps back to the reference position" begin
         ra, dec = AO.gnomonic_offsets_radec(ra_ref, dec_ref, 0.0, 0.0)
         @test ra  ≈ ra_ref  atol=1e-12
         @test dec ≈ dec_ref atol=1e-12
         @test !isnan(ra)                         # rho = 0 must not produce NaN
         @test !isnan(dec)
      end


      @testset "scalar method returns scalars" begin
         ra, dec = AO.gnomonic_offsets_radec(ra_ref, dec_ref, 10.0, -8.0)
         @test ra  isa Float64
         @test dec isa Float64
      end


      @testset "scalar and vector methods agree" begin
         x_as = [10.0, -5.0, 120.0, -300.0, 500.0]
         y_as = [-8.0, 15.0, -30.0,  200.0, -100.0]
         rav, decv = AO.gnomonic_offsets_radec(ra_ref, dec_ref, x_as, y_as)

         @test rav  isa Vector{Float64}
         @test decv isa Vector{Float64}
         @test length(rav)  == length(x_as)
         @test length(decv) == length(x_as)

         for i in eachindex(x_as)
            ras, decs = AO.gnomonic_offsets_radec(ra_ref, dec_ref, x_as[i], y_as[i])
            @test ras  ≈ rav[i]  rtol=1e-12
            @test decs ≈ decv[i] rtol=1e-12
         end
      end


      @testset "convention: west-positive x, north-positive y" begin
         # Positive x is west -> smaller RA; positive y is north -> larger Dec
         ra_w, _   = AO.gnomonic_offsets_radec(ra_ref, dec_ref, 10.0, 0.0)
         @test ra_w < ra_ref

         ra_e, _   = AO.gnomonic_offsets_radec(ra_ref, dec_ref, -10.0, 0.0)
         @test ra_e > ra_ref

         _, dec_n  = AO.gnomonic_offsets_radec(ra_ref, dec_ref, 0.0, 10.0)
         @test dec_n > dec_ref
      end
   end


   @testset "forward / inverse round-trips" begin
      # The two functions must be exact mutual inverses (to floating-point precision).

      @testset "RA/Dec -> offsets -> RA/Dec" begin
         # Spread of catalog points, including a few degrees from the reference
         ra_cat  = [150.05, 151.3, 148.7, 150.0, 153.0, 155.0, 145.0]
         dec_cat = [ 30.03, 31.2,  28.5,  32.0,  30.0,  35.0,  25.0]

         x, y   = AO.gnomonic_offsets_arcsec(ra_ref, dec_ref, ra_cat, dec_cat)
         ra2, dec2 = AO.gnomonic_offsets_radec(ra_ref, dec_ref, x, y)

         @test ra2  ≈ ra_cat  atol=1e-9
         @test dec2 ≈ dec_cat atol=1e-9
      end

      @testset "offsets -> RA/Dec -> offsets" begin
         x_as = [10.0, -5.0, 120.0, -300.0, 500.0]
         y_as = [-8.0, 15.0, -30.0,  200.0, -100.0]

         ra, dec = AO.gnomonic_offsets_radec(ra_ref, dec_ref, x_as, y_as)
         x2, y2  = AO.gnomonic_offsets_arcsec(ra_ref, dec_ref, ra, dec)

         @test x2 ≈ x_as atol=1e-6
         @test y2 ≈ y_as atol=1e-6
      end

      @testset "high-declination reference" begin
         # Convention and round-trip should hold well away from the celestial equator
         ra0, dec0 = 80.0, 85.0
         ra_cat  = [85.0, 75.0, 82.0]
         dec_cat = [85.5, 84.5, 86.0]

         x, y   = AO.gnomonic_offsets_arcsec(ra0, dec0, ra_cat, dec_cat)
         ra2, dec2 = AO.gnomonic_offsets_radec(ra0, dec0, x, y)

         @test ra2  ≈ ra_cat  atol=1e-9
         @test dec2 ≈ dec_cat atol=1e-9
      end
   end
end
