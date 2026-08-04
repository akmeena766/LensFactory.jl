#!!!!!!!!!!!!!! Testing a single square Pixel lens !!!!!!!!!!!!!!
# A PixelLens is a square patch of uniform convergence κ. Its potential is
#     ψ(θ) = (κ/π) ∫∫_pixel ln|θ - θ'| d²θ',
# evaluated in closed form. The golden numbers below come from an independent direct
# 2-D quadrature of that integral (agreement to machine precision); deflection and
# jacobian are the gradient and Hessian of ψ.
#
# Pixel lens quantities do not depend on the source distance, so (unlike the mass
# lenses) they are compared to get_* directly, without the adis factor.

@testset "Pixel lens" begin
   # Square patch centred at the origin, κ = 1, side = 1  ⇒  patch spans [-0.5, 0.5]²
   lens = Lenses.init_PixelLens(x_c=0.0, y_c=0.0, kappa=1.0, pixel_size=1.0)

   # ---- Point (xt1, yt1) = (1, 1) -----------------------------------------------------
   pot1 = Lenses.get_potential(lens, xt1, yt1)
   dex1 = Lenses.get_deflection(lens, xt1, yt1)
   jac1 = Lenses.get_jacobian(lens, xt1, yt1)
   @test pot1    ≈  0.10998270017740712  atol=1e-12
   @test dex1[1] ≈  0.15983234777274674  atol=1e-12
   @test dex1[2] ≈  0.15983234777274677  atol=1e-12
   @test jac1[1] ≈  0.0                  atol=1e-12
   @test jac1[2] ≈  0.0                  atol=1e-12
   @test jac1[3] ≈ -0.16260084616071638  atol=1e-12

   # ---- Point (xt2, yt2) = (1, 0) -----------------------------------------------------
   pot2 = Lenses.get_potential(lens, xt2, yt2)
   dex2 = Lenses.get_deflection(lens, xt2, yt2)
   jac2 = Lenses.get_jacobian(lens, xt2, yt2)
   @test pot2    ≈  0.0012751343922037528 atol=1e-12
   @test dex2[1] ≈  0.3133991464120883    atol=1e-12
   @test dex2[2] ≈  0.0                   atol=1e-12
   @test jac2[1] ≈ -0.2951672353008665    atol=1e-12
   @test jac2[2] ≈  0.2951672353008665    atol=1e-12
   @test jac2[3] ≈  0.0                   atol=1e-12

   # ---- Array call must reproduce the scalar calls, element by element -----------------
   potc = Lenses.get_potential(lens, [xt1, xt2], [yt1, yt2])
   dexc = Lenses.get_deflection(lens, [xt1, xt2], [yt1, yt2])
   jacc = Lenses.get_jacobian(lens, [xt1, xt2], [yt1, yt2])
   @test potc[1] ≈ pot1       atol=1e-13
   @test potc[2] ≈ pot2       atol=1e-13
   @test dexc[1][1] ≈ dex1[1] atol=1e-13
   @test dexc[2][1] ≈ dex1[2] atol=1e-13
   @test dexc[1][2] ≈ dex2[1] atol=1e-13
   @test dexc[2][2] ≈ dex2[2] atol=1e-13
   @test jacc[1][1] ≈ jac1[1] atol=1e-13
   @test jacc[2][1] ≈ jac1[2] atol=1e-13
   @test jacc[3][1] ≈ jac1[3] atol=1e-13
   @test jacc[1][2] ≈ jac2[1] atol=1e-13
   @test jacc[2][2] ≈ jac2[2] atol=1e-13
   @test jacc[3][2] ≈ jac2[3] atol=1e-13

   # ---- Int coordinates are promoted and give identical results -----------------------
   @test Lenses.get_potential(lens, 1, 1) ≈ pot1 atol=1e-13

   # ---- Deflection = ∇ψ  (central finite differences) ---------------------------------
   h = 1e-6
   for (x, y) in ((xt1, yt1), (xt2, yt2), (-0.9, 0.4))
      dx, dy = Lenses.get_deflection(lens, x, y)
      gx = (Lenses.get_potential(lens, x + h, y) - Lenses.get_potential(lens, x - h, y)) / (2h)
      gy = (Lenses.get_potential(lens, x, y + h) - Lenses.get_potential(lens, x, y - h)) / (2h)
      @test dx ≈ gx atol=1e-6
      @test dy ≈ gy atol=1e-6
   end

   # ---- Poisson equation: ψxx + ψyy = 2κ inside the patch, ≈ 0 far outside -------------
   for (x, y) in ((0.0, 0.0), (0.2, -0.3), (-0.45, 0.1))       # strictly inside
      j = Lenses.get_jacobian(lens, x, y)
      @test j[1] + j[2] ≈ 2.0 atol=1e-12
   end
   for (x, y) in ((12.0, 0.0), (0.0, -9.0), (7.0, 7.0))        # far outside
      j = Lenses.get_jacobian(lens, x, y)
      @test j[1] + j[2] ≈ 0.0 atol=1e-12
   end

   # ---- Convergence scales linearly with κ --------------------------------------------
   lens2 = Lenses.init_PixelLens(x_c=0.0, y_c=0.0, kappa=2.5, pixel_size=1.0)
   @test Lenses.get_potential(lens2, xt1, yt1) ≈ 2.5 * pot1 atol=1e-12
end
