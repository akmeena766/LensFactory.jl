#!!!!!!!!!!!!!! Testing a Multi-component Pixel lens !!!!!!!!!!!!!!
# A MultiPixelLens is a sum of independent square pixels, so it must equal a
# CompositeLens built from the same pixels (primary reference, mirroring
# test_Lenses_MultiPlummerLens.jl), backed by golden numbers from 2-D quadrature.

@testset "Multi-Pixel lens" begin
   xc1 = 0.0; yc1 = 0.0; k1 = 1.0; p1 = 1.0
   xc2 = 0.5; yc2 = 0.5; k2 = 0.7; p2 = 0.6

   # Composite of two single pixels
   lens1 = Lenses.init_CompositeLens([(lens=:PixelLens, x_c=xc1, y_c=yc1, kappa=k1, pixel_size=p1),
                                      (lens=:PixelLens, x_c=xc2, y_c=yc2, kappa=k2, pixel_size=p2)])
   # Equivalent MultiPixelLens
   lens2 = Lenses.init_MultiPixelLens(x_c=[xc1, xc2], y_c=[yc1, yc2], kappa=[k1, k2], pixel_size=[p1, p2])

   # ---- MultiPixelLens == CompositeLens at (1, 1) -------------------------------------
   pot11 = Lenses.get_potential(lens1, xt1, yt1)
   dex11 = Lenses.get_deflection(lens1, xt1, yt1)
   jac11 = Lenses.get_jacobian(lens1, xt1, yt1)

   pot21 = Lenses.get_potential(lens2, xt1, yt1)
   dex21 = Lenses.get_deflection(lens2, xt1, yt1)
   jac21 = Lenses.get_jacobian(lens2, xt1, yt1)

   @test pot11    ≈ pot21    atol=1e-13
   @test dex11[1] ≈ dex21[1] atol=1e-13
   @test dex11[2] ≈ dex21[2] atol=1e-13
   @test jac11[1] ≈ jac21[1] atol=1e-13
   @test jac11[2] ≈ jac21[2] atol=1e-13
   @test jac11[3] ≈ jac21[3] atol=1e-13

   # ---- MultiPixelLens == CompositeLens at (1, 0) -------------------------------------
   pot12 = Lenses.get_potential(lens1, xt2, yt2)
   dex12 = Lenses.get_deflection(lens1, xt2, yt2)
   jac12 = Lenses.get_jacobian(lens1, xt2, yt2)

   pot22 = Lenses.get_potential(lens2, xt2, yt2)
   dex22 = Lenses.get_deflection(lens2, xt2, yt2)
   jac22 = Lenses.get_jacobian(lens2, xt2, yt2)

   @test pot12    ≈ pot22    atol=1e-13
   @test dex12[1] ≈ dex22[1] atol=1e-13
   @test dex12[2] ≈ dex22[2] atol=1e-13
   @test jac12[1] ≈ jac22[1] atol=1e-13
   @test jac12[2] ≈ jac22[2] atol=1e-13
   @test jac12[3] ≈ jac22[3] atol=1e-13

   # ---- Golden values at (1, 1) (independent 2-D quadrature, summed over both pixels) --
   @test pot21    ≈  0.08200543982233413  atol=1e-12
   @test dex21[1] ≈  0.24077150926032356  atol=1e-12
   @test dex21[2] ≈  0.2407715092603236   atol=1e-12
   @test jac21[1] ≈  0.0                  atol=1e-12
   @test jac21[2] ≈  0.0                  atol=1e-12
   @test jac21[3] ≈ -0.33055395779679986  atol=1e-12

   # ---- Array call reproduces the scalar calls ----------------------------------------
   potc = Lenses.get_potential(lens2, [xt1, xt2], [yt1, yt2])
   dexc = Lenses.get_deflection(lens2, [xt1, xt2], [yt1, yt2])
   jacc = Lenses.get_jacobian(lens2, [xt1, xt2], [yt1, yt2])
   @test potc[1] ≈ pot21 atol=1e-13
   @test potc[2] ≈ pot22 atol=1e-13
   @test dexc[1][1] ≈ dex21[1] atol=1e-13
   @test dexc[2][1] ≈ dex21[2] atol=1e-13
   @test dexc[1][2] ≈ dex22[1] atol=1e-13
   @test dexc[2][2] ≈ dex22[2] atol=1e-13
   @test jacc[1][1] ≈ jac21[1] atol=1e-13
   @test jacc[2][1] ≈ jac21[2] atol=1e-13
   @test jacc[3][1] ≈ jac21[3] atol=1e-13
   @test jacc[1][2] ≈ jac22[1] atol=1e-13
   @test jacc[2][2] ≈ jac22[2] atol=1e-13
   @test jacc[3][2] ≈ jac22[3] atol=1e-13

   # ---- A single-component MultiPixelLens reduces to a plain PixelLens -----------------
   single = Lenses.init_PixelLens(x_c=0.3, y_c=-0.2, kappa=0.9, pixel_size=0.8)
   multi1 = Lenses.init_MultiPixelLens(x_c=[0.3], y_c=[-0.2], kappa=[0.9], pixel_size=[0.8])
   @test Lenses.get_potential(multi1, xt1, yt1) ≈ Lenses.get_potential(single, xt1, yt1) atol=1e-13

   # ---- Poisson equation over the overlap of the two pixels ---------------------------
   # (0.25, 0.25) lies inside pixel 1 ([-0.5,0.5]²) and pixel 2 ([0.2,0.8]²)
   j = Lenses.get_jacobian(lens2, 0.25, 0.25)
   @test j[1] + j[2] ≈ 2.0 * (k1 + k2) atol=1e-12

   # ---- get_kappa_gamma over the overlap: κ = Σκ_k, γ1 = 0 by diagonal symmetry --------
   κo, γ1o, γ2o = Lenses.get_kappa_gamma(lens2, 0.25, 0.25, 1.0)
   @test κo  ≈  k1 + k2               atol=1e-12
   @test γ1o ≈  0.0                   atol=1e-12
   @test γ2o ≈ -0.5442814604797956    atol=1e-12
   let jj = Lenses.get_jacobian(lens2, 0.25, 0.25)
      @test κo  ≈ 0.5 * (jj[1] + jj[2]) atol=1e-13
      @test γ2o ≈ jj[3]                 atol=1e-13
   end
end
