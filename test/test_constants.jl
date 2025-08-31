@testset "Cosmology" begin
    @testset "Distance units" begin
        @test DIST_AU  < DIST_PC
        @test DIST_PC  * 1e3 ≈ DIST_KPC atol=1e-15 rtol=1e-15
        @test DIST_KPC * 1e3 ≈ DIST_MPC atol=1e-15 rtol=1e-15
        @test DIST_MPC * 1e3 ≈ DIST_GPC atol=1e-15 rtol=1e-15
    end

    @testset "Angular units" begin
        @test ANGLE_MICROARCSEC * 1e3 ≈ ANGLE_MILLIARCSEC atol=1e-15 rtol=1e-15
        @test ANGLE_MILLIARCSEC * 1e3 ≈ ANGLE_ARCSEC atol=1e-15 rtol=1e-15
        @test ANGLE_ARCSEC * 60 ≈ ANGLE_ARCMIN atol=1e-15 rtol=1e-15
        @test ANGLE_ARCMIN * 60 ≈ ANGLE_DEGREE atol=1e-15 rtol=1e-15
    end

    @testset "Temporal units" begin
        @test YEAR2SECOND ≈ YEAR2DAY * DAY2SECOND atol=1e-15 rtol=1e-15
        @test YEAR2HOUR   ≈ YEAR2DAY * DAY2HOUR atol=1e-15 rtol=1e-15
        @test YEAR2DAY * DAY2HOUR ≈ YEAR2HOUR atol=1e-15 rtol=1e-15
    end
end