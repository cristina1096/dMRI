@testset "analysis helpers" begin
    x = [0.1, 0.2, 0.4, 0.8]
    y = 3.0 .* x .^ 1.5
    (a, p, se) = fit_power(x, y)
    @test a ≈ 3.0 rtol=1e-10
    @test p ≈ 1.5 rtol=1e-10
    @test se < 1e-8
    noisy = y .* [1.01, 0.99, 1.02, 0.98]
    (_, p2, se2) = fit_power(x, noisy)
    @test abs(p2 - 1.5) < 0.05 && se2 > 0
    @test threshold_at(x, y, 3.0 * 0.4^1.5) ≈ 0.4 rtol=1e-10
    @test threshold_at(x, y, 3.0 * 0.3^1.5) ≈ 0.3 rtol=1e-10          # log-log interpolation is exact for a power law
    @test isnan(threshold_at(x, y, 1e-6))
    @test threshold_at(x, y, 100.0) == 0.8                              # every point below: the largest scanned x
    @test rel_error(1.1, 1.0) ≈ 0.1
    @test_throws ArgumentError fit_power([0.1, 0.2], [1.0, -1.0])
end
