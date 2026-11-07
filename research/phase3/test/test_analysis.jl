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

@testset "analysis helpers: censoring, few points, reference se (Group A review)" begin
    x = [0.1, 0.2, 0.4, 0.8]
    y = 3.0 .* x .^ 1.5
    # threshold_at returns (value, status)
    @test threshold_status(x, y, 3.0 * 0.3^1.5) == (0.3, :crossed) || threshold_status(x, y, 3.0 * 0.3^1.5)[2] == :crossed
    @test threshold_status(x, y, 3.0 * 0.3^1.5)[1] ≈ 0.3 rtol=1e-10
    @test threshold_status(x, y, 100.0) == (0.8, :below_range)        # never reached: a lower bound, not a crossing
    @test threshold_status(x, y, 1e-6)[2] == :above_range              # first point already above
    # non-monotone: the threshold is where y first exceeds the level from small x (every smaller x within the level)
    @test threshold_status([0.1, 0.2, 0.4], [0.01, 0.05, 0.03], 0.02)[2] == :crossed
    @test threshold_status([0.1, 0.2, 0.4], [0.01, 0.05, 0.03], 0.02)[1] < 0.2
    # duplicate x: the larger y counts (conservative), independent of input order
    a = threshold_status([0.2, 0.2, 0.4], [0.01, 0.03, 0.05], 0.02)
    b = threshold_status([0.2, 0.2, 0.4], [0.03, 0.01, 0.05], 0.02)
    @test isequal(a, b) && a[2] == :above_range               # isequal: NaN == NaN is false
    # fit_power: few points
    @test isnan(fit_power([0.1, 0.2], [1.0, 2.0])[3])                  # n = 2: slope exact, no error estimate
    @test_throws ArgumentError fit_power([0.1], [1.0])                 # n = 1
    @test_throws ArgumentError fit_power([0.1, 0.1], [1.0, 2.0])       # no spread in x
    # reference se from an independent matrix least squares
    noisy = y .* [1.01, 0.99, 1.02, 0.98]
    X = hcat(ones(4), log.(x)); β = X \ log.(noisy); r = log.(noisy) .- X * β
    se_ref = sqrt(sum(r .^ 2) / 2 * inv(X' * X)[2, 2])
    (_, p, se) = fit_power(x, noisy)
    @test p ≈ β[2] rtol=1e-12
    @test se ≈ se_ref rtol=1e-10
end
