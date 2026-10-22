@testset "overlapping layers: profile realism" begin
    w, ρ = 2.0, 0.1
    xs = range(1e-9, w - 1e-9, length=201)          # a spin exactly on a wall counts both sides' layers (measure zero): sample inside the gap
    @testset "package agrees with closed form" begin
        for (shape, hs) in (("linear", (0.5, 0.6, 0.75, 0.9, 1.0) .* w), ("exponential", (0.05, 0.25, 0.5, 1.0) .* w))
            for h in hs, x in xs
                @test summed_rate_package(shape, ρ, h, w, x) ≈ summed_rate_closed(shape, ρ, h, w, x) rtol=1e-10 atol=1e-14
            end
        end
    end
    @testset "linear and exponential: never higher mid-gap than at the wall" begin
        for h in (0.25, 0.5, 0.6, 0.75, 0.9, 1.0) .* w
            @test is_realistic("linear", ρ, h, w)
        end
        for h in (0.05, 0.1, 0.25, 0.5, 0.75, 1.0) .* w
            @test is_realistic("exponential", ρ, h, w)
        end
    end
    @testset "closed forms of the overlap plateau and the uniform limit" begin
        # linear, w/2 < h ≤ w: flat plateau (ρ/h)(4 − 2w/h) on the overlap region [w − h, h]; at h = w uniform 2ρ/w
        h = 0.75w
        @test summed_rate_closed("linear", ρ, h, w, w / 2) ≈ ρ / h * (4 - 2w / h)
        @test summed_rate_closed("linear", ρ, h, w, w - h + 0.01) ≈ ρ / h * (4 - 2w / h)
        @test all(summed_rate_closed("linear", ρ, w, w, x) ≈ 2ρ / w for x in xs)
        # exponential: mid/wall = 2e^(−w/2h)/(1 + e^(−w/h))
        h = 0.5w
        m = overlap_metrics("exponential", h, w)
        @test m.mid_over_wall ≈ 2exp(-w / (2h)) / (1 + exp(-w / h)) rtol=1e-10
        @test overlap_metrics("linear", 0.75w, w).overlap_width ≈ 2 * 0.75w - w
        @test overlap_metrics("linear", 0.25w, w).overlap_width == 0.0
        @test overlap_metrics("linear", 0.25w, w).mid_over_wall == 0.0
    end
    @testset "counter-example: step with h > w/2 is unrealistic" begin
        @test !is_realistic("step", ρ, 0.75w, w)
        @test is_realistic("step", ρ, 0.5w, w)
    end
end
