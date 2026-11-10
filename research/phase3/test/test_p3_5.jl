@testset "P3.5 Task 1: rectangular brute force and overlap classification" begin
    outer = (shape="step", h=0.3, L=0.3, pre=1.0)
    # square: (sx, sy) form equals the single-s form
    let s = 2.5
        pts = [(s * (rand() - 0.5), s * (rand() - 0.5)) for _ in 1:100]
        @test all(brute_rate_profile(x, y; R=1.0, sx=s, sy=s, inner=nothing, outer=outer) ==
                  brute_rate_profile(x, y; R=1.0, s=s, inner=nothing, outer=outer) for (x, y) in pts)
    end
    # rectangular (2.5, 4.0): gap 0.5 in x (step L = 0.3 overlaps), gap 2.0 in y (no overlap at L = 0.5)
    @test brute_rate_profile(1.25, 0.0; R=1.0, sx=2.5, sy=4.0, inner=nothing, outer=outer) == 2.0
    @test brute_rate_profile(0.0, 2.0; R=1.0, sx=2.5, sy=4.0, inner=nothing, outer=(shape="step", h=0.5, L=0.5, pre=1.0)) == 0.0
    # h = gap: a point on the facing surface is covered by its own layer (d = 0) and the neighbour's (d = L, inclusive)
    @test brute_rate_profile(1.0, 0.0; R=1.0, s=2.5, inner=nothing, outer=(shape="step", h=0.5, L=0.5, pre=1.0)) == 2.0
    @test classify_overlap(2.5, 2.5, 1.0, 0.3) == (true, true, 2)
    @test classify_overlap(2.5, 4.0, 1.0, 0.3) == (true, false, 2)
    @test classify_overlap(4.0, 4.0, 1.0, 0.5) == (false, false, 1)
    @test classify_overlap(4.0, 4.0, 1.0, 2.0)[3] == 4
    @test max_overlap_count(4.0, 4.0, 1.0, 2.0) == 4
    @test max_overlap_count(2.5, 2.5, 1.0, 0.5) == 2
    @test max_overlap_count(4.0, 4.0, 1.0, 0.5) == 1
end
