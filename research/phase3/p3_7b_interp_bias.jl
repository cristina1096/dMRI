# P3.7b (Group C review I3, final review minor): how much does reading a fixed-ρ threshold off three decade points
# (ρR/D ∈ {0.1, 1, 10}, as the Monte Carlo maps do) bias it against a dense ρR/D grid? Curve-shape check on the s = 4
# annulus (R = 1, D = 3, no bulk; 41 readouts to ln 100/λ₁,conv); the annulus is used here BEYOND its validated range
# (ρR/D ≤ 0.1) only as a smooth test curve, not as a result.
# Bias = threshold(decade points) / threshold(dense, 49 log points 0.03–30) − 1, both with `threshold_status` at 0.01.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_7b_interp_bias.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_annulus.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
const OUT = phase3_outdir("p3_7b_interp_bias")
const R1, D3, S4 = 1.0, 3.0, 4.0
gap = S4 - 2R1
function maxdS(sh, hr, x)
    ρ = x * D3 / R1
    λc = annulus_rates(R=R1, R_out=krogh_radius(S4), D=D3, surface=ρ)[1]
    t = collect(range(0.0, log(100) / λc, length=41))
    Sc = annulus_signal(R=R1, R_out=krogh_radius(S4), D=D3, times=t, surface=ρ)
    Sl = annulus_signal(R=R1, R_out=krogh_radius(S4), D=D3, times=t, layer=(sh, ρ, hr, sh == "exponential" ? gap : hr))
    return maximum(abs.(Sl .- Sc))
end
dense_x = exp.(range(log(0.03), log(30.0), length=49))
jobs = [(sh, hr) for sh in ("step", "linear", "exponential") for hr in (0.1, 0.25, 0.5)]
rows = Vector{NamedTuple}(undef, length(jobs))
Threads.@threads for j in eachindex(jobs)
    (sh, hr) = jobs[j]
    (td, sd) = threshold_status(dense_x, [maxdS(sh, hr, x) for x in dense_x], 0.01)
    (t3, s3) = threshold_status([0.1, 1.0, 10.0], [maxdS(sh, hr, x) for x in (0.1, 1.0, 10.0)], 0.01)
    rows[j] = (shape=sh, hR=hr, threshold_dense=td, status_dense=String(sd), threshold_decade=t3, status_decade=String(s3),
               bias=(sd == :crossed && s3 == :crossed) ? t3 / td - 1 : NaN)
end
write_csv(joinpath(OUT, "interp_bias.csv"), rows)
foreach(r -> @printf("%-12s h/R %-5g dense %.4f (%s)  decade %.4f (%s)  bias %+.1f %%\n", r.shape, r.hR, r.threshold_dense, r.status_dense,
    r.threshold_decade, r.status_decade, 100r.bias), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "max_bias" => maximum(abs(r.bias) for r in rows if !isnan(r.bias)),
    "note" => "annulus at s = 4 used beyond its validated range as a smooth test curve only"))
println("saved to $OUT")
