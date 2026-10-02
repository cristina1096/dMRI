# P2.3.7 V11 — exact segment overlap vs endpoint rule, measured (no assumed sign of the bias).
#
# 1D toy (toy.jl), h = 0.1 µm, ΔR2(0) = 0.1 /ms, T = 50 ms, τ ∈ TAUS, seeds 1–10.
# (1) bias = S_endpoint − S_exact at T, in combined SEM, per τ (characterisation).
# (2) cross-check: toy :exact at τ = 1e-2 vs MCMR sweep (P2.3.4(b)) at h = 0.1, S(50) within 2 combined SEM.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_3_7_endpoint.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "fit.jl"))
include(joinpath(@__DIR__, "toy.jl"))

const OUT = phase2_outdir("p2_3_7_endpoint")
const H = 0.1
const TAUS = [1e-3, 1e-2, 4e-2]
const NSPINS = Dict(1e-3 => parse(Int, get(ENV, "NSPINS_1E3", "20000")), 1e-2 => 100_000, 4e-2 => 100_000)

seedmeans(τ, m) = [mean(toy_layer_signal(; h=H, rate=RATE_REF, tau=τ, nspins=NSPINS[τ], seed=s, method=m)) for s in 1:10]
rows = NamedTuple[]
for τ in TAUS
    ex = seedmeans(τ, :exact); ep = seedmeans(τ, :endpoint)
    se(v) = std(v) / sqrt(length(v))
    z = (mean(ep) - mean(ex)) / sqrt(se(ex)^2 + se(ep)^2)
    push!(rows, (tau_ms=τ, nspins=NSPINS[τ], S_exact=mean(ex), S_exact_sem=se(ex), S_endpoint=mean(ep), S_endpoint_sem=se(ep), bias=mean(ep) - mean(ex), bias_in_sem=z))
    @printf("τ = %.0e: S_exact %.5f ± %.5f, S_endpoint %.5f ± %.5f, bias %+.2e (%+.1f SEM)\n", τ, mean(ex), se(ex), mean(ep), se(ep), mean(ep) - mean(ex), z)
end
(hgrid, S_by_h) = read_sweep(joinpath(REPO_ROOT, "research", "results", "phase2", "p2_3_4_h_sweep", "sweep.csv"))
k = findfirst(==(H), hgrid)
mc = S_by_h[k][:, end]
toy = rows[findfirst(r -> r.tau_ms == 1e-2, rows)]
zc = (toy.S_exact - mean(mc)) / sqrt(toy.S_exact_sem^2 + (std(mc) / sqrt(length(mc)))^2)
cross_ok = abs(zc) <= 2
@printf("cross-check toy :exact vs MCMR at h = %.1f, τ = 1e-2: %.5f vs %.5f (%+.1f SEM) → %s\n", H, toy.S_exact, mean(mc), zc, cross_ok ? "PASS" : "FAIL")
write_csv(joinpath(OUT, "endpoint_vs_exact.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "h_um" => H, "rate" => RATE_REF, "T_ms" => 50.0,
    "rows" => rows, "cross_check_z" => zc, "cross_check_pass" => cross_ok,
    "note" => "toy and MCMR use different random streams; the cross-check is statistical"))
println("saved to $OUT")
