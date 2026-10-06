# P2.6.1 — MCMR at selected τ vs the toy scan (h = 0.1 µm, ΔR2(0) = 0.1, D = 3, T = 10 ms, seeds 1–10).
# Pass: |R_MCMR − R_toy| ≤ 2 combined SEM at every τ in TAUS (independent random streams, so the check is statistical).
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_6_1_mcmr_check.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_6_1_mcmr_check")
const TAUS = [1e-1, 1e-2, 1e-3, 1e-4]           # grid points of the V1 scan that divide T, so toy and MCMR simulate the same T
const T = 10.0
const NSPINS = parse(Int, get(ENV, "NSPINS", "100000"))
# scan.csv columns: tau_ms, R, R_sem, …; keyed by τ rounded to 6 significant digits so grid values match the literals above
scan = Dict{Float64, Tuple{Float64, Float64}}()
for l in readlines(joinpath(REPO_ROOT, "research", "results", "phase2", "p2_6_scan_V1_h0.1", "scan.csv"))[2:end]
    c = parse.(Float64, split(l, ",")[1:3])
    scan[round(c[1], sigdigits=6)] = (c[2], c[3])
end

rows = NamedTuple[]
for τ in TAUS
    r = run_walls(layer_walls(h=0.1, rate=0.1); timestep=τ, times=[0.0, T], nspins=NSPINS)
    Rs = -log.(r.S[:, end]) ./ T                  # r.S is seeds × readouts; the readout is at exactly T
    Rm, Rse = mean(Rs), std(Rs) / sqrt(length(Rs))
    (Rt, Rtse) = scan[round(τ, sigdigits=6)]
    z = (Rm - Rt) / sqrt(Rse^2 + Rtse^2)
    push!(rows, (tau_ms=τ, R_mcmr=Rm, R_mcmr_sem=Rse, R_toy=Rt, R_toy_sem=Rtse, z=z, runtime_s=r.runtime))
    @printf("τ = %.1e: MCMR %.6e ± %.1e, toy %.6e ± %.1e, z = %+.2f\n", τ, Rm, Rse, Rt, Rtse, z)
end
pass = all(abs(r.z) <= 2 for r in rows)
@printf("MCMR vs toy → %s\n", pass ? "PASS" : "FAIL")
write_csv(joinpath(OUT, "mcmr_vs_toy.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rows" => rows, "pass" => pass))
println("saved to $OUT")
