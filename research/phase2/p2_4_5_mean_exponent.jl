# P2.4.5 — mean-exponent identity for every profile, full size: E[−ln M⊥(t)] = 2ρt/w (R2_bulk = 0).
# Holds exactly in expectation for any profile, h and τ whenever each layer's support lies inside the gap:
# the uniform density is stationary under the reflected walk, also at every point along a straight piece.
# Pass: |z| ≤ 3 at t = 50 ms for every configuration, z = (mean over seeds − 2ρt/w) / SEM over seeds.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_4_5_mean_exponent.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_4_5_mean_exponent")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(NSPINS_REF)))
const T = 50.0
const CONFIGS = [("step", 0.2), ("linear", 0.2), ("linear", 1.5), ("exponential", 0.2), ("exponential", 0.5)]

rows = NamedTuple[]
for (shape, h) in CONFIGS, τ in (1e-2, 4e-2)
    ρ = rho_for(RATE_REF, h; shape=shape)
    r = run_walls(layer_walls(h=h, rate=RATE_REF, shape=shape); timestep=τ, times=[0.0, T], nspins=NSPINS, snapshots=true)
    per_seed = [mean(-log(s.orientations[1].transverse) for s in snap.spins) for snap in r.snapshots]
    m, se = mean(per_seed), std(per_seed) / sqrt(length(per_seed))
    expected = 2ρ * T / W_REF
    z = (m - expected) / se
    push!(rows, (shape=shape, h_um=h, tau_ms=τ, rho=ρ, expected=expected, mean=m, sem=se, z=z, pass=abs(z) <= 3))
    @printf("%-11s h = %.2f, τ = %.0e: mean −ln M = %.6f ± %.1e, 2ρt/w = %.6f, z = %+.2f → %s\n",
        shape, h, τ, m, se, expected, z, abs(z) <= 3 ? "PASS" : "FAIL")
end
write_csv(joinpath(OUT, "mean_exponent.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "nspins" => NSPINS, "T_ms" => T,
    "rate" => RATE_REF, "rows" => rows, "pass" => all(r.pass for r in rows)))
println("saved to $OUT")
