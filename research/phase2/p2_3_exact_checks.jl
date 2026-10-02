# P2.3.1 V0b, P2.3.4(a) exact end point, P2.3.11 additivity — full size (seeds 1–10, NSPINS spins).
#
# V0b       : rho = 0 (layer_h = 0.5) vs plain Walls(repeats=2), default τ (as B1/B2), R2_bulk = 0 and 1/80:
#             final positions and per-spin M⊥ bit-identical (==), every seed.
# h = w/2   : h = 1.0 µm, ΔR2(0) = 0.1, τ = 1e-2, R2_bulk = 0 and 1/80: per spin |M⊥/exp(-(ΔR2(0)+R2_bulk)·50) − 1| ≤ 1e-10;
#             ensemble at every readout within roundoff_bound.
# additivity: h = 0.2 µm, ΔR2(0) = 0.1, τ = 1e-2, R2_bulk = 0 vs 1/80, same seeds: positions ==, per spin
#             |M⊥_bulk / (M⊥_layer·exp(-50/80)) − 1| ≤ 1e-10; ensemble within roundoff_bound.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_3_exact_checks.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "expected.jl"))

const OUT = phase2_outdir("p2_3_exact_checks")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(NSPINS_REF)))
const TOL_SPIN = 1e-10
per_spin_final(r) = [[s.orientations[1].transverse for s in snap.spins] for snap in r.snapshots]
positions_final(r) = [mr.position.(snap) for snap in r.snapshots]

results = Dict{String, Any}()
function report(name, pass, extra)
    results[name] = merge(Dict{String, Any}("pass" => pass), extra)
    @printf("%-18s → %s  %s\n", name, pass ? "PASS" : "FAIL", string(extra))
end

for (label, R2) in (("bulk_off", 0.0), ("bulk_on", 1 / 80))
    a = run_walls(Walls(repeats=W_REF); R2_bulk=R2, timestep=nothing, nspins=NSPINS, snapshots=true)
    b = run_walls(layer_walls(h=0.5, rate=0.0); R2_bulk=R2, timestep=nothing, nspins=NSPINS, snapshots=true)
    same = positions_final(a) == positions_final(b) && per_spin_final(a) == per_spin_final(b) && a.S == b.S
    report("V0b_" * label, same, Dict("tau_ms" => a.tau, "R2_bulk" => R2))
end

bound = roundoff_bound(TIMES_REF, TAU_REF, NSPINS)
for (label, R2) in (("bulk_off", 0.0), ("bulk_on", 1 / 80))
    r = run_walls(layer_walls(h=1.0, rate=RATE_REF); R2_bulk=R2, nspins=NSPINS, snapshots=true)
    target = expected_half_gap(TIMES_REF, RATE_REF; R2_bulk=R2)
    spin_dev = maximum(maximum(abs.(m ./ target[end] .- 1)) for m in per_spin_final(r))
    ens_dev = maximum(abs.(r.S ./ target' .- 1))
    report("half_gap_" * label, spin_dev <= TOL_SPIN && ens_dev <= bound,
        Dict("max_rel_dev_per_spin" => spin_dev, "max_rel_dev_ensemble" => ens_dev, "bound" => bound, "R2_bulk" => R2))
end

layer_only = run_walls(layer_walls(h=0.2, rate=RATE_REF); nspins=NSPINS, snapshots=true)
with_bulk = run_walls(layer_walls(h=0.2, rate=RATE_REF); R2_bulk=1 / 80, nspins=NSPINS, snapshots=true)
same_pos = positions_final(layer_only) == positions_final(with_bulk)
spin_dev = maximum(maximum(abs.(mb ./ (ml .* exp(-TIMES_REF[end] / 80)) .- 1))
    for (mb, ml) in zip(per_spin_final(with_bulk), per_spin_final(layer_only)))
ens_dev = maximum(abs.(with_bulk.S ./ (layer_only.S .* exp.(-TIMES_REF' ./ 80)) .- 1))
report("additivity", same_pos && spin_dev <= TOL_SPIN && ens_dev <= bound,
    Dict("positions_identical" => same_pos, "max_rel_dev_per_spin" => spin_dev, "max_rel_dev_ensemble" => ens_dev,
        "bound" => bound, "S50_layer_only" => layer_only.mean[end]))

write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "nspins" => NSPINS, "results" => results))
println("saved to $OUT")
