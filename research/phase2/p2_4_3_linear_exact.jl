# P2.4.3 — exact checks for the linear profile, full size (seeds 1–10, NSPINS spins), τ = 1e-2 ms.
#
# V0b       : linear, ρ = 0 (h = 0.5) vs plain walls, R2_bulk = 1/80: positions and per-spin M⊥ bit-identical.
# uniform   : linear, h = w = 2 µm, both faces, ΔR2(0) = 0.1 → ΔR2(u) = 2ρ/w = 0.1 everywhere in the gap, so every spin
#             decays at exactly 0.1 (+ R2_bulk): per spin ≤ 1e-10, ensemble ≤ roundoff bound. Bulk off and on.
# additivity: linear, h = 0.2, R2_bulk 0 vs 1/80, same seeds: positions ==, per spin ≤ 1e-10, ensemble ≤ roundoff bound.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_4_3_linear_exact.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "expected.jl"))

const OUT = phase2_outdir("p2_4_3_linear_exact")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(NSPINS_REF)))
const TOL_SPIN = 1e-10
per_spin_final(r) = [[s.orientations[1].transverse for s in snap.spins] for snap in r.snapshots]
positions_final(r) = [mr.position.(snap) for snap in r.snapshots]

results = Dict{String, Any}()
function report(name, pass, extra)
    results[name] = merge(Dict{String, Any}("pass" => pass), extra)
    @printf("%-18s → %s  %s\n", name, pass ? "PASS" : "FAIL", string(extra))
end

a = run_walls(Walls(repeats=W_REF); R2_bulk=1 / 80, nspins=NSPINS, snapshots=true)
b = run_walls(layer_walls(h=0.5, rate=0.0, shape="linear"); R2_bulk=1 / 80, nspins=NSPINS, snapshots=true)
report("V0b", positions_final(a) == positions_final(b) && per_spin_final(a) == per_spin_final(b) && a.S == b.S, Dict("R2_bulk" => 1 / 80))

bound = roundoff_bound(TIMES_REF, TAU_REF, NSPINS)
for (label, R2) in (("uniform_bulk_off", 0.0), ("uniform_bulk_on", 1 / 80))
    r = run_walls(layer_walls(h=W_REF, rate=RATE_REF, shape="linear"); R2_bulk=R2, nspins=NSPINS, snapshots=true)
    target = expected_half_gap(TIMES_REF, RATE_REF; R2_bulk=R2)          # same formula: uniform rate RATE_REF + R2
    spin_dev = maximum(maximum(abs.(m ./ target[end] .- 1)) for m in per_spin_final(r))
    ens_dev = maximum(abs.(r.S ./ target' .- 1))
    report(label, spin_dev <= TOL_SPIN && ens_dev <= bound,
        Dict("max_rel_dev_per_spin" => spin_dev, "max_rel_dev_ensemble" => ens_dev, "bound" => bound, "R2_bulk" => R2))
end

lo = run_walls(layer_walls(h=0.2, rate=RATE_REF, shape="linear"); nspins=NSPINS, snapshots=true)
wb = run_walls(layer_walls(h=0.2, rate=RATE_REF, shape="linear"); R2_bulk=1 / 80, nspins=NSPINS, snapshots=true)
same_pos = positions_final(lo) == positions_final(wb)
spin_dev = maximum(maximum(abs.(x ./ (y .* exp(-TIMES_REF[end] / 80)) .- 1)) for (x, y) in zip(per_spin_final(wb), per_spin_final(lo)))
ens_dev = maximum(abs.(wb.S ./ (lo.S .* exp.(-TIMES_REF' ./ 80)) .- 1))
report("additivity", same_pos && spin_dev <= TOL_SPIN && ens_dev <= bound,
    Dict("positions_identical" => same_pos, "max_rel_dev_per_spin" => spin_dev, "max_rel_dev_ensemble" => ens_dev, "bound" => bound))

write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "shape" => "linear", "nspins" => NSPINS, "results" => results))
println("saved to $OUT")
