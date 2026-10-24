# P3.1.1 C-B1 (null run) and P3.1.2 C-B2 (bulk R2 only) on cylinders, old code (no layer).
# Configurations: single cylinder R = 1 (inside), square array s = 2.5, R = 1 (inside, outside). Default timestep
# (recorded: C-B1 asks for τ_max), seeds 1–10, 100000 spins per seed, readouts 0:5:50.
# C-B1: every spin keeps M⊥ = 1 exactly. C-B2 (R2 = 1/80): positions identical to C-B1 (same seeds); per spin
# |M⊥/e^(−t/80) − 1| ≤ 1e-10; ensemble ≤ roundoff_bound.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_1_1_null_bulk.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_1_1_null_bulk")
const R2B = 1 / 80

configs = [
    ("single_inside", Cylinders(radius=R_REF), :inside, mr.BoundingBox([-R_REF, -R_REF, -R_REF], [R_REF, R_REF, R_REF]), Inf),
    ("array_inside", Cylinders(radius=R_REF, repeats=[S_REF, S_REF]), :inside, 500, S_REF),
    ("array_outside", Cylinders(radius=R_REF, repeats=[S_REF, S_REF]), :outside, 500, S_REF),
]
results = Dict{String, Any}()
allpass = true
for (name, geometry, region, bb, s) in configs
    kw = (region=region, bounding_box=bb, sx=s, sy=s, timestep=nothing)
    null = run_cylinders(geometry; kw...)
    bulk = run_cylinders(geometry; R2_bulk=R2B, kw...)
    null_exact = all(all(sp.orientations[1].transverse == 1.0 for sp in f.spins) for f in null.finals)
    same_pos = all(mr.position.(a) == mr.position.(b) for (a, b) in zip(null.finals, bulk.finals))
    target = exp(-R2B * TIMES_REF[end])
    spin_dev = maximum(maximum(abs(sp.orientations[1].transverse / target - 1) for (k, sp) in enumerate(f.spins) if m[k])
                       for (f, m) in zip(bulk.finals, bulk.masks))
    ens_dev = maximum(abs.(bulk.mean ./ exp.(-R2B .* TIMES_REF) .- 1))
    bound = roundoff_bound(TIMES_REF, bulk.tau, NSPINS_REF)
    pass = null_exact && all(null.S .== 1.0) && same_pos && spin_dev <= 1e-10 && ens_dev <= bound
    global allpass &= pass
    results[name] = Dict("tau_default" => null.tau, "null_exact" => null_exact, "positions_identical" => same_pos,
        "max_dev_per_spin" => spin_dev, "max_dev_ensemble" => ens_dev, "bound" => bound,
        "nregion_mean" => mean(null.nregion), "runtime_s" => null.runtime + bulk.runtime, "pass" => pass)
    @printf("%-14s τ = %.4g ms: C-B1 exact %s; C-B2 positions identical %s, per spin %.1e, ensemble %.1e (bound %.1e) → %s\n",
        name, null.tau, null_exact, same_pos, spin_dev, ens_dev, bound, pass ? "PASS" : "FAIL")
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "results" => results, "pass" => allpass,
    "params" => Dict("R" => R_REF, "s" => S_REF, "D" => D_REF, "R2_bulk" => R2B, "seeds" => collect(SEEDS_REF), "nspins" => NSPINS_REF)))
println("saved to $OUT")
