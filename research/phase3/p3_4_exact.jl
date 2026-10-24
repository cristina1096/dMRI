# P3.2.5 V0b, P3.4.2 uniform inner layer h = R, P3.4.4 additivity, P3.4.7 side-specific. Exact in any diffusion regime.
# Seeds 1–10, 100000 spins, readouts 0:5:50, τ = 1e-2 (V0b: default τ).
# - V0b: ρ = 0, h = 0.3 vs no layer (square array, R2 = 1/80): positions and every M⊥ identical (==).
# - Uniform: single cylinder, inner step h = R = 1, ρ = 0.05 → κ = 1/2, ΔR2 = 2ρ/R = 0.1 /ms everywhere inside, so every
#   inside spin has M⊥ = e^(−(0.1 + R2)t): per spin ≤ 1e-10, ensemble ≤ roundoff_bound; bulk off and on (1/80).
# - Additivity: square array ρ = 0.02, h = 0.2, R2 0 vs 1/80: positions identical, per spin ratio e^(−t/80) ≤ 1e-10.
# - Side-specific: inner-only layer leaves every outside spin at M⊥ = 1, outer-only every inside spin, up to roundoff
#   (1 − M⊥ ≤ 1e-9: reflection points sit a roundoff distance across r = R, and the ≈ 1e-15 per reflection accumulates to
#   ≈ 1e-12 over 50 ms; the layer's own effect in the other compartment is 1 − S ≈ 0.3–0.6).
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_exact.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_4_exact")
transverse(f, m) = [sp.orientations[1].transverse for (k, sp) in enumerate(f.spins) if m[k]]
bb = mr.BoundingBox([-R_REF, -R_REF, -R_REF], [R_REF, R_REF, R_REF])
res = Dict{String, Any}()
report(name, pass, d) = (res[name] = merge(Dict{String, Any}(d), Dict{String, Any}("pass" => pass)); @printf("%-22s %s\n", name, pass ? "PASS" : "FAIL"))

# V0b (both regions)
for region in (:inside, :outside)
    kw = (region=region, sx=S_REF, sy=S_REF, R2_bulk=1 / 80, timestep=nothing)
    a = run_cylinders(Cylinders(radius=R_REF, repeats=[S_REF, S_REF]); kw...)
    b = run_cylinders(Cylinders(radius=R_REF, repeats=[S_REF, S_REF], layer_h=0.3); kw...)
    same = all(mr.position.(x) == mr.position.(y) && [s.orientations[1].transverse for s in x.spins] == [s.orientations[1].transverse for s in y.spins]
               for (x, y) in zip(a.finals, b.finals))
    report("V0b_$(region)", same, Dict("identical" => same))
end
# uniform inner layer
for R2 in (0.0, 1 / 80)
    r = run_cylinders(Cylinders(radius=R_REF, layer_rho_inner=0.05, layer_h_inner=R_REF); region=:inside, bounding_box=bb, R2_bulk=R2)
    target = exp.(-(0.1 + R2) .* TIMES_REF)
    spin_dev = maximum(maximum(abs.(transverse(f, m) ./ target[end] .- 1)) for (f, m) in zip(r.finals, r.masks))
    ens_dev = maximum(abs.(r.mean ./ target .- 1))
    bound = roundoff_bound(TIMES_REF, r.tau, NSPINS_REF)
    report("uniform_R2_$(R2)", spin_dev <= 1e-10 && ens_dev <= bound,
        Dict("max_dev_per_spin" => spin_dev, "max_dev_ensemble" => ens_dev, "bound" => bound))
end
# additivity
g = Cylinders(radius=R_REF, repeats=[S_REF, S_REF], layer_rho=0.02, layer_h=0.2)
for region in (:inside, :outside)
    a = run_cylinders(g; region=region, sx=S_REF, sy=S_REF)
    b = run_cylinders(g; region=region, sx=S_REF, sy=S_REF, R2_bulk=1 / 80)
    same = all(mr.position.(x) == mr.position.(y) for (x, y) in zip(a.finals, b.finals))
    dev = maximum(maximum(abs.(transverse(fb, m) ./ (transverse(fa, m) .* exp(-50 / 80)) .- 1)) for (fa, fb, m) in zip(a.finals, b.finals, a.masks))
    report("additivity_$(region)", same && dev <= 1e-10, Dict("positions_identical" => same, "max_dev_per_spin" => dev,
        "S50_layer_only" => a.mean[end]))
end
# side-specific
for (label, kw, untouched) in (("inner_only", (layer_rho_inner=0.02, layer_h_inner=0.2), :outside),
                               ("outer_only", (layer_rho_outer=0.02, layer_h_outer=0.2), :inside))
    gs = Cylinders(radius=R_REF, repeats=[S_REF, S_REF]; kw...)
    u = run_cylinders(gs; region=untouched, sx=S_REF, sy=S_REF)
    acted = run_cylinders(gs; region=(untouched == :inside ? :outside : :inside), sx=S_REF, sy=S_REF)
    leak = maximum(maximum(1 .- transverse(f, m)) for (f, m) in zip(u.finals, u.masks))
    exact = leak <= 1e-9          # roundoff at reflection points on r = R, accumulated over 50 ms (ledger rulings, Tasks 4 and 11)
    report("sides_$(label)", exact && acted.mean[end] < 1, Dict("untouched_max_leak" => leak, "untouched_within_bound" => exact, "S50_acted" => acted.mean[end]))
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "results" => res, "pass" => all(v["pass"] for v in values(res))))
println("saved to $OUT")
