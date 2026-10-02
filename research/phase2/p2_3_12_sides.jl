# P2.3.12 — side-specific layer.
# (a) Exact mirror check (pass/fail): a positive-side layer on walls with normal −x is geometrically the same as a
#     negative-side layer on walls with normal +x; with the same seeds the signals must be identical (==) and the
#     positions identical. This tests that the code treats both sides symmetrically, without statistics.
# (b) Statistical mirror comparison (reported, not graded): positive-only vs negative-only (h = 0.2 µm, ΔR2(0) = 0.1)
#     on +x walls, same seeds: paired per-seed differences d_s(t) = S₊ − S₋, z = mean(d)/SEM(d) per readout.
#     (The first version graded |z| ≤ 2 with unpaired SEMs at all 11 readouts; that criterion was mis-specified —
#     see phase2_record.md.)
# (c) Characterisation: both sides, one side, and asymmetric (positive h = 0.2, ΔR2(0) = 0.1; negative h = 0.4,
#     ΔR2(0) = 0.05 → same ρ = 0.02 µm/ms on both sides) — S(t) with SEM. No baseline exists for asymmetric surfaces.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_3_12_sides.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "fit.jl"))

const OUT = phase2_outdir("p2_3_12_sides")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(sweep_nspins())))

cases = Dict(
    "both" => layer_walls(h=0.2, rate=RATE_REF),
    "positive" => layer_walls(h=0.2, rate=RATE_REF, side=:positive),
    "negative" => layer_walls(h=0.2, rate=RATE_REF, side=:negative),
    "asymmetric" => layer_walls(h=0.2, rate=RATE_REF, side=:asymmetric, h_neg=0.4, rate_neg=0.05),
)
runs = Dict(k => run_walls(v; nspins=NSPINS) for (k, v) in cases)
flipped = run_walls(Walls(repeats=W_REF, rotation=[-1.0, 0.0, 0.0], layer_rho_positive=rho_for(RATE_REF, 0.2),
    layer_h_positive=0.2); nspins=NSPINS)
exact_mirror = flipped.S == runs["negative"].S
@printf("exact mirror (positive on −x walls == negative on +x walls): %s → %s\n", exact_mirror, exact_mirror ? "PASS" : "FAIL")
p, n = runs["positive"], runs["negative"]
d = p.S .- n.S
dmean = vec(mean(d, dims=1)); dsem = vec(std(d, dims=1)) ./ sqrt(size(d, 1))
z = [dsem[i] > 0 ? dmean[i] / dsem[i] : 0.0 for i in eachindex(dsem)]
@printf("statistical mirror (paired, reported only): max |z| = %.2f, S₊(50) − S₋(50) = %.2e ± %.1e\n", maximum(abs.(z)), dmean[end], dsem[end])
for k in ("both", "positive", "negative", "asymmetric")
    @printf("%-10s S(50) = %.5f ± %.5f\n", k, runs[k].mean[end], runs[k].sem[end])
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "nspins" => NSPINS, "times" => TIMES_REF,
    "exact_mirror_pass" => exact_mirror, "paired_diff_mean" => dmean, "paired_diff_sem" => dsem, "paired_z" => z, "paired_max_abs_z" => maximum(abs.(z)),
    "cases" => Dict(k => Dict("mean" => r.mean, "sem" => r.sem) for (k, r) in runs)))
println("saved to $OUT")
