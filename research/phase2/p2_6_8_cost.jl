# P2.6.8 runtime vs τ with the layer (vs B5, old code); P2.6.10 runtime vs h and vs ΔR2(0) at the timestep the
# new default constraint gives (τ ≤ 0.005/ΔR2(0), capped at the tortuosity default 0.04 ms); P2.6.11 noise cost at
# fixed budget; P2.6.12 accuracy–cost figure; and the consequence for the earlier τ = 1e-2 ms results
# (ΔR2(0) = 0.1, so ΔR2(0)·τ = 1e-3), using the measured bias at τ = 1e-2 in every toy scan.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_6_8_cost.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_6_8_cost")
const C_LAYER = mr.TimeSteps.DEFAULT_LAYER_SCALING
const TAU_DEFAULT = 0.03 * W_REF^2 / D_REF
R = joinpath(REPO_ROOT, "research", "results")
b5 = Dict(parse(Float64, split(l, ",")[2]) => parse(Float64, split(l, ",")[9])
    for l in readlines(joinpath(R, "baseline", "b5_runtime_vs_tau", "runtime.csv"))[2:end] if split(l, ",")[1] == "0")
us_per_step(r, nspins, T, τ) = r.runtime / (nspins * round(T / τ)) * 1e6
warmup() = run_walls(layer_walls(h=0.1, rate=0.1); timestep=1e-2, times=[0.0, 0.2], seeds=1:1, nspins=1000)

warmup()
rows = NamedTuple[]
# (a) runtime per spin-step with the layer (h = 0.1, ΔR2(0) = 0.1) vs B5 (no layer); 10 000 spins × 20 000 steps
for τ in (1e-1, 1e-2, 1e-3, 1e-4)
    r = run_walls(layer_walls(h=0.1, rate=0.1); timestep=τ, times=[0.0, 20_000τ], seeds=1:1, nspins=10_000)
    us = us_per_step(r, 10_000, 20_000τ, τ)
    push!(rows, (kind="tau", tau_ms=τ, h_um=0.1, rate=0.1, us_per_spin_step=us, b5_us_per_spin_step=b5[τ], overhead=us / b5[τ], hours_full_run=NaN))
    @printf("τ = %.0e: %.4f µs/spin-step with layer, %.4f without (B5) → ×%.2f\n", τ, us, b5[τ], us / b5[τ])
end
# (b) runtime vs h at τ = 1e-2 (ΔR2(0) = 0.1, constraint not binding) and (c) vs ΔR2(0) at the constrained τ (h = 0.1);
# 10 000 spins, T = 10 ms, 1 seed; full run = 1e5 spins × 50 ms × 10 seeds
full_hours(us, τ) = us * 1e5 * (50 / τ) * 10 / 3.6e9
for h in (0.01, 0.05, 0.1, 0.2, 0.5)
    r = run_walls(layer_walls(h=h, rate=0.1); timestep=1e-2, times=[0.0, 10.0], seeds=1:1, nspins=10_000)
    us = us_per_step(r, 10_000, 10.0, 1e-2)
    push!(rows, (kind="h", tau_ms=1e-2, h_um=h, rate=0.1, us_per_spin_step=us, b5_us_per_spin_step=b5[1e-2], overhead=us / b5[1e-2], hours_full_run=full_hours(us, 1e-2)))
    @printf("h = %.2f, τ = 1e-2: %.4f µs/spin-step (×%.2f vs B5), full run ≈ %.2f h\n", h, us, us / b5[1e-2], full_hours(us, 1e-2))
end
for rate in (0.1, 0.5, 1.0, 2.0, 10.0)
    τ = min(C_LAYER / rate, TAU_DEFAULT)
    r = run_walls(layer_walls(h=0.1, rate=rate); timestep=τ, times=[0.0, 10.0], seeds=1:1, nspins=10_000)
    us = us_per_step(r, 10_000, 10.0, τ)
    push!(rows, (kind="rate", tau_ms=τ, h_um=0.1, rate=rate, us_per_spin_step=us, b5_us_per_spin_step=NaN, overhead=NaN, hours_full_run=full_hours(us, τ)))
    @printf("ΔR2(0) = %5.1f: τ = %.2e ms, %.4f µs/spin-step, full run ≈ %.2f h\n", rate, τ, us, full_hours(us, τ))
end

# consequences: measured bias at τ = 1e-2 (and at the default 0.04-ish grid point 4.64e-2) in every toy scan
scanrow(l, τ) = (s = [split(x, ",") for x in readlines(joinpath(R, "phase2", "p2_6_scan_" * l, "scan.csv"))[2:end]];
    τs = parse.(Float64, first.(s)); Rv = [parse(Float64, c[2]) for c in s]; Rv[argmin(abs.(log.(τs ./ τ)))] / Rv[argmin(τs)] - 1)
cons = NamedTuple[]
for l in ("V2_h0.01", "V2_h0.05", "V1_h0.1", "V2_h0.2", "V3_D1", "V3_D6", "V2_rate0.5", "V2_rate1", "V2_rate2", "V2_h0.01_rate10")
    s = JSON.parsefile(joinpath(R, "phase2", "p2_6_scan_" * l, "summary.json"))
    push!(cons, (scan=l, h_um=s["h_um"], rate=s["rate"], D=s["D"], rate_tau_at_1e_2=s["rate"] * 1e-2,
        bias_at_1e_2=scanrow(l, 1e-2), bias_at_4p6e_2=scanrow(l, 4.6415888336127795e-2)))
    @printf("%-16s ΔR2(0)·τ(1e-2) = %.3f: bias at τ = 1e-2 %+.1e, at 4.6e-2 %+.1e\n", l, s["rate"] * 1e-2, cons[end].bias_at_1e_2, cons[end].bias_at_4p6e_2)
end
write_csv(joinpath(OUT, "runtime.csv"), rows)
write_csv(joinpath(OUT, "consequences.csv"), cons)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "c_layer" => C_LAYER, "tau_default" => TAU_DEFAULT,
    "rows" => [map(x -> x isa Float64 && isnan(x) ? nothing : x, r) for r in rows], "consequences" => cons))

using CairoMakie
fig = Figure(size=(600, 400))
ax = Axis(fig[1, 1], xscale=log10, yscale=log10, xlabel="runtime per spin per ms of sequence (µs)", ylabel="|R/R_ref − 1|",
    title="accuracy–cost (toy error, MCMR runtime with layer)")
us_step = mean(r.us_per_spin_step for r in rows if r.kind == "tau")
for l in ("V1_h0.1", "V2_rate1", "V2_rate2")
    s = [split(x, ",") for x in readlines(joinpath(R, "phase2", "p2_6_scan_" * l, "scan.csv"))[2:end]]
    τs = parse.(Float64, first.(s)); Rv = [parse(Float64, c[2]) for c in s]
    err = abs.(Rv ./ Rv[argmin(τs)] .- 1)
    keep = (err .> 0) .& (τs .>= 1e-4)
    scatterlines!(ax, (us_step ./ τs)[keep], err[keep], label=l)
end
hlines!(ax, [0.01, 0.001], color=:gray, linestyle=:dash)
axislegend(ax, position=:rt, labelsize=9)
save(joinpath(OUT, "cost.png"), fig)
println("saved to $OUT")
