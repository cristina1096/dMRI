# P2.3.4(b) attenuation vs h, P2.3.5 monotonicity, P2.3.6 decay shape.
#
# Profile from SHAPE (default step; linear uses H_GRID_LINEAR), both faces, ΔR2(0) = 0.1 /ms, h ∈ GRID, τ = 1e-2 ms (provisional, see plan), seeds 1–10,
# NSPINS = sweep_nspins() spins per seed. Records S(t) per seed (sweep.csv; also the input of the V4 fit).
# P2.3.5 pass: attenuation 1 − S(50) strictly increasing with h, each step by > 2 combined SEM.
# P2.3.6: curvature c of ln S(t) per seed (mean ± SEM) vs h — characterisation, no pass/fail.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_3_4_h_sweep.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "fit.jl"))

const SHAPE = get(ENV, "SHAPE", "step")
const GRID = SHAPE == "step" ? H_GRID : H_GRID_LINEAR
const OUT = phase2_outdir(SHAPE == "step" ? "p2_3_4_h_sweep" : "p2_3_4_h_sweep_" * SHAPE)
const NSPINS = parse(Int, get(ENV, "NSPINS", string(sweep_nspins())))

rows = NamedTuple[]
per_h = Dict{String, Any}[]
for h in GRID
    r = run_walls(layer_walls(h=h, rate=RATE_REF, shape=SHAPE); nspins=NSPINS)
    for s in 1:size(r.S, 1), (j, t) in enumerate(r.times)
        push!(rows, (h_um=h, seed=s, t_ms=t, signal=r.S[s, j]))
    end
    curv = [curvature(r.times, r.S[s, :]) for s in 1:size(r.S, 1)]
    push!(per_h, Dict("h_um" => h, "mean" => r.mean, "sem" => r.sem, "sigma" => r.sigma,
        "attenuation_50" => 1 - r.mean[end], "attenuation_50_sem" => r.sem[end],
        "curvature_mean" => mean(curv), "curvature_sem" => std(curv) / sqrt(length(curv)), "runtime_s" => r.runtime))
    @printf("h = %.3f: S(50) = %.5f ± %.5f (SEM), curvature %.2e ± %.1e, %.0f s\n",
        h, r.mean[end], r.sem[end], mean(curv), std(curv) / sqrt(length(curv)), r.runtime)
end
att = [p["attenuation_50"] for p in per_h]
sem = [p["attenuation_50_sem"] for p in per_h]
steps = increasing_steps(att, sem)
monotone = all(steps)
@printf("P2.3.5 monotonicity: %d of %d steps increase by > 2 SEM → %s\n", count(steps), length(steps), monotone ? "PASS" : "FAIL")
write_csv(joinpath(OUT, "sweep.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "nspins" => NSPINS, "rate" => RATE_REF,
    "tau_ms" => TAU_REF, "times" => TIMES_REF, "h_grid" => GRID, "shape" => SHAPE, "per_h" => per_h,
    "monotone_steps" => steps, "monotone_pass" => monotone))

using CairoMakie
fig = Figure(size=(1100, 340))
ax1 = Axis(fig[1, 1], xlabel="t (ms)", ylabel="S(t)", yscale=log10, title="S(t) per h, $(SHAPE) (ΔR₂(0) = 0.1 /ms)")
for p in per_h
    lines!(ax1, TIMES_REF, max.(p["mean"], 1e-6), label="h = $(p["h_um"])")
end
Legend(fig[1, 2], ax1, labelsize=8, rowgap=0)
ax2 = Axis(fig[1, 3], xlabel="h (µm)", ylabel="1 − S(50)", title="P2.3.5 attenuation vs h")
errorbars!(ax2, GRID, att, 2 .* sem); scatterlines!(ax2, GRID, att)
ax3 = Axis(fig[1, 4], xlabel="h (µm)", ylabel="curvature c of ln S", title="P2.3.6 decay shape")
cm = [p["curvature_mean"] for p in per_h]; cs = [p["curvature_sem"] for p in per_h]
errorbars!(ax3, GRID, cm, 2 .* cs); scatterlines!(ax3, GRID, cm); hlines!(ax3, [0.0], color=:gray, linestyle=:dash)
save(joinpath(OUT, "p2_3_4_h_sweep.png"), fig)
println("saved to $OUT")
