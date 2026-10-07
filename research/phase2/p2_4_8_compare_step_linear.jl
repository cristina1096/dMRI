# P2.4.8 (step vs linear) — profile comparison at fixed surface excess rate ΔR2(0) = 0.1 /ms. Characterisation only.
# (1) Attenuation 1 − S(50) vs h for both profiles (from the two sweeps).
# (2) V4: h*(θ) and ρ*(θ) = ΔR2(0)·h*/g(0) for both profiles, their ratios, and the max |residual| of each fit.
#
# Run: julia --project=research/baseline research/phase2/p2_4_8_compare_step_linear.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_4_8_compare_step_linear")
R = joinpath(REPO_ROOT, "research", "results", "phase2")
sweep(sfx) = JSON.parsefile(joinpath(R, "p2_3_4_h_sweep" * sfx, "summary.json"))["per_h"]
fits(sfx) = JSON.parsefile(joinpath(R, "p2_3_9_fit_h" * sfx, "summary.json"))["fits"]

ss, sl = sweep(""), sweep("_linear")
fs, fl = fits(""), fits("_linear")
rows = NamedTuple[]
for (a, b) in zip(fs, fl)
    @assert a["theta"] == b["theta"]
    hs, hl = a["h_final"], b["h_final"]
    ρs, ρl = RATE_REF * hs / 1, RATE_REF * hl / 2
    push!(rows, (theta=a["theta"], h_star_step=hs, h_star_linear=hl, h_ratio=hl / hs,
        rho_star_step=ρs, rho_star_linear=ρl, rho_ratio=ρl / ρs,
        max_abs_residual_step=maximum(abs.(a["residuals_final"])), max_abs_residual_linear=maximum(abs.(b["residuals_final"]))))
    @printf("θ = %.4f: h* step %.5f, linear %.5f (ratio %.4f); ρ* ratio %.4f; max |res| %.1f / %.1f\n",
        a["theta"], hs, hl, hl / hs, ρl / ρs, rows[end].max_abs_residual_step, rows[end].max_abs_residual_linear)
end
write_csv(joinpath(OUT, "comparison.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rate" => RATE_REF, "rows" => rows))

using CairoMakie
fig = Figure(size=(1000, 360))
ax1 = Axis(fig[1, 1], xlabel="h (µm)", ylabel="1 − S(50)", title="Attenuation vs h at ΔR₂(0) = 0.1 /ms")
scatterlines!(ax1, [p["h_um"] for p in ss], [p["attenuation_50"] for p in ss], label="step")
scatterlines!(ax1, [p["h_um"] for p in sl], [p["attenuation_50"] for p in sl], label="linear")
axislegend(ax1, position=:rb)
ax2 = Axis(fig[1, 2], xlabel="θ_relax", ylabel="ratio linear / step", xscale=log10, title="V4: h* and ρ* ratios")
scatterlines!(ax2, [r.theta for r in rows], [r.h_ratio for r in rows], label="h* ratio")
scatterlines!(ax2, [r.theta for r in rows], [r.rho_ratio for r in rows], label="ρ* ratio")
axislegend(ax2, position=:rt)
save(joinpath(OUT, "p2_4_8_compare.png"), fig)
println("saved to $OUT")
