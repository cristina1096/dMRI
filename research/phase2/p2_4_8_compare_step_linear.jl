# P2.4.8 (step vs linear, plus exponential attenuation) — profile comparison at fixed surface excess rate ΔR2(0) = 0.1 /ms. Characterisation only.
# (1) Attenuation 1 − S(50) vs h for both profiles (from the two sweeps).
# (2) V4: h*(θ) and ρ*(θ) = ΔR2(0)·h*/g(0) for both profiles, their ratios, and the max |residual| of each fit.
#     σ of the ρ* ratio: jackknife σ of each fit (parabola stage, `h_star_jackknife_sigma`) added in quadrature.
#     Both fits share seeds 1–10 and the B4 reference, so their errors are correlated and this σ is an
#     overestimate; z = (ratio − 1)/σ is therefore conservative.
#
# Run: julia --project=research/baseline research/phase2/p2_4_8_compare_step_linear.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_4_8_compare_step_linear")
R = joinpath(REPO_ROOT, "research", "results", "phase2")
sweep(sfx) = JSON.parsefile(joinpath(R, "p2_3_4_h_sweep" * sfx, "summary.json"))["per_h"]
fits(sfx) = JSON.parsefile(joinpath(R, "p2_3_9_fit_h" * sfx, "summary.json"))["fits"]

ss, sl, sx = sweep(""), sweep("_linear"), sweep("_exponential")
fs, fl = fits(""), fits("_linear")
rows = NamedTuple[]
for (a, b) in zip(fs, fl)
    @assert a["theta"] == b["theta"]
    hs, hl = a["h_final"], b["h_final"]
    ρs, ρl = RATE_REF * hs / 1, RATE_REF * hl / 2
    σ_ratio = (ρl / ρs) * hypot(a["h_star_jackknife_sigma"] / hs, b["h_star_jackknife_sigma"] / hl)
    push!(rows, (theta=a["theta"], h_star_step=hs, h_star_linear=hl, h_ratio=hl / hs,
        rho_star_step=ρs, rho_star_linear=ρl, rho_ratio=ρl / ρs, rho_ratio_sigma=σ_ratio, rho_ratio_z=(ρl / ρs - 1) / σ_ratio,
        max_abs_residual_step=maximum(abs.(a["residuals_final"])), max_abs_residual_linear=maximum(abs.(b["residuals_final"]))))
    @printf("θ = %.4f: h* step %.5f, linear %.5f (ratio %.4f); ρ* ratio − 1 = %+.1e ± %.0e (z = %+.1f); max |res| %.1f / %.1f\n",
        a["theta"], hs, hl, hl / hs, ρl / ρs - 1, σ_ratio, rows[end].rho_ratio_z, rows[end].max_abs_residual_step, rows[end].max_abs_residual_linear)
end
write_csv(joinpath(OUT, "comparison.csv"), rows)

# (3) Attenuation at fixed ΔR2(0) for all three profiles, at every h present in all three sweeps (no exponential V4).
byh(sw) = Dict(p["h_um"] => p for p in sw)
bs, bl, bx = byh(ss), byh(sl), byh(sx)
att = NamedTuple[]
for h in sort([h for h in keys(bs) if haskey(bl, h) && haskey(bx, h)])
    a, b, c = bs[h], bl[h], bx[h]
    push!(att, (h_um=h, S50_step=a["mean"][end], S50_step_sem=a["sem"][end], S50_linear=b["mean"][end], S50_linear_sem=b["sem"][end],
        S50_exponential=c["mean"][end], S50_exponential_sem=c["sem"][end],
        rho_step=rho_for(RATE_REF, h), rho_linear=rho_for(RATE_REF, h; shape="linear"), rho_exponential=rho_for(RATE_REF, h; shape="exponential")))
    @printf("h = %.3f: S(50) step %.5f, linear %.5f, exponential %.5f\n", h, a["mean"][end], b["mean"][end], c["mean"][end])
end
write_csv(joinpath(OUT, "attenuation.csv"), att)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rate" => RATE_REF, "rows" => rows, "attenuation" => att))

using CairoMakie
fig = Figure(size=(1000, 360))
ax1 = Axis(fig[1, 1], xlabel="h (µm)", ylabel="1 − S(50)", title="Attenuation vs h at ΔR₂(0) = 0.1 /ms")
scatterlines!(ax1, [p["h_um"] for p in ss], [p["attenuation_50"] for p in ss], label="step")
scatterlines!(ax1, [p["h_um"] for p in sl], [p["attenuation_50"] for p in sl], label="linear")
scatterlines!(ax1, [p["h_um"] for p in sx], [p["attenuation_50"] for p in sx], label="exponential")
axislegend(ax1, position=:rb)
ax2 = Axis(fig[1, 2], xlabel="θ_relax", ylabel="ratio linear / step", xscale=log10, title="V4: h* and ρ* ratios")
scatterlines!(ax2, [r.theta for r in rows], [r.h_ratio for r in rows], label="h* ratio")
scatterlines!(ax2, [r.theta for r in rows], [r.rho_ratio for r in rows], label="ρ* ratio")
axislegend(ax2, position=:rt)
save(joinpath(OUT, "p2_4_8_compare.png"), fig)
println("saved to $OUT")
