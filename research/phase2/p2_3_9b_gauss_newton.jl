# P2.3.9 V4, final step — one Gauss–Newton correction of h* after the parabola refinement (p2_3_9_fit_h.jl).
#
# The parabola vertex over h₀·(1 ± 5 %) left residuals of one sign at every θ, i.e. it was not yet the χ² minimum
# at this precision (SEM ≈ 1e-5). Here: Δh = Gauss–Newton step from the residuals of the direct run at the parabola
# h* and the slope dS/dh from the P2.3.4(b) sweep (ln S linear in h between neighbouring grid points); a final direct
# run at h_final = h* + Δh gives χ²/dof and the residual at every readout.
# The jackknife σ is the parabola-stage one (same seeds and spins). Residuals are reported, never pass/fail.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_3_9b_gauss_newton.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "expected.jl"))
include(joinpath(@__DIR__, "fit.jl"))

const SHAPE = get(ENV, "SHAPE", "step")
const SUFFIX = SHAPE == "step" ? "" : "_" * SHAPE
const OUT = phase2_outdir("p2_3_9_fit_h" * SUFFIX)
const SWEEP = joinpath(REPO_ROOT, "research", "results", "phase2", "p2_3_4_h_sweep" * SUFFIX, "sweep.csv")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(sweep_nspins())))

summary = JSON.parsefile(joinpath(OUT, "summary.json"))
(hgrid, S_by_h) = read_sweep(SWEEP)
model_mean = [vec(mean(S, dims=1)) for S in S_by_h]
residual(m, ms, r, rs) = (σ = sqrt.(ms .^ 2 .+ rs .^ 2); [σ[i] > 0 ? (m[i] - r[i]) / σ[i] : 0.0 for i in eachindex(σ)])

rows = NamedTuple[]
for f in summary["fits"]
    θ = f["theta"]
    hp = f["h_star"]
    Sm, Sms = Float64.(f["S_direct"]), Float64.(f["S_direct_sem"])
    Sr, Srs = Float64.(f["S_ref"]), Float64.(f["S_ref_sem"])
    J = sweep_slope(hgrid, model_mean, hp)
    h1 = hp + gauss_newton_step(Sm, Sms, Sr, Srs, J)
    final = run_walls(layer_walls(h=h1, rate=RATE_REF, shape=SHAPE); nspins=NSPINS)
    res = residual(final.mean, final.sem, Sr, Srs)
    χ2 = chi2_of(final.mean, final.sem, Sr, Srs)
    dof = count(sqrt.(final.sem .^ 2 .+ Srs .^ 2) .> 0) - 1
    rel = final.mean ./ Sr .- 1
    f["h_final"] = h1
    f["chi2_final"] = χ2
    f["dof_final"] = dof
    f["residuals_final"] = res
    f["rel_diff_final"] = rel
    f["S_final"] = final.mean
    f["S_final_sem"] = final.sem
    push!(rows, (theta=θ, h_parabola=hp, h_final=h1, h_sigma_jackknife=f["h_star_jackknife_sigma"],
        chi2_per_dof=χ2 / max(dof, 1), max_abs_residual=maximum(abs.(res)), max_abs_rel_diff=maximum(abs.(rel)),
        residual_5ms=res[2], residual_50ms=res[end], S50_ref=Sr[end], S50_model=final.mean[end]))
    @printf("θ = %.4f: h %.5f → %.5f µm, χ²/dof = %.2f, residuals 5 ms %+.1f … 50 ms %+.1f, max |rel diff| = %.1e\n",
        θ, hp, h1, χ2 / max(dof, 1), res[2], res[end], maximum(abs.(rel)))
end
summary["gauss_newton"] = "h_final = parabola h* + one Gauss–Newton step (p2_3_9b_gauss_newton.jl); see h_star_final.csv"
summary["provenance_gauss_newton"] = provenance()
write_csv(joinpath(OUT, "h_star_final.csv"), rows)
write_json(joinpath(OUT, "summary.json"), summary)

using CairoMakie
fig = Figure(size=(1000, 360))
ax1 = Axis(fig[1, 1], xlabel="θ_relax", ylabel="h* (µm)", xscale=log10, title="V4 ($(SHAPE)): h*(θ), ΔR₂(0) = 0.1 /ms, τ = 1e-2")
scatterlines!(ax1, [r.theta for r in rows], [r.h_final for r in rows])
ax2 = Axis(fig[1, 2], xlabel="t (ms)", ylabel="(model − reference) / σ", title="residuals at h* (final)")
for f in summary["fits"]
    lines!(ax2, TIMES_REF, Float64.(f["residuals_final"]), label="θ = $(f["theta"])")
end
hlines!(ax2, [-2.0, 2.0], color=:gray, linestyle=:dash)
Legend(fig[1, 3], ax2, labelsize=8, rowgap=0)
save(joinpath(OUT, "p2_3_9_fit_h.png"), fig)
println("saved to $OUT")
