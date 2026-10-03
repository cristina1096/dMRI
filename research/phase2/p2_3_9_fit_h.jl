# P2.3.9 V4 — fit h*(θ) of the step-profile layer (ΔR2(0) = 0.1 /ms, both faces) to the B4 baseline reference.
#
# (1) Start: h₀ from the P2.3.4(b) sweep, S(t; h) interpolated linearly in h (fit_h). The first pass used h₀ directly,
#     but its residuals were dominated by interpolation error (results kept in first_pass/), so:
# (2) Refine with direct runs: χ²(h) at h₀·(1 − δ), h₀, h₀·(1 + δ) (δ = 5 %), h* = vertex of the parabola through them;
#     final direct run at h* gives χ²/dof and the residual at every readout.
#     Jackknife σ(h*): the same three runs, leaving out one seed at a time.
# (3) Timestep: direct run at h*(θ) with τ = 1e-3 for θ ∈ TAU_CHECK; S(50) shift from τ = 1e-2, in combined SEM.
# (4) Bulk on: θ = 0.01, direct run at h* with R2_bulk = 1/80 vs reference × exp(−t/80) (B4 showed this factorises).
# The baseline is a reference, not the truth: residuals are reported, never pass/fail.
# Model and reference share seeds 1–10, so √(SEM₁² + SEM₂²) overestimates the noise of their difference.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_3_9_fit_h.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "expected.jl"))
include(joinpath(@__DIR__, "fit.jl"))

const OUT = phase2_outdir("p2_3_9_fit_h")
const SWEEP = joinpath(REPO_ROOT, "research", "results", "phase2", "p2_3_4_h_sweep", "sweep.csv")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(sweep_nspins())))
const TAU_CHECK = [0.002, 0.01, 0.05]
const DELTA = 0.05

ref = load_b4()
(hgrid, S_by_h) = read_sweep(SWEEP)
model_mean = [vec(mean(S, dims=1)) for S in S_by_h]
model_sem = [vec(std(S, dims=1)) ./ sqrt(size(S, 1)) for S in S_by_h]
residual(m, ms, r, rs) = (σ = sqrt.(ms .^ 2 .+ rs .^ 2); [σ[i] > 0 ? (m[i] - r[i]) / σ[i] : 0.0 for i in eachindex(σ)])
seed_stats_rows(S, rows) = (vec(mean(S[rows, :], dims=1)), vec(std(S[rows, :], dims=1)) ./ sqrt(length(rows)))

fits = Dict{String, Any}[]
rows = NamedTuple[]
for θ in sort([x for x in keys(ref) if x > 0])
    R = ref[θ]
    h0 = fit_h(hgrid, model_mean, model_sem, R.mean, R.sem).h
    hs = [h0 * (1 - DELTA), h0, h0 * (1 + DELTA)]
    trio = [run_walls(layer_walls(h=h, rate=RATE_REF); nspins=NSPINS) for h in hs]
    χs = [chi2_of(r.mean, r.sem, R.mean, R.sem) for r in trio]
    hstar = parabola_vertex(hs, χs)
    inside = hs[1] <= hstar <= hs[3]
    n = size(trio[1].S, 1)
    loo = map(1:n) do s
        keep = setdiff(1:n, s)
        parabola_vertex(hs, [chi2_of(seed_stats_rows(r.S, keep)..., R.mean, R.sem) for r in trio])
    end
    σjack = sqrt((n - 1) / n * sum((loo .- mean(loo)) .^ 2))
    final = run_walls(layer_walls(h=hstar, rate=RATE_REF); nspins=NSPINS)
    res = residual(final.mean, final.sem, R.mean, R.sem)
    χ2 = chi2_of(final.mean, final.sem, R.mean, R.sem)
    dof = count(sqrt.(final.sem .^ 2 .+ R.sem .^ 2) .> 0) - 1
    entry = Dict{String, Any}("theta" => θ, "h_start" => h0, "h_trio" => hs, "chi2_trio" => χs, "h_star" => hstar,
        "h_star_jackknife_sigma" => σjack, "vertex_inside_trio" => inside, "chi2" => χ2, "dof" => dof,
        "residuals_direct" => res, "S_direct" => final.mean, "S_direct_sem" => final.sem, "S_ref" => R.mean, "S_ref_sem" => R.sem,
        "rel_diff_direct" => final.mean ./ R.mean .- 1)
    if θ in TAU_CHECK
        fine = run_walls(layer_walls(h=hstar, rate=RATE_REF); timestep=1e-3, nspins=NSPINS)
        entry["S50_tau_1e-3"] = fine.mean[end]
        entry["S50_tau_1e-3_sem"] = fine.sem[end]
        entry["tau_shift_in_sem"] = (fine.mean[end] - final.mean[end]) / sqrt(fine.sem[end]^2 + final.sem[end]^2)
        entry["tau_shift_rel"] = fine.mean[end] / final.mean[end] - 1
    end
    if θ == 0.01
        bulk = run_walls(layer_walls(h=hstar, rate=RATE_REF); R2_bulk=1 / 80, nspins=NSPINS)
        f = exp.(-TIMES_REF ./ 80)
        entry["residuals_bulk_on"] = residual(bulk.mean, bulk.sem, R.mean .* f, R.sem .* f)
    end
    push!(fits, entry)
    push!(rows, (theta=θ, h_start=h0, h_star=hstar, h_star_sigma=σjack, vertex_inside=inside, chi2_per_dof=χ2 / max(dof, 1),
        max_abs_residual=maximum(abs.(res)), max_abs_rel_diff=maximum(abs.(final.mean ./ R.mean .- 1)),
        S50_ref=R.mean[end], S50_model=final.mean[end]))
    @printf("θ = %.4f: h₀ = %.4f → h* = %.5f ± %.5f µm%s, χ²/dof = %.2f, max |residual| = %.2f, max |rel diff| = %.1e\n",
        θ, h0, hstar, σjack, inside ? "" : " (vertex outside trio)", χ2 / max(dof, 1), maximum(abs.(res)), maximum(abs.(final.mean ./ R.mean .- 1)))
end
write_csv(joinpath(OUT, "h_star.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rate" => RATE_REF, "tau_ms" => TAU_REF,
    "nspins" => NSPINS, "delta" => DELTA, "h_grid" => hgrid, "fits" => fits,
    "note" => "model and reference share seeds 1–10 (same initial positions): combined SEM overestimates the noise of the difference"))

using CairoMakie
fig = Figure(size=(1000, 360))
ax1 = Axis(fig[1, 1], xlabel="θ_relax", ylabel="h* (µm)", xscale=log10, title="V4: h*(θ), ΔR₂(0) = 0.1 /ms, τ = 1e-2")
θs = [r.theta for r in rows]
errorbars!(ax1, θs, [r.h_star for r in rows], 2 .* [r.h_star_sigma for r in rows]); scatterlines!(ax1, θs, [r.h_star for r in rows])
ax2 = Axis(fig[1, 2], xlabel="t (ms)", ylabel="(model − reference) / σ", title="residuals of direct runs at h*")
for e in fits
    lines!(ax2, TIMES_REF, e["residuals_direct"], label="θ = $(e["theta"])")
end
hlines!(ax2, [-2.0, 2.0], color=:gray, linestyle=:dash)
Legend(fig[1, 3], ax2, labelsize=8, rowgap=0)
save(joinpath(OUT, "p2_3_9_fit_h.png"), fig)
println("saved to $OUT")
