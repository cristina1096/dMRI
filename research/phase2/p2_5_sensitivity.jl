# P2.5 — sensitivity to the spatial distribution of near-surface relaxation at fixed integrated relaxivity ρ.
#
# For every configuration of sweep SWEEP (h, w, D or rho; see sweep_configs): the layer (shape, h, ρ) and the
# conventional model with θ = ρ√(π/D) (P2.5.0), same walls, D, τ = tau_for(...), seeds 1–10, NSPINS spins, readouts
# 0:10:200 ms. Both models follow identical trajectories, so ΔS(t) = S_layer − S_base is a paired difference per seed.
# Per configuration: max |ΔS| and its time, paired SEM, |ΔS|/SEM, |ΔS| in units of 1/SNR for SNR 100 and 200,
# the best mono-exponential residual of the layer curve, and the baseline's curvature (mono-exponential check).
# Characterisation: nothing here is pass/fail.
#
# Run: SWEEP=h julia --project=research/baseline -t 8 research/phase2/p2_5_sensitivity.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "fit.jl"))
include(joinpath(@__DIR__, "sensitivity.jl"))

const SWEEP = get(ENV, "SWEEP", "h")
const NSPINS = parse(Int, get(ENV, "NSPINS", "50000"))
const TIMES = collect(0.0:10.0:200.0)
const OUT = phase2_outdir("p2_5_sensitivity_" * SWEEP)

base_cache = Dict{Tuple{Float64, Float64, Float64, Float64}, Any}()
function baseline_run(w, D, θ, τ)
    get!(base_cache, (w, D, θ, τ)) do
        run_walls(baseline_walls(theta=θ, w=w); D=D, timestep=τ, times=TIMES, nspins=NSPINS)
    end
end

rows = NamedTuple[]
curves = Dict{String, Any}[]
for c in sweep_configs(SWEEP)
    τ = tau_for(c.rho, c.h, c.shape; w=c.w)
    θ = theta_for(c.rho; D=c.D)
    lay = run_walls(rho_walls(rho=c.rho, h=c.h, shape=c.shape, w=c.w); D=c.D, timestep=τ, times=TIMES, nspins=NSPINS)
    base = baseline_run(c.w, c.D, θ, τ)
    d = paired_difference(lay.S, base.S)
    p = peak(TIMES, d.mean)
    fit = monoexp_residual(TIMES, lay.mean)
    cb = curvature(TIMES, base.mean)
    row = (sweep=SWEEP, shape=c.shape, h_um=c.h, w_um=c.w, D=c.D, rho=c.rho, theta=θ, tau_ms=τ,
        centroid_um=centroid(c.shape, c.h; cutoff=c.w),
        dS_peak=p.value, dS_peak_t_ms=p.t, dS_peak_sem=d.sem[p.index], dS_peak_z=p.value / max(d.sem[p.index], eps()),
        dS_in_snr100_units=abs(p.value) * 100, dS_in_snr200_units=abs(p.value) * 200,
        S_base_end=base.mean[end], S_layer_end=lay.mean[end],
        monoexp_R=fit.R, monoexp_residual=fit.residual, base_curvature=cb, runtime_s=lay.runtime)
    push!(rows, row)
    push!(curves, Dict("config" => Dict(pairs(c)), "tau_ms" => τ, "theta" => θ, "times" => TIMES,
        "S_layer" => lay.mean, "S_layer_sem" => lay.sem, "S_base" => base.mean, "S_base_sem" => base.sem,
        "dS" => d.mean, "dS_sem" => d.sem, "S_layer_per_seed" => [lay.S[i, :] for i in 1:size(lay.S, 1)]))
    @printf("%-11s h=%.3f w=%.1f D=%.0f ρ=%.3f τ=%.0e: max|ΔS| = %+.2e ± %.0e at %3.0f ms (%.0f SEM; ×SNR100 %.2f, ×SNR200 %.2f); mono-exp residual %.1e\n",
        c.shape, c.h, c.w, c.D, c.rho, τ, p.value, d.sem[p.index], p.t, abs(row.dS_peak_z), row.dS_in_snr100_units, row.dS_in_snr200_units, fit.residual)
    flush(stdout)
end
write_csv(joinpath(OUT, "rows.csv"), rows)
write_json(joinpath(OUT, "curves.json"), curves)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "sweep" => SWEEP, "nspins" => NSPINS,
    "times" => TIMES, "rows" => rows,
    "note" => "ΔS = S_layer − S_base on identical trajectories (same seeds); observable at SNR X if |ΔS| ≥ 2/X"))
println("saved to $OUT")
