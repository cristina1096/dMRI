# P2.5.6/P2.5.7 — apparent T2 (paper Figs. 9 and 8A analogues) of the conventional model and the three layer profiles at fixed ρ.
#
# SWEEP = te: Fig. 8A analogue, readouts 0:5:200, T2_app at every TE = 5:5:200 for the t2_configs("te") points.
# Sweep SWEEP (w, h, D or rho; see t2_configs) around T2_REF (w = 2, h = 0.2, D = 3, ρ = 1/70 → T2 ≈ 70 ms).
# Both models: R2_bulk = 1/2000, same walls, D, τ = tau_t2(...), seeds 1–10, NSPINS spins, readouts 0:5:100 ms.
# Conventional θ = ρ√(π/D). Per configuration and TE ∈ {20, 50, 100}: T2_app of layer and conventional (per seed),
# paired ΔT2 ± SEM, relative ΔT2/T2, and the two-point measurement noise σ_T2/T2 at SNR 100 and 200.
#
# MATCH = h (default; profiles share ρ and h) or centroid (profiles share ρ and the centroid; sweep "c" instead of "h").
# Run: SWEEP=w [MATCH=centroid] julia --project=research/baseline -t 8 research/phase2/p2_5_t2_apparent.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "sensitivity.jl"))

const SWEEP = get(ENV, "SWEEP", "w")
const NSPINS = parse(Int, get(ENV, "NSPINS", "50000"))
const TIMES = SWEEP == "te" ? collect(0.0:5.0:200.0) : collect(0.0:5.0:100.0)
const TES = SWEEP == "te" ? collect(5.0:5.0:200.0) : [20.0, 50.0, 100.0]     # Fig. 8A: every TE; Fig. 9: three TEs
const R2B = 1 / 2000
const MATCH = Symbol(get(ENV, "MATCH", "h"))                                  # :h (same h) or :centroid (same centroid)
const OUT = phase2_outdir("p2_5_t2_" * SWEEP * (MATCH == :centroid ? "_centroid" : ""))

base_cache = Dict{Tuple{Float64, Float64, Float64, Float64}, Any}()
baseline_run(w, D, θ, τ) = get!(() -> run_walls(baseline_walls(theta=θ, w=w); D=D, R2_bulk=R2B, timestep=τ, times=TIMES, nspins=NSPINS),
    base_cache, (w, D, θ, τ))

rows = NamedTuple[]
for c in t2_configs(SWEEP; match=MATCH)
    τ = tau_t2(c.rho, c.h, c.shape; w=c.w, D=c.D)
    θ = theta_for(c.rho; D=c.D)
    lay = run_walls(rho_walls(rho=c.rho, h=c.h, shape=c.shape, w=c.w); D=c.D, R2_bulk=R2B, timestep=τ, times=TIMES, nspins=NSPINS)
    base = baseline_run(c.w, c.D, θ, τ)
    for TE in TES
        p = paired_t2(lay.S, base.S, TIMES, TE)
        r = base.mean[findfirst(==(TE), TIMES)]
        n100, n200 = t2_noise_rel(r, 100), t2_noise_rel(r, 200)
        level = abs(p.rel) >= 2n100 ? "SNR 100" : abs(p.rel) >= 2n200 ? "SNR 200" : abs(p.dT2) >= 3p.dT2_sem ? "MC only" : "not resolved"
        push!(rows, (sweep=SWEEP, match=String(MATCH), cbar=centroid(c.shape, c.h; cutoff=c.w), shape=c.shape, h_um=c.h, w_um=c.w, D=c.D, rho=c.rho, theta=θ, tau_ms=τ, TE_ms=TE,
            centroid_um=centroid(c.shape, c.h; cutoff=c.w), T2_layer=p.T2a, T2_layer_sem=p.T2a_sem, T2_conv=p.T2b, T2_conv_sem=p.T2b_sem,
            dT2=p.dT2, dT2_sem=p.dT2_sem, rel=p.rel, rel_sem=p.rel_sem, S_conv_TE=r, noise_rel_snr100=n100, noise_rel_snr200=n200, level=level))
        (SWEEP != "te" || TE in (20.0, 50.0, 100.0, 200.0)) && @printf("%-11s h=%.3f w=%4.1f D=%.1f ρ=%.4f TE=%3.0f: T2 layer %.2f, conventional %.2f ms; ΔT2/T2 = %+.2e ± %.0e (2σ_SNR200 = %.1e) → %s\n",
            c.shape, c.h, c.w, c.D, c.rho, TE, p.T2a, p.T2b, p.rel, p.rel_sem, 2n200, level)
    end
    flush(stdout)
end
write_csv(joinpath(OUT, "rows.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "sweep" => SWEEP, "nspins" => NSPINS, "times" => TIMES,
    "TEs" => collect(TES), "R2_bulk" => R2B, "rows" => rows,
    "note" => "T2_app = −TE/ln(S(TE)/S(0)) per seed; ΔT2 paired on identical trajectories; observable at SNR X if |ΔT2/T2| ≥ 2σ_T2/T2"))
println("saved to $OUT")
