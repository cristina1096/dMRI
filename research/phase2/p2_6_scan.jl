# P2.6 τ scan on the 1D toy (exact layer rule, step profile, layer on both faces of walls every 2 µm).
#
# Environment: H (µm), RATE = ΔR2(0) (1/ms), D (µm²/ms), TAU_LO, TAU_HI (ms, 3 per decade), NSPINS, T (ms), LABEL.
# Observable: R = −ln S(T_eff)/T_eff per seed (seeds 1–10). Reference: the smallest τ.
# Reports τ_conv at relative tolerance 1 % and 0.1 % (interpolated in log τ), whether it is only a lower bound,
# whether the reference is too noisy for that tolerance, and the sign of the bias at the largest τ.
#
# Run (example, V1): H=0.1 RATE=0.1 D=3 TAU_LO=1e-5 TAU_HI=1e-1 NSPINS=100000 T=10 LABEL=V1_h0.1 \
#      julia --project=research/baseline -t 8 research/phase2/p2_6_scan.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "toy.jl"))
include(joinpath(@__DIR__, "timestep.jl"))

envf(k, d) = parse(Float64, get(ENV, k, string(d)))
const H, RATE, DIFF = envf("H", 0.1), envf("RATE", 0.1), envf("D", 3.0)
const TAUS = log_taus(envf("TAU_LO", 1e-5), envf("TAU_HI", 1e-1))
const NSPINS = parse(Int, get(ENV, "NSPINS", "100000"))
const T = envf("T", 10.0)
const LABEL = get(ENV, "LABEL", "h$(H)_rate$(RATE)_D$(DIFF)")
const OUT = phase2_outdir("p2_6_scan_" * LABEL)

rows = NamedTuple[]
for τ in TAUS
    t0 = time()
    r = toy_rate(; h=H, rate=RATE, D=DIFF, tau=τ, T=T, nspins=NSPINS, seeds=1:10)
    rt = time() - t0
    push!(rows, (tau_ms=τ, R=r.mean, R_sem=r.sem, T_eff=r.T_eff, runtime_s=rt,
        ns_per_spin_step=rt / (10 * NSPINS * round(T / τ)) * 1e9))
    @printf("τ = %.2e: R = %.6e ± %.1e (%.0f s)\n", τ, r.mean, r.sem, rt)
end
iref = argmin(TAUS)
R_ref, R_ref_sem = rows[iref].R, rows[iref].R_sem
Rs = [r.R for r in rows]
conv = Dict{String, Any}()
for tol in (0.01, 0.001)
    ci = tau_conv_interp(TAUS, Rs, R_ref; tol=tol)
    cg = tau_conv_grid(TAUS, Rs, R_ref; tol=tol)
    conv[string(tol)] = Dict("tau_conv" => ci.tau, "tau_conv_grid" => cg.tau, "lower_bound" => ci.lower_bound,
        "noise_limited" => noise_limited(R_ref, R_ref_sem; tol=tol))
    @printf("tol %.1f %%: τ_conv = %.3e ms%s%s\n", 100tol, ci.tau, ci.lower_bound ? " (lower bound)" : "",
        noise_limited(R_ref, R_ref_sem; tol=tol) ? " (noise-limited)" : "")
end
bias_large = (rows[end].R - R_ref) / R_ref
@printf("bias at τ = %.1e: %+.2e relative (%s attenuation than the reference)\n", TAUS[end], bias_large, bias_large < 0 ? "less" : "more")
write_csv(joinpath(OUT, "scan.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "h_um" => H, "rate" => RATE, "D" => DIFF,
    "T_ms" => T, "nspins" => NSPINS, "nseeds" => 10, "taus" => TAUS, "R_ref" => R_ref, "R_ref_sem" => R_ref_sem,
    "tau_conv" => conv, "bias_at_largest_tau" => bias_large,
    "h2_over_D" => H^2 / DIFF, "step_length_at_tau_conv_1pct" => sqrt(2DIFF * conv["0.01"]["tau_conv"])))

using CairoMakie
fig = Figure(size=(560, 380))
ax = Axis(fig[1, 1], xscale=log10, xlabel="τ (ms)", ylabel="R/R_ref − 1", title="toy τ scan: h = $H µm, ΔR₂(0) = $RATE, D = $DIFF")
errorbars!(ax, TAUS, Rs ./ R_ref .- 1, [r.R_sem for r in rows] ./ R_ref); scatterlines!(ax, TAUS, Rs ./ R_ref .- 1)
hlines!(ax, [-0.01, 0.01], color=:gray, linestyle=:dash); hlines!(ax, [-0.001, 0.001], color=:gray, linestyle=:dot)
save(joinpath(OUT, "scan.png"), fig)
println("saved to $OUT")
