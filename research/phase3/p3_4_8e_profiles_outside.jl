# P3.4.8 outside for linear and exponential (user request 2026-11-06): square array s = 3.0, R = 1 (gap 1.0), outer layer
# only, ρ = 0.02, D = 3, plus the conventional model at the same τ (paired: identical trajectories, checked).
#   linear: h = 0.1, 0.2, 0.3, 0.5, 0.7, 1.0 (support h; outer layers overlap for h > 0.5);
#   exponential: h = 0.05, 0.1, 0.2, 0.3, 0.5, 1.0 (decay length; support = default cutoff = gap 1.0, so outer supports
#   always reach the neighbouring surfaces and overlap, as the exponential on walls).
# Linear starts at h = 0.1: its surface rate 2ρ/(hκ) at h = 0.05 (0.79 /ms) would force τ < 1e-2 (user: avoid thin layers).
# τ_a = 1e-2 and τ_b = 2.5e-3 (Review C1: fixed-ρ differences checked at two τ). Seeds 1–10, NSPINS (default 20000), readouts
# 0:5:50 — the same seeds, spins and τ as the step run at s = 3.0 (p3_4_8_h_sweep_outside.jl), so all three profiles can be
# compared at the same h and the same centroid (each paired with its own conventional run).
# Monotonicity: neighbour step resolved when |paired mean| > 3 paired SEM at τ_b; sign recorded, not assumed.
# Run: SHAPE=linear|exponential julia --project=research/baseline -t 8 research/phase3/p3_4_8e_profiles_outside.jl

include(joinpath(@__DIR__, "harness.jl"))
const SHAPE = get(ENV, "SHAPE", "linear")
const OUT = phase3_outdir("p3_4_8e_profiles_outside_" * SHAPE)
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const RHO = 0.02
const S = 3.0
const HS = SHAPE == "linear" ? [0.1, 0.2, 0.3, 0.5, 0.7, 1.0] : [0.05, 0.1, 0.2, 0.3, 0.5, 1.0]
const θ = RHO * sqrt(π / D_REF)
const TAUS = [1e-2, 2.5e-3]
centroid(h) = SHAPE == "linear" ? h / 3 : (c = (S - 2R_REF) / h; h * (1 - (c + 1) * exp(-c)) / -expm1(-c))
paired(a, b) = (d = a.S[:, end] .- b.S[:, end]; (mean(d), std(d) / sqrt(length(d))))
models = vcat([("conventional", Cylinders(radius=R_REF, repeats=[S, S], surface_relaxation=θ))],
              [("$(SHAPE)_h$(h)", Cylinders(radius=R_REF, repeats=[S, S], layer_rho_outer=RHO, layer_h_outer=h, layer_shape=SHAPE)) for h in HS])
runs = Dict((n, τ) => run_cylinders(g; region=:outside, sx=S, sy=S, timestep=τ, nspins=NSPINS) for τ in TAUS for (n, g) in models)
identical = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs[("conventional", τ)].finals, runs[(n, τ)].finals))
                for τ in TAUS for (n, _) in models)
@printf("identical trajectories across models (per τ): %s\n", identical)
E50 = RHO * sv_outside(R_REF, S, S) * 50
rows = NamedTuple[]
for τ in TAUS, (k, (n, _)) in enumerate(models)
    r = runs[(n, τ)]
    (dm, ds) = paired(r, runs[("conventional", τ)])
    h = k == 1 ? 0.0 : HS[k - 1]
    push!(rows, (model=n, h=h, centroid=k == 1 ? 0.0 : centroid(h), tau=τ, S50=r.mean[end], sem=r.sem[end], J50=log(r.mean[end]) + E50,
        neglog50=r.neglog_mean[end], neglog_sem=r.neglog_sem[end], neglog_expected=E50, dS_vs_conv=dm, dS_vs_conv_sem=ds))
    @printf("%-18s τ = %.1e: S(50) %.6f ± %.1e, E[−ln M] %.5f ± %.1e (exp. %.5f), ΔS vs conv %+.2e ± %.1e\n",
        n, τ, r.mean[end], r.sem[end], r.neglog_mean[end], r.neglog_sem[end], E50, dm, ds)
end
names = [n for (n, _) in models]
steps = NamedTuple[]
for k in 1:length(names) - 1
    (ma, sa) = paired(runs[(names[k + 1], TAUS[1])], runs[(names[k], TAUS[1])])
    (mb, sb) = paired(runs[(names[k + 1], TAUS[2])], runs[(names[k], TAUS[2])])
    push!(steps, (from=names[k], to=names[k + 1], d_tau_a=ma, sem_tau_a=sa, d_tau_b=mb, sem_tau_b=sb, z_b=mb / sb,
        resolved=abs(mb / sb) > 3, sign=sign(mb), z_tau_shift=(ma - mb) / sqrt(sa^2 + sb^2)))
    @printf("%-18s → %-18s: τ_a %+.2e ± %.1e, τ_b %+.2e ± %.1e (z %+.1f), τ shift z %+.2f\n",
        names[k], names[k + 1], ma, sa, mb, sb, mb / sb, (ma - mb) / sqrt(sa^2 + sb^2))
end
layer_steps = steps[2:end]
signs = unique(s.sign for s in layer_steps if s.resolved)
write_csv(joinpath(OUT, "sweep.csv"), rows)
write_csv(joinpath(OUT, "steps.csv"), steps)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "shape" => SHAPE, "identical_trajectories" => identical,
    "taus" => TAUS, "nspins" => NSPINS, "n_layer_steps" => length(layer_steps), "n_resolved" => count(s -> s.resolved, layer_steps),
    "resolved_signs" => signs, "max_abs_tau_shift_z" => maximum(abs(s.z_tau_shift) for s in steps)))
println("saved to $OUT")
