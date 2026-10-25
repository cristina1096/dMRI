# P3.4.8 (outside): attenuation vs h at fixed ρ outside cylinders, outer layer only, plus the conventional model.
# Two square arrays (R = 1) so that h covers a large range while h_outer ≤ gap (validation):
#   s = 2.5 (gap 0.5): h = 0.05, 0.1, 0.2 (separated), 0.3, 0.4, 0.5 (outer layers overlap: 2h > gap);
#   s = 3.0 (gap 1.0): h = 0.05, 0.1, 0.2, 0.3, 0.5 (separated), 0.7, 1.0 (overlap).
# No h < 0.05 (user request 2026-10-25: thin layers force τ ∝ h at fixed ρ); with h ≥ 0.05, τ_a = 1e-2 for every run.
# No continuum exists outside (2D), so Monte Carlo with pairing: all h and the conventional model of one array share one
# τ, hence identical trajectories (asserted) and per-seed paired differences between neighbouring h.
# Two passes τ_a = min(1e-2, 0.005/ΔR2(0) of h = 0.05) and τ_b = τ_a/4; the difference between the paired estimates of
# the two passes measures the timestep bias of the h dependence (Review C1).
# Seeds 1–10, NSPINS (default 20000) spins per seed in the 500 µm box, readouts 0:5:50.
# Monotonicity: a neighbour step is resolved when |paired mean| > 3 paired SEM at τ_b; its sign is recorded, not assumed.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_8_h_sweep_outside.jl

include(joinpath(@__DIR__, "harness.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_4_8_h_sweep_outside")
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const RHO = 0.02
const ARRAYS = [(2.5, [0.05, 0.1, 0.2, 0.3, 0.4, 0.5]), (3.0, [0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0])]
const θ = RHO * sqrt(π / D_REF)
const τa = min(1e-2, 0.005 / rate_outside(RHO, 0.05, R_REF))
const TAUS = [τa, τa / 4]
paired(a, b) = (d = a.S[:, end] .- b.S[:, end]; (mean(d), std(d) / sqrt(length(d))))

rows = NamedTuple[]
steps = NamedTuple[]
allidentical = true
fig = Figure(size=(1100, 420))
for (ia, (s, HS)) in enumerate(ARRAYS)
    models = vcat([("conventional", Cylinders(radius=R_REF, repeats=[s, s], surface_relaxation=θ))],
                  [("step_h$(h)", Cylinders(radius=R_REF, repeats=[s, s], layer_rho_outer=RHO, layer_h_outer=h)) for h in HS])
    runs = Dict{Tuple{String, Float64}, Any}()
    for τ in TAUS, (name, g) in models
        runs[(name, τ)] = run_cylinders(g; region=:outside, sx=s, sy=s, timestep=τ, nspins=NSPINS)
        @assert runs[(name, τ)].tau == τ
    end
    identical = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs[("conventional", τ)].finals, runs[(name, τ)].finals))
                    for τ in TAUS for (name, _) in models)
    global allidentical &= identical
    @printf("s = %.1f: identical trajectories across models (per τ): %s\n", s, identical)
    E50 = RHO * sv_outside(R_REF, s, s) * 50
    for τ in TAUS, (name, _) in models
        r = runs[(name, τ)]
        push!(rows, (s=s, model=name, tau=τ, S50=r.mean[end], sem=r.sem[end], J50=log(r.mean[end]) + E50,
            neglog50=r.neglog_mean[end], neglog_sem=r.neglog_sem[end], neglog_expected=E50))
        @printf("s = %.1f %-14s τ = %.2e: S(50) %.6f ± %.1e, J %.3e, E[−ln M] %.5f ± %.1e (expected %.5f)\n", s, name, τ,
            r.mean[end], r.sem[end], log(r.mean[end]) + E50, r.neglog_mean[end], r.neglog_sem[end], E50)
    end
    names = [n for (n, _) in models]
    for k in 1:length(names) - 1
        (ma, sa) = paired(runs[(names[k + 1], TAUS[1])], runs[(names[k], TAUS[1])])
        (mb, sb) = paired(runs[(names[k + 1], TAUS[2])], runs[(names[k], TAUS[2])])
        push!(steps, (s=s, from=names[k], to=names[k + 1], d_tau_a=ma, sem_tau_a=sa, d_tau_b=mb, sem_tau_b=sb, z_b=mb / sb,
            resolved=abs(mb / sb) > 3, sign=sign(mb), z_tau_shift=(ma - mb) / sqrt(sa^2 + sb^2), layer_step=k > 1))
        @printf("s = %.1f %-14s → %-14s: ΔS(50) τ_a %+.2e ± %.1e, τ_b %+.2e ± %.1e (z %+.1f), τ shift z %+.2f\n",
            s, names[k], names[k + 1], ma, sa, mb, sb, mb / sb, (ma - mb) / sqrt(sa^2 + sb^2))
    end
    ax = Axis(fig[1, ia], xlabel="h (µm)", ylabel="S(50) − S_conv(50) (paired)", title="outside, s = $s, R = 1 (gap $(s - 2)), ρ = 0.02")
    for (τ, col) in zip(TAUS, (:gray, :black))
        d = [paired(runs[(n, τ)], runs[("conventional", τ)]) for n in names[2:end]]
        errorbars!(ax, HS, first.(d), last.(d), color=col)
        scatterlines!(ax, HS, first.(d), color=col, label="τ = $(round(τ, sigdigits=3)) ms")
    end
    vlines!(ax, [(s - 2) / 2], color=:red, linestyle=:dash, label="overlap starts (2h = gap)")
    hlines!(ax, [0.0], color=:black, linewidth=0.5)
    axislegend(ax, labelsize=9)
end
layer_steps = filter(st -> st.layer_step, steps)
signs = unique(st.sign for st in layer_steps if st.resolved)
write_csv(joinpath(OUT, "sweep.csv"), rows)
write_csv(joinpath(OUT, "steps.csv"), steps)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "identical_trajectories" => allidentical,
    "taus" => TAUS, "nspins" => NSPINS, "rho" => RHO, "n_layer_steps" => length(layer_steps),
    "n_layer_steps_resolved" => count(st -> st.resolved, layer_steps), "resolved_signs_layer_steps" => signs,
    "monotone" => length(signs) <= 1, "max_abs_tau_shift_z" => maximum(abs(st.z_tau_shift) for st in steps)))
save(joinpath(OUT, "p3_4_8_h_sweep_outside.png"), fig)
println("saved to $OUT")
