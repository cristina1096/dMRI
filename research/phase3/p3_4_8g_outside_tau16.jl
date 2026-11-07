# P3.4.8 outside, thin-end third timestep for linear and exponential (user decision 2026-11-07, option B): the Group 2
# check inside showed |ΔS| of surface-concentrated profiles still ≈ 10 % low at τ_b = 2.5e-3 (Review C1), and outside
# there is no continuum. For each shape: conventional and the three thinnest layers (linear h = 0.1, 0.2, 0.3;
# exponential h = 0.05, 0.1, 0.2; outer layer only, ρ = 0.02, square array s = 3.0) at τ_c = 6.25e-4 ms, seeds 1–10,
# NSPINS (default 10000) spins per seed (half the main pass; cost ∝ steps), identical trajectories checked.
# Reported per step (conventional → h₁ → h₂ → h₃) and per layer vs conventional: paired ΔS(50) at τ_c, its shift from
# τ_b (p3_4_8e_profiles_outside.jl results; runs at different τ treated as independent), and the linear-in-τ Richardson
# estimate ΔS(τ → 0) ≈ ΔS_c − (ΔS_b − ΔS_c)/3 (τ_b/τ_c = 4).
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_8g_outside_tau16.jl

include(joinpath(@__DIR__, "harness.jl"))
const NSPINS = parse(Int, get(ENV, "NSPINS", "10000"))
const RHO = 0.02
const S = 3.0
const θ = RHO * sqrt(π / D_REF)
const τc = 1e-2 / 16
const τb = 2.5e-3
paired(a, b) = (d = a.S[:, end] .- b.S[:, end]; (mean(d), std(d) / sqrt(length(d))))
readrows(f) = (ls = readlines(f); h = split(ls[1], ','); [Dict(zip(h, split(l, ','))) for l in ls[2:end]])
for (shape, HS) in (("linear", [0.1, 0.2, 0.3]), ("exponential", [0.05, 0.1, 0.2]))
    OUT = phase3_outdir("p3_4_8e_profiles_outside_" * shape)
    main_steps = readrows(joinpath(OUT, "steps.csv"))
    main_sweep = readrows(joinpath(OUT, "sweep.csv"))
    models = vcat([("conventional", Cylinders(radius=R_REF, repeats=[S, S], surface_relaxation=θ))],
                  [("$(shape)_h$(h)", Cylinders(radius=R_REF, repeats=[S, S], layer_rho_outer=RHO, layer_h_outer=h, layer_shape=shape)) for h in HS])
    runs = Dict(n => run_cylinders(g; region=:outside, sx=S, sy=S, timestep=τc, nspins=NSPINS) for (n, g) in models)
    identical = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs["conventional"].finals, runs[n].finals)) for (n, _) in models)
    @printf("%s: identical trajectories: %s\n", shape, identical)
    names = [n for (n, _) in models]
    steps = NamedTuple[]
    for k in 1:length(names) - 1
        (mc, sc) = paired(runs[names[k + 1]], runs[names[k]])
        b = only(r for r in main_steps if r["from"] == names[k] && r["to"] == names[k + 1])
        mb, sb = parse(Float64, b["d_tau_b"]), parse(Float64, b["sem_tau_b"])
        push!(steps, (from=names[k], to=names[k + 1], d_tau_c=mc, sem_tau_c=sc, z_c=mc / sc, resolved=abs(mc / sc) > 3,
            d_tau_b=mb, sem_tau_b=sb, z_shift=(mb - mc) / sqrt(sb^2 + sc^2), richardson=mc - (mb - mc) / 3))
        @printf("%s %-18s → %-18s: τ_c %+.2e ± %.1e (z %+.1f); τ_b %+.2e ± %.1e; shift z %+.2f; Richardson %+.2e\n",
            shape, names[k], names[k + 1], mc, sc, mc / sc, mb, sb, (mb - mc) / sqrt(sb^2 + sc^2), mc - (mb - mc) / 3)
    end
    vsconv = NamedTuple[]
    for n in names[2:end]
        (mc, sc) = paired(runs[n], runs["conventional"])
        b = only(r for r in main_sweep if r["model"] == n && parse(Float64, r["tau"]) == τb)
        mb, sb = parse(Float64, b["dS_vs_conv"]), parse(Float64, b["dS_vs_conv_sem"])
        push!(vsconv, (model=n, dS_tau_c=mc, sem_tau_c=sc, dS_tau_b=mb, sem_tau_b=sb, z_shift=(mb - mc) / sqrt(sb^2 + sc^2),
            richardson=mc - (mb - mc) / 3, rel_shift=(mc - mb) / mb))
        @printf("%s %-18s vs conv: τ_c %+.2e ± %.1e; τ_b %+.2e ± %.1e; shift z %+.2f (%+.0f %%); Richardson %+.2e\n",
            shape, n, mc, sc, mb, sb, (mb - mc) / sqrt(sb^2 + sc^2), 100 * (mc - mb) / mb, mc - (mb - mc) / 3)
    end
    write_csv(joinpath(OUT, "steps_tau16.csv"), steps)
    write_csv(joinpath(OUT, "vs_conv_tau16.csv"), vsconv)
    write_json(joinpath(OUT, "summary_tau16.json"), Dict("provenance" => provenance(), "shape" => shape, "tau_c" => τc,
        "nspins" => NSPINS, "identical_trajectories" => identical, "max_abs_shift_z" => maximum(abs(r.z_shift) for r in vsconv)))
end
println("done")
