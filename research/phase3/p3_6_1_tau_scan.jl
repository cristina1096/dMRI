# P3.6.1: Monte Carlo timestep scan inside a single cylinder against the exact radial continuum.
# SCAN = thin | profiles | radius | rate (table in the plan). For every (configuration, τ): all models share τ, so
# trajectories are identical (checked) and ΔS = S_layer − S_conv is paired per seed. Recorded per model and τ:
#   S(T) MC ± SEM, continuum, z; ΔS(T) MC ± paired SEM, continuum, z; R = −ln S(T)/T MC and continuum;
#   ΔR2(0)τ, √(2Dτ)/h, √(Dτ)/R. T = 50 ms, except the `rate` scan (T = 2R/(2ρ), so E[−ln M] = 2).
# Settings match P3.4.6 / P3.4.6d (seeds 1–10, 2e4 spins in the box [−R, R]³, ρ = 0.02, D = 3) so their rows are reused.
# Run: SCAN=thin julia --project=research/baseline -t 8 research/phase3/p3_6_1_tau_scan.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
const SCAN = get(ENV, "SCAN", "thin")
const OUT = phase3_outdir("p3_6_1_tau_scan_" * SCAN)
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
θ(ρ) = ρ * sqrt(π / D_REF)
κstep(h, R) = 1 - h / (2R)
cyl(R, ρ, shape, h) = shape == "conventional" ? Cylinders(radius=R, surface_relaxation=θ(ρ)) :
                       Cylinders(radius=R, layer_rho_inner=ρ, layer_h_inner=h, layer_shape=shape)
rate0(shape, ρ, h, R) = shape == "conventional" ? 0.0 : ρ / (h * radial_kappa(shape, h, R)) * radial_g(shape, 0.0, radial_support(shape, h, R) / h)

# configurations: (label, R, ρ, T, models::Vector{(name, shape, h)}, taus)
function configs(scan)
    m(shapes) = [(s == "conventional" ? "conventional" : "$(s)_h$(h)", s, h) for (s, h) in shapes]
    if scan == "thin"
        return [("R1", 1.0, 0.02, 50.0, m([("conventional", 0.0), ("step", 0.05), ("step", 0.1), ("step", 0.2)]), [5e-3, 1.25e-3]),
                ("R1_h0.2_tau_c", 1.0, 0.02, 50.0, m([("conventional", 0.0), ("step", 0.2)]), [6.25e-4]),
                ("R1_h0.05_tau_d", 1.0, 0.02, 50.0, m([("conventional", 0.0), ("step", 0.05)]), [3.125e-4])]
    elseif scan == "profiles"
        return [("R1", 1.0, 0.02, 50.0, m([("conventional", 0.0), ("linear", 0.3), ("exponential", 0.1)]), [1.25e-3, 6.25e-4])]
    elseif scan == "radius"
        return [("R$(R)", R, 0.02, 50.0, m([("conventional", 0.0), ("step", 0.1R)]), [1e-2, 2.5e-3, 6.25e-4]) for R in (0.5, 2.0)]
    elseif scan == "rate"
        out = []
        for r0 in (0.5, 2.0)
            ρ = r0 * 0.1 * κstep(0.1, 1.0)
            push!(out, ("rate$(r0)", 1.0, ρ, 2 * 1.0 / (2ρ), m([("step", 0.1)]), [c / r0 for c in (0.05, 0.02, 0.005, 0.00125)]))
        end
        return out
    end
    error("unknown SCAN $scan")
end

rows = NamedTuple[]
allidentical = true
for (label, R, ρ, T, models, taus) in configs(SCAN), τ in taus
    times = collect(0.0:T / 10:T)
    bb = mr.BoundingBox([-R, -R, -R], [R, R, R])
    runs = Dict(n => run_cylinders(cyl(R, ρ, s, h); region=:inside, R=R, bounding_box=bb, timestep=τ, times=times, nspins=NSPINS)
                for (n, s, h) in models)
    if haskey(runs, "conventional")
        ident = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs["conventional"].finals, runs[n].finals)) for (n, _, _) in models)
        global allidentical &= ident
    end
    Sconv = haskey(runs, "conventional") ? radial_signal(R=R, D=D_REF, times=[T], surface=ρ)[1] : NaN
    for (n, s, h) in models
        r = runs[n]
        Sc = s == "conventional" ? Sconv : radial_signal(R=R, D=D_REF, times=[T], layer=(s, ρ, h))[1]
        if s != "conventional" && haskey(runs, "conventional")
            d = r.S[:, end] .- runs["conventional"].S[:, end]
            dm, ds, dc = mean(d), std(d) / sqrt(length(d)), Sc - Sconv
        else
            dm, ds, dc = NaN, NaN, NaN
        end
        push!(rows, (scan=SCAN, config=label, model=n, shape=s, h=h, R=R, rho=ρ, T=T, tau=τ,
            S_mc=r.mean[end], S_sem=r.sem[end], S_cont=Sc, zS=(r.mean[end] - Sc) / r.sem[end],
            dS_mc=dm, dS_sem=ds, dS_cont=dc, zdS=(dm - dc) / ds,
            R_mc=-log(r.mean[end]) / T, R_cont=-log(Sc) / T,
            rate0_tau=rate0(s, ρ, h, R) * τ, step_over_h=s == "conventional" ? NaN : sqrt(2 * D_REF * τ) / h,
            step_over_R=sqrt(D_REF * τ) / R, nspins=NSPINS))
        @printf("%-14s %-16s τ = %.3e: S %.6f ± %.1e vs %.6f (z %+.2f); ΔS %+.2e ± %.1e vs %+.2e (z %+.2f)\n",
            label, n, τ, r.mean[end], r.sem[end], Sc, (r.mean[end] - Sc) / r.sem[end], dm, ds, dc, (dm - dc) / ds)
    end
end
write_csv(joinpath(OUT, "comparison.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "scan" => SCAN, "nspins" => NSPINS,
    "identical_trajectories" => allidentical, "settings" => Dict("D" => D_REF, "seeds" => collect(SEEDS_REF))))
println("saved to $OUT")
