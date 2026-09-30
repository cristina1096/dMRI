# P1.3.3 — R3: timestep plateau of the extra-axonal ADC (paper Supp. Fig. S1)
#
# Reproduces mcmr_paper_figures/Figure_S1/turtuosity.ipynb.
# Pass: plateau onset near τ ≈ 0.01 ms (the simulator's default tortuosity bound for r = 1 µm, D = 3).
#
# Variables (copied from the notebook):
#   D = 3 µm²/ms (simulator default); no relaxation, no sequence (Simulation([]))
#   evolve time t = 10 ms;  ADC = var(Δx ∪ Δy) / (2 t)   (x and y displacements pooled; cylinders along z)
#   spins: Snapshot(NSPINS, sim) in the default 1 mm box, keep only extra-axonal spins (inside = false)
#   timestep τ: fixed number, 10 .^ (-3:0.2:1) ms (21 values)
#   configurations:
#     ordered_high : Cylinders(radius=1, repeats=[2.1, 2.1])          gap 0.1 µm, density 0.712
#     ordered_low  : Cylinders(radius=1, repeats=[3.0, 3.0])          gap 1.0 µm, density 0.349
#     random_high  : random_positions_radii([50,50], 0.712, 2; mean=1, variance=0.01), repeats [50,50]
#     random_low   : random_positions_radii([50,50], 0.349, 2; mean=1, variance=0.01), repeats [50,50]
#     (random geometries use Random.seed!(GEOM_SEED) = 1 before generation — the notebook sets no seed)
#   N_spins: NSPINS (default 10_000, as in the notebook);  seeds 1…NSEEDS (default 10)
#
# Plateau analysis: reference = mean ADC over the three smallest τ (pooled seeds).
# Onset τ_2σ = largest τ such that every τ' ≤ τ lies within 2·σ_diff of the reference,
# with σ_diff = sqrt(σ(τ')² + σ_ref²) (σ = std over seeds). Also reported with a 1 % tolerance.
#
# Run: julia --project=research/baseline -t 8 research/baseline/r3_timestep_plateau.jl

include(joinpath(@__DIR__, "common.jl"))

const OUT = outdir("r3_timestep_plateau")
const NSPINS = parse(Int, get(ENV, "NSPINS", "10000"))
const NSEEDS = parse(Int, get(ENV, "NSEEDS", "10"))
const TIMESTEPS = 10 .^ collect(-3:0.2:1)
const T_EVOLVE = 10.0
const GEOM_SEED = 1

get_density(radius, dist) = π * radius^2 / (radius * 2 + dist)^2

function make_geometry(name)
    if name == "ordered_high"
        return Cylinders(radius=1., repeats=[2.1, 2.1])
    elseif name == "ordered_low"
        return Cylinders(radius=1., repeats=[3.0, 3.0])
    else
        density = name == "random_high" ? get_density(1, 0.1) : get_density(1, 1.)
        Random.seed!(GEOM_SEED)
        pos, _ = random_positions_radii([50., 50.], density, 2; mean=1., variance=0.1^2)
        # the notebook passes radius=1 (not r) for the random configuration; kept identical here
        return Cylinders(radius=1., position=pos, repeats=[50., 50.])
    end
end

function adc(sim, seed)
    Random.seed!(seed)
    start = mr.get_subset(Snapshot(NSPINS, sim), sim; inside=false)
    final = evolve(start, sim, T_EVOLVE)
    d = [f.position .- s.position for (s, f) in zip(start, final)]
    flat = vcat([x[1] for x in d], [x[2] for x in d])
    return var(flat) / (2 * T_EVOLVE), length(start)
end

configs = ["ordered_high", "ordered_low", "random_high", "random_low"]
results = Dict{String, Any}()
rows = NamedTuple[]
for name in configs
    geom = make_geometry(name)
    default_sim = Simulation([]; geometry=geom, verbose=false)
    default_tau = default_sim.timestep.max_timestep
    println("\n== $name: default τ_max = $default_tau ms (size_scale = $(mr.Geometries.Internal.size_scale(default_sim.geometry)) µm)")
    means = Float64[]; sigmas = Float64[]; runt = Float64[]; nextra = Float64[]
    for τ in TIMESTEPS
        sim = Simulation([]; geometry=geom, timestep=τ, verbose=false)
        vals = Float64[]; ns = Int[]
        t0 = time()
        for seed in 1:NSEEDS
            (v, n) = adc(sim, seed)
            push!(vals, v); push!(ns, n)
            push!(rows, (config=name, tau_ms=τ, seed=seed, ADC=v, n_extra=n))
        end
        push!(runt, time() - t0)
        st = seed_stats(vals)
        push!(means, st.mean); push!(sigmas, st.sigma); push!(nextra, mean(ns))
        @printf("  τ=%.3e ms  ADC=%.4f ± %.4f  (n_extra≈%d, %.1f s)\n", τ, st.mean, st.sigma, round(Int, mean(ns)), runt[end])
    end
    ref = mean(means[1:3]); σref = sqrt(mean(sigmas[1:3] .^ 2) / 3)
    # onset: largest τ such that all smaller τ' are within tolerance
    ok2σ = [abs(means[i] - ref) <= 2 * sqrt(sigmas[i]^2 + σref^2) for i in eachindex(TIMESTEPS)]
    ok1 = [abs(means[i] - ref) / ref <= 0.01 for i in eachindex(TIMESTEPS)]
    lastok(v) = (k = findfirst(!, v); isnothing(k) ? TIMESTEPS[end] : TIMESTEPS[max(k - 1, 1)])
    onset2σ = lastok(ok2σ); onset1pct = lastok(ok1)
    idx_default = argmin(abs.(log10.(TIMESTEPS) .- log10(default_tau)))
    @printf("  plateau ADC = %.4f;  onset (2σ) = %.3e ms;  onset (1%%) = %.3e ms;  ADC at τ≈default (%.3e) = %.4f\n",
        ref, onset2σ, onset1pct, TIMESTEPS[idx_default], means[idx_default])
    results[name] = Dict("default_tau_max_ms" => default_tau, "timesteps_ms" => TIMESTEPS,
        "ADC_mean" => means, "ADC_sigma" => sigmas, "n_extra_mean" => nextra, "runtime_s" => runt,
        "plateau_ADC" => ref, "plateau_sigma" => σref, "onset_2sigma_ms" => onset2σ, "onset_1pct_ms" => onset1pct,
        "ADC_at_default_tau" => means[idx_default])
end
write_csv(joinpath(OUT, "adc.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "results" => results,
    "params" => Dict("D" => 3.0, "t_evolve_ms" => T_EVOLVE, "nspins_initial" => NSPINS, "nseeds" => NSEEDS,
        "geometry_seed_random" => GEOM_SEED, "ADC_definition" => "var of pooled x,y displacement / (2 t)")))

using CairoMakie
fig = Figure(size=(620, 400))
ax = Axis(fig[1, 1], xscale=log10, xlabel="timestep τ (ms)", ylabel="ADC after 10 ms (µm²/ms)", title="Supp. Fig. S1 reproduction")
for name in configs
    r = results[name]
    lines!(ax, TIMESTEPS, r["ADC_mean"], label=name, linewidth=3)
    band!(ax, TIMESTEPS, r["ADC_mean"] .- r["ADC_sigma"], r["ADC_mean"] .+ r["ADC_sigma"], alpha=0.3)
end
vlines!(ax, [0.01], color=:black, linestyle=:dash)
ylims!(ax, 0, 3)
axislegend(ax, position=:lb)
save(joinpath(OUT, "r3_timestep_plateau.png"), fig)
println("saved to $OUT")
