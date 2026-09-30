# P1.3.4 — R4: timestep-correction schemes for surface relaxation, 1D slab (paper Supp. Fig. S2C)
#
# Part A reproduces the stand-alone toy model in mcmr_paper_figures/Figure_S2/toy_model.ipynb
# (cell 6: "Transverse signal attenuation" for three correction schemes).
# Part B runs the same physical set-up through the unmodified MCMRSimulator, which implements
# the logarithmic scheme (evolve.jl L555–561: M⊥ *= exp(-θ·sqrt(τ)) per collision).
#
# Variables:
#   geometry  : walls every 1 µm (toy: compartments [k, k+1); MCMR: Walls(repeats=1))
#   D         : 0.5 µm²/ms (toy steps are randn·sqrt(τ) → variance τ = 2Dτ)
#   durations : T = 10, 1, 0.1 ms;   θ_relax = 1/T  (µm⁻¹·ms^-1/2 units of the toy; "θ used" for R4)
#   N_steps   : [1, 2, 3, 5, 8, 10, 20, 30, 50, 80, 100, 300, 1000, 3000];  τ = T / N_steps
#   schemes   : linear      survival per hit = 1 − θ√τ
#               logarithmic survival per hit = exp(−θ√τ)                   ← MCMRSimulator
#               log+return  survival per hit = exp(−θ√τ)·I0(θ√τ)
#   observable: attenuation = 1 − ⟨survival^n_hits⟩ = 1 − S(T)/S(0)
#   spins     : toy NSPINS_TOY (default 1e5, as notebook); MCMR NSPINS_MCMR (default 1e4)
#   seeds     : 1…NSEEDS (default 10)
#
# Pass (tracker): logarithmic scheme flat to τ ≈ 1e-3 ms. Quantified per duration as the largest τ for
# which the log-scheme attenuation stays within 2σ (and within 5 %) of its smallest-τ value.
#
# Run: julia --project=research/baseline -t 8 research/baseline/r4_correction_schemes.jl

include(joinpath(@__DIR__, "common.jl"))
using MRIBuilder
import Bessels: besseli0

const OUT = outdir("r4_correction_schemes")
const NSEEDS = parse(Int, get(ENV, "NSEEDS", "10"))
const NSPINS_TOY = parse(Int, get(ENV, "NSPINS_TOY", "100000"))
const NSPINS_MCMR = parse(Int, get(ENV, "NSPINS_MCMR", "10000"))
const NSTEPS = [1, 2, 3, 5, 8, 10, 20, 30, 50, 80, 100, 300, 1000, 3000]
const DURATIONS = [10.0, 1.0, 0.1]
const SCHEMES = [:linear, :exp, :bessel]
const D = 0.5

# ---- Part A: toy model (verbatim logic from the notebook, collisions counted per spin) ----
function survival(rate, timestep, model)
    model == :linear && return 1 - sqrt(timestep) * rate
    model == :exp && return exp(-sqrt(timestep) * rate)
    model == :bessel && return exp(-sqrt(timestep) * rate) * besseli0(sqrt(timestep) * rate)
    error("unknown model")
end

"Number of wall collisions per spin for a 1D walk with walls at every integer (all reflective)."
function toy_collisions(Nsteps, Nspins, duration)
    timestep = duration / Nsteps
    pos = rand(Nspins) .* 1000
    comp = floor.(Int, pos)
    ncoll = zeros(Int, Nspins)
    Threads.@threads for i in 1:Nspins
        p = pos[i]; c = comp[i]; n = 0
        for _ in 1:Nsteps
            step = randn() * sqrt(timestep)
            while floor(Int, p + step) != c
                pos_step = step > 0
                n += 1
                wall = pos_step ? c + 1 : c
                step = abs(step) - abs(p - wall)
                p = wall
                step = pos_step ? -step : step
            end
            p += step
        end
        ncoll[i] = n
    end
    return ncoll
end

toy_rows = NamedTuple[]
toy = Dict{String, Any}()
for T in DURATIONS
    θ = 1 / T
    att = Dict(m => zeros(NSEEDS, length(NSTEPS)) for m in SCHEMES)
    for seed in 1:NSEEDS, (j, ns) in enumerate(NSTEPS)
        Random.seed!(1000 * seed + j)
        nc = toy_collisions(ns, NSPINS_TOY, T)
        τ = T / ns
        for m in SCHEMES
            a = mean(@. 1 - survival(θ, τ, m)^nc)
            att[m][seed, j] = a
            push!(toy_rows, (duration_ms=T, theta=θ, scheme=string(m), seed=seed, nsteps=ns, tau_ms=τ, attenuation=a))
        end
    end
    toy[string(T)] = Dict(string(m) => Dict("mean" => vec(mean(att[m], dims=1)), "sigma" => vec(std(att[m], dims=1))) for m in SCHEMES)
    @printf("toy T=%5.2f ms θ=%.2f: attenuation at τ_min (%.2e ms): linear %.4f  log %.4f  log+return %.4f\n", T, θ, T / NSTEPS[end],
        toy[string(T)]["linear"]["mean"][end], toy[string(T)]["exp"]["mean"][end], toy[string(T)]["bessel"]["mean"][end])
end
write_csv(joinpath(OUT, "toy_attenuation.csv"), toy_rows)

# ---- Part B: same set-up through MCMRSimulator (logarithmic scheme) --------------------
mcmr_rows = NamedTuple[]
mcmr = Dict{String, Any}()
for T in DURATIONS
    θ = 1 / T
    seq = GradientEcho(TE=T, TR=10 * T, scanner=MRIBuilder.Siemens_Prisma, excitation=(phase=90.,))
    att = zeros(NSEEDS, length(NSTEPS))
    default_tau = Simulation(seq; diffusivity=D, geometry=Walls(repeats=1., surface_relaxation=θ), verbose=false).timestep.max_timestep
    for (j, ns) in enumerate(NSTEPS)
        τ = T / ns
        sim = Simulation(seq; diffusivity=D, geometry=Walls(repeats=1., surface_relaxation=θ), timestep=τ, verbose=false)
        for seed in 1:NSEEDS
            Random.seed!(seed)
            s = readout(NSPINS_MCMR, sim, T; bounding_box=500)
            att[seed, j] = 1 - per_spin(s)
            push!(mcmr_rows, (duration_ms=T, theta=θ, seed=seed, nsteps=ns, tau_ms=τ, attenuation=att[seed, j]))
        end
    end
    mcmr[string(T)] = Dict("mean" => vec(mean(att, dims=1)), "sigma" => vec(std(att, dims=1)), "default_tau_max_ms" => default_tau)
    @printf("MCMR T=%5.2f ms θ=%.2f: attenuation %s  (default τ_max = %.3e ms)\n", T, θ,
        join([@sprintf("%.4f", a) for a in mcmr[string(T)]["mean"]], " "), default_tau)
end
write_csv(joinpath(OUT, "mcmr_attenuation.csv"), mcmr_rows)

# ---- flatness analysis ----------------------------------------------------------------
function flat_limit(mean_v, sig_v, taus; rel=0.05)
    # reference: smallest τ (last element). Walk from small τ to large τ.
    ref = mean_v[end]; σr = sig_v[end]
    lim2σ = taus[end]; lim_rel = taus[end]
    for k in length(taus)-1:-1:1
        within2σ = abs(mean_v[k] - ref) <= 2 * sqrt(sig_v[k]^2 + σr^2)
        within2σ && lim2σ == taus[k+1] && (lim2σ = taus[k])
        withinrel = abs(mean_v[k] - ref) / ref <= rel
        withinrel && lim_rel == taus[k+1] && (lim_rel = taus[k])
    end
    return (ref=ref, flat_2sigma_up_to_ms=lim2σ, flat_5pct_up_to_ms=lim_rel)
end

flat = Dict{String, Any}()
for T in DURATIONS
    taus = T ./ NSTEPS
    ft = Dict(string(m) => flat_limit(toy[string(T)][string(m)]["mean"], toy[string(T)][string(m)]["sigma"], taus) for m in SCHEMES)
    fm = flat_limit(mcmr[string(T)]["mean"], mcmr[string(T)]["sigma"], taus)
    flat[string(T)] = Dict("toy" => ft, "mcmr_log" => fm)
    @printf("T=%5.2f: flat(5%%) up to τ = linear %.2e | log %.2e | log+return %.2e | MCMR(log) %.2e ms;  flat(2σ): log %.2e, MCMR %.2e\n",
        T, ft["linear"].flat_5pct_up_to_ms, ft["exp"].flat_5pct_up_to_ms, ft["bessel"].flat_5pct_up_to_ms, fm.flat_5pct_up_to_ms,
        ft["exp"].flat_2sigma_up_to_ms, fm.flat_2sigma_up_to_ms)
end

write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "toy" => toy, "mcmr" => mcmr, "flatness" => flat,
    "params" => Dict("D" => D, "wall_spacing_um" => 1.0, "durations_ms" => DURATIONS, "theta_relax" => "1/duration",
        "nsteps" => NSTEPS, "nspins_toy" => NSPINS_TOY, "nspins_mcmr" => NSPINS_MCMR, "nseeds" => NSEEDS)))

using CairoMakie
fig = Figure(size=(1000, 330))
for (i, T) in enumerate(DURATIONS)
    taus = T ./ NSTEPS
    ax = Axis(fig[1, i], xscale=log10, xlabel="timestep τ (ms)", ylabel=i == 1 ? "transverse attenuation" : "",
        title="T = $T ms, θ = $(1/T)")
    for (m, lab) in zip(SCHEMES, ["linear 1−θ√τ", "log e^{−θ√τ}", "log+return"])
        r = toy[string(T)][string(m)]
        lines!(ax, taus, r["mean"], label="toy: $lab", linewidth=3)
    end
    r = mcmr[string(T)]
    scatter!(ax, taus, r["mean"], color=:black, marker=:x, markersize=12, label="MCMRSimulator")
    vlines!(ax, [0.01 * T], color=:black, linestyle=:dash)
    i == 1 && axislegend(ax, position=:lb, labelsize=10)
end
save(joinpath(OUT, "r4_correction_schemes.png"), fig)
println("saved to $OUT")
