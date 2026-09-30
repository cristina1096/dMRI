# P1.3.1 — R1: surface relaxation vs. the analytic rate (paper Eq. 17)
#
#     R2_surface = θ_relax · sqrt(D/π) · (S/V)        (fast-diffusion limit)
#
# The simulator applies exp(-θ·sqrt(τ)) to M⊥ at every collision (evolve.jl L555–561).
# The expected number of wall crossings per step for a uniform density is
# sqrt(Dτ/π)·(S/V), so the rate is θ·sqrt(D/π)·(S/V), independent of τ.
# We also report the exact Brownstein–Tarr lowest-mode rate, to show the fast-diffusion
# approximation error is far below the noise.
#
# Variables (all configurations):
#   θ_relax (surface_relaxation) = 0.01 µm⁻¹·ms^-1/2 ... see THETA below
#   D (diffusivity)              = 3 µm²/ms
#   R1 = R2 (global) = 0, off-resonance = 0, permeability = 0, no MT/sticking
#   sequence : GradientEcho, instant 90° excitation (phase 90°), TE = 200 ms (no gradients)
#   readouts : t = 10, 20, …, 100 ms
#   timestep : simulator default (tortuosity-bound, printed and recorded)
#   N_spins  : NSPINS (default 10_000) per seed;  seeds 1…NSEEDS (default 10)
#
# Configurations:
#   slab     : Walls(repeats=w), w = 2 µm → S/V = 2/w  (both faces of every wall)
#   cylinder : single Cylinders(radius=a), a = 1 µm; spins initialised in a box
#              around it and only the inside subset is read → S/V = 2/a
#
# Fitted R2: least-squares slope of -log(S(t)/N) vs t (with intercept).
# Pass: |R2_fit − R2_theory| ≤ 2σ (σ = std of R2_fit across seeds).
#
# Run: julia --project=research/baseline -t 8 research/baseline/r1_surface_relaxation.jl

include(joinpath(@__DIR__, "common.jl"))
using MRIBuilder
import SpecialFunctions: besselj0, besselj1

const OUT = outdir("r1_surface_relaxation")
const NSPINS = parse(Int, get(ENV, "NSPINS", "10000"))
const NSEEDS = parse(Int, get(ENV, "NSEEDS", "10"))
const THETA = 0.01
const D = 3.0
const TIMES = collect(10.0:10.0:100.0)

# --- analytic references -------------------------------------------------------------
eq17(θ, D, SV) = θ * sqrt(D / π) * SV

"Bisection root of f on [lo, hi]."
function bisect(f, lo, hi; tol=1e-14)
    flo = f(lo)
    for _ in 1:200
        mid = (lo + hi) / 2
        fm = f(mid)
        (sign(fm) == sign(flo)) ? (lo = mid; flo = fm) : (hi = mid)
        hi - lo < tol && break
    end
    return (lo + hi) / 2
end

"Brownstein–Tarr lowest mode. slab: ξ tan ξ = ρa/D, a = half-width. cylinder: ξ J1(ξ)/J0(ξ) = ρa/D, a = radius."
function brownstein_tarr(kind, ρ, a, D)
    κ = ρ * a / D
    ξ = kind == :slab ? bisect(x -> x * tan(x) - κ, 1e-12, π / 2 - 1e-9) :
                        bisect(x -> x * besselj1(x) / besselj0(x) - κ, 1e-12, 2.4)
    return D * ξ^2 / a^2
end

function fit_rate(times, ratio)
    y = -log.(ratio)
    X = hcat(ones(length(times)), times)
    coef = X \ y
    return (R2=coef[2], intercept=coef[1])
end

# --- runs ---------------------------------------------------------------------------
function run_config(name, geometry, SV, kind, a; bounding_box, subset)
    seq = GradientEcho(TE=200., TR=1000., scanner=MRIBuilder.Siemens_Prisma, excitation=(phase=90.,))
    sim = Simulation(seq; diffusivity=D, geometry=geometry, verbose=false)
    ts = sim.timestep.max_timestep
    ρ = THETA * sqrt(D / π)
    R2_eq17 = eq17(THETA, D, SV)
    R2_bt = brownstein_tarr(kind, ρ, a, D)
    println("\n== $name: S/V=$(SV) µm⁻¹, τ_max=$(ts) ms, R2_eq17=$(R2_eq17), R2_BT=$(R2_bt)")

    rows = NamedTuple[]
    fits = Float64[]
    signals = Matrix{Float64}(undef, NSEEDS, length(TIMES))
    nsp = Int[]
    runtime = 0.0
    for seed in 1:NSEEDS
        Random.seed!(seed)
        t0 = time()
        res = readout(NSPINS, sim, TIMES; bounding_box=bounding_box, subset=subset)
        runtime += time() - t0
        ratio = [per_spin(r) for r in res]
        signals[seed, :] = ratio
        push!(nsp, res[1].nspins)
        f = fit_rate(TIMES, ratio)
        push!(fits, f.R2)
        @printf("  seed %2d: nspins=%d  R2_fit=%.6f  intercept=%.2e  S(100)=%.5f  (%.1f s)\n", seed, res[1].nspins, f.R2, f.intercept, ratio[end], time() - t0)
        for (t, s) in zip(TIMES, ratio)
            push!(rows, (config=name, seed=seed, t_ms=t, signal=s, nspins=res[1].nspins))
        end
    end
    st = seed_stats(fits)
    dev = st.mean - R2_eq17
    pass = abs(dev) <= 2 * st.sigma
    dev_bt = st.mean - R2_bt
    pass_bt = abs(dev_bt) <= 2 * st.sigma
    @printf("  R2_fit = %.6f ± %.6f (σ over %d seeds, SEM %.2e); theory Eq17 %.6f; deviation %.2e = %.2fσ → %s\n",
        st.mean, st.sigma, st.n, st.sem, R2_eq17, dev, dev / st.sigma, pass ? "PASS" : "FAIL")
    @printf("  vs exact Brownstein–Tarr %.6f: deviation %.2e = %.2fσ → %s  (Eq17/BT − 1 = %.2e)\n",
        R2_bt, dev_bt, dev_bt / st.sigma, pass_bt ? "PASS" : "FAIL", R2_eq17 / R2_bt - 1)
    write_csv(joinpath(OUT, "signal_$(name).csv"), rows)
    return Dict(
        "config" => name, "S_over_V_per_um" => SV, "theta_relax" => THETA, "D" => D,
        "rho_um_per_ms" => ρ, "tau_max_ms" => ts, "binding_constraint" => "tortuosity (0.03·size_scale²/D)",
        "size_scale_um" => mr.Geometries.Internal.size_scale(sim.geometry),
        "nspins_per_seed_simulated" => NSPINS, "nspins_in_subset" => nsp, "nseeds" => NSEEDS,
        "readout_times_ms" => TIMES, "mean_signal" => vec(mean(signals, dims=1)), "sigma_signal" => vec(std(signals, dims=1)),
        "R2_fit_per_seed" => fits, "R2_fit_mean" => st.mean, "R2_fit_sigma" => st.sigma, "R2_fit_sem" => st.sem,
        "R2_eq17" => R2_eq17, "R2_brownstein_tarr" => R2_bt,
        "deviation_in_sigma" => dev / st.sigma, "pass_2sigma" => pass, "runtime_s" => runtime,
        "deviation_vs_BT_in_sigma" => dev_bt / st.sigma, "pass_2sigma_vs_BT" => pass_bt,
        "eq17_over_BT_minus_1" => R2_eq17 / R2_bt - 1,
        "apparent_T2_ms" => 1 / st.mean,
    )
end

results = Any[]
w = 2.0
push!(results, run_config("slab", Walls(repeats=w, surface_relaxation=THETA), 2 / w, :slab, w / 2;
    bounding_box=500, subset=Subset()))
a = 1.0
push!(results, run_config("cylinder", Cylinders(radius=a, surface_relaxation=THETA), 2 / a, :cylinder, a;
    bounding_box=mr.BoundingBox([-a, -a, -a], [a, a, a]), subset=Subset(inside=true)))

write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "results" => results,
    "sequence" => "GradientEcho(TE=200, TR=1000, scanner=Siemens_Prisma, excitation=(phase=90,)) — instant 90° pulse, no gradients"))

# figure
using CairoMakie
fig = Figure(size=(760, 330))
for (i, r) in enumerate(results)
    ax = Axis(fig[1, i], title="$(r["config"]): S/V = $(r["S_over_V_per_um"]) µm⁻¹", xlabel="t (ms)", ylabel="M⊥(t)/M⊥(0)", yscale=log10)
    tt = 0:1:100
    lines!(ax, tt, exp.(-r["R2_eq17"] .* tt), color=:black, label="Eq. 17")
    errorbars!(ax, TIMES, r["mean_signal"], r["sigma_signal"], color=:red)
    scatter!(ax, TIMES, r["mean_signal"], color=:red, label="MCMR (mean ± σ, $(NSEEDS) seeds)")
    axislegend(ax, position=:lb)
end
save(joinpath(OUT, "r1_surface_relaxation.png"), fig)
println("saved to $OUT")
