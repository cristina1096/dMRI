# B4 — baseline surface-relaxation reference curves S_θ(t) on walls (reference for fitting the R₂(d) model)
#
# The proposed model R2_total = R2_bulk + ΔR2(d) is fitted to the baseline θ_relax model: for each θ,
# Phase 2 finds the layer thickness h (fixed profile g, fixed surface excess rate ΔR2(0)) whose signal
# best matches S_θ(t). The baseline is a REFERENCE, not the true answer, so this run records the
# reference curves with their Monte Carlo σ and no analytic formula. It also checks that the reference
# is timestep-converged and that bulk R2 multiplies it exactly.
#
# Variables:
#   geometry : Walls(repeats=w, surface_relaxation=θ), w = 2 µm, normal along x; θ acts on BOTH faces of
#              every gap (the R₂(d) layer will also be put on both faces)
#   D = 3 µm²/ms; R1 = 0; R2_bulk ∈ {0, 1/80}; permeability / off-resonance / MT off
#   sequence : none (empty sequence, no RF); spins start transverse = 1, longitudinal = 0, phase = 0
#   readouts : t = 0, 5, …, 50 ms
#   θ        : THETAS (0 plus 12 values 0.001–0.05; pilot: S(50) ≈ 0.95 → 0.09)
#   τ        : default (0.04 ms, tortuosity) and 1e-2 ms for every θ (NSEEDS seeds);
#              1e-3 ms for THETAS_FINE (NSEEDS_FINE seeds) as a convergence check
#   N_spins  : NSPINS (default 100_000) per seed; seed k gives the same initial positions in every run
#
# Reported per (θ, τ, t): mean S, σ over seeds, SEM. Descriptive only: least-squares slope of −log S vs t.
# Checks:
#   timestep : |S_default − S_1e-2| and |S_1e-3 − S_1e-2| in units of the combined SEM, at every t
#   bulk     : S(θ, R2_bulk = 1/80) / S(θ, 0) = exp(−t/80), same seeds, to the roundoff bound (n_steps + N)·eps
# Reference curve = τ = 1e-2 ms (all θ, NSEEDS seeds), provided the 1e-3 check agrees within 2 SEM.
#
# Run: julia --project=research/baseline -t 8 research/baseline/b4_theta_reference.jl

include(joinpath(@__DIR__, "common.jl"))

const OUT = outdir("b4_theta_reference")
const NSPINS = parse(Int, get(ENV, "NSPINS", "100000"))
const NSEEDS = parse(Int, get(ENV, "NSEEDS", "10"))
const NSEEDS_FINE = parse(Int, get(ENV, "NSEEDS_FINE", "5"))
const D = 3.0
const W = 2.0
const R2_BULK = 1 / 80
const TIMES = collect(0.0:5.0:50.0)
const THETAS = [0.0, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05]
const THETAS_FINE = [0.001, 0.005, 0.02, 0.05]

make_sim(θ, τ, R2) = Simulation([mr.SequenceParts.empty_sequence()]; geometry=Walls(repeats=W, surface_relaxation=θ),
    diffusivity=D, R2=R2, verbose=false, (isnothing(τ) ? (;) : (; timestep=τ))...)

"Signal matrix (seed × time) and runtime for one configuration."
function run(θ, τ, R2, nseeds)
    sim = make_sim(θ, τ, R2)
    S = Matrix{Float64}(undef, nseeds, length(TIMES))
    rt = 0.0
    for seed in 1:nseeds
        Random.seed!(seed)
        snap = Snapshot(NSPINS, sim, 500; transverse=1.0, longitudinal=0.0)
        t0 = time()
        res = readout(snap, sim, TIMES)
        rt += time() - t0
        S[seed, :] = [per_spin(r) for r in res]
    end
    return (S=S, tau=sim.timestep.max_timestep, runtime=rt)
end

function fit_rate(S)    # decay rate of a straight-line fit to -ln S(t)
    X = hcat(ones(length(TIMES)), TIMES)
    return (X \ (-log.(S)))[2]
end

colmean(S) = vec(mean(S, dims=1))
colsem(S) = vec(std(S, dims=1)) ./ sqrt(size(S, 1))

rows = NamedTuple[]
checks_tau = NamedTuple[]
checks_bulk = NamedTuple[]
curves = Dict{Tuple{Float64, String}, Any}()
bulk_tol = (TIMES[end] / 0.04 + NSPINS) * eps()

open(joinpath(OUT, "run.log"), "w") do io
    for θ in THETAS
        tee(io, @sprintf("\n== θ = %g", θ))
        runs = Dict{String, Any}()
        runs["default"] = run(θ, nothing, 0.0, NSEEDS)
        runs["tau_1e-2"] = run(θ, 1e-2, 0.0, NSEEDS)
        θ in THETAS_FINE && (runs["tau_1e-3"] = run(θ, 1e-3, 0.0, NSEEDS_FINE))
        for (label, r) in runs
            m, s = colmean(r.S), vec(std(r.S, dims=1))
            curves[(θ, label)] = (mean=m, sigma=s, sem=colsem(r.S), tau=r.tau, nseeds=size(r.S, 1))
            for (it, t) in enumerate(TIMES)
                push!(rows, (theta=θ, config=label, tau_ms=r.tau, R2_bulk=0.0, t_ms=t, mean_signal=m[it], sigma=s[it],
                    sem=s[it] / sqrt(size(r.S, 1)), nseeds=size(r.S, 1), nspins=NSPINS))
            end
            rate = θ == 0 ? 0.0 : fit_rate(m)
            tee(io, @sprintf("  %-9s τ=%-7g S(25)=%.5f±%.5f  S(50)=%.5f±%.5f (σ)  descriptive rate %.6f ms⁻¹  runtime %.0f s",
                label, r.tau, m[6], s[6], m[end], s[end], rate, r.runtime))
        end
        # timestep checks against τ = 1e-2
        ref = curves[(θ, "tau_1e-2")]
        for label in ("default", "tau_1e-3")
            haskey(curves, (θ, label)) || continue
            c = curves[(θ, label)]
            z = (c.mean .- ref.mean) ./ sqrt.(c.sem .^ 2 .+ ref.sem .^ 2)
            z[1] = 0.0  # t = 0 identical by construction
            θ == 0 && (z .= 0.0)  # θ = 0: S ≡ 1 in both, z = 0/0
            zmax = θ == 0 ? 0.0 : maximum(abs.(z[2:end]))
            push!(checks_tau, (theta=θ, compare=label * " vs tau_1e-2", max_abs_z=zmax, z_at_50=z[end],
                rel_diff_at_50=c.mean[end] / ref.mean[end] - 1))
            tee(io, @sprintf("  τ check %-8s vs 1e-2: max |z| = %.2f, S(50) rel diff = %+.2e", label, zmax, c.mean[end] / ref.mean[end] - 1))
        end
        # bulk factorisation at default τ (same seeds)
        rb = run(θ, nothing, R2_BULK, NSEEDS)
        dev = maximum(abs.(rb.S ./ (runs["default"].S .* exp.(-R2_BULK .* TIMES')) .- 1))
        push!(checks_bulk, (theta=θ, R2_bulk=R2_BULK, max_abs_rel_dev=dev, tol=bulk_tol, pass=dev <= bulk_tol))
        tee(io, @sprintf("  bulk: max |S(θ,1/80)/(S(θ,0)·e^(−t/80)) − 1| = %.2e (tol %.1e) → %s", dev, bulk_tol, dev <= bulk_tol ? "PASS" : "FAIL"))
    end
end

write_csv(joinpath(OUT, "reference_curves.csv"), rows)
write_csv(joinpath(OUT, "timestep_check.csv"), checks_tau)
write_csv(joinpath(OUT, "bulk_check.csv"), checks_bulk)
tau_ok = all(c -> c.max_abs_z <= 2, filter(c -> startswith(c.compare, "tau_1e-3"), checks_tau))
write_json(joinpath(OUT, "summary.json"), Dict(
    "provenance" => provenance(),
    "geometry" => "Walls(repeats=$W, surface_relaxation=θ): planes x = k·$W, θ on both faces",
    "physics" => "D=$D, R1=0, R2_bulk ∈ {0, $R2_BULK}, permeability/off-resonance/MT off",
    "sequence" => "empty_sequence() — no RF; spins start transverse=1, longitudinal=0, phase=0",
    "thetas" => THETAS, "thetas_fine" => THETAS_FINE, "readout_times_ms" => TIMES,
    "nspins_per_seed" => NSPINS, "nseeds" => NSEEDS, "nseeds_fine" => NSEEDS_FINE,
    "reference_config" => "tau_1e-2", "tau_1e-3_agrees_within_2sem" => tau_ok,
    "bulk_all_pass" => all(c -> c.pass, checks_bulk), "bulk_tol" => bulk_tol,
    "timestep_check" => checks_tau, "bulk_check" => checks_bulk,
))

# figure
using CairoMakie
fig = Figure(size=(1150, 380))
ax1 = Axis(fig[1, 1], title="Reference curves S_θ(t) (τ = 1e-2 ms, mean ± σ)", xlabel="t (ms)", ylabel="M⊥(t)/M⊥(0)", yscale=log10)
cmap = cgrad(:viridis, length(THETAS), categorical=true)
for (i, θ) in enumerate(THETAS)
    c = curves[(θ, "tau_1e-2")]
    lines!(ax1, TIMES, c.mean, color=cmap[i], label="θ = $θ")
    errorbars!(ax1, TIMES, c.mean, c.sigma, color=cmap[i])
end
Legend(fig[1, 2], ax1, labelsize=9, rowgap=0)
ax2 = Axis(fig[1, 3], title="S(50 ms) vs θ", xlabel="θ_relax", ylabel="S(50)", xscale=log10, yscale=log10)
th = THETAS[2:end]
scatterlines!(ax2, th, [curves[(θ, "tau_1e-2")].mean[end] for θ in th], label="τ = 1e-2")
scatter!(ax2, th, [curves[(θ, "default")].mean[end] for θ in th], marker=:x, label="default τ = 0.04")
scatter!(ax2, THETAS_FINE, [curves[(θ, "tau_1e-3")].mean[end] for θ in THETAS_FINE], marker=:utriangle, label="τ = 1e-3")
axislegend(ax2, position=:lb, labelsize=9)
ax3 = Axis(fig[1, 4], title="Timestep check vs τ = 1e-2 (max |z| over t)", xlabel="θ_relax", ylabel="max |z|", xscale=log10, limits=(nothing, (0, 5)))
band!(ax3, [minimum(th), maximum(th)], [0, 0], [2, 2], color=(:gray, 0.25))
for (label, mk) in (("default", :x), ("tau_1e-3", :utriangle))
    cs = filter(c -> startswith(c.compare, label) && c.theta > 0, checks_tau)
    scatter!(ax3, [c.theta for c in cs], [c.max_abs_z for c in cs], marker=mk, label=label)
end
axislegend(ax3, position=:lt, labelsize=9)
save(joinpath(OUT, "b4_theta_reference.png"), fig)
println("saved to $OUT")
