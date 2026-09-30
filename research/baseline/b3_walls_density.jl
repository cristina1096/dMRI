# B3 — equilibrium spin density and layer occupancy between walls (baseline for P2.3.2 V9, P2.3.3)
#
# The R₂(d) layer relaxes a spin only while it is within λ (= h) of a wall face, so the relaxation
# it produces depends on the spin density near the wall. If the old reflection code distorts that
# density, every layer result inherits it. This run checks, with the UNMODIFIED simulator and all
# relaxation off, that
#   (a) the density profile across the gap is flat (histogram of x mod w), and
#   (b) the fraction of spins within λ of each face equals λ/w,
# at every readout time and at several timesteps. Reflection errors show up within about one
# step length √(2Dτ) of a wall; at the default τ = 0.04 ms that is 0.49 µm, i.e. longer than every
# λ tested, so any artefact found here is inherited by all Phase 2 layer results.
#
# Variables:
#   geometry : Walls(repeats=w), w = 2 µm, normal along x; u = mod(x, w) ∈ [0, w) is the position
#              within one gap; distance to the lower face is u, to the upper face w − u
#   D = 3 µm²/ms; all relaxation, permeability, off-resonance, MT off
#   sequence : none (empty sequence, no RF); positions only
#   timesteps: default (0.04 ms, tortuosity), 1e-2, 1e-3, 1e-4 ms (see CONFIGS for N_spins and times)
#   λ        : 0.05, 0.1, 0.2, 0.4 µm
#   seeds    : 1…NSEEDS (default 10); same seed → same initial positions for every τ
#   histogram: NBINS (default 100) bins of w/NBINS = 0.02 µm, pooled over seeds
#
# Pass criteria (per τ, per readout time):
#   occupancy : pooled count within λ of a face vs N·λ/w, z = (n − Np)/√(Np(1−p)); |z| ≤ 2
#   histogram : χ² vs uniform; z_χ² = (χ² − dof)/√(2·dof); |z_χ²| ≤ 3; and the two bins touching
#               each wall within 3σ (Poisson)
# With ~200 occupancy comparisons, about 5% are expected outside 2σ by chance; a systematic
# failure is one that repeats across seeds / readouts / λ with the same sign.
#
# Run: julia --project=research/baseline -t 8 research/baseline/b3_walls_density.jl

include(joinpath(@__DIR__, "common.jl"))

const OUT = outdir("b3_walls_density")
const NSEEDS = parse(Int, get(ENV, "NSEEDS", "10"))
const NBINS = parse(Int, get(ENV, "NBINS", "100"))
const D = 3.0
const W = 2.0
const LAMBDAS = [0.05, 0.1, 0.2, 0.4]
const FULL_TIMES = [0.0, 10.0, 20.0, 30.0, 40.0, 50.0]

# (label, timestep (nothing = default), N_spins per seed, readout times)
const CONFIGS = [
    ("default", nothing, parse(Int, get(ENV, "NSPINS", "100000")), FULL_TIMES),
    ("tau_1e-2", 1e-2, parse(Int, get(ENV, "NSPINS", "100000")), FULL_TIMES),
    ("tau_1e-3", 1e-3, parse(Int, get(ENV, "NSPINS_1E3", "100000")), parse.(Float64, split(get(ENV, "TIMES_1E3", "0,10,20,30,40,50"), ","))),
    ("tau_1e-4", 1e-4, parse(Int, get(ENV, "NSPINS_1E4", "100000")), parse.(Float64, split(get(ENV, "TIMES_1E4", "0,1,2,5"), ","))),
]

geometry = Walls(repeats=W)

function occupancy_stats(counts, ntot, p)
    e = ntot * p
    s = sqrt(ntot * p * (1 - p))
    return (expected=e, sigma_binom=s, z=(sum(counts) - e) / s)
end

function run_config(label, timestep, nspins, times, io)
    kw = isnothing(timestep) ? (;) : (; timestep=timestep)
    sim = Simulation([mr.SequenceParts.empty_sequence()]; geometry=geometry, diffusivity=D, verbose=false, kw...)
    τ = sim.timestep.max_timestep
    tee(io, @sprintf("\n== %s: τ = %.3g ms, step length √(2Dτ) = %.3g µm, N_spins = %d, seeds = %d, times = %s",
        label, τ, sqrt(2D * τ), nspins, NSEEDS, string(times)))

    nt = length(times)
    hist = zeros(Int, NBINS, nt)                              # pooled over seeds
    occ = zeros(Int, NSEEDS, nt, length(LAMBDAS), 2)          # (seed, time, λ, face: 1 = lower, 2 = upper)
    runtime = 0.0
    for seed in 1:NSEEDS
        Random.seed!(seed)
        snap = Snapshot(nspins, sim, 500)
        t0 = time()
        res = readout(snap, sim, times; return_snapshot=true)
        runtime += time() - t0
        for (it, s) in enumerate(res)
            u = [mod(p[1], W) for p in mr.position.(s)]
            for x in u
                hist[min(NBINS, floor(Int, x / W * NBINS) + 1), it] += 1
            end
            for (il, λ) in enumerate(LAMBDAS)
                occ[seed, it, il, 1] = count(<(λ), u)
                occ[seed, it, il, 2] = count(x -> W - x < λ, u)
            end
        end
        @printf("  seed %2d done (%.1f s)\n", seed, time() - t0)
    end
    ntot = nspins * NSEEDS

    occ_rows = NamedTuple[]
    n_out_2σ = 0
    n_occ = 0
    worst_occ = 0.0
    for (it, t) in enumerate(times), (il, λ) in enumerate(LAMBDAS), face in 1:2
        p = λ / W
        st = occupancy_stats(occ[:, it, il, face], ntot, p)
        frac = sum(occ[:, it, il, face]) / ntot
        seed_sigma = std(occ[:, it, il, face] ./ nspins)
        n_occ += 1
        abs(st.z) > 2 && (n_out_2σ += 1)
        worst_occ = abs(st.z) > abs(worst_occ) ? st.z : worst_occ
        push!(occ_rows, (config=label, tau_ms=τ, t_ms=t, lambda_um=λ, face=face == 1 ? "lower" : "upper",
            expected_fraction=p, measured_fraction=frac, rel_dev=frac / p - 1, z=st.z,
            sigma_frac_binomial=sqrt(p * (1 - p) / nspins), sigma_frac_seeds=seed_sigma))
    end

    hist_rows = NamedTuple[]
    chi_rows = NamedTuple[]
    e = ntot / NBINS
    hist_pass = true
    for (it, t) in enumerate(times)
        h = hist[:, it]
        χ² = sum((h .- e) .^ 2 ./ e)
        dof = NBINS - 1
        zχ = (χ² - dof) / sqrt(2dof)
        wall_z = [(h[1] - e) / sqrt(e), (h[end] - e) / sqrt(e)]
        ok = abs(zχ) <= 3 && all(abs.(wall_z) .<= 3)
        hist_pass &= ok
        push!(chi_rows, (config=label, tau_ms=τ, t_ms=t, chi2=χ², dof=dof, z_chi2=zχ,
            z_bin_lower_wall=wall_z[1], z_bin_upper_wall=wall_z[2], pass=ok))
        for b in 1:NBINS
            push!(hist_rows, (config=label, tau_ms=τ, t_ms=t, bin_lo_um=(b - 1) * W / NBINS, bin_hi_um=b * W / NBINS,
                count=h[b], expected=e, ratio=h[b] / e))
        end
        tee(io, @sprintf("  t = %5.1f ms: χ²/dof = %.3f (z = %+.2f); wall bins z = %+.2f / %+.2f → %s",
            t, χ² / dof, zχ, wall_z[1], wall_z[2], ok ? "flat" : "NOT flat"))
    end
    for r in occ_rows
        tee(io, @sprintf("    t=%5.1f λ=%.2f %-5s: fraction %.5f vs %.5f (rel %+.2e, z %+.2f)",
            r.t_ms, r.lambda_um, r.face, r.measured_fraction, r.expected_fraction, r.rel_dev, r.z))
    end
    occ_pass = n_out_2σ <= max(1, ceil(Int, 0.1 * n_occ))
    tee(io, @sprintf("  occupancy: %d of %d comparisons outside 2σ (≈%.1f expected by chance), worst z = %+.2f → %s",
        n_out_2σ, n_occ, 0.0455 * n_occ, worst_occ, occ_pass ? "PASS" : "CHECK"))
    tee(io, @sprintf("  histogram → %s;  runtime %.1f s", hist_pass ? "PASS" : "CHECK", runtime))

    return (
        rows_occ=occ_rows, rows_hist=hist_rows, rows_chi=chi_rows,
        summary=Dict("config" => label, "tau_ms" => τ, "step_length_um" => sqrt(2D * τ),
            "nspins_per_seed" => nspins, "nseeds" => NSEEDS, "readout_times_ms" => times,
            "occupancy_outside_2sigma" => n_out_2σ, "occupancy_comparisons" => n_occ,
            "occupancy_worst_z" => worst_occ, "occupancy_pass" => occ_pass,
            "histogram_pass" => hist_pass, "runtime_s" => runtime),
    )
end

all_occ = NamedTuple[]
all_hist = NamedTuple[]
all_chi = NamedTuple[]
summaries = Any[]
open(joinpath(OUT, "run.log"), "w") do io
    for (label, ts, n, times) in CONFIGS
        r = run_config(label, ts, n, times, io)
        append!(all_occ, r.rows_occ)
        append!(all_hist, r.rows_hist)
        append!(all_chi, r.rows_chi)
        push!(summaries, r.summary)
    end
end
write_csv(joinpath(OUT, "occupancy.csv"), all_occ)
write_csv(joinpath(OUT, "histogram.csv"), all_hist)
write_csv(joinpath(OUT, "histogram_chi2.csv"), all_chi)
write_json(joinpath(OUT, "summary.json"), Dict(
    "provenance" => provenance(), "results" => summaries,
    "geometry" => "Walls(repeats=$W): planes x = k·$W, normal along x",
    "physics" => "D=$D, all relaxation/permeability/off-resonance/MT off",
    "lambdas_um" => LAMBDAS, "nbins" => NBINS, "bin_width_um" => W / NBINS,
))

# figure: density profile (last readout, pooled) per τ, and occupancy z vs λ
using CairoMakie
fig = Figure(size=(1100, 620))
for (i, s) in enumerate(summaries)
    label = s["config"]
    ax = Axis(fig[1, i], title=@sprintf("%s (τ = %.0e ms)", label, s["tau_ms"]), xlabel="u = x mod w (µm)",
        ylabel=i == 1 ? "count / expected" : "", limits=(0, W, 0.9, 1.1))
    for t in s["readout_times_ms"]
        rows = filter(r -> r.config == label && r.t_ms == t, all_hist)
        lines!(ax, [(r.bin_lo_um + r.bin_hi_um) / 2 for r in rows], [r.ratio for r in rows], label="t = $(t) ms")
    end
    ne = s["nspins_per_seed"] * s["nseeds"] / NBINS
    band!(ax, [0, W], fill(1 - 2 / sqrt(ne), 2), fill(1 + 2 / sqrt(ne), 2), color=(:gray, 0.25))
    i == 1 && axislegend(ax, position=:cb, labelsize=9)

    ax2 = Axis(fig[2, i], xlabel="λ (µm)", ylabel=i == 1 ? "occupancy z" : "", limits=(0, 0.45, -4, 4))
    band!(ax2, [0, 0.45], [-2, -2], [2, 2], color=(:gray, 0.25))
    for face in ("lower", "upper")
        rows = filter(r -> r.config == label && r.face == face, all_occ)
        scatter!(ax2, [r.lambda_um for r in rows] .+ (face == "lower" ? -0.007 : 0.007), [r.z for r in rows],
            marker=face == "lower" ? :circle : :utriangle, label=face)
    end
    i == 1 && axislegend(ax2, position=:lb, labelsize=9)
end
save(joinpath(OUT, "b3_walls_density.png"), fig)
println("saved to $OUT")
