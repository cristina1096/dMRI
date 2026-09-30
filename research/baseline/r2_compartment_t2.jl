# P1.3.2 — R2: compartment T2 in myelinated white matter (paper §3.3, Fig. 8)
#
# Reproduces `se_grid_plot` from mcmr_paper_figures/Figure_8_9/gradient_spin_echo.ipynb
# (spin-echo panel with magnetisation transfer, MT = 5e-3).
# Reference values (tracker): intra-axonal ~80 ms, extra-axonal ~60 ms, myelin ~45 ms.
#
# Variables (copied from the notebook):
#   geometry : Random.seed!(123); random_positions_radii((60, 60), 0.7, 2; mean=1, variance=0.2, min_radius=0.3)
#              Annuli(inner = 0.7·r, outer = r, rotation = :y, repeats = (60, 60), myelin = MYELIN,
#                     relaxation_inner_surface = MT, relaxation_outer_surface = MT)
#   MT (θ_relax on both annulus surfaces) = 5e-3 µm⁻¹·ms^-1/2
#   D = 3 µm²/ms,  R1 = 1/4000 kHz,  R2 (global, CSF) = 1/2000 kHz
#   timestep = (size_scale = 1,) → τ_max = 0.03·1²/3 = 0.01 ms
#   sequences: SpinEcho(TE = 5:5:200 ms, TR = 1000, scanner = Siemens_Prisma, excitation = (phase = 90,)) — instant pulses
#   readout  : one readout per sequence at its TE
#   subsets  : inside = 0 (extra-axonal), 1 (myelin), 2 (intra-axonal)
#   N_spins  : NSPINS (default 1000, as in the notebook) per seed; seeds 1…NSEEDS (default 10)
#   myelin susceptibility: run with MYELIN = false (default) and MYELIN = true (env MYELIN=true)
#
# Apparent T2 per compartment:  T2(TE) = -TE / log(S(TE)/N)  (notebook definition),
# plus a mono-exponential fit of -log(S/N) vs TE over all TEs (slope → 1/T2).
#
# Run: MYELIN=false julia --project=research/baseline -t 8 research/baseline/r2_compartment_t2.jl

include(joinpath(@__DIR__, "common.jl"))
using MRIBuilder

const MYELIN = parse(Bool, get(ENV, "MYELIN", "false"))
const OUT = outdir("r2_compartment_t2_myelin_$(MYELIN)")
const NSPINS = parse(Int, get(ENV, "NSPINS", "1000"))
const NSEEDS = parse(Int, get(ENV, "NSEEDS", "10"))
const MT = 5e-3
const TES = collect(5.0:5.0:200.0)
const COMPARTMENTS = ["extra", "myelin", "intra"]   # inside = 0, 1, 2
const REFERENCE = Dict("intra" => 80.0, "extra" => 60.0, "myelin" => 45.0)

Random.seed!(123)
(positions, radii) = random_positions_radii((60., 60.), 0.7, 2; mean=1., variance=0.2, min_radius=0.3)
geometry = Annuli(; inner=0.7 .* radii, outer=radii, position=positions, rotation=:y, repeats=(60., 60.),
    myelin=MYELIN, relaxation_inner_surface=MT, relaxation_outer_surface=MT)
n_axons = length(radii)
density = sum(π .* radii .^ 2) / 60^2
println("geometry: $n_axons annuli (notebook: 651), outer-radius density $(round(density, digits=4)), myelin=$MYELIN")

sequences = [SpinEcho(TE=t, TR=1000., scanner=MRIBuilder.Siemens_Prisma, excitation=(phase=90.,)) for t in TES]
sim = Simulation(sequences; diffusivity=3., geometry=geometry, R1=1/4000, R2=1/2000, timestep=(size_scale=1.,), verbose=false)
println("τ_max = $(sim.timestep.max_timestep) ms")

# exact S/V per compartment (per unit length along the axon), for the record
Ai = sum(π .* (0.7 .* radii) .^ 2); Ao = sum(π .* radii .^ 2)
Pi = sum(2π .* 0.7 .* radii);      Po = sum(2π .* radii)
SV = Dict("intra" => Pi / Ai, "myelin" => (Pi + Po) / (Ao - Ai), "extra" => Po / (60^2 - Ao))
volfrac = Dict("intra" => Ai / 60^2, "myelin" => (Ao - Ai) / 60^2, "extra" => 1 - Ao / 60^2)
println("S/V (µm⁻¹): ", SV, "\nvolume fractions: ", volfrac)

signal = zeros(NSEEDS, length(TES), 3)
nspins = zeros(Int, NSEEDS, 3)
runtime = Float64[]
for seed in 1:NSEEDS
    Random.seed!(seed)
    t0 = time()
    res = readout(NSPINS, sim; subset=[Subset(inside=0), Subset(inside=1), Subset(inside=2)])
    push!(runtime, time() - t0)
    # res: (n_sequences, n_subsets)
    for c in 1:3, (i, _) in enumerate(TES)
        signal[seed, i, c] = per_spin(res[i, c])
    end
    nspins[seed, :] = [res[1, c].nspins for c in 1:3]
    @printf("seed %2d  (%.1f s)  nspins extra/myelin/intra = %s  T2(TE=50) = %s\n", seed, runtime[end], nspins[seed, :],
        join([@sprintf("%.1f", apparent_T2(50., signal[seed, 10, c])) for c in 1:3], " / "))
end

function fitT2(ratio)
    y = -log.(ratio); X = hcat(ones(length(TES)), TES)
    c = X \ y
    return 1 / c[2]
end

summary = Dict{String, Any}()
rows = NamedTuple[]
for (c, name) in enumerate(COMPARTMENTS)
    T2_TE = [apparent_T2.(TES, signal[s, :, c]) for s in 1:NSEEDS]      # per seed, per TE
    T2_TE_mat = reduce(hcat, T2_TE)'                                     # seeds × TEs
    T2_fit = [fitT2(signal[s, :, c]) for s in 1:NSEEDS]
    st_fit = seed_stats(T2_fit)
    st50 = seed_stats(T2_TE_mat[:, findfirst(==(50.), TES)])
    st100 = seed_stats(T2_TE_mat[:, findfirst(==(100.), TES)])
    ref = REFERENCE[name]
    # tracker criterion for values read off a figure: within 5 %
    within5 = abs(st_fit.mean - ref) / ref <= 0.05
    summary[name] = Dict(
        "reference_T2_ms" => ref, "S_over_V_per_um" => SV[name], "volume_fraction" => volfrac[name],
        "mean_nspins_per_seed" => mean(nspins[:, c]),
        "T2_fit_mean" => st_fit.mean, "T2_fit_sigma" => st_fit.sigma, "T2_fit_sem" => st_fit.sem,
        "T2_TE50_mean" => st50.mean, "T2_TE50_sigma" => st50.sigma,
        "T2_TE100_mean" => st100.mean, "T2_TE100_sigma" => st100.sigma,
        "T2_vs_TE_mean" => vec(mean(T2_TE_mat, dims=1)), "T2_vs_TE_sigma" => vec(std(T2_TE_mat, dims=1)),
        "signal_mean" => vec(mean(signal[:, :, c], dims=1)), "signal_sigma" => vec(std(signal[:, :, c], dims=1)),
        "within_5pct_of_reference" => within5,
        "within_2sigma_of_reference" => abs(st_fit.mean - ref) <= 2 * st_fit.sigma,
    )
    @printf("%-7s T2_fit = %6.1f ± %4.1f ms (σ)  T2(TE=50) = %6.1f ± %4.1f  T2(TE=100) = %6.1f ± %4.1f  ref ≈ %g  → %s\n",
        name, st_fit.mean, st_fit.sigma, st50.mean, st50.sigma, st100.mean, st100.sigma, ref, within5 ? "within 5%" : "outside 5%")
    for s in 1:NSEEDS, (i, te) in enumerate(TES)
        push!(rows, (compartment=name, seed=s, TE_ms=te, signal=signal[s, i, c], nspins=nspins[s, c]))
    end
end
write_csv(joinpath(OUT, "signal.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict(
    "provenance" => provenance(), "compartments" => summary,
    "params" => Dict("myelin" => MYELIN, "MT_theta_relax" => MT, "D" => 3.0, "R1" => 1/4000, "R2_global" => 1/2000,
        "tau_max_ms" => sim.timestep.max_timestep, "binding_constraint" => "tortuosity with size_scale=1 µm (user-set)",
        "n_axons" => n_axons, "outer_density" => density, "repeats_um" => [60, 60], "g_ratio" => 0.7,
        "geometry_seed" => 123, "nspins_per_seed" => NSPINS, "nseeds" => NSEEDS, "TE_ms" => TES,
        "sequence" => "SpinEcho(TE, TR=1000, scanner=Siemens_Prisma, excitation=(phase=90,)), instant pulses"),
    "runtime_s_per_seed" => runtime))

using CairoMakie
fig = Figure(size=(820, 330))
ax1 = Axis(fig[1, 1], xlabel="TE (ms)", ylabel="apparent T2 = -TE/log(S/N) (ms)", yscale=log10,
    title="spin echo, MT=5e-3, myelin=$MYELIN")
ax2 = Axis(fig[1, 2], xlabel="TE (ms)", ylabel="S(TE)/N", yscale=log10)
for name in COMPARTMENTS
    s = summary[name]
    lines!(ax1, TES, s["T2_vs_TE_mean"], label=name, linewidth=3)
    band!(ax1, TES, s["T2_vs_TE_mean"] .- s["T2_vs_TE_sigma"], s["T2_vs_TE_mean"] .+ s["T2_vs_TE_sigma"], alpha=0.3)
    hlines!(ax1, [REFERENCE[name]], linestyle=:dash, color=:gray)
    lines!(ax2, TES, s["signal_mean"], label=name, linewidth=3)
end
ax1.yticks = ([10, 30, 45, 60, 80, 100, 300, 1000], string.([10, 30, 45, 60, 80, 100, 300, 1000]))
ylims!(ax1, 8, 3000)
axislegend(ax1, position=:rt)
save(joinpath(OUT, "r2_compartment_t2.png"), fig)
println("saved to $OUT")
