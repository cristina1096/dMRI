# P2.3.2 V9 equilibrium density and P2.3.3 layer occupancy, with the layer on (h = 0.5 µm, ΔR2(0) = 0.1).
#
# For τ = default (0.04 ms) and 1e-2 ms, seeds 1–10, NSPINS spins:
#  (1) final positions with the layer == final positions without it (same seeds)  [layer never moves spins]
#  (2) histogram of x mod w (100 bins of 0.02 µm) pooled over seeds: χ² vs uniform, z = (χ² − dof)/√(2 dof), |z| ≤ 3
#  (3) occupancy within h of either face for h ∈ HS: z = (n − N·2h/w)/√(N·p(1−p)), |z| ≤ 2 for every h
# Pass (gate P2.3.3): (1) and (2) and (3) for both τ.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_3_2_density_occupancy.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_3_2_density_occupancy")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(NSPINS_REF)))
const NBINS = 100
const HS = [0.05, 0.1, 0.2, 0.4, 0.5]

summary = Dict{String, Any}()
rows = NamedTuple[]
all_pass = true
for (label, τ) in (("default", nothing), ("tau_1e-2", 1e-2))
    with = run_walls(layer_walls(h=0.5, rate=RATE_REF); timestep=τ, times=[0.0, 50.0], nspins=NSPINS, snapshots=true)
    without = run_walls(Walls(repeats=W_REF); timestep=τ, times=[0.0, 50.0], nspins=NSPINS, snapshots=true)
    identical = [mr.position.(a) for a in with.snapshots] == [mr.position.(b) for b in without.snapshots]
    u = vcat([[mod(p[1], W_REF) for p in mr.position.(s)] for s in with.snapshots]...)
    ntot = length(u)
    hist = zeros(Int, NBINS)
    for x in u
        hist[min(NBINS, floor(Int, x / W_REF * NBINS) + 1)] += 1
    end
    e = ntot / NBINS
    χ2 = sum((hist .- e) .^ 2 ./ e)
    zχ = (χ2 - (NBINS - 1)) / sqrt(2 * (NBINS - 1))
    occ_ok = true
    for h in HS
        p = 2h / W_REF
        n = count(x -> x <= h || W_REF - x <= h, u)
        z = (n - ntot * p) / sqrt(ntot * p * (1 - p))
        occ_ok &= abs(z) <= 2
        push!(rows, (config=label, tau_ms=with.tau, h_um=h, expected_fraction=p, measured_fraction=n / ntot, z=z))
        @printf("  %-8s h=%.2f: fraction %.5f vs %.5f (z %+.2f)\n", label, h, n / ntot, p, z)
    end
    pass = identical && abs(zχ) <= 3 && occ_ok
    global all_pass &= pass
    summary[label] = Dict("tau_ms" => with.tau, "positions_identical_to_no_layer" => identical,
        "histogram_chi2_per_dof" => χ2 / (NBINS - 1), "histogram_z_chi2" => zχ, "occupancy_all_within_2sigma" => occ_ok, "pass" => pass)
    @printf("%-8s: positions identical %s, χ²/dof %.3f (z %+.2f), occupancy ok %s → %s\n",
        label, identical, χ2 / (NBINS - 1), zχ, occ_ok, pass ? "PASS" : "FAIL")
end
write_csv(joinpath(OUT, "occupancy.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "nspins" => NSPINS, "nbins" => NBINS,
    "h_values" => HS, "results" => summary, "pass" => all_pass))
println("saved to $OUT")
