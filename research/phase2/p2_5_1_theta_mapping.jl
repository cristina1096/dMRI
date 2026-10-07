# P2.5 prerequisite — the conventional model's relaxivity is ρ = θ·√(D/π).
#
# MCMR multiplies M⊥ by e^(−θ√τ) at each wall collision; Gaussian steps hit a wall at n·√(D/(πτ)) per area and time,
# so E[−ln M_base(t)] = 2·θ·√(D/π)·t/w for a uniform start (the mean-exponent identity of the baseline). With
# θ = ρ·√(π/D) this equals the layer's 2ρt/w (P2.4.5). Checked at D = 1, 3, 6 and τ = 1e-2, 4e-2 (w = 2, ρ = 0.01,
# t = 50 ms, seeds 1–10). Pass: |z| ≤ 3 for the baseline at every (D, τ). Also reported: the paired difference of the
# mean exponents, baseline − exponential layer (h = 0.2) on identical trajectories.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_5_1_theta_mapping.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_5_1_theta_mapping")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(NSPINS_REF)))
const RHO, T = 0.01, 50.0
mean_expo(r) = [mean(-log(s.orientations[1].transverse) for s in snap.spins) for snap in r.snapshots]
positions(r) = [mr.position.(snap) for snap in r.snapshots]

rows = NamedTuple[]
for D in (1.0, 3.0, 6.0), τ in (1e-2, 4e-2)
    θ = theta_for(RHO; D=D)
    b = run_walls(baseline_walls(theta=θ); D=D, timestep=τ, times=[0.0, T], nspins=NSPINS, snapshots=true)
    l = run_walls(rho_walls(rho=RHO, h=0.2, shape="exponential"); D=D, timestep=τ, times=[0.0, T], nspins=NSPINS, snapshots=true)
    eb, el = mean_expo(b), mean_expo(l)
    expected = 2RHO * T / W_REF
    zb = (mean(eb) - expected) / (std(eb) / sqrt(length(eb)))
    dpair = eb .- el
    zp = mean(dpair) / (std(dpair) / sqrt(length(dpair)))
    same_paths = positions(b) == positions(l)
    push!(rows, (D=D, tau_ms=τ, theta=θ, expected=expected, base_mean=mean(eb), base_sem=std(eb) / sqrt(length(eb)), z_base=zb,
        layer_mean=mean(el), paired_diff=mean(dpair), paired_sem=std(dpair) / sqrt(length(dpair)), z_paired=zp,
        same_paths=same_paths, pass=abs(zb) <= 3))
    @printf("D = %.0f, τ = %.0e: θ = %.5f; baseline E[−ln M] = %.6f ± %.1e (2ρt/w = %.6f, z = %+.2f); paired base − layer = %+.1e ± %.0e (z = %+.1f); identical paths: %s → %s\n",
        D, τ, θ, mean(eb), std(eb) / sqrt(length(eb)), expected, zb, mean(dpair), std(dpair) / sqrt(length(dpair)), zp, same_paths, abs(zb) <= 3 ? "PASS" : "FAIL")
end
write_csv(joinpath(OUT, "mapping.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rho" => RHO, "T_ms" => T, "nspins" => NSPINS,
    "rows" => rows, "pass" => all(r.pass for r in rows), "paths_identical" => all(r.same_paths for r in rows)))
println("saved to $OUT")
