# P2.1.4 — noise floor of the layer model and the spin count for the h sweep.
#
# (1) σ and SEM of S(t) at every readout for h = 0.1 and 0.5 µm (ΔR2(0) = 0.1 /ms), 10 seeds × NSPINS.
#     At h = 0 σ is exactly 0 (B1), so it is not run here.
# (2) Pilot of the P2.3.4(b) sweep over H_GRID with NSEEDS_PILOT seeds: for each neighbouring pair of h, the
#     separation of 1 − S(50) in units of the combined SEM of 10-seed means (pilot σ / √10).
#     The tracker asks for > 10 SEM between neighbours; recommended_nspins = required_spins(...) for the worst pair.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_1_4_noise_floor.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "noise.jl"))

const OUT = phase2_outdir("p2_1_4_noise_floor")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(NSPINS_REF)))
const NSEEDS_PILOT = parse(Int, get(ENV, "NSEEDS_PILOT", "3"))

floor_runs = Dict{String, Any}()
for h in (0.1, 0.5)
    r = run_walls(layer_walls(h=h, rate=RATE_REF); nspins=NSPINS)
    floor_runs[string(h)] = Dict("mean" => r.mean, "sigma" => r.sigma, "sem" => r.sem, "runtime_s" => r.runtime)
    @printf("h = %.2f: S(50) = %.5f, σ(50) = %.2e, SEM(50) = %.2e, max σ over t = %.2e, runtime %.0f s\n",
        h, r.mean[end], r.sigma[end], r.sem[end], maximum(r.sigma), r.runtime)
end

pilot = [run_walls(layer_walls(h=h, rate=RATE_REF); nspins=NSPINS, seeds=1:NSEEDS_PILOT) for h in H_GRID]
att = [1 - p.mean[end] for p in pilot]
sem10 = [p.sigma[end] / sqrt(10) for p in pilot]
pairs = NamedTuple[]
for k in 1:length(H_GRID) - 1
    sep = abs(att[k + 1] - att[k]) / sqrt(sem10[k]^2 + sem10[k + 1]^2 + eps())
    push!(pairs, (h_low=H_GRID[k], h_high=H_GRID[k + 1], attenuation_low=att[k], attenuation_high=att[k + 1], separation_in_sem=sep))
    @printf("  h %.3f → %.3f: 1-S(50) %.5f → %.5f, separation %.1f SEM\n", H_GRID[k], H_GRID[k + 1], att[k], att[k + 1], sep)
end
worst = minimum(p.separation_in_sem for p in pairs)
recommended = required_spins(NSPINS, worst)
@printf("worst neighbouring separation %.1f SEM → recommended spins per seed: %d\n", worst, recommended)
write_csv(joinpath(OUT, "pilot_pairs.csv"), pairs)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "nspins" => NSPINS,
    "nseeds_pilot" => NSEEDS_PILOT, "h_grid" => H_GRID, "rate" => RATE_REF, "times" => TIMES_REF,
    "floor" => floor_runs, "pilot_S50" => [p.mean[end] for p in pilot],
    "worst_separation_in_sem" => worst, "recommended_nspins" => recommended))
println("saved to $OUT")
