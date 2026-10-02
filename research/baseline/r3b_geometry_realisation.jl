# P1.3.3 add-on — how much does the random-packing ADC plateau vary between geometry realisations?
#
# The paper notebook (Figure_S1/turtuosity.ipynb) generated its random packings without a seed, so its
# realisation differs from ours (GEOM_SEED = 1 in r3_timestep_plateau.jl). This script measures the
# geometry-to-geometry spread of the plateau ADC to judge the ~3 % offset against the notebook values.
#
# Variables: random_low configuration (density 0.349, radius 1, repeats [50, 50], variance 0.01),
#            geometry seeds 1…NGEOM (default 10), one spin seed per geometry,
#            τ = 1e-3 ms, t = 10 ms, NSPINS = 10_000 initial spins, D = 3.
#
# Run: julia --project=research/baseline -t 8 research/baseline/r3b_geometry_realisation.jl

include(joinpath(@__DIR__, "common.jl"))

const OUT = outdir("r3_timestep_plateau")
const NGEOM = parse(Int, get(ENV, "NGEOM", "10"))
const NSPINS = 10_000
const TAU = 1e-3
const T_EVOLVE = 10.0
const DENSITY = π / 3.0^2 # random_low with r = 1, gap 1

vals = Float64[]
for g in 1:NGEOM
    Random.seed!(g)
    pos, _ = random_positions_radii([50., 50.], DENSITY, 2; mean=1., variance=0.1^2)
    sim = Simulation([]; geometry=Cylinders(radius=1., position=pos, repeats=[50., 50.]), timestep=TAU, verbose=false)
    Random.seed!(1)
    start = mr.get_subset(Snapshot(NSPINS, sim), sim; inside=false)
    final = evolve(start, sim, T_EVOLVE)
    d = [f.position .- s.position for (s, f) in zip(start, final)]
    a = var(vcat([x[1] for x in d], [x[2] for x in d])) / (2T_EVOLVE)
    push!(vals, a)
    @printf("geometry seed %2d: ADC = %.4f\n", g, a)
end
st = seed_stats(vals)
@printf("random_low plateau over %d geometries: %.4f ± %.4f (σ, includes spin noise ≈0.027)\n", NGEOM, st.mean, st.sigma)
write_json(joinpath(OUT, "geometry_realisation.json"), Dict("provenance" => provenance(), "tau_ms" => TAU,
    "geometry_seeds" => collect(1:NGEOM), "ADC" => vals, "mean" => st.mean, "sigma" => st.sigma,
    "notebook_value_at_tau_1e-3" => 2.021457697959048))
