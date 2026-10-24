# P3.1.3 follow-up: is the near-surface occupancy deviation found in C-B3 (inside below, outside above the uniform
# expectation by ≈ 0.5 % at the default τ = 0.01 ms) a timestep effect of reflecting finite steps on a curved surface?
#
# Occupancy within h of the surface, relative deviation δ = fraction/p − 1 with σ_δ = √((1 − p)/(N p)) (binomial,
# pooled spins), for the single cylinder (inside, R = 1) and the square array (outside, s = 2.5), at t = 0 (control: the
# initial positions are uniform by construction) and t = 50 ms, τ = 4e-2, 1e-2, 2.5e-3 ms. Seeds 1–10, 100000 spins.
# A deviation that shrinks with τ (∝ √τ or τ) is a property of the discrete walk at finite τ, not of the layer code.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_1_3b_density_tau.jl

include(joinpath(@__DIR__, "harness.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_1_3b_density_tau")
const HS = [0.02, 0.05, 0.1, 0.2]
const TAUS = [4e-2, 1e-2, 2.5e-3]

radii(r, j; s) = vcat([[radial_distance(sp.position; sx=s, sy=s) for (k, sp) in enumerate(f.spins) if m[k]]
                       for (f, m) in zip(r.snaps[j], r.masks)]...)

"run_cylinders variant keeping the snapshots at t = 0 and t = 50 (index 1 and 2)."
function run_keep(geometry; region, s, bounding_box, timestep)
    sim = Simulation([mr.SequenceParts.empty_sequence()]; geometry=geometry, diffusivity=D_REF, verbose=false, timestep=timestep)
    snaps = [Snapshot[], Snapshot[]]
    masks = BitVector[]
    for seed in SEEDS_REF
        Random.seed!(seed)
        snap = Snapshot(NSPINS_REF, sim, bounding_box; transverse=1.0, longitudinal=0.0)
        r0 = [radial_distance(sp.position; sx=s, sy=s) for sp in snap.spins]
        push!(masks, BitVector(region == :inside ? r0 .< R_REF : r0 .> R_REF))
        res = readout(snap, sim, [0.0, 50.0]; return_snapshot=true)
        push!(snaps[1], res[1]); push!(snaps[2], res[2])
    end
    return (snaps=snaps, masks=masks)
end

rows = NamedTuple[]
bb = mr.BoundingBox([-R_REF, -R_REF, -R_REF], [R_REF, R_REF, R_REF])
for τ in TAUS
    for (region, g, s, b) in ((:inside, Cylinders(radius=R_REF), Inf, bb), (:outside, Cylinders(radius=R_REF, repeats=[S_REF, S_REF]), S_REF, 500))
        r = run_keep(g; region=region, s=s, bounding_box=b, timestep=τ)
        for (j, t) in ((1, 0.0), (2, 50.0))
            rr = radii(r, j; s=s)
            N = length(rr)
            for h in HS
                p = region == :inside ? shell_area_inside(R_REF, h) : shell_area_outside(R_REF, h, S_REF, S_REF)
                n = region == :inside ? count(x -> x >= R_REF - h, rr) : count(x -> x <= R_REF + h, rr)
                δ = n / (N * p) - 1
                σ = sqrt((1 - p) / (N * p))
                push!(rows, (region=String(region), tau=τ, t=t, h=h, p=p, N=N, delta=δ, sigma=σ, z=δ / σ))
                @printf("%-7s τ = %.1e t = %4.0f h = %.2f: δ = %+.2e ± %.1e (z %+.2f)\n", region, τ, t, h, δ, σ, δ / σ)
            end
        end
    end
end
write_csv(joinpath(OUT, "occupancy_tau.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rows" => rows))

fig = Figure(size=(900, 380))
for (k, region) in enumerate(("inside", "outside"))
    ax = Axis(fig[1, k], xlabel="τ (ms)", ylabel="occupancy / uniform − 1 (t = 50 ms)", xscale=log10, title=region)
    for h in HS
        q = filter(r -> r.region == region && r.t == 50.0 && r.h == h, rows)
        errorbars!(ax, [r.tau for r in q], [r.delta for r in q], [r.sigma for r in q])
        scatterlines!(ax, [r.tau for r in q], [r.delta for r in q], label="h = $h")
    end
    hlines!(ax, [0.0], color=:black)
    axislegend(ax, labelsize=9)
end
save(joinpath(OUT, "p3_1_3b_density_tau.png"), fig)
println("saved to $OUT")
