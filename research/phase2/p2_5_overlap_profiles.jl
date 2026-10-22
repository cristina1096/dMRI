# P2.5.9 — summed rate profiles across the gap for overlapping layers (linear, exponential; step h > w/2 as the
# counter-example), from the package (points) and closed form (lines); overlap metrics table.
#
# Run: julia --project=research/baseline research/phase2/p2_5_overlap_profiles.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "overlap.jl"))
using CairoMakie

const OUT = phase2_outdir("p2_5_overlap")
const W, RHO = 2.0, 0.1
grids = Dict("linear" => [0.25, 0.5, 0.6, 0.75, 0.9, 1.0], "exponential" => [0.05, 0.1, 0.25, 0.5, 0.75, 1.0], "step" => [0.5, 0.75])
rows = NamedTuple[]
fig = Figure(size=(1300, 380))
for (j, shape) in enumerate(("linear", "exponential", "step"))
    ax = Axis(fig[1, j], xlabel="x (µm), walls at 0 and $W", ylabel="ΔR₂(x) / (ρ/h)", title=shape == "step" ? "step (counter-example, h > w/2)" : shape)
    for hw in grids[shape]
        h = hw * W
        xs = range(0.0, W, length=400)                    # even count: no sample exactly at the touching point w/2
        lines!(ax, xs, [summed_rate_closed(shape, RHO, h, W, x) / (RHO / h) for x in xs], label="h/w = $hw")
        xp = range(1e-9, W - 1e-9, length=21)            # interior: a spin exactly on a wall counts both sides
        scatter!(ax, xp, [summed_rate_package(shape, RHO, h, W, x) / (RHO / h) for x in xp], markersize=4)
        m = overlap_metrics(shape, h, W)
        push!(rows, (shape=shape, h_over_w=hw, overlap_width_um=m.overlap_width, mid_over_wall=m.mid_over_wall,
            opposite_share=m.opposite_share, realistic=is_realistic(shape, RHO, h, W)))
        @printf("%-11s h/w=%.2f: overlap width %.2f µm, mid/wall %.3e, opposite share %.3e, realistic %s\n", shape, hw, m.overlap_width, m.mid_over_wall, m.opposite_share, rows[end].realistic)
    end
    axislegend(ax, position=:ct, labelsize=8)
end
save(joinpath(OUT, "p2_5_overlap_profiles.png"), fig)
write_csv(joinpath(OUT, "overlap_metrics.csv"), rows)
println("saved to $OUT")
