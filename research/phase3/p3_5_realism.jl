# P3.5.2: realism of the summed outer-layer rate in 2D arrays. Criterion (as on walls, P2.5.9): realistic ⟺ the
# maximum of the summed rate over the fluid region ≤ (1 + 1e-9)·ΔR₂(0), the single-cylinder surface value
# ΔR₂(0) = (ρ/(hκ_out))·g(0). Evaluated on an 801 × 801 grid of the unit cell plus exact special points: the gap
# midpoints, the cell corner (interstitial point), and the surface points facing a neighbour. Overlapping layers add.
# Geometries: square 2.5, rectangle 2.5 × 4.0, square 4.0 (R = 1); step/linear h ∈ {0.15, 0.25, 0.35, 0.5} (≤ gap);
# exponential support L = smallest gap, decay length h ∈ the same list. ρ = 1 (the ratio does not depend on ρ).
# Run: julia --project=research/baseline -t 8 research/phase3/p3_5_realism.jl

include(joinpath(@__DIR__, "harness.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_5_realism")
const R1, N = 1.0, 801
geoms = (("square 2.5", 2.5, 2.5), ("rect 2.5x4.0", 2.5, 4.0), ("square 4.0", 4.0, 4.0))
cases = NamedTuple[]
for (name, sx, sy) in geoms
    gap = min(sx, sy) - 2R1
    for sh in ("step", "linear", "exponential"), h in (0.15, 0.25, 0.35, 0.5)
        L = sh == "exponential" ? gap : h
        L <= gap + 1e-12 || continue
        push!(cases, (geometry=name, sx=sx, sy=sy, shape=sh, h=h, L=L))
    end
end
results = Vector{NamedTuple}(undef, length(cases))
maps = Dict{Int, Matrix{Float64}}()
lk = ReentrantLock()
Threads.@threads for k in eachindex(cases)
    c = cases[k]
    outer = (shape=c.shape, h=c.h, L=c.L, pre=side_prefactor(1.0, c.shape, c.h, c.L, R1, false))
    surface = outer.pre * g_profile(c.shape, 0.0, c.L / c.h)
    rate(x, y) = brute_rate_profile(x, y; R=R1, sx=c.sx, sy=c.sy, inner=nothing, outer=outer)
    fluid(x, y) = hypot(x - round(x / c.sx) * c.sx, y - round(y / c.sy) * c.sy) >= R1 - 1e-12
    xs, ys = range(-c.sx / 2, c.sx / 2, length=N), range(-c.sy / 2, c.sy / 2, length=N)
    M = [fluid(x, y) ? rate(x, y) / surface : NaN for x in xs, y in ys]
    special = [(c.sx / 2, 0.0), (0.0, c.sy / 2), (c.sx / 2, c.sy / 2), (R1, 0.0), (0.0, R1), (-R1, 0.0), (0.0, -R1)]
    sp = [rate(x, y) / surface for (x, y) in special if fluid(x, y)]
    gmax = maximum(filter(!isnan, M))
    imax = argmax(replace(M, NaN => -Inf))
    mx = max(gmax, maximum(sp))
    results[k] = merge(c, (surface_rate=surface, max_ratio=mx, grid_max_ratio=gmax, x_at_max=xs[imax[1]], y_at_max=ys[imax[2]],
        gap_mid_ratio=rate(c.sx / 2, 0.0) / surface, corner_ratio=fluid(c.sx / 2, c.sy / 2) ? rate(c.sx / 2, c.sy / 2) / surface : NaN,
        facing_surface_ratio=rate(R1, 0.0) / surface, realistic=mx <= 1 + 1e-9))
    lock(lk) do
        maps[k] = M
    end
end
write_csv(joinpath(OUT, "realism.csv"), results)
foreach(r -> @printf("%-13s %-12s h %.2f L %.2f  max/surface %.4f  gap mid %.4f  corner %s  facing surface %.4f  %s\n", r.geometry, r.shape,
    r.h, r.L, r.max_ratio, r.gap_mid_ratio, isnan(r.corner_ratio) ? "—" : @sprintf("%.4f", r.corner_ratio), r.facing_surface_ratio,
    r.realistic ? "realistic" : "NOT realistic"), results)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_cases" => length(results),
    "n_unrealistic" => count(r -> !r.realistic, results), "criterion" => "max summed rate ≤ (1 + 1e-9)·single surface value",
    "files" => ["realism.csv", "p3_5_realism.png"]))

# figure: ratio maps for the overlapping cases at the largest h (and the 4-fold exponential)
sel = [k for k in eachindex(results) if (results[k].h == 0.5 || (results[k].shape == "exponential" && results[k].geometry == "square 4.0" && results[k].h == 0.5))
       && results[k].max_ratio > 0 && !(results[k].geometry == "square 4.0" && results[k].shape != "exponential")]
fig = Figure(size=(360 * min(length(sel), 4), 340 * cld(length(sel), 4)))
for (i, k) in enumerate(sel)
    r = results[k]
    ax = Axis(fig[cld(i, 4), mod1(i, 4)], aspect=DataAspect(), title="$(r.geometry), $(r.shape), h = $(r.h)\nmax/surface = $(round(r.max_ratio, digits=3))")
    hm = heatmap!(ax, range(-r.sx / 2, r.sx / 2, length=N), range(-r.sy / 2, r.sy / 2, length=N), maps[k], colorrange=(0, max(1.0, r.max_ratio)),
        colormap=:viridis, nan_color=:white)
    scatter!(ax, [r.x_at_max], [r.y_at_max], color=:red, marker=:xcross, markersize=10)
    hidedecorations!(ax)
end
save(joinpath(OUT, "p3_5_realism.png"), fig)
println("saved to $OUT")
