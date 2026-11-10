# P3.5.2: realism of the summed outer-layer rate in 2D arrays. Criterion (tracker P3.5.2, as on walls P2.5.9 "never
# higher mid-gap"; ruling 2026-11-11 — the plan's first wording compared with the SINGLE-cylinder surface value): realistic
# ⟺ no interior maximum above the surface value. In 2D this is tested along rays (the faithful generalisation of the 1D
# walls test): from every surface point along its outward normal until another surface or the cell edge, the summed rate
# never exceeds its value at that surface point (excess ≤ 1e-9). A single global max/max comparison fails here, because a
# neighbour's support boundary (circle R + L around the neighbour) can cut into the gap, so a truncated profile switches on
# inside the fluid where the local surface value does not include it. Also reported (context): the global interior and
# surface maxima (surface value held over an arc of two adjacent samples), and the peak
# relative to the single-cylinder surface value ΔR₂(0) = (ρ/(hκ_out))·g(0) (how much summing raises the peak; relevant to
# the τ bound), and whether a violation is confined to a touching point (no overlap of positive area).
# Evaluated on an 801 × 801 grid of the unit cell, 3600 points on the surface, and exact special points: the gap
# midpoints and the cell corner (interstitial point). Overlapping layers add.
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
    special = [(c.sx / 2, 0.0), (0.0, c.sy / 2), (c.sx / 2, c.sy / 2)]
    interior(x, y) = hypot(x - round(x / c.sx) * c.sx, y - round(y / c.sy) * c.sy) > R1 + 1e-9
    Mi = [interior(x, y) ? M[a, b] : NaN for (a, x) in enumerate(xs), (b, y) in enumerate(ys)]
    sp = [rate(x, y) / surface for (x, y) in special if interior(x, y)]
    imax = argmax(replace(Mi, NaN => -Inf))
    interior_max = max(maximum(filter(!isnan, Mi)), maximum(sp; init=0.0))
    # surface value held over an arc (two adjacent samples): an isolated point where a neighbour's support ends exactly on
    # the surface (step/linear at h = gap, inclusive endpoint) does not count as the surface value
    sv = [rate((R1 + 1e-12) * cos(θ), (R1 + 1e-12) * sin(θ)) / surface for θ in range(0, 2π, length=3601)]
    surf_max = maximum(min(sv[i], sv[i + 1]) for i in 1:length(sv) - 1)
    mx = max(interior_max, surf_max)
    # ray test: 1800 surface points, outward normal rays in steps of 2e-4 µm up to 2·max(gap) or until another cylinder
    rmax_ray = 2 * max(c.sx, c.sy)
    excess, nviol, θworst = 0.0, 0, 0.0
    θs = range(0, 2π, length=1801)[1:end-1]
    for θ in θs
        u = (cos(θ), sin(θ))
        s0 = rate((R1 + 1e-12) * u[1], (R1 + 1e-12) * u[2])   # just outside: cos/sin rounding can put R·u at r < R
        best = s0
        for t in 2e-4:2e-4:rmax_ray
            x, y = (R1 + t) * u[1], (R1 + t) * u[2]
            fluid(x, y) || break
            v = rate(x, y)
            best = max(best, v)
        end
        e = s0 > 0 ? best / s0 - 1 : (best > 0 ? Inf : 0.0)
        if e > 1e-9
            nviol += 1
            e > excess && (θworst = θ)
        end
        excess = max(excess, e)
    end
    touching_only = interior_max > surf_max * (1 + 1e-9) &&
                    all(!(rate(x, y) / surface > surf_max * (1 + 1e-9)) for (a, x) in enumerate(xs), (b, y) in enumerate(ys) if interior(x, y) && !(abs(abs(x) - c.sx / 2) < 1e-9 || abs(abs(y) - c.sy / 2) < 1e-9))
    results[k] = merge(c, (surface_rate=surface, max_ratio=mx, interior_max_ratio=interior_max, surface_max_ratio=surf_max,
        x_at_max=xs[imax[1]], y_at_max=ys[imax[2]],
        gap_mid_ratio=rate(c.sx / 2, 0.0) / surface, corner_ratio=fluid(c.sx / 2, c.sy / 2) ? rate(c.sx / 2, c.sy / 2) / surface : NaN,
        facing_surface_ratio=rate(R1, 0.0) / surface, ray_max_excess=excess, ray_violation_fraction=nviol / length(θs),
        ray_worst_angle_deg=rad2deg(θworst), realistic=excess <= 1e-9,
        violation_touching_point_only=touching_only))
    lock(lk) do
        maps[k] = M
    end
end
write_csv(joinpath(OUT, "realism.csv"), results)
foreach(r -> @printf("%-13s %-12s h %.2f L %.2f  interior max %.4f  surface max %.4f (× single)  gap mid %.4f  corner %s  %s\n", r.geometry, r.shape,
    r.h, r.L, r.interior_max_ratio, r.surface_max_ratio, r.gap_mid_ratio, isnan(r.corner_ratio) ? "—" : @sprintf("%.4f", r.corner_ratio),
    r.realistic ? "realistic" : @sprintf("NOT realistic (ray excess %.3f on %.1f %% of the surface)", r.ray_max_excess, 100r.ray_violation_fraction)), results)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_cases" => length(results),
    "n_unrealistic" => count(r -> !r.realistic, results),
    "criterion" => "ray test: along every outward surface normal the summed rate never exceeds its surface value (excess ≤ 1e-9)",
    "files" => ["realism.csv", "p3_5_realism.png"]))

# figure: ratio maps for the overlapping cases at the largest h (and the 4-fold exponential)
sel = [k for k in eachindex(results) if (results[k].h == 0.5 || (results[k].shape == "exponential" && results[k].geometry == "square 4.0" && results[k].h == 0.5))
       && results[k].max_ratio > 0 && !(results[k].geometry == "square 4.0" && results[k].shape != "exponential")]
fig = Figure(size=(360 * min(length(sel), 4), 340 * cld(length(sel), 4)))
for (i, k) in enumerate(sel)
    r = results[k]
    ax = Axis(fig[cld(i, 4), mod1(i, 4)], aspect=DataAspect(), title="$(r.geometry), $(r.shape), h = $(r.h)\npeak/single surface $(round(r.max_ratio, digits=3)), ray excess $(round(r.ray_max_excess, digits=3))")
    hm = heatmap!(ax, range(-r.sx / 2, r.sx / 2, length=N), range(-r.sy / 2, r.sy / 2, length=N), maps[k], colorrange=(0, max(1.0, r.max_ratio)),
        colormap=:viridis, nan_color=:white)
    scatter!(ax, [r.x_at_max], [r.y_at_max], color=:red, marker=:xcross, markersize=10)
    hidedecorations!(ax)
end
save(joinpath(OUT, "p3_5_realism.png"), fig)
println("saved to $OUT")
