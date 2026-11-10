# P3.5.1: overlap configurations of outer layers in arrays (R = 1). For each geometry and outer support L: overlap per
# direction (L > gap/2), the fold at the cell corner (interstitial point), the maximum number of supports covering one
# fluid point (grid), and the fraction of the fluid area covered by ≥ 2 supports. Supports: step/linear L = h ∈ {0.15,
# 0.25, 0.35, 0.5} (h ≤ smallest gap, package rule); exponential L = default cutoff = smallest surface gap (D3.5).
# Run: julia --project=research/baseline research/phase3/p3_5_configs.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_5_configs")
const R1 = 1.0
geoms = (("square 2.5", 2.5, 2.5), ("rect 2.5x4.0", 2.5, 4.0), ("square 4.0", 4.0, 4.0))
"Fraction of the fluid area of the unit cell covered by ≥ 2 supports (n × n midpoint grid)."
function overlap_fraction(sx, sy, L; n=600)
    fluid = over = 0
    for a in 1:n, b in 1:n
        x, y = (a - 0.5) / n * sx - sx / 2, (b - 0.5) / n * sy - sy / 2
        cnt, inside = 0, false
        for i in -2:2, j in -2:2
            r = hypot(x - i * sx, y - j * sy)
            r < R1 && (inside = true; break)
            r - R1 <= L && (cnt += 1)
        end
        inside && continue
        fluid += 1
        cnt >= 2 && (over += 1)
    end
    return over / fluid
end
rows = NamedTuple[]
for (name, sx, sy) in geoms
    gap = min(sx, sy) - 2R1
    supports = vcat([(sh, h, h) for sh in ("step", "linear") for h in (0.15, 0.25, 0.35, 0.5) if h <= gap + 1e-12],
                    [("exponential", NaN, gap)])
    for (sh, h, L) in supports
        (xo, yo, fold) = classify_overlap(sx, sy, R1, L)
        config = fold == 4 ? "interstitial 4-fold" : (xo && yo) ? "two-direction" : (xo || yo) ? "one-direction" : "separated"
        push!(rows, (geometry=name, sx=sx, sy=sy, gap_x=sx - 2R1, gap_y=sy - 2R1, shape=sh, h=h, L=L, x_overlap=xo, y_overlap=yo,
            corner_fold=fold, configuration=config, max_overlap_count=max_overlap_count(sx, sy, R1, L),
            overlap_area_fraction=overlap_fraction(sx, sy, L)))
    end
end
write_csv(joinpath(OUT, "configs.csv"), rows)
foreach(r -> @printf("%-13s %-12s L %.2f  %-20s max count %d  overlap area %.3f\n", r.geometry, r.shape, r.L, r.configuration,
    r.max_overlap_count, r.overlap_area_fraction), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n" => length(rows), "files" => ["configs.csv"]))
println("saved to $OUT")
