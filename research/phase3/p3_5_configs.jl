# P3.5.1: overlap configurations of outer layers in arrays (R = 1). For each geometry and outer support L: overlap per
# direction (L > gap/2), the fold at the cell corner (interstitial point), the maximum number of supports covering one
# fluid point (grid), the fraction of the fluid area covered by ≥ 2 supports, the fraction covered by ≥ 1 support
# (coverage; Group B review I1), and the largest distance of a fluid point from every surface. Supports: step/linear L = h ∈ {0.15,
# 0.25, 0.35, 0.5} (h ≤ smallest gap, package rule); exponential L = default cutoff = smallest surface gap (D3.5).
# Run: julia --project=research/baseline research/phase3/p3_5_configs.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_5_configs")
const R1 = 1.0
geoms = (("square 2.5", 2.5, 2.5), ("rect 2.5x4.0", 2.5, 4.0), ("square 4.0", 4.0, 4.0))
"Fractions of the fluid area of the unit cell covered by ≥ 2 and by ≥ 1 supports, and the largest distance of a fluid
point from every surface (n × n midpoint grid)."
function coverage(sx, sy, L; n=600)
    fluid = over = cov = 0
    dfar = 0.0
    for a in 1:n, b in 1:n
        x, y = (a - 0.5) / n * sx - sx / 2, (b - 0.5) / n * sy - sy / 2
        cnt, inside, dmin = 0, false, Inf
        for i in -2:2, j in -2:2
            r = hypot(x - i * sx, y - j * sy)
            r < R1 && (inside = true; break)
            r - R1 <= L && (cnt += 1)
            dmin = min(dmin, r - R1)
        end
        inside && continue
        fluid += 1
        cnt >= 2 && (over += 1)
        cnt >= 1 && (cov += 1)
        dfar = max(dfar, dmin)
    end
    return (over=over / fluid, covered=cov / fluid, dfar=dfar)
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
            (cv = coverage(sx, sy, L); (overlap_area_fraction=cv.over, covered_fraction=cv.covered, farthest_from_surface=cv.dfar))...))
    end
end
write_csv(joinpath(OUT, "configs.csv"), rows)
foreach(r -> @printf("%-13s %-12s L %.2f  %-20s max count %d  overlap area %.3f  covered %.3f  farthest %.3f\n", r.geometry, r.shape,
    r.L, r.configuration, r.max_overlap_count, r.overlap_area_fraction, r.covered_fraction, r.farthest_from_surface), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n" => length(rows), "files" => ["configs.csv"]))
println("saved to $OUT")
