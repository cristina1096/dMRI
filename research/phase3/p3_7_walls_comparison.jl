# P3.7.4 inside part (final review C1): walls (Phase 2 Review C1 follow-up map, w = 2, D = 3) against cylinders inside
# (P3.7a map, R = 1, D = 3) with ONE threshold convention: `threshold_status` (log-log interpolation of the first
# crossing of 0.01) on each continuum map, applied to max|ΔS| at fixed ρ and to the best-fit residual among valid fits
# (walls: not at the grid edge; cylinders: not at_edge and not over the ceiling). Thresholds are converted to physical
# ρ (µm/ms): walls ρ = x·D/w, cylinders ρ = y·D/R, and compared at equal absolute h (µm) where both maps have it.
# Sources: results/phase2/p2_5_continuum/map.csv, results/phase3/p3_7_map_analysis/continuum_map.csv.
# Run: julia --project=research/baseline research/phase3/p3_7_walls_comparison.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
include(joinpath(@__DIR__, "p3_7_shared.jl"))
const OUT = phase3_outdir("p3_7_walls_comparison")
const W, R1, D3, LEVEL = 2.0, 1.0, 3.0, 0.01
walls = readrows(joinpath(REPO_ROOT, "research", "results", "phase2", "p2_5_continuum", "map.csv"))
cyl = readrows(joinpath(REPO_ROOT, "research", "results", "phase3", "p3_7_map_analysis", "continuum_map.csv"))

"Interpolated threshold in the dimensionless strength x, as (value, status) — `> x_max` reported as (x_max, below_range)."
function thr(rows, xkey, ykey, valid)
    sel = sort([r for r in rows if valid(r)], by=r -> num(r, xkey))
    length(sel) < 2 && return (NaN, :no_valid_points)
    return threshold_status([num(r, xkey) for r in sel], [num(r, ykey) for r in sel], LEVEL)
end
fmt(v, st, scale) = st == :crossed ? @sprintf("%.2f", v * scale) : st == :below_range ? @sprintf("> %.2f", v * scale) :
                    st == :above_range ? "< range" : "—"

out = NamedTuple[]
for s in ("step", "linear", "exponential"), h in (0.1, 0.5, 1.0)
    hw, hr = h / W, h / R1
    wr = [r for r in walls if r["shape"] == s && isapprox(num(r, "h_over_w"), hw)]
    cr = [r for r in cyl if r["shape"] == s && isapprox(num(r, "hR"), hr)]
    (isempty(wr) || isempty(cr)) && continue
    (wf, wfs) = thr(wr, "rho_w_over_D", "max_abs_dS_fixed_rho", r -> true)
    (wb, wbs) = thr(wr, "rho_w_over_D", "best_fit_residual", r -> r["fit_at_edge"] == "false")
    (cf, cfs) = thr(cr, "rhoRD", "maxdS", r -> true)
    (cb, cbs) = thr(cr, "rhoRD", "residual", r -> r["at_edge"] == "false" && r["ceiling"] == "false")
    push!(out, (shape=s, h_um=h, walls_fixed_rho=wf * D3 / W, walls_fixed_status=String(wfs), cyl_fixed_rho=cf * D3 / R1,
        cyl_fixed_status=String(cfs), walls_bestfit_rho=wb * D3 / W, walls_bestfit_status=String(wbs),
        cyl_bestfit_rho=cb * D3 / R1, cyl_bestfit_status=String(cbs)))
    @printf("%-12s h = %.1f µm  fixed ρ: walls %-8s cyl %-8s   best fit: walls %-8s cyl %-8s (µm/ms)\n", s, h,
        fmt(wf, wfs, D3 / W), fmt(cf, cfs, D3 / R1), fmt(wb, wbs, D3 / W), fmt(cb, cbs, D3 / R1))
end
write_csv(joinpath(OUT, "thresholds_rho.csv"), out)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "level" => LEVEL, "w" => W, "R" => R1, "D" => D3,
    "convention" => "threshold_status (log-log interpolated first crossing) on both continuum maps; values in µm/ms"))
println("saved to $OUT")
