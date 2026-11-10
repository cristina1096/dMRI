# P3.5.4 (plan 2026-11-11-overlapping-layers-arrays.md, Task 4): signal against the overlap degree at fixed ρ.
# MC only (p3_5_run.jl TASK=overlap): square 2.5 (two-direction overlap) and rectangle 2.5 × 4.0 (one-direction),
# R = 1, D = 3, no bulk, ρR/D = 1; conventional + step/linear/exponential at h ∈ {0.15, 0.25, 0.35, 0.5}
# (h/gap_x = 0.3 separated, 0.5 touching, 0.7, 1.0 overlapping; exponential support fixed at the gap 0.5).
# Quantities: max_t |ΔS(t)|, ΔS = S_layer − S_conv (paired, SEM at the arg-max readout), against h/gap and against the
# fluid-area fraction covered by ≥ 2 supports (P3.5.1, configs.csv). Monotonicity in h per profile: steps between
# neighbouring h with σ = √(σ₁² + σ₂²) — conservative, because the h models share trajectories but only per-model SEMs
# are stored. Separated reference: the 3.7b map at s = 4, ρR/D = 1 (h/R 0.25, 0.5; no overlap for step/linear; its
# exponential support is 2, overlapping) — a different geometry (S/V 0.49 vs 2.02 µm⁻¹), shown for orientation only.
# Checks: identical trajectories, ΔR₂,max·τ = 0.005 (summed bound), conventional θ√τ.
# Run: julia --project=research/baseline research/phase3/p3_5_overlap_analysis.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "p3_7_shared.jl"))
using CairoMakie
const RES = joinpath(REPO_ROOT, "research", "results", "phase3")
const OUT = phase3_outdir("p3_5_overlap_analysis")
const SHAPES = ("step", "linear", "exponential")

conds = readrows(joinpath(RES, "p3_5_overlap", "conditions.csv"))
curves = readrows(joinpath(RES, "p3_5_overlap", "curves.csv"))
cfg = readrows(joinpath(RES, "p3_5_configs", "configs.csv"))
geomname(sx, sy) = sx == sy ? "square $(sx)" : "rect $(sx)x$(sy)"
function overlap_fraction(sx, sy, sh, h)
    L = sh == "exponential" ? min(sx, sy) - 2.0 : h
    r = [c for c in cfg if c["geometry"] == geomname(sx, sy) && c["shape"] == sh && isapprox(num(c, "L"), L)]
    return isempty(r) ? NaN : num(first(r), "overlap_area_fraction")
end
rows = NamedTuple[]
for c in conds
    label, sx, sy, ρ, τ = c["label"], num(c, "sx"), num(c, "sy"), num(c, "rho"), num(c, "tau")
    gx = sx - 2.0
    for m in unique(r["model"] for r in curves if r["label"] == label)
        m == "conventional" && continue
        p = sort([r for r in curves if r["label"] == label && r["model"] == m], by=r -> num(r, "t"))
        dS = [num(r, "dS_mean") for r in p]; dsem = [num(r, "dS_sem") for r in p]
        k = argmax(abs.(dS))
        sh, h = p[1]["shape"], num(p[1], "h")
        push!(rows, (label=label, sx=sx, sy=sy, model=m, shape=sh, h=h, h_over_gap=h / gx, overlap_fraction=overlap_fraction(sx, sy, sh, h),
            maxdS=abs(dS[k]), maxdS_sem=dsem[k], sign=sign(dS[k]), t_max=num(p[k], "t"), tau=τ,
            max_rate_tau=num(c, "max_rate_tau"), theta_sqrt_tau=ρ * sqrt(π / num(c, "D")) * sqrt(τ),
            identical=c["identical_trajectories"] == "true"))
    end
end
write_csv(joinpath(OUT, "overlap.csv"), rows)
steps = NamedTuple[]
for l in unique(r.label for r in rows), sh in SHAPES
    p = sort([r for r in rows if r.label == l && r.shape == sh], by=r -> r.h)
    for i in 1:length(p) - 1
        d = p[i + 1].maxdS - p[i].maxdS
        σ = sqrt(p[i].maxdS_sem^2 + p[i + 1].maxdS_sem^2)
        push!(steps, (label=l, shape=sh, h_from=p[i].h, h_to=p[i + 1].h, change=d, sigma_unpaired=σ, z=d / σ,
            relative_change=d / p[i].maxdS))
    end
end
write_csv(joinpath(OUT, "steps.csv"), steps)
for r in rows
    @printf("%-20s %-18s h/gap %.2f overlap area %.3f  max|ΔS| %.4f ± %.4f (sign %+d)\n", r.label, r.model, r.h_over_gap, r.overlap_fraction,
        r.maxdS, r.maxdS_sem, Int(r.sign))
end
for s in steps
    @printf("step %-20s %-12s h %.2f → %.2f: Δ max|ΔS| %+.4f ± %.4f (z %+.1f, %+.0f %%)\n", s.label, s.shape, s.h_from, s.h_to,
        s.change, s.sigma_unpaired, s.z, 100s.relative_change)
end
# separated reference (3.7b map, s = 4, ρR/D = 1)
ref = NamedTuple[]
mref = joinpath(RES, "p3_7b_map_analysis", "mc_map.csv")
if isfile(mref)
    for r in readrows(mref)
        (r["label"] == "s4.0_rhoRD1.0" && num(r, "hR") in (0.25, 0.5)) || continue
        push!(ref, (shape=r["shape"], h=num(r, "hR"), maxdS=num(r, "maxdS_mc"), maxdS_sem=num(r, "maxdS_sem")))
    end
end
write_csv(joinpath(OUT, "separated_reference_s4.csv"), isempty(ref) ? [(note="missing",)] : ref)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "all_identical" => all(r -> r.identical, rows),
    "max_rate_tau" => maximum(r.max_rate_tau for r in rows), "files" => ["overlap.csv", "steps.csv", "separated_reference_s4.csv", "p3_5_overlap.png"]))

# figure
cols = Dict("step" => Makie.wong_colors()[1], "linear" => Makie.wong_colors()[2], "exponential" => Makie.wong_colors()[3])
fig = Figure(size=(1300, 520))
for (j, (xk, xl)) in enumerate(((:h_over_gap, "h / gap_x"), (:overlap_fraction, "fluid fraction covered by ≥ 2 supports")))
    ax = Axis(fig[1, j], xlabel=xl, ylabel="max_t |ΔS| (fixed ρ, ρR/D = 1)",
        title=j == 1 ? "(a) against h/gap (dashed: touching, h = gap/2)" : "(b) against the overlap area (P3.5.1)")
    j == 1 && vlines!(ax, [0.5], color=:gray, linestyle=:dash)
    for l in unique(r.label for r in rows), sh in SHAPES
        p = sort([r for r in rows if r.label == l && r.shape == sh], by=r -> r.h)
        x = [getproperty(r, xk) for r in p]
        errorbars!(ax, x, [r.maxdS for r in p], [2r.maxdS_sem for r in p], color=cols[sh])
        scatterlines!(ax, x, [r.maxdS for r in p], color=cols[sh], linestyle=startswith(l, "sq") ? :solid : :dash,
            marker=startswith(l, "sq") ? :circle : :utriangle, label="$(sh), $(startswith(l, "sq") ? "square 2.5" : "rect 2.5×4")")
    end
    j == 1 && Legend(fig[1, 3], ax, labelsize=9, merge=true)
end
save(joinpath(OUT, "p3_5_overlap.png"), fig)
println("saved to $OUT")
