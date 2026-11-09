# P3.7b Task 3: validation of the annulus continuum (Krogh cell, R_out = s/√π) against paired Monte Carlo outside
# cylinders in square arrays (p3_7b_run.jl TASK=validate: s ∈ {2.5, 4.0}, ρ ∈ {0.02, 0.2}, D = 3, no bulk; conventional
# + step/linear/exponential at the same centroid c̄ = 0.1, exponential support = gap).
# Quantities:
#   e_S(t) = (S_mc − S_ann)/S_ann and z_S = (S_mc − S_ann)/SEM_S per model and readout (t > 0); max_t |e_S|;
#   ΔS ratio ΔS_mc/ΔS_ann at the readout of max_t |ΔS_ann|, σ = SEM_ΔS/|ΔS_ann| (paired), and z_ΔS;
#   effective rates λ_eff = −ln S(t_end)/t_end, MC against annulus; the kept fraction of the layer in the annulus.
# Decision rule (plan, fixed before the run):
#   full    — at every condition max|e_S| ≤ 1 % for all models AND every ΔS ratio within [0.8, 1.2] within 2σ
#             (ratio + 2σ ≥ 0.8 and ratio − 2σ ≤ 1.2): annulus gives dense curves, TE*, windows, fit family, ceiling;
#   partial — only the conventional model meets max|e_S| ≤ 1 %: TE*, windows, ceiling (approximate) only;
#   failure — otherwise: TE and windows only; every outside result MC-only.
# Run: julia --project=research/baseline research/phase3/p3_7b_validate_analysis.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "p3_7_shared.jl"))
using CairoMakie
const SRC = joinpath(REPO_ROOT, "research", "results", "phase3", "p3_7b_validate")
const OUT = phase3_outdir("p3_7b_validate_analysis")
const ETOL, RLO, RHI = 0.01, 0.8, 1.2

conds = readrows(joinpath(SRC, "conditions.csv"))
curves = readrows(joinpath(SRC, "curves.csv"))
rows = NamedTuple[]
series = Dict{Tuple{String, String}, NamedTuple}()
for c in conds
    label = c["label"]
    sel = [r for r in curves if r["label"] == label]
    for m in unique(r["model"] for r in sel)
        p = sort([r for r in sel if r["model"] == m], by=r -> num(r, "t"))
        t = [num(r, "t") for r in p]
        S, sem, Sa = [num(r, "S_mean") for r in p], [num(r, "S_sem") for r in p], [num(r, "S_sup") for r in p]
        eS = (S .- Sa) ./ Sa
        zS = (S .- Sa) ./ max.(sem, eps())
        k = argmax(abs.(eS[2:end])) + 1
        layer = m != "conventional"
        if layer
            dS, dsem, dSa = [num(r, "dS_mean") for r in p], [num(r, "dS_sem") for r in p], [num(r, "dS_sup") for r in p]
            j = argmax(abs.(dSa))
            ratio, rσ, zd = dS[j] / dSa[j], dsem[j] / abs(dSa[j]), (dS[j] - dSa[j]) / dsem[j]
        else
            dS = dsem = dSa = fill(NaN, length(t)); ratio = rσ = zd = NaN; j = 1
        end
        push!(rows, (label=label, s=num(c, "s"), rho=num(c, "rho"), tau=num(c, "tau"), model=m, shape=p[1]["shape"], h=num(p[1], "h"),
            cbar=num(p[1], "cbar"), kept=num(p[1], "kept"), max_abs_eS=abs(eS[k]), eS_at=eS[k], zS_at=zS[k], t_max_eS=t[k],
            max_abs_zS=maximum(abs.(zS[2:end])), lambda_eff_mc=-log(S[end]) / t[end], lambda_eff_ann=-log(Sa[end]) / t[end],
            dS_ratio=ratio, dS_ratio_sigma=rσ, z_dS=zd, t_max_dS=t[j], maxdS_ann=layer ? abs(dSa[j]) : NaN,
            identical=c["identical_trajectories"] == "true"))
        series[(label, m)] = (t=t, eS=eS, dS=dS, dsem=dsem, dSa=dSa)
    end
end
write_csv(joinpath(OUT, "validation.csv"), rows)
for r in rows
    @printf("%-12s %-18s kept %.3f  max|e_S| %.2e (z %+5.1f)  λ_eff mc/ann %.4f/%.4f  ΔS ratio %s\n", r.label, r.model, r.kept,
        r.max_abs_eS, r.zS_at, r.lambda_eff_mc, r.lambda_eff_ann,
        isnan(r.dS_ratio) ? "—" : @sprintf("%.3f ± %.3f (z %+4.1f, max|ΔS_ann| %.1e)", r.dS_ratio, r.dS_ratio_sigma, r.z_dS, r.maxdS_ann))
end

# ---- decision ----
conv_ok = all(r.max_abs_eS <= ETOL for r in rows if r.model == "conventional")
layer_S_ok = all(r.max_abs_eS <= ETOL for r in rows if r.model != "conventional")
ratio_ok = all(r.dS_ratio + 2r.dS_ratio_sigma >= RLO && r.dS_ratio - 2r.dS_ratio_sigma <= RHI for r in rows if r.model != "conventional")
decision = conv_ok && layer_S_ok && ratio_ok ? "full" : conv_ok ? "partial" : "failure"
@printf("decision: %s (conventional max|e_S| ≤ 1 %%: %s; layers: %s; ΔS ratios in [0.8, 1.2] within 2σ: %s)\n",
    decision, conv_ok, layer_S_ok, ratio_ok)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "decision" => decision,
    "conventional_eS_ok" => conv_ok, "layer_eS_ok" => layer_S_ok, "dS_ratio_ok" => ratio_ok,
    "max_abs_eS" => maximum(r.max_abs_eS for r in rows), "all_identical" => all(r -> r.identical, rows),
    "rule" => "full: max|e_S| ≤ 0.01 all models and ΔS ratio in [0.8,1.2] within 2σ; partial: conventional only; else failure",
    "files" => ["validation.csv", "p3_7b_validate.png"]))

# ---- figure: ΔS(t) MC vs annulus (top), e_S(t) per model (bottom), one column per condition ----
cols = Dict("conventional" => :black, "step" => Makie.wong_colors()[1], "linear" => Makie.wong_colors()[2],
    "exponential" => Makie.wong_colors()[3])
labels = [c["label"] for c in conds]
fig = Figure(size=(380 * length(labels), 760))
for (i, l) in enumerate(labels)
    ax = Axis(fig[1, i], xlabel="t (ms)", ylabel="ΔS = S_layer − S_conv", title="$(l): ΔS, MC (points) vs annulus (lines)")
    ax2 = Axis(fig[2, i], xlabel="t (ms)", ylabel="(S_mc − S_ann)/S_ann", title="$(l): e_S")
    for r in (q for q in rows if q.label == l)
        sr = series[(l, r.model)]
        if r.model != "conventional"
            lines!(ax, sr.t, sr.dSa, color=cols[r.shape], label=r.shape)
            errorbars!(ax, sr.t, sr.dS, 2 .* sr.dsem, color=cols[r.shape])
            scatter!(ax, sr.t, sr.dS, color=cols[r.shape], markersize=6)
        end
        scatterlines!(ax2, sr.t[2:end], sr.eS[2:end], color=cols[r.shape], markersize=5, label=r.model == "conventional" ? "conventional" : r.shape)
    end
    hspan!(ax2, -ETOL, ETOL, color=(:gray, 0.2))
    i == 1 && axislegend(ax, position=:lb, labelsize=9)
    i == 1 && axislegend(ax2, position=:lb, labelsize=9)
end
Label(fig[3, :], "decision: $(decision) (rule fixed before the run; grey band ±1 %; bars 2σ paired)", fontsize=12)
save(joinpath(OUT, "p3_7b_validate.png"), fig)
println("saved to $OUT")
