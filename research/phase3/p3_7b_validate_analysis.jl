# P3.7b Task 3: validation of the annulus continuum (Krogh cell, R_out = s/√π) against paired Monte Carlo outside
# cylinders in square arrays (p3_7b_run.jl TASK=validate: s ∈ {2.5, 4.0}, ρ ∈ {0.02, 0.2}, D = 3, no bulk; conventional
# + step/linear/exponential at the same centroid c̄ = 0.1, exponential support = gap).
# Quantities:
#   e_S(t) = (S_mc − S_ann)/S_ann and z_S = (S_mc − S_ann)/SEM_S per model and readout (t > 0); max_t |e_S|;
#   ΔS ratio ΔS_mc/ΔS_ann at the readout of max_t |ΔS_ann|, σ = SEM_ΔS/|ΔS_ann| (paired), and z_ΔS;
#   effective rates λ_eff = −ln S(t_end)/t_end (endpoint, not the plan's log-linear fit; reported only), MC against
#   annulus; the kept fraction of the layer in the annulus.
# Decision rule (plan, fixed before the run):
#   (S criterion, ruling before any data: |S_mc − S_ann| ≤ 1 %·S_ann + 2·SEM_S at every readout — the plan's bare
#    max|e_S| ≤ 1 % would fail on MC noise alone where S ≈ 0.01 late in the ρ = 0.2 windows)
#   full    — at every condition the S criterion holds for all models AND every ΔS ratio within [0.8, 1.2] within 2σ
#             (ratio + 2σ ≥ 0.8 and ratio − 2σ ≤ 1.2): annulus gives dense curves, TE*, windows, fit family, ceiling;
#   partial — only the conventional model meets the S criterion: TE*, windows, ceiling (approximate) only;
#   failure — otherwise: TE and windows only; every outside result MC-only.
# After the run (user decision 2026-11-09): the annulus layer is renormalised by 1/kept (continuum_annulus.jl), recomputed
# here (the curves' S_sup column holds the cut, unscaled annulus of the run and is kept as `*_cut` for transparency);
# the same rule is applied PER SPACING: full support where it passes, TE/readout windows only where it fails.
# Run: julia --project=research/baseline research/phase3/p3_7b_validate_analysis.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_annulus.jl"))
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
        s_, ρ_, sh, h_ = num(c, "s"), num(c, "rho"), p[1]["shape"], num(p[1], "h")
        ann(layer) = annulus_signal(R=1.0, R_out=krogh_radius(s_), D=num(c, "D"), times=t, R2_bulk=num(c, "R2_bulk"); layer...)
        Sconv_a = ann((surface=ρ_,))
        S, sem = [num(r, "S_mean") for r in p], [num(r, "S_sem") for r in p]
        Sa = m == "conventional" ? Sconv_a : ann((layer=(sh, ρ_, h_, sh == "exponential" ? s_ - 2 : h_),))
        Sa_cut = [num(r, "S_sup") for r in p]
        eS = (S .- Sa) ./ Sa
        zS = (S .- Sa) ./ max.(sem, eps())
        k = argmax(abs.(eS[2:end])) + 1
        layer = m != "conventional"
        if layer
            dS, dsem, dSa = [num(r, "dS_mean") for r in p], [num(r, "dS_sem") for r in p], Sa .- Sconv_a
            dSa_cut = [num(r, "dS_sup") for r in p]
            j = argmax(abs.(dSa))
            ratio_cut = dS[argmax(abs.(dSa_cut))] / dSa_cut[argmax(abs.(dSa_cut))]
            ratio, rσ, zd = dS[j] / dSa[j], dsem[j] / abs(dSa[j]), (dS[j] - dSa[j]) / dsem[j]
        else
            dS = dsem = dSa = fill(NaN, length(t)); ratio = rσ = zd = ratio_cut = NaN; j = 1
        end
        push!(rows, (label=label, s=num(c, "s"), rho=num(c, "rho"), tau=num(c, "tau"), model=m, shape=p[1]["shape"], h=num(p[1], "h"),
            cbar=num(p[1], "cbar"), kept=num(p[1], "kept"), max_abs_eS=abs(eS[k]), eS_at=eS[k], zS_at=zS[k], t_max_eS=t[k],
            max_abs_zS=maximum(abs.(zS[2:end])), max_abs_dS_abs=maximum(abs.(S[2:end] .- Sa[2:end])),
            max_abs_eS_cut=maximum(abs.((S[2:end] .- Sa_cut[2:end]) ./ Sa_cut[2:end])), dS_ratio_cut=ratio_cut, lambda_eff_mc=-log(S[end]) / t[end], lambda_eff_ann=-log(Sa[end]) / t[end],
            dS_ratio=ratio, dS_ratio_sigma=rσ, z_dS=zd, t_max_dS=t[j], maxdS_ann=layer ? abs(dSa[j]) : NaN,
            S_ok=all(abs.(S[2:end] .- Sa[2:end]) .<= ETOL .* Sa[2:end] .+ 2 .* sem[2:end]),
            identical=c["identical_trajectories"] == "true"))
        series[(label, m)] = (t=t, eS=eS, dS=dS, dsem=dsem, dSa=dSa)
    end
end
write_csv(joinpath(OUT, "validation.csv"), rows)
for r in rows
    @printf("%-12s %-18s kept %.3f  S ok %-5s max|e_S| %.2e (z %+5.1f)  λ_eff mc/ann %.4f/%.4f  ΔS ratio %s\n", r.label, r.model, r.kept, r.S_ok,
        r.max_abs_eS, r.zS_at, r.lambda_eff_mc, r.lambda_eff_ann,
        isnan(r.dS_ratio) ? "—" : @sprintf("%.3f ± %.3f (z %+4.1f, max|ΔS_ann| %.1e)", r.dS_ratio, r.dS_ratio_sigma, r.z_dS, r.maxdS_ann))
end

# ---- decision: the pre-registered rule, applied per spacing (user decision 2026-11-09) ----
ratio_in(r) = r.dS_ratio + 2r.dS_ratio_sigma >= RLO && r.dS_ratio - 2r.dS_ratio_sigma <= RHI
dec = Dict{String, Any}()
for sp in sort(unique(r.s for r in rows))
    rs = [r for r in rows if r.s == sp]
    conv_ok = all(r.S_ok for r in rs if r.model == "conventional")
    layer_S_ok = all(r.S_ok for r in rs if r.model != "conventional")
    ratio_ok = all(ratio_in(r) for r in rs if r.model != "conventional")
    d = conv_ok && layer_S_ok && ratio_ok ? "full" : conv_ok ? "partial" : "failure"
    dec["s$(sp)"] = Dict("decision" => d, "conventional_S_ok" => conv_ok, "layer_S_ok" => layer_S_ok, "dS_ratio_ok" => ratio_ok,
                         "max_abs_S_diff" => maximum(r.max_abs_dS_abs for r in rs))
    @printf("s = %.1f: decision %s (conventional S ok %s; layers S ok %s; ΔS ratios ok %s; max |S_mc − S_ann| %.1e)\n",
        sp, d, conv_ok, layer_S_ok, ratio_ok, maximum(r.max_abs_dS_abs for r in rs))
end
decision = join(["$(k): $(v["decision"])" for (k, v) in sort(collect(dec), by=first)], ", ")
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "decision_per_spacing" => dec,
    "annulus" => "renormalised by 1/kept (user decision 2026-11-09); *_cut columns: the cut, unscaled annulus of the run",
    "all_identical" => all(r -> r.identical, rows),
    "rule" => "per spacing — full: |S_mc − S_ann| ≤ 0.01 S_ann + 2 SEM at every readout for all models and ΔS ratio in [0.8,1.2] within 2σ; partial: conventional only; else failure (TE/windows only)",
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
