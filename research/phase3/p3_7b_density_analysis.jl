# P3.7b Task 4 (P3.7.1 arrays, P3.7.2): apparent T₂ against packing density — the paper's Fig. 9B analogue — for the
# intra-axonal (inside), extra-axonal (outside) and summed signals. Monte Carlo primary: p3_7b_run.jl TASK=inside (one
# inside run at R = 1, shared by every density) and TASK=density (outside, square arrays s ∈ {2.3, 2.5, 3.0, 4.0},
# f = πR²/s² = 0.594, 0.503, 0.349, 0.196). ρ = 1/140, D = 3, bulk 1/2000; same centroid c̄ = 0.1 on both sides
# (inside exponential support R, outside support = gap).
# TE = TE_sum*(s): the summed conventional signal f·S_in + (1 − f)·S_out = e⁻¹ (continuum inside, annulus outside);
# all three compartments are read at that TE (exact readout 11 of 0:TE/10:1.5 TE).
# Quantities per density, compartment (intra, extra, sum) and model:
#   T₂_app = −TE/ln S(TE), σ = TE·SEM/(S ln²S);
#   ΔT₂/T₂ = (T₂_layer − T₂_conv)/T₂_conv, σ ≈ TE·SEM_ΔS/(S_layer ln² S_layer)/T₂_conv (paired per compartment; sum via
#   `sum_compartments`: S = f S_in + (1 − f) S_out, SEMs in quadrature — inside and outside runs are independent);
#   2σ two-point noise at SNR 100/200 with r = S_conv(TE) of that compartment;
#   z against support: intra vs the radial continuum (exact); extra and sum vs the annulus ONLY where the annulus has
#   full support (s ≥ 4, validated range; Task 3 + Group A review I1); elsewhere MC-only (annulus gives TE only).
# Dense support curves: intra continuum over s ∈ [2.3, 6] (exact at any s; its TE uses the annulus, allowed); extra/sum
# annulus for s ∈ [4, 6] only — EXTRAPOLATED beyond the validated s = 4, justified by ρR_out/D ≤ 0.008 ≪ 0.067 here
# (Group B review I5).
# Note: the inside run is shared by all densities, so the intra columns and the intra part of the sums are correlated
# across densities.
# Run: julia --project=research/baseline research/phase3/p3_7b_density_analysis.jl  (works on partial outside results)

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
include(joinpath(@__DIR__, "continuum_annulus.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
include(joinpath(@__DIR__, "p3_7_shared.jl"))
include(joinpath(@__DIR__, "p3_7b_mc.jl"))
using CairoMakie
const RES = joinpath(REPO_ROOT, "research", "results", "phase3")
const OUT = phase3_outdir("p3_7b_density_analysis")
const R1, D3, RHO, BULK, CB = 1.0, 3.0, 1 / 140, 1 / 2000, 0.1
const S_FULL = 4.0                                          # annulus full support from this spacing on (validated)
const SHAPES = ("step", "linear", "exponential")
const MODELS = vcat(["conventional"], ["$(s)_c$(CB)" for s in SHAPES])

"Per-model readout tables of one condition label: Dicts by model of vectors over `times` (sorted)."
function tables(curves, label, times)
    sel = [r for r in curves if r["label"] == label]
    pick(m, k) = (p = Dict(num(r, "t") => num(r, k) for r in sel if r["model"] == m); [p[t] for t in times])
    D(k) = Dict(m => pick(m, k) for m in MODELS)
    dD(k) = Dict(m => pick(m, k) for m in MODELS[2:end])
    return (times=times, S_mean=D("S_mean"), S_sem=D("S_sem"), S_sup=D("S_sup"), dS_mean=dD("dS_mean"), dS_sem=dD("dS_sem"),
            dS_sup=dD("dS_sup"), shape=Dict(m => first(r["shape"] for r in sel if r["model"] == m) for m in MODELS))
end
t2(S, te) = apparent_t2(S, te)
t2σ(S, sem, te) = te * sem / (S * log(S)^2)

ins_c = readrows(joinpath(RES, "p3_7b_inside", "conditions.csv"))
ins_k = readrows(joinpath(RES, "p3_7b_inside", "curves.csv"))
den_c = readrows(joinpath(RES, "p3_7b_density", "conditions.csv"))
den_k = readrows(joinpath(RES, "p3_7b_density", "curves.csv"))
inlabel = only(c["label"] for c in ins_c)
rows = NamedTuple[]
for c in den_c
    label, s, f, te = c["label"], num(c, "s"), num(c, "f"), num(c, "te")
    times = sort(unique(num(r, "t") for r in den_k if r["label"] == label))
    k = findfirst(t -> abs(t - te) <= 1e-9 * te, times)
    isnothing(k) && error("$label: no readout at TE")
    out = tables(den_k, label, times)
    ins = tables(ins_k, inlabel, times)                    # inside run holds the union of all density grids
    full = s >= S_FULL
    for comp in ("intra", "extra", "sum"), m in MODELS
        if comp == "sum"
            cm = sum_compartments(f, ins, out, m)
            cc = sum_compartments(f, ins, out, "conventional")
            S, sem, Sc = cm.S[k], cm.S_sem[k], cc.S[k]
            dS, dsem = cm.dS[k], cm.dS_sem[k]
            sup = full ? f * ins.S_sup[m][k] + (1 - f) * out.S_sup[m][k] : NaN
            dsup = (full && m != "conventional") ? f * ins.dS_sup[m][k] + (1 - f) * out.dS_sup[m][k] : NaN
        else
            T = comp == "intra" ? ins : out
            S, sem, Sc = T.S_mean[m][k], T.S_sem[m][k], T.S_mean["conventional"][k]
            dS = m == "conventional" ? NaN : T.dS_mean[m][k]
            dsem = m == "conventional" ? NaN : T.dS_sem[m][k]
            ok = comp == "intra" || full
            sup = ok ? T.S_sup[m][k] : NaN
            dsup = (ok && m != "conventional") ? T.dS_sup[m][k] : NaN
        end
        T2, T2c = t2(S, te), t2(Sc, te)
        layer = m != "conventional"
        push!(rows, (label=label, s=s, f=f, te=te, compartment=comp, model=m, shape=comp == "intra" ? ins.shape[m] : out.shape[m],
            S_mc=S, S_sem=sem, T2_mc=T2, T2_sigma=t2σ(S, sem, te),
            dT2_rel=layer ? (T2 - T2c) / T2c : NaN, dT2_rel_sigma=layer ? t2σ(S, dsem, te) / T2c : NaN,
            noise2_snr200=2 * t2_noise_rel(Sc, 200), noise2_snr100=2 * t2_noise_rel(Sc, 100),
            detectable_snr200=layer && abs((T2 - T2c) / T2c) > 2 * t2_noise_rel(Sc, 200),
            S_sup=sup, z_S=(S - sup) / sem, dS_mc=dS, dS_sem=dsem, dS_sup=dsup, z_dS=(dS - dsup) / dsem,
            T2_sup=isnan(sup) ? NaN : t2(sup, te), support=comp == "intra" ? "radial continuum" : full ? "annulus (full)" : "none (MC only)"))
    end
end
write_csv(joinpath(OUT, "at_te.csv"), rows)
for r in rows
    @printf("%-5s f %.3f %-5s %-18s T₂ %7.2f ± %.2f  ΔT₂/T₂ %s  z_S %s  [%s]\n", r.label, r.f, r.compartment, r.model, r.T2_mc,
        r.T2_sigma, r.model == "conventional" ? "        —        " : @sprintf("%+.4f%% ± %.4f%%", 100r.dT2_rel, 100r.dT2_rel_sigma),
        isnan(r.z_S) ? "—" : @sprintf("%+.1f", r.z_S), r.support)
end

# ---- dense support curves: intra over s ∈ [2.3, 6]; extra and sum on the annulus for s ∈ [4, 6] only ----
dense = NamedTuple[]
for s in sort(unique(vcat(collect(range(2.3, 6.0, length=19)), [S_FULL])))
    f = π * R1^2 / s^2
    Ro = krogh_radius(s); gap = s - 2R1
    Sin = Dict("conventional" => signal_fun(R=R1, D=D3, surface=RHO, bulk=BULK))
    Sout = Dict("conventional" => annulus_signal_fun(R=R1, R_out=Ro, D=D3, surface=RHO, bulk=BULK))
    for sh in SHAPES
        Sin["$(sh)_c$(CB)"] = signal_fun(R=R1, D=D3, layer=(sh, RHO, h_for_cbar(sh, CB, R1)), bulk=BULK)
        ho = h_for_cbar_out(sh, CB, gap)
        Sout["$(sh)_c$(CB)"] = annulus_signal_fun(R=R1, R_out=Ro, D=D3, layer=(sh, RHO, ho, sh == "exponential" ? gap : ho), bulk=BULK)
    end
    te = te_star(t -> f * Sin["conventional"](t) + (1 - f) * Sout["conventional"](t); t_hi=50 / (BULK + RHO * 2π * R1 / (s^2 - π)))
    for m in MODELS
        vals = Dict("intra" => Sin[m](te), "extra" => Sout[m](te), "sum" => f * Sin[m](te) + (1 - f) * Sout[m](te))
        cvals = Dict("intra" => Sin["conventional"](te), "extra" => Sout["conventional"](te),
                     "sum" => f * Sin["conventional"](te) + (1 - f) * Sout["conventional"](te))
        for comp in (s >= S_FULL ? ("intra", "extra", "sum") : ("intra",))
            T2, T2c = t2(vals[comp], te), t2(cvals[comp], te)
            push!(dense, (s=s, f=f, te=te, compartment=comp, model=m, T2=T2, dT2_rel=(T2 - T2c) / T2c))
        end
    end
end
write_csv(joinpath(OUT, "support_dense.csv"), dense)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_densities" => length(den_c),
    "all_identical" => all(c -> c["identical_trajectories"] == "true", vcat(den_c, ins_c)),
    "max_abs_dT2_rel" => maximum((abs(r.dT2_rel) for r in rows if r.model != "conventional"); init=0.0),
    "any_detectable_snr200" => any(r -> r.detectable_snr200, rows), "annulus_full_support_from_s" => S_FULL,
    "files" => ["at_te.csv", "support_dense.csv", "p3_7b_density.png"]))

# ---- figure: Fig. 9B analogue ----
cols = Dict("conventional" => :black, "step" => Makie.wong_colors()[1], "linear" => Makie.wong_colors()[2],
    "exponential" => Makie.wong_colors()[3])
mk = Dict("intra" => :circle, "extra" => :utriangle, "sum" => :rect)
sty = Dict("intra" => :dot, "extra" => :dash, "sum" => :solid)
fig = Figure(size=(1500, 820))
axa = Axis(fig[1, 1], xlabel="intra-axonal volume fraction f = πR²/s²", ylabel="T₂,app at TE_sum* (ms)",
    title="(a) apparent T₂ vs packing, conventional (layers coincide to ≤ 0.03 %, see b)")
axb = Axis(fig[1, 2], xlabel="f", ylabel="|ΔT₂/T₂| vs conventional (%)", yscale=log10,
    title="(b) |layer − conventional|, paired, 2σ bars, lower ends clipped at y/10 (filled ΔT₂ < 0, open > 0)")
axc = Axis(fig[2, 1:2], xlabel="density / compartment / model", ylabel="z = (MC − support)/SEM",
    title="(c) Monte Carlo against support (intra: exact continuum; extra, sum: annulus only at s ≥ 4)")
for comp in ("intra", "extra", "sum"), m in MODELS
    sh = m == "conventional" ? "conventional" : split(m, "_")[1]
    d = sort([q for q in dense if q.compartment == comp && q.model == m], by=q -> q.f)
    m == "conventional" && lines!(axa, [q.f for q in d], [q.T2 for q in d], color=:black, linestyle=sty[comp])
    m == "conventional" || lines!(axb, [q.f for q in d], [max(100abs(q.dT2_rel), 1e-6) for q in d], color=cols[sh], linestyle=sty[comp])
    p = sort([r for r in rows if r.compartment == comp && r.model == m], by=r -> r.f)
    if m == "conventional"
        errorbars!(axa, [r.f for r in p], [r.T2_mc for r in p], [2r.T2_sigma for r in p], color=:black)
        scatter!(axa, [r.f for r in p], [r.T2_mc for r in p], color=:black, marker=mk[comp], markersize=9, label="MC, $(comp)")
    end
    if m != "conventional"
        y = [100abs(r.dT2_rel) for r in p]; e = [200r.dT2_rel_sigma for r in p]
        rangebars!(axb, [r.f for r in p], max.(y .- e, y ./ 10), y .+ e, color=cols[sh])
        scatter!(axb, [r.f for r in p], y, color=[r.dT2_rel < 0 ? cols[sh] : :white for r in p], strokecolor=cols[sh],
            strokewidth=1.5, marker=mk[comp], markersize=9, label=comp == "sum" ? sh : nothing)
    end
end
n200 = 100 * 2 * t2_noise_rel(exp(-1), 200)
hlines!(axb, [n200], color=:gray, linestyle=:dash, label="2σ noise SNR 200 (sum, r = e⁻¹)")
axislegend(axa, position=:lt, labelsize=9)
Legend(fig[1, 3], axb, labelsize=9)
Label(fig[1, 3, Bottom()], "markers: ● intra, ▲ extra, ■ sum (MC)\nlines: dotted intra (continuum, any s),\ndashed extra, solid sum (annulus: validated\nat s = 4, extrapolated to s = 6)", fontsize=10)
zr = [r for r in rows if !isnan(r.z_S)]
scatter!(axc, 1:length(zr), [r.z_S for r in zr], color=[cols[r.model == "conventional" ? "conventional" : split(r.model, "_")[1]] for r in zr],
    marker=[mk[r.compartment] for r in zr], markersize=8)
hlines!(axc, [-2, 2], color=:gray, linestyle=:dash)
axc.xticks = (1:length(zr), ["$(r.label) $(r.compartment) $(r.model == "conventional" ? "conv" : split(r.model, "_")[1])" for r in zr])
axc.xticklabelrotation = π / 2.5
axc.xticklabelsize = 7
save(joinpath(OUT, "p3_7b_density.png"), fig)
println("saved to $OUT")
