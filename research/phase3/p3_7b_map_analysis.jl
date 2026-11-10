# P3.7b Task 5 (P3.7.3 outside): dimensionless map OUTSIDE cylinders in square arrays, s ∈ {2.5, 4.0}; R = 1, D = 3,
# no bulk; ρR/D ∈ {0.1, 1, 10} × h/R ∈ {0.1, 0.25, 0.5} × outer step/linear/exponential (+ conventional) per condition;
# 41 readouts to ln 100/λ₁ (annulus conventional). Monte Carlo primary (p3_7b_run.jl TASK=map); 5 seeds at s = 4,
# ρR/D ∈ {0.1, 10}, 10 elsewhere (user decision 2026-11-10).
# Support (user decisions 2026-11-09/10):
#   s = 2.5 — annulus failed validation: MC only (no dense map, no best fit, no ceiling);
#   s = 4   — "map MC decides": each condition validates the annulus with the pre-registered rule (every model:
#             |S_mc − S_ann| ≤ 1 %·S_ann + 2·SEM at every readout; every layer: ΔS ratio at the readout of max|ΔS_ann| in
#             [0.8, 1.2] within 2σ). Pass → the annulus gives the dense curves, the best-fit family and the ceiling at
#             that ρR/D. The dense map is shown only up to the largest ρR/D below which every condition passed (Group A
#             validated ρR/D ≤ 0.067).
# Quantities per layer model: fixed ρ — max_t |ΔS| (paired MC, SEM at the arg-max readout);
#   best fit (s = 4, passing conditions) — ρ_c = argmin max_t |S_layer,MC(t) − S_conv,ann(t; ρ_c)|; flags at_edge and
#   ceiling (λ₁ of the layer in the annulus > `annulus_ceiling`). Flagged fits report no ρ_c and no residual; their
#   residual at the search edge is reported separately as an UPPER BOUND on the distance (3.7a lesson I1). Noise floor:
#   max_t SEM(S_layer). Thresholds: `threshold_status` at 0.01 in ρR/D (one convention, 3.7a lesson C1).
# Checks: conventional θ√τ and MC − annulus z (s = 4); annulus grid at h/R = 0.1, ρR/D = 1 and 30 (N 2000 vs 4000); fit
# family N = 1000 vs 2000 at the most demanding VALID fit (3.7a lesson I4). c̄_min = 0.033 (linear h/R = 0.1) < 0.05:
# the τ rule is extrapolated there (flagged).
# Run: julia --project=research/baseline -t 8 research/phase3/p3_7b_map_analysis.jl  (works on partial results)

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_annulus.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
include(joinpath(@__DIR__, "p3_7_shared.jl"))
using CairoMakie
const SRC = joinpath(REPO_ROOT, "research", "results", "phase3", "p3_7b_map")
const OUT = phase3_outdir("p3_7b_map_analysis")
const R1, D3, LEVEL, ETOL = 1.0, 3.0, 0.01, 0.01
const SHAPES = ("step", "linear", "exponential")
const NFAM = 1000
supp(sh, h, s) = sh == "exponential" ? s - 2R1 : h
ann_fun(s; kw...) = annulus_signal_fun(R=R1, R_out=krogh_radius(s), D=D3; kw...)
conv_family(s, times; N=NFAM) = ρc -> ann_fun(s; surface=ρc, N=N).(times)
function fit_flags(S_layer, times, ρ, sh, h, s; N=NFAM)
    (ρc, res, edge) = best_fit_conventional(S_layer, times, conv_family(s, times; N=N); rho0=ρ)
    λ1 = annulus_rates(R=R1, R_out=krogh_radius(s), D=D3, layer=(sh, ρ, h, supp(sh, h, s)))[1]
    ceil = λ1 > annulus_ceiling(R1, krogh_radius(s), D3)
    ok = !edge && !ceil
    return (rho_c=ok ? ρc : NaN, residual=ok ? res : NaN, at_edge=edge, ceiling=ceil, lambda1=λ1,
            residual_edge_bound=ok ? NaN : res, rho_c_raw=ρc)
end

conds = readrows(joinpath(SRC, "conditions.csv"))
curves = readrows(joinpath(SRC, "curves.csv"))
series(label, m, k) = (p = sort([r for r in curves if r["label"] == label && r["model"] == m], by=r -> num(r, "t")); [num(r, k) for r in p])

# ---- per condition: validation (s = 4) and MC quantities ----
mc = NamedTuple[]
condrows = NamedTuple[]
for c in conds
    label, s, ρ, τ = c["label"], num(c, "s"), num(c, "rho"), num(c, "tau")
    x = ρ * R1 / D3
    times = sort(unique(num(r, "t") for r in curves if r["label"] == label))
    models = unique(r["model"] for r in curves if r["label"] == label)
    S_ok = Dict{String, Bool}()
    ratio = Dict{String, Tuple{Float64, Float64}}()
    for m in models
        S, sem, Sa = series(label, m, "S_mean"), series(label, m, "S_sem"), series(label, m, "S_sup")
        S_ok[m] = all(abs.(S[2:end] .- Sa[2:end]) .<= ETOL .* Sa[2:end] .+ 2 .* sem[2:end])
        if m != "conventional"
            dS, dsem, dSa = series(label, m, "dS_mean"), series(label, m, "dS_sem"), series(label, m, "dS_sup")
            j = argmax(abs.(dSa))
            ratio[m] = (dS[j] / dSa[j], dsem[j] / abs(dSa[j]))
        end
    end
    passed = s >= 4.0 && all(values(S_ok)) && all(r + 2σ >= 0.8 && r - 2σ <= 1.2 for (r, σ) in values(ratio))
    Sc, semc, Sca = series(label, "conventional", "S_mean"), series(label, "conventional", "S_sem"), series(label, "conventional", "S_sup")
    θ = ρ * sqrt(π / D3)
    push!(condrows, (label=label, s=s, rhoRD=x, tau=τ, nseeds=Int(num(c, "nseeds")), cbar_min=num(c, "cbar_min"),
        tau_rule_extrapolated=num(c, "cbar_min") < 0.05, theta_sqrt_tau=θ * sqrt(τ),
        conv_max_abs_z_vs_annulus=s >= 4.0 ? maximum(abs.((Sc[2:end] .- Sca[2:end]) ./ semc[2:end])) : NaN,
        annulus_validated=passed, conv_S_ok=S_ok["conventional"], layers_S_ok=all(S_ok[m] for m in models if m != "conventional"),
        ratios_ok=all(r + 2σ >= 0.8 && r - 2σ <= 1.2 for (r, σ) in values(ratio)), identical=c["identical_trajectories"] == "true"))
    for m in models
        m == "conventional" && continue
        r1 = first(r for r in curves if r["label"] == label && r["model"] == m)
        sh, h = r1["shape"], num(r1, "h")
        dS, dsem = series(label, m, "dS_mean"), series(label, m, "dS_sem")
        k = argmax(abs.(dS))
        S = series(label, m, "S_mean")
        fit = passed ? fit_flags(S, times, ρ, sh, h, s) :
              (rho_c=NaN, residual=NaN, at_edge=false, ceiling=false, lambda1=NaN, residual_edge_bound=NaN, rho_c_raw=NaN)
        push!(mc, (label=label, s=s, rhoRD=x, rho=ρ, model=m, shape=sh, hR=h / R1, kept=num(r1, "kept"),
            maxdS_mc=abs(dS[k]), maxdS_sem=dsem[k], t_max=times[k], dS_ratio_ann=ratio[m][1], dS_ratio_ann_sigma=ratio[m][2],
            annulus_validated=passed, noise_floor=maximum(series(label, m, "S_sem")), fit...))
    end
end
write_csv(joinpath(OUT, "conditions_check.csv"), condrows)
write_csv(joinpath(OUT, "mc_map.csv"), mc)
for c in condrows
    @printf("%-14s τ %.2e seeds %d θ√τ %.3f  annulus validated %s (conv S %s, layers S %s, ratios %s)%s\n", c.label, c.tau, c.nseeds,
        c.theta_sqrt_tau, c.annulus_validated, c.conv_S_ok, c.layers_S_ok, c.ratios_ok, c.tau_rule_extrapolated ? "  [τ extrapolated: c̄_min < 0.05]" : "")
end
for r in mc
    @printf("%-14s %-16s max|ΔS| %.4f ± %.4f  ann ratio %.3f ± %.3f  ρc/ρ %s  residual %s %s\n", r.label, r.model, r.maxdS_mc, r.maxdS_sem,
        r.dS_ratio_ann, r.dS_ratio_ann_sigma, isnan(r.rho_c) ? "—" : @sprintf("%.3f", r.rho_c / r.rho),
        isnan(r.residual) ? "—" : @sprintf("%.2e", r.residual),
        r.at_edge || r.ceiling ? @sprintf("[flagged: edge %s ceiling %s; edge residual ≤ %.3f]", r.at_edge, r.ceiling, r.residual_edge_bound) : "")
end

# ---- thresholds (MC, fixed ρ, both spacings; best fit only at s = 4 where valid) ----
ths = NamedTuple[]
for s in sort(unique(r.s for r in mc)), sh in SHAPES, hr in sort(unique(r.hR for r in mc))
    sel = sort([r for r in mc if r.s == s && r.shape == sh && r.hR == hr], by=r -> r.rhoRD)
    length(sel) < 2 && continue
    (vd, sd) = threshold_status([r.rhoRD for r in sel], [r.maxdS_mc for r in sel], LEVEL)
    ok = [r for r in sel if !isnan(r.residual)]
    (vb, sb) = length(ok) >= 2 ? threshold_status([r.rhoRD for r in ok], [r.residual for r in ok], LEVEL) : (NaN, :no_valid_fits)
    push!(ths, (source="mc", s=s, shape=sh, hR=hr, fixed_rho_threshold=vd, fixed_rho_status=String(sd), bestfit_threshold=vb,
        bestfit_status=String(sb), n_valid_fits=length(ok)))
end

# ---- dense annulus map at s = 4, only up to the validated ρR/D ----
s4 = sort([c for c in condrows if c.s == 4.0], by=c -> c.rhoRD)
xvalid = 0.067
for c in s4
    c.annulus_validated || break
    global xvalid = max(xvalid, c.rhoRD)
end
xs = [x for x in (0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0) if x <= xvalid * (1 + 1e-9)]
hrs = [0.1, 0.25, 0.5]
jobs = [(sh, x, hr) for sh in SHAPES for x in xs for hr in hrs]
dense = Vector{NamedTuple}(undef, length(jobs))
Threads.@threads for j in eachindex(jobs)
    (sh, x, hr) = jobs[j]
    ρ = x * D3 / R1
    λc = annulus_rates(R=R1, R_out=krogh_radius(4.0), D=D3, surface=ρ)[1]
    times = collect(range(0.0, log(100) / λc, length=41))
    Sc = ann_fun(4.0; surface=ρ).(times)
    Sl = ann_fun(4.0; layer=(sh, ρ, hr * R1, supp(sh, hr * R1, 4.0))).(times)
    dense[j] = (shape=sh, rhoRD=x, rho=ρ, hR=hr, maxdS=maximum(abs.(Sl .- Sc)), fit_flags(Sl, times, ρ, sh, hr * R1, 4.0)...)
end
write_csv(joinpath(OUT, "annulus_map_s4.csv"), isempty(dense) ? [(note="no validated range",)] : dense)
for sh in SHAPES, hr in hrs
    sel = sort([d for d in dense if d.shape == sh && d.hR == hr], by=d -> d.rhoRD)
    length(sel) < 2 && continue
    (vd, sd) = threshold_status([d.rhoRD for d in sel], [d.maxdS for d in sel], LEVEL)
    ok = [d for d in sel if !isnan(d.residual)]
    (vb, sb) = length(ok) >= 2 ? threshold_status([d.rhoRD for d in ok], [d.residual for d in ok], LEVEL) : (NaN, :no_valid_fits)
    push!(ths, (source="annulus s4 (≤ ρR/D $(xvalid))", s=4.0, shape=sh, hR=hr, fixed_rho_threshold=vd, fixed_rho_status=String(sd),
        bestfit_threshold=vb, bestfit_status=String(sb), n_valid_fits=length(ok)))
end
write_csv(joinpath(OUT, "thresholds.csv"), ths)
foreach(t -> @printf("threshold %-26s s %.1f %-12s h/R %-5g fixed-ρ %s (%s)  best-fit %s (%s)\n", t.source, t.s, t.shape, t.hR,
    string(round(t.fixed_rho_threshold, sigdigits=3)), t.fixed_rho_status, string(round(t.bestfit_threshold, sigdigits=3)), t.bestfit_status), ths)

# ---- checks: annulus grid (s = 4, h/R = 0.1) and fit family at the most demanding valid fit ----
gc = NamedTuple[]
for x in (1.0, 30.0), sh in SHAPES
    ρ = x * D3 / R1
    λc = annulus_rates(R=R1, R_out=krogh_radius(4.0), D=D3, surface=ρ)[1]
    t = collect(range(0.0, log(100) / λc, length=41))
    a = ann_fun(4.0; layer=(sh, ρ, 0.1, supp(sh, 0.1, 4.0)), N=2000).(t)
    b = ann_fun(4.0; layer=(sh, ρ, 0.1, supp(sh, 0.1, 4.0)), N=4000).(t)
    push!(gc, (rhoRD=x, shape=sh, hR=0.1, max_abs_diff_S=maximum(abs.(a .- b))))
end
write_csv(joinpath(OUT, "grid_convergence.csv"), gc)
foreach(g -> @printf("annulus grid s = 4, h/R 0.1, ρR/D %-4g %-12s |S2000 − S4000| %.2e\n", g.rhoRD, g.shape, g.max_abs_diff_S), gc)
valid = [r for r in mc if !isnan(r.residual)]
famchk = Dict{String, Any}("note" => "no valid MC fit")
if !isempty(valid)
    chk = first(sort(valid, by=r -> (-r.rhoRD, -r.residual)))
    t = sort(unique(num(r, "t") for r in curves if r["label"] == chk.label))
    f2 = best_fit_conventional(series(chk.label, chk.model, "S_mean"), t, conv_family(4.0, t; N=2000); rho0=chk.rho)
    famchk = Dict("label" => chk.label, "model" => chk.model, "rho_c_N1000" => chk.rho_c_raw, "rho_c_N2000" => f2[1],
                  "residual_N1000" => chk.residual, "residual_N2000" => f2[2])
    @printf("family check %s %s: ρc %.6g vs %.6g, residual %.4e vs %.4e\n", chk.label, chk.model, chk.rho_c_raw, f2[1], chk.residual, f2[2])
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_conditions" => length(conds),
    "all_identical" => all(c -> c.identical, condrows), "annulus_valid_up_to_rhoRD_s4" => xvalid,
    "annulus_validated" => Dict(c.label => c.annulus_validated for c in condrows), "family_check" => famchk,
    "tau_rule_extrapolated" => any(c -> c.tau_rule_extrapolated, condrows),
    "files" => ["conditions_check.csv", "mc_map.csv", "thresholds.csv", "annulus_map_s4.csv", "grid_convergence.csv", "p3_7b_map.png"]))

# ---- figure: max|ΔS| vs ρR/D (rows: s = 2.5, s = 4), best-fit residual at s = 4 (row 3); columns: profiles ----
cols = Dict(0.1 => Makie.wong_colors()[1], 0.25 => Makie.wong_colors()[2], 0.5 => Makie.wong_colors()[3])
fig = Figure(size=(1400, 1100))
for (j, sh) in enumerate(SHAPES)
    for (i, s) in enumerate((2.5, 4.0))
        ax = Axis(fig[i, j], xlabel="ρR/D", ylabel="max_t |ΔS| (fixed ρ)", xscale=log10, yscale=log10,
            title="s = $(s), $(sh)" * (s == 2.5 ? " (MC only)" : " (lines: annulus ≤ ρR/D $(round(xvalid, sigdigits=2)))"))
        for hr in hrs
            p = sort([r for r in mc if r.s == s && r.shape == sh && r.hR == hr], by=r -> r.rhoRD)
            isempty(p) && continue
            errorbars!(ax, [r.rhoRD for r in p], [r.maxdS_mc for r in p], [2r.maxdS_sem for r in p], color=cols[hr])
            scatterlines!(ax, [r.rhoRD for r in p], [r.maxdS_mc for r in p], color=cols[hr], label="h/R = $(hr)")
            if s == 4.0
                d = sort([q for q in dense if q.shape == sh && q.hR == hr], by=q -> q.rhoRD)
                isempty(d) || lines!(ax, [q.rhoRD for q in d], [q.maxdS for q in d], color=cols[hr], linestyle=:dash)
            end
        end
        hlines!(ax, [LEVEL], color=:gray, linestyle=:dot)
        (i == 1 && j == 1) && axislegend(ax, position=:lt, labelsize=9)
    end
    ax3 = Axis(fig[3, j], xlabel="ρR/D", ylabel="best-fit residual (s = 4)", xscale=log10, yscale=log10,
        title="s = 4, $(sh): best fit to the annulus family (× flagged: edge-residual upper bound)")
    for hr in hrs
        p = sort([r for r in mc if r.s == 4.0 && r.shape == sh && r.hR == hr], by=r -> r.rhoRD)
        v = [r for r in p if !isnan(r.residual)]
        isempty(v) || scatterlines!(ax3, [r.rhoRD for r in v], [r.residual for r in v], color=cols[hr])
        fl = [r for r in p if !isnan(r.residual_edge_bound)]
        isempty(fl) || scatter!(ax3, [r.rhoRD for r in fl], [r.residual_edge_bound for r in fl], color=cols[hr], marker=:xcross, markersize=12)
    end
    hlines!(ax3, [LEVEL], color=:gray, linestyle=:dot)
end
Label(fig[4, :], "dotted: 0.01 (2σ at SNR 200); 2σ bars; best fit only where the annulus is validated (s = 4, per condition); c̄_min = 0.033 < 0.05: τ rule extrapolated", fontsize=11)
save(joinpath(OUT, "p3_7b_map.png"), fig)
println("saved to $OUT")
