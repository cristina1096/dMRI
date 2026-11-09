# P3.7.3 inside (plan 2026-11-08-comparison-inside-cylinders.md, Task 5): dimensionless map over ρR/D and h/R inside a
# single cylinder, Monte Carlo primary (p3_7_run.jl TASK=map: R = 1, D = 3, no bulk, ρR/D ∈ {0.1, 1, 10} × h/R ∈ {0.1,
# 0.25, 0.5} × step/linear/exponential, 41 readouts to t_end = ln 100/λ₁,conv), the exact radial continuum as support.
# Quantities per layer model:
#   fixed ρ:   max_t |ΔS(t)|, ΔS = S_layer − S_conv (paired MC; SEM at the arg-max readout; continuum alongside, z);
#   best fit:  ρ_c = argmin max_t |S_layer(t) − S_conv(t; ρ_c)| over the conventional continuum family (minimax,
#              `best_fit_conventional`); residual against the noise floor max_t SEM(S_layer);
#   flags:     at_edge (fit leaves the searched range or sits on a plateau reaching it; Group A review I1) and ceiling
#              (the layer's lowest continuum rate λ₁ exceeds the fastest conventional decay D·j₀₁²/R²). For flagged
#              fits ρ_c and the residual are not reported (NaN): the conventional model cannot represent the curve, and
#              the residual of an undetermined fit is not a distance;
#   thresholds: first ρR/D at which max|ΔS| and the best-fit residual reach 0.01 (`threshold_status`, per profile and h/R).
# Also reported (Group A review M2): θ√τ of the conventional run (θ = ρ√(π/D)) and its MC − continuum z.
# Support: continuum dense map ρR/D ∈ {0.003 … 30} × h/R ∈ {0.01 … 1} (same quantities, both from the continuum);
# grid convergence of the layer curve at h/R = 0.01, ρR/D = 1 (N = 2000 vs 4000). The conventional fit family uses
# N = 1000 (no layer to resolve; second-order Robin face), checked against N = 2000 at one MC point.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_7_map_analysis.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
include(joinpath(@__DIR__, "p3_7_shared.jl"))
using CairoMakie
const SRC = joinpath(REPO_ROOT, "research", "results", "phase3", "p3_7_map")
const OUT = phase3_outdir("p3_7_map_analysis")
const R1, D3 = 1.0, 3.0
const SHAPES = ("step", "linear", "exponential")
const NFAM = 1000
const LEVEL = 0.01
const CEIL = conventional_ceiling(R1, D3)

"Conventional continuum family on `times` as a function of ρ_c (N = NFAM)."
conv_family(times; N=NFAM) = ρc -> (f = signal_fun(R=R1, D=D3, surface=ρc, N=N); f.(times))
lambda1_layer(s, ρ, h) = radial_rates(R=R1, D=D3, layer=(s, ρ, h))[1]
function fit_flags(S_layer, times, ρ, s, h; N=NFAM)
    (ρc, res, edge) = best_fit_conventional(S_layer, times, conv_family(times; N=N); rho0=ρ)
    λ1 = lambda1_layer(s, ρ, h)
    ceil = λ1 > CEIL
    ok = !edge && !ceil
    return (rho_c=ok ? ρc : NaN, residual=ok ? res : NaN, at_edge=edge, ceiling=ceil, lambda1=λ1, rho_c_raw=ρc, residual_raw=res)
end

# ---- Monte Carlo ----
conds = Dict(c["label"] => c for c in readrows(joinpath(SRC, "conditions.csv")))
curves = readrows(joinpath(SRC, "curves.csv"))
mc = NamedTuple[]
convrows = NamedTuple[]
for (label, c) in sort(collect(conds), by=first)
    ρ, τ = num(c, "rho"), num(c, "tau")
    x = ρ * R1 / D3
    sel = [r for r in curves if r["label"] == label]
    times = sort(unique(num(r, "t") for r in sel))
    get(m, k) = [num(r, k) for r in sort([r for r in sel if r["model"] == m], by=r -> num(r, "t"))]
    zc = (get("conventional", "S_mean") .- get("conventional", "S_cont")) ./ max.(get("conventional", "S_sem"), eps())
    θ = ρ * sqrt(π / D3)
    push!(convrows, (label=label, rhoRD=x, tau=τ, theta_sqrt_tau=θ * sqrt(τ), max_abs_z_S=maximum(abs.(zc[2:end])),
        identical=c["identical_trajectories"] == "true"))
    for m in unique(r["model"] for r in sel)
        m == "conventional" && continue
        r1 = first(r for r in sel if r["model"] == m)
        s, h = r1["shape"], num(r1, "h")
        dS, dSsem, dSc = get(m, "dS_mean"), get(m, "dS_sem"), get(m, "dS_cont")
        k = argmax(abs.(dS))
        kc = argmax(abs.(dSc))
        Sl = get(m, "S_mean")
        fit = fit_flags(Sl, times, ρ, s, h)
        push!(mc, (label=label, rhoRD=x, rho=ρ, model=m, shape=s, hR=h / R1, tau=τ,
            maxdS_mc=abs(dS[k]), maxdS_sem=dSsem[k], t_max_mc=times[k], maxdS_cont=abs(dSc[kc]), t_max_cont=times[kc],
            z_dS=(dS[k] - dSc[k]) / dSsem[k], ratio=dS[k] / dSc[k],
            noise_floor=maximum(get(m, "S_sem")), fit...))
    end
end
write_csv(joinpath(OUT, "mc_map.csv"), mc)
write_csv(joinpath(OUT, "conventional_check.csv"), convrows)
for r in mc
    @printf("ρR/D %-5g %-18s max|ΔS| %.4f ± %.4f (cont %.4f, z %+5.1f)  ρc/ρ %6.3f  residual %.2e (floor %.1e) edge %s ceiling %s\n",
        r.rhoRD, r.model, r.maxdS_mc, r.maxdS_sem, r.maxdS_cont, r.z_dS, r.rho_c / r.rho, r.residual, r.noise_floor, r.at_edge, r.ceiling)
end
for r in convrows
    @printf("conventional ρR/D %-5g τ %.2e θ√τ %.3f max|z_S| %.1f\n", r.rhoRD, r.tau, r.theta_sqrt_tau, r.max_abs_z_S)
end
# family resolution check at the most demanding VALID fit (final review I4: a flagged point sits at the grid end for both
# N and proves nothing): the non-flagged MC fit with the largest ρR/D, then the largest residual
chk = first(sort([r for r in mc if !r.at_edge && !r.ceiling], by=r -> (-r.rhoRD, -r.residual)))
chk_times = sort(unique(num(r, "t") for r in curves if r["label"] == chk.label))
chk_S = [num(r, "S_mean") for r in sort([r for r in curves if r["label"] == chk.label && r["model"] == chk.model], by=r -> num(r, "t"))]
chk2000 = best_fit_conventional(chk_S, chk_times, conv_family(chk_times; N=2000); rho0=chk.rho)
@printf("family check (%s, %s): ρc N=1000 %.6g, N=2000 %.6g; residual %.3e vs %.3e\n", chk.label, chk.model, chk.rho_c_raw,
    chk2000[1], chk.residual_raw, chk2000[2])

# ---- continuum dense map ----
xs = [0.003, 0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]
hrs = [0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1.0]
jobs = [(s, x, hr) for s in SHAPES for x in xs for hr in hrs]
dense = Vector{NamedTuple}(undef, length(jobs))
Threads.@threads for j in eachindex(jobs)
    (s, x, hr) = jobs[j]
    ρ = x * D3 / R1
    λc = radial_rates(R=R1, D=D3, surface=ρ)[1]
    times = collect(range(0.0, log(100) / λc, length=41))
    Sc = signal_fun(R=R1, D=D3, surface=ρ).(times)
    Sl = signal_fun(R=R1, D=D3, layer=(s, ρ, hr * R1)).(times)
    fit = fit_flags(Sl, times, ρ, s, hr * R1)
    dense[j] = (shape=s, rhoRD=x, rho=ρ, hR=hr, maxdS=maximum(abs.(Sl .- Sc)), fit...)
end
write_csv(joinpath(OUT, "continuum_map.csv"), dense)
# grid convergence of the layer curve at the thinnest layer
gc = NamedTuple[]
for s in SHAPES
    ρ = 1.0 * D3 / R1
    λc = radial_rates(R=R1, D=D3, surface=ρ)[1]
    times = collect(range(0.0, log(100) / λc, length=41))
    a = signal_fun(R=R1, D=D3, layer=(s, ρ, 0.01), N=2000).(times)
    b = signal_fun(R=R1, D=D3, layer=(s, ρ, 0.01), N=4000).(times)
    Sc = signal_fun(R=R1, D=D3, surface=ρ).(times)
    push!(gc, (shape=s, hR=0.01, rhoRD=1.0, max_abs_diff_S=maximum(abs.(a .- b)), maxdS_N2000=maximum(abs.(a .- Sc)),
        maxdS_N4000=maximum(abs.(b .- Sc))))
end
write_csv(joinpath(OUT, "grid_convergence.csv"), gc)
foreach(g -> @printf("grid h/R = 0.01 %-12s |S2000 − S4000| %.2e, max|ΔS| %.5f vs %.5f\n", g.shape, g.max_abs_diff_S, g.maxdS_N2000, g.maxdS_N4000), gc)

# ---- thresholds ----
ths = NamedTuple[]
for (src, rows) in (("continuum", dense), ("mc", mc)), s in SHAPES, hr in sort(unique(r.hR for r in rows))
    sel = sort([r for r in rows if r.shape == s && r.hR == hr], by=r -> r.rhoRD)
    length(sel) < 2 && continue
    xsel = [r.rhoRD for r in sel]
    yd = src == "mc" ? [r.maxdS_mc for r in sel] : [r.maxdS for r in sel]
    (vd, sd) = threshold_status(xsel, yd, LEVEL)
    okr = [!isnan(r.residual) for r in sel]
    (vr, sr) = count(okr) >= 2 ? threshold_status(xsel[okr], [r.residual for r in sel][okr], LEVEL) : (NaN, :no_valid_fits)
    push!(ths, (source=src, shape=s, hR=hr, level=LEVEL, fixed_rho_threshold=vd, fixed_rho_status=String(sd),
        bestfit_threshold=vr, bestfit_status=String(sr), n_flagged=count(.!okr)))
end
write_csv(joinpath(OUT, "thresholds.csv"), ths)
foreach(t -> @printf("threshold %-9s %-12s h/R %-5g fixed-ρ %s (%s)  best-fit %s (%s)  flagged %d\n", t.source, t.shape, t.hR,
    string(round(t.fixed_rho_threshold, sigdigits=3)), t.fixed_rho_status, string(round(t.bestfit_threshold, sigdigits=3)),
    t.bestfit_status, t.n_flagged), ths)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_mc" => length(mc), "ceiling" => CEIL,
    "all_identical" => all(r -> r.identical, convrows), "max_abs_z_dS" => maximum(abs(r.z_dS) for r in mc),
    "family_check" => Dict("rho_c_N1000" => chk.rho_c_raw, "rho_c_N2000" => chk2000[1]),
    "files" => ["mc_map.csv", "conventional_check.csv", "continuum_map.csv", "grid_convergence.csv", "thresholds.csv", "p3_7_map.png"]))

# ---- figure: rows = fixed-ρ max|ΔS|, best-fit residual, ρ_c/ρ; columns = profiles ----
fig = Figure(size=(1400, 1150))
lx = log10.(xs); lh = log10.(hrs)
cell(rows, s, f) = [f(only(r for r in rows if r.shape == s && r.rhoRD == x && r.hR == hr)) for x in xs, hr in hrs]
for (j, s) in enumerate(SHAPES)
    for (i, (name, f, fmc, cr)) in enumerate((
            ("log₁₀ max|ΔS| (fixed ρ)", r -> log10(max(r.maxdS, 1e-8)), r -> r.maxdS_mc, (-6, 0)),
            ("log₁₀ best-fit residual", r -> isnan(r.residual) ? NaN : log10(max(r.residual, 1e-8)), r -> r.residual, (-6, 0)),
            ("log₁₀ ρ_c/ρ", r -> isnan(r.rho_c) ? NaN : log10(r.rho_c / r.rho), r -> r.rho_c / r.rho, (-1.5, 1.5))))
        ax = Axis(fig[i, j], xlabel="log₁₀ ρR/D", ylabel="log₁₀ h/R", title="$(s): $(name)")
        Z = cell(dense, s, f)
        hm = heatmap!(ax, lx, lh, Z, colormap=i == 3 ? :balance : :viridis, colorrange=cr, nan_color=(:gray, 0.35))
        i < 3 && contour!(ax, lx, lh, Z, levels=[log10(LEVEL)], color=:white, linewidth=2)
        fl = [r for r in dense if r.shape == s && (r.at_edge || r.ceiling)]
        scatter!(ax, log10.([r.rhoRD for r in fl]), log10.([r.hR for r in fl]), marker=:xcross, color=:red, markersize=9)
        for r in (m for m in mc if m.shape == s)
            v = fmc(r)
            txt = isnan(v) ? (r.ceiling ? "ceil" : "edge") : (i == 3 ? @sprintf("%.2f", v) : @sprintf("%.0e", v))
            scatter!(ax, [log10(r.rhoRD)], [log10(r.hR)], color=:white, strokecolor=:black, strokewidth=1, markersize=11)
            text!(ax, log10(r.rhoRD), log10(r.hR), text=txt, fontsize=9, offset=(6, 4), color=:white)
        end
        j == 3 && Colorbar(fig[i, 4], hm, label=name)
    end
end
Label(fig[4, :], "heatmap: continuum; white contour: 0.01; grey cells / red ×: best fit flagged (edge or λ₁ > D·j₀₁²/R², ρ_c and residual not reported); circles: Monte Carlo (value printed)", fontsize=11)
save(joinpath(OUT, "p3_7_map.png"), fig)
println("saved to $OUT")
