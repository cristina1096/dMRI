# P3.6.1 / P3.6.2: analysis of the Monte Carlo timestep error inside single cylinders, against the exact radial
# continuum (recomputed here with the current solver for every row, so all rows share one reference).
# Sources (same settings: seeds 1–10, 2e4 spins in the box [−R, R]³, D = 3, readouts 0:T/10:T):
#   P3.4.6 comparison.csv + comparison_tau16.csv (step, conventional; R = 1, ρ = 0.02, T = 50);
#   P3.4.6d comparison.csv (linear, exponential); P3.6.1 scans thin, profiles, radius, rate.
# Rules (plan Task 3 amendments after the Group A review):
#   - duplicates (same shape, h, R, ρ, T, τ) must be bit-identical in S_mc (reuse check) and count once;
#   - rate scan: τ_eff = (T/10)/ceil((T/10)/τ);
#   - errors: e_R = (R_mc − R_cont)/R_cont with R = −ln S(T)/T, σ = SEM_S/(S_mc |ln S_cont|);
#             e_ΔS = (ΔS_mc − ΔS_cont)/|ΔS_cont| (> 0: |ΔS| underestimated), σ = SEM_ΔS/|ΔS_cont|;
#             conventional e_S = (S_mc − S_cont)/S_cont, σ = SEM_S/S_cont (reported, not fitted);
#   - fits: per series, on points with |e| > 2σ (sensitivity: 3σ), exponent ± SE (NaN for 2 points, none for 1);
#   - thresholds: per series, `threshold_status` on all points (ruling: the 2σ filter would remove exactly the
#     sub-threshold points that place a threshold), with the significance of the bracketing points;
#     only the 10 % ΔS level is resolvable (5 %, 2 % reported for completeness with their status).
# Run: julia --project=research/baseline research/phase3/p3_6_2_analysis.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_6_2_analysis")
res(p...) = joinpath(REPO_ROOT, "research", "results", "phase3", p...)
readrows(f) = (ls = readlines(f); h = split(ls[1], ','); [Dict(zip(h, split(l, ','))) for l in ls[2:end]])
num(r, k) = parse(Float64, r[k])
function shape_h(model)
    model == "conventional" && return ("conventional", 0.0)
    (s, h) = split(model, "_h")
    return (String(s), parse(Float64, h))
end

# ---- assemble ----
raw = NamedTuple[]
for (src, f) in (("P3.4.6", res("p3_4_6_continuum_check", "comparison.csv")), ("P3.4.6c", res("p3_4_6_continuum_check", "comparison_tau16.csv")),
                 ("P3.4.6d", res("p3_4_6d_profiles", "comparison.csv")))
    for r in readrows(f)
        (s, h) = shape_h(r["model"])
        push!(raw, (source=src, shape=s, h=h, R=1.0, rho=0.02, T=50.0, tau=num(r, "tau"), S_mc=num(r, "S50_mc"), S_sem=num(r, "S50_sem"),
            dS_mc=s == "conventional" ? NaN : num(r, "dS50_mc"), dS_sem=s == "conventional" ? NaN : num(r, "dS50_sem")))
    end
end
for scan in ("thin", "profiles", "radius", "rate")
    for r in readrows(res("p3_6_1_tau_scan_" * scan, "comparison.csv"))
        push!(raw, (source="P3.6.1 " * scan, shape=r["shape"], h=num(r, "h"), R=num(r, "R"), rho=num(r, "rho"), T=num(r, "T"),
            tau=num(r, "tau"), S_mc=num(r, "S_mc"), S_sem=num(r, "S_sem"), dS_mc=num(r, "dS_mc"), dS_sem=num(r, "dS_sem")))
    end
end
# dedupe + reuse check
key(r) = (r.shape, r.h, r.R, r.rho, r.T, r.tau)
uniq = Dict{Any, NamedTuple}()
ndup = 0
for r in raw
    k = key(r)
    if haskey(uniq, k)
        uniq[k].S_mc == r.S_mc || error("reuse check failed: $(k) differs between $(uniq[k].source) and $(r.source)")
        global ndup += 1
    else
        uniq[k] = r
    end
end
@printf("rows: %d raw, %d unique (%d bit-identical duplicates)\n", length(raw), length(uniq), ndup)

# ---- continuum and errors ----
cont = Dict{Any, Float64}()
Scont(shape, h, R, ρ, T) = get!(cont, (shape, h, R, ρ, T)) do
    shape == "conventional" ? radial_signal(R=R, D=D_REF, times=[T], surface=ρ)[1] :
                              radial_signal(R=R, D=D_REF, times=[T], layer=(shape, ρ, h))[1]
end
rate0(shape, ρ, h, R) = ρ / (h * radial_kappa(shape, h, R)) * radial_g(shape, 0.0, radial_support(shape, h, R) / h)
τeff(τ, T) = (T / 10) / ceil((T / 10) / τ - 1e-9)
rows = NamedTuple[]
for r in values(uniq)
    te = τeff(r.tau, r.T)
    Sc = Scont(r.shape, r.h, r.R, r.rho, r.T)
    Rmc, Rc = -log(r.S_mc) / r.T, -log(Sc) / r.T
    eR, σR = (Rmc - Rc) / Rc, r.S_sem / (r.S_mc * abs(log(Sc)))
    if r.shape != "conventional" && !isnan(r.dS_mc)
        dSc = Sc - Scont("conventional", 0.0, r.R, r.rho, r.T)
        edS, σdS = (r.dS_mc - dSc) / abs(dSc), r.dS_sem / abs(dSc)
    else
        dSc, edS, σdS = NaN, NaN, NaN
    end
    push!(rows, (source=r.source, shape=r.shape, h=r.h, R=r.R, rho=r.rho, T=r.T, tau=r.tau, tau_eff=te,
        rate0_tau=r.shape == "conventional" ? 0.0 : rate0(r.shape, r.rho, r.h, r.R) * te,
        step_over_h=r.shape == "conventional" ? NaN : sqrt(2 * D_REF * te) / r.h, step_over_R=sqrt(D_REF * te) / r.R,
        S_mc=r.S_mc, S_sem=r.S_sem, S_cont=Sc, e_S=(r.S_mc - Sc) / Sc, e_S_sigma=r.S_sem / Sc,
        e_R=eR, e_R_sigma=σR, dS_mc=r.dS_mc, dS_sem=r.dS_sem, dS_cont=dSc, e_dS=edS, e_dS_sigma=σdS))
end
sort!(rows, by=r -> (r.shape, r.R, r.h, r.rho, r.tau))
write_csv(joinpath(OUT, "errors.csv"), rows)

# ---- series, fits, thresholds ----
series(f) = unique([(r.shape, r.h, r.R, r.rho, r.T) for r in rows if f(r)])
fits = NamedTuple[]
function addfit!(measure, sel, xk, ek, sk, label; nsig=2)
    pts = [r for r in sel if !isnan(getproperty(r, ek)) && abs(getproperty(r, ek)) > nsig * getproperty(r, sk)]
    n = length(unique(getproperty(r, xk) for r in pts))
    if n >= 2
        (a, p, se) = fit_power([getproperty(r, xk) for r in pts], [abs(getproperty(r, ek)) for r in pts])
    else
        (a, p, se) = (NaN, NaN, NaN)
    end
    push!(fits, (measure=measure, series=label, variable=String(xk), nsig=nsig, n_points=n, a=a, p=p, se_p=se))
end
# (a) rate error vs ΔR2(0)τ: rate scan series (single model, short T)
for (s, h, R, ρ, T) in series(r -> startswith(r.source, "P3.6.1 rate"))
    sel = [r for r in rows if (r.shape, r.h, r.R, r.rho, r.T) == (s, h, R, ρ, T)]
    for nsig in (2, 3)
        addfit!("e_R", sel, :rate0_tau, :e_R, :e_R_sigma, "rate ΔR2(0) = $(round(rate0(s, ρ, h, R), digits=3))"; nsig=nsig)
    end
end
# (b), (c) difference error vs √(2Dτ)/h and vs ΔR2(0)τ, per series
for (s, h, R, ρ, T) in series(r -> r.shape != "conventional" && !isnan(r.e_dS))
    sel = [r for r in rows if (r.shape, r.h, r.R, r.rho, r.T) == (s, h, R, ρ, T)]
    for nsig in (2, 3), xk in (:step_over_h, :rate0_tau)
        addfit!("e_dS", sel, xk, :e_dS, :e_dS_sigma, "$(s) h=$(h) R=$(R)"; nsig=nsig)
    end
end
write_csv(joinpath(OUT, "fits.csv"), fits)

ths = NamedTuple[]
for (s, h, R, ρ, T) in series(r -> r.shape != "conventional" && !isnan(r.e_dS))
    sel = sort([r for r in rows if (r.shape, r.h, r.R, r.rho, r.T) == (s, h, R, ρ, T)], by=r -> r.step_over_h)
    x = [r.step_over_h for r in sel]; y = [abs(r.e_dS) for r in sel]
    for level in (0.10, 0.05, 0.02)
        (v, st) = threshold_status(x, y, level)
        # significance of the first point above the level (is the exceedance real?) and of the last point at/below it
        above = findfirst(>(level), y)
        sig_above = isnothing(above) ? NaN : y[above] / sel[above].e_dS_sigma
        push!(ths, (series="$(s) h=$(h) R=$(R)", measure="e_dS", level=level, variable="step_over_h", value=v, status=String(st),
            z_first_above=sig_above, k_rule=st == :crossed ? v^2 : NaN))
    end
end
for (s, h, R, ρ, T) in series(r -> startswith(r.source, "P3.6.1 rate"))
    sel = sort([r for r in rows if (r.shape, r.h, r.R, r.rho, r.T) == (s, h, R, ρ, T)], by=r -> r.rate0_tau)
    x = [r.rate0_tau for r in sel]; y = [abs(r.e_R) for r in sel]
    for level in (0.01, 0.001)
        (v, st) = threshold_status(x, y, level)
        above = findfirst(>(level), y)
        push!(ths, (series="rate ΔR2(0) = $(round(rate0(s, ρ, h, R), digits=3))", measure="e_R", level=level, variable="rate0_tau",
            value=v, status=String(st), z_first_above=isnothing(above) ? NaN : y[above] / sel[above].e_R_sigma, k_rule=NaN))
    end
end
write_csv(joinpath(OUT, "thresholds.csv"), ths)

for f in fits
    f.nsig == 2 && @printf("fit %-5s %-26s vs %-12s n = %d: p = %.2f ± %.2f\n", f.measure, f.series, f.variable, f.n_points, f.p, f.se_p)
end
for t in ths
    @printf("threshold %-5s %-26s level %.3f: %s %s (z first above %.1f)\n", t.measure, t.series, t.level, string(t.value), t.status, t.z_first_above)
end
conv = [r for r in rows if r.shape == "conventional"]
for r in conv
    @printf("conventional R = %.1f τ = %.3e: e_S %+.2e ± %.1e (z %+.1f)\n", r.R, r.tau, r.e_S, r.e_S_sigma, r.e_S / r.e_S_sigma)
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_raw" => length(raw), "n_unique" => length(uniq),
    "n_duplicates_bit_identical" => ndup, "files" => ["errors.csv", "fits.csv", "thresholds.csv", "p3_6_2_errors.png"]))   # NaN-bearing tables stay in CSV

# ---- figure ----
fig = Figure(size=(1100, 860))
axa = Axis(fig[1, 1], xlabel="ΔR₂(0)·τ_eff", ylabel="|e_R|", xscale=log10, yscale=log10, title="(a) rate error (rate scan, single model)")
axb = Axis(fig[1, 2], xlabel="√(2Dτ)/h", ylabel="|e_ΔS|", xscale=log10, yscale=log10, title="(b) difference error vs step/h (R = 1 unless noted)")
axc = Axis(fig[2, 1], xlabel="ΔR₂(0)·τ", ylabel="|e_ΔS|", xscale=log10, yscale=log10, title="(c) difference error vs per-step exponent")
axd = Axis(fig[2, 2], xlabel="√(Dτ)/R", ylabel="|e_S| conventional", xscale=log10, yscale=log10, title="(d) conventional S error (reported, not fitted)")
plotpts!(ax, xs, ys, σs, col, lab) = (ok = ys .> 0; scatterlines!(ax, xs[ok], ys[ok], color=col, label=lab);
                                     sig = ok .& (ys .> 2 .* σs); scatter!(ax, xs[sig], ys[sig], color=col, marker=:star5, markersize=12))
cols = Makie.wong_colors()
cols20 = Makie.to_colormap(:tab20)                       # 11 difference-error series need distinct colours
for (k, (s, h, R, ρ, T)) in enumerate(series(r -> startswith(r.source, "P3.6.1 rate")))
    sel = sort([r for r in rows if (r.shape, r.h, r.R, r.rho, r.T) == (s, h, R, ρ, T)], by=r -> r.rate0_tau)
    plotpts!(axa, [r.rate0_tau for r in sel], [abs(r.e_R) for r in sel], [r.e_R_sigma for r in sel], cols[k], "ΔR₂(0) = $(round(rate0(s, ρ, h, R), digits=2))")
end
vlines!(axa, [0.005], color=:black, linestyle=:dash, label="c = 0.005")
axislegend(axa, position=:lt, labelsize=9)
for (k, (s, h, R, ρ, T)) in enumerate(series(r -> r.shape != "conventional" && !isnan(r.e_dS)))
    sel = sort([r for r in rows if (r.shape, r.h, r.R, r.rho, r.T) == (s, h, R, ρ, T)], by=r -> r.tau)
    lab = "$(s) h=$(h)" * (R == 1.0 ? "" : " R=$(R)")
    plotpts!(axb, [r.step_over_h for r in sel], [abs(r.e_dS) for r in sel], [r.e_dS_sigma for r in sel], cols20[mod1(2k - 1, 20)], lab)
    plotpts!(axc, [r.rate0_tau for r in sel], [abs(r.e_dS) for r in sel], [r.e_dS_sigma for r in sel], cols20[mod1(2k - 1, 20)], lab)
end
hlines!(axb, [0.1], color=:black, linestyle=:dash)
axislegend(axb, position=:lt, labelsize=8)
for (k, R) in enumerate(sort(unique(r.R for r in conv)))
    sel = sort([r for r in conv if r.R == R], by=r -> r.tau)
    plotpts!(axd, [r.step_over_R for r in sel], abs.([r.e_S for r in sel]), [r.e_S_sigma for r in sel], cols[k], "R = $(R)")
end
axislegend(axd, position=:lt, labelsize=9)
Label(fig[3, :], "stars: |e| > 2σ (used in fits); open points: not significant", fontsize=11)
save(joinpath(OUT, "p3_6_2_errors.png"), fig)
println("saved to $OUT")
