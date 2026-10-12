# Review C1 follow-up — P2.5 and P2.5.6/7 recomputed with the exact 1D continuum benchmark (continuum.jl), which has no
# Monte Carlo noise and no timestep bias. The Monte Carlo tables at the adopted τ underestimate |ΔS| and flip its sign
# for thin layers (review check C1); these continuum values replace them, with the C1 Monte Carlo runs at small τ as
# confirmation.
#
# (A) Validation: continuum ΔS(90) vs Monte Carlo at τ = 6.25e-4 (p2_5_tau_check) and vs P2.5 rows with h ≥ 0.5.
# (B) P2.5 sweeps (h, w, D, ρ; readouts 0:10:200): peak ΔS = S_layer − S_conv at fixed ρ, profile pairs, centroid.
# (C) P2.5.6/7: apparent T2 (TE = 50 and 5:5:200; R2_bulk = 1/2000) for same h and same centroid; ΔT2/T2 vs the
#     two-point T2 noise at SNR 100/200.
# (D) Wider map: ρw/D from 0.003 to 30 (w = 2, D = 3), h/w ∈ {0.01, 0.05, 0.1, 0.25, 0.5}, three shapes, no bulk R2,
#     41 readouts from 0 until the conventional signal reaches ≈ 1 % (log(100)/its slowest rate): max |ΔS| at fixed ρ,
#     and the minimax residual after fitting the conventional model's ρ freely (coarse scan over e^−4…e^8 × ρ, then
#     golden section; edge minima flagged). Both are compared with the single-point thresholds 2/SNR; the max over 41
#     correlated readouts makes this an indicator, not a calibrated detection test.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_5_continuum.jl   (figures: p2_5_continuum_figs.jl)

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "sensitivity.jl"))
include(joinpath(@__DIR__, "continuum.jl"))

const OUT = phase2_outdir("p2_5_continuum")
const R = joinpath(REPO_ROOT, "research", "results", "phase2")
const T25 = collect(0.0:10.0:200.0)
csvrows(f) = (L = readlines(f); hdr = split(L[1], ","); [Dict(zip(hdr, split(l, ","))) for l in L[2:end]])
num(d, k) = parse(Float64, d[k])

conv_cache = Dict{Any, Vector{Float64}}()
Sconv(w, D, ρ, times; R2=0.0) = get!(() -> continuum_signal(w=w, D=D, times=times, surface=ρ, R2_bulk=R2), conv_cache, (w, D, ρ, times, R2))
Slay(shape, ρ, h, w, D, times; R2=0.0) = continuum_signal(w=w, D=D, times=times, layer=(shape, ρ, h), R2_bulk=R2)

# (A) validation
valrows = NamedTuple[]
for r in csvrows(joinpath(R, "p2_5_tau_check", "tau_check.csv"))
    num(r, "tau_ms") == 6.25e-4 || continue
    c = Slay(r["shape"], 0.01, num(r, "h_um"), 2.0, 3.0, [90.0])[1] - Sconv(2.0, 3.0, 0.01, [90.0])[1]
    push!(valrows, (source="tau_check τ=6.25e-4", shape=r["shape"], h_um=num(r, "h_um"), mc=num(r, "dS_90"), mc_sem=num(r, "dS_90_sem"), continuum=c,
        z=(num(r, "dS_90") - c) / num(r, "dS_90_sem")))
end
for r in csvrows(joinpath(R, "p2_5_sensitivity_h", "rows.csv"))
    num(r, "h_um") >= 0.5 || continue
    tp = num(r, "dS_peak_t_ms")
    c = Slay(r["shape"], 0.01, num(r, "h_um"), 2.0, 3.0, [tp])[1] - Sconv(2.0, 3.0, 0.01, [tp])[1]
    push!(valrows, (source="P2.5 h sweep τ=1e-2", shape=r["shape"], h_um=num(r, "h_um"), mc=num(r, "dS_peak"), mc_sem=num(r, "dS_peak_sem"), continuum=c,
        z=(num(r, "dS_peak") - c) / num(r, "dS_peak_sem")))
end
for v in valrows
    @printf("(A) %-20s %-11s h=%.2f: MC %+.2e ± %.0e, continuum %+.2e, z = %+.1f\n", v.source, v.shape, v.h_um, v.mc, v.mc_sem, v.continuum, v.z)
end
write_csv(joinpath(OUT, "validation.csv"), valrows)

# (B) P2.5 sweeps
p25 = NamedTuple[]
for s in ("h", "w", "D", "rho"), c in sweep_configs(s)
    d = Slay(c.shape, c.rho, c.h, c.w, c.D, T25) .- Sconv(c.w, c.D, c.rho, T25)
    k = argmax(abs.(d))
    level = abs(d[k]) >= 0.02 ? "SNR 100" : abs(d[k]) >= 0.01 ? "SNR 200" : "below SNR 200"
    push!(p25, (sweep=s, shape=c.shape, h_um=c.h, w_um=c.w, D=c.D, rho=c.rho, centroid_um=centroid(c.shape, c.h; cutoff=c.w),
        dS_peak=d[k], dS_peak_t_ms=T25[k], level=level))
end
write_csv(joinpath(OUT, "p2_5_continuum.csv"), p25)
uniq25 = Dict((r.shape, r.h_um, r.w_um, r.D, r.rho) => r for r in p25)
pairs25 = NamedTuple[]
for (k, r) in uniq25, (k2, r2) in uniq25
    (k[2:end] == k2[2:end] && r.shape < r2.shape) || continue
    push!(pairs25, (pair=r.shape * "−" * r2.shape, h_um=k[2], w_um=k[3], D=k[4], rho=k[5], diff_peaks=r.dS_peak - r2.dS_peak))
end
write_csv(joinpath(OUT, "p2_5_continuum_pairs.csv"), pairs25)
mx = argmax(r -> abs(r.dS_peak), collect(values(uniq25)))
@printf("(B) P2.5: %d distinct configurations, all negative: %s; largest |ΔS| = %.2e (%s, h=%.2f w=%.1f D=%.0f ρ=%.3f)\n",
    length(uniq25), all(r -> r.dS_peak < 0, values(uniq25)), abs(mx.dS_peak), mx.shape, mx.h_um, mx.w_um, mx.D, mx.rho)

# (C) apparent T2, same h and same centroid
t2rows = NamedTuple[]
for match in (:h, :centroid), s in (match == :h ? ("te", "w", "h", "D", "rho") : ("te", "w", "c", "D", "rho"))
    TEs = s == "te" ? collect(5.0:5.0:200.0) : [20.0, 50.0, 100.0]
    times = vcat(0.0, TEs)
    for c in t2_configs(s; match=match)
        Sl = Slay(c.shape, c.rho, c.h, c.w, c.D, times; R2=1 / 2000)
        Sc = Sconv(c.w, c.D, c.rho, times; R2=1 / 2000)
        for (j, TE) in enumerate(TEs)
            Tl, Tc = apparent_t2(Sl[1], Sl[j + 1], TE), apparent_t2(Sc[1], Sc[j + 1], TE)
            rel = (Tl - Tc) / Tc
            r = Sc[j + 1] / Sc[1]
            level = abs(rel) >= 2t2_noise_rel(r, 100) ? "SNR 100" : abs(rel) >= 2t2_noise_rel(r, 200) ? "SNR 200" : "below SNR 200"
            push!(t2rows, (match=String(match), sweep=s, shape=c.shape, h_um=c.h, w_um=c.w, D=c.D, rho=c.rho,
                centroid_um=centroid(c.shape, c.h; cutoff=c.w), TE_ms=TE, T2_layer=Tl, T2_conv=Tc, rel=rel,
                noise_rel_snr200=t2_noise_rel(r, 200), level=level))
        end
    end
end
write_csv(joinpath(OUT, "t2_continuum.csv"), t2rows)
for match in ("h", "centroid")
    sel = [r for r in t2rows if r.match == match && r.TE_ms == 50.0 && r.sweep != "te"]
    m = argmax(r -> abs(r.rel), sel)
    @printf("(C) T2 (%s), TE = 50: all ΔT2/T2 < 0: %s; largest |ΔT2/T2| = %.2e (%s, %s sweep); observable at SNR 200: %d\n",
        match, all(r -> r.rel < 0, sel), abs(m.rel), m.shape, m.sweep, count(r -> r.level != "below SNR 200", sel))
end

# profile spread at TE = 50: computed from t2_continuum.csv in p2_5_continuum_figs.jl (distinct configurations)

# (D) wider map, fixed ρ and best-fit conventional (single condition)
function best_fit_residual(Sl, w, D, ρ0, times)
    # coarse scan of log ρ_fit over ρ0·e^(−4) … ρ0·e^(8), then golden-section refinement around the best scan point;
    # `at_edge` flags a minimum on the scan boundary (would make the residual only an upper bound)
    f(lr) = maximum(abs.(Sl .- continuum_signal(w=w, D=D, times=times, surface=exp(lr), N=1500)))
    grid = log(ρ0) .+ (-4.0:0.25:8.0)
    vals = f.(grid)
    k = argmin(vals)
    at_edge = k == 1 || k == length(grid)
    a, b = grid[max(k - 1, 1)], grid[min(k + 1, length(grid))]
    g = (sqrt(5) - 1) / 2
    x1, x2 = b - g * (b - a), a + g * (b - a)
    f1, f2 = f(x1), f(x2)
    for _ in 1:30
        if f1 < f2
            b, x2, f2 = x2, x1, f1; x1 = b - g * (b - a); f1 = f(x1)
        else
            a, x1, f1 = x1, x2, f2; x2 = a + g * (b - a); f2 = f(x2)
        end
    end
    best = min(f1, f2, vals[k])
    xbest = best == f1 ? x1 : best == f2 ? x2 : grid[k]
    return (resid=best, rho_fit=exp(xbest), at_edge=at_edge)
end
mapr = NamedTuple[]
w, D = 2.0, 3.0
for ρwD in 10 .^ range(log10(0.003), log10(30), length=9), hw in (0.01, 0.05, 0.1, 0.25, 0.5), shape in ("step", "linear", "exponential")
    shape == "step" && hw > 0.5 && continue
    ρ = ρwD * D / w
    h = hw * w
    tmax = log(100) / continuum_rates(w=w, D=D, surface=ρ, N=1500)[1]   # window: until the conventional signal ≈ 1 %
    times = collect(range(0.0, tmax, length=41))
    Sl = continuum_signal(w=w, D=D, times=times, layer=(shape, ρ, h), N=1500)
    Sc = continuum_signal(w=w, D=D, times=times, surface=ρ, N=1500)
    dS = maximum(abs.(Sl .- Sc))
    bf = best_fit_residual(Sl, w, D, ρ, times)
    push!(mapr, (rho_w_over_D=ρwD, h_over_w=hw, shape=shape, rho=ρ, tmax_ms=tmax, max_abs_dS_fixed_rho=dS, best_fit_residual=bf.resid,
        rho_fit_over_rho=bf.rho_fit / ρ, fit_at_edge=bf.at_edge))
    @printf("(D) ρw/D=%.3g h/w=%.2f %-11s: max|ΔS| fixed ρ %.2e, best-fit residual %.2e (ρ_fit/ρ = %.3f%s)\n", ρwD, hw, shape, dS, bf.resid, bf.rho_fit / ρ, bf.at_edge ? ", AT EDGE" : "")
    flush(stdout)
end
write_csv(joinpath(OUT, "map.csv"), mapr)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "solver" => "continuum.jl, N = 3000 (map N = 1500)",
    "validation" => valrows, "p25_distinct" => length(uniq25), "p25_all_negative" => all(r -> r.dS_peak < 0, values(uniq25)),
    "p25_max_abs_dS" => abs(mx.dS_peak), "map_fit_at_edge" => count(r -> r.fit_at_edge, mapr)))

println("saved to $OUT (figures: p2_5_continuum_figs.jl)")
