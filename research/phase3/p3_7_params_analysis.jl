# P3.7.1 (plan 2026-11-08-comparison-inside-cylinders.md, Task 4): parameter sweeps at R = 1 — one parameter at a time
# from the reference ρ = 1/140, D = 3, c̄ = 0.1 (same-centroid matching; R2_bulk = 1/2000). Monte Carlo primary
# (p3_7_run.jl TASK=params; the reference point is fig9a's R1.0_c0.1), the exact radial continuum as support.
# Sweeps: D ∈ {0.5, 1, 2, 3}; ρ ∈ {1/280, 1/140, 1/70, 1/35, 0.1}; c̄ ∈ {0.05, 0.1, 0.25}.
# Quantities at TE* = continuum conventional T₂ of each condition: as in p3_7_fig9a_analysis.jl (T₂_app, paired ΔT₂/T₂
# with σ, 2σ two-point noise at SNR 100/200, z values) and ΔS_mc/ΔS_cont (Group A review I2).
# Support: continuum dense curves over each sweep (D 0.25–4, ρ 1/560–0.2, c̄ 0.01–0.3, log-spaced); the continuum-only
# points ρ = 1/560 and c̄ ∈ {0.01, 0.025} (below the MC range; c̄ < 0.05 would need τ < 3.4e-4) are marked separately.
# Run: julia --project=research/baseline research/phase3/p3_7_params_analysis.jl  (works on partial MC results)

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
include(joinpath(@__DIR__, "p3_7_shared.jl"))
using CairoMakie
const RES = joinpath(REPO_ROOT, "research", "results", "phase3")
const OUT = phase3_outdir("p3_7_params_analysis")
const RHO, D3, BULK, CB, R1 = 1 / 140, 3.0, 1 / 2000, 0.1, 1.0
const SHAPES = ("step", "linear", "exponential")

# ---- Monte Carlo at TE* (params + the fig9a reference) ----
ref = [r for r in at_te_rows(joinpath(RES, "p3_7_fig9a")) if r.label == "R1.0_c0.1"]
rows = vcat(ref, at_te_rows(joinpath(RES, "p3_7_params")))
cond = Dict(c["label"] => c for c in vcat(readrows(joinpath(RES, "p3_7_fig9a", "conditions.csv")),
                                          readrows(joinpath(RES, "p3_7_params", "conditions.csv"))))
# sweep membership and the swept value of each MC condition
function sweeps_of(label)
    c = cond[label]
    label == "R1.0_c0.1" && return [("D", D3), ("rho", RHO), ("cbar", CB)]
    startswith(label, "D") && return [("D", num(c, "D"))]
    startswith(label, "rho") && return [("rho", num(c, "rho"))]
    startswith(label, "c") && return [("cbar", parse(Float64, label[2:end]))]
    error("unknown label $label")
end
mc = NamedTuple[]
for r in rows, (sw, x) in sweeps_of(r.label)
    push!(mc, merge((sweep=sw, x=x), r))
end
write_csv(joinpath(OUT, "at_te.csv"), mc)
bias = dS_ratio_rows(rows)
write_csv(joinpath(OUT, "dS_ratio.csv"), bias)
for r in rows
    @printf("%-12s %-20s T₂ %7.2f ± %.2f (cont %7.2f, z_S %+5.1f)  ΔT₂/T₂ %+8.3f%% ± %.3f%% (cont %+8.3f%%)  ΔS ratio %s\n",
        r.label, r.model, r.T2_mc, r.T2_sigma, r.T2_cont, r.z_S, 100r.dT2_rel, 100r.dT2_rel_sigma, 100r.dT2_rel_cont,
        r.model == "conventional" ? "—" : @sprintf("%.3f ± %.3f", r.dS_mc / r.dS_cont, r.dS_sem / abs(r.dS_cont)))
end

# ---- continuum dense curves and continuum-only points ----
function cont_point(; D=D3, ρ=RHO, cb=CB)
    Sconv = signal_fun(R=R1, D=D, surface=ρ, bulk=BULK)
    te = te_star(Sconv; t_hi=50 / (BULK + 2ρ / R1))
    T2c = apparent_t2(Sconv(te), te)
    out = NamedTuple[]
    for s in SHAPES
        h = h_for_cbar(s, cb, R1)
        radial_support(s, h, R1) <= R1 || continue
        T2 = apparent_t2(signal_fun(R=R1, D=D, layer=(s, ρ, h), bulk=BULK)(te), te)
        push!(out, (shape=s, h=h, te=te, T2_conv=T2c, T2=T2, dT2_rel=(T2 - T2c) / T2c, noise2_snr200=2 * t2_noise_rel(exp(-1), 200)))
    end
    return out
end
logrange(a, b, n) = exp.(range(log(a), log(b), length=n))
dense = NamedTuple[]
for (sw, xs, kw) in (("D", logrange(0.25, 4.0, 25), x -> (D=x,)), ("rho", logrange(1 / 560, 0.2, 25), x -> (ρ=x,)),
                     ("cbar", logrange(0.01, 0.3, 25), x -> (cb=x,)))
    for x in xs, p in cont_point(; kw(x)...)
        push!(dense, merge((sweep=sw, x=x), p))
    end
end
write_csv(joinpath(OUT, "continuum_dense.csv"), dense)
conly = NamedTuple[]
for (sw, x, kw) in (("rho", 1 / 560, (ρ=1 / 560,)), ("cbar", 0.01, (cb=0.01,)), ("cbar", 0.025, (cb=0.025,))), p in cont_point(; kw...)
    push!(conly, merge((sweep=sw, x=x, continuum_only=true), p))
end
write_csv(joinpath(OUT, "continuum_only.csv"), conly)
for p in conly
    @printf("continuum only %-4s = %.4g %-12s ΔT₂/T₂ %+8.3f%%\n", p.sweep, p.x, p.shape, 100p.dT2_rel)
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_mc_conditions" => length(unique(r.label for r in rows)),
    "all_identical" => all(r -> r.identical, rows), "max_abs_z_S" => maximum(abs(r.z_S) for r in rows),
    "max_abs_dT2_rel_mc" => maximum((abs(r.dT2_rel) for r in rows if r.model != "conventional"); init=0.0),
    "noise2_snr200" => 2 * t2_noise_rel(exp(-1), 200), "noise2_snr100" => 2 * t2_noise_rel(exp(-1), 100),
    "files" => ["at_te.csv", "dS_ratio.csv", "continuum_dense.csv", "continuum_only.csv", "p3_7_params.png"]))

# ---- figure: |ΔT₂/T₂| (top) and ΔS_mc/ΔS_cont (bottom) per sweep ----
cols = Dict("step" => Makie.wong_colors()[1], "linear" => Makie.wong_colors()[2], "exponential" => Makie.wong_colors()[3])
xlab = Dict("D" => "D (µm²/ms)", "rho" => "ρ (µm/ms)", "cbar" => "c̄ (µm)")
n200, n100 = 100 * 2 * t2_noise_rel(exp(-1), 200), 100 * 2 * t2_noise_rel(exp(-1), 100)
xt = Dict("D" => [0.25, 0.5, 1, 2, 4], "rho" => [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2], "cbar" => [0.01, 0.02, 0.05, 0.1, 0.2])
xr = Dict("D" => (0.2, 5.0), "rho" => (1.5e-3, 0.25), "cbar" => (8e-3, 0.4))
yt = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1, 3, 10]
fig = Figure(size=(1400, 820))
for (j, sw) in enumerate(("D", "rho", "cbar"))
    ax = Axis(fig[1, j], xlabel=xlab[sw], ylabel="|ΔT₂/T₂| vs conventional (%)", xscale=log10, yscale=log10,
        title="($( 'a' + j - 1)) sweep in $(xlab[sw]) (filled ΔT₂ < 0, open > 0)",
        xticks=(xt[sw], string.(xt[sw])), yticks=(yt, string.(yt)))
    xlims!(ax, xr[sw]...)
    for s in SHAPES
        sel = [d for d in dense if d.sweep == sw && d.shape == s]
        lines!(ax, [d.x for d in sel], [max(100abs(d.dT2_rel), 1e-6) for d in sel], color=cols[s], label="$(s) continuum")
        co = [p for p in conly if p.sweep == sw && p.shape == s]
        isempty(co) || scatter!(ax, [p.x for p in co], [100abs(p.dT2_rel) for p in co], color=:white, strokecolor=cols[s],
            strokewidth=1.5, marker=:diamond, markersize=10)
        for r in (m for m in mc if m.sweep == sw && m.shape == s)
            y, e = 100abs(r.dT2_rel), 200r.dT2_rel_sigma
            rangebars!(ax, [r.x], [max(y - e, y / 10)], [y + e], color=cols[s])
            scatter!(ax, [r.x], [y], color=r.dT2_rel < 0 ? cols[s] : :white, strokecolor=cols[s], strokewidth=1.5, markersize=9)
        end
    end
    hlines!(ax, [n200], color=:gray, linestyle=:dash, label="2σ noise SNR 200")
    hlines!(ax, [n100], color=:gray, linestyle=:dot, label="2σ noise SNR 100")
    j == 3 && Legend(fig[1, 4], ax, labelsize=9)
    axr = Axis(fig[2, j], xlabel=xlab[sw], ylabel="ΔS_mc / ΔS_cont", xscale=log10,
        title="($( 'd' + j - 1)) Monte Carlo against the continuum (paired ΔS)", xticks=(xt[sw], string.(xt[sw])))
    xlims!(axr, xr[sw]...)
    hspan!(axr, 0.9, 1.1, color=(:gray, 0.2))
    hlines!(axr, [1.0], color=:black)
    for (k, s) in enumerate(SHAPES), r in (m for m in mc if m.sweep == sw && m.shape == s)
        x = r.x * (1 + 0.04 * (k - 2))                                  # small horizontal offset per profile
        rt, σ = r.dS_mc / r.dS_cont, r.dS_sem / abs(r.dS_cont)
        errorbars!(axr, [x], [rt], [2σ], color=cols[s])
        scatter!(axr, [x], [rt], color=cols[s], markersize=9)
    end
end
Label(fig[3, :], "points: Monte Carlo (2σ bars); lines: continuum; open diamonds: continuum-only points; grey band (d–f): ±10 % (comparison timestep rule)", fontsize=11)
save(joinpath(OUT, "p3_7_params.png"), fig)
println("saved to $OUT")
