# P3.7.1 (plan 2026-11-08-comparison-inside-cylinders.md, Task 3): Fig. 9A analogue inside single cylinders — apparent
# T₂ against axon radius for the conventional model and the three layer profiles, Monte Carlo primary (p3_7_run.jl,
# TASK=fig9a) with the exact radial continuum as support.
# Settings: ρ = 1/140, D = 3, R2_bulk = 1/2000; matchings: same centroid c̄ = 0.1 (step h 0.2, linear 0.3, exponential by
# bisection with support R) at R ∈ {0.5, 1, 2}, and same h = 0.2 at R = 1.
# Quantities at TE* = continuum conventional T₂ (readout 11 of 0:TE*/10:1.5 TE*, exact):
#   T₂_app = −TE*/ln S(TE*);  σ(T₂_app) = TE*·SEM_S/(S ln²S)  (linearised);
#   ΔT₂/T₂ = (T₂_layer − T₂_conv)/T₂_conv, σ from the paired SEM of ΔS: σ ≈ TE*·SEM_ΔS/(S_layer ln²S_layer)/T₂_conv
#     (the common trajectory noise cancels in the pair; linearised about the layer signal);
#   2σ two-point noise at SNR 100, 200: 2·√(1 + r²)/(SNR·r·|ln r|), r = S_conv(TE*) = e^(−1) by construction;
#   z(S) = (S_mc − S_cont)/SEM_S,  z(ΔS) = (ΔS_mc − ΔS_cont)/SEM_ΔS.
# Continuum dense curves over R ∈ [0.25, 3] for every profile in both matchings; a point is skipped when the layer
# support exceeds R (linear h = 0.3 needs R ≥ 0.3, h = 0.2 needs R ≥ 0.2).
# Run: julia --project=research/baseline research/phase3/p3_7_fig9a_analysis.jl  (works on partial MC results)

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
using CairoMakie
const SRC = joinpath(REPO_ROOT, "research", "results", "phase3", "p3_7_fig9a")
const OUT = phase3_outdir("p3_7_fig9a_analysis")
const RHO, D3, BULK, CB, HSAME = 1 / 140, 3.0, 1 / 2000, 0.1, 0.2
const SHAPES = ("step", "linear", "exponential")
readrows(f) = (ls = readlines(f); h = split(ls[1], ','); [Dict(zip(h, split(l, ','))) for l in ls[2:end]])
num(r, k) = parse(Float64, r[k])

cbar(shape, h, R) = shape == "step" ? h / 2 : shape == "linear" ? h / 3 : (c = R / h; h * (1 - (c + 1) * exp(-c)) / -expm1(-c))
function h_for_cbar(shape, cb, R)
    shape == "step" && return 2cb
    shape == "linear" && return 3cb
    lo, hi = 1e-6, 100R
    for _ in 1:200
        mid = sqrt(lo * hi); cbar("exponential", mid, R) < cb ? (lo = mid) : (hi = mid)
    end
    return sqrt(lo * hi)
end
"S(t) closure from one eigendecomposition of the radial operator (same discretisation as radial_signal)."
function signal_fun(; R, D, layer=nothing, surface=0.0, bulk=0.0, N=2000)
    (M, A) = radial_operator(; R=R, D=D, layer=layer, surface=surface, N=N)
    F = eigen(M)
    c2 = (F.vectors' * sqrt.(A)) .^ 2 ./ sum(A)
    return t -> sum(c2 .* exp.(-F.values .* t)) * exp(-bulk * t)
end
t2app_sigma(S, sem, te) = te * sem / (S * log(S)^2)

# ---- Monte Carlo at TE* ----
conds = readrows(joinpath(SRC, "conditions.csv"))
curves = readrows(joinpath(SRC, "curves.csv"))
rows = NamedTuple[]
for c in conds
    label, R, te = c["label"], num(c, "R"), num(c, "te")
    sel = [r for r in curves if r["label"] == label]
    ts = sort(unique(num(r, "t") for r in sel))
    k = argmin(abs.(ts .- te))
    abs(ts[k] - te) <= 1e-9 * te || error("$label: no readout at TE*")
    at(m) = only(r for r in sel if r["model"] == m && num(r, "t") == ts[k])
    conv = at("conventional")
    Sc = num(conv, "S_mean")
    T2c = apparent_t2(Sc, te)
    for m in unique(r["model"] for r in sel)
        r = at(m)
        S, sem, Scont = num(r, "S_mean"), num(r, "S_sem"), num(r, "S_cont")
        T2 = apparent_t2(S, te)
        layer = m != "conventional"
        dT2 = layer ? (T2 - T2c) / T2c : NaN
        dT2σ = layer ? t2app_sigma(S, num(r, "dS_sem"), te) / T2c : NaN
        dT2cont = layer ? (apparent_t2(Scont, te) - apparent_t2(num(conv, "S_cont"), te)) / apparent_t2(num(conv, "S_cont"), te) : NaN
        push!(rows, (label=label, matching=occursin("_h", label) ? "same_h" : "same_centroid", R=R, te=te, model=m, shape=r["shape"],
            h=num(r, "h"), cbar=num(r, "cbar"), S_mc=S, S_sem=sem, S_cont=Scont, z_S=(S - Scont) / sem,
            T2_mc=T2, T2_sigma=t2app_sigma(S, sem, te), T2_cont=apparent_t2(Scont, te),
            dT2_rel=dT2, dT2_rel_sigma=dT2σ, dT2_rel_cont=dT2cont,
            noise2_snr100=2 * t2_noise_rel(Sc, 100), noise2_snr200=2 * t2_noise_rel(Sc, 200),
            detectable_snr200=layer ? abs(dT2) > 2 * t2_noise_rel(Sc, 200) : false,
            dS_mc=layer ? num(r, "dS_mean") : NaN, dS_sem=layer ? num(r, "dS_sem") : NaN, dS_cont=layer ? num(r, "dS_cont") : NaN,
            z_dS=layer ? (num(r, "dS_mean") - num(r, "dS_cont")) / num(r, "dS_sem") : NaN,
            tau=num(c, "tau"), identical=c["identical_trajectories"] == "true"))
    end
end
write_csv(joinpath(OUT, "at_te.csv"), rows)
for r in rows
    @printf("%-11s %-24s T₂ %7.2f ± %.2f (cont %7.2f, z_S %+5.1f)  ΔT₂/T₂ %+7.3f%% ± %.3f%% (cont %+7.3f%%, z_ΔS %+5.1f)\n",
        r.label, r.model, r.T2_mc, r.T2_sigma, r.T2_cont, r.z_S, 100r.dT2_rel, 100r.dT2_rel_sigma, 100r.dT2_rel_cont, r.z_dS)
end

# Group A review I2: the per-point |z| hides a same-sign pattern, so report ΔS_mc/ΔS_cont ± SEM_ΔS/|ΔS_cont| per point and
# per condition (the profiles of one condition share trajectories, so their ratios are correlated, not independent)
bias = [(label=r.label, model=r.model, tau=r.tau, ratio=r.dS_mc / r.dS_cont, ratio_sigma=r.dS_sem / abs(r.dS_cont), z_dS=r.z_dS)
        for r in rows if r.model != "conventional"]
write_csv(joinpath(OUT, "dS_ratio.csv"), bias)
for l in unique(b.label for b in bias)
    sel = [b for b in bias if b.label == l]
    @printf("%-11s ΔS_mc/ΔS_cont: %s  (z_ΔS signs %s)\n", l, join([@sprintf("%.3f ± %.3f", b.ratio, b.ratio_sigma) for b in sel], ", "),
        join([b.z_dS > 0 ? "+" : "−" for b in sel]))
end

# ---- continuum dense curves ----
Rs = collect(range(0.25, 3.0, length=45))
dense = NamedTuple[]
for R in Rs
    Sconv = signal_fun(R=R, D=D3, surface=RHO, bulk=BULK)
    te = te_star(Sconv; t_hi=50 / (BULK + 2RHO / R))
    T2c = apparent_t2(Sconv(te), te)
    push!(dense, (R=R, matching="none", model="conventional", shape="conventional", h=0.0, te=te, T2=T2c, dT2_rel=0.0))
    for (matching, hf) in (("same_centroid", s -> h_for_cbar(s, CB, R)), ("same_h", s -> HSAME)), s in SHAPES
        h = hf(s)
        radial_support(s, h, R) <= R || continue
        S = signal_fun(R=R, D=D3, layer=(s, RHO, h), bulk=BULK)(te)
        T2 = apparent_t2(S, te)
        push!(dense, (R=R, matching=matching, model=s, shape=s, h=h, te=te, T2=T2, dT2_rel=(T2 - T2c) / T2c))
    end
end
write_csv(joinpath(OUT, "continuum_dense.csv"), dense)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_conditions" => length(conds),
    "all_identical" => all(r -> r.identical, rows), "max_abs_z_S" => maximum(abs(r.z_S) for r in rows),
    "max_abs_z_dS" => maximum((abs(r.z_dS) for r in rows if !isnan(r.z_dS)); init=0.0),
    "noise2_snr200" => 2 * t2_noise_rel(exp(-1), 200), "noise2_snr100" => 2 * t2_noise_rel(exp(-1), 100),
    "files" => ["at_te.csv", "dS_ratio.csv", "continuum_dense.csv", "p3_7_fig9a.png"]))

# ---- figure ----
cols = Dict("conventional" => :black, "step" => Makie.wong_colors()[1], "linear" => Makie.wong_colors()[2],
    "exponential" => Makie.wong_colors()[3])
sty = Dict("same_centroid" => :solid, "same_h" => :dash, "none" => :solid)
mk = Dict("same_centroid" => :circle, "same_h" => :utriangle)
fig = Figure(size=(1150, 820))
axa = Axis(fig[1, 1], xlabel="R (µm)", ylabel="T₂,app at TE* (ms)", title="(a) apparent T₂ vs radius (ρ = 1/140, D = 3)")
axb = Axis(fig[1, 2], xlabel="R (µm)", ylabel="|ΔT₂/T₂| vs conventional (%)", yscale=log10,
    title="(b) |layer − conventional|, paired (filled: ΔT₂ < 0, open: > 0)")
axc = Axis(fig[2, 1:3], xlabel="condition / model", ylabel="z (MC − continuum)/SEM", title="(c) Monte Carlo against the continuum")
for (matching, s) in vcat([("none", "conventional")], [(m, s) for m in ("same_centroid", "same_h") for s in SHAPES])
    sel = [d for d in dense if d.matching == matching && d.shape == s]
    isempty(sel) && continue
    lab = s == "conventional" ? "conventional" : "$(s) " * (matching == "same_h" ? "h = 0.2" : "c̄ = 0.1")
    lines!(axa, [d.R for d in sel], [d.T2 for d in sel], color=cols[s], linestyle=sty[matching], label=lab)
    s == "conventional" || lines!(axb, [d.R for d in sel], [max(100abs(d.dT2_rel), 1e-6) for d in sel], color=cols[s], linestyle=sty[matching], label=lab)
end
for r in rows
    m = r.model == "conventional" ? :rect : mk[r.matching]
    errorbars!(axa, [r.R], [r.T2_mc], [2r.T2_sigma], color=cols[r.shape])
    scatter!(axa, [r.R], [r.T2_mc], color=cols[r.shape], marker=m, markersize=9)
    if r.model != "conventional"
        y, e = 100abs(r.dT2_rel), 200r.dT2_rel_sigma
        rangebars!(axb, [r.R], [max(y - e, y / 10)], [y + e], color=cols[r.shape])
        scatter!(axb, [r.R], [y], color=r.dT2_rel < 0 ? cols[r.shape] : :white, strokecolor=cols[r.shape], strokewidth=1.5,
            marker=m, markersize=9)
    end
end
n200, n100 = 100 * 2 * t2_noise_rel(exp(-1), 200), 100 * 2 * t2_noise_rel(exp(-1), 100)
hlines!(axb, [n200], color=:gray, linestyle=:dash, label="2σ noise SNR 200")
hlines!(axb, [n100], color=:gray, linestyle=:dot, label="2σ noise SNR 100")
axislegend(axa, position=:lt, labelsize=9)
Legend(fig[1, 3], axb, labelsize=9)
xs = 1:length(rows)
scatter!(axc, xs, [r.z_S for r in rows], color=[cols[r.shape] for r in rows], marker=:circle, label="z(S)")
zd = [r.z_dS for r in rows]
ok = .!isnan.(zd)
scatter!(axc, xs[ok], zd[ok], color=[cols[r.shape] for r in rows[ok]], marker=:xcross, markersize=12, label="z(ΔS)")
hlines!(axc, [-2, 2], color=:gray, linestyle=:dash)
axc.xticks = (xs, ["$(r.label)\n$(r.shape)" for r in rows])
axc.xticklabelrotation = π / 3
axc.xticklabelsize = 8
axislegend(axc, position=:lt, labelsize=9)
save(joinpath(OUT, "p3_7_fig9a.png"), fig)
println("saved to $OUT")
