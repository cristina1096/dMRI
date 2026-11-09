# P3.7.3 inside (plan 2026-11-08-comparison-inside-cylinders.md, Task 6): D-dependence signature inside a single
# cylinder. Monte Carlo primary (p3_7_run.jl TASK=ddep: R = 1, bulk 1/2000, ρ ∈ {0.1, 1}, D ∈ {0.5, 1, 2, 3},
# conventional + step/linear/exponential at h/R ∈ {0.25, 0.5, 1}), the exact radial continuum as support.
# TE per ρ is fixed across D: the continuum conventional T₂ at D = 3 (readout 11 of 0:TE/10:1.5 TE, exact).
# Quantities per model: T₂_app(D) = −TE/ln S(TE); the D-signature q = T₂(0.5)/T₂(3) − 1, with σ_q from the two
# linearised T₂ errors added in quadrature (the D conditions share seeds and scaled draws, so this overstates σ_q);
# the layer-specific part q_layer − q_conv (what a D-dependence measurement could attribute to the layer).
# Detection band for a two-measurement ratio: 2·√(n(S₀.₅)² + n(S₃)²) with n = two-point relative noise at SNR 200
# (√2 × the single 2σ when both S are equal). Continuum-only: ρ = 0.01.
# Run: julia --project=research/baseline research/phase3/p3_7_ddep_analysis.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
include(joinpath(@__DIR__, "p3_7_shared.jl"))
using CairoMakie
const SRC = joinpath(REPO_ROOT, "research", "results", "phase3", "p3_7_ddep")
const OUT = phase3_outdir("p3_7_ddep_analysis")
const R1, BULK, SNR = 1.0, 1 / 2000, 200
const DS = (0.5, 1.0, 2.0, 3.0)
const SHAPES = ("step", "linear", "exponential")
const HRS = (0.25, 0.5, 1.0)

# ---- Monte Carlo at TE ----
rows = at_te_rows(SRC)
cond = Dict(c["label"] => c for c in readrows(joinpath(SRC, "conditions.csv")))
mc = [merge((rho=num(cond[r.label], "rho"), D=num(cond[r.label], "D")), r) for r in rows]
write_csv(joinpath(OUT, "at_te.csv"), mc)
sig = NamedTuple[]
for ρ in sort(unique(r.rho for r in mc)), m in unique(r.model for r in mc)
    p = [r for r in mc if r.rho == ρ && r.model == m]
    a = filter(r -> r.D == 0.5, p); b = filter(r -> r.D == 3.0, p)
    (isempty(a) || isempty(b)) && continue
    a, b = only(a), only(b)
    q = a.T2_mc / b.T2_mc - 1
    σq = (a.T2_mc / b.T2_mc) * sqrt((a.T2_sigma / a.T2_mc)^2 + (b.T2_sigma / b.T2_mc)^2)
    qc = a.T2_cont / b.T2_cont - 1
    band = 2 * sqrt(t2_noise_rel(a.S_mc, SNR)^2 + t2_noise_rel(b.S_mc, SNR)^2)
    push!(sig, (rho=ρ, model=m, shape=a.shape, hR=a.h / R1, q_mc=q, q_sigma=σq, q_cont=qc, z_q=σq > 1e-12 ? (q - qc) / σq : NaN,   # step at h = R: uniform R₂, T₂ independent of D, σ_q = 0
        band2_snr200=band))
end
for (i, s) in enumerate(sig)
    conv = only(c for c in sig if c.rho == s.rho && c.model == "conventional")
    sig[i] = merge(s, (q_layer_minus_conv_mc=s.q_mc - conv.q_mc, q_layer_minus_conv_cont=s.q_cont - conv.q_cont,
        detectable_snr200=abs(s.q_mc - conv.q_mc) > s.band2_snr200))
end
write_csv(joinpath(OUT, "signature.csv"), sig)
for s in sig
    @printf("ρ %-4g %-20s q = T₂(0.5)/T₂(3) − 1: %+.4f ± %.4f (cont %+.4f, z %+4.1f); layer − conv %+.4f (cont %+.4f); band %.3f %s\n",
        s.rho, s.model, s.q_mc, s.q_sigma, s.q_cont, s.z_q, s.q_layer_minus_conv_mc, s.q_layer_minus_conv_cont, s.band2_snr200,
        s.model == "conventional" ? "" : (s.detectable_snr200 ? "DETECTABLE" : "below band"))
end

# ---- continuum: every MC point plus ρ = 0.01 ----
tefor(ρ) = te_star(signal_fun(R=R1, D=3.0, surface=ρ, bulk=BULK); t_hi=50 / (BULK + 2ρ / R1))
cont = NamedTuple[]
for ρ in (0.01, 0.1, 1.0)
    te = tefor(ρ)
    for D in DS
        T2c = apparent_t2(signal_fun(R=R1, D=D, surface=ρ, bulk=BULK)(te), te)
        push!(cont, (rho=ρ, D=D, model="conventional", shape="conventional", hR=0.0, te=te, T2=T2c, continuum_only=ρ == 0.01))
        for s in SHAPES, hr in HRS
            T2 = apparent_t2(signal_fun(R=R1, D=D, layer=(s, ρ, hr * R1), bulk=BULK)(te), te)
            push!(cont, (rho=ρ, D=D, model="$(s)_hR$(hr)", shape=s, hR=hr, te=te, T2=T2, continuum_only=ρ == 0.01))
        end
    end
end
write_csv(joinpath(OUT, "continuum.csv"), cont)
for ρ in (0.01,)
    conv = Dict(c.D => c.T2 for c in cont if c.rho == ρ && c.model == "conventional")
    qc = conv[0.5] / conv[3.0] - 1
    for m in unique(c.model for c in cont if c.rho == ρ)
        T = Dict(c.D => c.T2 for c in cont if c.rho == ρ && c.model == m)
        @printf("continuum only ρ 0.01 %-20s q %+.5f  layer − conv %+.5f\n", m, T[0.5] / T[3.0] - 1, T[0.5] / T[3.0] - 1 - qc)
    end
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_conditions" => length(cond),
    "all_identical" => all(r -> r.identical, rows), "max_abs_z_S" => maximum(abs(r.z_S) for r in rows),
    "max_abs_z_q" => maximum((abs(s.z_q) for s in sig if !isnan(s.z_q)); init=0.0),
    "files" => ["at_te.csv", "signature.csv", "continuum.csv", "p3_7_ddep.png"]))

# ---- figure: T₂(D)/T₂(3) per ρ (top) and the layer − conventional signature against the band (bottom) ----
cols = Dict("conventional" => :black, "step" => Makie.wong_colors()[1], "linear" => Makie.wong_colors()[2],
    "exponential" => Makie.wong_colors()[3])
sty = Dict(0.25 => :solid, 0.5 => :dash, 1.0 => :dot, 0.0 => :solid)
fig = Figure(size=(1400, 820))
for (j, ρ) in enumerate((0.01, 0.1, 1.0))
    ax = Axis(fig[1, j], xlabel="D (µm²/ms)", ylabel="T₂,app(D) / T₂,app(3)", title="ρ = $(ρ)" * (ρ == 0.01 ? " (continuum only)" : ""))
    for m in unique(c.model for c in cont if c.rho == ρ)
        sel = sort([c for c in cont if c.rho == ρ && c.model == m], by=c -> c.D)
        ref = only(c for c in sel if c.D == 3.0).T2
        lines!(ax, [c.D for c in sel], [c.T2 / ref for c in sel], color=cols[sel[1].shape], linestyle=sty[sel[1].hR],
            label=m == "conventional" ? "conventional" : "$(sel[1].shape) h/R = $(sel[1].hR)")
        p = sort([r for r in mc if r.rho == ρ && r.model == m], by=r -> r.D)
        isempty(p) && continue
        refm = only(r for r in p if r.D == 3.0).T2_mc
        scatter!(ax, [r.D for r in p], [r.T2_mc / refm for r in p], color=cols[sel[1].shape], markersize=7)
    end
    j == 3 && Legend(fig[1, 4], ax, labelsize=9, merge=true)
    ax2 = Axis(fig[2, j], xlabel="h/R", ylabel="q_layer − q_conv", xscale=log10, xticks=([0.25, 0.5, 1.0], ["0.25", "0.5", "1"]),
        title="ρ = $(ρ): q_layer − q_conv, q = T₂(0.5)/T₂(3) − 1")
    if ρ == 0.01
        conv = Dict(c.D => c.T2 for c in cont if c.rho == ρ && c.model == "conventional")
        for s in SHAPES
            q = [begin T = Dict(c.D => c.T2 for c in cont if c.rho == ρ && c.model == "$(s)_hR$(hr)"); T[0.5] / T[3.0] - conv[0.5] / conv[3.0] end for hr in HRS]
            lines!(ax2, collect(HRS), q, color=cols[s])
        end
    else
        for s in SHAPES
            sel = sort([g for g in sig if g.rho == ρ && g.shape == s], by=g -> g.hR)
            lines!(ax2, [g.hR for g in sel], [g.q_layer_minus_conv_cont for g in sel], color=cols[s])
            errorbars!(ax2, [g.hR for g in sel], [g.q_layer_minus_conv_mc for g in sel], [2g.q_sigma for g in sel], color=cols[s])
            scatter!(ax2, [g.hR for g in sel], [g.q_layer_minus_conv_mc for g in sel], color=cols[s], markersize=8)
        end
        b = maximum(g.band2_snr200 for g in sig if g.rho == ρ)
        hspan!(ax2, -b, b, color=(:gray, 0.2))
    end
    hlines!(ax2, [0.0], color=:black)
end
Label(fig[3, :], "lines: continuum; points: Monte Carlo (2σ bars, σ overstated: D conditions share seeds); grey band: two-measurement 2σ noise at SNR 200", fontsize=11)
save(joinpath(OUT, "p3_7_ddep.png"), fig)
println("saved to $OUT")
