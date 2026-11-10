# P3.5.4 (Task 5): D dependence of overlapping outer layers, square 2.5 (gap 0.5), ρ = 1, bulk 1/2000, D 0.5 vs 3.
# Sources (independent seeds; same TE convention, the annulus conventional T₂ at D = 3, read from each conditions.csv):
#   p3_5_ddep (h ∈ {0.15, 0.35}) and p3_7b_ddep (h ∈ {0.25, 0.5}, rows D = 0.5 and 3 of ρ = 1). Each source's layer
#   models are compared with ITS OWN conventional run (paired).
# q = T₂(0.5)/T₂(3) − 1; Δq = q_layer − q_conv = (1 + q_conv)[(1 + δ₀.₅)/(1 + δ₃) − 1] with δ the paired ΔT₂/T₂ per D;
# σ_Δq ≈ (1 + q_conv)·√(σ_δ₀.₅² + σ_δ₃²); band = 2·√(n(S₀.₅)² + n(S₃)²) at SNR 200 (n at the conventional S, as 3.7b).
# Walls (P2.5.9, ρ = 1): thick or overlapping layers nearly D-independent (q_layer ≤ 4.6 % for h/w ≥ 0.5 vs q_conv 47 %).
# Run: julia --project=research/baseline research/phase3/p3_5_ddep_analysis.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
include(joinpath(@__DIR__, "p3_7_shared.jl"))
using CairoMakie
const RES = joinpath(REPO_ROOT, "research", "results", "phase3")
const OUT = phase3_outdir("p3_5_ddep_analysis")
const GAP = 0.5

function at_te(src, labelfilter)
    conds = [c for c in readrows(joinpath(RES, src, "conditions.csv")) if labelfilter(c["label"])]
    curves = readrows(joinpath(RES, src, "curves.csv"))
    out = NamedTuple[]
    for c in conds
        te, D = num(c, "te"), num(c, "D")
        sel = [r for r in curves if r["label"] == c["label"]]
        ts = sort(unique(num(r, "t") for r in sel))
        k = findfirst(t -> abs(t - te) <= 1e-9 * te, ts)
        isnothing(k) && error("$(c["label"]): no readout at TE")
        at(m) = only(r for r in sel if r["model"] == m && num(r, "t") == ts[k])
        Sc = num(at("conventional"), "S_mean"); T2c = apparent_t2(Sc, te)
        for m in unique(r["model"] for r in sel)
            r = at(m); S = num(r, "S_mean"); T2 = apparent_t2(S, te)
            lay = m != "conventional"
            push!(out, (source=src, D=D, te=te, model=m, shape=r["shape"], h=num(r, "h"), T2=T2,
                T2_sigma=te * num(r, "S_sem") / (S * log(S)^2), delta=lay ? (T2 - T2c) / T2c : 0.0,
                delta_sigma=lay ? te * num(r, "dS_sem") / (S * log(S)^2) / T2c : 0.0, noise_rel=t2_noise_rel(Sc, 200)))
        end
    end
    return out
end
rows = vcat(at_te("p3_5_ddep", l -> true), at_te("p3_7b_ddep", l -> startswith(l, "s2.5_rho1.0_D") && (endswith(l, "D0.5") || endswith(l, "D3.0"))))
write_csv(joinpath(OUT, "at_te.csv"), rows)
sig = NamedTuple[]
for src in unique(r.source for r in rows), m in unique(r.model for r in rows if r.source == src)
    m == "conventional" && continue
    g(D, mm) = only(r for r in rows if r.source == src && r.model == mm && r.D == D)
    a, b, ca, cb = g(0.5, m), g(3.0, m), g(0.5, "conventional"), g(3.0, "conventional")
    qc = ca.T2 / cb.T2 - 1
    q = a.T2 / b.T2 - 1
    Δq = (1 + qc) * ((1 + a.delta) / (1 + b.delta) - 1)
    σ = (1 + qc) * sqrt(a.delta_sigma^2 + b.delta_sigma^2)
    band = 2 * sqrt(a.noise_rel^2 + b.noise_rel^2)
    push!(sig, (source=src, model=m, shape=a.shape, h=a.h, h_over_gap=a.h / GAP, q=q, q_conv=qc, dq=Δq, dq_sigma=σ,
        band2_snr200=band, ratio_to_band=abs(Δq) / band, detectable=abs(Δq) > band, q_over_qconv=q / qc))
end
sort!(sig, by=s -> (s.shape, s.h))
write_csv(joinpath(OUT, "signature.csv"), sig)
foreach(s -> @printf("%-11s %-12s h/gap %.2f  q %+.4f (q_conv %.4f; q/q_conv %.3f)  Δq %+.4f ± %.4f  band %.4f (%.2f×)\n", s.source, s.shape,
    s.h_over_gap, s.q, s.q_conv, s.q_over_qconv, s.dq, s.dq_sigma, s.band2_snr200, s.ratio_to_band), sig)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "files" => ["at_te.csv", "signature.csv", "p3_5_ddep.png"]))

cols = Dict("step" => Makie.wong_colors()[1], "linear" => Makie.wong_colors()[2], "exponential" => Makie.wong_colors()[3])
fig = Figure(size=(1100, 450))
ax = Axis(fig[1, 1], xlabel="h / gap", ylabel="q_layer / q_conv", title="(a) remaining D dependence (q = T₂(0.5)/T₂(3) − 1)")
ax2 = Axis(fig[1, 2], xlabel="h / gap", ylabel="Δq = q_layer − q_conv", title="(b) D signature vs two-measurement 2σ band (SNR 200)")
b = maximum(s.band2_snr200 for s in sig)
hspan!(ax2, -b, b, color=(:gray, 0.2)); hlines!(ax2, [0.0], color=:black)
vlines!(ax, [0.5], color=:gray, linestyle=:dash); vlines!(ax2, [0.5], color=:gray, linestyle=:dash)
for sh in ("step", "linear", "exponential")
    p = sort([s for s in sig if s.shape == sh], by=s -> s.h)
    scatterlines!(ax, [s.h_over_gap for s in p], [s.q_over_qconv for s in p], color=cols[sh], label=sh)
    errorbars!(ax2, [s.h_over_gap for s in p], [s.dq for s in p], [2s.dq_sigma for s in p], color=cols[sh])
    scatterlines!(ax2, [s.h_over_gap for s in p], [s.dq for s in p], color=cols[sh])
end
axislegend(ax, position=:rt, labelsize=9)
save(joinpath(OUT, "p3_5_ddep.png"), fig)
println("saved to $OUT")
