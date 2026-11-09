# P3.7b Task 6 (P3.7.3 outside): D-dependence signature OUTSIDE cylinders, square array s = 2.5 (gap 0.5, f = 0.503).
# Monte Carlo only (p3_7b_run.jl TASK=ddep): the annulus failed validation at s = 2.5 and sets only TE (user decision
# 2026-11-09). R = 1, bulk 1/2000, ρ ∈ {0.1, 1}, D ∈ {0.5, 1, 2, 3}, conventional + outer step/linear/exponential at
# h/R ∈ {0.25, 0.5} (h/R = 0.1 dropped, user decision 2026-11-10). TE per ρ = annulus conventional T₂ at D = 3, the same
# for every D (exact readout 11 of 0:TE/10:1.5 TE). Every D condition has its own seeds (independent).
# Quantities per model: T₂_app(D) = −TE/ln S(TE), σ = TE·SEM/(S ln²S);
#   q = T₂(0.5)/T₂(3) − 1, σ_q = (1 + q)·√((σ_T₂/T₂)₀.₅² + (σ_T₂/T₂)₃²) (independent D conditions);
#   the layer-specific part Δq = q_layer − q_conv. Within a D condition the layer and conventional runs are paired, so
#   with δ(D) = (T₂_layer − T₂_conv)/T₂_conv (paired σ_δ from SEM_ΔS):
#   Δq = (1 + q_conv)·[(1 + δ₀.₅)/(1 + δ₃) − 1] ≈ (1 + q_conv)(δ₀.₅ − δ₃), σ_Δq ≈ (1 + q_conv)·√(σ_δ₀.₅² + σ_δ₃²);
#   two-measurement band at SNR 200: 2·√(n(S₀.₅)² + n(S₃)²), n = two-point relative noise of the conventional S.
#   Δq is compared with the band; cases at the band edge (|Δq| ≥ 0.9·band but ≤ band) are named explicitly (3.7a lesson I3).
# The comparison is at the same ρ for every D (ρ independent of D; MCMR θ = ρ√(π/D) changes with D).
# Run: julia --project=research/baseline research/phase3/p3_7b_ddep_analysis.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
include(joinpath(@__DIR__, "p3_7_shared.jl"))
using CairoMakie
const SRC = joinpath(REPO_ROOT, "research", "results", "phase3", "p3_7b_ddep")
const OUT = phase3_outdir("p3_7b_ddep_analysis")
const SNR = 200

conds = readrows(joinpath(SRC, "conditions.csv"))
curves = readrows(joinpath(SRC, "curves.csv"))
rows = NamedTuple[]
for c in conds
    label, te, ρ, D = c["label"], num(c, "te"), num(c, "rho"), num(c, "D")
    sel = [r for r in curves if r["label"] == label]
    ts = sort(unique(num(r, "t") for r in sel))
    k = findfirst(t -> abs(t - te) <= 1e-9 * te, ts)
    isnothing(k) && error("$label: no readout at TE")
    at(m) = only(r for r in sel if r["model"] == m && num(r, "t") == ts[k])
    conv = at("conventional")
    Sc = num(conv, "S_mean")
    T2c = apparent_t2(Sc, te)
    for m in unique(r["model"] for r in sel)
        r = at(m)
        S, sem = num(r, "S_mean"), num(r, "S_sem")
        T2 = apparent_t2(S, te)
        layer = m != "conventional"
        push!(rows, (label=label, rho=ρ, D=D, te=te, model=m, shape=r["shape"], hR=num(r, "h"), S_mc=S, S_sem=sem,
            T2=T2, T2_sigma=te * sem / (S * log(S)^2),
            delta=layer ? (T2 - T2c) / T2c : 0.0, delta_sigma=layer ? te * num(r, "dS_sem") / (S * log(S)^2) / T2c : 0.0,
            noise_rel=t2_noise_rel(Sc, SNR), tau=num(c, "tau"), cbar_min=num(c, "cbar_min"), nseeds=Int(num(c, "nseeds")),
            identical=c["identical_trajectories"] == "true"))
    end
end
write_csv(joinpath(OUT, "at_te.csv"), rows)

sig = NamedTuple[]
for ρ in sort(unique(r.rho for r in rows)), m in unique(r.model for r in rows)
    get(D) = (q = filter(r -> r.rho == ρ && r.model == m && r.D == D, rows); isempty(q) ? nothing : only(q))
    a, b = get(0.5), get(3.0)
    (isnothing(a) || isnothing(b)) && continue
    ca = only(filter(r -> r.rho == ρ && r.model == "conventional" && r.D == 0.5, rows))
    cb = only(filter(r -> r.rho == ρ && r.model == "conventional" && r.D == 3.0, rows))
    q = a.T2 / b.T2 - 1
    σq = (1 + q) * sqrt((a.T2_sigma / a.T2)^2 + (b.T2_sigma / b.T2)^2)
    qc = ca.T2 / cb.T2 - 1
    Δq = (1 + qc) * ((1 + a.delta) / (1 + b.delta) - 1)
    σΔq = (1 + qc) * sqrt(a.delta_sigma^2 + b.delta_sigma^2)
    band = 2 * sqrt(a.noise_rel^2 + b.noise_rel^2)
    status = m == "conventional" ? "reference" : abs(Δq) > band ? "detectable" : abs(Δq) >= 0.9band ? "band edge" : "below band"
    push!(sig, (rho=ρ, model=m, shape=a.shape, hR=a.hR, q=q, q_sigma=σq, q_conv=qc, dq=m == "conventional" ? 0.0 : Δq,
        dq_sigma=m == "conventional" ? 0.0 : σΔq, band2_snr200=band, ratio_to_band=abs(Δq) / band, status=status))
end
write_csv(joinpath(OUT, "signature.csv"), sig)
for s in sig
    @printf("ρ %-4g %-16s q = %+.4f ± %.4f  Δq = %+.4f ± %.4f  band %.4f  (%.2f×)  %s\n", s.rho, s.model, s.q, s.q_sigma, s.dq,
        s.dq_sigma, s.band2_snr200, s.ratio_to_band, s.status)
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_conditions" => length(conds),
    "all_identical" => all(r -> r.identical, rows), "support" => "none at s = 2.5 (annulus failed validation; TE only)",
    "cbar_min" => minimum(r.cbar_min for r in rows), "cbar_min_flag" => minimum(r.cbar_min for r in rows) < 0.05 ? "τ rule extrapolated (c̄ < 0.05)" : "within the measured rule",
    "files" => ["at_te.csv", "signature.csv", "p3_7b_ddep.png"]))

# ---- figure ----
cols = Dict("conventional" => :black, "step" => Makie.wong_colors()[1], "linear" => Makie.wong_colors()[2],
    "exponential" => Makie.wong_colors()[3])
sty = Dict(0.0 => :solid, 0.25 => :solid, 0.5 => :dash)
fig = Figure(size=(1300, 800))
for (j, ρ) in enumerate(sort(unique(r.rho for r in rows)))
    ax = Axis(fig[1, j], xlabel="D (µm²/ms)", ylabel="T₂,app(D) / T₂,app(3)", title="outside, s = 2.5, ρ = $(ρ) (MC only)")
    for m in unique(r.model for r in rows)
        p = sort([r for r in rows if r.rho == ρ && r.model == m], by=r -> r.D)
        ref = only(r for r in p if r.D == 3.0).T2
        scatterlines!(ax, [r.D for r in p], [r.T2 / ref for r in p], color=cols[p[1].shape], linestyle=sty[p[1].hR],
            label=m == "conventional" ? "conventional" : "$(p[1].shape) h/R = $(p[1].hR)")
    end
    j == 2 && Legend(fig[1, 3], ax, labelsize=9)
    ax2 = Axis(fig[2, j], xlabel="h/R", ylabel="Δq = q_layer − q_conv", xticks=[0.25, 0.5],
        title="ρ = $(ρ): Δq (2σ bars), grey: two-measurement 2σ noise, SNR 200")
    b = maximum(g.band2_snr200 for g in sig if g.rho == ρ)
    hspan!(ax2, -b, b, color=(:gray, 0.2))
    hlines!(ax2, [0.0], color=:black)
    for (k, sh) in enumerate(("step", "linear", "exponential"))
        g = sort([x for x in sig if x.rho == ρ && x.shape == sh], by=x -> x.hR)
        xs = [x.hR for x in g] .+ 0.01 * (k - 2)
        errorbars!(ax2, xs, [x.dq for x in g], [2x.dq_sigma for x in g], color=cols[sh])
        scatterlines!(ax2, xs, [x.dq for x in g], color=cols[sh])
    end
end
save(joinpath(OUT, "p3_7b_ddep.png"), fig)
println("saved to $OUT")
