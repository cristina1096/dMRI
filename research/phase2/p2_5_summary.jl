# P2.5 summary — observability of the profile (spatial distribution) at fixed ρ, from the four sweeps.
#
# (1) Per sweep: |ΔS| at its peak vs the swept parameter, per shape, with 3·SEM and the MRI thresholds 2/SNR
#     (SNR 100: 0.02, SNR 200: 0.01). (2) h sweep replotted against the profile centroid (P2.5.4). (3) ΔS(t) curves
#     for the h sweep. (4) observability.csv: one row per configuration.
#
# Run: julia --project=research/baseline research/phase2/p2_5_summary.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_5_summary")
const SWEEPS = ("h", "w", "D", "rho")
R = joinpath(REPO_ROOT, "research", "results", "phase2")
load_rows(s) = [Dict(zip(split(readlines(joinpath(R, "p2_5_sensitivity_" * s, "rows.csv"))[1], ","), split(l, ",")))
                for l in readlines(joinpath(R, "p2_5_sensitivity_" * s, "rows.csv"))[2:end]]
num(d, k) = parse(Float64, d[k])
classify(a, sem) = a >= 0.02 ? "SNR 100" : a >= 0.01 ? "SNR 200" : a >= 3sem ? "MC only" : "not resolved"

obs = NamedTuple[]
for s in SWEEPS, r in load_rows(s)
    a, se = abs(num(r, "dS_peak")), num(r, "dS_peak_sem")
    push!(obs, (sweep=s, shape=r["shape"], h_um=num(r, "h_um"), w_um=num(r, "w_um"), D=num(r, "D"), rho=num(r, "rho"),
        tau_ms=num(r, "tau_ms"), centroid_um=num(r, "centroid_um"), dS_peak=num(r, "dS_peak"), dS_peak_sem=se,
        dS_peak_t_ms=num(r, "dS_peak_t_ms"), level=classify(a, se), monoexp_residual=num(r, "monoexp_residual"),
        base_curvature=num(r, "base_curvature")))
    o = obs[end]
    @printf("%-3s %-11s h=%.3f w=%4.1f D=%.0f ρ=%.3f: ΔS = %+.2e ± %.0e at %3.0f ms → %s\n",
        s, o.shape, o.h_um, o.w_um, o.D, o.rho, o.dS_peak, se, o.dS_peak_t_ms, o.level)
end
write_csv(joinpath(OUT, "observability.csv"), obs)
counts = Dict(l => count(o -> o.level == l, obs) for l in ("SNR 100", "SNR 200", "MC only", "not resolved"))
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rows" => obs, "counts" => counts,
    "thresholds" => Dict("SNR 100" => 0.02, "SNR 200" => 0.01, "MC" => "3·paired SEM")))
println("counts: ", counts)

using CairoMakie
param = Dict("h" => :h_um, "w" => :w_um, "D" => :D, "rho" => :rho)
fig = Figure(size=(1200, 330))
for (j, s) in enumerate(SWEEPS)
    ax = Axis(fig[1, j], xscale=log10, yscale=log10, xlabel=string(param[s]), ylabel="|ΔS| at peak", title="sweep $s")
    for shape in ("step", "linear", "exponential")
        sel = [o for o in obs if o.sweep == s && o.shape == shape]
        isempty(sel) && continue
        x = [getproperty(o, param[s]) for o in sel]
        scatterlines!(ax, x, max.(abs.([o.dS_peak for o in sel]), 1e-7), label=shape)
    end
    hlines!(ax, [0.02, 0.01], color=:gray, linestyle=:dash)
    j == 1 && axislegend(ax, position=:lt, labelsize=8)
end
save(joinpath(OUT, "p2_5_sweeps.png"), fig)

fig2 = Figure(size=(560, 380))
ax = Axis(fig2[1, 1], xscale=log10, yscale=log10, xlabel="centroid of the relaxation (µm)", ylabel="|ΔS| at peak", title="h sweep vs centroid (P2.5.4)")
for shape in ("step", "linear", "exponential")
    sel = [o for o in obs if o.sweep == "h" && o.shape == shape]
    scatterlines!(ax, [o.centroid_um for o in sel], max.(abs.([o.dS_peak for o in sel]), 1e-7), label=shape)
end
hlines!(ax, [0.02, 0.01], color=:gray, linestyle=:dash); axislegend(ax, position=:lt)
save(joinpath(OUT, "p2_5_centroid.png"), fig2)

curves = JSON.parsefile(joinpath(R, "p2_5_sensitivity_h", "curves.json"))
fig3 = Figure(size=(560, 380))
ax = Axis(fig3[1, 1], xlabel="t (ms)", ylabel="ΔS = S_layer − S_base", title="h sweep: ΔS(t) at fixed ρ = 0.01 µm/ms")
for c in curves
    cfg = c["config"]
    lines!(ax, Float64.(c["times"]), Float64.(c["dS"]), label="$(cfg["shape"]) h=$(cfg["h"])",
        linestyle=cfg["shape"] == "step" ? :dash : cfg["shape"] == "linear" ? :solid : :dot)
end
hlines!(ax, [0.0], color=:gray)          # MRI thresholds (±0.01, ±0.02) are 25× outside this range; see p2_5_sweeps.png
Legend(fig3[1, 2], ax, labelsize=7, rowgap=0)
save(joinpath(OUT, "p2_5_curves.png"), fig3)
println("saved to $OUT")
