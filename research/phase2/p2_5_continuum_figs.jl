# Figures for the continuum recomputation (p2_5_continuum.jl), read from its CSVs.
# (1) P2.5 sweeps: continuum peak ΔS (lines) vs the Monte Carlo at τ = 1e-2 (dots). (2) Wider map in log10 coordinates:
# max |ΔS| at fixed ρ and the best-fit-conventional residual, with the 2σ thresholds of SNR 200 (0.01) and 100 (0.02).
# (3) Apparent T2 (TE = 50): ΔT2/T2 continuum vs Monte Carlo at τ = 1e-2, same h and same centroid.
#
# Run: julia --project=research/baseline research/phase2/p2_5_continuum_figs.jl

include(joinpath(@__DIR__, "harness.jl"))
using CairoMakie

const OUT = phase2_outdir("p2_5_continuum")
const R = joinpath(REPO_ROOT, "research", "results", "phase2")
csvrows(f) = (L = readlines(f); hdr = split(L[1], ","); [Dict(zip(hdr, split(l, ","))) for l in L[2:end]])
num(d, k) = parse(Float64, d[k])
colors = Dict("step" => :steelblue, "linear" => :darkorange, "exponential" => :seagreen)

p25 = csvrows(joinpath(OUT, "p2_5_continuum.csv"))
key = Dict("h" => "h_um", "w" => "w_um", "D" => "D", "rho" => "rho")
fig = Figure(size=(1300, 380))
for (j, s) in enumerate(("h", "w", "D", "rho"))
    ax = Axis(fig[1, j], xscale=log10, xlabel=key[s], ylabel="ΔS at peak", title="sweep $s: continuum (lines), MC τ=1e-2 (dots)")
    mc = csvrows(joinpath(R, "p2_5_sensitivity_" * s, "rows.csv"))
    for shape in ("step", "linear", "exponential")
        q = sort([r for r in p25 if r["sweep"] == s && r["shape"] == shape], by=r -> num(r, key[s]))
        lines!(ax, [num(r, key[s]) for r in q], [num(r, "dS_peak") for r in q], color=colors[shape], label=shape)
        m = sort([r for r in mc if r["shape"] == shape], by=r -> num(r, key[s]))
        scatter!(ax, [num(r, key[s]) for r in m], [num(r, "dS_peak") for r in m], color=colors[shape], markersize=7)
    end
    hlines!(ax, [0.0], color=:black, linewidth=0.5)
    j == 1 && axislegend(ax, position=:lb, labelsize=8)
end
save(joinpath(OUT, "p2_5_continuum_sweeps.png"), fig)

mp = csvrows(joinpath(OUT, "map.csv"))
fig2 = Figure(size=(1300, 720))
for (j, shape) in enumerate(("step", "linear", "exponential")), (i, (k, ttl)) in enumerate((("max_abs_dS_fixed_rho", "max |ΔS|, fixed ρ"), ("best_fit_residual", "residual, best-fit conventional")))
    q = [r for r in mp if r["shape"] == shape]
    xs = sort(unique([num(r, "rho_w_over_D") for r in q])); ys = sort(unique([num(r, "h_over_w") for r in q]))
    Z = [log10(max(num(only([r for r in q if num(r, "rho_w_over_D") == x && num(r, "h_over_w") == y]), k), 1e-10)) for x in xs, y in ys]
    ax = Axis(fig2[i, j], xlabel="log10(ρw/D)", ylabel="log10(h/w)", title="$shape: $ttl")
    hm = heatmap!(ax, log10.(xs), log10.(ys), Z, colorrange=(-8, 0))
    contour!(ax, log10.(xs), log10.(ys), Z, levels=[log10(0.01), log10(0.02)], color=:white, linewidth=2)
    j == 3 && Colorbar(fig2[i, 4], hm, label="log10")
end
Label(fig2[0, :], "Continuum (w = 2, D = 3), 5 mean decay times: white contours = 2σ single-measurement thresholds at SNR 200 (0.01) and 100 (0.02)")
save(joinpath(OUT, "p2_5_continuum_map.png"), fig2)

t2 = csvrows(joinpath(OUT, "t2_continuum.csv"))
fig3 = Figure(size=(1300, 700))
for (i, match) in enumerate(("h", "centroid")), (j, s) in enumerate(("w", match == "h" ? "h" : "c", "D", "rho"))
    pk = Dict("w" => "w_um", "h" => "h_um", "c" => "centroid_um", "D" => "D", "rho" => "rho")[s]
    ax = Axis(fig3[i, j], xscale=log10, xlabel=pk, ylabel="ΔT₂/T₂ at TE = 50", title="same $match, sweep $s: continuum (lines), MC τ=1e-2 (dots)")
    mc = csvrows(joinpath(R, "p2_5_t2_" * s * (match == "centroid" ? "_centroid" : ""), "rows.csv"))
    for shape in ("step", "linear", "exponential")
        q = sort([r for r in t2 if r["match"] == match && r["sweep"] == s && r["shape"] == shape && num(r, "TE_ms") == 50.0], by=r -> num(r, pk))
        lines!(ax, [num(r, pk) for r in q], [num(r, "rel") for r in q], color=colors[shape])
        m = sort([r for r in mc if r["shape"] == shape && num(r, "TE_ms") == 50.0], by=r -> num(r, pk))
        scatter!(ax, [num(r, pk) for r in m], [num(r, "rel") for r in m], color=colors[shape], markersize=6)
    end
    hlines!(ax, [0.0], color=:black, linewidth=0.5)
end
save(joinpath(OUT, "p2_5_continuum_t2.png"), fig3)
println("saved to $OUT")
