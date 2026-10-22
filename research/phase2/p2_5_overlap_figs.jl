# P2.5.9 figures from grid.csv / ddep.csv: (1) vs overlap (h/w, and mid/wall for the exponential): max |ΔS|, ΔT2/T2 with
# the SNR-200 2σ band, best-fit residual — one line per ρw/D, panels per profile; (2) D dependence: T2(D)/T2(D = 3) − 1
# for the conventional model and the layer at each overlap, ρ = 0.1 and 1.
#
# Run: julia --project=research/baseline research/phase2/p2_5_overlap_figs.jl

include(joinpath(@__DIR__, "harness.jl"))
using CairoMakie
const OUT = phase2_outdir("p2_5_overlap")
rd(f) = (L = readlines(joinpath(OUT, f)); hdr = split(L[1], ","); [Dict(zip(hdr, split(l, ","))) for l in L[2:end]])
num(d, k) = parse(Float64, d[k])
g = rd("grid.csv"); dd = rd("ddep.csv")
rwds = sort(unique([num(r, "rho_w_over_D") for r in g]))
fig = Figure(size=(1300, 700))
for (i, shape) in enumerate(("linear", "exponential")), (j, (k, lab)) in enumerate((("max_abs_dS", "max |ΔS| (fixed ρ)"), ("rel_T2", "ΔT₂/T₂ at TE* (conventional T₂)"), ("best_fit_residual", "best-fit conventional residual")))
    ax = Axis(fig[i, j], xlabel="h/w", ylabel=lab, yscale=k == "rel_T2" ? identity : log10, title="$shape: $lab")
    for rwd in rwds
        q = sort([r for r in g if r["shape"] == shape && num(r, "rho_w_over_D") == rwd], by=r -> num(r, "h_over_w"))
        y = [k == "rel_T2" ? num(r, k) : max(num(r, k), 1e-9) for r in q]
        scatterlines!(ax, [num(r, "h_over_w") for r in q], y, label="ρw/D = $rwd")
    end
    k == "max_abs_dS" || k == "best_fit_residual" ? hlines!(ax, [0.01, 0.02], color=:gray, linestyle=:dash) : nothing
    shape == "linear" && vlines!(ax, [0.5], color=:black, linestyle=:dot)        # layers touch at h/w = 0.5
    (i, j) == (1, 3) && axislegend(ax, position=:rb, labelsize=8)
end
save(joinpath(OUT, "p2_5_overlap_signal.png"), fig)

fig2 = Figure(size=(900, 380))
for (j, ρ) in enumerate((0.1, 1.0))
    ax = Axis(fig2[1, j], xlabel="h/w", ylabel="T₂(D = 0.5)/T₂(D = 3) − 1", title="D dependence of apparent T₂, ρ = $ρ")
    conv = only([r for r in dd if num(r, "rho") == ρ && r["model"] == "conventional"])
    hlines!(ax, [num(conv, "rel_change_D05_vs_D3")], color=:black, linestyle=:dash, label="conventional")
    for shape in ("linear", "exponential")
        q = sort([r for r in dd if num(r, "rho") == ρ && r["model"] == shape], by=r -> num(r, "h_over_w"))
        scatterlines!(ax, [num(r, "h_over_w") for r in q], [num(r, "rel_change_D05_vs_D3") for r in q], label=shape)
    end
    axislegend(ax, position=:rt, labelsize=8)
end
save(joinpath(OUT, "p2_5_overlap_ddep.png"), fig2)
println("saved to $OUT")
