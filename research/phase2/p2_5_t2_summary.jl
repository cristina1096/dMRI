# P2.5.6/P2.5.7 summary — apparent T2 vs microstructure parameters (Fig. 9 analogue) and vs echo time (Fig. 8A
# analogue, (4) below) of Cottaar et al. 2026, walls.
#
# (1) Top row: T2_app(TE = 50 ms) vs w, h, D, ρ for the conventional model and the three profiles (step, linear,
#     exponential). (2) Bottom row: ΔT2/T2 (layer − conventional) ± paired SEM, with the two-point measurement
#     noise bands ±2σ_T2/T2 at SNR 100 and 200. (3) TE dependence: ΔT2/T2 at TE = 20, 50, 100 for the w sweep.
# Distinct configurations are counted once (the reference point appears in every sweep).
#
# Run: julia --project=research/baseline research/phase2/p2_5_t2_summary.jl

include(joinpath(@__DIR__, "harness.jl"))

const MATCH = get(ENV, "MATCH", "h")                       # "h" (same h) or "centroid" (same centroid): separate figures
const SFX = MATCH == "centroid" ? "_centroid" : ""
const OUT = phase2_outdir("p2_5_t2_summary")
const SWEEPS = MATCH == "centroid" ? ("w", "c", "D", "rho") : ("w", "h", "D", "rho")
const PARAM = Dict("w" => "w_um", "h" => "h_um", "c" => "centroid_um", "D" => "D", "rho" => "rho")
const LABEL = Dict("w" => "wall spacing w (µm)", "h" => "layer thickness h (µm)", "c" => "centroid c̄ (µm)", "D" => "D (µm²/ms)", "rho" => "ρ (µm/ms)")
R = joinpath(REPO_ROOT, "research", "results", "phase2")
function load(s)
    L = readlines(joinpath(R, "p2_5_t2_" * s * SFX, "rows.csv"))
    hdr = split(L[1], ",")
    return [Dict(zip(hdr, split(l, ","))) for l in L[2:end]]
end
num(d, k) = parse(Float64, d[k])
rows = Dict(s => load(s) for s in SWEEPS)

tab = NamedTuple[]
for s in SWEEPS, r in rows[s]
    push!(tab, (sweep=s, shape=r["shape"], param=num(r, PARAM[s]), TE_ms=num(r, "TE_ms"), T2_layer=num(r, "T2_layer"),
        T2_conv=num(r, "T2_conv"), rel=num(r, "rel"), rel_sem=num(r, "rel_sem"), noise200=num(r, "noise_rel_snr200"),
        noise100=num(r, "noise_rel_snr100"), level=r["level"], key=(r["shape"], num(r, "h_um"), num(r, "w_um"), num(r, "D"), num(r, "rho"), num(r, "TE_ms"))))
end
write_csv(joinpath(OUT, "t2_table$(SFX).csv"), [Base.structdiff(t, NamedTuple{(:key,)}) for t in tab])
uniq = Dict(t.key => t for t in tab)
te50 = [t for t in values(uniq) if t.TE_ms == 50.0]
counts = Dict(l => count(t -> t.level == l, te50) for l in ("SNR 100", "SNR 200", "MC only", "not resolved"))
worst = te50[argmax([abs(t.rel) for t in te50])]
@printf("TE = 50: %d distinct configurations; levels %s; largest |ΔT2/T2| = %.2e (%s, %s sweep)\n", length(te50), counts, abs(worst.rel), worst.shape, worst.sweep)
write_json(joinpath(OUT, "summary$(SFX).json"), Dict("provenance" => provenance(), "match" => MATCH, "counts_TE50" => counts, "n_distinct_TE50" => length(te50),
    "largest_rel_TE50" => Dict("rel" => worst.rel, "shape" => worst.shape, "sweep" => worst.sweep, "param" => worst.param)))

using CairoMakie
colors = Dict("step" => :steelblue, "linear" => :darkorange, "exponential" => :seagreen)
fig = Figure(size=(1400, 900))
for (j, s) in enumerate(SWEEPS)
    sel = [t for t in tab if t.sweep == s && t.TE_ms == 50.0]
    ax = Axis(fig[1, j], xscale=log10, xlabel=LABEL[s], ylabel="apparent T₂ at TE = 50 ms (ms)", title=("A", "B", "C", "D")[j] * ": sweep " * s)
    conv = sort(unique([(t.param, t.T2_conv) for t in sel]))
    scatterlines!(ax, first.(conv), last.(conv), color=:black, linestyle=:dash, label="conventional")
    for shape in ("step", "linear", "exponential")
        q = sort([t for t in sel if t.shape == shape], by=t -> t.param)
        isempty(q) || scatterlines!(ax, [t.param for t in q], [t.T2_layer for t in q], color=colors[shape], label=shape)
    end
    j == 1 && axislegend(ax, position=:lt, labelsize=9)
    bx = Axis(fig[2, j], xscale=log10, xlabel=LABEL[s], ylabel="ΔT₂/T₂ (layer − conventional)")
    for shape in ("step", "linear", "exponential")
        q = sort([t for t in sel if t.shape == shape], by=t -> t.param)
        isempty(q) && continue
        errorbars!(bx, [t.param for t in q], [t.rel for t in q], 3 .* [t.rel_sem for t in q], color=colors[shape])
        scatterlines!(bx, [t.param for t in q], [t.rel for t in q], color=colors[shape])
    end
    cq = sort(unique([(t.param, t.noise200, t.noise100) for t in sel]))
    band!(bx, first.(cq), -2 .* getindex.(cq, 2), 2 .* getindex.(cq, 2), color=(:gray, 0.15))
    lines!(bx, first.(cq), 2 .* getindex.(cq, 3), color=:gray, linestyle=:dash); lines!(bx, first.(cq), -2 .* getindex.(cq, 3), color=:gray, linestyle=:dash)
    hlines!(bx, [0.0], color=:black, linewidth=0.5)
    # zoomed: the profile differences are far inside the noise bands, so show them on their own scale (±3·SEM)
    zx = Axis(fig[3, j], xscale=log10, xlabel=LABEL[s], ylabel="ΔT₂/T₂, zoomed (±3·SEM)")
    for shape in ("step", "linear", "exponential")
        q = sort([t for t in sel if t.shape == shape], by=t -> t.param)
        isempty(q) && continue
        errorbars!(zx, [t.param for t in q], [t.rel for t in q], 3 .* [t.rel_sem for t in q], color=colors[shape])
        scatterlines!(zx, [t.param for t in q], [t.rel for t in q], color=colors[shape])
    end
    hlines!(zx, [0.0], color=:black, linewidth=0.5)
end
Label(fig[0, :], "Apparent T₂ vs microstructure at fixed ρ, profiles matched by $(MATCH == "centroid" ? "centroid" : "h") (cf. Cottaar et al. 2026, Fig. 9); shaded: ±2σ T₂ noise at SNR 200, dashed: SNR 100")
save(joinpath(OUT, "p2_5_t2_fig9$(SFX).png"), fig)

fig2 = Figure(size=(560, 380))
ax = Axis(fig2[1, 1], xscale=log10, xlabel="wall spacing w (µm)", ylabel="ΔT₂/T₂", title="w sweep: TE dependence")
for shape in ("step", "linear", "exponential"), (TE, ls) in ((20.0, :dot), (50.0, :solid), (100.0, :dash))
    q = sort([t for t in tab if t.sweep == "w" && t.shape == shape && t.TE_ms == TE], by=t -> t.param)
    lines!(ax, [t.param for t in q], [t.rel for t in q], color=colors[shape], linestyle=ls, label="$shape TE=$(Int(TE))")
end
Legend(fig2[1, 2], ax, labelsize=8, rowgap=0)
save(joinpath(OUT, "p2_5_t2_te$(SFX).png"), fig2)

# Fig. 8A analogue: T2_app vs TE (5–200 ms) per configuration (w = 1, 2, 5 at h = 0.2; h = 0.05, 0.5 at w = 2).
te = load("te")
# columns: one per (w, matching value), where the matching value is h (same h) or the centroid (same centroid)
mkey = MATCH == "centroid" ? "cbar" : "h_um"
mref = MATCH == "centroid" ? 0.1 : 0.2
pts = sort(unique([(num(r, "w_um"), round(num(r, mkey), sigdigits=6)) for r in te]), by=p -> (p[2] == mref ? 0 : 1, p[1], p[2]))
fig3 = Figure(size=(300 * length(pts), 900))
for (j, (w, m)) in enumerate(pts)
    sel = [r for r in te if num(r, "w_um") == w && round(num(r, mkey), sigdigits=6) == m]
    ax = Axis(fig3[1, j], xlabel="TE (ms)", ylabel="apparent T₂ (ms)", title=MATCH == "centroid" ? "w = $w µm, c̄ = $m µm" : "w = $w µm, h = $m µm")
    conv = sort(unique([(num(r, "TE_ms"), num(r, "T2_conv"), num(r, "T2_conv_sem")) for r in sel]))
    band!(ax, first.(conv), getindex.(conv, 2) .- 2 .* getindex.(conv, 3), getindex.(conv, 2) .+ 2 .* getindex.(conv, 3), color=(:black, 0.1))
    lines!(ax, first.(conv), getindex.(conv, 2), color=:black, linestyle=:dash, label="conventional")
    for shape in ("step", "linear", "exponential")
        q = sort([r for r in sel if r["shape"] == shape], by=r -> num(r, "TE_ms"))
        lines!(ax, [num(r, "TE_ms") for r in q], [num(r, "T2_layer") for r in q], color=colors[shape], label=shape)
    end
    j == 1 && axislegend(ax, position=:rb, labelsize=8)
    bx = Axis(fig3[2, j], xlabel="TE (ms)", ylabel="ΔT₂/T₂ (layer − conventional)")
    for shape in ("step", "linear", "exponential")
        q = sort([r for r in sel if r["shape"] == shape], by=r -> num(r, "TE_ms"))
        x = [num(r, "TE_ms") for r in q]; y = [num(r, "rel") for r in q]; e3 = 3 .* [num(r, "rel_sem") for r in q]
        band!(bx, x, y .- e3, y .+ e3, color=(colors[shape], 0.2)); lines!(bx, x, y, color=colors[shape])
    end
    n = sort(unique([(num(r, "TE_ms"), num(r, "noise_rel_snr200")) for r in sel]))
    lines!(bx, first.(n), 2 .* last.(n), color=:gray, linestyle=:dash); lines!(bx, first.(n), -2 .* last.(n), color=:gray, linestyle=:dash)
    hlines!(bx, [0.0], color=:black, linewidth=0.5)
    zx = Axis(fig3[3, j], xlabel="TE (ms)", ylabel="ΔT₂/T₂, zoomed (±3·SEM)")
    for shape in ("step", "linear", "exponential")
        q = sort([r for r in sel if r["shape"] == shape], by=r -> num(r, "TE_ms"))
        x = [num(r, "TE_ms") for r in q]; y = [num(r, "rel") for r in q]; e3 = 3 .* [num(r, "rel_sem") for r in q]
        band!(zx, x, y .- e3, y .+ e3, color=(colors[shape], 0.2)); lines!(zx, x, y, color=colors[shape])
    end
    hlines!(zx, [0.0], color=:black, linewidth=0.5)
end
Label(fig3[0, :], "Apparent T₂ vs echo time at fixed ρ, profiles matched by $(MATCH == "centroid" ? "centroid" : "h") (cf. Cottaar et al. 2026, Fig. 8A; no off-resonance). Without surface relaxation T₂ = 2000 ms. Bottom: ±3·SEM bands, dashed ±2σ T₂ noise at SNR 200")
save(joinpath(OUT, "p2_5_t2_fig8a$(SFX).png"), fig3)

# Direct comparison of the two matchings (written on the centroid run, once both result sets exist):
# ΔT2/T2 at TE = 50 vs w, per shape, same h (solid) vs same centroid (dashed).
if MATCH == "centroid" && isfile(joinpath(R, "p2_5_t2_w", "rows.csv"))
    loadh(s) = (L = readlines(joinpath(R, "p2_5_t2_" * s, "rows.csv")); hdr = split(L[1], ","); [Dict(zip(hdr, split(l, ","))) for l in L[2:end]])
    fig4 = Figure(size=(620, 400))
    ax = Axis(fig4[1, 1], xscale=log10, xlabel="wall spacing w (µm)", ylabel="ΔT₂/T₂ at TE = 50 ms", title="same h (solid) vs same centroid (dashed)")
    for shape in ("step", "linear", "exponential"), (src, ls) in ((loadh("w"), :solid), (load("w"), :dash))
        q = sort([r for r in src if r["shape"] == shape && num(r, "TE_ms") == 50.0], by=r -> num(r, "w_um"))
        isempty(q) || scatterlines!(ax, [num(r, "w_um") for r in q], [num(r, "rel") for r in q], color=colors[shape], linestyle=ls, label="$shape $(ls == :solid ? "same h" : "same c̄")")
    end
    hlines!(ax, [0.0], color=:black, linewidth=0.5)
    Legend(fig4[1, 2], ax, labelsize=8, rowgap=0)
    save(joinpath(OUT, "p2_5_t2_match_compare.png"), fig4)
end
println("saved to $OUT")
