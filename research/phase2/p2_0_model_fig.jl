# Phase 2 record, section "Model, notation and estimators" — illustration only, no new simulation.
#
# (A) The three profiles g(u), ∫g = 1 (exponential truncated at c = w/h = 10 for h = 0.2, w = 2).
# (B) ΔR2(x) across a gap of w = 2 µm, both faces, at the same ρ = 0.01 µm/ms and h = 0.2 µm:
#     ΔR2(x) = (ρ/h)[g(x/h) + g((w − x)/h)]. The conventional model puts the same ρ on the walls themselves.
# (C) Jensen gap J = ln S(t) + 2ρt/w at fixed ρ = 0.02 µm/ms (E[X] = 1 at t = 50 ms), w = 2, D = 3, from the exact
#     continuum benchmark (continuum.jl; no timestep bias, see Review C1), vs h for the three profiles and for the
#     conventional surface model (the h → 0 limit). With X = −ln M⊥ of one spin, E[X] = 2ρt/w (mean-exponent identity,
#     P2.4.5) and S = E[e^(−X)], so J = ln E[e^(−X)] + E[X] ≥ 0 (Jensen), ≈ Var[X]/2 to second order; J = 0 for a
#     uniform rate. A smaller J means a less variable exponent and so a lower signal at the same ρ.
#
# Run: julia --project=research/baseline research/phase2/p2_0_model_fig.jl

include(joinpath(@__DIR__, "..", "baseline", "common.jl"))
using CairoMakie
include(joinpath(@__DIR__, "continuum.jl"))

const OUT = mkpath(joinpath(REPO_ROOT, "research", "results", "phase2", "p2_0_model"))
const W = 2.0
const T = 50.0

g(shape, u, c) = shape == "step" ? (0 <= u < 1 ? 1.0 : 0.0) :
                 shape == "linear" ? (0 <= u <= 1 ? 2 * (1 - u) : 0.0) :
                 (0 <= u <= c ? exp(-u) / (1 - exp(-c)) : 0.0)
centroid(shape, h) = shape == "step" ? h / 2 : shape == "linear" ? h / 3 :
                     (c = W / h; h * (1 - (c + 1) * exp(-c)) / (1 - exp(-c)))

colors = Dict("step" => :steelblue, "linear" => :darkorange, "exponential" => :seagreen)
shapes = ("step", "linear", "exponential")
fig = Figure(size=(1350, 400))

ax1 = Axis(fig[1, 1], xlabel="u = d/h", ylabel="g(u)", title="(A) profiles, ∫g du = 1")
us = range(0, 3, length=600)
for s in shapes
    lines!(ax1, us, [g(s, u, 10.0) for u in us], color=colors[s], label=s)
end
axislegend(ax1, labelsize=9)

ax2 = Axis(fig[1, 2], xlabel="x (µm), walls at 0 and w = 2", ylabel="ΔR₂(x) (ms⁻¹)",
    title="(B) both faces, same ρ = 0.01 µm/ms, h = 0.2 µm")
h, ρ = 0.2, 0.01
xs = range(0, W, length=2000)
for s in shapes
    lines!(ax2, xs, [ρ / h * (g(s, x / h, W / h) + g(s, (W - x) / h, W / h)) for x in xs], color=colors[s], label=s)
end
hlines!(ax2, [2ρ / W], color=:gray, linestyle=:dash, label="uniform 2ρ/w (same ρ)")
axislegend(ax2, labelsize=9, position=:ct)

ax3 = Axis(fig[1, 3], xlabel="h/w", ylabel="J = ln S(50) + 2ρ·50/w",
    title="(C) Jensen gap at fixed ρ = 0.02 (continuum)", xscale=log10, yscale=log10)
rows = NamedTuple[]
ρC, D = 0.02, 3.0
Jof(S) = log(S) + 2ρC * T / W
Jconv = Jof(continuum_signal(w=W, D=D, times=[T], surface=ρC)[1])
for s in shapes
    hs = s == "step" ? [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9] : [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 1.9]
    Js = [Jof(continuum_signal(w=W, D=D, times=[T], layer=(s, ρC, hh))[1]) for hh in hs]
    append!(rows, [(model=s, h_um=hh, centroid_um=centroid(s, hh), jensen_gap=J) for (hh, J) in zip(hs, Js)])
    keep = Js .> 1e-12          # log axis; the uniform limit has J = 0 up to roundoff
    scatterlines!(ax3, hs[keep] ./ W, Js[keep], color=colors[s], label=s)
end
push!(rows, (model="conventional", h_um=0.0, centroid_um=0.0, jensen_gap=Jconv))
hlines!(ax3, [Jconv], color=:black, linestyle=:dash, label="conventional (surface)")
axislegend(ax3, labelsize=9, position=:lb)
write_csv(joinpath(OUT, "jensen_gap.csv"), rows)
save(joinpath(OUT, "p2_0_model.png"), fig)
for r in rows
    @printf("%-12s h = %.3f: c̄ = %.4f, J = %.3e\n", r.model, r.h_um, r.centroid_um, r.jensen_gap)
end
println("saved to $OUT")
