# P3.4.8 outside: same-centroid comparison of the three profiles in the square array s = 3.0 (R = 1, gap 1.0, outer layer
# only, ρ = 0.02), from existing results only (no new simulation):
#   step: p3_4_8_h_sweep_outside (s = 3.0 rows; τ_b) and its τ_c thin-end pass (steps_tau16.csv, h ≤ 0.2);
#   linear, exponential: p3_4_8e_profiles_outside_<shape> (τ_b) and the τ_c thin-end pass (vs_conv_tau16.csv).
# ΔS(50) vs the conventional model of the same run (identical trajectories, so a difference of means is the paired mean).
# At the thin end (step h ≤ 0.2, linear h ≤ 0.3, exponential h ≤ 0.2) the τ_c = 6.25e-4 value is used, otherwise τ_b.
# Centroids: step h/2, linear h/3, exponential (support L = gap = 1) h(1 − (c + 1)e^(−c))/(1 − e^(−c)), c = L/h.
# ΔS is interpolated linearly in c̄ to c̄ = 0.05, 0.1, 0.15, 0.2, 0.3; spread = max pairwise |a − b|/max(|a|, |b|).
# Run: julia --project=research/baseline research/phase3/p3_4_8h_outside_centroid.jl

include(joinpath(@__DIR__, "harness.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_4_8h_outside_centroid")
readrows(f) = (ls = readlines(f); h = split(ls[1], ','); [Dict(zip(h, split(l, ','))) for l in ls[2:end]])
num(r, k) = parse(Float64, r[k])
res(name) = joinpath(REPO_ROOT, "research", "results", "phase3", name)
cexp(h; L=1.0) = (c = L / h; h * (1 - (c + 1) * exp(-c)) / -expm1(-c))

# step (s = 3.0): τ_b values from S50 differences; thin end (h ≤ 0.2) from the τ_c paired steps, accumulated from conventional
sw = filter(r -> num(r, "s") == 3.0 && num(r, "tau") == 2.5e-3, readrows(joinpath(res("p3_4_8_h_sweep_outside"), "sweep.csv")))
Sconv = num(only(r for r in sw if r["model"] == "conventional"), "S50")
step = Dict(parse(Float64, replace(r["model"], "step_h" => "")) => num(r, "S50") - Sconv for r in sw if r["model"] != "conventional")
st16 = filter(r -> num(r, "s") == 3.0, readrows(joinpath(res("p3_4_8_h_sweep_outside"), "steps_tau16.csv")))
acc = 0.0
for r in st16                                                     # conventional → 0.05 → 0.1 → 0.2, in order
    global acc += num(r, "d_tau_c")
    step[parse(Float64, replace(r["to"], "step_h" => ""))] = acc
end
data = Dict("step" => sort([(h / 2, d) for (h, d) in step]))
for (shape, cen) in (("linear", h -> h / 3), ("exponential", h -> cexp(h)))
    d = Dict{Float64, Float64}()
    for r in readrows(joinpath(res("p3_4_8e_profiles_outside_" * shape), "sweep.csv"))
        (r["model"] == "conventional" || num(r, "tau") != 2.5e-3) && continue
        d[num(r, "h")] = num(r, "dS_vs_conv")
    end
    for r in readrows(joinpath(res("p3_4_8e_profiles_outside_" * shape), "vs_conv_tau16.csv"))
        d[parse(Float64, replace(r["model"], shape * "_h" => ""))] = num(r, "dS_tau_c")
    end
    data[shape] = sort([(cen(h), v) for (h, v) in d])
end
function interp(pts, x)
    for k in 1:length(pts) - 1
        (x0, y0), (x1, y1) = pts[k], pts[k + 1]
        x0 <= x <= x1 && return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    end
    return NaN
end
rows = NamedTuple[]
for cb in (0.05, 0.1, 0.15, 0.2, 0.3)
    v = Dict(s => interp(data[s], cb) for s in ("step", "linear", "exponential"))
    ok = [s for s in keys(v) if !isnan(v[s])]
    spread = maximum((abs(v[a] - v[b]) / max(abs(v[a]), abs(v[b])) for a in ok for b in ok if a < b); init=0.0)
    push!(rows, (centroid=cb, dS_step=v["step"], dS_linear=v["linear"], dS_exponential=v["exponential"], spread=spread))
    @printf("c̄ = %.2f: step %+.3e linear %+.3e exponential %+.3e, spread %.3f\n", cb, v["step"], v["linear"], v["exponential"], spread)
end
write_csv(joinpath(OUT, "centroid_outside.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "points" => Dict(k => collect.(v) for (k, v) in data), "rows" => rows))
fig = Figure(size=(600, 420))
ax = Axis(fig[1, 1], xlabel="centroid c̄ (µm)", ylabel="ΔS(50) = S − S_conv", title="outside, square array s = 3.0, R = 1, ρ = 0.02")
for (s, col) in zip(("step", "linear", "exponential"), (:steelblue, :darkorange, :seagreen))
    scatterlines!(ax, first.(data[s]), last.(data[s]), color=col, label=s)
end
axislegend(ax, position=:lb, labelsize=9)
save(joinpath(OUT, "p3_4_8h_outside_centroid.png"), fig)
println("saved to $OUT")
