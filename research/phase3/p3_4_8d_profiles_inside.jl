# P3.4.8 inside for all profiles, from the radial continuum (validated for step in P3.4.6, linear/exponential in Task 7).
# Fixed ρ = 0.02, D = 3, t = 50 ms, R = 0.5, 1, 2.
# (1) h sweep: step and linear h/R = 0.005 … 1; exponential h/R = 0.005 … 5 (support R). S(50), J(50) = ln S + ρ(2/R)·50,
#     ΔS(50) vs the conventional model; monotonicity: sign of neighbouring steps (resolved if > 10 × grid error).
# (2) Same centroid (Phase 2 P2.5.6 found profiles nearly coincide at the same centroid on walls): for c̄/R = 0.02, 0.05,
#     0.1, 0.2, choose h per profile with that centroid (step h = 2c̄, linear h = 3c̄, exponential by bisection with
#     support R), and report ΔS(50) per profile and the spread max|a − b|/max(|a|, |b|) over profile pairs; also the same
#     comparison at the same h.
# Run: julia --project=research/baseline -t 4 research/phase3/p3_4_8d_profiles_inside.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_4_8d_profiles_inside")
const RHO = 0.02
const T = 50.0
HR = Dict("step" => [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0],
          "linear" => [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0],
          "exponential" => [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0])
centroid_of(shape, h, R) = shape == "step" ? h / 2 : shape == "linear" ? h / 3 :
                           (c = R / h; h * (1 - (c + 1) * exp(-c)) / -expm1(-c))
function h_for_centroid(shape, cbar, R)
    shape == "step" && return 2cbar
    shape == "linear" && return 3cbar
    lo, hi = 1e-6, 100R                     # exponential centroid increases monotonically with h, from 0 to R/2
    for _ in 1:200
        mid = sqrt(lo * hi)
        centroid_of("exponential", mid, R) < cbar ? (lo = mid) : (hi = mid)
    end
    return sqrt(lo * hi)
end
S50(shape, h, R) = radial_signal(R=R, D=D_REF, times=[T], layer=(shape, RHO, h))[1]

rows = NamedTuple[]
steps = NamedTuple[]
for R in (0.5, 1.0, 2.0)
    Sc = radial_signal(R=R, D=D_REF, times=[T], surface=RHO)[1]
    E = RHO * (2 / R) * T
    for shape in ("step", "linear", "exponential")
        prev = nothing
        for hr in HR[shape]
            h = hr * R
            S = S50(shape, h, R)
            ge = abs(S - radial_signal(R=R, D=D_REF, times=[T], layer=(shape, RHO, h), N=4000)[1])
            row = (R=R, shape=shape, h=h, h_over_R=hr, centroid_over_R=centroid_of(shape, h, R) / R, S50=S, J50=log(S) + E,
                   dS50=S - Sc, grid_err=ge)
            push!(rows, row)
            if !isnothing(prev)
                d = S - prev.S50
                push!(steps, (R=R, shape=shape, h_from=prev.h, h_to=h, dS=d, resolved=abs(d) > 10 * max(ge, prev.grid_err), sign=sign(d)))
            end
            prev = row
        end
    end
end
cen = NamedTuple[]
for R in (0.5, 1.0, 2.0), cb in (0.02, 0.05, 0.1, 0.2)
    Sc = radial_signal(R=R, D=D_REF, times=[T], surface=RHO)[1]
    vals = Dict(sh => S50(sh, h_for_centroid(sh, cb * R, R), R) - Sc for sh in ("step", "linear", "exponential"))
    hvals = Dict(sh => S50(sh, 2cb * R, R) - Sc for sh in ("step", "linear", "exponential"))     # same h = 2c̄ (step's)
    spread(v) = maximum(abs(v[a] - v[b]) / max(abs(v[a]), abs(v[b])) for a in keys(v) for b in keys(v) if a < b)
    push!(cen, (R=R, centroid_over_R=cb, dS_step=vals["step"], dS_linear=vals["linear"], dS_exponential=vals["exponential"],
        spread_same_centroid=spread(vals), spread_same_h=spread(hvals)))
    @printf("R = %.1f c̄/R = %.2f: ΔS step %+.3e linear %+.3e exponential %+.3e; spread same c̄ %.3f, same h %.3f\n",
        R, cb, vals["step"], vals["linear"], vals["exponential"], spread(vals), spread(hvals))
end
signs = unique(s.sign for s in steps if s.resolved)
@printf("h steps resolved: %d of %d; signs: %s\n", count(s -> s.resolved, steps), length(steps), signs)
write_csv(joinpath(OUT, "sweep.csv"), rows)
write_csv(joinpath(OUT, "steps.csv"), steps)
write_csv(joinpath(OUT, "centroid.csv"), cen)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rho" => RHO, "D" => D_REF, "t" => T,
    "n_steps" => length(steps), "n_resolved" => count(s -> s.resolved, steps), "resolved_signs" => signs,
    "max_spread_same_centroid" => maximum(c.spread_same_centroid for c in cen), "max_spread_same_h" => maximum(c.spread_same_h for c in cen)))

fig = Figure(size=(1100, 400))
ax1 = Axis(fig[1, 1], xlabel="h/R", ylabel="ΔS(50) = S − S_conv", xscale=log10, title="R = 1, fixed ρ = 0.02: same h")
ax2 = Axis(fig[1, 2], xlabel="c̄/R", ylabel="ΔS(50) = S − S_conv", xscale=log10, title="R = 1: same centroid")
for (sh, col) in zip(("step", "linear", "exponential"), (:steelblue, :darkorange, :seagreen))
    q = filter(r -> r.R == 1.0 && r.shape == sh, rows)
    scatterlines!(ax1, [r.h_over_R for r in q], [r.dS50 for r in q], color=col, label=sh)
    scatterlines!(ax2, [r.centroid_over_R for r in q], [r.dS50 for r in q], color=col, label=sh)
end
axislegend(ax1, position=:lb, labelsize=9)
save(joinpath(OUT, "p3_4_8d_profiles_inside.png"), fig)
println("saved to $OUT")
