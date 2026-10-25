# P3.4.8 (inside): attenuation vs h at fixed ρ inside single cylinders, from the radial continuum (validated in P3.4.6).
# R = 0.5, 1, 2; D = 3; ρ = 0.02 (so E[−ln M] = ρ(2/R)t is the same for every h); step; h/R = 0.005 … 1.
# Quantities at readouts 0:5:50:
#   S(50); Jensen gap J(50) = ln S(50) + ρ(2/R)·50 ≥ 0 (Phase 2 record, M3: S = e^(−E[X]) e^J);
#   ΔS(50) = S − S_conv (conventional: Robin, same ρ); curvature c of ln S = a + bt + ct² (least squares);
#   grid error = |S_N=2000 − S_N=4000| (layer), |S_N=1000 − S_N=2000| (conventional; Robin on the face since the
#   final review, I3 — the first version put the sink on the ring centre, a first-order error of ≈ 5e-7).
# Monotonicity: the sign of S(50; h_{k+1}) − S(50; h_k) is recorded for every neighbouring pair; a step counts as resolved
# when |difference| > 10 × grid error. h → 0 approaches the conventional model; h = R is uniform (J = 0).
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_8_h_sweep_inside.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_4_8_h_sweep_inside")
const RHO = 0.02
const HR = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
curv(ts, S) = (X = hcat(ones(length(ts)), ts, ts .^ 2); (X \ log.(S))[3])

rows = NamedTuple[]
steps = NamedTuple[]
conv = Dict{Float64, Float64}()
conv_grid = Dict{Float64, Float64}()
for R in (0.5, 1.0, 2.0)
    E50 = RHO * (2 / R) * 50
    Sc = radial_signal(R=R, D=D_REF, times=TIMES_REF, surface=RHO)
    conv[R] = Sc[end]
    conv_grid[R] = abs(Sc[end] - radial_signal(R=R, D=D_REF, times=[50.0], surface=RHO, N=1000)[1])   # Robin on the face (I3)
    prev = nothing
    for hr in HR
        h = hr * R
        S = radial_signal(R=R, D=D_REF, times=TIMES_REF, layer=("step", RHO, h))
        S4 = radial_signal(R=R, D=D_REF, times=[50.0], layer=("step", RHO, h), N=4000)[1]
        row = (R=R, h=h, h_over_R=hr, S50=S[end], J50=log(S[end]) + E50, dS50_vs_conv=S[end] - Sc[end],
               curvature=curv(TIMES_REF, S), grid_err=abs(S[end] - S4))
        push!(rows, row)
        if !isnothing(prev)
            d = row.S50 - prev.S50
            push!(steps, (R=R, h_from=prev.h, h_to=h, dS50=d, resolved=abs(d) > 10 * max(row.grid_err, prev.grid_err), sign=sign(d)))
        end
        prev = row
        @printf("R = %.1f h/R = %.3f: S(50) %.6f  J %.3e  ΔS vs conv %+.2e  curvature %+.1e  grid %.1e\n",
            R, hr, row.S50, row.J50, row.dS50_vs_conv, row.curvature, row.grid_err)
    end
    @printf("R = %.1f conventional: S(50) %.9f  J %.4e  grid (N 1000 vs 2000) %.1e\n", R, Sc[end], log(Sc[end]) + E50, conv_grid[R])
end
write_csv(joinpath(OUT, "sweep.csv"), rows)
write_csv(joinpath(OUT, "steps.csv"), steps)
signs = unique(s.sign for s in steps if s.resolved)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rho" => RHO, "D" => D_REF,
    "conventional_S50" => conv, "conventional_grid_err" => conv_grid, "n_steps" => length(steps), "n_resolved" => count(s -> s.resolved, steps),
    "resolved_signs" => signs, "monotone" => length(signs) == 1))
@printf("neighbour steps resolved: %d of %d; signs of resolved steps: %s\n", count(s -> s.resolved, steps), length(steps), signs)

fig = Figure(size=(1100, 400))
ax1 = Axis(fig[1, 1], xlabel="h/R", ylabel="J(50) = ln S + ρ(S/V)t", xscale=log10, title="Jensen gap (fixed ρ = 0.02)")
ax2 = Axis(fig[1, 2], xlabel="h/R", ylabel="ΔS(50) = S − S_conv", xscale=log10, title="difference from the conventional model")
for (R, col) in zip((0.5, 1.0, 2.0), (:steelblue, :darkorange, :seagreen))
    q = filter(r -> r.R == R, rows)
    scatterlines!(ax1, [r.h_over_R for r in q], [r.J50 for r in q], color=col, label="R = $R")
    hlines!(ax1, [log(conv[R]) + RHO * (2 / R) * 50], color=col, linestyle=:dash)
    scatterlines!(ax2, [r.h_over_R for r in q], [r.dS50_vs_conv for r in q], color=col, label="R = $R")
end
hlines!(ax2, [0.0], color=:black)
axislegend(ax1, labelsize=9)
save(joinpath(OUT, "p3_4_8_h_sweep_inside.png"), fig)
println("saved to $OUT")
