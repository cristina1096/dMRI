# P2.6.4 — profile dependence of the timestep error: linear vs step toy scans at matching (h, ΔR2(0), D).
#
# For every scan: bias of R at τ = 0.1 ms, τ_conv at 1 % and 0.1 % (lower-bound / noise flags), ΔR2(0)·τ_conv, and the
# error of R at the default constraint τ = 0.005/ΔR2(0) (P2.6.5), interpolated linearly in log τ between grid points
# (capped at 0.1 ms, the largest τ scanned). Question answered: does the constraint ΔR2(0)·τ ≤ c, with the profile's own
# ΔR2(0) = (ρ/h)·g(0), hold for the linear profile with the same c?
#
# Run: julia --project=research/baseline research/phase2/p2_6_4_profile.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_6_4_profile")
const C_LAYER = mr.TimeSteps.DEFAULT_LAYER_SCALING
const PAIRS = [("V1_h0.1", "L_h0.1"), ("V2_rate0.5", "L_rate0.5"), ("V2_rate1", "L_rate1"), ("V2_rate2", "L_rate2"),
    ("V2_h0.4_rate1", "L_h0.4_rate1"), ("V3_D1_rate1", "L_D1_rate1"), ("V2_h0.01_rate10", "L_h0.01_rate10"), ("", "L_h0.2_rate1")]
R = joinpath(REPO_ROOT, "research", "results", "phase2")
summ(l) = JSON.parsefile(joinpath(R, "p2_6_scan_" * l, "summary.json"))
function scan(l)
    s = [split(x, ",") for x in readlines(joinpath(R, "p2_6_scan_" * l, "scan.csv"))[2:end]]
    τ = [parse(Float64, c[1]) for c in s]; Rv = [parse(Float64, c[2]) for c in s]
    o = sortperm(τ)
    return τ[o], Rv[o] ./ Rv[o[1]] .- 1
end
"Relative error of R at τ, linear in log τ between grid points (τ capped to the scanned range)."
function err_at(l, τ)
    (τs, e) = scan(l)
    τ = clamp(τ, τs[1], τs[end])
    k = clamp(searchsortedlast(τs, τ), 1, length(τs) - 1)
    f = (log(τ) - log(τs[k])) / (log(τs[k + 1]) - log(τs[k]))
    return e[k] + f * (e[k + 1] - e[k])
end
function describe(l)
    isempty(l) && return nothing
    s = summ(l)
    c(tol) = s["tau_conv"][string(tol)]
    flag(tol) = c(tol)["lower_bound"] ? "≥" : (c(tol)["noise_limited"] ? "noisy" : "")
    τc = C_LAYER / s["rate"]
    return (label=l, shape=get(s, "shape", "step"), h_um=s["h_um"], rate=s["rate"], D=s["D"], bias_0p1ms=s["bias_at_largest_tau"],
        rate_tau_conv_1pct=s["rate"] * c(0.01)["tau_conv"], flag_1pct=flag(0.01),
        rate_tau_conv_0p1pct=s["rate"] * c(0.001)["tau_conv"], flag_0p1pct=flag(0.001),
        tau_default=min(τc, 0.1), err_at_default=err_at(l, τc))
end

rows = NamedTuple[]
for (a, b) in PAIRS
    for d in (describe(a), describe(b))
        isnothing(d) && continue
        push!(rows, d)
        @printf("%-16s %-6s h=%.2f ΔR2(0)=%5.1f D=%.0f: bias(0.1 ms) %+.2e; ΔR2(0)·τ_conv 1%% %s%.3g, 0.1%% %s%.3g; error at τ = 0.005/ΔR2(0) %+.1e\n",
            d.label, d.shape, d.h_um, d.rate, d.D, d.bias_0p1ms, d.flag_1pct, d.rate_tau_conv_1pct, d.flag_0p1pct, d.rate_tau_conv_0p1pct, d.err_at_default)
    end
end
lin = [r for r in rows if r.shape == "linear"]
worst = maximum(abs(r.err_at_default) for r in lin)
@printf("linear: max |error| at the default constraint (c = %.3f) = %.1e\n", C_LAYER, worst)
write_csv(joinpath(OUT, "profile.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "c_layer" => C_LAYER, "rows" => rows,
    "linear_max_abs_error_at_default" => worst,
    "note" => "error at the default is interpolated in log τ between scan points; τ_conv flags: ≥ lower bound, noisy = noise-limited"))

using CairoMakie
fig = Figure(size=(1000, 380))
ax1 = Axis(fig[1, 1], xscale=log10, xlabel="ΔR₂(0)·τ", ylabel="R/R_ref − 1", title="timestep error vs ΔR₂(0)·τ (step: dashed, linear: solid)")
for (a, b) in PAIRS, (l, ls) in ((a, :dash), (b, :solid))
    isempty(l) && continue
    (τs, e) = scan(l)
    lines!(ax1, summ(l)["rate"] .* τs, e, linestyle=ls, label=l)
end
vlines!(ax1, [C_LAYER], color=:black, linestyle=:dot)
hlines!(ax1, [-0.01, -0.001], color=:gray, linestyle=:dash)
Legend(fig[1, 2], ax1, labelsize=8, rowgap=0)
save(joinpath(OUT, "p2_6_4_profile.png"), fig)
println("saved to $OUT")
