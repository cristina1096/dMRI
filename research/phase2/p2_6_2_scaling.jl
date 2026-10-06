# P2.6.2 V2 (h² scaling), the tracker's decision point (rate dependence), P2.6.3 V3 (1/D scaling).
#
# Reads the toy scans (p2_6_scan.jl) and tests τ_conv = c_h·h²/D at tolerances 1 % and 0.1 %.
# A τ_conv that is only a lower bound (no τ ≤ 0.1 ms fails) or noise-limited does not enter a fit; it is listed.
# V2 pass: τ_conv·D/h² equal within ±30 % across h = 0.01, 0.05, 0.1, 0.2 (rate 0.1, D = 3).
# V3 pass: τ_conv·D equal within ±30 % across D = 1, 3, 6 at h = 0.1.
# Decision point (planned when the probes contradicted the h² premise; ledger ruling): τ_conv against ΔR2(0) and
# against ρ = ΔR2(0)·h (h = 0.1 at ΔR2(0) = 0.1, 0.5, 1, 2; h = 0.01 at ΔR2(0) = 10 has the same ρ as h = 0.1 at 1).
# The bias at the largest τ (0.1 ms) is reported for every scan as the size of the effect, independent of any tolerance.
#
# Run: julia --project=research/baseline research/phase2/p2_6_2_scaling.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "timestep.jl"))

const OUT = phase2_outdir("p2_6_2_scaling")
const LABELS = ["V2_h0.01", "V2_h0.05", "V1_h0.1", "V2_h0.2", "V2_rate0.5", "V2_rate1", "V2_rate2", "V2_h0.01_rate10", "V3_D1", "V3_D6", "V2_h0.4_rate1", "V3_D1_rate1"]
load(label) = JSON.parsefile(joinpath(REPO_ROOT, "research", "results", "phase2", "p2_6_scan_" * label, "summary.json"))
scans = Dict(l => load(l) for l in LABELS)
conv(s, tol) = s["tau_conv"][string(tol)]
usable(s, tol) = !conv(s, tol)["lower_bound"] && !conv(s, tol)["noise_limited"]
tc(s, tol) = conv(s, tol)["tau_conv"]

function ratio_test(labels, f, tol)
    ok = [l for l in labels if usable(scans[l], tol)]
    length(ok) < 2 && return (n=length(ok), values=Float64[], max_dev=nothing, pass=nothing)
    v = [f(scans[l]) * tc(scans[l], tol) for l in ok]
    dev = maximum(abs.(v ./ mean(v) .- 1))
    return (n=length(ok), labels=ok, values=v, max_dev=dev, pass=dev <= 0.3)
end

rows = NamedTuple[]
for l in LABELS
    s = scans[l]
    push!(rows, (label=l, h_um=s["h_um"], rate=s["rate"], rho=s["rate"] * s["h_um"], D=s["D"],
        R_ref=s["R_ref"], bias_at_0p1ms=s["bias_at_largest_tau"],
        tau_conv_1pct=tc(s, 0.01), lower_bound_1pct=conv(s, 0.01)["lower_bound"], noise_1pct=conv(s, 0.01)["noise_limited"],
        tau_conv_0p1pct=tc(s, 0.001), lower_bound_0p1pct=conv(s, 0.001)["lower_bound"], noise_0p1pct=conv(s, 0.001)["noise_limited"]))
    @printf("%-16s h=%.2f ΔR2(0)=%5.2f ρ=%.3f D=%.0f: bias(0.1 ms) %+.2e; τ_conv 1%% %.2e%s; 0.1%% %.2e%s\n", l, s["h_um"], s["rate"],
        s["rate"] * s["h_um"], s["D"], s["bias_at_largest_tau"],
        tc(s, 0.01), usable(s, 0.01) ? "" : (conv(s, 0.01)["lower_bound"] ? " (≥, lower bound)" : " (noise-limited)"),
        tc(s, 0.001), usable(s, 0.001) ? "" : (conv(s, 0.001)["lower_bound"] ? " (≥, lower bound)" : " (noise-limited)"))
end

h_set = ["V2_h0.01", "V2_h0.05", "V1_h0.1", "V2_h0.2"]
D_set = ["V3_D1", "V1_h0.1", "V3_D6"]
rate_set = ["V1_h0.1", "V2_rate0.5", "V2_rate1", "V2_rate2"]
# coverage of c = ΔR2(0)·τ_conv away from h = 0.1, D = 3 (review I3): all at ΔR2(0) ≥ 1
cover_set = ["V2_rate1", "V2_h0.4_rate1", "V3_D1_rate1", "V2_h0.01_rate10"]
tests = Dict{String, Any}()
for tol in (0.01, 0.001)
    t = Dict(
        "V2_tau_D_over_h2" => ratio_test(h_set, s -> s["D"] / s["h_um"]^2, tol),
        "V3_tau_D" => ratio_test(D_set, s -> s["D"], tol),
        "rate_tau_times_rate" => ratio_test(rate_set, s -> s["rate"], tol),
        "coverage_tau_times_rate" => ratio_test(cover_set, s -> s["rate"], tol))
    tests[string(tol)] = t
    for (k, r) in sort(collect(t), by=first)
        @printf("tol %.1f %%, %-20s: %d usable scans, values %s, max deviation %s → %s\n", 100tol, k, r.n,
            round.(r.values, sigdigits=3), isnothing(r.max_dev) ? "–" : @sprintf("%.0f %%", 100r.max_dev),
            isnothing(r.pass) ? "not testable" : (r.pass ? "PASS" : "FAIL"))
    end
end
same_rho = (scans["V2_rate1"]["bias_at_largest_tau"], scans["V2_h0.01_rate10"]["bias_at_largest_tau"])
@printf("same ρ = 0.1 µm/ms: bias(0.1 ms) h=0.1/ΔR2(0)=1: %+.2e, h=0.01/ΔR2(0)=10: %+.2e\n", same_rho...)

write_csv(joinpath(OUT, "scaling.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rows" => rows, "tests" => tests,
    "same_rho_bias" => same_rho,
    "note" => "τ_conv that is a lower bound or noise-limited is not used in any ratio test; no c_h is adopted here (ledger ruling, Task 2)"))

using CairoMakie
fig = Figure(size=(1000, 380))
ax1 = Axis(fig[1, 1], xscale=log10, xlabel="τ (ms)", ylabel="R/R_ref − 1", title="all toy scans")
for l in LABELS
    s = [split(x, ",") for x in readlines(joinpath(REPO_ROOT, "research", "results", "phase2", "p2_6_scan_" * l, "scan.csv"))[2:end]]
    τ = [parse(Float64, c[1]) for c in s]; Rv = [parse(Float64, c[2]) for c in s]
    scatterlines!(ax1, τ, Rv ./ Rv[argmin(τ)] .- 1, label=l, markersize=5)
end
hlines!(ax1, [-0.01, 0.01], color=:gray, linestyle=:dash); hlines!(ax1, [-0.001, 0.001], color=:gray, linestyle=:dot)
Legend(fig[1, 2], ax1, labelsize=8, rowgap=0)
ax2 = Axis(fig[1, 3], xscale=log10, xlabel="ΔR₂(0) (1/ms)", ylabel="|bias at τ = 0.1 ms|", yscale=log10, title="size of the τ effect vs layer rate (h = 0.1)")
b = [abs(scans[l]["bias_at_largest_tau"]) for l in rate_set]
scatterlines!(ax2, [scans[l]["rate"] for l in rate_set], max.(b, 1e-6))
save(joinpath(OUT, "scaling.png"), fig)
println("saved to $OUT")
