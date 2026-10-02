# Phase 2 exact expectations and the B4 reference curves (P2.1.5).
#
# Exact cases (tracker, "Exact expectations"): they hold in any diffusion regime.
# Everything else is measured and compared with the B4 baseline θ_relax curves, which are a reference, not the truth.

"h = 0 (no layer), R2_bulk = 0: the signal stays 1."
expected_null(times) = ones(length(times))

"Bulk R2 only: exp(-R2_bulk·t)."
expected_bulk(times, R2_bulk) = exp.(-R2_bulk .* times)

"Step profile, h = w/2, both faces: every point of the gap is in exactly one layer, so S = exp(-(ΔR2(0) + R2_bulk)·t)."
expected_half_gap(times, rate; R2_bulk=0.0) = exp.(-(rate + R2_bulk) .* times)

"Bulk on top of a layer (same seeds): S_bulk+layer = S_layer·exp(-R2_bulk·t)."
expected_additive(S_layer, times, R2_bulk) = S_layer .* exp.(-R2_bulk .* times)

"Floating-point bound for an ensemble mean after maximum(times)/τ multiplications per spin and a sum over nspins."
roundoff_bound(times, tau, nspins) = (maximum(times) / tau + nspins) * eps()

const B4_CSV = joinpath(REPO_ROOT, "research", "results", "baseline", "b4_theta_reference", "reference_curves.csv")

"""
    load_b4(; config="tau_1e-2")

B4 reference curves (bulk off) for one timestep configuration ("default", "tau_1e-2", "tau_1e-3").
Returns Dict θ => (times, mean, sigma, sem), readouts sorted by time.
"""
function load_b4(; config="tau_1e-2")
    lines = readlines(B4_CSV)
    header = split(lines[1], ",")
    col(name) = findfirst(==(name), header)
    rows = [split(l, ",") for l in lines[2:end]]
    num(r, name) = parse(Float64, r[col(name)])
    rows = filter(r -> r[col("config")] == config && num(r, "R2_bulk") == 0.0, rows)
    out = Dict{Float64, NamedTuple{(:times, :mean, :sigma, :sem), NTuple{4, Vector{Float64}}}}()
    for θ in sort(unique(num.(rows, "theta")))
        rs = sort(filter(r -> num(r, "theta") == θ, rows), by=r -> num(r, "t_ms"))
        f(name) = [num(r, name) for r in rs]
        out[θ] = (times=f("t_ms"), mean=f("mean_signal"), sigma=f("sigma"), sem=f("sem"))
    end
    return out
end

"""
    theta_for_signal(ref, S_target; t=50.0)

θ_relax of the baseline whose signal at time `t` equals `S_target`: log S interpolated linearly in log θ
between grid values (θ > 0). Returns NaN if `S_target` is outside the range spanned by the grid.
"""
function theta_for_signal(ref, S_target; t=50.0)
    θs = sort([θ for θ in keys(ref) if θ > 0])
    k = findfirst(==(t), ref[θs[1]].times)
    Ss = [ref[θ].mean[k] for θ in θs]                       # decreasing in θ
    (S_target > Ss[1] || S_target < Ss[end]) && return NaN
    j = findfirst(i -> Ss[i] >= S_target >= Ss[i + 1], 1:length(θs) - 1)
    S_target == Ss[j] && return θs[j]
    f = (log(S_target) - log(Ss[j])) / (log(Ss[j + 1]) - log(Ss[j]))
    return exp(log(θs[j]) + f * (log(θs[j + 1]) - log(θs[j])))
end

"""
    reference_at(ref, θ)

Baseline curve S_θ(t) at every readout for any θ inside the grid: the stored curve on a grid point,
otherwise log S interpolated linearly in log θ between the two neighbouring grid values.
"""
function reference_at(ref, θ)
    haskey(ref, θ) && return ref[θ].mean
    θs = sort([x for x in keys(ref) if x > 0])
    (θ < θs[1] || θ > θs[end]) && error("θ = $θ is outside the B4 grid [$(θs[1]), $(θs[end])].")
    j = searchsortedlast(θs, θ)
    f = (log(θ) - log(θs[j])) / (log(θs[j + 1]) - log(θs[j]))
    return exp.((1 - f) .* log.(ref[θs[j]].mean) .+ f .* log.(ref[θs[j + 1]].mean))
end
