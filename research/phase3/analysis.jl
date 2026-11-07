# Phase 3 analysis helpers (P3.6): power-law fits and thresholds on scanned error curves.

"""
    fit_power(x, y; w=nothing)

Least squares of log y = log a + p log x (weights `w` on the log points); returns (a, p, se_p), with se_p the standard
error of the slope from the residual variance (n − 2 degrees of freedom; NaN for n = 2). Requires x, y > 0, at least
two points and two distinct x (ArgumentError otherwise).
"""
function fit_power(x, y; w=nothing)
    (all(>(0), x) && all(>(0), y)) || throw(ArgumentError("fit_power needs x, y > 0"))
    length(x) >= 2 || throw(ArgumentError("fit_power needs at least 2 points"))
    lx, ly = log.(x), log.(y)
    W = isnothing(w) ? ones(length(x)) : collect(Float64, w)
    sw = sum(W)
    mx, my = sum(W .* lx) / sw, sum(W .* ly) / sw
    sxx = sum(W .* (lx .- mx) .^ 2)
    sxx > 0 || throw(ArgumentError("fit_power needs at least two distinct x"))
    p = sum(W .* (lx .- mx) .* (ly .- my)) / sxx
    a = exp(my - p * mx)
    n = length(x)
    res = ly .- (log(a) .+ p .* lx)
    se = n > 2 ? sqrt(sum(W .* res .^ 2) / (n - 2) / sxx) : NaN       # n = 2: exact line, no error estimate
    return (a, p, se)
end

"""
    threshold_status(x, y, level)

Threshold of a scanned error curve y(x), assumed to grow with x, as (value, status):
- `:crossed`: x where y first exceeds `level` coming from small x (every smaller scanned x has y ≤ level), by log-log
  interpolation between the last point at or below and the first point above `level`;
- `:below_range`: y ≤ level at every scanned x; value = the largest scanned x, a lower bound on the threshold;
- `:above_range`: y > level already at the smallest x; value = NaN (threshold below the scanned range).
Duplicate x are merged keeping the largest y (conservative, independent of input order). Non-monotone noisy curves give
the first exceedance, the convention of τ_conv in Phase 2 (all smaller τ within the level).
"""
function threshold_status(x, y, level)
    ux = sort(unique(x))
    uy = [maximum(y[i] for i in eachindex(x) if x[i] == v) for v in ux]
    uy[1] > level && return (NaN, :above_range)
    for k in 1:length(ux) - 1
        if uy[k + 1] > level
            (uy[k] > 0 && uy[k + 1] > 0) || return (ux[k], :crossed)
            f = (log(level) - log(uy[k])) / (log(uy[k + 1]) - log(uy[k]))
            return (exp(log(ux[k]) + f * (log(ux[k + 1]) - log(ux[k]))), :crossed)
        end
    end
    return (ux[end], :below_range)
end

"Value of [`threshold_status`](@ref) (NaN when the threshold is below the scanned range)."
threshold_at(x, y, level) = threshold_status(x, y, level)[1]

"Relative error (mc − cont)/cont."
rel_error(mc, cont) = (mc - cont) / cont
