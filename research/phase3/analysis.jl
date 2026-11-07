# Phase 3 analysis helpers (P3.6): power-law fits and thresholds on scanned error curves.

"""
    fit_power(x, y; w=nothing)

Least squares of log y = log a + p log x (weights `w` on the log points); returns (a, p, se_p), with se_p the standard
error of the slope from the residual variance (n − 2 degrees of freedom; 0 for n = 2). Requires x, y > 0.
"""
function fit_power(x, y; w=nothing)
    (all(>(0), x) && all(>(0), y)) || throw(ArgumentError("fit_power needs x, y > 0"))
    lx, ly = log.(x), log.(y)
    W = isnothing(w) ? ones(length(x)) : collect(Float64, w)
    sw = sum(W)
    mx, my = sum(W .* lx) / sw, sum(W .* ly) / sw
    sxx = sum(W .* (lx .- mx) .^ 2)
    p = sum(W .* (lx .- mx) .* (ly .- my)) / sxx
    a = exp(my - p * mx)
    n = length(x)
    res = ly .- (log(a) .+ p .* lx)
    se = n > 2 ? sqrt(sum(W .* res .^ 2) / (n - 2) / sxx) : 0.0
    return (a, p, se)
end

"""
    threshold_at(x, y, level)

Largest x in the scanned range with y ≤ level, assuming y increases with x: the last point at or below `level`, then
log-log interpolation towards the next point if that one is above. NaN if the first point is already above `level`.
"""
function threshold_at(x, y, level)
    o = sortperm(x)
    xs, ys = x[o], y[o]
    ys[1] > level && return NaN
    for k in 1:length(xs) - 1
        if ys[k + 1] > level
            f = (log(level) - log(ys[k])) / (log(ys[k + 1]) - log(ys[k]))
            return exp(log(xs[k]) + f * (log(xs[k + 1]) - log(xs[k])))
        end
    end
    return xs[end]
end

"Relative error (mc − cont)/cont."
rel_error(mc, cont) = (mc - cont) / cont
