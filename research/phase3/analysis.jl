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
"Two-point apparent T₂ = −TE/ln(S(TE)/S(0)) for a signal normalised to S(0) = 1."
apparent_t2(S_TE, TE) = -TE / log(S_TE)

"Relative noise of the two-point T₂ at SNR `snr`: σ_T₂/T₂ = √(1 + r²)/(SNR·r·|ln r|), r = S(TE)/S(0)."
t2_noise_rel(r, snr) = sqrt(1 + r^2) / (snr * r * abs(log(r)))

"Fastest decay of any surface sink inside a cylinder (absorbing wall): D·j₀₁²/R², j₀₁ = 2.404825557695773."
conventional_ceiling(R, D) = D * (2.404825557695773 / R)^2

"The t in [0, t_hi] with Sfun(t) = e^(−1) for a monotonically decreasing Sfun (bisection, 200 steps)."
function te_star(Sfun; t_hi)
    lo, hi = 0.0, Float64(t_hi)
    Sfun(hi) > exp(-1) && error("te_star: S(t_hi) > e^(−1); increase t_hi")
    for _ in 1:200
        mid = (lo + hi) / 2
        Sfun(mid) > exp(-1) ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

"""
    best_fit_conventional(S_layer, times, Sconv_of_rho; rho0, lo=exp(-4), hi=exp(8), ngrid=61, iters=60)

Minimax fit of the conventional model's ρ_c to a layer curve: ρ_c minimising max_t |S_layer(t) − S_conv(t; ρ_c)|.
A log-grid scan over [lo·rho0, hi·rho0], then golden-section search between the neighbours of the best grid point.
Returns (rho_c, residual, at_edge); at_edge = true when the best grid point is the first or last (the fit wants to
leave the searched range, e.g. beyond the conventional ceiling). A plateau that ties the minimum within 1e-9 up to
either end of the grid (a family that has stopped changing, e.g. near D·j₀₁²/R²) also counts as at the edge: its ρ_c is
not determined and callers must not report ρ_c or the residual for it.
"""
function best_fit_conventional(S_layer, times, Sconv_of_rho; rho0, lo=exp(-4), hi=exp(8), ngrid=61, iters=60)
    obj(ρ) = maximum(abs.(S_layer .- Sconv_of_rho(ρ)))
    grid = rho0 .* exp.(range(log(lo), log(hi), length=ngrid))
    vals = obj.(grid)
    k = argmin(vals)
    at_edge = vals[1] - vals[k] <= 1e-9 || vals[end] - vals[k] <= 1e-9
    a, b = log(grid[max(k - 1, 1)]), log(grid[min(k + 1, ngrid)])
    g = (sqrt(5) - 1) / 2
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = obj(exp(c)), obj(exp(d))
    for _ in 1:iters
        if fc < fd
            b, d, fd = d, c, fc
            c = b - g * (b - a); fc = obj(exp(c))
        else
            a, c, fc = c, d, fd
            d = a + g * (b - a); fd = obj(exp(d))
        end
    end
    ρc = exp((a + b) / 2)
    return (ρc, obj(ρc), at_edge)
end
