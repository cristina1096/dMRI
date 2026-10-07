# Timestep-study tools (P2.6): toy-based rate per τ, converged timestep, c_h.
#
# Observable: R = −ln S(T_eff)/T_eff, the effective decay rate after T_eff = round(T/τ)·τ.
# τ_conv(tol): the largest τ such that every τ' ≤ τ in the grid has |R(τ') − R_ref| / R_ref ≤ tol,
# with R_ref = R at the smallest τ; interpolated linearly in log τ to where the error crosses tol.

"Log-spaced timesteps from lo to hi (ms), `per_decade` points per decade, both ends included."
log_taus(lo, hi; per_decade=3) = 10 .^ range(log10(lo), log10(hi), length=round(Int, per_decade * log10(hi / lo)) + 1)

"""
    toy_rate(; h, rate, D, tau, T, nspins, seeds, method=:exact, shape="step")

Effective rate R = −ln S/T_eff per seed from the 1D toy (walls every W_REF, step layer on both faces),
with T_eff = round(T/τ)·τ the time actually simulated. Returns (mean, sem, per_seed, T_eff).
"""
function toy_rate(; h, rate, D, tau, T, nspins, seeds, method=:exact, shape="step")
    T_eff = round(Int, T / tau) * tau
    Rs = [-log(mean(toy_layer_signal(; h=h, rate=rate, tau=tau, T=T, nspins=nspins, seed=s, method=method, D=D, shape=shape))) / T_eff
          for s in seeds]
    return (mean=mean(Rs), sem=length(Rs) > 1 ? std(Rs) / sqrt(length(Rs)) : NaN, per_seed=Rs, T_eff=T_eff)
end

"Largest grid τ such that all τ' ≤ τ have relative error ≤ tol. `all_pass` is true when no τ fails."
function tau_conv_grid(taus, R, R_ref; tol)
    order = sortperm(taus)
    best = taus[order[1]]
    for k in order
        abs(R[k] - R_ref) / R_ref <= tol || return (tau=best, all_pass=false)
        best = taus[k]
    end
    return (tau=best, all_pass=true)
end

"""
    tau_conv_interp(taus, R, R_ref; tol)

As `tau_conv_grid`, but interpolated: linear in log τ between the last passing and the first failing grid point,
to where the relative error reaches `tol`. If no τ fails, returns the largest τ with `lower_bound = true`.
"""
function tau_conv_interp(taus, R, R_ref; tol)
    order = sortperm(taus)
    err(k) = abs(R[k] - R_ref) / R_ref
    for (i, k) in enumerate(order)
        err(k) <= tol && continue
        i == 1 && return (tau=taus[k], lower_bound=false)
        kp = order[i - 1]
        f = (tol - err(kp)) / (err(k) - err(kp))
        return (tau=exp(log(taus[kp]) + f * (log(taus[k]) - log(taus[kp]))), lower_bound=false)
    end
    return (tau=taus[order[end]], lower_bound=true)
end

"True when the reference's relative SEM is more than a third of the tolerance (τ_conv not resolvable)."
noise_limited(R_ref, R_ref_sem; tol) = R_ref_sem / R_ref > tol / 3

"c_h from τ_conv = c_h·h²/D: mean of the ratios τ_conv·D/h², with their relative spread (std/mean)."
function fit_c_h(tau_convs, hs, Ds)
    ratios = tau_convs .* Ds ./ hs .^ 2
    return (c_h=mean(ratios), spread=length(ratios) > 1 ? std(ratios) / mean(ratios) : 0.0, ratios=ratios)
end
