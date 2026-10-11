# Sensitivity-study helpers (P2.5): paired differences, peaks, mono-exponential residual, profile centroids.

"Per-readout mean and SEM over seeds of the paired difference Sa − Sb (rows = seeds, same seeds in both)."
function paired_difference(Sa::AbstractMatrix, Sb::AbstractMatrix)
    D = Sa .- Sb
    n = size(D, 1)
    return (mean=vec(mean(D, dims=1)), sem=n > 1 ? vec(std(D, dims=1)) ./ sqrt(n) : zeros(size(D, 2)), per_seed=D)
end

"Index, time and value of the largest |m|."
function peak(times, m)
    k = argmax(abs.(m))
    return (index=k, t=times[k], value=m[k])
end

"""
    monoexp_residual(times, S)

Best mono-exponential S ≈ e^(−R t) by least squares on ln S through the origin (t > 0), and the largest
|S − e^(−R t)| over the readouts. Stands for the best the conventional model can do when it is mono-exponential.
"""
function monoexp_residual(times, S)
    k = findall(t -> t > 0, times)
    R = -sum(times[k] .* log.(S[k])) / sum(times[k] .^ 2)
    return (R=R, residual=maximum(abs.(S .- exp.(-R .* times))))
end

"""
    centroid(shape, h; cutoff=Inf)

Mean distance from the wall of the relaxation, ∫ d·g(d/h) dd / ∫ g(d/h) dd: step h/2, linear h/3,
exponential truncated at `cutoff`: h·(1 − (c + 1)e^(−c))/(1 − e^(−c)), c = cutoff/h.
"""
function centroid(shape, h; cutoff=Inf)
    shape == "step" && return h / 2
    shape == "linear" && return h / 3
    shape == "exponential" || error("Unknown shape $shape")
    c = cutoff / h
    return isinf(c) ? h : h * (1 - (c + 1) * exp(-c)) / -expm1(-c)
end

"Reference point of the sensitivity study (Global Constraints of the plan)."
const SENS_REF = (h=0.2, w=2.0, D=3.0, rho=0.01)
const SENS_SHAPES = ("step", "linear", "exponential")

"""
    sweep_configs(sweep)

One-at-a-time configurations around SENS_REF for sweep "h", "w", "D" or "rho": every shape at every value,
except the step profile where h > w/2 (overlapping step layers relax the middle of the gap most; excluded).
"""
function sweep_configs(sweep::AbstractString)
    values = sweep == "h" ? [0.05, 0.1, 0.2, 0.5, 1.0, 2.0] :
             sweep == "w" ? [1.0, 2.0, 5.0, 10.0] :
             sweep == "D" ? [1.0, 3.0, 6.0] :
             sweep == "rho" ? [0.003, 0.01, 0.03, 0.1] : error("Unknown sweep $sweep; use h, w, D or rho.")
    key = Symbol(sweep)
    out = NamedTuple[]
    for v in values, shape in SENS_SHAPES
        c = merge(SENS_REF, NamedTuple{(key,)}((v,)))
        (shape == "step" && c.h > c.w / 2) && continue
        push!(out, (shape=shape, h=c.h, w=c.w, D=c.D, rho=c.rho))
    end
    return out
end

"Apparent T2 (ms) from the ratio S(TE)/S(0), as in Cottaar et al. 2026 (Figs. 8–9)."
apparent_t2(S0, STE, TE) = -TE / log(STE / S0)

"""
    t2_noise_rel(r, snr)

Relative σ of T2 from a two-point ratio r = S(TE)/S(0) with independent Gaussian noise of σ = 1/snr (in units of S(0))
on both points: σ_T2/T2 = √(1 + r²)/(snr·r·|ln r|).
"""
t2_noise_rel(r, snr) = sqrt(1 + r^2) / (snr * r * abs(log(r)))

"""
    paired_t2(Sa, Sb, times, TE)

Per-seed apparent T2 at `TE` for two runs on the same seeds (rows), and their paired difference:
(T2a, T2a_sem, T2b, T2b_sem, dT2, dT2_sem, rel, rel_sem) with rel = (T2a − T2b)/T2b per seed.
"""
function paired_t2(Sa::AbstractMatrix, Sb::AbstractMatrix, times, TE)
    k = findfirst(==(TE), times)
    isnothing(k) && error("TE = $TE ms is not one of the readout times")
    i0 = findfirst(==(0.0), times)
    ta = apparent_t2.(Sa[:, i0], Sa[:, k], TE)
    tb = apparent_t2.(Sb[:, i0], Sb[:, k], TE)
    se(v) = length(v) > 1 ? std(v) / sqrt(length(v)) : 0.0
    rel = (ta .- tb) ./ tb
    return (T2a=mean(ta), T2a_sem=se(ta), T2b=mean(tb), T2b_sem=se(tb), dT2=mean(ta .- tb), dT2_sem=se(ta .- tb),
        rel=mean(rel), rel_sem=se(rel))
end

"Timestep for an apparent-T2 configuration: layer constraint, 1e-2 cap and MCMR's tortuosity default 0.03·w²/D."
tau_t2(rho, h, shape; w, D) = min(tau_for(rho, h, shape; w=w), 0.03 * w^2 / D)

"Reference point of the apparent-T2 study: ρ matched to the paper's tissue T2 of 70 ms at w = 2 (2ρ/w = 1/70)."
const T2_REF = (h=0.2, w=2.0, D=3.0, rho=1 / 70)

"""
    h_for_centroid(shape, cbar; w)

Thickness h giving centroid `cbar` (µm) for a profile in a gap of width `w` (exponential truncated at w):
step 2c̄, linear 3c̄, exponential by bisection (its centroid rises with h, is ≤ h, and stays below w/2).
Returns `nothing` when the profile cannot reach c̄ with a valid h (step h ≤ w/2; linear and exponential h ≤ w).
"""
function h_for_centroid(shape, cbar; w)
    if shape == "step"
        h = 2cbar
        return h <= w / 2 ? h : nothing
    elseif shape == "linear"
        h = 3cbar
        return h <= w ? h : nothing
    elseif shape == "exponential"
        cbar < w / 2 || return nothing
        f(h) = centroid("exponential", h; cutoff=w) - cbar
        lo, hi = cbar, cbar
        while f(hi) < 0
            hi *= 2
            hi > 1e8 * w && return nothing
        end
        for _ in 1:200
            mid = (lo + hi) / 2
            f(mid) < 0 ? (lo = mid) : (hi = mid)
            hi - lo <= 1e-14 * hi && break
        end
        h = (lo + hi) / 2
        return h <= w ? h : nothing
    end
    error("Unknown shape $shape")
end

"""
    t2_configs(sweep)

One-at-a-time configurations around T2_REF for sweep "w", "h", "D" or "rho": every layer shape at every value,
except the step profile where h > w/2. Sweep "te" gives the Fig. 8A configurations (all have h ≤ w/2).
`match = :centroid`: same sweeps but the profiles share the centroid c̄ instead of h (h from `h_for_centroid`;
field `cbar` added); the h sweep is replaced by the centroid sweep "c".
"""
function t2_configs(sweep::AbstractString; match::Symbol=:h)
    match in (:h, :centroid) || error("match must be :h or :centroid, got $match")
    if match == :centroid
        sweep == "h" && error("with match = :centroid use the centroid sweep \"c\" instead of \"h\"")
        cref, R0 = centroid("step", T2_REF.h), T2_REF                # c̄ = 0.1 µm
        pts = sweep == "te" ? [(cref, 1.0, R0.D, R0.rho), (cref, 2.0, R0.D, R0.rho), (cref, 5.0, R0.D, R0.rho),
                               (0.025, 2.0, R0.D, R0.rho), (0.25, 2.0, R0.D, R0.rho)] :
              sweep == "c" ? [(c, R0.w, R0.D, R0.rho) for c in [0.01, 0.025, 0.05, 0.1, 0.25, 0.5]] :
              sweep == "w" ? [(cref, w, R0.D, R0.rho) for w in [0.5, 1.0, 2.0, 3.0, 5.0, 10.0]] :
              sweep == "D" ? [(cref, R0.w, D, R0.rho) for D in [0.5, 1.0, 2.0, 3.0]] :
              sweep == "rho" ? [(cref, R0.w, R0.D, r) for r in [0.005, 1 / 70, 0.03, 0.06]] :
              error("Unknown sweep $sweep; use te, c, w, D or rho with match = :centroid.")
        out = NamedTuple[]
        for (cbar, w, D, rho) in pts, shape in SENS_SHAPES
            h = h_for_centroid(shape, cbar; w=w)
            isnothing(h) && continue
            push!(out, (shape=shape, h=h, w=w, D=D, rho=rho, cbar=cbar))
        end
        return out
    end
    sweep == "c" && error("the centroid sweep \"c\" needs match = :centroid")
    if sweep == "te"            # Fig. 8A analogue: T2 vs TE for several S/V (w) and two extra h at w = 2
        pts = [(T2_REF.h, 1.0), (T2_REF.h, 2.0), (T2_REF.h, 5.0), (0.05, 2.0), (0.5, 2.0)]
        return [(shape=shape, h=h, w=w, D=T2_REF.D, rho=T2_REF.rho) for (h, w) in pts for shape in SENS_SHAPES]
    end
    values = sweep == "w" ? [0.5, 1.0, 2.0, 3.0, 5.0, 10.0] :
             sweep == "h" ? [0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0] :
             sweep == "D" ? [0.5, 1.0, 2.0, 3.0] :
             sweep == "rho" ? [0.005, 1 / 70, 0.03, 0.06] : error("Unknown sweep $sweep; use w, h, D or rho.")
    key = Symbol(sweep)
    out = NamedTuple[]
    for v in values, shape in SENS_SHAPES
        c = merge(T2_REF, NamedTuple{(key,)}((v,)))
        (shape == "step" && c.h > c.w / 2) && continue
        push!(out, (shape=shape, h=c.h, w=c.w, D=c.D, rho=c.rho))
    end
    return out
end
