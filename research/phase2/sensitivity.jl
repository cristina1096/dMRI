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
