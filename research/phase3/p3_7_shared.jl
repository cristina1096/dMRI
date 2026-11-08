# Phase 3.7a shared analysis helpers (Tasks 3–6): profile centroids, continuum S(t) closures and the Monte Carlo
# quantities at TE* read back from p3_7_run.jl outputs. Requires harness.jl, continuum_radial.jl and analysis.jl.

readrows(f) = (ls = readlines(f); h = split(ls[1], ','); [Dict(zip(h, split(l, ','))) for l in ls[2:end]])
num(r, k) = parse(Float64, r[k])

cbar(shape, h, R) = shape == "step" ? h / 2 : shape == "linear" ? h / 3 : (c = R / h; h * (1 - (c + 1) * exp(-c)) / -expm1(-c))
function h_for_cbar(shape, cb, R)
    shape == "step" && return 2cb
    shape == "linear" && return 3cb
    lo, hi = 1e-6, 100R
    for _ in 1:200
        mid = sqrt(lo * hi); cbar("exponential", mid, R) < cb ? (lo = mid) : (hi = mid)
    end
    return sqrt(lo * hi)
end
"S(t) closure from one eigendecomposition of the radial operator (same discretisation as radial_signal)."
function signal_fun(; R, D, layer=nothing, surface=0.0, bulk=0.0, N=2000)
    (M, A) = radial_operator(; R=R, D=D, layer=layer, surface=surface, N=N)
    F = eigen(M)
    c2 = (F.vectors' * sqrt.(A)) .^ 2 ./ sum(A)
    return t -> sum(c2 .* exp.(-F.values .* t)) * exp(-bulk * t)
end
t2app_sigma(S, sem, te) = te * sem / (S * log(S)^2)

"""
    at_te_rows(src)

Per condition and model of a p3_7_run.jl output directory: S, apparent T₂ and paired ΔT₂/T₂ at the TE* readout (exact
readout, checked), with the MC − continuum z values (formulas: header of p3_7_fig9a_analysis.jl).
"""
function at_te_rows(src)
    conds = readrows(joinpath(src, "conditions.csv"))
    curves = readrows(joinpath(src, "curves.csv"))
    rows = NamedTuple[]
    for c in conds
        label, R, te = c["label"], num(c, "R"), num(c, "te")
        sel = [r for r in curves if r["label"] == label]
        ts = sort(unique(num(r, "t") for r in sel))
        k = argmin(abs.(ts .- te))
        abs(ts[k] - te) <= 1e-9 * te || error("$label: no readout at TE*")
        at(m) = only(r for r in sel if r["model"] == m && num(r, "t") == ts[k])
        conv = at("conventional")
        Sc = num(conv, "S_mean")
        T2c = apparent_t2(Sc, te)
        for m in unique(r["model"] for r in sel)
            r = at(m)
            S, sem, Scont = num(r, "S_mean"), num(r, "S_sem"), num(r, "S_cont")
            T2 = apparent_t2(S, te)
            layer = m != "conventional"
            dT2 = layer ? (T2 - T2c) / T2c : NaN
            dT2σ = layer ? t2app_sigma(S, num(r, "dS_sem"), te) / T2c : NaN
            dT2cont = layer ? (apparent_t2(Scont, te) - apparent_t2(num(conv, "S_cont"), te)) / apparent_t2(num(conv, "S_cont"), te) : NaN
            push!(rows, (label=label, matching=occursin("_h", label) ? "same_h" : "same_centroid", R=R, te=te, model=m, shape=r["shape"],
                h=num(r, "h"), cbar=num(r, "cbar"), S_mc=S, S_sem=sem, S_cont=Scont, z_S=(S - Scont) / sem,
                T2_mc=T2, T2_sigma=t2app_sigma(S, sem, te), T2_cont=apparent_t2(Scont, te),
                dT2_rel=dT2, dT2_rel_sigma=dT2σ, dT2_rel_cont=dT2cont,
                noise2_snr100=2 * t2_noise_rel(Sc, 100), noise2_snr200=2 * t2_noise_rel(Sc, 200),
                detectable_snr200=layer ? abs(dT2) > 2 * t2_noise_rel(Sc, 200) : false,
                dS_mc=layer ? num(r, "dS_mean") : NaN, dS_sem=layer ? num(r, "dS_sem") : NaN, dS_cont=layer ? num(r, "dS_cont") : NaN,
                z_dS=layer ? (num(r, "dS_mean") - num(r, "dS_cont")) / num(r, "dS_sem") : NaN,
                tau=num(c, "tau"), identical=c["identical_trajectories"] == "true"))
        end
    end
    return rows
end

"ΔS_mc/ΔS_cont ± SEM_ΔS/|ΔS_cont| per layer row of `at_te_rows` (Group A review I2)."
dS_ratio_rows(rows) = [(label=r.label, model=r.model, tau=r.tau, ratio=r.dS_mc / r.dS_cont, ratio_sigma=r.dS_sem / abs(r.dS_cont), z_dS=r.z_dS)
        for r in rows if r.model != "conventional"]
