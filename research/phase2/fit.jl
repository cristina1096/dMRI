# Phase 2 analysis helpers: sweep spin count, decay-shape curvature, monotonicity, h* fit (Task 4 adds the fit).

"Spins per seed for the h sweep: `recommended_nspins` from P2.1.4, else NSPINS_REF."
function sweep_nspins()
    f = joinpath(REPO_ROOT, "research", "results", "phase2", "p2_1_4_noise_floor", "summary.json")
    isfile(f) || return NSPINS_REF
    return Int(JSON.parsefile(f)["recommended_nspins"])
end

"c in the least-squares fit ln S = a + b·t + c·t² (curvature of the log-decay; 0 for a single exponential)."
function curvature(times, S)
    X = hcat(ones(length(times)), times, times .^ 2)
    return (X \ log.(S))[3]
end

"Entry k: attenuation rises from step k to k+1 by more than 2 combined SEM."
increasing_steps(att, sem) = [att[k + 1] - att[k] > 2 * sqrt(sem[k]^2 + sem[k + 1]^2) for k in 1:length(att) - 1]

"Per-readout linear interpolation in h between grid curves; clamped to the grid ends."
function interp_curve(hgrid, curves, h)
    h <= hgrid[1] && return curves[1]
    h >= hgrid[end] && return curves[end]
    j = searchsortedlast(hgrid, h)
    f = (h - hgrid[j]) / (hgrid[j + 1] - hgrid[j])
    return (1 - f) .* curves[j] .+ f .* curves[j + 1]
end

"""
    fit_h(hgrid, model_mean, model_sem, ref_mean, ref_sem; npts=10_001)

h minimising χ²(h) = Σ_t [(S_model(t; h) − S_ref(t)) / σ_t]², σ_t = √(sem_model² + sem_ref²), over a dense
h grid with the model curve interpolated in h. Readouts with σ_t = 0 are excluded. `at_edge` is true when
the best h is the first or last grid value (the reference is not reached inside the grid).
"""
function fit_h(hgrid, model_mean, model_sem, ref_mean, ref_sem; npts=10_001)
    best_h, best_chi2 = NaN, Inf
    for h in range(hgrid[1], hgrid[end], length=npts)
        m = interp_curve(hgrid, model_mean, h)
        σ = sqrt.(interp_curve(hgrid, model_sem, h) .^ 2 .+ ref_sem .^ 2)
        ok = σ .> 0
        χ2 = sum(((m .- ref_mean)[ok] ./ σ[ok]) .^ 2)
        χ2 < best_chi2 && ((best_h, best_chi2) = (h, χ2))
    end
    m = interp_curve(hgrid, model_mean, best_h)
    σ = sqrt.(interp_curve(hgrid, model_sem, best_h) .^ 2 .+ ref_sem .^ 2)
    ok = σ .> 0
    residuals = [ok[i] ? (m[i] - ref_mean[i]) / σ[i] : 0.0 for i in eachindex(m)]
    return (h=best_h, chi2=best_chi2, dof=count(ok) - 1, residuals=residuals, model=m,
        at_edge=(best_h == hgrid[1] || best_h == hgrid[end]))
end

"""
    jackknife_h(hgrid, S_by_h, ref_mean, ref_sem)

h* from the mean over all seeds, and its jackknife standard error: refit leaving out one seed at a time,
σ = √((n−1)/n · Σ (h₋ₛ − mean h₋)²). `S_by_h[i]` is the seeds × times matrix for hgrid[i].
"""
function jackknife_h(hgrid, S_by_h, ref_mean, ref_sem)
    n = size(S_by_h[1], 1)
    fit_rows(rows) = fit_h(hgrid,
        [vec(mean(S[rows, :], dims=1)) for S in S_by_h],
        [vec(std(S[rows, :], dims=1)) ./ sqrt(length(rows)) for S in S_by_h],
        ref_mean, ref_sem).h
    full = fit_rows(1:n)
    loo = [fit_rows(setdiff(1:n, s)) for s in 1:n]
    return (h=full, sigma=sqrt((n - 1) / n * sum((loo .- mean(loo)) .^ 2)))
end

"Reads sweep.csv (h_um, seed, t_ms, signal) into (hgrid, S_by_h)."
function read_sweep(path)
    lines = readlines(path)
    header = split(lines[1], ",")
    col(name) = findfirst(==(name), header)
    data = [parse.(Float64, split(l, ",")) for l in lines[2:end]]
    hgrid = sort(unique(r[col("h_um")] for r in data))
    S_by_h = map(hgrid) do h
        rs = filter(r -> r[col("h_um")] == h, data)
        seeds = sort(unique(Int(r[col("seed")]) for r in rs))
        times = sort(unique(r[col("t_ms")] for r in rs))
        M = zeros(length(seeds), length(times))
        for r in rs
            M[findfirst(==(Int(r[col("seed")])), seeds), findfirst(==(r[col("t_ms")]), times)] = r[col("signal")]
        end
        M
    end
    return (hgrid, S_by_h)
end
