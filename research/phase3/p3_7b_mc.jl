# Phase 3.7b: paired Monte Carlo outside cylinders in square arrays, summed intra + extra signals, and crash-safe
# driver I/O. Requires harness.jl, continuum_annulus.jl and p3_7_mc.jl (for `write_csv` via common.jl) to be included.
# Outside runs set only the OUTER layer (the cylinders are impermeable, so the inner layer cannot act on outside spins;
# user 2026-11-09). Every model of a condition shares seeds and τ, so trajectories are identical (checked) and ΔS is
# paired per seed. The annulus continuum (R_out = s/√π, exponential support = gap) runs alongside as support.

"""
    h_for_cbar_out(shape, cb, gap)

Layer thickness h giving centroid `cb` for an OUTER profile in a gap `gap` (µm): step 2cb, linear 3cb, exponential by
bisection with support L = gap. Errors if the support would exceed the gap or the exponential cannot reach `cb`
(its centroid is < gap/2 for every h).
"""
function h_for_cbar_out(shape, cb, gap)
    if shape == "step" || shape == "linear"
        h = shape == "step" ? 2cb : 3cb
        h <= gap + 1e-12 * max(1.0, gap) || error("h_for_cbar_out: $shape support h = $h exceeds the gap $gap")
        return h
    end
    shape == "exponential" || error("h_for_cbar_out: profile $shape is not supported")
    cb < gap / 2 || error("h_for_cbar_out: exponential centroid $cb is unreachable with support $gap (max < $(gap / 2))")
    lo, hi = 1e-8, 1e6 * gap
    for _ in 1:300
        mid = sqrt(lo * hi)
        outer_centroid("exponential", mid, gap) < cb ? (lo = mid) : (hi = mid)
    end
    return sqrt(lo * hi)
end

"""
    run_condition_outside(; R=1.0, s=nothing, sx=s, sy=s, D, rho, R2_bulk, models, times, tau, nspins=20_000, seeds)

`models`: (name, shape, h) with shape "conventional" (MCMR θ = ρ√(π/D), P3.1.6) or an outer profile (outer layer only).
Square array of spacing `s`, or rectangular (`sx`, `sy`; P3.5), outside spins. Returns the `run_condition` fields with the
annulus continuum (equal area: π R_out² = sx·sy; exponential support = smallest gap) in place of the radial one
(`S_ann`, `dS_ann`) and `kept`, the fraction of the layer inside the annulus.
"""
function run_condition_outside(; R=1.0, s=nothing, sx=s, sy=s, D, rho, R2_bulk, models, times, tau, nspins=20_000, seeds)
    (isnothing(sx) || isnothing(sy)) && error("run_condition_outside: give s or both sx and sy")
    gap = min(sx, sy) - 2R
    Ro = sqrt(sx * sy / π)
    support(shape, h) = shape == "exponential" ? gap : h
    geom(shape, h) = shape == "conventional" ?
        Cylinders(radius=R, repeats=[sx, sy], surface_relaxation=rho * sqrt(π / D)) :
        Cylinders(radius=R, repeats=[sx, sy], layer_rho_outer=rho, layer_h_outer=h, layer_shape=shape)
    runs = Dict(name => run_cylinders(geom(shape, h); region=:outside, R=R, sx=sx, sy=sy, R2_bulk=R2_bulk, timestep=tau,
                                      times=times, seeds=seeds, nspins=nspins, D=D) for (name, shape, h) in models)
    ref = any(m -> m[2] == "conventional", models) ? "conventional" : first(models)[1]
    identical = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs[ref].finals, runs[n].finals)) for (n, _, _) in models)
    n = length(seeds)
    S = Dict(name => runs[name].S for (name, _, _) in models)
    S_mean = Dict(name => runs[name].mean for (name, _, _) in models)
    S_sem = Dict(name => runs[name].sem for (name, _, _) in models)
    S_ann = Dict(name => (shape == "conventional" ?
                          annulus_signal(R=R, R_out=Ro, D=D, times=times, R2_bulk=R2_bulk, surface=rho) :
                          annulus_signal(R=R, R_out=Ro, D=D, times=times, R2_bulk=R2_bulk, layer=(shape, rho, h, support(shape, h))))
                 for (name, shape, h) in models)
    kept = Dict(name => (shape == "conventional" ? 1.0 : annulus_layer_kept(shape, h, support(shape, h), R, Ro)) for (name, shape, h) in models)
    dS_mean, dS_sem, dS_ann = Dict{String, Vector{Float64}}(), Dict{String, Vector{Float64}}(), Dict{String, Vector{Float64}}()
    if haskey(S, "conventional")
        for (name, shape, _) in models
            shape == "conventional" && continue
            d = S[name] .- S["conventional"]
            dS_mean[name] = vec(mean(d, dims=1))
            dS_sem[name] = n > 1 ? vec(std(d, dims=1)) ./ sqrt(n) : zeros(length(times))
            dS_ann[name] = S_ann[name] .- S_ann["conventional"]
        end
    end
    return (times=collect(Float64, times), tau=runs[ref].tau, identical_trajectories=identical, S=S, S_mean=S_mean, S_sem=S_sem,
            S_ann=S_ann, dS_mean=dS_mean, dS_sem=dS_sem, dS_ann=dS_ann, kept=kept, nregion=runs[ref].nregion)
end

"""
    sum_compartments(f, inside, outside, name)

Summed signal of model `name`: S = f S_in + (1 − f) S_out, ΔS likewise; SEMs in quadrature (independent runs).
`inside`/`outside` need `times`, `S_mean`, `S_sem`, `dS_mean`, `dS_sem` (Dicts by model name); their times must agree.
"""
function sum_compartments(f, inside, outside, name)
    inside.times == outside.times || error("sum_compartments: inside and outside readout times differ")
    q(a, b) = sqrt.(f^2 .* a .^ 2 .+ (1 - f)^2 .* b .^ 2)
    hasd = haskey(inside.dS_mean, name) && haskey(outside.dS_mean, name)
    return (times=inside.times, S=f .* inside.S_mean[name] .+ (1 - f) .* outside.S_mean[name],
            S_sem=q(inside.S_sem[name], outside.S_sem[name]),
            dS=hasd ? f .* inside.dS_mean[name] .+ (1 - f) .* outside.dS_mean[name] : fill(NaN, length(inside.times)),
            dS_sem=hasd ? q(inside.dS_sem[name], outside.dS_sem[name]) : fill(NaN, length(inside.times)))
end

"`write_csv` to a temporary file, then rename over `path` (no partial file on a crash; 3.7a Group B minor M1)."
function write_atomic(path, rows)
    write_csv(path * ".tmp", rows)
    mv(path * ".tmp", path; force=true)
    return path
end

"""
    read_done_both(cfile, kfile; cheader, kheader)

Labels written to BOTH conditions (`cfile`) and curves (`kfile`): a label in only one of them (crash between the two
writes) is not done. Errors if an existing file's header differs from the expected one (checked before any new run).
Returns (labels, cprev, kprev) with `cprev`/`kprev` as from `read_done`, restricted to the done labels.
"""
function read_done_both(cfile, kfile; cheader, kheader)
    c, k = read_done(cfile), read_done(kfile)
    isempty(c.header) || c.header == cheader || error("read_done_both: header of $cfile differs from the expected one")
    isempty(k.header) || k.header == kheader || error("read_done_both: header of $kfile differs from the expected one")
    done = intersect(c.labels, k.labels)
    keep(prev) = (header=prev.header, lines=[l for l in prev.lines if split(l, ',')[findfirst(==("label"), split(prev.header, ','))] in done],
                  labels=done)
    return (labels=done, cprev=isempty(c.header) ? c : keep(c), kprev=isempty(k.header) ? k : keep(k))
end
