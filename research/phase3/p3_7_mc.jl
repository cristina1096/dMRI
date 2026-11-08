# Phase 3.7a: shared paired Monte Carlo runner. One condition (R, D, ρ, bulk R2, τ) runs the conventional model (MCMR θ =
# ρ√(π/D), P3.1.6) and any layer models inside a single cylinder at the same τ and seeds, so their trajectories are
# identical (checked) and differences are paired per seed. Every model's exact radial-continuum S(t) is returned
# alongside (support method). Requires harness.jl and continuum_radial.jl to be included first.

"""
    run_condition(; R, D, rho, R2_bulk, models, times, tau, nspins=20_000, seeds=1:10)

`models`: vector of (name, shape, h), shape ∈ ("conventional", "step", "linear", "exponential"); layers are inner-only
(ρ per unit surface area, D3.1; exponential support R, D3.5). Spins start uniformly in the box [−R, R]³ and only
those inside the cylinder are used. Returns a NamedTuple with `times`, `tau`, `identical_trajectories` and, per model
name: `S` (seeds × times), `S_mean`, `S_sem`, `S_cont`, and for layer models the paired difference against
"conventional": `dS_mean`, `dS_sem` (per-seed differences), `dS_cont`.
"""
function run_condition(; R, D, rho, R2_bulk, models, times, tau, nspins=20_000, seeds=1:10)
    bb = mr.BoundingBox([-R, -R, -R], [R, R, R])
    geom(shape, h) = shape == "conventional" ? Cylinders(radius=R, surface_relaxation=rho * sqrt(π / D)) :
                     Cylinders(radius=R, layer_rho_inner=rho, layer_h_inner=h, layer_shape=shape)
    runs = Dict(name => run_cylinders(geom(shape, h); region=:inside, R=R, bounding_box=bb, R2_bulk=R2_bulk, timestep=tau,
                                      times=times, seeds=seeds, nspins=nspins, D=D) for (name, shape, h) in models)
    ref = any(m -> m[2] == "conventional", models) ? "conventional" : first(models)[1]
    identical = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs[ref].finals, runs[n].finals)) for (n, _, _) in models)
    n = length(seeds)
    S = Dict(name => runs[name].S for (name, _, _) in models)
    S_mean = Dict(name => runs[name].mean for (name, _, _) in models)
    S_sem = Dict(name => runs[name].sem for (name, _, _) in models)
    S_cont = Dict(name => (shape == "conventional" ? radial_signal(R=R, D=D, times=times, R2_bulk=R2_bulk, surface=rho) :
                                                     radial_signal(R=R, D=D, times=times, R2_bulk=R2_bulk, layer=(shape, rho, h)))
                  for (name, shape, h) in models)
    dS_mean, dS_sem, dS_cont = Dict{String, Vector{Float64}}(), Dict{String, Vector{Float64}}(), Dict{String, Vector{Float64}}()
    if haskey(S, "conventional")
        for (name, shape, _) in models
            shape == "conventional" && continue
            d = S[name] .- S["conventional"]
            dS_mean[name] = vec(mean(d, dims=1))
            dS_sem[name] = n > 1 ? vec(std(d, dims=1)) ./ sqrt(n) : zeros(length(times))
            dS_cont[name] = S_cont[name] .- S_cont["conventional"]
        end
    end
    return (times=collect(Float64, times), tau=runs[ref].tau, identical_trajectories=identical, S=S, S_mean=S_mean, S_sem=S_sem,
            S_cont=S_cont, dS_mean=dS_mean, dS_sem=dS_sem, dS_cont=dS_cont)
end

# Resume support for the long drivers (Group A review M3): rows already written are kept verbatim and their labels skipped.
"Rows already in a CSV written by `write_csv`: (header, lines, labels of the `label` column); empty when absent."
function read_done(path)
    isfile(path) || return (header="", lines=String[], labels=Set{String}())
    ls = readlines(path)
    k = findfirst(==("label"), split(ls[1], ','))
    return (header=ls[1], lines=ls[2:end], labels=Set(String(split(l, ',')[k]) for l in ls[2:end]))
end

"Write the previous rows of `read_done` followed by `rows` (NamedTuples, `write_csv` format); headers must agree."
function write_merged(path, prev, rows)
    if isempty(rows)
        isempty(prev.header) || open(io -> foreach(l -> println(io, l), vcat(prev.header, prev.lines)), path, "w")
        return path
    end
    tmp = write_csv(path * ".new", rows)
    new = readlines(tmp)
    rm(tmp)
    isempty(prev.header) || new[1] == prev.header || error("write_merged: header mismatch in $path")
    open(io -> foreach(l -> println(io, l), vcat(new[1], prev.lines, new[2:end])), path, "w")
    return path
end
