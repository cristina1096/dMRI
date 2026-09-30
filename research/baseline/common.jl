# Shared helpers for the Phase 1 baseline scripts.
#
# Every script `include`s this file. It provides:
#   - `provenance()`     : software versions, git commit, hardware, thread count
#   - `outdir(name)`     : research/results/baseline/<name>/, created on demand
#   - `write_json`, `write_csv` : minimal writers with no extra dependencies
#   - `seed_stats`       : mean / σ / SEM over seeds
#   - `apparent_T2`      : -t / log(S(t)/S(0))
#
# Nothing here modifies the simulator; all simulator calls go through its public API
# (plus a few read-only internals that are named explicitly where used).

using MCMRSimulator
import MCMRSimulator as mr
using Statistics, Random, Printf, Dates, LinearAlgebra
import JSON
import Pkg

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const RESULTS_ROOT = joinpath(REPO_ROOT, "research", "results", "baseline")

function outdir(name::AbstractString)
    d = joinpath(RESULTS_ROOT, name)
    mkpath(d)
    return d
end

gitcmd(args...) = strip(read(Cmd(`git -C $REPO_ROOT $(collect(args))`), String))

function provenance()
    deps = Pkg.dependencies()
    pkgver(name) = begin
        for (_, info) in deps
            if info.name == name
                src = isnothing(info.git_revision) ? "" : " @ $(info.git_revision)"
                return string(info.version) * src
            end
        end
        return "not loaded"
    end
    dirty = !isempty(gitcmd("status", "--porcelain", "--", "src"))
    return Dict(
        "date" => string(now()),
        "julia" => string(VERSION),
        "MCMRSimulator" => pkgver("MCMRSimulator"),
        "MRIBuilder" => pkgver("MRIBuilder"),
        "git_branch" => gitcmd("rev-parse", "--abbrev-ref", "HEAD"),
        "git_commit" => gitcmd("rev-parse", "HEAD"),
        "src_modified" => dirty,
        "cpu" => Sys.cpu_info()[1].model,
        "ncpu" => Sys.CPU_THREADS,
        "memory_GB" => round(Sys.total_memory() / 2^30, digits=1),
        "os" => string(Sys.KERNEL, " ", Sys.MACHINE),
        "julia_threads" => Threads.nthreads(),
        "script" => basename(PROGRAM_FILE),
    )
end

function write_json(path, data)
    open(path, "w") do io
        JSON.print(io, data, 2)
    end
    return path
end

"Write a vector of NamedTuples (all with the same keys) as CSV."
function write_csv(path, rows::AbstractVector{<:NamedTuple})
    open(path, "w") do io
        println(io, join(string.(keys(rows[1])), ","))
        for r in rows
            println(io, join([v isa AbstractFloat ? @sprintf("%.10g", v) : string(v) for v in values(r)], ","))
        end
    end
    return path
end

"Mean, standard deviation across seeds (σ) and standard error of the mean."
function seed_stats(x::AbstractVector{<:Real})
    n = length(x)
    m = mean(x)
    s = n > 1 ? std(x) : NaN
    return (mean=m, sigma=s, sem=s / sqrt(n), n=n)
end

apparent_T2(t, ratio) = -t / log(ratio)

"Transverse magnetisation per spin from a SpinOrientationSum."
per_spin(s) = mr.transverse(s) / s.nspins

"Print and also log a line to an IO (e.g. a log file)."
function tee(io, args...)
    println(args...)
    println(io, args...)
    flush(io)
end
