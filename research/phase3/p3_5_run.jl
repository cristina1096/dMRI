# Phase 3.5 driver (plan 2026-11-11-overlapping-layers-arrays.md, Tasks 4–5): paired Monte Carlo outside cylinders with
# overlapping outer layers. TASK = overlap | ddep.
#   overlap: square 2.5 (two-direction) and rectangle 2.5 × 4.0 (one-direction); R = 1, D = 3, no bulk, ρR/D = 1 (ρ = 3);
#            conventional + step/linear/exponential at h ∈ {0.15, 0.25, 0.35, 0.5} (h/gap_x 0.3 … 1.0; exponential
#            support = smallest gap 0.5); 41 readouts to ln 100/λ₁ of the equal-area annulus (window only).
#   ddep:    square 2.5, ρ = 1, bulk 1/2000, D ∈ {0.5, 3}; h ∈ {0.15, 0.35} (3.7b ddep has h 0.25, 0.5 with the same TE
#            convention); TE = annulus conventional T₂ at D = 3 (exact readout 11 of 0:TE/10:1.5 TE).
# τ = min(0.135·c̄_min², 0.005/ΔR₂,max) with ΔR₂,max = max over the condition's models of the package's SUMMED bound
# `max_layer_rate` (overlap included; 3.7 final-review lesson). Seeds 1000·(offset + k) + (1:10) per condition.
# Output results/phase3/p3_5_<task>/ via curves-then-conditions atomic writes; resume skips labels in both files.
# Run: TASK=overlap julia --project=research/baseline -t 8 research/phase3/p3_5_run.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
include(joinpath(@__DIR__, "continuum_annulus.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
include(joinpath(@__DIR__, "p3_7_mc.jl"))
include(joinpath(@__DIR__, "p3_7_shared.jl"))
include(joinpath(@__DIR__, "p3_7b_mc.jl"))
const TASK = get(ENV, "TASK", "overlap")
const OUT = phase3_outdir("p3_5_" * TASK)
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const NSEEDS = parse(Int, get(ENV, "NSEEDS", "10"))
const R1, D3, R2B = 1.0, 3.0, 1 / 2000
const SHAPES = ("step", "linear", "exponential")
const OFFSET = Dict("overlap" => 500, "ddep" => 600)
const Layers = mr.Geometries.Internal.Layers

Ro(sx, sy) = sqrt(sx * sy / π)
gapof(sx, sy) = min(sx, sy) - 2R1
cbar5(sh, h, sx, sy) = outer_centroid(sh, h, sh == "exponential" ? gapof(sx, sy) : h)
"Package summed peak-rate bound (overlap included) of one outer-layer model."
function pkg_max_rate(sh, h, ρ, sx, sy, D)
    g = Cylinders(radius=R1, repeats=[sx, sy], layer_rho_outer=ρ, layer_h_outer=h, layer_shape=sh)
    return Layers.max_layer_rate(Simulation([]; geometry=g, diffusivity=D, verbose=false).geometry)
end
models_h(hs) = [("$(sh)_h$(h)", sh, h) for sh in SHAPES for h in hs]

# conditions: (label, sx, sy, D, ρ, bulk, layer models, times, te)
function conditions(task)
    c = []
    if task == "overlap"
        for (name, sx, sy) in (("sq2.5", 2.5, 2.5), ("rect2.5x4.0", 2.5, 4.0))
            ρ = 1.0 * D3 / R1
            λ1 = annulus_rates(R=R1, R_out=Ro(sx, sy), D=D3, surface=ρ)[1]
            push!(c, ("$(name)_rhoRD1", sx, sy, D3, ρ, 0.0, models_h((0.15, 0.25, 0.35, 0.5)),
                      collect(range(0.0, log(100) / λ1, length=41)), NaN))
        end
    elseif task == "ddep"
        sx = sy = 2.5
        ρ = 1.0
        te = te_star(annulus_signal_fun(R=R1, R_out=krogh_radius(2.5), D=D3, surface=ρ, bulk=R2B);
                     t_hi=50 / (R2B + ρ * 2π * R1 / (2.5^2 - π)))
        for D in (0.5, 3.0)
            push!(c, ("sq2.5_rho1.0_D$(D)", sx, sy, D, ρ, R2B, models_h((0.15, 0.35)), collect(range(0.0, 1.5te, length=16)), te))
        end
    else
        error("unknown TASK $task")
    end
    return c
end

const CKEYS = (:task, :label, :sx, :sy, :gap, :D, :rho, :R2_bulk, :tau, :cbar_min, :max_rate_summed, :max_rate_tau, :te, :t_end,
               :identical_trajectories, :runtime_s, :nspins, :nseeds, :seed_first)
const KKEYS = (:label, :model, :shape, :h, :cbar, :kept, :t, :S_mean, :S_sem, :S_sup, :dS_mean, :dS_sem, :dS_sup)
const CFILE, KFILE = joinpath(OUT, "conditions.csv"), joinpath(OUT, "curves.csv")
const PREV = get(ENV, "RESUME", "1") != "0" ? read_done_both(CFILE, KFILE; cheader=join(CKEYS, ","), kheader=join(KKEYS, ",")) :
                                               read_done_both("", ""; cheader="", kheader="")
isempty(PREV.labels) || println("resuming: skipping ", join(sort(collect(PREV.labels)), ", "))
function merge_atomic(path, prev, rows)
    write_merged(path * ".tmp", prev, rows)
    isfile(path * ".tmp") && mv(path * ".tmp", path; force=true)
end

conds = NamedTuple[]
curves = NamedTuple[]
for (k, (label, sx, sy, D, ρ, bulk, layers, times, te)) in enumerate(conditions(TASK))
    label in PREV.labels && continue
    seeds = 1000 * (OFFSET[TASK] + k) .+ (1:NSEEDS)
    cmin = minimum(cbar5(sh, h, sx, sy) for (_, sh, h) in layers)
    rmax = maximum(pkg_max_rate(sh, h, ρ, sx, sy, D) for (_, sh, h) in layers)
    τ = min(0.135 * cmin^2, 0.005 / rmax)
    models = vcat([("conventional", "conventional", 0.0)], layers)
    t0 = time()
    r = run_condition_outside(R=R1, sx=sx, sy=sy, D=D, rho=ρ, R2_bulk=bulk, models=models, times=times, tau=τ, nspins=NSPINS, seeds=seeds)
    push!(conds, NamedTuple{CKEYS}((TASK, label, sx, sy, gapof(sx, sy), D, ρ, bulk, r.tau, cmin, rmax, rmax * r.tau, te, times[end],
        r.identical_trajectories, time() - t0, NSPINS, NSEEDS, first(seeds))))
    for (name, shape, h) in models, (j, t) in enumerate(times)
        lay = shape != "conventional"
        push!(curves, NamedTuple{KKEYS}((label, name, shape, h, lay ? cbar5(shape, h, sx, sy) : 0.0, r.kept[name], t,
            r.S_mean[name][j], r.S_sem[name][j], r.S_ann[name][j],
            lay ? r.dS_mean[name][j] : NaN, lay ? r.dS_sem[name][j] : NaN, lay ? r.dS_ann[name][j] : NaN)))
    end
    merge_atomic(KFILE, PREV.kprev, curves)
    merge_atomic(CFILE, PREV.cprev, conds)
    @printf("%s %-20s τ = %.3e (c̄_min %.3f, ΔR₂max·τ %.4f), %d readouts to %.4g ms, identical %s, %.0f s\n", TASK, label, r.tau, cmin,
        rmax * r.tau, length(times), times[end], r.identical_trajectories, time() - t0)
    flush(stdout)
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "task" => TASK,
    "n_conditions" => length(conds) + length(PREV.labels), "n_resumed" => length(PREV.labels),
    "all_identical_new" => all(c -> c.identical_trajectories, conds), "nspins" => NSPINS, "nseeds" => NSEEDS))
println("saved to $OUT")
