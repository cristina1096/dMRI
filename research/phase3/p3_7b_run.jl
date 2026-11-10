# Phase 3.7b driver: paired Monte Carlo outside cylinders in square arrays (primary) with the annulus continuum
# alongside (support), plus the one inside run that the summed signal needs.
# TASK = validate | inside | density | map | ddep (plan 2026-11-09-comparison-arrays-outside.md, Tasks 3–6).
# Per condition: τ = min(0.135·c̄_min², 0.005/ΔR₂(0)_max) over its layer models (comparison rule, D = 3 value for all D);
# seeds = 1000·k + (1:10) with k = task offset + condition index (conditions independent; 3.7a Group B review I1);
# 2·10⁴ spins per seed. R = 1 throughout; outside runs set only the outer layer.
# Group C trim (user decision 2026-11-10): ddep uses h/R ∈ {0.25, 0.5} (not 0.1, whose linear c̄ = 0.033 forced τ = 1.5e-4);
# the map at s = 4, ρR/D ∈ {0.1, 10} uses 5 seeds (the longest windows); the seed count is recorded per condition.
# Output (results/phase3/p3_7b_<task>/): curves.csv then conditions.csv, both via temporary file + rename after every
# condition; resume skips labels present in both (RESUME=0 starts afresh); summary.json at the end.
# Run: TASK=validate julia --project=research/baseline -t 8 research/phase3/p3_7b_run.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
include(joinpath(@__DIR__, "continuum_annulus.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
include(joinpath(@__DIR__, "p3_7_mc.jl"))
include(joinpath(@__DIR__, "p3_7_shared.jl"))
include(joinpath(@__DIR__, "p3_7b_mc.jl"))
const TASK = get(ENV, "TASK", "validate")
const OUT = phase3_outdir("p3_7b_" * TASK)
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const NSEEDS = parse(Int, get(ENV, "NSEEDS", "10"))
const R1, D3, RHO_REF, R2B, CB = 1.0, 3.0, 1 / 140, 1 / 2000, 0.1
const SHAPES = ("step", "linear", "exponential")
const OFFSET = Dict("validate" => 0, "inside" => 100, "density" => 200, "map" => 300, "ddep" => 400)

gapof(s) = s - 2R1
fof(s) = π * R1^2 / s^2
"ΔR₂ at the wall (d = 0) of an outer layer: (ρ/(hκ_out))·g(0)."
rate0_out(shape, ρ, h, L) = ρ / (h * (1 + outer_centroid(shape, h, L) / R1)) * (shape == "step" ? 1.0 : shape == "linear" ? 2.0 : 1 / -expm1(-L / h))
rate0_in(shape, ρ, h) = ρ / (h * radial_kappa(shape, h, R1)) * radial_g(shape, 0.0, radial_support(shape, h, R1) / h)
support_out(shape, h, s) = shape == "exponential" ? gapof(s) : h
cbar_out(shape, h, s) = outer_centroid(shape, h, support_out(shape, h, s))
outer_same_centroid(s) = [("$(sh)_c$(CB)", sh, h_for_cbar_out(sh, CB, gapof(s))) for sh in SHAPES]
inner_same_centroid() = [("$(sh)_c$(CB)", sh, h_for_cbar(sh, CB, R1)) for sh in SHAPES]
outer_hR(hrs) = [("$(sh)_hR$(hr)", sh, hr * R1) for sh in SHAPES for hr in hrs]
area_out(s) = π * (krogh_radius(s)^2 - R1^2)
"Annulus conventional T₂ (S = e⁻¹) outside, with bulk."
te_out(s, D, ρ, bulk) = te_star(annulus_signal_fun(R=R1, R_out=krogh_radius(s), D=D, surface=ρ, bulk=bulk);
                                t_hi=50 / (bulk + ρ * 2π * R1 / area_out(s)))
"TE where the summed conventional signal f·S_in(t) + (1 − f)·S_out(t) = e⁻¹ (radial continuum inside, annulus outside)."
function te_sum(s, D, ρ, bulk)
    f = fof(s)
    Sin = signal_fun(R=R1, D=D, surface=ρ, bulk=bulk)
    Sout = annulus_signal_fun(R=R1, R_out=krogh_radius(s), D=D, surface=ρ, bulk=bulk)
    return te_star(t -> f * Sin(t) + (1 - f) * Sout(t); t_hi=50 / (bulk + ρ * 2π * R1 / area_out(s)))
end
const DENSITY_S = (2.3, 2.5, 3.0, 4.0)
density_grid(s) = collect(range(0.0, 1.5 * te_sum(s, D3, RHO_REF, R2B), length=16))

# conditions: (label, compartment, s, D, ρ, bulk, layer models, times, te)
function conditions(task)
    c = []
    if task == "validate"
        for s in (2.5, 4.0), ρ in (0.02, 0.2)
            λ1 = annulus_rates(R=R1, R_out=krogh_radius(s), D=D3, surface=ρ)[1]
            times = ρ == 0.02 ? collect(0.0:5.0:50.0) : collect(range(0.0, log(100) / λ1, length=41))
            push!(c, ("s$(s)_rho$(ρ)", :outside, s, D3, ρ, 0.0, outer_same_centroid(s), times, NaN))
        end
    elseif task == "inside"
        times = sort(unique(vcat([density_grid(s) for s in DENSITY_S]...)))
        push!(c, ("inside_R1_c$(CB)", :inside, NaN, D3, RHO_REF, R2B, inner_same_centroid(), times, NaN))
    elseif task == "density"
        for s in DENSITY_S
            g = density_grid(s)
            push!(c, ("s$(s)", :outside, s, D3, RHO_REF, R2B, outer_same_centroid(s), g, g[11]))
        end
    elseif task == "map"
        for s in (2.5, 4.0), x in (0.1, 1.0, 10.0)
            ρ = x * D3 / R1
            λ1 = annulus_rates(R=R1, R_out=krogh_radius(s), D=D3, surface=ρ)[1]
            push!(c, ("s$(s)_rhoRD$(x)", :outside, s, D3, ρ, 0.0, outer_hR((0.1, 0.25, 0.5)),
                      collect(range(0.0, log(100) / λ1, length=41)), NaN))
        end
    elseif task == "ddep"
        s = 2.5
        for ρ in (0.1, 1.0)
            te = te_out(s, D3, ρ, R2B)
            for D in (0.5, 1.0, 2.0, 3.0)
                push!(c, ("s$(s)_rho$(ρ)_D$(D)", :outside, s, D, ρ, R2B, outer_hR((0.25, 0.5)),
                          collect(range(0.0, 1.5te, length=16)), te))
            end
        end
    else
        error("unknown TASK $task")
    end
    return c
end

const CKEYS = (:task, :label, :compartment, :s, :f, :gap, :D, :rho, :R2_bulk, :tau, :cbar_min, :rate0_max, :te, :t_end,
               :identical_trajectories, :runtime_s, :nspins, :nseeds, :seed_first)
const KKEYS = (:label, :model, :shape, :h, :cbar, :kept, :t, :S_mean, :S_sem, :S_sup, :dS_mean, :dS_sem, :dS_sup)
const CFILE, KFILE = joinpath(OUT, "conditions.csv"), joinpath(OUT, "curves.csv")
const RESUME = get(ENV, "RESUME", "1") != "0"
const PREV = RESUME ? read_done_both(CFILE, KFILE; cheader=join(CKEYS, ","), kheader=join(KKEYS, ",")) :
                      read_done_both("", ""; cheader="", kheader="")
isempty(PREV.labels) || println("resuming: skipping ", join(sort(collect(PREV.labels)), ", "))
"Previous rows plus new rows, via a temporary file and rename."
function merge_atomic(path, prev, rows)
    write_merged(path * ".tmp", prev, rows)
    isfile(path * ".tmp") && mv(path * ".tmp", path; force=true)
end

conds = NamedTuple[]
curves = NamedTuple[]
for (k, (label, comp, s, D, ρ, bulk, layers, times, te)) in enumerate(conditions(TASK))
    label in PREV.labels && continue
    nseeds = (TASK == "map" && s == 4.0 && any(x -> isapprox(ρ * R1 / D3, x; rtol=1e-9), (0.1, 10.0))) ? min(5, NSEEDS) : NSEEDS   # isapprox: ρR/D = 0.1000…02 in floating point
    seeds = 1000 * (OFFSET[TASK] + k) .+ (1:nseeds)
    cb(sh, h) = comp == :inside ? cbar(sh, h, R1) : cbar_out(sh, h, s)
    cmin = minimum(cb(sh, h) for (_, sh, h) in layers)
    rmax = maximum(comp == :inside ? rate0_in(sh, ρ, h) : rate0_out(sh, ρ, h, support_out(sh, h, s)) for (_, sh, h) in layers)
    τ = min(0.135 * cmin^2, 0.005 / rmax)
    models = vcat([("conventional", "conventional", 0.0)], layers)
    t0 = time()
    if comp == :inside
        r = run_condition(R=R1, D=D, rho=ρ, R2_bulk=bulk, models=models, times=times, tau=τ, nspins=NSPINS, seeds=seeds)
        Ssup, dSsup, kept = r.S_cont, r.dS_cont, Dict(n => 1.0 for (n, _, _) in models)
    else
        r = run_condition_outside(R=R1, s=s, D=D, rho=ρ, R2_bulk=bulk, models=models, times=times, tau=τ, nspins=NSPINS, seeds=seeds)
        Ssup, dSsup, kept = r.S_ann, r.dS_ann, r.kept
    end
    push!(conds, NamedTuple{CKEYS}((TASK, label, String(comp), s, comp == :inside ? NaN : fof(s), comp == :inside ? NaN : gapof(s),
        D, ρ, bulk, r.tau, cmin, rmax, te, times[end], r.identical_trajectories, time() - t0, NSPINS, nseeds, first(seeds))))
    for (name, shape, h) in models, (j, t) in enumerate(times)
        lay = shape != "conventional"
        push!(curves, NamedTuple{KKEYS}((label, name, shape, h, lay ? cb(shape, h) : 0.0, kept[name], t,
            r.S_mean[name][j], r.S_sem[name][j], Ssup[name][j],
            lay ? r.dS_mean[name][j] : NaN, lay ? r.dS_sem[name][j] : NaN, lay ? dSsup[name][j] : NaN)))
    end
    merge_atomic(KFILE, PREV.kprev, curves)          # curves first: a crash before conditions leaves the label not done
    merge_atomic(CFILE, PREV.cprev, conds)
    @printf("%s %-18s τ = %.3e, %d readouts to %.4g ms, identical %s, %.0f s\n", TASK, label, r.tau, length(times), times[end],
        r.identical_trajectories, time() - t0)
    flush(stdout)
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "task" => TASK,
    "n_conditions" => length(conds) + length(PREV.labels), "n_resumed" => length(PREV.labels),
    "all_identical_new" => all(c -> c.identical_trajectories, conds), "nspins" => NSPINS, "nseeds" => NSEEDS))
println("saved to $OUT")
