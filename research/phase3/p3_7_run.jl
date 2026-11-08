# Phase 3.7a driver: paired Monte Carlo conditions inside single cylinders (primary), with the exact radial continuum
# alongside (support). TASK = fig9a | params | map | ddep (plan 2026-11-08-comparison-inside-cylinders.md, Tasks 3–6).
# Per condition: all models share τ = min(0.135·c̄_min², 0.005/ΔR2(0)_max) (comparison timestep rule; c̄_min over the
# layer models), seeds 1–10, 2e4 spins in the box [−R, R]³. Readouts: 0:TE*/10:1.5·TE* with TE* = continuum conventional
# T₂ of the condition (ddep: of D = 3); map: 41 readouts to t_end = ln 100/λ₁,conv.
# Output (rewritten after every condition, so partial results survive): results/phase3/p3_7_<task>/conditions.csv,
# curves.csv, summary.json.
# Resume: conditions whose label is already in conditions.csv are skipped and their rows kept (RESUME=0 starts afresh).
# Run: TASK=fig9a julia --project=research/baseline -t 8 research/phase3/p3_7_run.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
include(joinpath(@__DIR__, "analysis.jl"))
include(joinpath(@__DIR__, "p3_7_mc.jl"))
const TASK = get(ENV, "TASK", "fig9a")
const OUT = phase3_outdir("p3_7_" * TASK)
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const RHO_REF, D_REF3, R2B = 1 / 140, 3.0, 1 / 2000

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
rate0(shape, ρ, h, R) = ρ / (h * radial_kappa(shape, h, R)) * radial_g(shape, 0.0, radial_support(shape, h, R) / h)
same_centroid(cb, R) = [("$(s)_c$(cb)", s, h_for_cbar(s, cb, R)) for s in ("step", "linear", "exponential")]
same_h(h) = [("$(s)_h$(h)", s, h) for s in ("step", "linear", "exponential")]
tecont(R, D, ρ, bulk) = te_star(t -> radial_signal(R=R, D=D, times=[t], R2_bulk=bulk, surface=ρ)[1]; t_hi=50 / (bulk + 2ρ / R))

# conditions: (label, R, D, ρ, bulk, layer models, readout mode, TE for readouts or nothing)
function conditions(task)
    c = []
    if task == "fig9a"
        for R in (0.5, 1.0, 2.0)
            push!(c, ("R$(R)_c0.1", R, D_REF3, RHO_REF, R2B, same_centroid(0.1, R), :te, nothing))
        end
        push!(c, ("R1.0_h0.2", 1.0, D_REF3, RHO_REF, R2B, same_h(0.2), :te, nothing))
    elseif task == "params"
        for D in (0.5, 1.0, 2.0)
            push!(c, ("D$(D)", 1.0, D, RHO_REF, R2B, same_centroid(0.1, 1.0), :te, nothing))
        end
        for ρ in (1 / 280, 1 / 70, 1 / 35, 0.1)
            push!(c, ("rho$(round(ρ, sigdigits=4))", 1.0, D_REF3, ρ, R2B, same_centroid(0.1, 1.0), :te, nothing))
        end
        for cb in (0.05, 0.25)
            push!(c, ("c$(cb)", 1.0, D_REF3, RHO_REF, R2B, same_centroid(cb, 1.0), :te, nothing))
        end
    elseif task == "map"
        for x in (0.1, 1.0, 10.0)                                    # ρR/D at R = 1, D = 3
            models = [("$(s)_hR$(hr)", s, hr) for s in ("step", "linear", "exponential") for hr in (0.1, 0.25, 0.5)]
            push!(c, ("rhoRD$(x)", 1.0, D_REF3, x * D_REF3, 0.0, models, :map, nothing))
        end
    elseif task == "ddep"
        for ρ in (0.1, 1.0), D in (0.5, 1.0, 2.0, 3.0)
            models = [("$(s)_hR$(hr)", s, hr) for s in ("step", "linear", "exponential") for hr in (0.25, 0.5, 1.0)]
            push!(c, ("rho$(ρ)_D$(D)", 1.0, D, ρ, R2B, models, :te, tecont(1.0, D_REF3, ρ, R2B)))
        end
    else
        error("unknown TASK $task")
    end
    return c
end

const CFILE, KFILE = joinpath(OUT, "conditions.csv"), joinpath(OUT, "curves.csv")
const RESUME = get(ENV, "RESUME", "1") != "0"
const PREV_C = RESUME ? read_done(CFILE) : read_done("")
const PREV_K = RESUME ? read_done(KFILE) : read_done("")
isempty(PREV_C.labels) || println("resuming: skipping ", join(sort(collect(PREV_C.labels)), ", "))
conds = NamedTuple[]
curves = NamedTuple[]
for (label, R, D, ρ, bulk, layers, mode, te_fixed) in conditions(TASK)
    label in PREV_C.labels && continue
    cmin = minimum(cbar(s, h, R) for (_, s, h) in layers)
    rmax = maximum(rate0(s, ρ, h, R) for (_, s, h) in layers)
    τ = min(0.135 * cmin^2, 0.005 / rmax)
    if mode == :map
        λ1 = radial_rates(R=R, D=D, surface=ρ)[1]
        te = NaN
        times = collect(range(0.0, log(100) / λ1, length=41))
    else
        te = isnothing(te_fixed) ? tecont(R, D, ρ, bulk) : te_fixed
        times = collect(range(0.0, 1.5te, length=16))           # step te/10; te itself is the 11th readout
    end
    models = vcat([("conventional", "conventional", 0.0)], layers)
    t0 = time()
    r = run_condition(R=R, D=D, rho=ρ, R2_bulk=bulk, models=models, times=times, tau=τ, nspins=NSPINS)
    push!(conds, (task=TASK, label=label, R=R, D=D, rho=ρ, R2_bulk=bulk, tau=r.tau, cbar_min=cmin, rate0_max=rmax, te=te,
        t_end=times[end], identical_trajectories=r.identical_trajectories, runtime_s=time() - t0, nspins=NSPINS))
    for (name, shape, h) in models, (k, t) in enumerate(times)
        haslayer = shape != "conventional"
        push!(curves, (label=label, model=name, shape=shape, h=h, cbar=haslayer ? cbar(shape, h, R) : 0.0, t=t,
            S_mean=r.S_mean[name][k], S_sem=r.S_sem[name][k], S_cont=r.S_cont[name][k],
            dS_mean=haslayer ? r.dS_mean[name][k] : NaN, dS_sem=haslayer ? r.dS_sem[name][k] : NaN, dS_cont=haslayer ? r.dS_cont[name][k] : NaN))
    end
    write_merged(CFILE, PREV_C, conds)
    write_merged(KFILE, PREV_K, curves)
    @printf("%s %-12s τ = %.3e, %d readouts to %.3g ms, identical %s, %.0f s\n", TASK, label, r.tau, length(times), times[end],
        r.identical_trajectories, time() - t0)
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "task" => TASK, "n_conditions" => length(conds) + length(PREV_C.labels),
    "n_resumed" => length(PREV_C.labels), "all_identical" => all(c -> c.identical_trajectories, conds) &&
    all(l -> occursin(",true,", l), PREV_C.lines), "nspins" => NSPINS))
println("saved to $OUT")
