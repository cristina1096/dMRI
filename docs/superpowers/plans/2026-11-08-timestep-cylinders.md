# Timestep on cylinders (Phase 3.6: P3.6.1, P3.6.2, P3.6.4) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure how the Monte Carlo timestep biases results on cylinders, with the exact radial continuum as the reference inside a cylinder. Find what controls the bias: the per-step exponent ΔR₂(0)τ, the step length relative to the layer √(2Dτ)/h, or curvature √(Dτ)/R. Separate the error in a single model's decay from the error in the small difference between models at fixed ρ (the Review C1 problem). Then give a documented timestep rule (P3.6.2) for the user to decide on.

**Architecture:**
- **Research code only;** the package is unchanged.
- **Scans:** one parametrised scan script (`SCAN = thin | profiles | radius | rate`) runs paired Monte Carlo inside a single cylinder at a list of τ.
  - Every model of a scan shares τ, so trajectories are identical and ΔS = S_layer − S_conv is paired.
  - Results go to CSVs with the same columns as P3.4.6, so earlier runs at the same settings are reused rather than repeated.
- **Analysis:** an analysis script fits error laws with tested helpers (`research/phase3/analysis.jl`) and writes figures and a threshold table.
- **Decision:** whether to add a package option is the user's decision after the record (P3.6.2), not part of this plan.

**Tech Stack:** Julia 1.12, research environment `research/baseline` (MCMRSimulator, CairoMakie), the radial continuum (`research/phase3/continuum_radial.jl`).

**Spec:**
- `research/project-progress-tracker.md` 3.6: P3.6.1, P3.6.2, P3.6.4.
- Evidence so far:
  - P3.4.6 (step): h = 0.05 still fails at τ = 6.25e-4;
  - P3.4.6d (linear, exponential): |ΔS| 2.5–12 % low at τ = 2.5e-3;
  - P3.4.8 outside: thin-end shifts of about 3σ;
  - Phase 2 P2.6: the rate error is controlled by ΔR₂(0)τ, and the default is c = 0.005;
  - Review C1 (Phase 2 record): the difference error matters more than the rate error.

## Global Constraints

- **Branch.** Work only on branch `cc/near-surface-r2`. Never modify or commit on other branches. Do not push; the user pushes.
- **Commits.** Every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **No package change.** No change to `src/`. If a constraint should change, the record ends with options for the user (P3.6.2).
- **Reference.** Inside single cylinders the radial continuum is exact, so errors are measured against it. Outside there is no reference, so earlier outside data are reported only as τ trends.
- **Settings for reuse.** R = 1 unless scanned, D = 3, readouts 0:5:50 ms (`rate` scan: shorter, see Task 2), seeds 1–10, 2·10⁴ spins per seed in the box [−R, R]³, ρ = 0.02 unless scanned, conventional θ = ρ√(π/D).
- **Provenance.** Commit scripts **before** running them, so `git_commit` in the results matches.
- **Long runs.** Run detached (nohup chain plus a waiter) because the background-tool limit is 2 h. Wait for a full result before committing results.
- **Records.** Design / Method / Result / Conclusion with formulas and figures. Signs and laws are recorded as found, not assumed.
- **Commands.**
  - Research tests: `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`
  - Scans: `SCAN=<name> julia --project=research/baseline -t 8 research/phase3/p3_6_1_tau_scan.jl`
- **Plan end.** A fresh subagent reviewer (most capable model) reviews the whole branch (user rule).

## Review Focus

1. **Pairing across τ.** Runs at different τ are independent: different trajectories, though the same seeds. Errors are measured per τ against the continuum, never as paired differences across τ.
2. **Mixing reused rows with new rows.** Reused rows must have identical settings (R, ρ, spins, seeds, readouts). The analysis checks this from each source's summary and errors on a mismatch.
3. **Small denominators.** Relative errors of ΔS are meaningless where |ΔS_cont| is comparable to its MC SEM. Such rows are flagged and excluded from fits, not silently kept.
4. **Power-law fits on few points.** Report the fitted exponent with its standard error, and do not extrapolate beyond the scanned range.
5. **The `rate` scan uses a shorter T.** R = −ln S(T)/T must use the same T in the MC and the continuum.

---

### Task 1: Analysis helpers with tests

**Files:**
- Create: `research/phase3/analysis.jl`
- Create: `research/phase3/test/test_analysis.jl`
- Modify: `research/phase3/test/runtests.jl` (include `analysis.jl`)

**Interfaces — Produces:**
- `fit_power(x, y; w=nothing)` → `(a, p, se_p)` for y ≈ a·x^p, by least squares on log y vs log x (weights optional; requires y > 0);
- `threshold_at(x, y, level)` → the largest x in the scanned range with y ≤ level, by log-log interpolation between neighbours; `NaN` if none;
- `rel_error(mc, cont)` = (mc − cont)/cont.

- [ ] **Step 1: Write the failing tests** (`research/phase3/test/test_analysis.jl`)

```julia
@testset "analysis helpers" begin
    x = [0.1, 0.2, 0.4, 0.8]
    y = 3.0 .* x .^ 1.5
    (a, p, se) = fit_power(x, y)
    @test a ≈ 3.0 rtol=1e-10
    @test p ≈ 1.5 rtol=1e-10
    @test se < 1e-8
    noisy = y .* [1.01, 0.99, 1.02, 0.98]
    (_, p2, se2) = fit_power(x, noisy)
    @test abs(p2 - 1.5) < 0.05 && se2 > 0
    @test threshold_at(x, y, 3.0 * 0.4^1.5) ≈ 0.4 rtol=1e-10
    @test threshold_at(x, y, 3.0 * 0.3^1.5) ≈ 0.3 rtol=1e-10          # log-log interpolation is exact for a power law
    @test isnan(threshold_at(x, y, 1e-6))
    @test threshold_at(x, y, 100.0) == 0.8                              # every point below: the largest scanned x
    @test rel_error(1.1, 1.0) ≈ 0.1
    @test_throws ArgumentError fit_power([0.1, 0.2], [1.0, -1.0])
end
```

Add `include(joinpath(@__DIR__, "..", "analysis.jl"))` after the continuum include in `runtests.jl`.

- [ ] **Step 2: Run to verify they fail** — `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl` → FAIL (`could not open file … analysis.jl`).

- [ ] **Step 3: Implement** (`research/phase3/analysis.jl`)

```julia
# Phase 3 analysis helpers (P3.6): power-law fits and thresholds on scanned error curves.

"""
    fit_power(x, y; w=nothing)

Least squares of log y = log a + p log x (weights `w` on the log points); returns (a, p, se_p), with se_p the standard
error of the slope from the residual variance (n − 2 degrees of freedom; 0 for n = 2). Requires x, y > 0.
"""
function fit_power(x, y; w=nothing)
    (all(>(0), x) && all(>(0), y)) || throw(ArgumentError("fit_power needs x, y > 0"))
    lx, ly = log.(x), log.(y)
    W = isnothing(w) ? ones(length(x)) : collect(Float64, w)
    sw = sum(W)
    mx, my = sum(W .* lx) / sw, sum(W .* ly) / sw
    sxx = sum(W .* (lx .- mx) .^ 2)
    p = sum(W .* (lx .- mx) .* (ly .- my)) / sxx
    a = exp(my - p * mx)
    n = length(x)
    res = ly .- (log(a) .+ p .* lx)
    se = n > 2 ? sqrt(sum(W .* res .^ 2) / (n - 2) / sxx) : 0.0
    return (a, p, se)
end

"""
    threshold_at(x, y, level)

Largest x in the scanned range with y ≤ level, assuming y increases with x: the last point at or below `level`, then
log-log interpolation towards the next point if that one is above. NaN if the first point is already above `level`.
"""
function threshold_at(x, y, level)
    o = sortperm(x)
    xs, ys = x[o], y[o]
    ys[1] > level && return NaN
    for k in 1:length(xs) - 1
        if ys[k + 1] > level
            f = (log(level) - log(ys[k])) / (log(ys[k + 1]) - log(ys[k]))
            return exp(log(xs[k]) + f * (log(xs[k + 1]) - log(xs[k])))
        end
    end
    return xs[end]
end

"Relative error (mc − cont)/cont."
rel_error(mc, cont) = (mc - cont) / cont
```

- [ ] **Step 4: Run to verify they pass** — research suite PASS.

- [ ] **Step 5: Commit**

```bash
git add research/phase3/analysis.jl research/phase3/test/test_analysis.jl research/phase3/test/runtests.jl
git commit -m "Phase 3 analysis helpers: power-law fit with slope error, log-log threshold, relative error

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: The τ scan script (P3.6.1 data)

**Files:**
- Create: `research/phase3/p3_6_1_tau_scan.jl` → `research/results/phase3/p3_6_1_tau_scan_<scan>/` (`comparison.csv`, `summary.json`)

**Scans** (inside a single cylinder; all models of one scan share τ, so they are paired):

| SCAN | R | models (ρ = 0.02 unless stated) | τ (ms) run here | reused from |
|---|---|---|---|---|
| `thin` | 1 | conventional; step h = 0.05, 0.1, 0.2 | 5e-3, 1.25e-3 (all four); 6.25e-4 (conventional, h = 0.2); 3.125e-4 (conventional, h = 0.05) | P3.4.6: 1e-2, 2.5e-3 (all), 6.25e-4 (conventional, 0.05, 0.1) |
| `profiles` | 1 | conventional; linear h = 0.3; exponential h = 0.1 | 1.25e-3, 6.25e-4 | P3.4.6d: 1e-2, 2.5e-3 |
| `radius` | 0.5 and 2 | conventional; step h = 0.1R | 1e-2, 2.5e-3, 6.25e-4 | — |
| `rate` | 1 | step h = 0.1 at ΔR₂(0) = 0.5 and 2 ms⁻¹ (ρ = ΔR₂(0)·h·κ, κ = 0.95); T = 2R/(2ρ) ms, so E[X] = 2 | τ with ΔR₂(0)τ = 0.05, 0.02, 0.005, 0.00125 | — |

The `rate` scan measures the error of a single model's decay, R = −ln S(T)/T against the continuum. The others measure both S and the paired ΔS against the continuum.

- [ ] **Step 1: Write the script**

```julia
# P3.6.1: Monte Carlo timestep scan inside a single cylinder against the exact radial continuum.
# SCAN = thin | profiles | radius | rate (table in the plan). For every (configuration, τ): all models share τ, so
# trajectories are identical (checked) and ΔS = S_layer − S_conv is paired per seed. Recorded per model and τ:
#   S(T) MC ± SEM, continuum, z; ΔS(T) MC ± paired SEM, continuum, z; R = −ln S(T)/T MC and continuum;
#   ΔR2(0)τ, √(2Dτ)/h, √(Dτ)/R. T = 50 ms, except the `rate` scan (T = 2R/(2ρ), so E[−ln M] = 2).
# Settings match P3.4.6 / P3.4.6d (seeds 1–10, 2e4 spins in the box [−R, R]³, ρ = 0.02, D = 3) so their rows are reused.
# Run: SCAN=thin julia --project=research/baseline -t 8 research/phase3/p3_6_1_tau_scan.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
const SCAN = get(ENV, "SCAN", "thin")
const OUT = phase3_outdir("p3_6_1_tau_scan_" * SCAN)
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
θ(ρ) = ρ * sqrt(π / D_REF)
κstep(h, R) = 1 - h / (2R)
cyl(R, ρ, shape, h) = shape == "conventional" ? Cylinders(radius=R, surface_relaxation=θ(ρ)) :
                       Cylinders(radius=R, layer_rho_inner=ρ, layer_h_inner=h, layer_shape=shape)
rate0(shape, ρ, h, R) = shape == "conventional" ? 0.0 : ρ / (h * radial_kappa(shape, h, R)) * radial_g(shape, 0.0, radial_support(shape, h, R) / h)

# configurations: (label, R, ρ, T, models::Vector{(name, shape, h)}, taus)
function configs(scan)
    m(shapes) = [(s == "conventional" ? "conventional" : "$(s)_h$(h)", s, h) for (s, h) in shapes]
    if scan == "thin"
        return [("R1", 1.0, 0.02, 50.0, m([("conventional", 0.0), ("step", 0.05), ("step", 0.1), ("step", 0.2)]), [5e-3, 1.25e-3]),
                ("R1_h0.2_tau_c", 1.0, 0.02, 50.0, m([("conventional", 0.0), ("step", 0.2)]), [6.25e-4]),
                ("R1_h0.05_tau_d", 1.0, 0.02, 50.0, m([("conventional", 0.0), ("step", 0.05)]), [3.125e-4])]
    elseif scan == "profiles"
        return [("R1", 1.0, 0.02, 50.0, m([("conventional", 0.0), ("linear", 0.3), ("exponential", 0.1)]), [1.25e-3, 6.25e-4])]
    elseif scan == "radius"
        return [("R$(R)", R, 0.02, 50.0, m([("conventional", 0.0), ("step", 0.1R)]), [1e-2, 2.5e-3, 6.25e-4]) for R in (0.5, 2.0)]
    elseif scan == "rate"
        out = []
        for r0 in (0.5, 2.0)
            ρ = r0 * 0.1 * κstep(0.1, 1.0)
            push!(out, ("rate$(r0)", 1.0, ρ, 2 * 1.0 / (2ρ), m([("step", 0.1)]), [c / r0 for c in (0.05, 0.02, 0.005, 0.00125)]))
        end
        return out
    end
    error("unknown SCAN $scan")
end

rows = NamedTuple[]
allidentical = true
for (label, R, ρ, T, models, taus) in configs(SCAN), τ in taus
    times = collect(0.0:T / 10:T)
    bb = mr.BoundingBox([-R, -R, -R], [R, R, R])
    runs = Dict(n => run_cylinders(cyl(R, ρ, s, h); region=:inside, R=R, bounding_box=bb, timestep=τ, times=times, nspins=NSPINS)
                for (n, s, h) in models)
    if haskey(runs, "conventional")
        ident = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs["conventional"].finals, runs[n].finals)) for (n, _, _) in models)
        global allidentical &= ident
    end
    Sconv = haskey(runs, "conventional") ? radial_signal(R=R, D=D_REF, times=[T], surface=ρ)[1] : NaN
    for (n, s, h) in models
        r = runs[n]
        Sc = s == "conventional" ? Sconv : radial_signal(R=R, D=D_REF, times=[T], layer=(s, ρ, h))[1]
        if s != "conventional" && haskey(runs, "conventional")
            d = r.S[:, end] .- runs["conventional"].S[:, end]
            dm, ds, dc = mean(d), std(d) / sqrt(length(d)), Sc - Sconv
        else
            dm, ds, dc = NaN, NaN, NaN
        end
        push!(rows, (scan=SCAN, config=label, model=n, shape=s, h=h, R=R, rho=ρ, T=T, tau=τ,
            S_mc=r.mean[end], S_sem=r.sem[end], S_cont=Sc, zS=(r.mean[end] - Sc) / r.sem[end],
            dS_mc=dm, dS_sem=ds, dS_cont=dc, zdS=(dm - dc) / ds,
            R_mc=-log(r.mean[end]) / T, R_cont=-log(Sc) / T,
            rate0_tau=rate0(s, ρ, h, R) * τ, step_over_h=s == "conventional" ? NaN : sqrt(2 * D_REF * τ) / h,
            step_over_R=sqrt(D_REF * τ) / R, nspins=NSPINS))
        @printf("%-14s %-16s τ = %.3e: S %.6f ± %.1e vs %.6f (z %+.2f); ΔS %+.2e ± %.1e vs %+.2e (z %+.2f)\n",
            label, n, τ, r.mean[end], r.sem[end], Sc, (r.mean[end] - Sc) / r.sem[end], dm, ds, dc, (dm - dc) / ds)
    end
end
write_csv(joinpath(OUT, "comparison.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "scan" => SCAN, "nspins" => NSPINS,
    "identical_trajectories" => allidentical, "settings" => Dict("D" => D_REF, "seeds" => collect(SEEDS_REF))))
println("saved to $OUT")
```

- [ ] **Step 2: Smoke test** (seconds; not committed): `SCAN=rate NSPINS=200 julia --project=research/baseline -t 8 research/phase3/p3_6_1_tau_scan.jl`. Expected: 8 lines, no error. Then delete `research/results/phase3/p3_6_1_tau_scan_rate`.

- [ ] **Step 3: Commit the script before running**

```bash
git add research/phase3/p3_6_1_tau_scan.jl
git commit -m "P3.6.1: parametrised τ scan inside a cylinder against the radial continuum

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 4: Run all scans detached** (≈ 6–7 h: `thin` ≈ 3 h, `profiles` ≈ 2 h, `radius` ≈ 1.5 h, `rate` ≈ 0.3 h), as a nohup chain writing logs and a done-marker, with a waiter `until [ -f done ]; do sleep 120; done` (restart the waiter if it hits its 2 h limit).

- [ ] **Step 5: Commit results** (only after all four logs show `saved to` and exit 0)

```bash
git add research/results/phase3/p3_6_1_tau_scan_*
git commit -m "P3.6.1: τ scans inside cylinders (thin, profiles, radius, rate) vs the radial continuum

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Analysis of the timestep error (P3.6.1 laws, P3.6.2 thresholds)

**Files:**
- Create: `research/phase3/p3_6_2_analysis.jl` → `research/results/phase3/p3_6_2_analysis/` (`errors.csv`, `fits.csv`, `thresholds.csv`, `summary.json`, `p3_6_2_errors.png`)

**Method:**
- **Assemble:** the new scan rows, plus reused rows:
  - P3.4.6: `comparison.csv` (τ_a, τ_b) and `comparison_tau16.csv` (τ_c), for conventional and step h = 0.05–1.0;
  - P3.4.6d: `comparison.csv`.

  Reused rows are mapped to the same columns. The continuum values are recomputed with the current solver (face-centred Robin) for every row, so all rows share one reference. The analysis checks that every source has R = 1, ρ = 0.02, 2·10⁴ spins and seeds 1–10, and errors otherwise (Review Focus 2).
- **Error measures, per model and τ:**
  - **Rate error:** e_R = (R_mc − R_cont)/R_cont with R = −ln S(T)/T.
  - **Difference error:** e_ΔS = (ΔS_mc − ΔS_cont)/|ΔS_cont|, kept only where |ΔS_cont| ≥ 5·SEM(ΔS) (Review Focus 3).
- **Fits:**
  - e_R against ΔR₂(0)τ, from the `rate` scan and the step rows;
  - |e_ΔS| against √(2Dτ)/h and against ΔR₂(0)τ, separately for step, linear and exponential;
  - the conventional model's S error against √(Dτ)/R (thin and radius scans).

  Use `fit_power` on points with |e| > 2σ_e (σ_e from the MC SEM); report the exponent ± SE (Review Focus 4).
- **Thresholds:** for |e_ΔS| = 10 %, 5 % and 2 %, `threshold_at` gives the largest √(2Dτ)/h, converted to a rule τ ≤ k·h²/(2D) with k = threshold². The same for e_R at 1 % and 0.1 % in terms of ΔR₂(0)τ, to compare with c = 0.005.

- [ ] **Step 1: Write the script** (uses `analysis.jl`, `continuum_radial.jl`, `harness.jl`).

  It reads the CSVs above and builds one `errors.csv` with columns `source, scan, model, shape, h, R, rho, T, tau, rate0_tau, step_over_h, step_over_R, e_R, e_dS, e_dS_sigma, used_in_fit`. It writes:
  - `fits.csv`: one row per (measure, profile, variable): `a, p, se_p, n_points`;
  - `thresholds.csv`: level, variable, threshold, rule;
  - `p3_6_2_errors.png`, 2 × 2 panels on log-log axes, with the fitted lines over the scanned range only:
    - (a) e_R against ΔR₂(0)τ;
    - (b) |e_ΔS| against √(2Dτ)/h for each profile;
    - (c) |e_ΔS| against ΔR₂(0)τ;
    - (d) the conventional |S error| against √(Dτ)/R.

  (The executor writes the code from this specification. Tested helpers do all the fitting and thresholding, and a unit test in Task 1 covers each.)

- [ ] **Step 2: Run** (seconds to minutes) → inspect the figure. Expected:
  - e_R scales roughly linearly with ΔR₂(0)τ, as on walls; record the exponent as found;
  - the threshold table is filled.

- [ ] **Step 3: Commit** the script first, then rerun, then commit the results (provenance).

---

### Task 4: Runtime summary (P3.6.4, no new runs)

- [ ] Assemble a table from existing results:
  - C-B5 (`p3_1_5_runtime`, including `runtime_many.csv`);
  - P3.3.4 (`p3_3_quadrature/runtime.csv`);
  - the walls cost (Phase 2 P2.6.8 table in `research/phase2_record.md`).

  Report µs per spin-step: no layer; step, linear and exponential; 1 cylinder against 400 per cell; walls against cylinders. Record the deferred cell-list optimisation as P3.6.4's open item. Record only, no script needed; the numbers are quoted from the files with their paths.

---

### Task 5: Record, P3.6.2 options, tracker

- [ ] **Step 1:** Write a record section "Timestep on cylinders (P3.6.1, P3.6.2, P3.6.4)" with Design / Method / Result / Conclusion. It includes the figure, the fits with errors, the threshold table, and the runtime table.
- [ ] **Step 2:** End the section with **P3.6.2 options for the user**, each with its cost from the measured runtime:
  - (a) keep c = 0.005 for single-model runs, if e_R ≤ 0.1 % holds;
  - (b) for comparison studies, a documented rule τ ≤ k·h_min²/(2D) from the thresholds;
  - (c) optionally, a package timestep option implementing (b) (a later plan);
  - (d) use the continuum instead of Monte Carlo wherever it exists.

  State which the data support. Do **not** implement (c).
- [ ] **Step 3:** Tracker:
  - P3.6.1 ☑ if the scans and fits are complete;
  - P3.6.2 ◐ ("options documented, user decision pending");
  - P3.6.4 ☑ for the table, with the cell list noted as open.

  Update the done count.
- [ ] **Step 4:** Research suite PASS, then commit.

---

## After the last task

Final whole-branch review by a fresh subagent on the most capable model (user rule). Give it this plan, the spec pointers, the Review Focus and the ledger. Fix Critical and Important findings, and ledger the minors.
