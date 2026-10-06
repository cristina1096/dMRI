# Layer Timestep Study (V1–V3) and Cost Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Measure how the layer model's signal depends on the timestep τ (V1), and how the converged timestep τ_conv scales with h and D (V2, V3), which gives c_h in τ ≤ c_h·h²/D (proposal Eq. 11). Then build that constraint into the simulator's timestep choice, measure the cost, and state which earlier τ = 1e-2 ms results are converged.

**Architecture:** The dense τ scans run on the 1D toy walk (`research/phase2/toy.jl`). It reproduces MCMR's motion between walls exactly (V11 cross-check, +1.8 SEM) and costs far less per step, so a scan down to τ = 1e-5 ms is affordable. MCMR itself is run at selected τ to confirm the toy curve. The observable is the effective rate R = −ln S(T)/T after T = 10 ms. τ_conv is the largest τ whose relative error in R, against the smallest-τ reference, stays within a tolerance, interpolated on a log scale between grid points. The constraint goes into `src/timesteps.jl` as a sixth option, `layer`, using the thinnest active layer.

**Tech Stack:** Julia 1.12, MCMRSimulator (this repo), `research/phase2/` (harness, toy), CairoMakie, `Test` stdlib.

**Spec:** `research/project-progress-tracker.md` section 2.6 (P2.6.1–P2.6.5, P2.6.8, P2.6.10–P2.6.12, and the decision point under V3); proposal §3.2 V (Eq. 11: τ_max = min(τ_current, c_h·h²/D); straight steps miss excursions into the layer, so large τ under-estimates layer time); the user's framing (no fast-diffusion arguments).

**Depends on:** the walls layer, the harness/trace and the step-verification plans (all done); B5 (old-code runtime vs τ).

## Global Constraints

- Work only on branch `cc/near-surface-r2`.
- **Configuration:** walls every w = 2 µm, layer on both faces, step profile, ΔR₂(0) = 0.1 ms⁻¹ unless a task says otherwise. Toy: T = 10 ms, 10 seeds; spins per seed set by the pilot (Task 2).
- **Error measure:** relative error of R = −ln S(T)/T against the reference (smallest τ). Tolerances **1%** (primary) and **0.1%** (reported; marked "noise-limited" where the reference's relative SEM exceeds 0.03%).
- The **sign of the large-τ bias** is measured and reported, not assumed. The proposal (§3.2 V) expects too little attenuation at large τ; the plan records whether that holds.
- **Out of scope** (recorded in the tracker as such):
  - P2.6.4 (profile dependence): needs the linear plan first. It will reuse this plan's tools.
  - P2.6.6, P2.6.7, P2.6.9 (integration tolerance, quadrature vs closed form): there's no numerical-integration path yet, because every implemented profile has a closed form.
  - P2.6.13 (`rd-walls-v1`): needs 2.4 and 2.5.
- Research tests: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`. Package tests: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `research/phase2/timestep.jl` | Create | `toy_rate`, `tau_conv_grid`, `tau_conv_interp`, `fit_c_h`, `log_taus` |
| `research/phase2/test/test_timestep.jl` | Create | unit tests for the above |
| `research/phase2/test/runtests.jl` | Modify | include `timestep.jl` when it exists |
| `research/phase2/p2_6_scan.jl` | Create | one τ scan on the toy for a given (h, ΔR₂(0), D); used by V1, V2, V3 |
| `research/phase2/p2_6_1_mcmr_check.jl` | Create | MCMR at selected τ vs the toy curve (h = 0.1) |
| `research/phase2/p2_6_2_scaling.jl` | Create | combines the scans: τ_conv vs h²/D, c_h, rate dependence, D dependence |
| `research/phase2/p2_6_8_cost.jl` | Create | runtime vs τ and vs h with the layer, noise cost, accuracy–cost figure |
| `src/geometries/internal/layers.jl`, `internal.jl` | Modify | `min_layer_h(geometry)` |
| `src/timesteps.jl` | Modify | new option `layer` (c_h) in `TimeStep` |
| `test/test_layer.jl` | Modify | timestep-constraint tests |
| `research/phase2_record.md`, `research/project-progress-tracker.md` | Modify | records, ticks |

## Review Focus

1. **Toy steps must fit T exactly.** With log-spaced τ, T/τ is not an integer. The toy runs round(T/τ) steps, so R must be computed with the actual time T_eff = round(T/τ)·τ, not T. Test: Task 1's `toy_rate` test with a τ that doesn't divide T.
2. **No τ fails:** if no τ in the grid exceeds the tolerance, τ_conv is "at least the largest τ", and must be reported as a lower bound, not as a value. Test: Task 1.
3. **A noisy reference:** when the reference's relative SEM is comparable to the tolerance, τ_conv is meaningless. Test: Task 1 (the `noise_limited` flag).
4. **The constraint must ignore sides with ρ = 0,** and must give Inf for a geometry without a layer, so the old timestep choice is unchanged. Test: Task 5.
5. **The verbose message** must say when the layer constraint is the binding one. Test: Task 5 (`@test_logs` with the expected text).

---

### Task 1: Timestep-analysis tools

**Files:**
- Create: `research/phase2/timestep.jl`, `research/phase2/test/test_timestep.jl`
- Modify: `research/phase2/test/runtests.jl`

**Interfaces:**
- Consumes: `toy_layer_signal(; h, rate, tau, T, nspins, seed, method, w, D)` (toy.jl).
- Produces:
  - `log_taus(lo, hi; per_decade=3)::Vector{Float64}`
  - `toy_rate(; h, rate, D, tau, T, nspins, seeds, method=:exact)`, which returns `(mean, sem, per_seed, T_eff)`, with R = −ln S(T_eff)/T_eff per seed
  - `tau_conv_grid(taus, R, R_ref; tol)`, which returns `(tau, all_pass)`
  - `tau_conv_interp(taus, R, R_ref; tol)`, which returns `(tau, lower_bound::Bool)`
  - `noise_limited(R_ref, R_ref_sem; tol)::Bool`, true when R_ref_sem/R_ref > tol/3
  - `fit_c_h(tau_convs, hs, Ds)`, which returns `(c_h, spread, ratios)`, with ratios τ_conv·D/h²

- [ ] **Step 1: Write the failing tests**

`research/phase2/test/test_timestep.jl`:

```julia
@testset "timestep tools" begin
    @testset "log_taus" begin
        t = log_taus(1e-5, 1e-1)
        @test t[1] ≈ 1e-5 && t[end] ≈ 1e-1
        @test length(t) == 13
        @test all(diff(log10.(t)) .≈ 1 / 3)
    end

    @testset "toy_rate uses the actual simulated time" begin
        # h = w/2: every spin decays at exactly `rate`, so R = rate for any τ, even when τ does not divide T
        r = toy_rate(; h=1.0, rate=0.1, D=3.0, tau=0.003, T=1.0, nspins=50, seeds=1:2)
        @test r.T_eff ≈ round(1.0 / 0.003) * 0.003
        @test all(isapprox.(r.per_seed, 0.1; rtol=1e-10))
    end

    @testset "tau_conv on synthetic data" begin
        taus = [1e-4, 1e-3, 1e-2, 1e-1]
        R = [1.000, 1.001, 1.02, 1.2]                 # relative errors 0, 0.1%, 2%, 20%
        g = tau_conv_grid(taus, R, 1.0; tol=0.01)
        @test g.tau == 1e-3 && !g.all_pass
        i = tau_conv_interp(taus, R, 1.0; tol=0.01)
        @test 1e-3 < i.tau < 1e-2 && !i.lower_bound
        # where |error| crosses 1% between 0.1% (at 1e-3) and 2% (at 1e-2), linearly in log τ
        @test log10(i.tau) ≈ -3 + (0.01 - 0.001) / (0.02 - 0.001) atol=1e-12
        n = tau_conv_interp(taus, [1.0, 1.0, 1.0, 1.005], 1.0; tol=0.01)
        @test n.tau == 1e-1 && n.lower_bound      # never fails: lower bound only
    end

    @testset "noise_limited and fit_c_h" begin
        @test noise_limited(1.0, 0.001; tol=0.001)
        @test !noise_limited(1.0, 0.0001; tol=0.001)
        f = fit_c_h([0.01, 0.04], [0.1, 0.2], [3.0, 3.0])
        @test f.c_h ≈ 3.0 && f.spread ≈ 0.0
    end
end
```

Add to `research/phase2/test/runtests.jl`, after the `toy.jl` line:

```julia
isfile(joinpath(@__DIR__, "..", "timestep.jl")) && include(joinpath(@__DIR__, "..", "timestep.jl"))
```

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: FAIL with `UndefVarError: log_taus not defined`.

- [ ] **Step 3: Create `research/phase2/timestep.jl`**

```julia
# Timestep-study tools (P2.6): toy-based rate per τ, converged timestep, c_h.
#
# Observable: R = −ln S(T_eff)/T_eff, the effective decay rate after T_eff = round(T/τ)·τ.
# τ_conv(tol): the largest τ such that every τ' ≤ τ in the grid has |R(τ') − R_ref| / R_ref ≤ tol,
# with R_ref = R at the smallest τ; interpolated linearly in log τ to where the error crosses tol.

"Log-spaced timesteps from lo to hi (ms), `per_decade` points per decade, both ends included."
log_taus(lo, hi; per_decade=3) = 10 .^ range(log10(lo), log10(hi), length=round(Int, per_decade * log10(hi / lo)) + 1)

"""
    toy_rate(; h, rate, D, tau, T, nspins, seeds, method=:exact)

Effective rate R = −ln S/T_eff per seed from the 1D toy (walls every W_REF, step layer on both faces),
with T_eff = round(T/τ)·τ the time actually simulated. Returns (mean, sem, per_seed, T_eff).
"""
function toy_rate(; h, rate, D, tau, T, nspins, seeds, method=:exact)
    T_eff = round(Int, T / tau) * tau
    Rs = [-log(mean(toy_layer_signal(; h=h, rate=rate, tau=tau, T=T, nspins=nspins, seed=s, method=method, D=D))) / T_eff
          for s in seeds]
    return (mean=mean(Rs), sem=length(Rs) > 1 ? std(Rs) / sqrt(length(Rs)) : NaN, per_seed=Rs, T_eff=T_eff)
end

"Largest grid τ such that all τ' ≤ τ have relative error ≤ tol. `all_pass` is true when no τ fails."
function tau_conv_grid(taus, R, R_ref; tol)
    order = sortperm(taus)
    best = taus[order[1]]
    for k in order
        abs(R[k] - R_ref) / R_ref <= tol || return (tau=best, all_pass=false)
        best = taus[k]
    end
    return (tau=best, all_pass=true)
end

"""
    tau_conv_interp(taus, R, R_ref; tol)

As `tau_conv_grid`, but interpolated: linear in log τ between the last passing and the first failing grid point,
to where the relative error reaches `tol`. If no τ fails, returns the largest τ with `lower_bound = true`.
"""
function tau_conv_interp(taus, R, R_ref; tol)
    order = sortperm(taus)
    err(k) = abs(R[k] - R_ref) / R_ref
    for (i, k) in enumerate(order)
        err(k) <= tol && continue
        i == 1 && return (tau=taus[k], lower_bound=false)
        kp = order[i - 1]
        f = (tol - err(kp)) / (err(k) - err(kp))
        return (tau=exp(log(taus[kp]) + f * (log(taus[k]) - log(taus[kp]))), lower_bound=false)
    end
    return (tau=taus[order[end]], lower_bound=true)
end

"True when the reference's relative SEM is more than a third of the tolerance (τ_conv not resolvable)."
noise_limited(R_ref, R_ref_sem; tol) = R_ref_sem / R_ref > tol / 3

"c_h from τ_conv = c_h·h²/D: mean of the ratios τ_conv·D/h², with their relative spread (std/mean)."
function fit_c_h(tau_convs, hs, Ds)
    ratios = tau_convs .* Ds ./ hs .^ 2
    return (c_h=mean(ratios), spread=length(ratios) > 1 ? std(ratios) / mean(ratios) : 0.0, ratios=ratios)
end
```

- [ ] **Step 4: Run the tests and check they pass**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: PASS (testset "timestep tools").

- [ ] **Step 5: Commit**

```bash
git add research/phase2/timestep.jl research/phase2/test/test_timestep.jl research/phase2/test/runtests.jl
git commit -m "Add timestep-study tools (toy rate, tau_conv, c_h)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: The scan script and pilot

**Files:**
- Create: `research/phase2/p2_6_scan.jl`

**Interfaces:**
- Consumes: `toy_rate`, `log_taus`, `tau_conv_interp`, `tau_conv_grid`, `noise_limited` (Task 1).
- Produces: `results/phase2/p2_6_scan_<label>/{scan.csv, summary.json, scan.png}` for a configuration given by environment variables: `H` (µm), `RATE` (1/ms), `D` (µm²/ms), `TAU_LO`, `TAU_HI` (ms), `NSPINS`, `T` (ms), `LABEL`. `summary.json` contains, per tolerance (0.01, 0.001), `tau_conv`, `lower_bound`, `noise_limited`, plus the large-τ bias sign and the timing per spin-step.

- [ ] **Step 1: Create the script**

```julia
# P2.6 τ scan on the 1D toy (exact layer rule, step profile, layer on both faces of walls every 2 µm).
#
# Environment: H (µm), RATE = ΔR2(0) (1/ms), D (µm²/ms), TAU_LO, TAU_HI (ms, 3 per decade), NSPINS, T (ms), LABEL.
# Observable: R = −ln S(T_eff)/T_eff per seed (seeds 1–10). Reference: the smallest τ.
# Reports τ_conv at relative tolerance 1 % and 0.1 % (interpolated in log τ), whether it is only a lower bound,
# whether the reference is too noisy for that tolerance, and the sign of the bias at the largest τ.
#
# Run (example, V1): H=0.1 RATE=0.1 D=3 TAU_LO=1e-5 TAU_HI=1e-1 NSPINS=100000 T=10 LABEL=V1_h0.1 \
#      julia --project=research/baseline -t 8 research/phase2/p2_6_scan.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "toy.jl"))
include(joinpath(@__DIR__, "timestep.jl"))

envf(k, d) = parse(Float64, get(ENV, k, string(d)))
const H, RATE, DIFF = envf("H", 0.1), envf("RATE", 0.1), envf("D", 3.0)
const TAUS = log_taus(envf("TAU_LO", 1e-5), envf("TAU_HI", 1e-1))
const NSPINS = parse(Int, get(ENV, "NSPINS", "100000"))
const T = envf("T", 10.0)
const LABEL = get(ENV, "LABEL", "h$(H)_rate$(RATE)_D$(DIFF)")
const OUT = phase2_outdir("p2_6_scan_" * LABEL)

rows = NamedTuple[]
for τ in TAUS
    t0 = time()
    r = toy_rate(; h=H, rate=RATE, D=DIFF, tau=τ, T=T, nspins=NSPINS, seeds=1:10)
    rt = time() - t0
    push!(rows, (tau_ms=τ, R=r.mean, R_sem=r.sem, T_eff=r.T_eff, runtime_s=rt,
        ns_per_spin_step=rt / (10 * NSPINS * round(T / τ)) * 1e9))
    @printf("τ = %.2e: R = %.6e ± %.1e (%.0f s)\n", τ, r.mean, r.sem, rt)
end
iref = argmin(TAUS)
R_ref, R_ref_sem = rows[iref].R, rows[iref].R_sem
Rs = [r.R for r in rows]
conv = Dict{String, Any}()
for tol in (0.01, 0.001)
    ci = tau_conv_interp(TAUS, Rs, R_ref; tol=tol)
    cg = tau_conv_grid(TAUS, Rs, R_ref; tol=tol)
    conv[string(tol)] = Dict("tau_conv" => ci.tau, "tau_conv_grid" => cg.tau, "lower_bound" => ci.lower_bound,
        "noise_limited" => noise_limited(R_ref, R_ref_sem; tol=tol))
    @printf("tol %.1f %%: τ_conv = %.3e ms%s%s\n", 100tol, ci.tau, ci.lower_bound ? " (lower bound)" : "",
        noise_limited(R_ref, R_ref_sem; tol=tol) ? " (noise-limited)" : "")
end
bias_large = (rows[end].R - R_ref) / R_ref
@printf("bias at τ = %.1e: %+.2e relative (%s attenuation than the reference)\n", TAUS[end], bias_large, bias_large < 0 ? "less" : "more")
write_csv(joinpath(OUT, "scan.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "h_um" => H, "rate" => RATE, "D" => DIFF,
    "T_ms" => T, "nspins" => NSPINS, "nseeds" => 10, "taus" => TAUS, "R_ref" => R_ref, "R_ref_sem" => R_ref_sem,
    "tau_conv" => conv, "bias_at_largest_tau" => bias_large,
    "h2_over_D" => H^2 / DIFF, "step_length_at_tau_conv_1pct" => sqrt(2DIFF * conv["0.01"]["tau_conv"])))

using CairoMakie
fig = Figure(size=(560, 380))
ax = Axis(fig[1, 1], xscale=log10, xlabel="τ (ms)", ylabel="R/R_ref − 1", title="toy τ scan: h = $H µm, ΔR₂(0) = $RATE, D = $DIFF")
errorbars!(ax, TAUS, Rs ./ R_ref .- 1, [r.R_sem for r in rows] ./ R_ref); scatterlines!(ax, TAUS, Rs ./ R_ref .- 1)
hlines!(ax, [-0.01, 0.01], color=:gray, linestyle=:dash); hlines!(ax, [-0.001, 0.001], color=:gray, linestyle=:dot)
save(joinpath(OUT, "scan.png"), fig)
println("saved to $OUT")
```

- [ ] **Step 2: Run a pilot to set the spin count and confirm the cost**

Run: `H=0.1 RATE=0.1 D=3 TAU_LO=1e-4 TAU_HI=1e-1 NSPINS=20000 T=10 LABEL=pilot julia --project=research/baseline -t 8 research/phase2/p2_6_scan.jl`
Expected: 10 τ lines and two τ_conv lines. Then read `ns_per_spin_step` and `R_ref_sem/R_ref` from the results:
- **Cost:** a full V1 scan down to 1e-5 with N spins costs about 10 · N · (T/1e-5) · 1.5 · ns_per_spin_step. Choose N so that one V1 scan stays under 1 h, but no lower than 20 000.
- **Precision:** if R_ref_sem/R_ref at the chosen N exceeds 0.003 (tol/3 at 1%), increase N or T until it doesn't, and record the choice as a ruling.
- **If both can't hold** (a full scan to 1e-5 would take far over 1 h at N = 20 000): raise `TAU_LO` for all scans to one decade below the pilot's 1% τ_conv (scaled by h²/D for the V2/V3 scans). The reference only has to be converged, not at 1e-5. Ledger it as a ruling.

Delete `results/phase2/p2_6_scan_pilot/` afterwards. It's only for sizing.

- [ ] **Step 3: Commit the script**

```bash
git add research/phase2/p2_6_scan.jl
git commit -m "Add the P2.6 toy timestep scan script

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: V1 timestep plateau and MCMR confirmation (P2.6.1)

**Files:**
- Create: `research/phase2/p2_6_1_mcmr_check.jl`
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `p2_6_scan.jl` (Task 2); `run_walls`, `layer_walls` (harness).
- Produces: `results/phase2/p2_6_scan_V1_h0.1/` and `results/phase2/p2_6_1_mcmr_check/summary.json`, with per τ the MCMR rate, the toy rate and their z.

- [ ] **Step 1: Run the V1 scan**

Run: `H=0.1 RATE=0.1 D=3 TAU_LO=1e-5 TAU_HI=1e-1 NSPINS=<N from Task 2> T=10 LABEL=V1_h0.1 julia --project=research/baseline -t 8 research/phase2/p2_6_scan.jl`
Expected: 13 τ lines, τ_conv at 1% (and 0.1%, possibly noise-limited), and the large-τ bias with its sign.

- [ ] **Step 2: Create the MCMR confirmation script**

```julia
# P2.6.1 — MCMR at selected τ vs the toy scan (h = 0.1 µm, ΔR2(0) = 0.1, D = 3, T = 10 ms, seeds 1–10).
# Pass: |R_MCMR − R_toy| ≤ 2 combined SEM at every τ in TAUS (independent random streams, so the check is statistical).
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_6_1_mcmr_check.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_6_1_mcmr_check")
const TAUS = [1e-1, 1e-2, 1e-3, 1e-4]           # grid points of the V1 scan that divide T, so toy and MCMR simulate the same T
const T = 10.0
const NSPINS = parse(Int, get(ENV, "NSPINS", "100000"))
# scan.csv columns: tau_ms, R, R_sem, …; keyed by τ rounded to 6 significant digits so grid values match the literals above
scan = Dict{Float64, Tuple{Float64, Float64}}()
for l in readlines(joinpath(REPO_ROOT, "research", "results", "phase2", "p2_6_scan_V1_h0.1", "scan.csv"))[2:end]
    c = parse.(Float64, split(l, ",")[1:3])
    scan[round(c[1], sigdigits=6)] = (c[2], c[3])
end

rows = NamedTuple[]
for τ in TAUS
    r = run_walls(layer_walls(h=0.1, rate=0.1); timestep=τ, times=[0.0, T], nspins=NSPINS)
    Rs = -log.(r.S[:, end]) ./ T                  # r.S is seeds × readouts; the readout is at exactly T
    Rm, Rse = mean(Rs), std(Rs) / sqrt(length(Rs))
    (Rt, Rtse) = scan[round(τ, sigdigits=6)]
    z = (Rm - Rt) / sqrt(Rse^2 + Rtse^2)
    push!(rows, (tau_ms=τ, R_mcmr=Rm, R_mcmr_sem=Rse, R_toy=Rt, R_toy_sem=Rtse, z=z, runtime_s=r.runtime))
    @printf("τ = %.1e: MCMR %.6e ± %.1e, toy %.6e ± %.1e, z = %+.2f\n", τ, Rm, Rse, Rt, Rtse, z)
end
pass = all(abs(r.z) <= 2 for r in rows)
@printf("MCMR vs toy → %s\n", pass ? "PASS" : "FAIL")
write_csv(joinpath(OUT, "mcmr_vs_toy.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rows" => rows, "pass" => pass))
println("saved to $OUT")
```

The values in `TAUS` are grid points of the V1 scan (every third point), and each divides T = 10 ms, so the toy's T_eff equals T.

- [ ] **Step 3: Run it**

Run: `julia --project=research/baseline -t 8 research/phase2/p2_6_1_mcmr_check.jl`
Expected: 4 lines with |z| ≤ 2 and `MCMR vs toy → PASS`. About 40 min, mostly τ = 1e-4. A FAIL means the toy doesn't represent MCMR at that τ. Stop and investigate before using the toy for V2/V3.

- [ ] **Step 4: Record and commit**

Add a summary row and a section "P2.6.1: V1 timestep plateau" to `research/phase2_record.md`. It should contain the scan table (τ, R ± SEM, relative error), τ_conv at 1% and 0.1% (with lower-bound and noise-limited flags), the large-τ bias with its sign compared with the proposal's expectation, the MCMR confirmation table, and the scan figure. The interpretation must be numbers only.

```bash
git add research/phase2/p2_6_1_mcmr_check.jl research/results/phase2/p2_6_scan_V1_h0.1 research/results/phase2/p2_6_1_mcmr_check research/phase2_record.md
git commit -m "Run V1 timestep plateau (toy scan, MCMR confirmation)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: V2 h² scaling, rate dependence, V3 1/D scaling (P2.6.2, P2.6.3)

**Files:**
- Create: `research/phase2/p2_6_2_scaling.jl`
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `p2_6_scan.jl`, `fit_c_h`.
- Produces: scans `p2_6_scan_V2_h0.05`, `p2_6_scan_V2_h0.2`, `p2_6_scan_V2_rate0.5`, `p2_6_scan_V3_D1`, `p2_6_scan_V3_D6`; and `results/phase2/p2_6_2_scaling/summary.json` with `c_h` (1%), its spread, the per-scan ratios τ_conv·D/h², the rate test, and `c_h_recommended` = c_h/10 (one-decade margin, P2.6.5).

- [ ] **Step 1: Run the five scans**

Each one takes the same spin count N as V1. The τ range is chosen so the expected τ_conv lies well inside it; extend `TAU_LO` by one decade if a scan reports `lower_bound` = false at its smallest τ (it fails immediately).

```bash
H=0.05 RATE=0.1 D=3 TAU_LO=1e-6 TAU_HI=1e-1 NSPINS=<N> T=10 LABEL=V2_h0.05 julia --project=research/baseline -t 8 research/phase2/p2_6_scan.jl
H=0.2  RATE=0.1 D=3 TAU_LO=1e-5 TAU_HI=1e-1 NSPINS=<N> T=10 LABEL=V2_h0.2  julia --project=research/baseline -t 8 research/phase2/p2_6_scan.jl
H=0.1  RATE=0.5 D=3 TAU_LO=1e-5 TAU_HI=1e-1 NSPINS=<N> T=10 LABEL=V2_rate0.5 julia --project=research/baseline -t 8 research/phase2/p2_6_scan.jl
H=0.1  RATE=0.1 D=1 TAU_LO=1e-5 TAU_HI=1e-1 NSPINS=<N> T=10 LABEL=V3_D1   julia --project=research/baseline -t 8 research/phase2/p2_6_scan.jl
H=0.1  RATE=0.1 D=6 TAU_LO=1e-6 TAU_HI=1e-1 NSPINS=<N> T=10 LABEL=V3_D6   julia --project=research/baseline -t 8 research/phase2/p2_6_scan.jl
```

Expected: each prints τ_conv at 1% that is not a lower bound and not noise-limited. If the h = 0.05 or D = 6 scan is too slow at τ = 1e-6, reduce its `NSPINS` by up to 4× and record the ruling. Its precision flag will show whether that's still enough.

- [ ] **Step 2: Create the combining script**

```julia
# P2.6.2 V2 (h² scaling), rate dependence (tracker decision point), P2.6.3 V3 (1/D scaling).
# Reads the six toy scans and fits τ_conv = c_h·h²/D at tolerance 1 %.
# V2 pass: τ_conv·D/h² equal within ±30 % across h = 0.05, 0.1, 0.2 (a line through the origin).
# V3 pass: τ_conv·D equal within ±30 % across D = 1, 3, 6 at h = 0.1.
# Rate test (reported): τ_conv at ΔR2(0) = 0.5 vs 0.1 at h = 0.1. If V2 fails, the tracker's decision point applies
# (test τ_conv against 1/ΔR2(0)).
#
# Run: julia --project=research/baseline research/phase2/p2_6_2_scaling.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "timestep.jl"))

const OUT = phase2_outdir("p2_6_2_scaling")
load(label) = JSON.parsefile(joinpath(REPO_ROOT, "research", "results", "phase2", "p2_6_scan_" * label, "summary.json"))
scans = Dict(l => load(l) for l in ("V2_h0.05", "V1_h0.1", "V2_h0.2", "V2_rate0.5", "V3_D1", "V3_D6"))
tc(s) = s["tau_conv"]["0.01"]["tau_conv"]
flags(s) = (s["tau_conv"]["0.01"]["lower_bound"], s["tau_conv"]["0.01"]["noise_limited"])

v2 = ["V2_h0.05", "V1_h0.1", "V2_h0.2"]
f2 = fit_c_h([tc(scans[l]) for l in v2], [scans[l]["h_um"] for l in v2], [scans[l]["D"] for l in v2])
v2_pass = maximum(abs.(f2.ratios ./ f2.c_h .- 1)) <= 0.3
v3 = ["V3_D1", "V1_h0.1", "V3_D6"]
d_ratios = [tc(scans[l]) * scans[l]["D"] for l in v3]
v3_pass = maximum(abs.(d_ratios ./ mean(d_ratios) .- 1)) <= 0.3
rate_ratio = tc(scans["V2_rate0.5"]) / tc(scans["V1_h0.1"])
all_c = fit_c_h([tc(scans[l]) for l in [v2; v3[[1, 3]]]], [scans[l]["h_um"] for l in [v2; v3[[1, 3]]]], [scans[l]["D"] for l in [v2; v3[[1, 3]]]])

for l in sort(collect(keys(scans)))
    s = scans[l]
    @printf("%-11s h=%.3f rate=%.2f D=%.1f: τ_conv(1%%) = %.3e ms, τ_conv·D/h² = %.4f  flags(lower bound, noise) = %s\n",
        l, s["h_um"], s["rate"], s["D"], tc(s), tc(s) * s["D"] / s["h_um"]^2, flags(s))
end
@printf("V2: c_h = %.4f (ratios %s, max deviation %.0f %%) → %s\n", f2.c_h, round.(f2.ratios, digits=4), 100maximum(abs.(f2.ratios ./ f2.c_h .- 1)), v2_pass ? "PASS" : "FAIL")
@printf("V3: τ_conv·D = %s (max deviation %.0f %%) → %s\n", round.(d_ratios, sigdigits=4), 100maximum(abs.(d_ratios ./ mean(d_ratios) .- 1)), v3_pass ? "PASS" : "FAIL")
@printf("rate: τ_conv(ΔR2(0) = 0.5) / τ_conv(0.1) = %.3f\n", rate_ratio)
@printf("c_h over all h and D scans = %.4f (spread %.0f %%); recommended constraint c_h/10 = %.5f\n", all_c.c_h, 100all_c.spread, all_c.c_h / 10)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(),
    "tau_conv_1pct" => Dict(l => tc(s) for (l, s) in scans), "flags" => Dict(l => flags(s) for (l, s) in scans),
    "v2_c_h" => f2.c_h, "v2_ratios" => f2.ratios, "v2_pass" => v2_pass,
    "v3_tau_conv_times_D" => d_ratios, "v3_pass" => v3_pass, "rate_ratio" => rate_ratio,
    "c_h" => all_c.c_h, "c_h_spread" => all_c.spread, "c_h_recommended" => all_c.c_h / 10))

using CairoMakie
fig = Figure(size=(560, 400))
ax = Axis(fig[1, 1], xscale=log10, yscale=log10, xlabel="h²/D (ms)", ylabel="τ_conv at 1 % (ms)", title="V2/V3: τ_conv vs h²/D")
for l in sort(collect(keys(scans)))
    s = scans[l]
    scatter!(ax, [s["h_um"]^2 / s["D"]], [tc(s)], label=l)
end
x = 10 .^ range(-4, -1, length=50)
lines!(ax, x, all_c.c_h .* x, color=:gray, linestyle=:dash, label="c_h·h²/D")
axislegend(ax, position=:lt, labelsize=8)
save(joinpath(OUT, "scaling.png"), fig)
println("saved to $OUT")
```

- [ ] **Step 3: Run it**

Run: `julia --project=research/baseline research/phase2/p2_6_2_scaling.jl`
Expected: six scan lines, V2 and V3 verdicts, the rate ratio, and c_h with c_h/10. If V2 FAILs, stop and follow the tracker's decision point: compare τ_conv with 1/ΔR₂(0) using the rate scan and record the result, before doing Task 5. The constraint's form then has to be decided by you, not by the plan.

- [ ] **Step 4: Record and commit**

Add a section "P2.6.2, P2.6.3: scaling of τ_conv" to `research/phase2_record.md`, with the table of the six scans (h, ΔR₂(0), D, τ_conv at 1% and 0.1%, flags, τ_conv·D/h², step length √(2Dτ_conv)), the V2/V3 verdicts, the rate ratio, c_h with its spread, the recommended c_h/10, and the figure.

```bash
git add research/phase2/p2_6_2_scaling.jl research/results/phase2/p2_6_scan_V2_h0.05 research/results/phase2/p2_6_scan_V2_h0.2 research/results/phase2/p2_6_scan_V2_rate0.5 research/results/phase2/p2_6_scan_V3_D1 research/results/phase2/p2_6_scan_V3_D6 research/results/phase2/p2_6_2_scaling research/phase2_record.md
git commit -m "Run V2/V3 scaling of the converged timestep and fit c_h

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Layer constraint in the simulator's timestep (P2.6.5)

**Files:**
- Modify: `src/geometries/internal/layers.jl`, `src/geometries/internal/internal.jl`, `src/timesteps.jl`
- Modify: `test/test_layer.jl`

**Interfaces:**
- Consumes: `c_h_recommended` from `results/phase2/p2_6_2_scaling/summary.json` (Task 4).
- Produces:
  - `Internal.min_layer_h(geometry)::Float64`: the smallest h over sides with ρ > 0 (Inf without a layer)
  - `TimeStep(; …, layer=DEFAULT_LAYER_SCALING)` with the new option `layer · min_layer_h² / D`
  - `Simulation(...; timestep=(layer=…,))` to override it

- [ ] **Step 1: Append the failing tests to `test/test_layer.jl`**

```julia
@testset "test_layer.jl: layer timestep constraint" begin
    Layers = mr.Geometries.Internal.Layers
    geom(walls) = mr.Simulation([]; geometry=walls, verbose=false).geometry
    @test Layers.min_layer_h(geom(mr.Walls(repeats=2.))) == Inf
    @test Layers.min_layer_h(geom(mr.Walls(repeats=2., layer_h=0.5))) == Inf              # ρ = 0 → no layer
    @test Layers.min_layer_h(geom(mr.Walls(repeats=2., layer_rho_positive=0.01, layer_h_positive=0.1,
        layer_rho_negative=0.05, layer_h_negative=0.5))) == 0.1
    # constraint binds for a thin layer, is absent without a layer
    s = mr.Simulation([]; geometry=mr.Walls(repeats=2., layer_rho=0.01, layer_h=0.1), diffusivity=3., verbose=false, timestep=(layer=0.5,))
    @test s.timestep.max_timestep ≈ 0.5 * 0.1^2 / 3
    s0 = mr.Simulation([]; geometry=mr.Walls(repeats=2.), diffusivity=3., verbose=false, timestep=(layer=0.5,))
    @test s0.timestep.max_timestep ≈ 0.04                                                   # tortuosity, unchanged
    # runtests.jl calls Logging.disable_logging(Logging.Info) globally; lift it for this one check and restore it
    Logging.disable_logging(Logging.Debug)
    try
        @test_logs (:info, r"near-surface layer") match_mode=:any mr.Simulation([]; geometry=mr.Walls(repeats=2., layer_rho=0.01, layer_h=0.1),
            diffusivity=3., verbose=true, timestep=(layer=0.5,))
    finally
        Logging.disable_logging(Logging.Info)
    end
end
```

`FixedGeometry` is an `NTuple` of groups (`fixed_obstruction_groups.jl:104`), so `min_layer_h` iterates it with `init=Inf`, which also covers an empty geometry.

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`
Expected: FAIL with `UndefVarError: min_layer_h not defined`.

- [ ] **Step 3: Add `min_layer_h` to `src/geometries/internal/layers.jl`**

Before the final `end` of the module:

```julia
"""
    min_layer_h(geometry)

Smallest layer thickness h (um) over all sides with ρ > 0 in the geometry; Inf if there is no layer.
"""
min_layer_h(::Nothing) = Inf
min_layer_h(layers::Vector{WallLayer}) = minimum((s.h for l in layers for s in (l.positive, l.negative) if s.rho > 0); init=Inf)
min_layer_h(geometry::FixedGeometry) = minimum((min_layer_h(g.layer) for g in geometry); init=Inf)
```

In `src/geometries/internal/internal.jl`, extend the import: `import .Layers: LayerSide, WallLayer, layer_exponent, has_layer, min_layer_h`.

- [ ] **Step 4: Add the option to `src/timesteps.jl`**

At the top of the module, after the imports, add (the value comes from `c_h_recommended` in `results/phase2/p2_6_2_scaling/summary.json`, rounded down to 2 significant digits):

```julia
"Default `layer` scaling c_h in τ ≤ c_h·h²/D (P2.6.5): the measured c_h for a 1 % rate error divided by 10."
const DEFAULT_LAYER_SCALING = <c_h_recommended, 2 significant digits>
```

Replace the keyword list and the options tuple:

```julia
function TimeStep(; 
    diffusivity, geometry, size_scale=nothing, 
    tortuosity=3e-2, gradient=1e-4, verbose=true,
    permeability=0.5, surface_relaxation=0.01,
    transfer=0.01, dwell_time=0.1, layer=DEFAULT_LAYER_SCALING,
    )
```

```julia
    options = (
            tortuosity * use_size_scale^2 / diffusivity,
            Internal.max_timestep_sticking(geometry, diffusivity, transfer),
            max_timestep_permeability(geometry, permeability),
            max_timestep_surface_relaxation(geometry, surface_relaxation),
            Internal.min_dwell_time(geometry) * dwell_time,
            layer * Internal.min_layer_h(geometry)^2 / diffusivity,
        )
```

After the `elseif idx == 5 … end` block's last branch, add:

```julia
        elseif idx == 6
            push!(lines, "Maximum timestep set by the near-surface layer (τ ≤ c_h·h²/D with the thinnest layer h = $(Internal.min_layer_h(geometry)) um) to $(options[6]) ms.")
            push!(lines, "You can alter it by changing `timestep=(layer=...)` from its current value of $(layer).")
```

Also add a line 7 to the docstring list: "7. `layer` · h_min² / D, where h_min is the thinnest near-surface layer (P2.6.5)."

- [ ] **Step 5: Run the tests and check they pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer", "collisions", "permeability", "transfer"])'`
Expected: PASS. Without a layer every timestep is unchanged, because `min_layer_h` = Inf gives an Inf option.

- [ ] **Step 6: Commit**

```bash
git add src/geometries/internal/layers.jl src/geometries/internal/internal.jl src/timesteps.jl test/test_layer.jl
git commit -m "Add the near-surface layer timestep constraint tau <= c_h h^2/D

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Cost and consequences for earlier results (P2.6.8, P2.6.10–P2.6.12)

**Files:**
- Create: `research/phase2/p2_6_8_cost.jl`
- Modify: `research/phase2_record.md`, `research/project-progress-tracker.md`

**Interfaces:**
- Consumes: `run_walls`, `layer_walls`, B5 (`results/baseline/b5_runtime_vs_tau/runtime.csv`), the V1 scan, c_h from Task 4, V4 h\* values (`p2_3_9_fit_h/summary.json`).
- Produces: `results/phase2/p2_6_8_cost/{runtime.csv, consequences.csv, summary.json, cost.png}`.

- [ ] **Step 1: Create the script**

```julia
# P2.6.8 runtime vs τ with the layer (vs B5, old code); P2.6.10 runtime vs h at the converged τ;
# P2.6.11 noise cost at fixed budget; P2.6.12 accuracy–cost summary; and the consequence for the earlier
# τ = 1e-2 ms results: τ_conv(h) = c_h·h²/D for every h used in the sweep and in V4.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_6_8_cost.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_6_8_cost")
R = joinpath(REPO_ROOT, "research", "results")
c_h = JSON.parsefile(joinpath(R, "phase2", "p2_6_2_scaling", "summary.json"))["c_h"]
v1 = JSON.parsefile(joinpath(R, "phase2", "p2_6_scan_V1_h0.1", "summary.json"))
b5 = Dict(parse(Float64, split(l, ",")[2]) => parse(Float64, split(l, ",")[9])
    for l in readlines(joinpath(R, "baseline", "b5_runtime_vs_tau", "runtime.csv"))[2:end] if split(l, ",")[1] == "0")

# runtime per spin-step with the layer (h = 0.1), same workload as B5: 10 000 spins × 20 000 steps
rows = NamedTuple[]
for τ in (1e-1, 1e-2, 1e-3, 1e-4)
    T = 20_000 * τ
    run_walls(layer_walls(h=0.1, rate=0.1); timestep=τ, times=[0.0, 20τ], seeds=1:1, nspins=1000)   # warm-up
    r = run_walls(layer_walls(h=0.1, rate=0.1); timestep=τ, times=[0.0, T], seeds=1:1, nspins=10_000)
    us = r.runtime / (10_000 * 20_000) * 1e6
    push!(rows, (kind="tau", tau_ms=τ, h_um=0.1, us_per_spin_step=us, b5_us_per_spin_step=get(b5, τ, NaN), overhead=us / get(b5, τ, NaN)))
    @printf("τ = %.0e: %.4f µs/spin-step with layer, %.4f without (B5) → ×%.2f\n", τ, us, get(b5, τ, NaN), us / get(b5, τ, NaN))
end
# runtime vs h at the timestep the new default constraint gives (c_h/10, capped at 1e-2), 1e4 spins, T = 10 ms, 1 seed
c_rec = JSON.parsefile(joinpath(R, "phase2", "p2_6_2_scaling", "summary.json"))["c_h_recommended"]
for h in (0.01, 0.05, 0.1, 0.2, 0.5)
    τc = min(c_rec * h^2 / D_REF, 1e-2)
    r = run_walls(layer_walls(h=h, rate=0.1); timestep=τc, times=[0.0, 10.0], seeds=1:1, nspins=10_000)
    push!(rows, (kind="h", tau_ms=τc, h_um=h, us_per_spin_step=r.runtime / (10_000 * round(10 / τc)) * 1e6,
        b5_us_per_spin_step=NaN, overhead=NaN))
    @printf("h = %.2f: τ = %.2e ms, 1e5 spins × 50 ms × 10 seeds ≈ %.1f h\n", h, τc, rows[end].us_per_spin_step * 1e5 * (50 / τc) * 10 / 3.6e9)
end
# consequence for earlier results: which h used at τ = 1e-2 were converged (τ_conv ≥ 1e-2)?
fits = JSON.parsefile(joinpath(R, "phase2", "p2_3_9_fit_h", "summary.json"))["fits"]
cons = NamedTuple[]
for h in sort(unique(vcat(H_GRID[2:end], [f["h_final"] for f in fits])))
    τc = c_h * h^2 / D_REF
    push!(cons, (h_um=h, tau_conv_1pct=τc, converged_at_1e_2=τc >= 1e-2))
end
write_csv(joinpath(OUT, "runtime.csv"), rows)
write_csv(joinpath(OUT, "consequences.csv"), cons)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "c_h" => c_h, "rows" => rows,
    "consequences" => cons, "v1_tau_conv_1pct" => v1["tau_conv"]["0.01"]["tau_conv"]))

using CairoMakie
fig = Figure(size=(560, 380))
ax = Axis(fig[1, 1], xscale=log10, yscale=log10, xlabel="runtime per spin per ms of sequence (µs)", ylabel="|relative error of R|",
    title="Accuracy–cost, h = 0.1 µm (V1 scan)")
scan = [split(l, ",") for l in readlines(joinpath(R, "phase2", "p2_6_scan_V1_h0.1", "scan.csv"))[2:end]]
taus = [parse(Float64, s[1]) for s in scan]; Rv = [parse(Float64, s[2]) for s in scan]
err = abs.(Rv ./ Rv[argmin(taus)] .- 1)
us_step = mean(r.us_per_spin_step for r in rows if r.kind == "tau")
keep = err .> 0
scatter!(ax, (us_step ./ taus)[keep], err[keep])
for (x, y, t) in zip((us_step ./ taus)[keep], err[keep], taus[keep])
    text!(ax, x, y, text=@sprintf("%.0e", t), fontsize=8)
end
save(joinpath(OUT, "cost.png"), fig)
println("saved to $OUT")
```

The B5 CSV columns are `theta, tau_ms, sim_time_ms, nspins, nsteps, runtime_median_s, runtime_min_s, runtime_max_s, us_per_spin_step, …` (checked when the plan was written). The parser takes column 2 (τ) and column 9 (µs per spin-step) for θ = 0. B5 has exactly τ = 0.1, 0.01, 0.001, 1e-4, 1e-5, so the `get(b5, τ, NaN)` lookups hit.

- [ ] **Step 2: Run it**

Run: `julia --project=research/baseline -t 8 research/phase2/p2_6_8_cost.jl`
Expected:
- **Overhead:** four lines with the layer's cost relative to B5. An overhead well above ×2 is worth a note in the record.
- **Runtime vs h:** five lines, with run-time estimates for a full reference run at the converged τ.
- **Outputs:** `consequences.csv` and the figure.

About 30 min.

- [ ] **Step 3: Record, tick the tracker, commit**

Add a section "P2.6.5, P2.6.8, P2.6.10–P2.6.12: constraint, cost and consequences" to `research/phase2_record.md`. It should contain:
- c_h and the adopted default (c_h/10), and which constraint binds in the reference configuration for h = 0.01, 0.1 and 0.5 µm
- the runtime table (with layer vs B5)
- the runtime-vs-h table with full-run estimates
- the noise-cost statement (at a fixed budget, cutting τ by a factor k means √k fewer spins' worth of precision)
- the accuracy–cost figure
- **the consequences table:** which h in the sweep and the V4 fits were converged at τ = 1e-2 ms (τ_conv ≥ 1e-2) and which were not

For the unconverged ones, state which results (P2.3.4(b), P2.3.9, and the linear plan if it was run) need a rerun at τ ≤ τ_conv. That rerun is a follow-up decision for you, not part of this plan.

In `research/project-progress-tracker.md`:
- **P2.6.1, P2.6.2, P2.6.3, P2.6.5, P2.6.8, P2.6.10, P2.6.11, P2.6.12:** ☑, with outputs.
- **P2.6.4:** "after the linear plan (tools ready: `p2_6_scan.jl` with a shape option to add)".
- **P2.6.6, P2.6.7, P2.6.9:** "not applicable yet: every implemented profile has a closed form; revisit when a profile or geometry needs numerical integration".
- **P2.6.13:** stays ☐ (needs 2.4 and 2.5).
- Increase the Phase 2 "Done" count by 8.

```bash
git add research/phase2/p2_6_8_cost.jl research/results/phase2/p2_6_8_cost research/phase2_record.md research/project-progress-tracker.md
git commit -m "Measure layer cost, record timestep constraint and consequences

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
