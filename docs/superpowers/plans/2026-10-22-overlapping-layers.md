# Overlapping Near-Surface Layers (Linear, Exponential) vs the Conventional Model Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Characterise overlapping near-surface layers (where the layers of the two faces of a gap overlap), for the linear and exponential profiles only, and compare them with the conventional surface-relaxation model at the same ρ. First verify that the summed rate profile is realistic: the rate must never be higher mid-gap than near the walls. Then measure how the signal, the apparent T₂ and the D dependence change from small to large overlap.

**Architecture:**
- **Profile check (Task 1).** Evaluate the summed rate ΔR₂(x) across the gap in two independent ways: the package's implementation (`Layers.layer_exponent` for a stationary spin) and closed forms. Test that ΔR₂ does not increase from the wall to the middle, and that the two agree. Report overlap measures: the width of the overlap region, the mid/wall ratio, and the share of each face's integral that lies beyond the gap middle. The step profile with h > w/2 is shown only as the counter-example that motivated the exclusion.
- **Signal comparison (Task 2).** Uses the exact continuum benchmark (`research/phase2/continuum.jl`; no MC noise or timestep bias) over an overlap grid and several ρw/D:
  - ΔS(t) at fixed ρ;
  - ΔT₂/T₂ at TE = 50;
  - the best-fit conventional residual (single curve);
  - the D dependence of the decay (sensitivity of the apparent T₂ to D) for the layer and for the conventional model.
- **Monte Carlo spot check (Task 3).** Confirms the continuum on a few overlapping configurations. Thick layers have a low ΔR₂(0), so τ = 1e-2 suffices (C1: h ≥ 0.5 agrees within 1σ). The check uses τ/4 anyway.

**Tech Stack:** Julia 1.12, MCMRSimulator (this repo), `research/phase2/` (continuum.jl, harness.jl, sensitivity.jl), CairoMakie, `Test` stdlib.

**Spec:**
- The user (2026-10-22): "experiment about overlapping layers, only for linear and exponential profile, compare with conventional model. before the test, … check whether the overlapping area has unrealistic shape, like the value in middle is higher than near wall … overlap small area, large area … one of the character that different from conventional model".
- Tracker §2.5 scope decision: the step profile is excluded for h > w/2. Overlapping linear and exponential layers are included because their summed rate is never higher mid-gap than at the walls; that was derived and is now to be tested.
- Joint fit: left to future work (user, 2026-10-22).
- The user's framing: the conventional model is a reference, not the truth; no fast-diffusion assumptions.

## Global Constraints

- Work only on branch `cc/near-surface-r2`.
- **Overlap grid (w = 2 µm):**
  - **linear:** h/w ∈ {0.25, 0.5, 0.6, 0.75, 0.9, 1.0}. 0.25 means no overlap, 0.5 means the two layers touch, 1.0 means full overlap (uniform).
  - **exponential (cutoff = w):** h/w ∈ {0.05, 0.1, 0.25, 0.5, 0.75, 1.0}. Its tail always reaches the opposite wall, so the overlap is graded by the mid/wall ratio rather than switched on.
- **Strength:** ρw/D ∈ {0.0067 (the P2.5 reference, ρ = 0.01), 0.1, 0.3, 1, 3} at w = 2, D = 3.
- **D dependence:** D ∈ {0.5, 1, 2, 3} at ρ = 0.1 (ρw/D = 0.067 to 0.4) and at ρ = 1 (0.67 to 4).
- **Comparison:** the conventional model at the same ρ (θ = ρ√(π/D), P2.5.0). R₂ᵇ = 1/2000 for the T₂ quantities, and 0 for ΔS.
- **Observability:** as before. |ΔS| is compared with 2/SNR, |ΔT₂/T₂| with the two-point 2σ_T₂/T₂ at SNR 100 and 200, and the best-fit residual is reported with its caveats (single condition; maximum over 41 readouts). All are characterisation, not pass/fail. Only the Task 1 profile checks and the Task 3 spot check have pass/fail.
- **Solver:** N = 3000 (N = 1500 inside best-fit searches, within 4 % at the extremes).
- **Final review:** done by a fresh subagent reviewer (user's standing rule).
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `research/phase2/overlap.jl` | Create | `summed_rate_package`, `summed_rate_closed`, `overlap_metrics` |
| `research/phase2/test/test_overlap.jl`, `runtests.jl` | Create / Modify | profile realism and agreement tests |
| `research/phase2/p2_5_overlap_profiles.jl` | Create | profile figure and overlap-metric table |
| `research/phase2/p2_5_overlap.jl` | Create | continuum comparison grid and D dependence |
| `research/phase2/p2_5_overlap_mc_check.jl` | Create | Monte Carlo spot check |
| `research/phase2/p2_5_overlap_figs.jl` | Create | figures |
| `research/phase2_record.md`, `research/project-progress-tracker.md` | Modify | record, tracker row P2.5.9 |

## Review Focus

1. **The check must test what the simulator actually uses**, not only the formula. Test: Task 1 compares the package's `layer_exponent` (stationary spin, dt = 1) with the closed form at every x.
2. **Monotonicity from wall to middle must hold for every grid point**, including the edges h/w = 0.5 and 1.0. At h/w = 1 the profile is flat, which counts as non-increasing with tolerance. Test: Task 1.
3. **The counter-example must fail the test.** Step with h > w/2 must be flagged as unrealistic (mid > wall). This shows the check can detect the problem. Test: Task 1.
4. **The timestep for the Monte Carlo check must use the largest summed rate.** For the exponential near full overlap the wall rate is up to (1 + e^(−w/h)) times ΔR₂(0) (review M1), so τ must be set from the summed maximum. Test: Task 3 computes τ from `maximum(summed_rate_closed)`.
5. **The D-dependence result must compare like with like.** Same ρ for both models at each D, θ recomputed per D. Test: Task 2 uses `theta_for(ρ; D)` (continuum: the surface ρ is D-independent).

---

### Task 1: Is the overlapping profile realistic? (profile check)

**Files:**
- Create: `research/phase2/overlap.jl`, `research/phase2/test/test_overlap.jl`
- Modify: `research/phase2/test/runtests.jl`

**Interfaces:**
- Produces:
  - `summed_rate_package(shape, ρ, h, w, x)`: ΔR₂ at position x from the package (stationary spin, both faces, all images)
  - `summed_rate_closed(shape, ρ, h, w, x)`: closed form
  - `overlap_metrics(shape, h, w)` → (overlap_width, mid_over_wall, opposite_share), where opposite_share is the fraction of one face's ∫ΔR₂ that lies beyond the gap middle
  - `is_realistic(shape, ρ, h, w; n=401)`: ΔR₂ non-increasing from wall to middle (relative tolerance 1e-12)

- [ ] **Step 1: Write the failing tests**

`research/phase2/test/test_overlap.jl`:

```julia
@testset "overlapping layers: profile realism" begin
    w, ρ = 2.0, 0.1
    xs = range(0.0, w, length=201)
    @testset "package agrees with closed form" begin
        for (shape, hs) in (("linear", (0.5, 0.6, 0.75, 0.9, 1.0) .* w), ("exponential", (0.05, 0.25, 0.5, 1.0) .* w))
            for h in hs, x in xs
                @test summed_rate_package(shape, ρ, h, w, x) ≈ summed_rate_closed(shape, ρ, h, w, x) rtol=1e-10 atol=1e-14
            end
        end
    end
    @testset "linear and exponential: never higher mid-gap than at the wall" begin
        for h in (0.25, 0.5, 0.6, 0.75, 0.9, 1.0) .* w
            @test is_realistic("linear", ρ, h, w)
        end
        for h in (0.05, 0.1, 0.25, 0.5, 0.75, 1.0) .* w
            @test is_realistic("exponential", ρ, h, w)
        end
    end
    @testset "closed forms of the overlap plateau and the uniform limit" begin
        # linear, w/2 < h ≤ w: flat plateau (ρ/h)(4 − 2w/h) on the overlap region [w − h, h]; at h = w uniform 2ρ/w
        h = 0.75w
        @test summed_rate_closed("linear", ρ, h, w, w / 2) ≈ ρ / h * (4 - 2w / h)
        @test summed_rate_closed("linear", ρ, h, w, w - h + 0.01) ≈ ρ / h * (4 - 2w / h)
        @test all(summed_rate_closed("linear", ρ, w, w, x) ≈ 2ρ / w for x in xs)
        # exponential: mid/wall = 2e^(−w/2h)/(1 + e^(−w/h))
        h = 0.5w
        m = overlap_metrics("exponential", h, w)
        @test m.mid_over_wall ≈ 2exp(-w / (2h)) / (1 + exp(-w / h)) rtol=1e-10
        @test overlap_metrics("linear", 0.75w, w).overlap_width ≈ 2 * 0.75w - w
        @test overlap_metrics("linear", 0.25w, w).overlap_width == 0.0
        @test overlap_metrics("linear", 0.25w, w).mid_over_wall == 0.0
    end
    @testset "counter-example: step with h > w/2 is unrealistic" begin
        @test !is_realistic("step", ρ, 0.75w, w)
        @test is_realistic("step", ρ, 0.5w, w)
    end
end
```

Add to `research/phase2/test/runtests.jl`, after the `jointfit.jl` (or `continuum.jl`) include line:

```julia
isfile(joinpath(@__DIR__, "..", "overlap.jl")) && include(joinpath(@__DIR__, "..", "overlap.jl"))
```

- [ ] **Step 2: Run the research tests and check they fail**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: FAIL with `UndefVarError: summed_rate_package`.

- [ ] **Step 3: Implement `research/phase2/overlap.jl`**

```julia
# Overlapping near-surface layers (P2.5.9): the summed rate ΔR2(x) across the gap from both faces, from the package and in
# closed form; realism check (never higher mid-gap than at the wall); overlap measures.

using StaticArrays: SVector
const _LY = mr.Geometries.Internal.Layers

"ΔR2(x) as the simulator applies it: layer exponent of a stationary spin at x for dt = 1 (both faces, all wall images)."
function summed_rate_package(shape, ρ, h, w, x)
    g = Simulation([]; geometry=Walls(repeats=w, layer_rho=ρ, layer_h=h, layer_shape=shape), verbose=false).geometry
    p = SVector{3, Float64}(x, 0.3, -1.2)
    return _LY.layer_exponent(g, p, p, 1.0)
end

"Closed-form ΔR2(x) = (ρ/h)[g(x/h) + g((w − x)/h)]; exponential truncated at w and renormalised."
function summed_rate_closed(shape, ρ, h, w, x)
    c = w / h
    g(u) = shape == "step" ? (0 <= u <= 1 ? 1.0 : 0.0) :
           shape == "linear" ? (0 <= u <= 1 ? 2 * (1 - u) : 0.0) :
           shape == "exponential" ? (0 <= u <= c ? exp(-u) / (1 - exp(-c)) : 0.0) : error("Unknown shape $shape")
    return ρ / h * (g(x / h) + g((w - x) / h))
end

"""
    overlap_metrics(shape, h, w)

overlap_width: length of the region where both faces contribute (linear/step: max(0, 2h − w); exponential: w).
mid_over_wall: ΔR2(w/2) / ΔR2(0). opposite_share: fraction of one face's integral ∫(ρ/h)g(d/h)dd lying beyond the
middle of the gap (d > w/2).
"""
function overlap_metrics(shape, h, w)
    width = shape == "exponential" ? w : max(0.0, 2h - w)
    r(x) = summed_rate_closed(shape, 1.0, h, w, x)
    if shape == "exponential"                                    # one face, truncated at w: share beyond w/2 in closed form
        c = w / h
        share = (exp(-c / 2) - exp(-c)) / (1 - exp(-c))
    else
        share = shape == "linear" ? (h <= w / 2 ? 0.0 : (1 - w / (2h))^2) : max(0.0, 1 - w / (2h))
    end
    return (overlap_width=width, mid_over_wall=r(w / 2) / r(0.0), opposite_share=share)
end

"ΔR2 non-increasing from the wall to the middle of the gap (n points), relative tolerance 1e-12."
function is_realistic(shape, ρ, h, w; n=401)
    v = [summed_rate_closed(shape, ρ, h, w, x) for x in range(0.0, w / 2, length=n)]
    return all(v[k + 1] <= v[k] * (1 + 1e-12) + 1e-15 for k in 1:n - 1)
end
```

- [ ] **Step 4: Run the research tests and check they pass**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: PASS. If package vs closed form fails at a support edge (x = w − h, or exactly at the cutoff), check whether it is a measure-zero boundary value (both sides of a jump). In that case move the sample by 1e-9 rather than loosening the tolerance, and record a ruling.

- [ ] **Step 5: Profile figure and metric table**

Create `research/phase2/p2_5_overlap_profiles.jl`:

```julia
# P2.5.9 — summed rate profiles across the gap for overlapping layers (linear, exponential; step h > w/2 as the
# counter-example), from the package (points) and closed form (lines); overlap metrics table.
#
# Run: julia --project=research/baseline research/phase2/p2_5_overlap_profiles.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "overlap.jl"))
using CairoMakie

const OUT = phase2_outdir("p2_5_overlap")
const W, RHO = 2.0, 0.1
grids = Dict("linear" => [0.25, 0.5, 0.6, 0.75, 0.9, 1.0], "exponential" => [0.05, 0.1, 0.25, 0.5, 0.75, 1.0], "step" => [0.5, 0.75])
rows = NamedTuple[]
fig = Figure(size=(1300, 380))
for (j, shape) in enumerate(("linear", "exponential", "step"))
    ax = Axis(fig[1, j], xlabel="x (µm), walls at 0 and $W", ylabel="ΔR₂(x) / (ρ/h)", title=shape == "step" ? "step (counter-example, h > w/2)" : shape)
    for hw in grids[shape]
        h = hw * W
        xs = range(0.0, W, length=401)
        lines!(ax, xs, [summed_rate_closed(shape, RHO, h, W, x) / (RHO / h) for x in xs], label="h/w = $hw")
        xp = range(0.0, W, length=21)
        scatter!(ax, xp, [summed_rate_package(shape, RHO, h, W, x) / (RHO / h) for x in xp], markersize=4)
        m = overlap_metrics(shape, h, W)
        push!(rows, (shape=shape, h_over_w=hw, overlap_width_um=m.overlap_width, mid_over_wall=m.mid_over_wall,
            opposite_share=m.opposite_share, realistic=is_realistic(shape, RHO, h, W)))
        @printf("%-11s h/w=%.2f: overlap width %.2f µm, mid/wall %.3e, opposite share %.3e, realistic %s\n", shape, hw, m.overlap_width, m.mid_over_wall, m.opposite_share, rows[end].realistic)
    end
    axislegend(ax, position=:ct, labelsize=8)
end
save(joinpath(OUT, "p2_5_overlap_profiles.png"), fig)
write_csv(joinpath(OUT, "overlap_metrics.csv"), rows)
println("saved to $OUT")
```

Run: `julia --project=research/baseline research/phase2/p2_5_overlap_profiles.jl`
Expected:
- 14 lines;
- every linear and exponential row has `realistic true`;
- step h/w = 0.75 has `realistic false` (mid/wall = 2);
- in the figure the package points lie on the lines.

- [ ] **Step 6: Commit**

```bash
git add research/phase2/overlap.jl research/phase2/test/test_overlap.jl research/phase2/test/runtests.jl research/phase2/p2_5_overlap_profiles.jl research/results/phase2/p2_5_overlap
git commit -m "P2.5.9: overlapping-layer profiles — realism check (package vs closed form) and overlap metrics

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Overlapping layers vs the conventional model (continuum)

**Files:**
- Create: `research/phase2/p2_5_overlap.jl`

**Interfaces:**
- Consumes: `continuum_signal`, `continuum_rates`, `overlap_metrics`, `apparent_t2`, `t2_noise_rel`, `centroid`.
- Produces: `results/phase2/p2_5_overlap/{grid.csv, ddep.csv, summary.json}`.

- [ ] **Step 1: Create the script**

```julia
# P2.5.9 — overlapping layers (linear, exponential) vs the conventional model at the same ρ, exact continuum benchmark.
# (1) Overlap grid × ρw/D: max |ΔS| over the window (0 until the conventional signal ≈ 1 %, 41 readouts, R2_bulk = 0),
#     ΔT2/T2 at TE = 50 (R2_bulk = 1/2000) with the two-point noise at SNR 100/200, best-fit conventional residual
#     (ρ_c free, scan e^−4…e^8 × ρ + golden; single condition), and the layer's mono-exponential residual.
# (2) D dependence: apparent T2 (TE = 50) vs D ∈ {0.5, 1, 2, 3} at ρ = 0.1 and 1 for the layer at each overlap and for
#     the conventional model; relative change T2(D = 0.5)/T2(D = 3) − 1 per model.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_5_overlap.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "sensitivity.jl"))
include(joinpath(@__DIR__, "continuum.jl"))
include(joinpath(@__DIR__, "overlap.jl"))

const OUT = phase2_outdir("p2_5_overlap")
const W, D0 = 2.0, 3.0
const HW = Dict("linear" => [0.25, 0.5, 0.6, 0.75, 0.9, 1.0], "exponential" => [0.05, 0.1, 0.25, 0.5, 0.75, 1.0])
const RWD = [0.0067, 0.1, 0.3, 1.0, 3.0]

function best_fit(Sl, w, D, ρ0, times)
    f(lr) = maximum(abs.(Sl .- continuum_signal(w=w, D=D, times=times, surface=exp(lr), N=1500)))
    grid = log(ρ0) .+ (-4.0:0.25:8.0)
    v = f.(grid); k = argmin(v)
    a, b = grid[max(k - 1, 1)], grid[min(k + 1, length(grid))]
    gr = (sqrt(5) - 1) / 2
    x1, x2 = b - gr * (b - a), a + gr * (b - a)
    f1, f2 = f(x1), f(x2)
    for _ in 1:30
        if f1 < f2
            b, x2, f2 = x2, x1, f1
            x1 = b - gr * (b - a); f1 = f(x1)
        else
            a, x1, f1 = x1, x2, f2
            x2 = a + gr * (b - a); f2 = f(x2)
        end
    end
    return (resid=min(f1, f2, v[k]), at_edge=k == 1 || k == length(grid))
end
monoexp(t, S) = (k = findall(>(0), t); R = -sum(t[k] .* log.(S[k])) / sum(t[k] .^ 2); maximum(abs.(S .- exp.(-R .* t))))

rows = NamedTuple[]
for shape in ("linear", "exponential"), hw in HW[shape], rwd in RWD
    ρ, h = rwd * D0 / W, hw * W
    tmax = log(100) / continuum_rates(w=W, D=D0, surface=ρ, N=1500)[1]
    t = collect(range(0.0, tmax, length=41))
    Sl = continuum_signal(w=W, D=D0, times=t, layer=(shape, ρ, h))
    Sc = continuum_signal(w=W, D=D0, times=t, surface=ρ)
    bf = best_fit(continuum_signal(w=W, D=D0, times=t, layer=(shape, ρ, h), N=1500), W, D0, ρ, t)
    St = [continuum_signal(w=W, D=D0, times=[0.0, 50.0], layer=(shape, ρ, h), R2_bulk=1 / 2000), continuum_signal(w=W, D=D0, times=[0.0, 50.0], surface=ρ, R2_bulk=1 / 2000)]
    T2l, T2c = apparent_t2(St[1][1], St[1][2], 50.0), apparent_t2(St[2][1], St[2][2], 50.0)
    r = St[2][2]
    m = overlap_metrics(shape, h, W)
    push!(rows, (shape=shape, h_over_w=hw, rho_w_over_D=rwd, rho=ρ, overlap_width_um=m.overlap_width, mid_over_wall=m.mid_over_wall,
        centroid_um=centroid(shape, h; cutoff=W), max_abs_dS=maximum(abs.(Sl .- Sc)), dS_sign=sign(Sl[argmax(abs.(Sl .- Sc))] - Sc[argmax(abs.(Sl .- Sc))]),
        rel_T2_TE50=(T2l - T2c) / T2c, two_sigma_T2_snr200=2t2_noise_rel(r, 200), two_sigma_T2_snr100=2t2_noise_rel(r, 100),
        best_fit_residual=bf.resid, best_fit_at_edge=bf.at_edge, layer_monoexp_residual=monoexp(t, Sl), conv_monoexp_residual=monoexp(t, Sc)))
    q = rows[end]
    @printf("%-11s h/w=%.2f ρw/D=%-6g: mid/wall %.2e, max|ΔS| %.2e, ΔT2/T2 %+.2e (2σ200 %.3f), best-fit %.2e%s, mono-exp layer %.1e conv %.1e\n",
        shape, hw, rwd, q.mid_over_wall, q.max_abs_dS, q.rel_T2_TE50, q.two_sigma_T2_snr200, q.best_fit_residual, q.best_fit_at_edge ? " (edge)" : "", q.layer_monoexp_residual, q.conv_monoexp_residual)
end
write_csv(joinpath(OUT, "grid.csv"), rows)

ddep = NamedTuple[]
for ρ in (0.1, 1.0)
    T2(model, D) = (S = model == "conventional" ? continuum_signal(w=W, D=D, times=[0.0, 50.0], surface=ρ, R2_bulk=1 / 2000) :
                        continuum_signal(w=W, D=D, times=[0.0, 50.0], layer=(model[1], ρ, model[2] * W), R2_bulk=1 / 2000); apparent_t2(S[1], S[2], 50.0))
    models = vcat(["conventional"], [(s, hw) for s in ("linear", "exponential") for hw in HW[s]])
    for mdl in models
        Ts = [T2(mdl, D) for D in (0.5, 1.0, 2.0, 3.0)]
        push!(ddep, (rho=ρ, model=mdl == "conventional" ? "conventional" : mdl[1], h_over_w=mdl == "conventional" ? NaN : mdl[2],
            T2_D05=Ts[1], T2_D1=Ts[2], T2_D2=Ts[3], T2_D3=Ts[4], rel_change_D05_vs_D3=Ts[1] / Ts[4] - 1))
        @printf("D-dependence ρ=%.1f %-12s h/w=%-5s: T2(D=0.5..3) = %s ms; T2(0.5)/T2(3) − 1 = %+.2e\n", ρ, ddep[end].model,
            mdl == "conventional" ? "-" : string(mdl[2]), join(round.(Ts, digits=3), ", "), ddep[end].rel_change_D05_vs_D3)
    end
end
write_csv(joinpath(OUT, "ddep.csv"), ddep)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "w" => W, "D" => D0, "rows" => rows, "ddep" => ddep))
println("saved to $OUT")
```

- [ ] **Step 2: Run it (detached, about 20–30 min)**

Run it detached and wait for `saved to`.
Expected:
- 60 grid lines and 26 D-dependence lines;
- linear h/w = 1.0 has a layer mono-exponential residual ≈ 0, since full overlap is a uniform rate;
- the D-dependence of linear h/w = 1.0 is 0 up to bulk-free numerics, because 2ρ/w doesn't depend on D. The conventional model must show a non-zero D dependence that grows with ρ.

If linear h/w = 1.0 is not D-independent to better than 1e-6, there is a bug: stop and debug.

- [ ] **Step 3: Commit**

```bash
git add research/phase2/p2_5_overlap.jl research/results/phase2/p2_5_overlap
git commit -m "P2.5.9: overlapping layers vs the conventional model (continuum): signal, T2, best-fit, D dependence

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Monte Carlo spot check of overlapping layers

**Files:**
- Create: `research/phase2/p2_5_overlap_mc_check.jl`

**Interfaces:**
- Consumes: `run_walls`, `rho_walls`, `theta_for`, `continuum_signal`, `summed_rate_closed`.
- Produces: `results/phase2/p2_5_overlap_mc_check/summary.json` with `pass`.

- [ ] **Step 1: Create the script**

```julia
# P2.5.9 spot check — Monte Carlo vs continuum for overlapping layers: linear h/w = 0.75 and 1.0, exponential h/w = 0.5
# and 1.0, at ρ = 0.1 and 1.0 (w = 2, D = 3, R2_bulk = 0). τ = (0.005 / max_x ΔR2_summed(x)) / 4, capped at 1e-2, so
# the summed (overlapping) maximum sets the step (review M1). Readouts: 10 points up to the time where the conventional
# continuum signal is 0.05. Seeds 1–10, 50000 spins. Pass: |z| ≤ 3 at every readout (z = (MC − continuum)/SEM).
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_5_overlap_mc_check.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum.jl"))
include(joinpath(@__DIR__, "overlap.jl"))

const OUT = phase2_outdir("p2_5_overlap_mc_check")
const W, D = 2.0, 3.0
rows = NamedTuple[]
for ρ in (0.1, 1.0), (shape, hw) in (("linear", 0.75), ("linear", 1.0), ("exponential", 0.5), ("exponential", 1.0))
    h = hw * W
    rmax = maximum(summed_rate_closed(shape, ρ, h, W, x) for x in range(0.0, W, length=2001))
    τ = min(0.005 / rmax / 4, 1e-2)
    tend = log(20) / continuum_rates(w=W, D=D, surface=ρ, N=1500)[1]
    T = collect(range(0.0, tend, length=10))
    cont = continuum_signal(w=W, D=D, times=T, layer=(shape, ρ, h))
    r = run_walls(rho_walls(rho=ρ, h=h, shape=shape, w=W); D=D, timestep=τ, times=T, nspins=50_000)
    z = (r.mean[2:end] .- cont[2:end]) ./ max.(r.sem[2:end], 1e-9)   # floor: the uniform case has SEM ≈ 1e-16 (exact)
    push!(rows, (rho=ρ, shape=shape, h_over_w=hw, tau_ms=τ, max_abs_z=maximum(abs.(z)), max_abs_diff=maximum(abs.(r.mean .- cont)), pass=maximum(abs.(z)) <= 3))
    @printf("ρ=%.1f %-11s h/w=%.2f τ=%.1e: max |z| = %.2f, max |MC − continuum| = %.1e → %s\n", ρ, shape, hw, τ, rows[end].max_abs_z, rows[end].max_abs_diff, rows[end].pass ? "PASS" : "FAIL")
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rows" => rows, "pass" => all(r.pass for r in rows)))
println("saved to $OUT")
```

- [ ] **Step 2: Run it**

Run: `julia --project=research/baseline -t 8 research/phase2/p2_5_overlap_mc_check.jl`
Expected: 8 lines, all `→ PASS`. Linear h/w = 1.0 is the uniform exact case, so there the Monte Carlo should equal e^(−2ρt/w) to roundoff and z ≈ 0. A FAIL elsewhere means the overlap/image handling or the continuum disagrees at that point. Stop and debug.

- [ ] **Step 3: Commit**

```bash
git add research/phase2/p2_5_overlap_mc_check.jl research/results/phase2/p2_5_overlap_mc_check
git commit -m "P2.5.9: Monte Carlo spot check of overlapping layers against the continuum

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Figures, record, tracker

**Files:**
- Create: `research/phase2/p2_5_overlap_figs.jl`
- Modify: `research/phase2_record.md`, `research/project-progress-tracker.md`

- [ ] **Step 1: Create the figures script**

```julia
# P2.5.9 figures from grid.csv / ddep.csv: (1) vs overlap (h/w, and mid/wall for the exponential): max |ΔS|, ΔT2/T2 with
# the SNR-200 2σ band, best-fit residual — one line per ρw/D, panels per profile; (2) D dependence: T2(D)/T2(D = 3) − 1
# for the conventional model and the layer at each overlap, ρ = 0.1 and 1.
#
# Run: julia --project=research/baseline research/phase2/p2_5_overlap_figs.jl

include(joinpath(@__DIR__, "harness.jl"))
using CairoMakie
const OUT = phase2_outdir("p2_5_overlap")
rd(f) = (L = readlines(joinpath(OUT, f)); hdr = split(L[1], ","); [Dict(zip(hdr, split(l, ","))) for l in L[2:end]])
num(d, k) = parse(Float64, d[k])
g = rd("grid.csv"); dd = rd("ddep.csv")
rwds = sort(unique([num(r, "rho_w_over_D") for r in g]))
fig = Figure(size=(1300, 700))
for (i, shape) in enumerate(("linear", "exponential")), (j, (k, lab)) in enumerate((("max_abs_dS", "max |ΔS| (fixed ρ)"), ("rel_T2_TE50", "ΔT₂/T₂ at TE = 50"), ("best_fit_residual", "best-fit conventional residual")))
    ax = Axis(fig[i, j], xlabel="h/w", ylabel=lab, yscale=k == "rel_T2_TE50" ? identity : log10, title="$shape: $lab")
    for rwd in rwds
        q = sort([r for r in g if r["shape"] == shape && num(r, "rho_w_over_D") == rwd], by=r -> num(r, "h_over_w"))
        y = [k == "rel_T2_TE50" ? num(r, k) : max(num(r, k), 1e-9) for r in q]
        scatterlines!(ax, [num(r, "h_over_w") for r in q], y, label="ρw/D = $rwd")
    end
    k == "max_abs_dS" || k == "best_fit_residual" ? hlines!(ax, [0.01, 0.02], color=:gray, linestyle=:dash) : nothing
    shape == "linear" && vlines!(ax, [0.5], color=:black, linestyle=:dot)        # layers touch at h/w = 0.5
    (i, j) == (1, 3) && axislegend(ax, position=:rb, labelsize=8)
end
save(joinpath(OUT, "p2_5_overlap_signal.png"), fig)

fig2 = Figure(size=(900, 380))
for (j, ρ) in enumerate((0.1, 1.0))
    ax = Axis(fig2[1, j], xlabel="h/w", ylabel="T₂(D = 0.5)/T₂(D = 3) − 1", title="D dependence of apparent T₂, ρ = $ρ")
    conv = only([r for r in dd if num(r, "rho") == ρ && r["model"] == "conventional"])
    hlines!(ax, [num(conv, "rel_change_D05_vs_D3")], color=:black, linestyle=:dash, label="conventional")
    for shape in ("linear", "exponential")
        q = sort([r for r in dd if num(r, "rho") == ρ && r["model"] == shape], by=r -> num(r, "h_over_w"))
        scatterlines!(ax, [num(r, "h_over_w") for r in q], [num(r, "rel_change_D05_vs_D3") for r in q], label=shape)
    end
    axislegend(ax, position=:rt, labelsize=8)
end
save(joinpath(OUT, "p2_5_overlap_ddep.png"), fig2)
println("saved to $OUT")
```

- [ ] **Step 2: Run it and look at the figures**

Run: `julia --project=research/baseline research/phase2/p2_5_overlap_figs.jl`
Expected: two figures. If a panel is unreadable, adjust it and record a ruling.

- [ ] **Step 3: Record, tracker, commit**

Add a section "P2.5.9: overlapping layers (linear, exponential) vs the conventional model" to `research/phase2_record.md`. It should contain:
- **Profile check:** the closed forms (linear plateau (ρ/h)(4 − 2w/h) over [w − h, h]; uniform 2ρ/w at h = w; exponential 2(ρ/h)N·e^(−w/2h)cosh((x − w/2)/h), which is lowest mid-gap), the realism result for every grid point, package = closed form, and the step counter-example. Include the profile figure and the overlap-metric table.
- **Monte Carlo spot check:** the table.
- **Comparison tables per profile:** h/w and ρw/D with mid/wall, max |ΔS|, ΔT₂/T₂ against 2σ (SNR 200), the best-fit residual (edge flags) and the mono-exponential residuals.
- **D-dependence table and figure.**
- **A numbers-only conclusion:**
  - does the overlap change the difference from the conventional model (with ρ fixed and with ρ fitted);
  - from which ρw/D the overlap becomes observable (single measurement);
  - whether the D-independence of the fully overlapping (uniform) layer, against the D-dependent conventional decay, is a distinguishing feature, and how large it is at ρ = 0.1 and 1 compared with the T₂ noise.
- **Limitations:**
  - walls only, free decay;
  - overlapping requires h ≳ w/2, so it is relevant for small compartments;
  - D is assumed not to change ρ;
  - single-measurement thresholds.

In `research/project-progress-tracker.md`, add under §2.5 the row below, and increase the Phase 2 total and "Done" counts by 1:

```markdown
| P2.5.9 | ☑ | Overlapping layers (linear, exponential): realistic summed profile (never higher mid-gap), and comparison with the conventional model vs overlap and ρw/D, incl. D dependence | profile check pass/fail; comparison as characterisation | `p2_5_overlap/`, `p2_5_overlap_mc_check/` |
```

Then commit:

```bash
git add research/phase2/p2_5_overlap_figs.jl research/results/phase2/p2_5_overlap research/phase2_record.md research/project-progress-tracker.md
git commit -m "P2.5.9: overlapping-layer results, figures and record

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
