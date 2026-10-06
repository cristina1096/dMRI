# Linear Near-Surface Profile on Walls Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the linear profile g(u) = 2(1 − u) for 0 ≤ u ≤ 1 (proposal §3.2 III) to the wall layer, verify it the same way as the step profile, fit h\*(θ) to the B4 baseline, and compare it with the step profile at fixed ΔR₂(0).

**Architecture:** `LayerSide` gains a `shape` (`:step` or `:linear`). On a straight piece the distance d changes linearly, so g(d/h) of the linear profile is also linear along the piece. Its average over the part inside the layer is therefore exactly its value at the midpoint of that part, which keeps `side_exponent` a closed form with no extra integration. `Walls` gets `layer_shape` (both sides) and `layer_shape_positive` / `layer_shape_negative` (one side), as strings. The research scripts take `SHAPE=linear` from the environment and write to `…_linear` result folders, so the step results stay untouched.

**Tech Stack:** Julia 1.12, MCMRSimulator (this repo), the Phase 2 research code in `research/phase2/`, `Test` stdlib.

**Spec:** proposal (Y. Shi, 2026-09-14) §3.2 III (ΔR₂(d) = (ρ/h)·g(d/h), ∫₀¹ g = 1; linear g(u) = 2(1 − u)); `research/project-progress-tracker.md` P2.4.1, P2.4.2, P2.4.3, P2.4.7, P2.4.8; the user's framing (baseline = reference, not truth; no fast-diffusion arguments).

**Depends on:** plans `2026-10-01-near-surface-layer-walls.md`, `2026-10-03-phase2-harness-and-trace.md` and `2026-10-03-phase2-step-verification.md` (all done; tag `rd-step-v1`).

## Global Constraints

- Work only on branch `cc/near-surface-r2`.
- **Timestep constraint (added after the timestep study, 2026-10-07):** `Layers.max_layer_rate` (used by `TimeStep`'s `layer` option, τ ≤ c/ΔR₂(0)_max) currently returns `surface_rate(side)` = ρ/h, the step-profile ΔR₂(0). In Task 1 make it return the profile's ΔR₂(0) = (ρ/h)·g(0) (= 2ρ/h for linear), with a test: linear side ρ = 0.01, h = 0.1 → 0.2.
- **Step-profile behaviour must stay bit-identical.** All existing package tests (`test_args=["layer"]`) and research tests must still pass unchanged.
- Units as before: ρ in µm/ms, h in µm, ΔR₂ in 1/ms. **ΔR₂(0) = (ρ/h)·g(0)**: g(0) = 1 for step, **2 for linear**. So at fixed ΔR₂(0) the linear profile has **half the ρ** of the step profile: ρ = ΔR₂(0)·h/2.
- Reference configuration as in Phase 2: walls every 2 µm, D = 3 µm²/ms, no RF, readouts 0:5:50 ms, τ = 1e-2 ms (provisional until P2.6), seeds 1–10, 10⁵ spins, ΔR₂(0) = 0.1 ms⁻¹, layer on both faces.
- Uncertainty of a mean is the SEM. The baseline is a reference: V4 residuals are reported, never pass/fail.
- Package tests: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`. Research tests: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## The exact check for the linear profile

The step profile's exact case (h = w/2: every point in exactly one layer) does **not** hold for the linear profile, because ΔR₂ varies inside the layer. The linear profile has its own exact case instead: **h = w with layers on both faces**. At position u in the gap (walls at u = 0 and u = w):

  ΔR₂(u) = (2ρ/w)·(1 − u/w) + (2ρ/w)·(1 − (w − u)/w) = 2ρ/w, the same everywhere.

Since ΔR₂(0) = (ρ/h)·g(0) = 2ρ/w when h = w, the uniform rate **equals ΔR₂(0)**: every spin decays at exactly ΔR₂(0), whatever its path. With ΔR₂(0) = 0.1 ms⁻¹ (ρ = 0.1 µm/ms): S(t) = e^(−0.1·t). This case uses overlapping layers (h > w/2, allowed up to h = spacing), so it also tests the overlap sum (Eq. 10) end to end.

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `src/geometries/internal/layers.jl` | Modify | `LayerSide.shape`, `profile_g`, shape-aware `surface_rate` and `side_exponent` |
| `src/geometries/user/obstructions/obstructions.jl` | Modify | `layer_shape`, `layer_shape_positive`, `layer_shape_negative` fields on `Walls` |
| `src/geometries/user/fix.jl` | Modify | resolve and validate the shape in `wall_layers` |
| `src/cli/geometry.jl` | Modify | CLI parsing of string-valued fields |
| `test/test_layer.jl`, `test/test_cli.jl` | Modify | package tests |
| `research/phase2/harness.jl` | Modify | `rho_for(...; shape)`, `layer_walls(...; shape)`, `H_GRID_LINEAR` |
| `research/phase2/p2_2_9_layer_trace_lib.jl` | Modify | `sampled_profile_average` (brute force, any shape) |
| `research/phase2/test/test_harness.jl`, `test_trace.jl` | Modify | research tests |
| `research/phase2/p2_2_9_layer_trace.jl` | Modify | `SHAPE` env, profile-aware prediction |
| `research/phase2/p2_4_3_linear_exact.jl` | Create | V0b, exact h = w case, additivity for the linear profile |
| `research/phase2/p2_3_4_h_sweep.jl`, `p2_3_9_fit_h.jl`, `p2_3_9b_gauss_newton.jl` | Modify | `SHAPE` env, grid and output folder per shape |
| `research/phase2/p2_4_8_compare_step_linear.jl` | Create | step vs linear comparison at fixed ΔR₂(0) |
| `research/phase2_record.md`, `research/project-progress-tracker.md` | Modify | records, ticks |

## Review Focus

1. **The step profile must not change.** The refactor of `side_exponent` must give the same numbers for `:step`. Test: Task 1 keeps every existing step assertion, and the full `layer` package suite must stay green.
2. **The linear profile's edges:** g = 2 at the wall (d = 0) and g = 0 at d = h. A spin sitting exactly at d = 0 must get 2ρ/h; one at d = h must get 0. Test: Task 1.
3. **Different shapes on the two sides** (positive linear, negative step) must each use their own formula. Test: Task 2's `layer_exponent` case.
4. **A string-valued field** must survive JSON and the CLI. Test: Task 2 (JSON round trip and a `geometry create walls … --layer_shape linear` CLI call).
5. **Overlap up to h = w** for the linear profile must sum to a uniform rate. Test: Task 2 (`layer_exponent`) and Task 4 (full size, per spin).

---

### Task 1: Linear profile in the layer maths (P2.4.1, P2.4.2)

**Files:**
- Modify: `src/geometries/internal/layers.jl`
- Modify: `test/test_layer.jl` (append)

**Interfaces:**
- Produces:
  - `Layers.LayerSide(rho, h, shape::Symbol)`. The existing `LayerSide(rho, h)` still works and means `:step`.
  - `Layers.profile_g(shape::Symbol, u)`: the profile value, 0 outside 0 ≤ u ≤ 1
  - `Layers.surface_rate(side)` = (ρ/h)·g(0)
  - `Layers.side_exponent(side, d0, d1, dt)` for `:step` and `:linear`

- [ ] **Step 1: Append the failing tests to `test/test_layer.jl`**

```julia
@testset "test_layer.jl: linear profile maths" begin
    Layers = mr.Geometries.Internal.Layers

    @testset "profile values and normalisation (P2.4.2)" begin
        @test Layers.profile_g(:linear, 0.0) == 2.0
        @test Layers.profile_g(:linear, 0.5) == 1.0
        @test Layers.profile_g(:linear, 1.0) == 0.0
        @test Layers.profile_g(:linear, 1.5) == 0.0
        @test Layers.profile_g(:step, 0.3) == 1.0
        @test Layers.profile_g(:step, 1.2) == 0.0
        n = 100_000
        for shape in (:step, :linear)
            @test sum(Layers.profile_g(shape, (k - 0.5) / n) for k in 1:n) / n ≈ 1.0 atol=1e-10
        end
    end

    @testset "step is unchanged" begin
        side = Layers.LayerSide(0.05, 0.5)
        @test side.shape == :step
        @test Layers.side_exponent(side, 0.1, 0.3, 0.2) ≈ 0.1 * 0.2
        @test Layers.side_exponent(side, 0.4, 0.6, 1.0) ≈ 0.1 * 1.0 * 0.5
    end

    @testset "linear: surface rate and edges" begin
        side = Layers.LayerSide(0.05, 0.5, :linear)            # ρ/h = 0.1, ΔR2(0) = 0.2
        @test Layers.surface_rate(side) ≈ 0.2
        @test Layers.side_exponent(side, 0.0, 0.0, 1.0) ≈ 0.2  # sitting on the wall
        @test Layers.side_exponent(side, 0.5, 0.5, 1.0) ≈ 0.0  # sitting at d = h
        @test Layers.side_exponent(side, 0.7, 0.7, 1.0) == 0.0 # outside
    end

    @testset "linear: pieces" begin
        side = Layers.LayerSide(0.05, 0.5, :linear)            # ρ/h = 0.1
        # whole layer crossed, 0 → h: average g = 1 → exponent = (ρ/h)·dt
        @test Layers.side_exponent(side, 0.0, 0.5, 1.0) ≈ 0.1
        # 0.8 → 0.2: inside for the second half (d 0.5 → 0.2, mean d = 0.35, g = 2(1 − 0.7) = 0.6)
        @test Layers.side_exponent(side, 0.8, 0.2, 1.0) ≈ 0.1 * 0.6 * 0.5
        # matches brute-force sampling on random pieces
        Random.seed!(5)
        for _ in 1:200
            d0, d1, h = 1.2 * rand() - 0.1, 1.2 * rand() - 0.1, 0.1 + 0.9 * rand()
            s = Layers.LayerSide(0.3 * h, h, :linear)
            n = 200_000
            brute = sum(Layers.profile_g(:linear, (d0 + (d1 - d0) * (k - 0.5) / n) / h) for k in 1:n) / n * (s.rho / h)
            @test Layers.side_exponent(s, d0, d1, 1.0) ≈ brute atol=2e-5
        end
    end
end
```

`profile_g` returns 0 for u < 0, because points with d < 0 are on the other side of the wall. That makes the brute-force sum correct for pieces that cross d = 0.

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`
Expected: FAIL with `UndefVarError: profile_g not defined`.

- [ ] **Step 3: Update `src/geometries/internal/layers.jl`**

Replace the `LayerSide` struct and `surface_rate`:

```julia
"""
    LayerSide(rho, h[, shape])

Near-surface layer on one side of a surface: integrated relaxivity `rho` (um/ms), length scale `h` (um) and
profile `shape` (`:step` or `:linear`; default `:step`). A side with `rho == 0` has no layer.
"""
struct LayerSide
    rho :: Float64
    h :: Float64
    shape :: Symbol
end
LayerSide(rho, h) = LayerSide(rho, h, :step)

"""
    profile_g(shape, u)

Normalised profile g(u), ∫₀¹ g(u) du = 1 (proposal §3.2 III). Zero outside 0 ≤ u ≤ 1.
`:step`: g = 1. `:linear`: g = 2(1 − u).
"""
function profile_g(shape::Symbol, u::Real)
    (0 <= u <= 1) || return 0.0
    shape == :step && return 1.0
    shape == :linear && return 2 * (1 - u)
    error("Unknown layer profile $shape")
end

"""
    surface_rate(side)

Excess transverse relaxation rate at the surface, ΔR2(0) = (ρ/h)·g(0), in 1/ms.
"""
surface_rate(side::LayerSide) = iszero(side.rho) ? 0.0 : side.rho / side.h * profile_g(side.shape, 0.0)
```

Replace `side_exponent`:

```julia
"""
    side_exponent(side, d0, d1, dt)

∫ΔR2(d(t)) dt over a straight segment of duration `dt` (ms), where `d0` and `d1` are the distances from the
surface at the start and end, positive into this side. d is linear along the segment, so for both profiles
the integral is exact in closed form:
- `:step`: ΔR2(0)·(time spent with 0 ≤ d ≤ h);
- `:linear`: g(d/h) is linear along the segment, so its average over the part inside the layer is its value
  at the midpoint of that part.
"""
function side_exponent(side::LayerSide, d0::Float64, d1::Float64, dt::Float64)
    iszero(side.rho) && return 0.0
    (s_enter, s_exit) = segment_overlap(d0, d1, 0.0, side.h)
    fraction = max(0.0, s_exit - s_enter)
    iszero(fraction) && return 0.0
    side.shape == :step && return surface_rate(side) * dt * fraction
    d_mid = d0 + (d1 - d0) * (s_enter + s_exit) / 2
    return side.rho / side.h * profile_g(side.shape, d_mid / side.h) * dt * fraction
end
```

For `:step` this returns exactly the same numbers as before: `surface_rate` is ρ/h·1.0, and a zero fraction still gives 0.0. A stationary spin inside the layer gets `(s_enter, s_exit) = (0, 1)` from `segment_overlap` and d_mid = d0.

- [ ] **Step 4: Run the tests and check they pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`
Expected: PASS. All earlier `test_layer.jl` testsets plus "linear profile maths".

- [ ] **Step 5: Commit**

```bash
git add src/geometries/internal/layers.jl test/test_layer.jl
git commit -m "Add the linear near-surface profile to the layer maths

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: `layer_shape` on `Walls`, validation, CLI and JSON (P2.4.1)

**Files:**
- Modify: `src/geometries/user/obstructions/obstructions.jl` (the `:Wall` fields), `src/geometries/user/fix.jl` (`wall_layers`), `src/cli/geometry.jl` (string parsing)
- Modify: `test/test_layer.jl`, `test/test_cli.jl`

**Interfaces:**
- Consumes: `Layers.LayerSide(rho, h, shape)` (Task 1).
- Produces: `Walls` fields `layer_shape::String` (default `"step"`), `layer_shape_positive`, `layer_shape_negative` (default `nothing` = inherit). Allowed values are `"step"` and `"linear"`.

- [ ] **Step 1: Append the failing tests**

To `test/test_layer.jl`:

```julia
@testset "test_layer.jl: layer_shape on Walls" begin
    Layers = mr.Geometries.Internal.Layers
    geom(walls) = mr.Simulation([]; geometry=walls, verbose=false).geometry
    v(x) = SVector{3, Float64}(x, 0.3, -1.2)

    @testset "resolution" begin
        l = geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5))[1].layer[1]
        @test (l.positive.shape, l.negative.shape) == (:step, :step)
        l = geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5, layer_shape="linear"))[1].layer[1]
        @test (l.positive.shape, l.negative.shape) == (:linear, :linear)
        l = geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5, layer_shape_negative="linear"))[1].layer[1]
        @test (l.positive.shape, l.negative.shape) == (:step, :linear)
        @test_throws ErrorException geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5, layer_shape="gaussian"))
    end

    @testset "mixed shapes use their own formula" begin
        # positive side step, negative side linear; ρ/h = 0.1
        g = geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5, layer_shape_negative="linear"))
        @test Layers.layer_exponent(g, v(0.1), v(0.1), 1.0) ≈ 0.1            # positive (step) side of wall 0: g = 1
        @test Layers.layer_exponent(g, v(1.9), v(1.9), 1.0) ≈ 0.16           # negative (linear) side of wall 2, d = 0.1: g = 1.6
        @test Layers.layer_exponent(g, v(0.0), v(0.0), 1.0) ≈ 0.1 + 0.2      # on wall 0: inside both of its layers (step 1 + linear 2)
    end

    @testset "linear, h = w: uniform rate" begin
        g = geom(mr.Walls(repeats=2., layer_rho=0.1, layer_h=2.0, layer_shape="linear"))   # 2ρ/w = 0.1
        for (a, b) in ((0.05, 0.3), (0.8, 1.2), (1.0, 1.9), (0.01, 1.99), (0.3, 0.3))
            @test Layers.layer_exponent(g, v(a), v(b), 1.0) ≈ 0.1 rtol=1e-12
        end
    end

    @testset "linear, h = w: every spin decays at exactly 0.1 /ms" begin
        walls = mr.Walls(repeats=2., layer_rho=0.1, layer_h=2.0, layer_shape="linear")
        sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=walls, diffusivity=3., verbose=false)
        Random.seed!(1)
        snap = mr.Snapshot(2000, sim, 500; transverse=1., longitudinal=0.)
        res = mr.readout(snap, sim, [10., 50.]; return_snapshot=true)
        for (t, s) in zip((10., 50.), res)
            @test all(isapprox.([x.orientations[1].transverse for x in s.spins], exp(-0.1 * t); rtol=1e-10))
        end
    end

    @testset "JSON round trip keeps layer_shape" begin
        io = IOBuffer()
        mr.write_geometry(io, mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5, layer_shape="linear"))
        back = mr.read_geometry_json(String(take!(io)))
        @test back.layer_shape.value == "linear"
        @test isnothing(back.layer_shape_positive.value)
    end
end
```

To `test/test_cli.jl`, inside `@testset "mcmr geometry create"`, after the `"create walls"` testset:

```julia
        @testset "create walls with a linear layer" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create walls 1 test.json --repeats 2 --layer_rho 0.05 --layer_h 0.5 --layer_shape linear")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["layer_shape"] == "linear"
                @test result["layer_rho"] == 0.05
            end
        end
```

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer", "cli"])'`
Expected: FAIL. `Walls` has no keyword `layer_shape`, and the CLI test fails on the unknown option.

- [ ] **Step 3: Add the fields**

In `src/geometries/user/obstructions/obstructions.jl`, inside the `ObstructionType(:Wall; …, fields=[…])` list, after the `layer_h_negative` field:

```julia
        Field{String}(:layer_shape, "Near-surface layer: profile g on both sides of the wall, unless a side-specific value is set. \"step\" or \"linear\".", "step"),
        Field{String}(:layer_shape_positive, "Near-surface layer: profile on the positive side. Overrides `layer_shape`."),
        Field{String}(:layer_shape_negative, "Near-surface layer: profile on the negative side. Overrides `layer_shape`."),
```

- [ ] **Step 4: Resolve and validate the shape in `src/geometries/user/fix.jl`**

In `wall_layers`, add a helper next to `resolve` and use it when building each side. Replace

```julia
            Internal.LayerSide(rho, h)
```

with

```julia
            Internal.LayerSide(rho, h, resolve_shape(side, i))
```

and add, directly after the `resolve` helper inside `wall_layers`:

```julia
    function resolve_shape(side::String, i::Int)
        specific = getproperty(walls, Symbol("layer_shape_" * side))[i]
        shape = isnothing(specific) ? walls.layer_shape[i] : specific
        shape = isnothing(shape) ? "step" : shape
        shape in ("step", "linear") || error("Wall $i, $side side: layer_shape must be \"step\" or \"linear\", got \"$shape\".")
        return Symbol(shape)
    end
```

- [ ] **Step 5: Parse string fields in the CLI (`src/cli/geometry.jl`)**

After the existing `ArgParse.parse_item(::Type{FieldParser{T}}, …) where {T<:Number}` method, add:

```julia
function ArgParse.parse_item(::Type{FieldParser{String}}, text::AbstractString)
    parts = split(text, ',')
    return length(parts) == 1 ? FieldParser{String}(String(text)) : FieldParser{String}(String.(parts))
end
```

- [ ] **Step 6: Run the tests and check they pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer", "cli", "collisions", "transfer"])'`
Expected: PASS. If the CLI test still fails because ArgParse cannot build a default for `FieldParser{String}`, record a ruling and pass the default through `FieldParser{String}("step")` explicitly. Do not change the field type.

- [ ] **Step 7: Commit**

```bash
git add src/geometries/user/obstructions/obstructions.jl src/geometries/user/fix.jl src/cli/geometry.jl test/test_layer.jl test/test_cli.jl
git commit -m "Add layer_shape (step, linear) to Walls with CLI and JSON support

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Research harness, trace and the linear trace run (P2.4.3, part 1)

**Files:**
- Modify: `research/phase2/harness.jl`, `research/phase2/p2_2_9_layer_trace_lib.jl`, `research/phase2/p2_2_9_layer_trace.jl`
- Modify: `research/phase2/test/test_harness.jl`, `research/phase2/test/test_trace.jl`
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `Walls(…; layer_shape)` (Task 2).
- Produces:
  - `rho_for(rate, h; shape="step")` = rate·h / g(0), g(0) = 1 (step) or 2 (linear)
  - `layer_walls(; h, rate, side, h_neg, rate_neg, shape="step")`
  - `H_GRID_LINEAR`
  - `sampled_profile_average(a, b, h; shape=:step, n=10_000)`: the brute-force time-average of g(d/h) summed over both faces of one gap
  - `p2_2_9_layer_trace.jl` reads `SHAPE` (default `"step"`) and writes to `p2_2_9_layer_trace` (step) or `p2_2_9_layer_trace_linear`

- [ ] **Step 1: Append the failing research tests**

To `research/phase2/test/test_harness.jl`:

```julia
@testset "harness: shapes" begin
    @test rho_for(0.1, 0.5; shape="linear") == 0.025
    @test rho_for(0.1, 0.5) == 0.05
    l = Simulation([]; geometry=layer_walls(h=0.5, rate=0.1, shape="linear"), verbose=false).geometry[1].layer[1]
    @test (l.positive.shape, l.positive.rho) == (:linear, 0.025)
    @test H_GRID_LINEAR[end] == 2.0
    r = run_walls(layer_walls(h=2.0, rate=0.1, shape="linear"); times=[0.0, 50.0], seeds=1:2, nspins=1000)
    @test all(isapprox.(r.S, exp.(-0.1 .* [0.0 50.0]); rtol=1e-10))            # exact h = w case
end
```

To `research/phase2/test/test_trace.jl`:

```julia
@testset "trace: brute-force profile average" begin
    @test sampled_profile_average(0.1, 0.3, 0.5) ≈ sampled_layer_fraction(0.1, 0.3, 0.5)      # step: same as before
    @test sampled_profile_average(0.0, 0.5, 0.5; shape=:linear) ≈ 1.0 atol=1e-3              # g averages to 1 across the layer
    @test sampled_profile_average(0.3, 1.7, 2.0; shape=:linear) ≈ 2.0 atol=1e-12             # h = w: g(lower) + g(upper) = 2 everywhere
end
```

- [ ] **Step 2: Run the research tests and check they fail**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: FAIL with a `MethodError` for `rho_for` with keyword `shape`, or `UndefVarError: sampled_profile_average`.

- [ ] **Step 3: Update `research/phase2/harness.jl`**

Replace `rho_for` and `layer_walls`, and add `H_GRID_LINEAR` after `H_GRID`:

```julia
"Layer thicknesses (µm) for the linear-profile sweep: as H_GRID, extended to h = w (the linear exact case)."
const H_GRID_LINEAR = [0.0, 0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0]

"g(0) of a profile: 1 for step, 2 for linear."
profile_g0(shape) = shape == "step" ? 1.0 : shape == "linear" ? 2.0 : error("Unknown shape $shape")

"ρ (µm/ms) giving surface excess rate ΔR2(0) = `rate` (1/ms) for thickness `h` (µm): ρ = rate·h / g(0)."
rho_for(rate, h; shape="step") = Float64(rate * h / profile_g0(shape))

"""
    layer_walls(; h=0.0, rate=0.0, side=:both, h_neg=h, rate_neg=rate, shape="step")

Reference walls with a layer of profile `shape` ("step" or "linear") and surface excess rate `rate`.
`side` is `:both`, `:positive`, `:negative`, or `:asymmetric` (positive side `h`, `rate`; negative side `h_neg`, `rate_neg`).
"""
function layer_walls(; h=0.0, rate=0.0, side=:both, h_neg=h, rate_neg=rate, shape="step")
    ρ(r, hh) = rho_for(r, hh; shape=shape)
    side == :both && return Walls(repeats=W_REF, layer_rho=ρ(rate, h), layer_h=h, layer_shape=shape)
    side == :positive && return Walls(repeats=W_REF, layer_rho_positive=ρ(rate, h), layer_h_positive=h, layer_shape=shape)
    side == :negative && return Walls(repeats=W_REF, layer_rho_negative=ρ(rate, h), layer_h_negative=h, layer_shape=shape)
    side == :asymmetric && return Walls(repeats=W_REF,
        layer_rho_positive=ρ(rate, h), layer_h_positive=h,
        layer_rho_negative=ρ(rate_neg, h_neg), layer_h_negative=h_neg, layer_shape=shape)
    error("Unknown side $side; use :both, :positive, :negative or :asymmetric.")
end
```

- [ ] **Step 4: Add `sampled_profile_average` to `research/phase2/p2_2_9_layer_trace_lib.jl`**

```julia
"""
    sampled_profile_average(a, b, h; shape=:step, n=10_000)

Time-average over a straight piece x = a → b of g(d/h) summed over both faces of the gap (walls at multiples
of W_REF), by brute-force sampling of `n` midpoints. Independent of MCMR's `Layers` code. Valid for h ≤ W_REF
(only the two faces of the gap can reach it).
"""
function sampled_profile_average(a::Float64, b::Float64, h::Float64; shape=:step, n=10_000)
    g(u) = 0 <= u <= 1 ? (shape == :step ? 1.0 : 2 * (1 - u)) : 0.0
    total = 0.0
    for k in 1:n
        u = mod(a + (b - a) * (k - 0.5) / n, W_REF)
        total += g(u / h) + g((W_REF - u) / h)
    end
    return total / n
end
```

- [ ] **Step 5: Run the research tests and check they pass**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: PASS.

- [ ] **Step 6: Make the trace script shape-aware**

In `research/phase2/p2_2_9_layer_trace.jl`, replace

```julia
const OUT = phase2_outdir("p2_2_9_layer_trace")
```

with

```julia
const SHAPE = get(ENV, "SHAPE", "step")
const OUT = phase2_outdir(SHAPE == "step" ? "p2_2_9_layer_trace" : "p2_2_9_layer_trace_" * SHAPE)
```

Replace

```julia
sim = Simulation(seq; geometry=layer_walls(h=H, rate=RATE), diffusivity=D_REF, R2=R2B, verbose=false)
```

with

```julia
sim = Simulation(seq; geometry=layer_walls(h=H, rate=RATE, shape=SHAPE), diffusivity=D_REF, R2=R2B, verbose=false)
const RHO_OVER_H = rho_for(RATE, H; shape=SHAPE) / H
```

and replace

```julia
            expo += RATE * dt * sampled_layer_fraction(path[k][1], path[k + 1][1], H)
```

with

```julia
            expo += RHO_OVER_H * dt * sampled_profile_average(path[k][1], path[k + 1][1], H; shape=Symbol(SHAPE))
```

In the `write_json` params add `"shape" => SHAPE`. For the step profile, RHO_OVER_H·average = RATE·fraction, which is the same prediction as before.

- [ ] **Step 7: Run the linear trace**

Run: `SHAPE=linear julia --project=research/baseline -t 1 research/phase2/p2_2_9_layer_trace.jl`
Expected: `→ PASS`, with max error ≤ 1e-6 and a positive count of steps with a reflection inside the layer. About 15 s.

- [ ] **Step 8: Record and commit**

Append a summary row and a short section to `research/phase2_record.md`, filled from `results/phase2/p2_2_9_layer_trace_linear/summary.json`:

```markdown
| P2.4.3 (trace) | Single-spin trace, linear profile | every step: \|M⊥_sim/M⊥_pred − 1\| ≤ 1e-6; ≥ 1 reflection inside the layer | <PASS/FAIL>: max error <max_rel_error>; <n> steps with a reflection inside the layer |
```

```markdown
## P2.4.3: linear profile — single-spin trace

Script [p2_2_9_layer_trace.jl](phase2/p2_2_9_layer_trace.jl) with `SHAPE=linear` → [results/phase2/p2_2_9_layer_trace_linear/](results/phase2/p2_2_9_layer_trace_linear/). Same configuration as P2.2.9 (h = 0.5 µm, ΔR₂(0) = 0.1 ms⁻¹ so ρ = 0.025 µm/ms, R₂_bulk = 1/80, τ = 0.01 ms, 100 spins × 500 steps). The prediction samples g(d/h) = 2(1 − d/h) along every piece, independently of the `Layers` code.

| Max relative error per step | Steps with a reflection inside the layer | Outcome |
|---|---|---|
| <max_rel_error> | <n> | <PASS/FAIL> |
```

```bash
git add research/phase2/harness.jl research/phase2/p2_2_9_layer_trace_lib.jl research/phase2/p2_2_9_layer_trace.jl research/phase2/test/test_harness.jl research/phase2/test/test_trace.jl research/results/phase2/p2_2_9_layer_trace_linear research/phase2_record.md
git commit -m "Add linear-profile support to the Phase 2 harness and trace

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Exact checks for the linear profile (P2.4.3, part 2)

**Files:**
- Create: `research/phase2/p2_4_3_linear_exact.jl`
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `run_walls(...; snapshots=true)`, `layer_walls(...; shape)`, `expected_half_gap`, `roundoff_bound` (the exact uniform-rate formula applies to the h = w case too).
- Produces: `results/phase2/p2_4_3_linear_exact/summary.json` with `V0b`, `uniform_bulk_off`, `uniform_bulk_on`, `additivity`, each `{pass, …}`.

- [ ] **Step 1: Create the script**

```julia
# P2.4.3 — exact checks for the linear profile, full size (seeds 1–10, NSPINS spins), τ = 1e-2 ms.
#
# V0b       : linear, ρ = 0 (h = 0.5) vs plain walls, R2_bulk = 1/80: positions and per-spin M⊥ bit-identical.
# uniform   : linear, h = w = 2 µm, both faces, ΔR2(0) = 0.1 → ΔR2(u) = 2ρ/w = 0.1 everywhere in the gap, so every spin
#             decays at exactly 0.1 (+ R2_bulk): per spin ≤ 1e-10, ensemble ≤ roundoff bound. Bulk off and on.
# additivity: linear, h = 0.2, R2_bulk 0 vs 1/80, same seeds: positions ==, per spin ≤ 1e-10, ensemble ≤ roundoff bound.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_4_3_linear_exact.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "expected.jl"))

const OUT = phase2_outdir("p2_4_3_linear_exact")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(NSPINS_REF)))
const TOL_SPIN = 1e-10
per_spin_final(r) = [[s.orientations[1].transverse for s in snap.spins] for snap in r.snapshots]
positions_final(r) = [mr.position.(snap) for snap in r.snapshots]

results = Dict{String, Any}()
function report(name, pass, extra)
    results[name] = merge(Dict{String, Any}("pass" => pass), extra)
    @printf("%-18s → %s  %s\n", name, pass ? "PASS" : "FAIL", string(extra))
end

a = run_walls(Walls(repeats=W_REF); R2_bulk=1 / 80, nspins=NSPINS, snapshots=true)
b = run_walls(layer_walls(h=0.5, rate=0.0, shape="linear"); R2_bulk=1 / 80, nspins=NSPINS, snapshots=true)
report("V0b", positions_final(a) == positions_final(b) && per_spin_final(a) == per_spin_final(b) && a.S == b.S, Dict("R2_bulk" => 1 / 80))

bound = roundoff_bound(TIMES_REF, TAU_REF, NSPINS)
for (label, R2) in (("uniform_bulk_off", 0.0), ("uniform_bulk_on", 1 / 80))
    r = run_walls(layer_walls(h=W_REF, rate=RATE_REF, shape="linear"); R2_bulk=R2, nspins=NSPINS, snapshots=true)
    target = expected_half_gap(TIMES_REF, RATE_REF; R2_bulk=R2)          # same formula: uniform rate RATE_REF + R2
    spin_dev = maximum(maximum(abs.(m ./ target[end] .- 1)) for m in per_spin_final(r))
    ens_dev = maximum(abs.(r.S ./ target' .- 1))
    report(label, spin_dev <= TOL_SPIN && ens_dev <= bound,
        Dict("max_rel_dev_per_spin" => spin_dev, "max_rel_dev_ensemble" => ens_dev, "bound" => bound, "R2_bulk" => R2))
end

lo = run_walls(layer_walls(h=0.2, rate=RATE_REF, shape="linear"); nspins=NSPINS, snapshots=true)
wb = run_walls(layer_walls(h=0.2, rate=RATE_REF, shape="linear"); R2_bulk=1 / 80, nspins=NSPINS, snapshots=true)
same_pos = positions_final(lo) == positions_final(wb)
spin_dev = maximum(maximum(abs.(x ./ (y .* exp(-TIMES_REF[end] / 80)) .- 1)) for (x, y) in zip(per_spin_final(wb), per_spin_final(lo)))
ens_dev = maximum(abs.(wb.S ./ (lo.S .* exp.(-TIMES_REF' ./ 80)) .- 1))
report("additivity", same_pos && spin_dev <= TOL_SPIN && ens_dev <= bound,
    Dict("positions_identical" => same_pos, "max_rel_dev_per_spin" => spin_dev, "max_rel_dev_ensemble" => ens_dev, "bound" => bound))

write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "shape" => "linear", "nspins" => NSPINS, "results" => results))
println("saved to $OUT")
```

- [ ] **Step 2: Run it**

Run: `julia --project=research/baseline -t 8 research/phase2/p2_4_3_linear_exact.jl`
Expected: four lines, all `→ PASS`. About 10 min. A FAIL is a code bug (these are exact identities): use superpowers:systematic-debugging.

- [ ] **Step 3: Record and commit**

Summary rows and a section in `research/phase2_record.md`, from `summary.json`, in the same layout as "P2.3.1, P2.3.4(a), P2.3.11: exact checks at full size". The section text must state the exact case: with h = w on both faces, ΔR₂(u) = 2ρ/w is uniform, so S = e^(−ΔR₂(0)·t) for every spin.

```bash
git add research/phase2/p2_4_3_linear_exact.jl research/results/phase2/p2_4_3_linear_exact research/phase2_record.md
git commit -m "Add exact checks for the linear profile

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Linear h sweep, monotonicity and decay shape (P2.4.3, part 3)

**Files:**
- Modify: `research/phase2/p2_3_4_h_sweep.jl`
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `layer_walls(...; shape)`, `H_GRID`, `H_GRID_LINEAR`, `curvature`, `increasing_steps`, `sweep_nspins`.
- Produces: `results/phase2/p2_3_4_h_sweep_linear/{sweep.csv, summary.json, p2_3_4_h_sweep.png}`. The step folder is unchanged.

- [ ] **Step 1: Make the sweep script shape-aware**

In `research/phase2/p2_3_4_h_sweep.jl`, replace

```julia
const OUT = phase2_outdir("p2_3_4_h_sweep")
```

with

```julia
const SHAPE = get(ENV, "SHAPE", "step")
const GRID = SHAPE == "step" ? H_GRID : H_GRID_LINEAR
const OUT = phase2_outdir(SHAPE == "step" ? "p2_3_4_h_sweep" : "p2_3_4_h_sweep_" * SHAPE)
```

Then replace every remaining `H_GRID` in the script with `GRID`, and replace `layer_walls(h=h, rate=RATE_REF)` with `layer_walls(h=h, rate=RATE_REF, shape=SHAPE)`. In the `write_json` call, change `"h_grid" => H_GRID` to `"h_grid" => GRID` and add `"shape" => SHAPE`. In the first axis title, replace `"S(t) per h (ΔR₂(0) = 0.1 /ms)"` with `"S(t) per h, $(SHAPE) (ΔR₂(0) = 0.1 /ms)"`.

- [ ] **Step 2: Run the linear sweep**

Run: `SHAPE=linear julia --project=research/baseline -t 8 research/phase2/p2_3_4_h_sweep.jl`
Expected: 16 `h = …` lines. S(50) goes from 1.0 (h = 0) to about e⁻⁵ = 0.0067 (h = 2.0, the exact case). The run should print `P2.3.5 monotonicity: 15 of 15 … → PASS`. About 40–60 min. Also check that 1 − S(50) spans the B4 range 0.048–0.912 (needed for Task 6). If it doesn't, record a ruling and add grid points before Task 6.

- [ ] **Step 3: Record and commit**

Summary row and a section in `research/phase2_record.md`, from `results/phase2/p2_3_4_h_sweep_linear/summary.json`, with the same columns as the step sweep table (h, S(25), S(50) ± SEM, 1 − S(50), curvature ± SEM). The interpretation paragraph must be numbers only: monotonicity, whether the curvature is zero within noise, and the coverage of the B4 range.

```bash
git add research/phase2/p2_3_4_h_sweep.jl research/results/phase2/p2_3_4_h_sweep_linear research/phase2_record.md
git commit -m "Run the h sweep for the linear profile

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: V4 for the linear profile: fit h\*(θ) to B4 (P2.4.7, linear)

**Files:**
- Modify: `research/phase2/p2_3_9_fit_h.jl`, `research/phase2/p2_3_9b_gauss_newton.jl`
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `results/phase2/p2_3_4_h_sweep_linear/sweep.csv` (Task 5); the three-stage fit (`fit_h` → parabola → Gauss–Newton) from the step V4.
- Produces: `results/phase2/p2_3_9_fit_h_linear/{summary.json, h_star.csv, h_star_final.csv, p2_3_9_fit_h.png}`.

- [ ] **Step 1: Make both fit scripts shape-aware**

In **both** `p2_3_9_fit_h.jl` and `p2_3_9b_gauss_newton.jl`, replace

```julia
const OUT = phase2_outdir("p2_3_9_fit_h")
const SWEEP = joinpath(REPO_ROOT, "research", "results", "phase2", "p2_3_4_h_sweep", "sweep.csv")
```

with

```julia
const SHAPE = get(ENV, "SHAPE", "step")
const SUFFIX = SHAPE == "step" ? "" : "_" * SHAPE
const OUT = phase2_outdir("p2_3_9_fit_h" * SUFFIX)
const SWEEP = joinpath(REPO_ROOT, "research", "results", "phase2", "p2_3_4_h_sweep" * SUFFIX, "sweep.csv")
```

and replace every `layer_walls(h=…, rate=RATE_REF)` call with the same call plus `, shape=SHAPE`. In `p2_3_9_fit_h.jl` add `"shape" => SHAPE` to the `write_json` dictionary, and in the first axis title replace `"V4: h*(θ), ΔR₂(0) = 0.1 /ms, τ = 1e-2"` with `"V4 ($(SHAPE)): h*(θ), ΔR₂(0) = 0.1 /ms, τ = 1e-2"`. Make the same title change in `p2_3_9b_gauss_newton.jl`.

- [ ] **Step 2: Check that the step defaults are unchanged**

Run: `julia --project=research/baseline -e 'include("research/phase2/harness.jl"); println(get(ENV, "SHAPE", "step"))'`
Expected: `step`. With `SHAPE` unset, both scripts use the old folders. Don't rerun the step fit.

- [ ] **Step 3: Run the linear fit (both stages)**

Run: `SHAPE=linear julia --project=research/baseline -t 8 research/phase2/p2_3_9_fit_h.jl`, then `SHAPE=linear julia --project=research/baseline -t 8 research/phase2/p2_3_9b_gauss_newton.jl`
Expected: 12 `θ = …` lines from each. After the Gauss–Newton step, the residuals should be mixed in sign and small, as for the step profile. About 2.5–3 h in total. If any θ's residuals still share one sign after the Gauss–Newton step, run `p2_3_9b_gauss_newton.jl` once more (it starts from the current `h_star`) and record that as a ruling.

- [ ] **Step 4: Record and commit**

Summary row and a section in `research/phase2_record.md`, from `results/phase2/p2_3_9_fit_h_linear/`, with the same columns as the step V4 table (θ, h\*, σ_jack, h\*/θ, χ²/dof, max |residual|, max |S_model/S_ref − 1|, S(50) reference and model). Add the timestep and bulk-on lines in the same form as the step record. The interpretation must be numbers only, with no formula-based expectation.

```bash
git add research/phase2/p2_3_9_fit_h.jl research/phase2/p2_3_9b_gauss_newton.jl research/results/phase2/p2_3_9_fit_h_linear research/phase2_record.md
git commit -m "Fit h*(theta) of the linear profile to the B4 reference

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Step vs linear at fixed ΔR₂(0), tracker (P2.4.8, linear part)

**Files:**
- Create: `research/phase2/p2_4_8_compare_step_linear.jl`
- Modify: `research/phase2_record.md`, `research/project-progress-tracker.md`

**Interfaces:**
- Consumes: the `summary.json` files of `p2_3_4_h_sweep`, `p2_3_4_h_sweep_linear`, `p2_3_9_fit_h`, `p2_3_9_fit_h_linear`.
- Produces: `results/phase2/p2_4_8_compare_step_linear/{comparison.csv, summary.json, p2_4_8_compare.png}`.

- [ ] **Step 1: Create the comparison script**

```julia
# P2.4.8 (step vs linear) — profile comparison at fixed surface excess rate ΔR2(0) = 0.1 /ms. Characterisation only.
# (1) Attenuation 1 − S(50) vs h for both profiles (from the two sweeps).
# (2) V4: h*(θ) and ρ*(θ) = ΔR2(0)·h*/g(0) for both profiles, their ratios, and the max |residual| of each fit.
#
# Run: julia --project=research/baseline research/phase2/p2_4_8_compare_step_linear.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_4_8_compare_step_linear")
R = joinpath(REPO_ROOT, "research", "results", "phase2")
sweep(sfx) = JSON.parsefile(joinpath(R, "p2_3_4_h_sweep" * sfx, "summary.json"))["per_h"]
fits(sfx) = JSON.parsefile(joinpath(R, "p2_3_9_fit_h" * sfx, "summary.json"))["fits"]

ss, sl = sweep(""), sweep("_linear")
fs, fl = fits(""), fits("_linear")
rows = NamedTuple[]
for (a, b) in zip(fs, fl)
    @assert a["theta"] == b["theta"]
    hs, hl = a["h_final"], b["h_final"]
    ρs, ρl = RATE_REF * hs / 1, RATE_REF * hl / 2
    push!(rows, (theta=a["theta"], h_star_step=hs, h_star_linear=hl, h_ratio=hl / hs,
        rho_star_step=ρs, rho_star_linear=ρl, rho_ratio=ρl / ρs,
        max_abs_residual_step=maximum(abs.(a["residuals_final"])), max_abs_residual_linear=maximum(abs.(b["residuals_final"]))))
    @printf("θ = %.4f: h* step %.5f, linear %.5f (ratio %.4f); ρ* ratio %.4f; max |res| %.1f / %.1f\n",
        a["theta"], hs, hl, hl / hs, ρl / ρs, rows[end].max_abs_residual_step, rows[end].max_abs_residual_linear)
end
write_csv(joinpath(OUT, "comparison.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rate" => RATE_REF, "rows" => rows))

using CairoMakie
fig = Figure(size=(1000, 360))
ax1 = Axis(fig[1, 1], xlabel="h (µm)", ylabel="1 − S(50)", title="Attenuation vs h at ΔR₂(0) = 0.1 /ms")
scatterlines!(ax1, [p["h_um"] for p in ss], [p["attenuation_50"] for p in ss], label="step")
scatterlines!(ax1, [p["h_um"] for p in sl], [p["attenuation_50"] for p in sl], label="linear")
axislegend(ax1, position=:rb)
ax2 = Axis(fig[1, 2], xlabel="θ_relax", ylabel="ratio linear / step", xscale=log10, title="V4: h* and ρ* ratios")
scatterlines!(ax2, [r.theta for r in rows], [r.h_ratio for r in rows], label="h* ratio")
scatterlines!(ax2, [r.theta for r in rows], [r.rho_ratio for r in rows], label="ρ* ratio")
axislegend(ax2, position=:rt)
save(joinpath(OUT, "p2_4_8_compare.png"), fig)
println("saved to $OUT")
```

- [ ] **Step 2: Run it**

Run: `julia --project=research/baseline research/phase2/p2_4_8_compare_step_linear.jl`
Expected: 12 `θ = …` lines and a figure. It takes seconds, since it only reads saved results.

- [ ] **Step 3: Record, tick the tracker and commit**

Add a section "P2.4.8 (step vs linear): profile comparison at fixed ΔR₂(0)" to `research/phase2_record.md`. It should contain the comparison table (θ, h\* step, h\* linear, h\* ratio, ρ\* ratio, max |residual| for each) and a numbers-only interpretation: how the h\* ratio and ρ\* ratio behave across θ, and whether either profile fits the reference curves better.

In `research/project-progress-tracker.md`:
- **P2.4.1:** ◐, "step and linear in `layers.jl`; polynomial and exponential to come"
- **P2.4.2:** ◐, "step, linear: ∫g = 1 to 1e-10 (`test_layer.jl`)"
- **P2.4.3:** ☑, with the output folders of Tasks 3–5
- **P2.4.7:** ◐, "linear: h\*(θ) for all θ (`p2_3_9_fit_h_linear/`)"
- **P2.4.8:** ◐, "step vs linear (`p2_4_8_compare_step_linear/`)"
- Increase the Phase 2 "Done" count by 1 (P2.4.3).

```bash
git add research/phase2/p2_4_8_compare_step_linear.jl research/results/phase2/p2_4_8_compare_step_linear research/phase2_record.md research/project-progress-tracker.md
git commit -m "Compare step and linear profiles at fixed surface rate; update tracker

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
