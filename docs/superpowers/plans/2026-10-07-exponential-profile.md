# Exponential Near-Surface Profile on Walls Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a truncated, renormalised exponential profile g(u) ∝ e^(−u) to the wall layer and verify it like the step and linear profiles. Add an exact statistical check that holds for every profile: the mean of −ln M equals 2ρt/w. Also check the timestep constraint and compare all three profiles at fixed ΔR₂(0). (No V4 fit for the exponential: user decision, 2026-10-07.) This prepares the fixed-ρ sensitivity study.

**Architecture:**
- `LayerSide` gains `cutoff`, the support of the profile in µm measured from the wall. For step and linear it equals h, so their code paths and results are unchanged.
- For `:exponential`, g(u) = e^(−u)/(1 − e^(−c)) on 0 ≤ u ≤ c = cutoff/h, and 0 beyond. Renormalising keeps ∫g = 1, so ρ stays exactly the integrated relaxivity, which the fixed-ρ study needs.
- On a straight piece d changes linearly, so the mean of e^(−d/h) over the part inside the support has a closed form, h(e^(−d_a/h) − e^(−d_b/h))/(d_b − d_a). It is computed with `expm1` for accuracy.
- The `Walls` field `layer_cutoff` defaults to the spacing between walls, so an exponential layer reaches the opposite wall and the two faces' layers overlap. Non-repeating walls must set `layer_cutoff`.
- The research scripts take `SHAPE=exponential`.

**Tech Stack:** Julia 1.12, MCMRSimulator (this repo), `research/phase2/`, `Test` stdlib.

**Spec:**
- Proposal §3.2 III: ΔR₂(d) = (ρ/h)·g(d/h), ∫g = 1; exponential profile.
- `research/project-progress-tracker.md`: P2.4.1, P2.4.2, P2.4.5 (exponential end-to-end, truncation recorded), P2.4.8, P2.6.4. P2.4.7 (V4) is skipped for the exponential by user decision (2026-10-07).
- User decisions (2026-10-07):
  - exponential before the sensitivity study;
  - the fixed-ρ comparison is primary;
  - overlapping layers matter for decaying profiles;
  - the step profile is excluded where h > w/2.
- User framing: the baseline is a reference, not the truth; no fast-diffusion arguments.

**Depends on:** the linear plan (done) and P2.6.4 (done).

## Global Constraints

- Work only on branch `cc/near-surface-r2`.
- **Step and linear must stay bit-identical.** `cutoff == h` for both, and every existing package and research test must pass unchanged.
- **Units:** ρ in µm/ms; h and cutoff in µm; ΔR₂ in 1/ms.
- **Surface rate of the exponential:** ΔR₂(0) = (ρ/h)·g(0) = (ρ/h)/(1 − e^(−cutoff/h)). So at fixed ΔR₂(0), ρ = ΔR₂(0)·h·(1 − e^(−cutoff/h)).
- **Reference configuration:**
  - walls every 2 µm, D = 3 µm²/ms, no RF;
  - readouts 0:5:50 ms, τ = 1e-2 ms;
  - seeds 1–10, 10⁵ spins;
  - ΔR₂(0) = 0.1 ms⁻¹, layers on both faces;
  - exponential cutoff = 2 µm (the default).
- The uncertainty of a mean is its SEM. Characterisation results are reported, never pass/fail.
- **Package tests:** `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer", "cli"])'`. **Research tests:** `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`.
- **Long runs:** start them detached (`nohup … & disown`) and wait on them with a background `until` loop. A background command is killed after 2 h.
- At the end, the final review goes to a fresh subagent reviewer (user's standing rule).
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## The mean-exponent identity (exact for every profile)

Between reflecting walls, the uniform spin density is stationary under MCMR's step-and-mirror walk. Every point along a straight piece, x₀ + s·Δx folded by reflection, is therefore uniform too. Each spin's layer exponent is ∫ΔR₂ dt along its pieces, so its expectation is t·⟨ΔR₂⟩_volume.

If each layer's support lies within the gap (true for every configuration here), then ⟨ΔR₂⟩ = 2·(1/w)∫(ρ/h)g(d/h) dd = 2ρ/w, for any profile, any h and any τ. So:

  **E[−ln M⊥(t)] = 2ρt/w** (R₂_bulk = 0).

This holds whatever the regime or timestep, and it tests the normalisation, the truncation, the overlap sum and the image loop all at once. It is a statistical check: |z| ≤ 3 with z = (mean − 2ρt/w)/SEM over seeds. It is not a statement about the signal S = ⟨M⊥⟩. That depends on the higher moments, which is where profiles can differ.

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `src/geometries/internal/layers.jl` | Modify | `LayerSide.cutoff`, `SHAPES`, exponential `profile_g`, `surface_rate`, `side_exponent`, reach by cutoff, module docstring |
| `src/geometries/user/obstructions/obstructions.jl` | Modify | `layer_cutoff` field; shape description lists "exponential" |
| `src/geometries/user/fix.jl` | Modify | shape list from `Layers.SHAPES`, cutoff resolution and validation |
| `test/test_layer.jl`, `test/test_cli.jl` | Modify | package tests |
| `research/phase2/harness.jl` | Modify | `profile_g0(shape, h)`, `rho_for` for exponential |
| `research/phase2/toy.jl`, `research/phase2/p2_2_9_layer_trace_lib.jl` | Modify | exponential in the toy and in the brute-force profile average |
| `research/phase2/test/test_harness.jl`, `test_trace.jl`, `test_toy.jl` | Modify | research tests |
| `research/phase2/p2_4_5_mean_exponent.jl` | Create | mean-exponent identity, all three profiles, full size |
| `research/phase2/p2_6_4_profile.jl` | Modify | add the exponential scans |
| `research/phase2/p2_4_8_compare_step_linear.jl` | Modify | three profiles, ρ\* via `rho_for` |
| `research/phase2_record.md`, `research/project-progress-tracker.md` | Modify | records, ticks |

## Review Focus

1. **Step and linear unchanged.** With `cutoff == h`, the `segment_overlap` bounds and the image reach are identical. Test: all existing `layer` tests pass unchanged (Task 1), and the full suite at Task 2.
2. **Exponential near-degenerate pieces.** A piece with d_b − d_a ≈ 0 (stationary spin, or a tiny piece) must use the series branch and not divide 0 by 0. Test: Task 1 (stationary, and a piece of length 1e-12).
3. **Reach of the wall images.** The exponential layer reaches farther than h. The image loop must use `cutoff`, or contributions from the opposite wall are lost. Test: Task 2 (`layer_exponent` at x = 0.5 includes the far wall's term).
4. **Cutoff validation.** cutoff ≤ 0, cutoff > spacing, and a non-repeating wall without `layer_cutoff` must each raise a clear error. Test: Task 2.
5. **Normalisation and overlap end to end.** The mean-exponent identity at full size for every profile, including overlapping linear (h = 1.5) and exponential layers. Test: Task 3.

---

### Task 1: Exponential profile in the layer maths (P2.4.1, P2.4.2)

**Files:**
- Modify: `src/geometries/internal/layers.jl`
- Modify: `test/test_layer.jl` (append)

**Interfaces:**
- Produces:
  - `Layers.SHAPES = (:step, :linear, :exponential)`
  - `Layers.LayerSide(rho, h, shape, cutoff)`. `LayerSide(rho, h)` and `LayerSide(rho, h, shape)` still work (cutoff = h); for `:exponential` the three-argument form raises an error.
  - `Layers.profile_g(shape, u, c=Inf)`: c is the support in units of h, used only by `:exponential`.
  - `Layers.surface_rate(side)` = (ρ/h)·g(0; cutoff/h)
  - `Layers.side_exponent` for all three shapes
  - the reach of `layer_exponent` is `max(cutoff)`

- [ ] **Step 1: Append the failing tests to `test/test_layer.jl`**

```julia
@testset "test_layer.jl: exponential profile maths" begin
    Layers = mr.Geometries.Internal.Layers

    @testset "values and normalisation" begin
        c = 4.0
        N = 1 / (1 - exp(-c))
        @test Layers.profile_g(:exponential, 0.0, c) ≈ N
        @test Layers.profile_g(:exponential, 1.0, c) ≈ N * exp(-1)
        @test Layers.profile_g(:exponential, c + 1e-9, c) == 0.0
        @test Layers.profile_g(:exponential, -0.1, c) == 0.0
        @test Layers.profile_g(:exponential, 2.0) ≈ exp(-2.0)              # c = Inf: untruncated, ∫ = 1
        n = 400_000
        @test sum(Layers.profile_g(:exponential, (k - 0.5) / n * c, c) for k in 1:n) * c / n ≈ 1.0 atol=1e-9
        @test Layers.profile_g(:step, 0.5, 3.0) == 1.0                     # c ignored for step and linear
        @test Layers.profile_g(:linear, 0.5, 3.0) == 1.0
        @test Layers.SHAPES == (:step, :linear, :exponential)
    end

    @testset "constructor, cutoff and surface rate" begin
        @test_throws ErrorException Layers.LayerSide(0.1, 0.2, :exponential)   # needs a cutoff
        s = Layers.LayerSide(0.1, 0.2, :exponential, 2.0)
        @test s.cutoff == 2.0
        @test Layers.surface_rate(s) ≈ 0.1 / 0.2 / (1 - exp(-10))
        @test Layers.LayerSide(0.05, 0.5).cutoff == 0.5
        @test Layers.LayerSide(0.05, 0.5, :linear).cutoff == 0.5
        @test Layers.max_layer_rate([Layers.WallLayer(0.0, s, Layers.LayerSide(0.0, 0.0))]) ≈ 0.5 / (1 - exp(-10))
    end

    @testset "pieces vs brute force" begin
        Random.seed!(9)
        for _ in 1:200
            h = 0.05 + 0.5 * rand()
            cut = h + (2.0 - h) * rand()
            d0, d1 = 2.4 * rand() - 0.2, 2.4 * rand() - 0.2
            s = Layers.LayerSide(0.3 * h, h, :exponential, cut)
            n = 200_000
            brute = sum(Layers.profile_g(:exponential, (d0 + (d1 - d0) * (k - 0.5) / n) / h, cut / h) for k in 1:n) / n * (s.rho / h)
            @test Layers.side_exponent(s, d0, d1, 1.0) ≈ brute rtol=1e-4 atol=5e-6
        end
        s = Layers.LayerSide(0.02, 0.2, :exponential, 2.0)                 # ρ/h = 0.1
        r0 = Layers.surface_rate(s)
        @test Layers.side_exponent(s, 0.3, 0.3, 0.5) ≈ 0.5 * r0 * exp(-1.5)             # stationary
        @test Layers.side_exponent(s, 0.3, 0.3 + 1e-12, 0.5) ≈ 0.5 * r0 * exp(-1.5)     # tiny piece: series branch
        @test Layers.side_exponent(s, 0.0, 2.0, 1.0) ≈ 0.02 / 2.0 rtol=1e-12             # whole support: ρ/cutoff
        @test Layers.side_exponent(s, 2.5, 2.6, 1.0) == 0.0                              # beyond the cutoff
        @test Layers.side_exponent(s, -0.5, -0.1, 1.0) == 0.0                            # other side of the wall
    end
end
```

The "whole support" case: crossing the full support uniformly in time averages g over [0, c], which is 1/c. So the exponent is (ρ/h)·(1/c)·dt = ρ/cutoff = 0.02/2.0 = 0.01.

The brute-force tolerance `atol=5e-6` covers midpoint sampling across the jumps of g at d = 0 and d = cutoff (jump/n ≈ 1.5e-6).

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`
Expected: FAIL. `profile_g` has no 3-argument method; `LayerSide` has no field `cutoff`; `SHAPES` is not defined.

- [ ] **Step 3: Update `src/geometries/internal/layers.jl`**

Module docstring, replace the "Implemented:" line with:

```julia
Implemented: planar walls; profiles step g(u) = 1 (0 ≤ u ≤ 1), linear g(u) = 2(1 − u) (0 ≤ u ≤ 1) and
exponential g(u) = e^(−u)/(1 − e^(−c)) (0 ≤ u ≤ c = cutoff/h, truncated and renormalised); no finite RF pulses.
```

Replace the `LayerSide` docstring, struct and constructor with:

```julia
"Profiles implemented in `profile_g`."
const SHAPES = (:step, :linear, :exponential)

"""
    LayerSide(rho, h[, shape[, cutoff]])

Near-surface layer on one side of a surface: integrated relaxivity `rho` (um/ms), length scale `h` (um), profile
`shape` (one of `SHAPES`; default `:step`) and `cutoff` (um), the distance from the surface beyond which g = 0.
For `:step` and `:linear` the cutoff is `h`. `:exponential` needs an explicit cutoff. A side with `rho == 0` has no layer.
"""
struct LayerSide
    rho :: Float64
    h :: Float64
    shape :: Symbol
    cutoff :: Float64
end
LayerSide(rho, h) = LayerSide(rho, h, :step)
function LayerSide(rho, h, shape::Symbol)
    shape == :exponential && error("The exponential profile needs a cutoff: LayerSide(rho, h, :exponential, cutoff).")
    return LayerSide(rho, h, shape, h)
end
```

Replace `profile_g` with:

```julia
"""
    profile_g(shape, u, c=Inf)

Normalised profile g(u), ∫ g(u) du = 1 over its support (proposal §3.2 III). Zero outside the support.
`:step`: g = 1 on [0, 1]. `:linear`: g = 2(1 − u) on [0, 1]. `:exponential`: g = e^(−u)/(1 − e^(−c)) on [0, c],
with c = cutoff/h (c = Inf: untruncated). `c` is ignored by step and linear.
"""
function profile_g(shape::Symbol, u::Real, c::Real=Inf)
    shape == :exponential && return (0 <= u <= c) ? exp(-u) / -expm1(-c) : 0.0
    (0 <= u <= 1) || return 0.0
    shape == :step && return 1.0
    shape == :linear && return 2 * (1 - u)
    error("Unknown layer profile $shape")
end
```

Replace `surface_rate` with:

```julia
surface_rate(side::LayerSide) = iszero(side.rho) ? 0.0 : side.rho / side.h * profile_g(side.shape, 0.0, side.cutoff / side.h)
```

Replace `side_exponent` (docstring and body) with:

```julia
"""
    side_exponent(side, d0, d1, dt)

∫ΔR2(d(t)) dt over a straight segment of duration `dt` (ms), where `d0` and `d1` are the distances from the
surface at the start and end, positive into this side. d is linear along the segment, so the integral is exact
in closed form over the part of the segment inside the support 0 ≤ d ≤ cutoff:
- `:step`: ΔR2(0)·(time inside);
- `:linear`: g is linear along the segment, so its average over the part inside is its value at the midpoint;
- `:exponential`: the mean of e^(−d/h) for d uniform on [d_a, d_b] is e^(−d_a/h)·(1 − e^(−x))/x with
  x = (d_b − d_a)/h, evaluated with `expm1` (series 1 − x/2 for |x| < 1e-8).
"""
function side_exponent(side::LayerSide, d0::Float64, d1::Float64, dt::Float64)
    iszero(side.rho) && return 0.0
    (s_enter, s_exit) = segment_overlap(d0, d1, 0.0, side.cutoff)
    fraction = max(0.0, s_exit - s_enter)
    iszero(fraction) && return 0.0
    side.shape == :step && return surface_rate(side) * dt * fraction
    if side.shape == :exponential
        d_a = d0 + (d1 - d0) * s_enter
        x = (d1 - d0) * fraction / side.h
        mean_e = exp(-d_a / side.h) * (abs(x) < 1e-8 ? 1 - x / 2 : -expm1(-x) / x)
        return side.rho / side.h / -expm1(-side.cutoff / side.h) * mean_e * dt * fraction
    end
    d_mid = d0 + (d1 - d0) * (s_enter + s_exit) / 2
    return side.rho / side.h * profile_g(side.shape, d_mid / side.h) * dt * fraction
end
```

In `layer_exponent(layers::Vector{WallLayer}, …)`, replace

```julia
            reach = max(layer.positive.h, layer.negative.h)
```

with

```julia
            reach = max(layer.positive.cutoff, layer.negative.cutoff)
```

For step and linear, `cutoff == h`, so `segment_overlap(d0, d1, 0.0, side.cutoff)` and the reach are the same floating-point values as before. `surface_rate` gets `profile_g(:step, 0.0, 1.0) == 1.0` and `profile_g(:linear, 0.0, 1.0) == 2.0`, the same as before.

- [ ] **Step 4: Run the tests and check they pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`
Expected: PASS, with every earlier testset unchanged and "exponential profile maths" added.

- [ ] **Step 5: Commit**

```bash
git add src/geometries/internal/layers.jl test/test_layer.jl
git commit -m "Add the truncated exponential near-surface profile to the layer maths

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: `layer_shape = "exponential"` and `layer_cutoff` on `Walls` (P2.4.1, P2.4.5)

**Files:**
- Modify: `src/geometries/user/obstructions/obstructions.jl`, `src/geometries/user/fix.jl`
- Modify: `test/test_layer.jl`, `test/test_cli.jl`

**Interfaces:**
- Consumes: `Layers.SHAPES` and `LayerSide(rho, h, shape, cutoff)` (Task 1).
- Produces: the `Walls` field `layer_cutoff::Float64` (default `nothing`, meaning the spacing between walls; it applies to both sides and to every exponential side). `layer_shape` also accepts `"exponential"`.

- [ ] **Step 1: Append the failing tests**

To `test/test_layer.jl`:

```julia
@testset "test_layer.jl: exponential on Walls" begin
    Layers = mr.Geometries.Internal.Layers
    geom(walls) = mr.Simulation([]; geometry=walls, verbose=false).geometry
    v(x) = SVector{3, Float64}(x, 0.3, -1.2)

    @testset "resolution and validation" begin
        l = geom(mr.Walls(repeats=2., layer_rho=0.02, layer_h=0.2, layer_shape="exponential"))[1].layer[1]
        @test (l.positive.shape, l.positive.cutoff, l.negative.cutoff) == (:exponential, 2.0, 2.0)
        l = geom(mr.Walls(repeats=2., layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=1.0))[1].layer[1]
        @test l.positive.cutoff == 1.0
        l = geom(mr.Walls(repeats=2., layer_rho=0.02, layer_h=0.2, layer_shape_negative="exponential", layer_cutoff=1.0))[1].layer[1]
        @test (l.positive.shape, l.positive.cutoff, l.negative.shape, l.negative.cutoff) == (:step, 0.2, :exponential, 1.0)
        @test_throws ErrorException geom(mr.Walls(repeats=2., layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=2.5))
        @test_throws ErrorException geom(mr.Walls(repeats=2., layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=0.0))
        @test_throws ErrorException geom(mr.Walls(layer_rho=0.02, layer_h=0.2, layer_shape="exponential"))   # single wall: no spacing
        @test geom(mr.Walls(layer_shape="exponential"))[1].layer === nothing                                  # ρ = 0: inert, no error
        l = geom(mr.Walls(layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=3.0))[1].layer[1]
        @test l.positive.cutoff == 3.0
    end

    @testset "reach includes the opposite wall" begin
        g = geom(mr.Walls(repeats=2., layer_rho=0.02, layer_h=0.2, layer_shape="exponential"))   # ρ/h = 0.1, cutoff 2
        N = 1 / (1 - exp(-10))
        # stationary at x = 0.5: d = 0.5 from wall 0 (positive side), 1.5 from wall 2 (negative side)
        @test Layers.layer_exponent(g, v(0.5), v(0.5), 1.0) ≈ 0.1 * N * (exp(-2.5) + exp(-7.5)) rtol=1e-12
    end

    @testset "JSON round trip keeps layer_cutoff" begin
        io = IOBuffer()
        mr.write_geometry(io, mr.Walls(repeats=2., layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=1.5))
        back = mr.read_geometry_json(String(take!(io)))
        @test back.layer_shape.value == "exponential"
        @test back.layer_cutoff.value == 1.5
    end
end
```

To `test/test_cli.jl`, after the `"create walls with a linear layer"` testset:

```julia
        @testset "create walls with an exponential layer" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create walls 1 test.json --repeats 2 --layer_rho 0.02 --layer_h 0.2 --layer_shape exponential --layer_cutoff 1.5")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["layer_shape"] == "exponential"
                @test result["layer_cutoff"] == 1.5
            end
        end
```

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer", "cli"])'`
Expected: FAIL, with `layer_shape must be "step" or "linear"` for "exponential", `KeyError: layer_cutoff`, and the CLI reporting an unrecognized option.

- [ ] **Step 3: Fields**

In `src/geometries/user/obstructions/obstructions.jl`, change the `layer_shape` description to
`"Near-surface layer: profile g on both sides of the wall, unless a side-specific value is set. One of step, linear, exponential."`
(no double quotes inside). After the `layer_shape_negative` field, add:

```julia
        Field{Float64}(:layer_cutoff, "Near-surface layer: support (um) of the exponential profile, measured from the wall; g is truncated there and renormalised so that its integral is 1. Default: the spacing between walls. Ignored by step and linear."),
```

- [ ] **Step 4: Resolution and validation in `src/geometries/user/fix.jl` (`wall_layers`)**

In `resolve_shape`, replace the check line with:

```julia
        Symbol(shape) in Internal.Layers.SHAPES || error("Wall $i, $side side: layer_shape must be one of $(join(string.(Internal.Layers.SHAPES), ", ")), got \"$shape\".")
```

Replace

```julia
            Internal.LayerSide(rho, h, resolve_shape(side, i))
```

with

```julia
            shape = resolve_shape(side, i)
            shape == :exponential || return Internal.LayerSide(rho, h, shape)
            iszero(rho) && return Internal.LayerSide(rho, h, shape, h)          # inert side: no cutoff needed
            given = walls.layer_cutoff[i]
            cutoff = isnothing(given) ? spacing : Float64(given)
            isfinite(cutoff) || error("Wall $i, $side side: an exponential layer on a non-repeating wall needs `layer_cutoff`.")
            cutoff > 0 || error("Wall $i, $side side: layer_cutoff must be > 0, got $cutoff.")
            cutoff <= spacing || error("Wall $i, $side side: layer_cutoff = $cutoff um exceeds the spacing between walls ($spacing um).")
            Internal.LayerSide(rho, h, shape, cutoff)
```

`return` inside the `do` block returns from that block, i.e. it gives this side's value. A side with ρ = 0 is inert, so it skips the cutoff checks; otherwise `Walls(layer_shape="exponential")` on a single wall would fail with no layer at all. `spacing` is `Inf` for non-repeating walls, so an explicit cutoff passes the `≤ spacing` check there.

- [ ] **Step 5: Run the tests and check they pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer", "cli", "collisions", "transfer"])'`
Expected: PASS.

- [ ] **Step 6: Commit, then run the full suite**

```bash
git add src/geometries/user/obstructions/obstructions.jl src/geometries/user/fix.jl test/test_layer.jl test/test_cli.jl
git commit -m "Add exponential layer_shape and layer_cutoff to Walls

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["no-plots"])'`
Expected: `Testing MCMRSimulator tests passed`.

---

### Task 3: Research support and the mean-exponent identity (P2.4.5)

**Files:**
- Modify: `research/phase2/harness.jl`, `research/phase2/toy.jl`, `research/phase2/p2_2_9_layer_trace_lib.jl`
- Modify: `research/phase2/test/test_harness.jl`, `test_trace.jl`, `test_toy.jl`
- Create: `research/phase2/p2_4_5_mean_exponent.jl`
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `Walls(…; layer_shape="exponential")` (Task 2) and `Layers.profile_g(shape, u, c)` (Task 1).
- Produces:
  - `profile_g0(shape, h)` = g(0) with cutoff W_REF
  - `rho_for(rate, h; shape)` for all three shapes
  - `toy_layer_signal(…; shape="exponential")` (cutoff = w)
  - `sampled_profile_average(…; shape=:exponential)` (cutoff W_REF)
  - `results/phase2/p2_4_5_mean_exponent/summary.json`
  - `results/phase2/p2_2_9_layer_trace_exponential/`

- [ ] **Step 1: Append the failing research tests**

To `research/phase2/test/test_harness.jl`:

```julia
@testset "harness: exponential" begin
    @test profile_g0("step", 0.3) == 1.0 && profile_g0("linear", 0.3) == 2.0
    @test profile_g0("exponential", 0.5) ≈ 1 / (1 - exp(-4))
    @test rho_for(0.1, 0.5; shape="exponential") ≈ 0.05 * (1 - exp(-4))
    @test rho_for(0.1, 0.0; shape="exponential") == 0.0
    l = Simulation([]; geometry=layer_walls(h=0.5, rate=0.1, shape="exponential"), verbose=false).geometry[1].layer[1]
    @test (l.positive.shape, l.positive.cutoff) == (:exponential, 2.0)
    @test mr.Geometries.Internal.Layers.surface_rate(l.positive) ≈ 0.1
end
```

To `research/phase2/test/test_trace.jl`:

```julia
@testset "trace: exponential profile average" begin
    h = 0.4
    N = 1 / (1 - exp(-W_REF / h))
    @test sampled_profile_average(0.5, 0.5, h; shape=:exponential) ≈ N * (exp(-0.5 / h) + exp(-1.5 / h))
    # straight piece 0.2 → 0.7: closed-form mean of the two faces' terms
    mean_face(a, b) = h * (exp(-a / h) - exp(-b / h)) / (b - a)
    @test sampled_profile_average(0.2, 0.7, h; shape=:exponential) ≈ N * (mean_face(0.2, 0.7) + mean_face(1.3, 1.8)) rtol=1e-6
end
```

To `research/phase2/test/test_toy.jl`:

```julia
@testset "toy: exponential profile and the mean-exponent identity" begin
    # E[−ln M] = 2ρT/w for any profile (uniform density is stationary); exponential, overlapping both faces
    ρ = rho_for(0.5, 0.3; shape="exponential")
    s = toy_layer_signal(; h=0.3, rate=0.5, tau=0.01, T=2.0, nspins=20_000, seed=1, method=:exact, shape="exponential")
    m = mean(-log.(s)); se = std(-log.(s)) / sqrt(length(s))
    @test abs(m - 2ρ * 2.0 / W_REF) <= 4se
    # linear overlapping (h = 1.5 > w/2)
    ρl = rho_for(0.5, 1.5; shape="linear")
    s = toy_layer_signal(; h=1.5, rate=0.5, tau=0.01, T=2.0, nspins=20_000, seed=2, method=:exact, shape="linear")
    m = mean(-log.(s)); se = std(-log.(s)) / sqrt(length(s))
    @test abs(m - 2ρl * 2.0 / W_REF) <= 4se
end
```

- [ ] **Step 2: Run the research tests and check they fail**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: FAIL, with `profile_g0` having no 2-argument method, and the toy rejecting `shape="exponential"`.

- [ ] **Step 3: Implement**

In `research/phase2/harness.jl`, replace `profile_g0` and `rho_for` with:

```julia
"g(0) of a profile for thickness `h` (µm); the exponential is truncated at W_REF (the default `layer_cutoff`)."
profile_g0(shape, h) = mr.Geometries.Internal.Layers.profile_g(Symbol(shape), 0.0, W_REF / h)

"ρ (µm/ms) giving surface excess rate ΔR2(0) = `rate` (1/ms) for thickness `h` (µm): ρ = rate·h / g(0)."
rho_for(rate, h; shape="step") = iszero(h) ? 0.0 : Float64(rate * h / profile_g0(shape, h))
```

For step and linear, `profile_g0` returns exactly 1.0 and 2.0, so `rho_for` is unchanged.

In `research/phase2/toy.jl`, replace the `shape in …` check and the `side = …` line with:

```julia
    shape in ("step", "linear", "exponential") || error("shape must be step, linear or exponential, got $shape")
    side = _Layers.LayerSide(rho_for(rate, h; shape=shape), h, Symbol(shape), shape == "exponential" ? w : h)
```

Replace the linear endpoint line:

```julia
                expo += tau * side.rho / h * (_Layers.profile_g(:linear, x / h) + _Layers.profile_g(:linear, (w - x) / h))
```

with

```julia
                expo += tau * side.rho / h * (_Layers.profile_g(side.shape, x / h, side.cutoff / h) +
                                              _Layers.profile_g(side.shape, (w - x) / h, side.cutoff / h))
```

The `else` branch of `:exact` already calls `side_exponent`, so it covers the exponential. Update the docstring to say "linear or exponential (cutoff = w) use MCMR's `side_exponent`". `rho_for` assumes cutoff = W_REF, so add a line to the docstring: "the toy's w must equal W_REF for the exponential".

In `research/phase2/p2_2_9_layer_trace_lib.jl`, replace the `g(u) = …` line of `sampled_profile_average` with:

```julia
    cut = W_REF / h
    g(u) = shape == :exponential ? ((0 <= u <= cut) ? exp(-u) / -expm1(-cut) : 0.0) :
           (0 <= u <= 1 ? (shape == :step ? 1.0 : 2 * (1 - u)) : 0.0)
```

Add to its docstring: "The exponential is truncated at W_REF (support of each face = the gap)."

- [ ] **Step 4: Run the research tests and check they pass**

Run: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
Expected: PASS.

- [ ] **Step 5: Exponential single-spin trace**

Run: `SHAPE=exponential julia --project=research/baseline -t 1 research/phase2/p2_2_9_layer_trace.jl`
Expected: `→ PASS` (max error ≤ 1e-6) and a positive count of reflections inside the layer. The output goes to `results/phase2/p2_2_9_layer_trace_exponential/`.

- [ ] **Step 6: Create the mean-exponent script**

`research/phase2/p2_4_5_mean_exponent.jl`:

```julia
# P2.4.5 — mean-exponent identity for every profile, full size: E[−ln M⊥(t)] = 2ρt/w (R2_bulk = 0).
# Holds exactly in expectation for any profile, h and τ whenever each layer's support lies inside the gap:
# the uniform density is stationary under the reflected walk, also at every point along a straight piece.
# Pass: |z| ≤ 3 at t = 50 ms for every configuration, z = (mean over seeds − 2ρt/w) / SEM over seeds.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_4_5_mean_exponent.jl

include(joinpath(@__DIR__, "harness.jl"))

const OUT = phase2_outdir("p2_4_5_mean_exponent")
const NSPINS = parse(Int, get(ENV, "NSPINS", string(NSPINS_REF)))
const T = 50.0
const CONFIGS = [("step", 0.2), ("linear", 0.2), ("linear", 1.5), ("exponential", 0.2), ("exponential", 0.5)]

rows = NamedTuple[]
for (shape, h) in CONFIGS, τ in (1e-2, 4e-2)
    ρ = rho_for(RATE_REF, h; shape=shape)
    r = run_walls(layer_walls(h=h, rate=RATE_REF, shape=shape); timestep=τ, times=[0.0, T], nspins=NSPINS, snapshots=true)
    per_seed = [mean(-log(s.orientations[1].transverse) for s in snap.spins) for snap in r.snapshots]
    m, se = mean(per_seed), std(per_seed) / sqrt(length(per_seed))
    expected = 2ρ * T / W_REF
    z = (m - expected) / se
    push!(rows, (shape=shape, h_um=h, tau_ms=τ, rho=ρ, expected=expected, mean=m, sem=se, z=z, pass=abs(z) <= 3))
    @printf("%-11s h = %.2f, τ = %.0e: mean −ln M = %.6f ± %.1e, 2ρt/w = %.6f, z = %+.2f → %s\n",
        shape, h, τ, m, se, expected, z, abs(z) <= 3 ? "PASS" : "FAIL")
end
write_csv(joinpath(OUT, "mean_exponent.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "nspins" => NSPINS, "T_ms" => T,
    "rate" => RATE_REF, "rows" => rows, "pass" => all(r.pass for r in rows)))
println("saved to $OUT")
```

- [ ] **Step 7: Run it**

Run: `julia --project=research/baseline -t 8 research/phase2/p2_4_5_mean_exponent.jl`
Expected: 10 lines, all `→ PASS`. About 15 min. A FAIL means a normalisation, truncation, overlap or image error. Use superpowers:systematic-debugging, starting with the `layer_exponent` reach.

- [ ] **Step 8: Record and commit**

Add summary rows and sections to `research/phase2_record.md`:
- "P2.4.5: exponential profile — single-spin trace" (same layout as the linear trace section);
- "P2.4.5: mean-exponent identity (all profiles)", containing the derivation paragraph from this plan, a table (shape, h, τ, 2ρt/w, mean ± SEM, z, outcome), and one sentence on what the identity does and does not test.

```bash
git add research/phase2/harness.jl research/phase2/toy.jl research/phase2/p2_2_9_layer_trace_lib.jl research/phase2/test/test_harness.jl research/phase2/test/test_trace.jl research/phase2/test/test_toy.jl research/phase2/p2_4_5_mean_exponent.jl research/results/phase2/p2_2_9_layer_trace_exponential research/results/phase2/p2_4_5_mean_exponent research/phase2_record.md
git commit -m "Add exponential support to the research code; mean-exponent identity for all profiles

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Exponential h sweep (P2.4.5)

**Files:**
- Modify: `research/phase2_record.md`

**Interfaces:**
- Consumes: `p2_3_4_h_sweep.jl` with `SHAPE` (it uses `H_GRID_LINEAR` for any non-step shape, up to h = 2.0 = w, which is a valid decay length with cutoff w).
- Produces: `results/phase2/p2_3_4_h_sweep_exponential/`.

- [ ] **Step 1: Run the sweep (detached)**

```bash
nohup env SHAPE=exponential julia --project=research/baseline -t 8 research/phase2/p2_3_4_h_sweep.jl > /tmp/exp_sweep.log 2>&1 & disown
```

Wait with a background `until grep -q "saved to" /tmp/exp_sweep.log || ! pgrep -f p2_3_4_h_sweep; do sleep 60; done`.
Expected: 16 `h = …` lines and `15 of 15 … → PASS` for monotonicity. Record the range of 1 − S(50) (characterisation only; no V4 fit follows). About 30 min.

- [ ] **Step 2: Record and commit**

Add a summary row and a section "P2.4.5: exponential profile — h sweep" to `research/phase2_record.md`, with the same columns as the linear sweep table (h, S(25), S(50) ± SEM, 1 − S(50), curvature ± SEM). Add one more column, ρ = ΔR₂(0)·h·(1 − e^(−w/h)). The interpretation must be numbers only.

```bash
git add research/results/phase2/p2_3_4_h_sweep_exponential research/phase2_record.md
git commit -m "Run the h sweep for the exponential profile

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Timestep constraint for the exponential profile (P2.6.4)

**Files:**
- Modify: `research/phase2/p2_6_4_profile.jl`, `research/phase2_record.md`, `research/project-progress-tracker.md`

**Interfaces:**
- Consumes: `p2_6_scan.jl` with `SHAPE=exponential`; `p2_6_1_mcmr_check.jl` with `SCAN`.
- Produces: `results/phase2/p2_6_scan_E_*`, `p2_6_1_mcmr_check_E_rate1/`, and an updated `p2_6_4_profile/`.

- [ ] **Step 1: Run the scans and the MCMR check (detached chain)**

The settings are the same as the linear set. The scan's `RATE` is ΔR₂(0) of the exponential (its `rho_for` takes care of g(0)):

```bash
cat > /tmp/exp_ts.sh <<'EOF'
#!/bin/zsh
set -e
cd /Users/ishimiyabisaki/Documents/mcmrsimulator.jl-main-cc-near-surface-r2
J="julia --project=research/baseline -t 8"
scan() { env SHAPE=exponential H=$1 RATE=$2 D=$3 TAU_LO=$4 TAU_HI=1e-1 NSPINS=100000 T=10 LABEL=$5 ${=J} research/phase2/p2_6_scan.jl > /tmp/$5.log 2>&1; echo "$5 done"; }
scan 0.1 1.0 3 1e-4 E_rate1
SCAN=E_rate1 TAUS=1e-1,1e-2,1e-3 ${=J} research/phase2/p2_6_1_mcmr_check.jl > /tmp/E_mcmr.log 2>&1; echo "mcmr done"
scan 0.1 0.1 3 1e-5 E_h0.1
scan 0.1 0.5 3 1e-4 E_rate0.5
scan 0.1 2.0 3 1e-4 E_rate2
scan 0.4 1.0 3 1e-4 E_h0.4_rate1
scan 0.1 1.0 1 1e-4 E_D1_rate1
scan 0.01 10.0 3 1e-5 E_h0.01_rate10
echo CHAIN DONE
EOF
chmod +x /tmp/exp_ts.sh; nohup /tmp/exp_ts.sh > /tmp/exp_ts.out 2>&1 & disown
```

Wait for `CHAIN DONE`, re-arming the watcher at 2 h. Expected: about 2 h. The MCMR check should print `→ PASS`. A FAIL means the toy doesn't represent MCMR for this profile, so stop and investigate.

- [ ] **Step 2: Extend the comparison**

In `research/phase2/p2_6_4_profile.jl`, replace `PAIRS` with triples, and loop over all members:

```julia
const GROUPS = [("V1_h0.1", "L_h0.1", "E_h0.1"), ("V2_rate0.5", "L_rate0.5", "E_rate0.5"), ("V2_rate1", "L_rate1", "E_rate1"),
    ("V2_rate2", "L_rate2", "E_rate2"), ("V2_h0.4_rate1", "L_h0.4_rate1", "E_h0.4_rate1"), ("V3_D1_rate1", "L_D1_rate1", "E_D1_rate1"),
    ("V2_h0.01_rate10", "L_h0.01_rate10", "E_h0.01_rate10"), ("", "L_h0.2_rate1", "")]
```

Replace every `for (a, b) in PAIRS` with `for grp in GROUPS`. In the table loop, use `for d in describe.(grp)`. In the figure loop, use `for (l, ls) in zip(grp, (:dash, :solid, :dot))`. Add `worst_exp` for `shape == "exponential"`, mirroring `worst_step`. Write it to the summary and print it. Update the title to "(step: dashed, linear: solid, exponential: dotted)".

Run: `julia --project=research/baseline research/phase2/p2_6_4_profile.jl`
Expected: 22 lines and the three maxima at the default constraint, each printed with ± σ.

- [ ] **Step 3: Record, update the tracker, commit**

- **Record:** add the exponential rows and the MCMR check line to the P2.6.4 section of `research/phase2_record.md`, and extend the interpretation (numbers with σ).
- **Docstring:** if the exponential maximum at c = 0.005 exceeds the current docstring statement, update the `DEFAULT_LAYER_SCALING` docstring in `src/timesteps.jl`. Its wording must stay true for all three profiles.
- **Tracker:** update the P2.6.4 note to include the exponential.

```bash
git add research/phase2/p2_6_4_profile.jl research/results/phase2/p2_6_scan_E_* research/results/phase2/p2_6_1_mcmr_check_E_rate1 research/results/phase2/p2_6_4_profile research/phase2_record.md research/project-progress-tracker.md src/timesteps.jl
git commit -m "P2.6.4: timestep error of the exponential profile

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Three-profile comparison at fixed ΔR₂(0), tracker (P2.4.8)

**Files:**
- Modify: `research/phase2/p2_4_8_compare_step_linear.jl`, `research/phase2_record.md`, `research/project-progress-tracker.md`

**Interfaces:**
- Consumes: the sweeps of all three profiles, and the V4 fits of step and linear (no exponential V4).
- Produces: `results/phase2/p2_4_8_compare_step_linear/` with the exponential added to the attenuation comparison (the folder name is kept for continuity).

- [ ] **Step 1: Extend the script**

- Load `sx = sweep("_exponential")`.
- The V4 table (step vs linear, per θ) is unchanged.
- Add a second table, written to `attenuation.csv`: for every h present in all three sweeps, h, S(50) ± SEM for step, linear and exponential, and ρ for each, `rho_for(RATE_REF, h; shape)`. Step rows with h > w/2 are not present (the step grid ends at 1.0 = w/2). Print one line per h.
- Add the exponential to the attenuation panel (left) of the figure. The V4 ratio panel (right) stays step vs linear.

Run: `julia --project=research/baseline research/phase2/p2_4_8_compare_step_linear.jl`
Expected: the 12 unchanged θ lines, one line per common h, and a figure.

- [ ] **Step 2: Record, tick the tracker, commit**

- **Record:** extend the P2.4.8 section with the attenuation table (three profiles at fixed ΔR₂(0)) and a numbers-only interpretation with SEM.
- **Tracker:**
  - **P2.4.1:** ◐, "step, linear, exponential; polynomial to come"
  - **P2.4.2:** ◐, "step, linear, exponential (truncated, renormalised): ∫g = 1 to 1e-9"
  - **P2.4.5:** ☑, with the folders of Tasks 3–4, and the truncation stated as "cutoff = w, renormalised; no loss"
  - **P2.4.7:** ◐ unchanged, add "exponential skipped (user decision 2026-10-07)"
  - **P2.4.8:** ◐, "step, linear, exponential"
  - Increase the Phase 2 "Done" count by 1 (P2.4.5).

```bash
git add research/phase2/p2_4_8_compare_step_linear.jl research/results/phase2/p2_4_8_compare_step_linear research/phase2_record.md research/project-progress-tracker.md
git commit -m "Add the exponential profile to the fixed-rate comparison; update tracker

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
