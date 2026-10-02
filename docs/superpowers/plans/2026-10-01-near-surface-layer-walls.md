# Near-Surface R₂ Layer on Walls (Step Profile, Side-Specific) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a distance-dependent transverse relaxation layer, R₂_total(d) = R₂_bulk + ΔR₂(d) with ΔR₂(d) = (ρ/h)·g(d/h), to planar walls in MCMRSimulator. Each side of a wall gets its own (ρ, h), and the step profile g is used. There are no finite RF pulses.

**Architecture:** A new internal module `Layers` holds the per-side parameters and the exact segment integral. Along a straight free segment the distance to a wall changes linearly, so the time spent inside a step-profile layer has a closed form, `segment_overlap`. User-facing fields on `Walls` are resolved into a `Vector{WallLayer}` at `fix` time and stored on the `FixedObstructionGroup`. The group stores `nothing` when no wall has a layer, so the existing code path is untouched. In `draw_step!`, right after the existing `relax!` call for each free segment, `apply_layer!` multiplies the transverse magnetisation by exp(−∫ΔR₂ dt). The bulk R₂ factor from `relax!` and the layer factor multiply, which is exactly R₂_total = R₂_bulk + ΔR₂.

**Tech Stack:** Julia 1.12, MCMRSimulator (this repository), `Test` stdlib, StaticArrays.

**Spec:** Proposal (Y. Shi, 2026-09-14) §3.2 I–IV, Eqs. 5–10; `research/project-progress-tracker.md` P2.1.1, P2.1.2, P2.2.1–P2.2.8; side-specific design agreed 2026-10-01 (each side of a surface has its own ρ and h; "both sides equal" and "one side only" are settings).

## Global Constraints

- Work only on branch `cc/near-surface-r2`; never commit to other branches.
- Units follow MCMR: lengths in µm, times in ms, rates in 1/ms (MCMR calls this kHz). ρ in µm/ms, h in µm, ΔR₂(0) = ρ/h in 1/ms.
- With no layer set, results must be **bit-identical** to the unmodified simulator (`baseline-v1`): same positions, same magnetisation.
- Scope: `Walls` only; step profile only (g(u) = 1 for 0 ≤ u ≤ 1); no finite RF pulses (instantaneous pulses are fine).
- Side naming for walls: **positive** = the side where the wall's local coordinate is larger than the wall position (the +x side for the default `rotation=:x`); **negative** = the other side.
- Overlapping layers add (proposal Eq. 10).
- Run tests with: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

**Out of scope (later plans):** profiles other than step (P2.4), the timestep constraint τ ≤ c_h·h²/D (P2.6), cylinders/spheres (Phase 3), finite RF (Phase 3), the "rigid-lattice ceiling" warning in P2.1.2 (needs a physical value first).

## How the three designs are configured

| Design | User code |
|---|---|
| Same layer on both sides (used to fit to B4) | `Walls(repeats=2, layer_rho=ρ, layer_h=h)` |
| One side only | `Walls(repeats=2, layer_rho_positive=ρ, layer_h_positive=h)` |
| Different on each side | `Walls(repeats=2, layer_rho_positive=ρ₁, layer_h_positive=h₁, layer_rho_negative=ρ₂, layer_h_negative=h₂)` |

Side-specific fields override the both-sides fields; an unset side-specific field falls back to `layer_rho` / `layer_h` (default 0 = no layer).

## Map to the tracker

| Tracker | Where |
|---|---|
| P2.1.1 parameter interface | Task 2 |
| P2.1.2 validation | Task 2 |
| P2.2.1 distance evaluation | Task 3 (`x − c` on the positive side, `c − x` on the negative side) |
| P2.2.2 / P2.2.3 entry / exit crossings | Task 1 (`segment_overlap` returns the entry and exit points) |
| P2.2.4 crossing with reflection | Task 4 (reflection test) |
| P2.2.5 segmentation | Existing `draw_step!` already splits a step at every reflection. Within a straight piece, the step-profile integral is exact without further splitting (Task 1). |
| P2.2.6 local relaxation integral | Task 1 (`side_exponent`) |
| P2.2.7 transverse update | Task 4 (`apply_layer!`) |
| P2.2.8 λ = 0 short-circuit | Task 4 (bit-identical test) and Task 5 (runtime) |

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `src/geometries/internal/layers.jl` | Create | `LayerSide`, `WallLayer`, segment maths, `layer_exponent` for a group and a geometry |
| `src/geometries/internal/internal.jl` | Modify | include `layers.jl`, import its names |
| `src/geometries/internal/fixed_obstruction_groups.jl` | Modify | add `layer` field to `FixedObstructionGroup` |
| `src/geometries/user/obstructions/obstructions.jl` | Modify | six layer fields on `Walls` |
| `src/geometries/user/fix.jl` | Modify | `wall_layers` (resolve + validate), pass `layer` to the group |
| `src/evolve.jl` | Modify | `apply_layer!`, call it after `relax!` in `draw_step!` |
| `test/test_layer.jl` | Create | all tests for this feature |
| `test/runtests.jl` | Modify | register `"layer"` |

## Review Focus

1. **A segment longer than the wall spacing** (large τ or permeable walls): every wall image the segment passes must be counted. Test: Task 3, segment x = 0.1 → 4.1.
2. **A spin exactly on a wall or exactly at d = h**: the boundary counts as inside the layer (closed interval) and must not double count or vanish. Test: Task 1, stationary points at d = 0 and d = h.
3. **A finite RF pulse in a simulation with a layer** must stop with a clear error, not silently ignore the layer's interaction with the pulse. Test: Task 4.
4. **Surface sticking (`density > 0`) together with a layer**: the layer is not applied while a spin is stuck; the user must be warned. Test: Task 2.
5. **Per-wall values** (a vector of `layer_rho` over several walls) must reach the right wall. Test: Task 2 and Task 3.

---

### Task 1: Layer types and exact segment integral

**Files:**
- Create: `src/geometries/internal/layers.jl`
- Modify: `src/geometries/internal/internal.jl`
- Create: `test/test_layer.jl`
- Modify: `test/runtests.jl:15-30` (the `all_tests` list)

**Interfaces:**
- Produces:
  - `Layers.LayerSide(rho::Float64, h::Float64)`
  - `Layers.WallLayer(position::Float64, positive::LayerSide, negative::LayerSide)`
  - `Layers.surface_rate(side::LayerSide)::Float64`, which returns ΔR₂(0) = ρ/h
  - `Layers.segment_overlap(a0, a1, lo, hi)::Tuple{Float64, Float64}`, which returns `(s_enter, s_exit)`
  - `Layers.segment_fraction(a0, a1, lo, hi)::Float64`
  - `Layers.side_exponent(side::LayerSide, d0, d1, dt)::Float64`

- [ ] **Step 1: Register the test file and write the failing tests**

In `test/runtests.jl`, add `"layer",` to `all_tests` right after `"subsets",`:

```julia
all_tests = [
    "collisions",
    "evolve",
    "known_sequences",
    "meshes",
    "offresonance",
    "transfer",
    "permeability",
    "radio_frequency",
    "hierarchical_mri",
    "various",
    "subsets",
    "layer",
    "swc",
    "plots",
    "cli",
]
```

Create `test/test_layer.jl`:

```julia
@testset "test_layer.jl: segment maths" begin
    Layers = mr.Geometries.Internal.Layers

    @testset "segment_overlap: entry and exit points" begin
        # enters the layer [0, 0.5]: starts at 1.0 (outside), ends at 0.2 (inside)
        (s_in, s_out) = Layers.segment_overlap(1.0, 0.2, 0.0, 0.5)
        @test s_in ≈ 0.625
        @test s_out ≈ 1.0
        # exits the layer
        (s_in, s_out) = Layers.segment_overlap(0.2, 1.0, 0.0, 0.5)
        @test s_in ≈ 0.0
        @test s_out ≈ 0.375
    end

    @testset "segment_fraction" begin
        # neither endpoint inside, the segment passes through [0.25, 0.5]
        @test Layers.segment_fraction(-1.0, 2.0, 0.25, 0.5) ≈ 0.25 / 3
        # no overlap
        @test Layers.segment_fraction(0.6, 1.4, 0.0, 0.5) == 0.0
        # stationary inside / outside
        @test Layers.segment_fraction(0.3, 0.3, 0.0, 0.5) == 1.0
        @test Layers.segment_fraction(0.7, 0.7, 0.0, 0.5) == 0.0
        # stationary exactly on the wall (d = 0) and on the layer edge (d = h): inside
        @test Layers.segment_fraction(0.0, 0.0, 0.0, 0.5) == 1.0
        @test Layers.segment_fraction(0.5, 0.5, 0.0, 0.5) == 1.0
    end

    @testset "side_exponent: step profile" begin
        side = Layers.LayerSide(0.05, 0.5)          # ΔR2(0) = 0.05 / 0.5 = 0.1 /ms
        @test Layers.surface_rate(side) ≈ 0.1
        @test Layers.side_exponent(side, 0.1, 0.3, 0.2) ≈ 0.1 * 0.2         # fully inside
        @test Layers.side_exponent(side, 0.4, 0.6, 1.0) ≈ 0.1 * 1.0 * 0.5   # half inside
        @test Layers.side_exponent(Layers.LayerSide(0.0, 0.0), 0.1, 0.3, 0.2) == 0.0
    end
end
```

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`
Expected: FAIL with `UndefVarError: Layers not defined` (or `type Internal has no field Layers`).

- [ ] **Step 3: Create `src/geometries/internal/layers.jl`**

```julia
"""
Near-surface transverse relaxation layer (proposal Eqs. 5–10):

    R2_total(d) = R2_bulk + ΔR2(d),    ΔR2(d) = (ρ/h)·g(d/h)

Implemented: planar walls, step profile g(u) = 1 for 0 ≤ u ≤ 1, no finite RF pulses.
Each side of a wall has its own (ρ, h).
"""
module Layers

"""
    LayerSide(rho, h)

Near-surface layer on one side of a surface: integrated relaxivity `rho` (um/ms) and length scale `h` (um).
A side with `rho == 0` has no layer.
"""
struct LayerSide
    rho :: Float64
    h :: Float64
end

"""
    surface_rate(side)

Excess transverse relaxation rate at the surface, ΔR2(0) = (ρ/h)·g(0), in 1/ms (g(0) = 1 for the step profile).
"""
surface_rate(side::LayerSide) = iszero(side.rho) ? 0.0 : side.rho / side.h

"""
    WallLayer(position, positive, negative)

Layers on the two sides of one wall at local coordinate `position` (um).
`positive` is the side where the local coordinate is larger than `position`.
"""
struct WallLayer
    position :: Float64
    positive :: LayerSide
    negative :: LayerSide
end

"""
    segment_overlap(a0, a1, lo, hi)

A coordinate changes linearly from `a0` (s = 0) to `a1` (s = 1) along a straight segment.
Returns `(s_enter, s_exit)`, the part of s ∈ [0, 1] where lo ≤ a(s) ≤ hi.
There is no overlap when `s_exit <= s_enter`.
"""
function segment_overlap(a0::Float64, a1::Float64, lo::Float64, hi::Float64)
    if a0 == a1
        return (lo <= a0 <= hi) ? (0.0, 1.0) : (0.0, 0.0)
    end
    s_lo = (lo - a0) / (a1 - a0)
    s_hi = (hi - a0) / (a1 - a0)
    return (max(0.0, min(s_lo, s_hi)), min(1.0, max(s_lo, s_hi)))
end

"""
    segment_fraction(a0, a1, lo, hi)

Fraction of the segment (see [`segment_overlap`](@ref)) spent with lo ≤ a(s) ≤ hi.
"""
function segment_fraction(a0::Float64, a1::Float64, lo::Float64, hi::Float64)
    (s_enter, s_exit) = segment_overlap(a0, a1, lo, hi)
    return max(0.0, s_exit - s_enter)
end

"""
    side_exponent(side, d0, d1, dt)

∫ΔR2(d(t)) dt over a straight segment of duration `dt` (ms).
`d0` and `d1` are the distances from the surface at the start and end, positive into this side.
For the step profile this is ΔR2(0)·(time spent with 0 ≤ d ≤ h).
"""
function side_exponent(side::LayerSide, d0::Float64, d1::Float64, dt::Float64)
    iszero(side.rho) && return 0.0
    return surface_rate(side) * dt * segment_fraction(d0, d1, 0.0, side.h)
end

end
```

- [ ] **Step 4: Include the module in `src/geometries/internal/internal.jl`**

Add the include right after `include("fixed_obstruction_groups.jl")`, and the import after the existing `import .FixedObstructionGroups: ...` line:

```julia
include("fixed_obstruction_groups.jl")
include("layers.jl")
include("properties.jl")
```

```julia
import .Layers: LayerSide, WallLayer
```

- [ ] **Step 5: Run the tests and check they pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`
Expected: PASS (all tests in "test_layer.jl: segment maths").

- [ ] **Step 6: Commit**

```bash
git add src/geometries/internal/layers.jl src/geometries/internal/internal.jl test/test_layer.jl test/runtests.jl
git commit -m "Add near-surface layer types and exact step-profile segment integral

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: User parameters on `Walls`, validation, storage on the group

**Files:**
- Modify: `src/geometries/user/obstructions/obstructions.jl:29` (the `ObstructionType(:Wall; ...)` line)
- Modify: `src/geometries/internal/fixed_obstruction_groups.jl:45-90` (struct and constructor)
- Modify: `src/geometries/user/fix.jl:49-52` (`fix_type(::Walls)`), `fix.jl:139` (`apply_properties` signature), `fix.jl:197-208` (constructor call)
- Modify: `test/test_layer.jl` (append)

**Interfaces:**
- Consumes: `Internal.LayerSide`, `Internal.WallLayer` (Task 1).
- Produces:
  - `Walls` fields `layer_rho`, `layer_h` (Float64, default 0), and `layer_rho_positive`, `layer_rho_negative`, `layer_h_positive`, `layer_h_negative` (Float64, default `nothing`).
  - `FixedObstructionGroup.layer`, which is `nothing` or a `Vector{WallLayer}` (one per wall, in wall order).
  - `Fix.wall_layers(walls; density=0.)`, which returns `nothing` or a `Vector{WallLayer}`.

- [ ] **Step 1: Append the failing tests to `test/test_layer.jl`**

```julia
@testset "test_layer.jl: wall layer parameters" begin
    geom(walls) = mr.Simulation([]; geometry=walls, verbose=false).geometry

    @testset "no layer gives nothing" begin
        @test geom(mr.Walls(repeats=2.))[1].layer === nothing
        @test geom(mr.Walls(repeats=2., layer_h=0.5))[1].layer === nothing
    end

    @testset "same layer on both sides" begin
        l = geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5))[1].layer[1]
        @test l.position == 0.0
        @test (l.positive.rho, l.positive.h) == (0.05, 0.5)
        @test (l.negative.rho, l.negative.h) == (0.05, 0.5)
    end

    @testset "one side only" begin
        l = geom(mr.Walls(repeats=2., layer_rho_positive=0.05, layer_h_positive=0.5))[1].layer[1]
        @test (l.positive.rho, l.positive.h) == (0.05, 0.5)
        @test l.negative.rho == 0.0
    end

    @testset "side-specific values override both-sides values" begin
        l = geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5,
            layer_rho_negative=0.2, layer_h_negative=0.4))[1].layer[1]
        @test (l.positive.rho, l.positive.h) == (0.05, 0.5)
        @test (l.negative.rho, l.negative.h) == (0.2, 0.4)
    end

    @testset "per-wall values reach the right wall" begin
        layers = geom(mr.Walls(position=[0., 1.], layer_rho=[0.05, 0.], layer_h=0.4))[1].layer
        @test [l.position for l in layers] == [0.0, 1.0]
        @test layers[1].positive.rho == 0.05
        @test layers[2].positive.rho == 0.0
    end

    @testset "validation" begin
        @test_throws ErrorException geom(mr.Walls(repeats=2., layer_rho=-0.1, layer_h=0.5))
        @test_throws ErrorException geom(mr.Walls(repeats=2., layer_rho=0.1, layer_h=-0.5))
        @test_throws ErrorException geom(mr.Walls(repeats=2., layer_rho=0.1))              # h = 0
        @test_throws ErrorException geom(mr.Walls(repeats=2., layer_rho=0.1, layer_h=2.5)) # h > spacing 2
        @test_throws ErrorException geom(mr.Walls(position=[0., 1.], layer_rho=0.1, layer_h=1.5))
    end

    @testset "warning when spins can stick to the wall" begin
        @test_logs (:warn, r"not applied while a spin is stuck") geom(
            mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5, density=1., dwell_time=1.))
    end

    @testset "JSON round trip keeps the layer fields" begin
        walls = mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5, layer_rho_negative=0.)
        io = IOBuffer()
        mr.write_geometry(io, walls)
        back = mr.read_geometry_json(String(take!(io)))
        @test back.layer_rho.value == 0.05
        @test back.layer_h.value == 0.5
        @test back.layer_rho_negative.value == 0.0
        @test isnothing(back.layer_h_positive.value)
    end
end
```

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`
Expected: FAIL with an error that `Walls` has no keyword/property `layer_h`.

- [ ] **Step 3: Add the fields to `Walls`**

In `src/geometries/user/obstructions/obstructions.jl`, replace

```julia
    ObstructionType(:Wall; ndim=1, volumes=[]),
```

with

```julia
    ObstructionType(:Wall; ndim=1, volumes=[], fields=[
        Field{Float64}(:layer_rho, "Near-surface layer: integrated relaxivity ρ = ∫ΔR2(d)dd (um/ms) on both sides of the wall, unless a side-specific value is set. Zero means no layer.", 0.),
        Field{Float64}(:layer_h, "Near-surface layer: length scale h (um) on both sides of the wall, unless a side-specific value is set.", 0.),
        Field{Float64}(:layer_rho_positive, "Near-surface layer: ρ (um/ms) on the positive side (local coordinate larger than the wall position). Overrides `layer_rho`."),
        Field{Float64}(:layer_rho_negative, "Near-surface layer: ρ (um/ms) on the negative side. Overrides `layer_rho`."),
        Field{Float64}(:layer_h_positive, "Near-surface layer: h (um) on the positive side. Overrides `layer_h`."),
        Field{Float64}(:layer_h_negative, "Near-surface layer: h (um) on the negative side. Overrides `layer_h`."),
    ]),
```

- [ ] **Step 4: Add the `layer` field to `FixedObstructionGroup`**

In `src/geometries/internal/fixed_obstruction_groups.jl`:

Type parameters: replace

```julia
    A <: NamedTuple,
    K
    }
```

with

```julia
    A <: NamedTuple,
    K,
    L
    }
```

Fields: replace

```julia
    # Additional arguments that should be passed around for the `obstructions` to work
    args :: A
```

with

```julia
    # Additional arguments that should be passed around for the `obstructions` to work
    args :: A

    # Near-surface R2 layer: `nothing`, or one `WallLayer` per obstruction (walls only)
    layer :: L
```

Constructor: replace

```julia
    function FixedObstructionGroup(obstructions, repeats, parent_index, original_index, rotation, grid, volume, surface, size_scale, args)
        N = size(rotation, 2)
        repeats = isnothing(repeats) ? nothing : SVector{N, Float64}(repeats)
        new{
            N, typeof(repeats), eltype(obstructions), typeof(grid),
            typeof(volume), typeof(surface), typeof(args), 3 * size(rotation, 2)
        }(
            repeats, parent_index, original_index, 
            rotation, transpose(rotation), 
            grid, volume, surface,
            size_scale, args,
        )
    end
```

with

```julia
    function FixedObstructionGroup(obstructions, repeats, parent_index, original_index, rotation, grid, volume, surface, size_scale, args; layer=nothing)
        N = size(rotation, 2)
        repeats = isnothing(repeats) ? nothing : SVector{N, Float64}(repeats)
        new{
            N, typeof(repeats), eltype(obstructions), typeof(grid),
            typeof(volume), typeof(surface), typeof(args), 3 * size(rotation, 2), typeof(layer)
        }(
            repeats, parent_index, original_index, 
            rotation, transpose(rotation), 
            grid, volume, surface,
            size_scale, args, layer,
        )
    end
```

Also add to the docstring property list above the struct:

```
- `layer`: near-surface R2 layer, `nothing` or one `WallLayer` per obstruction (walls only).
```

- [ ] **Step 5: Resolve, validate and pass the layer in `src/geometries/user/fix.jl`**

Replace

```julia
function fix_type(walls::Walls, index::Int, original_index::Int; kwargs...)
    base_obstructions = fill(Internal.Wall(), length(walls))
    apply_properties(walls, base_obstructions, index, original_index; surface="surface", kwargs...)
end
```

with

```julia
function fix_type(walls::Walls, index::Int, original_index::Int; kwargs...)
    base_obstructions = fill(Internal.Wall(), length(walls))
    layer = wall_layers(walls; density=kwargs[:density])
    apply_properties(walls, base_obstructions, index, original_index; surface="surface", layer=layer, kwargs...)
end

"""
    wall_layers(walls; density=0.)

Resolves the near-surface layer parameters of every wall into `Internal.WallLayer` objects.
Side-specific fields (`layer_rho_positive`, ...) override the both-sides fields (`layer_rho`, `layer_h`).

Returns `nothing` if no wall has a layer, so that simulations without a layer run the unmodified code.
`density` is the global surface density passed to `Simulation` (used only for the stuck-spin warning).
"""
function wall_layers(walls::Walls; density=0.)
    function resolve(quantity::String, side::String, i::Int)
        specific = getproperty(walls, Symbol("layer_" * quantity * "_" * side))[i]
        isnothing(specific) || return Float64(specific)
        both = getproperty(walls, Symbol("layer_" * quantity))[i]
        return isnothing(both) ? 0.0 : Float64(both)
    end
    spacing = size_scale(walls; ignore_user_value=true)
    positions = value_as_vector(walls.position)
    layers = Internal.WallLayer[]
    for i in 1:length(walls)
        sides = map(("positive", "negative")) do side
            rho = resolve("rho", side, i)
            h = resolve("h", side, i)
            rho < 0 && error("Wall $i, $side side: layer_rho must be >= 0, got $rho.")
            h < 0 && error("Wall $i, $side side: layer_h must be >= 0, got $h.")
            rho > 0 && iszero(h) && error("Wall $i, $side side: layer_h must be > 0 when layer_rho > 0.")
            rho > 0 && h > spacing && error("Wall $i, $side side: layer_h = $h um exceeds the spacing between walls ($spacing um), so the layer would pass through a neighbouring wall.")
            Internal.LayerSide(rho, h)
        end
        push!(layers, Internal.WallLayer(Float64(positions[i]), sides...))
    end
    if all(l -> iszero(l.positive.rho) && iszero(l.negative.rho), layers)
        return nothing
    end
    stuck_possible = density > 0 || any(d -> !isnothing(d) && d > 0, value_as_vector(walls.density))
    stuck_possible && @warn "The near-surface layer is not applied while a spin is stuck to a wall (surface density > 0)."
    return layers
end
```

In `apply_properties`, add `layer=nothing` to the keyword list. Replace

```julia
function apply_properties(user_obstructions::ObstructionGroup, internal_obstructions::Vector{<:Internal.FixedObstruction}, index::Int, original_index::Int; surface=nothing, volume=nothing, apply_shift=true, kwargs...)
```

with

```julia
function apply_properties(user_obstructions::ObstructionGroup, internal_obstructions::Vector{<:Internal.FixedObstruction}, index::Int, original_index::Int; surface=nothing, volume=nothing, apply_shift=true, layer=nothing, kwargs...)
```

and pass it to the constructor. Replace

```julia
        size_scale(user_obstructions),
        args
    )
    return result
```

with

```julia
        size_scale(user_obstructions),
        args;
        layer=layer,
    )
    return result
```

- [ ] **Step 6: Run the tests and check they pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`
Expected: PASS (both testsets).

- [ ] **Step 7: Check that existing geometry and CLI tests still pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["collisions", "permeability", "transfer", "cli"])'`
Expected: PASS. The new `Walls` fields appear in JSON output and as CLI options, but must not change any existing behaviour.

- [ ] **Step 8: Commit**

```bash
git add src/geometries/user/obstructions/obstructions.jl src/geometries/internal/fixed_obstruction_groups.jl src/geometries/user/fix.jl test/test_layer.jl
git commit -m "Add side-specific near-surface layer parameters to Walls

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Layer exponent for a wall group and a whole geometry

**Files:**
- Modify: `src/geometries/internal/layers.jl` (append before the final `end`)
- Modify: `src/geometries/internal/internal.jl` (extend the `import .Layers` line)
- Modify: `test/test_layer.jl` (append)

**Interfaces:**
- Consumes: `FixedObstructionGroup.layer` (Task 2); `FixedObstructionGroups.repeating`, `FixedObstructionGroups.rotate_from_global`.
- Produces:
  - `Layers.layer_exponent(group::FixedObstructionGroup, start::SVector{3,Float64}, dest::SVector{3,Float64}, dt::Float64)::Float64`
  - `Layers.layer_exponent(geometry::FixedGeometry, start, dest, dt)::Float64`, the sum over groups (Eq. 10)
  - `Layers.has_layer(geometry::FixedGeometry)::Bool`

- [ ] **Step 1: Append the failing tests to `test/test_layer.jl`**

All expected values below use the wall at x = 0 repeating every w = 2 µm, ΔR₂(0) = ρ/h, and exponent = ΔR₂(0) · dt · (fraction of the segment inside the layer).

```julia
@testset "test_layer.jl: layer exponent on walls" begin
    Layers = mr.Geometries.Internal.Layers
    geom(walls) = mr.Simulation([]; geometry=walls, verbose=false).geometry
    v(x) = SVector{3, Float64}(x, 0.3, -1.2)
    both = geom(mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5))   # ΔR2(0) = 0.1 /ms on both sides

    @testset "inside, outside and partly inside" begin
        @test Layers.layer_exponent(both, v(0.1), v(0.3), 0.2) ≈ 0.1 * 0.2        # positive side of wall 0
        @test Layers.layer_exponent(both, v(1.7), v(1.9), 0.2) ≈ 0.1 * 0.2        # negative side of wall 2
        @test Layers.layer_exponent(both, v(0.6), v(1.4), 1.0) == 0.0             # middle of the gap
        @test Layers.layer_exponent(both, v(0.4), v(0.6), 1.0) ≈ 0.1 * 1.0 * 0.5
    end

    @testset "neither endpoint inside a layer" begin
        # 0.6 → -0.6 crosses both layers of wall 0: 1.0 of 1.2 um inside
        @test Layers.layer_exponent(both, v(0.6), v(-0.6), 1.2) ≈ 0.1 * 1.2 * (1.0 / 1.2)
    end

    @testset "repeated wall images" begin
        @test Layers.layer_exponent(both, v(10.1), v(10.3), 0.2) ≈ 0.1 * 0.2
        @test Layers.layer_exponent(both, v(-7.9), v(-7.7), 0.2) ≈ 0.1 * 0.2
        # segment longer than the spacing: 0.1 → 4.1 is inside layers for 2.0 of its 4.0 um
        @test Layers.layer_exponent(both, v(0.1), v(4.1), 1.0) ≈ 0.1 * 1.0 * 0.5
    end

    @testset "one side only" begin
        pos = geom(mr.Walls(repeats=2., layer_rho_positive=0.05, layer_h_positive=0.5))
        @test Layers.layer_exponent(pos, v(0.1), v(0.3), 0.2) ≈ 0.1 * 0.2
        @test Layers.layer_exponent(pos, v(1.7), v(1.9), 0.2) == 0.0
    end

    @testset "different on each side" begin
        asym = geom(mr.Walls(repeats=2., layer_rho_positive=0.05, layer_h_positive=0.5,
            layer_rho_negative=0.2, layer_h_negative=0.4))                       # 0.1 /ms and 0.5 /ms
        @test Layers.layer_exponent(asym, v(0.1), v(0.3), 1.0) ≈ 0.1
        @test Layers.layer_exponent(asym, v(1.7), v(1.9), 1.0) ≈ 0.5
    end

    @testset "overlapping layers add (Eq. 10)" begin
        wide = geom(mr.Walls(repeats=2., layer_rho=0.15, layer_h=1.5))           # 0.1 /ms, h > w/2
        @test Layers.layer_exponent(wide, v(0.9), v(1.1), 1.0) ≈ 2 * 0.1         # inside both layers
    end

    @testset "h = w/2 covers every point exactly once" begin
        half = geom(mr.Walls(repeats=2., layer_rho=0.1, layer_h=1.0))            # 0.1 /ms
        for (a, b) in ((0.05, 0.3), (0.8, 1.2), (1.0, 1.9), (0.01, 1.99))
            @test Layers.layer_exponent(half, v(a), v(b), 1.0) ≈ 0.1 rtol=1e-14
        end
    end

    @testset "rotation and shifted position" begin
        rot = geom(mr.Walls(repeats=2., rotation=:y, layer_rho=0.05, layer_h=0.5))
        @test Layers.layer_exponent(rot, SVector(5.0, 0.1, 7.0), SVector(9.0, 0.3, 7.0), 0.2) ≈ 0.1 * 0.2
        shifted = geom(mr.Walls(repeats=2., position=0.5, layer_rho=0.05, layer_h=0.5))
        @test Layers.layer_exponent(shifted, v(0.6), v(0.8), 0.2) ≈ 0.1 * 0.2
    end

    @testset "per-wall values" begin
        two = geom(mr.Walls(position=[0., 1.], layer_rho=[0.05, 0.], layer_h=0.4))  # only wall at 0
        @test Layers.layer_exponent(two, v(0.1), v(0.3), 1.0) ≈ 0.05 / 0.4
        @test Layers.layer_exponent(two, v(0.7), v(0.9), 1.0) == 0.0
    end

    @testset "no layer" begin
        none = geom(mr.Walls(repeats=2.))
        @test !Layers.has_layer(none)
        @test Layers.layer_exponent(none, v(0.1), v(0.3), 0.2) == 0.0
        @test Layers.has_layer(both)
    end
end
```

`test/runtests.jl` already has `using StaticArrays`, so `SVector` is available.

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`
Expected: FAIL with `UndefVarError: layer_exponent not defined` (in the new testset).

- [ ] **Step 3: Append the group and geometry functions to `src/geometries/internal/layers.jl`**

Add near the top of the module, right after `module Layers`:

```julia
import StaticArrays: SVector
import ..FixedObstructionGroups: FixedObstructionGroup, FixedGeometry, repeating, rotate_from_global
```

Append before the final `end` of the module:

```julia
"""
    layer_exponent(group/geometry, start, dest, dt)

∫ΔR2(d(t)) dt (dimensionless) for a spin moving in a straight line from `start` to `dest` (global coordinates, um)
in `dt` ms. For a geometry, the contributions of all groups and all walls are added (proposal Eq. 10).
Returns 0 when no layer is defined.
"""
layer_exponent(g::FixedObstructionGroup, start::SVector{3, Float64}, dest::SVector{3, Float64}, dt::Float64) = layer_exponent(g.layer, g, start, dest, dt)

layer_exponent(::Nothing, ::FixedObstructionGroup, ::SVector{3, Float64}, ::SVector{3, Float64}, ::Float64) = 0.0

function layer_exponent(layers::Vector{WallLayer}, g::FixedObstructionGroup, start::SVector{3, Float64}, dest::SVector{3, Float64}, dt::Float64)
    x0 = rotate_from_global(g, start)[1]
    x1 = rotate_from_global(g, dest)[1]
    lo = min(x0, x1)
    hi = max(x0, x1)
    spacing = repeating(g) ? g.repeats[1] : 0.0
    total = 0.0
    for layer in layers
        if repeating(g)
            # wall images c = position + k·spacing whose layers can reach [lo, hi]
            reach = max(layer.positive.h, layer.negative.h)
            kmin = ceil(Int, (lo - reach - layer.position) / spacing)
            kmax = floor(Int, (hi + reach - layer.position) / spacing)
        else
            kmin = kmax = 0
        end
        for k in kmin:kmax
            c = layer.position + k * spacing
            total += side_exponent(layer.positive, x0 - c, x1 - c, dt)
            total += side_exponent(layer.negative, c - x0, c - x1, dt)
        end
    end
    return total
end

function layer_exponent(geometry::FixedGeometry, start::SVector{3, Float64}, dest::SVector{3, Float64}, dt::Float64)
    return sum(g -> layer_exponent(g, start, dest, dt), geometry; init=0.0)
end

"""
    has_layer(geometry)

Whether any group in the [`FixedGeometry`](@ref) has a near-surface layer.
"""
has_layer(geometry::FixedGeometry) = any(g -> !isnothing(g.layer), geometry)
```

- [ ] **Step 4: Export the new names from `Internal`**

In `src/geometries/internal/internal.jl`, replace

```julia
import .Layers: LayerSide, WallLayer
```

with

```julia
import .Layers: LayerSide, WallLayer, layer_exponent, has_layer
```

- [ ] **Step 5: Run the tests and check they pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`
Expected: PASS (all three testsets).

- [ ] **Step 6: Commit**

```bash
git add src/geometries/internal/layers.jl src/geometries/internal/internal.jl test/test_layer.jl
git commit -m "Compute the near-surface layer exponent for wall groups

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Apply the layer in the spin step

**Files:**
- Modify: `src/evolve.jl:13` and `:20` (imports), `:481` (new `apply_layer!` before `draw_step!`), `:547` (call)
- Modify: `test/test_layer.jl` (append)

**Interfaces:**
- Consumes: `Internal.layer_exponent`, `Internal.has_layer` (Task 3).
- Produces: `Evolve.apply_layer!(spin, geometry, start, dest, dt, parts)`, which multiplies every orientation's `transverse` by exp(−exponent) and changes nothing else.

- [ ] **Step 1: Append the failing tests to `test/test_layer.jl`**

```julia
@testset "test_layer.jl: layer applied during the simulation" begin
    per_spin_transverse(snap) = [s.orientations[1].transverse for s in snap.spins]

    @testset "crossing with reflection in one step (P2.2.4)" begin
        # wall at 0, layer [0, 0.5] on the positive side, ΔR2(0) = 0.1 /ms
        walls = mr.Walls(position=0., layer_rho_positive=0.05, layer_h_positive=0.5)
        seq = build_sequence([1.2, :readout])
        sim = mr.Simulation(seq; geometry=walls, diffusivity=3., verbose=false)
        part = mr.parts([seq], 0., mr.TimeStep(1.2, Inf))[1]
        @test part.duration ≈ 1.2
        spin = mr.Spin(position=[0.8, 0., 0.], transverse=1., longitudinal=0.)
        # proposed move 0.8 → -0.4 reflects at 0 and ends at 0.4.
        # in layer: last 0.5 of the 0.8 um before the wall (time 0.5) + all 0.4 um after (time 0.4) = 0.9 of 1.2 ms
        mr.Evolve.draw_step!(spin, sim, part, [3.], [-0.4, 0., 0.])
        @test spin.position[1] ≈ 0.4 atol=1e-12
        @test spin.orientations[1].transverse ≈ exp(-0.1 * 0.9) rtol=1e-12
        @test spin.orientations[1].phase == 0.0
    end

    @testset "no layer: bit-identical to the unmodified code (P2.2.8)" begin
        function run(walls)
            sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=walls, diffusivity=3., R2=1/80, verbose=false)
            Random.seed!(7)
            snap = mr.Snapshot(2000, sim, 500; transverse=1., longitudinal=0.)
            return mr.readout(snap, sim, [50.]; return_snapshot=true)[1]
        end
        a = run(mr.Walls(repeats=2.))
        b = run(mr.Walls(repeats=2., layer_h=0.5))          # ρ = 0: no layer
        @test mr.position.(a) == mr.position.(b)
        @test per_spin_transverse(a) == per_spin_transverse(b)
    end

    @testset "h = w/2: every spin decays at exactly ΔR2(0)" begin
        walls = mr.Walls(repeats=2., layer_rho=0.1, layer_h=1.0)                # ΔR2(0) = 0.1 /ms
        sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=walls, diffusivity=3., verbose=false)
        Random.seed!(1)
        snap = mr.Snapshot(2000, sim, 500; transverse=1., longitudinal=0.)
        res = mr.readout(snap, sim, [10., 50.]; return_snapshot=true)
        for (t, s) in zip((10., 50.), res)
            @test all(isapprox.(per_spin_transverse(s), exp(-0.1 * t); rtol=1e-10))
        end
    end

    @testset "bulk R2 and the layer multiply exactly (R2_total = R2_bulk + ΔR2)" begin
        walls = mr.Walls(repeats=2., layer_rho=0.03, layer_h=0.3)
        function run(R2)
            sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=walls, diffusivity=3., R2=R2, verbose=false)
            Random.seed!(3)
            snap = mr.Snapshot(2000, sim, 500; transverse=1., longitudinal=0.)
            return mr.readout(snap, sim, [50.]; return_snapshot=true)[1]
        end
        layer_only = run(0.)
        with_bulk = run(1/80)
        @test mr.position.(layer_only) == mr.position.(with_bulk)
        @test all(isapprox.(per_spin_transverse(with_bulk), per_spin_transverse(layer_only) .* exp(-50 / 80); rtol=1e-10))
        @test minimum(per_spin_transverse(layer_only)) < 1.0                     # the layer did act
    end

    @testset "positive-only and negative-only layers are mirror images" begin
        function mean_signal(walls)
            sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=walls, diffusivity=3., verbose=false)
            Random.seed!(11)
            snap = mr.Snapshot(20000, sim, 500; transverse=1., longitudinal=0.)
            m = per_spin_transverse(mr.readout(snap, sim, [50.]; return_snapshot=true)[1])
            return (mean(m), std(m) / sqrt(length(m)))
        end
        (p, ep) = mean_signal(mr.Walls(repeats=2., layer_rho_positive=0.05, layer_h_positive=0.5))
        (n, en) = mean_signal(mr.Walls(repeats=2., layer_rho_negative=0.05, layer_h_negative=0.5))
        @test abs(p - n) < 5 * sqrt(ep^2 + en^2)
        @test p < 0.99                                                          # the layer did act
    end

    @testset "finite RF pulse with a layer is rejected" begin
        seq = mr.SequenceParts.SequenceWaveform(
            (([], []), ([], []), ([], [])),
            [(0., 9., [mr.SequenceParts.ConstantPulse(0.25 / 9, 0., 0.)])],
            [], [10.], 10.,
        )
        sim = mr.Simulation(seq; geometry=mr.Walls(repeats=2., layer_rho=0.05, layer_h=0.5), verbose=false)
        @test_throws Exception mr.readout(100, sim, [10.])
    end
end
```

`Random`, `Statistics` and `build_sequence` come from `test/runtests.jl`.

- [ ] **Step 2: Run the tests and check they fail**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`
Expected: FAIL. The reflection test gives transverse = 1.0 instead of exp(−0.09), the h = w/2 test gives 1.0, and the finite-RF test does not throw.

- [ ] **Step 3: Add the imports in `src/evolve.jl`**

Replace line 13

```julia
import ..SequenceParts: SequencePart, MultSequencePart, InstantSequencePart, get_readouts, IndexedReadout, empty_sequence, GradientEvent, PulseEvent, parts, repetition_time
```

with

```julia
import ..SequenceParts: SequencePart, MultSequencePart, InstantSequencePart, get_readouts, IndexedReadout, empty_sequence, GradientEvent, PulseEvent, parts, repetition_time, PulsePart
```

Replace line 20

```julia
import ..Geometries.Internal: Reflection, detect_intersection, empty_intersection, has_intersection, surface_relaxation, permeability, surface_density, direction, previous_hit, dwell_time, empty_reflection, FixedGeometry
```

with

```julia
import ..Geometries.Internal: Reflection, detect_intersection, empty_intersection, has_intersection, surface_relaxation, permeability, surface_density, direction, previous_hit, dwell_time, empty_reflection, FixedGeometry, layer_exponent, has_layer
```

- [ ] **Step 4: Add `apply_layer!` just above `function draw_step!(spin::Spin{N}, ...)` (line 482)**

```julia
"""
    apply_layer!(spin, geometry, start, dest, dt, parts)

Multiplies the transverse magnetisation by exp(-∫ΔR2(d(t))dt) for the straight free segment `start → dest`
of duration `dt` (ms). Bulk R2 is applied separately by [`relax!`](@ref); the two factors multiply,
which gives R2_total = R2_bulk + ΔR2(d). Phase and longitudinal magnetisation are not changed.
"""
function apply_layer!(spin::Spin, geometry::FixedGeometry, start::SVector{3, Float64}, dest::SVector{3, Float64}, dt::Float64, parts::MultSequencePart)
    has_layer(geometry) || return
    if any(p -> p isa PulsePart, parts.parts)
        error("The near-surface R2 layer does not support finite RF pulses yet. Use instantaneous pulses.")
    end
    exponent = layer_exponent(geometry, start, dest, dt)
    iszero(exponent) && return
    attenuation = exp(-exponent)
    for orientation in spin.orientations
        orientation.transverse *= attenuation
    end
end

```

- [ ] **Step 5: Call it after `relax!` in `draw_step!`**

Replace

```julia
            # spin relaxation
            relax!(spin, collision_pos, simulation, parts, fraction_timestep, next_fraction_timestep, B0s)
```

with

```julia
            # spin relaxation
            relax!(spin, collision_pos, simulation, parts, fraction_timestep, next_fraction_timestep, B0s)
            apply_layer!(spin, simulation.geometry, current_pos, collision_pos, (next_fraction_timestep - fraction_timestep) * timestep, parts)
```

`current_pos` is the start of this free segment, `collision_pos` its end (`new_pos` when nothing is hit), and `(next_fraction_timestep - fraction_timestep) * timestep` its duration in ms. These are the same three quantities `relax!` uses.

- [ ] **Step 6: Run the tests and check they pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer"])'`
Expected: PASS (all four testsets).

- [ ] **Step 7: Commit**

```bash
git add src/evolve.jl test/test_layer.jl
git commit -m "Apply the near-surface R2 layer along each free segment

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Regression, runtime check and records

**Files:**
- Modify: `research/project-progress-tracker.md` (P2.1.1, P2.1.2, P2.2.1–P2.2.8 status)
- Modify: `research/phase1_code_log.md` (§7: note that the insertion points are implemented)

**Interfaces:**
- Consumes: everything above. Produces no code.

- [ ] **Step 1: Run the whole test suite except plots**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["no-plots"])'`
Expected: every test passes. Baseline B0 had 442 passing in its five suites, and the totals now include the new `layer` tests. If a pre-existing test fails, stop and investigate before continuing.

- [ ] **Step 2: Check the no-layer runtime against baseline B1**

Run:

```bash
julia --project=research/baseline -t 8 -e '
import MCMRSimulator as mr; using Random
sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=mr.Walls(repeats=2.), diffusivity=3., verbose=false)
times = Float64[]
for seed in 1:4
    Random.seed!(seed); snap = mr.Snapshot(100000, sim, 500; transverse=1., longitudinal=0.)
    t0 = time(); mr.readout(snap, sim, [0., 10., 20., 30., 40., 50.]); push!(times, time() - t0)
end
println("µs per spin-step: ", sum(times[2:end]) / 3 / (100000 * 50 / 0.04) * 1e6)'
```

Expected: about 0.013 µs per spin-step (baseline B1: 0.0130). A value more than 10% higher means the no-layer path is slower. Check that `has_layer` returns `false` early.

- [ ] **Step 3: Update the tracker**

In `research/project-progress-tracker.md`, set P2.1.1, P2.1.2 and P2.2.1–P2.2.8 to ☑. The rows already describe this design and point to this plan, so only the status changes. In P2.2.8, add the measured µs per spin-step from Step 2. Leave P2.2.9 (single-spin trace script) unticked.

- [ ] **Step 4: Note the implementation in the code log**

In `research/phase1_code_log.md` §7, add below the table:

```markdown
**Status (2026-10-01):** 1a, 1b, 1d, 2a, 2c, 3a and 5 are implemented for walls (step profile, per-side ρ and h) in `src/geometries/internal/layers.jl`, `src/geometries/user/fix.jl` (`wall_layers`) and `src/evolve.jl` (`apply_layer!`). 2c uses the exact linear overlap (`segment_overlap`) instead of the hit grid. 3b and 4 are not yet implemented.
```

- [ ] **Step 5: Commit**

```bash
git add research/project-progress-tracker.md research/phase1_code_log.md
git commit -m "Record near-surface layer implementation status

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
