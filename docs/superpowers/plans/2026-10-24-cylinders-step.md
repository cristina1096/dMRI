# Cylinders, step profile (Phase 3.1, 3.2, 3.4-step) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Put the step-profile near-surface layer on cylinders (single and square/rectangular arrays, inner and outer side, overlapping outer layers summed), verify it, and record the cylinder baselines and the conventional mapping on curved surfaces.

**Architecture:** The layer exponent of a straight piece is computed in the cylinder's transverse plane. Along the piece the squared distance from the axis is a quadratic in s, r(s)² = a s² + b s + q. The time spent in a shell [r_lo, r_hi] is therefore exact from quadratic roots. Each side has the D3.1 normalisation ΔR₂ = (ρ/h)·g/κ, with κ = 1 ∓ h/(2R) for the step profile. Contributions of all cylinder images within reach are summed (Eq. 10).

New code:
- `src/geometries/internal/layers.jl` (internal maths);
- `src/geometries/user/fix.jl` (user fields → `CylinderLayer`);
- `research/phase3/` (harness and scripts).

**Tech Stack:** Julia 1.12, MCMRSimulator.jl (this repo), StaticArrays, CairoMakie (figures), research environment `research/baseline`.

**Spec:**
- `research/phase3_record.md`: design decisions D3.1–D3.6, binding;
- `research/project-progress-tracker.md`: Phase 3 sections 3.1, 3.2, 3.4 and the risks table.

## Global Constraints

- **Branch.** Work only on branch `cc/near-surface-r2`. Never modify or commit on other branches. Do not push; the user pushes.
- **Commits.** Every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **No RF.** Finite RF pulses with a layer stay rejected; this already holds for any geometry via `apply_layer!`.
- **Step profile only.** On cylinders `layer_shape` must be `"step"` (D3.2). Linear and exponential come in a later plan and are rejected with an error until then.
- **ρ per unit surface area (D3.1).** ΔR₂(d) = (ρ/h)·g(d/h)/κ, κ = 1 − c̄/R inside, 1 + c̄/R outside, c̄ = h/2 for the step profile.
- **Sides (D3.4).** Inner layer: 0 ≤ R − r ≤ h_inner. Outer layer: 0 ≤ r − R ≤ h_outer. Side-specific fields override both-sides fields. Both sides are symmetric by default.
- **Validation.** h_inner ≤ R. h_outer ≤ the smallest surface gap |c_i − c_j| − R_i − R_j to any other cylinder or image; this is the wall rule h ≤ spacing.
- **Arrays (D3.6).** Square and rectangular via `Cylinders(repeats=[sx, sy])`.
- **Model framing.**
  - Baselines (MCMR θ) are references, not truth.
  - No fast-diffusion or Brownstein–Tarr assumption anywhere.
  - Records give Design / Method / Result / Conclusion with formulas, and figures where useful.
- **Package tests.**
  - `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders"])'`
  - Full suite without plots: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["no-plots"])'`
- **Research tests.** `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`
- **Research scripts.** `julia --project=research/baseline -t 8 research/phase3/<script>.jl`. Each writes `research/results/phase3/<name>/summary.json` with `provenance()`.
- **Plan end.** A fresh subagent reviewer reviews the whole branch (user rule), not a self-review.

## Review Focus

1. **Nearly axial motion.** A piece with a tiny but non-zero transverse displacement (|Δ⊥|² ≈ 1e-30) must give the same shell fraction as a stationary point, with no NaN or Inf. Test in Task 1.
2. **Piece starting on the surface.** After a reflection the next piece starts exactly at r = R. The inner and outer fractions must still be exact (the boundary has measure zero). Test in Task 1.
3. **Rotated cylinder axis.** A cylinder along x must measure r in the (y, z) plane; motion along the axis must not count. Test in Task 3.
4. **Several cylinders per unit cell in a repeating array.** Every image of every cylinder within reach must be summed. Test in Task 3 (two cylinders per cell vs brute force).
5. **Permeable cylinders with a layer.** Spins cross the surface, so a piece may lie on either side. The run must work and give finite M⊥ ≤ 1 for every spin. Test in Task 4.

---

## File Structure

- **Modify** `src/geometries/internal/layers.jl`:
  - `CylinderLayer`, `centroid`, `kappa`, `cylinder_surface_rate`;
  - `disk_length`, `shell_fraction`, `cylinder_exponent`;
  - `image_centres`, `layer_exponent(::Vector{CylinderLayer}, ...)`;
  - `max_layer_rate` for cylinder groups.
- **Modify** `src/geometries/internal/internal.jl`: export `CylinderLayer`.
- **Modify** `src/geometries/user/obstructions/obstructions.jl`: cylinder layer fields.
- **Modify** `src/geometries/user/fix.jl`:
  - `obstruction_shifts` (refactored out of `apply_properties`);
  - `cylinder_layers`, `check_outer_reach`;
  - cylinder `fix_type` passes `layer`.
- **Create** `test/test_layer_cylinders.jl`. **Modify** `test/runtests.jl` (add `"layer_cylinders"`) and `test/test_cli.jl`.
- **Create** research code:
  - harness: `research/phase3/harness.jl`;
  - tests: `research/phase3/test/runtests.jl`, `research/phase3/test/test_harness.jl`;
  - scripts: `research/phase3/p3_1_1_null_bulk.jl`, `p3_1_3_density.jl`, `p3_1_4_theta_reference.jl`, `p3_1_5_runtime.jl`, `p3_1_6_mapping.jl`, `p3_4_1_trace.jl`, `p3_4_exact.jl`, `p3_4_3_identity.jl`.
- **Modify** `research/phase3_record.md` and `research/project-progress-tracker.md`.

---

### Task 1: Shell geometry along a straight piece (P3.2.2)

**Files:**
- Modify: `src/geometries/internal/layers.jl` (append after `side_exponent`; add `import LinearAlgebra: norm` next to the StaticArrays import)
- Create: `test/test_layer_cylinders.jl`
- Modify: `test/runtests.jl` (add `"layer_cylinders"` after `"layer"` in `all_tests`)

**Interfaces:**
- Consumes: `LayerSide(rho, h)`, `surface_rate(side)` (existing).
- Produces:
  - `CylinderLayer(centre::SVector{2,Float64}, radius::Float64, inner::LayerSide, outer::LayerSide)`;
  - `centroid(side::LayerSide)::Float64`;
  - `kappa(side, radius, inside::Bool)::Float64`;
  - `cylinder_surface_rate(side, radius, inside::Bool)::Float64`;
  - `disk_length(a, b, q, X)::Float64`;
  - `shell_fraction(p0, p1, centre, r_lo, r_hi)::Float64` (2D SVectors);
  - `cylinder_exponent(layer, centre, p0, p1, dt)::Float64`.

- [ ] **Step 1: Write the failing tests**

Create `test/test_layer_cylinders.jl`:

```julia
@testset "test_layer_cylinders.jl: shell maths" begin
    Layers = mr.Geometries.Internal.Layers
    P(x, y) = SVector{2, Float64}(x, y)
    O = P(0., 0.)

    @testset "normalisation κ and surface rate (D3.1)" begin
        side = Layers.LayerSide(0.02, 0.2)
        @test Layers.centroid(side) == 0.1
        @test Layers.kappa(side, 1.0, true) ≈ 0.9
        @test Layers.kappa(side, 1.0, false) ≈ 1.1
        @test Layers.cylinder_surface_rate(side, 1.0, true) ≈ 0.02 / (0.2 * 0.9)
        @test Layers.cylinder_surface_rate(side, 1.0, false) ≈ 0.02 / (0.2 * 1.1)
        @test Layers.cylinder_surface_rate(Layers.LayerSide(0.0, 0.0), 1.0, true) == 0.0
        @test_throws ErrorException Layers.centroid(Layers.LayerSide(0.02, 0.2, :linear))
    end

    @testset "disk_length" begin
        @test Layers.disk_length(1.0, 0.0, 0.0, 0.25) ≈ 0.5          # r² = s², s ≤ 0.5
        @test Layers.disk_length(1.0, 0.0, 2.0, 1.0) == 0.0           # never inside
        @test Layers.disk_length(0.0, 0.0, 0.5, 1.0) == 1.0           # stationary inside
        @test Layers.disk_length(0.0, 0.0, 1.5, 1.0) == 0.0           # stationary outside
    end

    @testset "shell_fraction: exact cases" begin
        # radial through the centre, shell 0.8 ≤ r ≤ 1: |x| ∈ [0.8, 1] → 0.4 of a length-2 piece
        @test Layers.shell_fraction(P(-1., 0.), P(1., 0.), O, 0.8, 1.0) ≈ 0.2
        # chord at y = 0.5 from x = −2 to 2: inside r ≤ 1 for |x| ≤ √0.75, inside r ≤ 0.8 for |x| ≤ √0.39
        @test Layers.shell_fraction(P(-2., 0.5), P(2., 0.5), O, 0.8, 1.0) ≈ (2sqrt(0.75) - 2sqrt(0.39)) / 4
        # tangent to the outer circle r = 1: measure zero
        @test Layers.shell_fraction(P(-2., 1.0), P(2., 1.0), O, 0.8, 1.0) == 0.0
        # tangent to the inner circle r = 0.8: the whole chord |x| ≤ 0.6 is in the shell
        @test Layers.shell_fraction(P(-2., 0.8), P(2., 0.8), O, 0.8, 1.0) ≈ 0.3
        # outer shell 1 ≤ r ≤ 1.3, radial piece 1 → 2
        @test Layers.shell_fraction(P(1., 0.), P(2., 0.), O, 1.0, 1.3) ≈ 0.3
        # r_lo = 0 (shell = whole disk, h = R)
        @test Layers.shell_fraction(P(-0.5, 0.), P(0.5, 0.), O, 0.0, 1.0) == 1.0
        # stationary
        @test Layers.shell_fraction(P(0.9, 0.), P(0.9, 0.), O, 0.8, 1.0) == 1.0
        @test Layers.shell_fraction(P(0.5, 0.), P(0.5, 0.), O, 0.8, 1.0) == 0.0
        # shifted centre
        @test Layers.shell_fraction(P(4., 3.), P(6., 3.), P(5., 3.), 0.8, 1.0) ≈ 0.2
    end

    @testset "Review Focus 1: nearly axial motion equals stationary" begin
        f = Layers.shell_fraction(P(0.9, 0.), P(0.9, 1e-15), O, 0.8, 1.0)
        @test isfinite(f)
        @test f ≈ 1.0
        @test Layers.shell_fraction(P(0.5, 0.), P(0.5 + 1e-15, 0.), O, 0.8, 1.0) == 0.0
    end

    @testset "Review Focus 2: piece starting exactly on the surface" begin
        # after a reflection at r = R = 1, moving inward to r = 0.7: inner shell [0.8, 1] → 0.2/0.3
        @test Layers.shell_fraction(P(1., 0.), P(0.7, 0.), O, 0.8, 1.0) ≈ 2 / 3
        @test Layers.shell_fraction(P(1., 0.), P(0.7, 0.), O, 1.0, 1.3) ≈ 0.0 atol=1e-12
        # moving outward from the surface to r = 1.5: outer shell [1, 1.3] → 0.3/0.5
        θ = 0.7
        u = P(cos(θ), sin(θ))
        @test Layers.shell_fraction(u, 1.5 .* u, O, 1.0, 1.3) ≈ 0.6
        @test Layers.shell_fraction(u, 1.5 .* u, O, 0.8, 1.0) ≈ 0.0 atol=1e-12
    end

    @testset "shell_fraction vs brute-force sampling" begin
        Random.seed!(2026)
        K = 100_000
        maxerr = 0.0
        for _ in 1:1000
            p0 = P((4 .* rand(2) .- 2)...)
            p1 = P((4 .* rand(2) .- 2)...)
            for (lo, hi) in ((0.8, 1.0), (1.0, 1.3), (0.0, 1.0))
                brute = count(k -> lo <= norm(p0 + (k - 0.5) / K * (p1 - p0)) <= hi, 1:K) / K
                maxerr = max(maxerr, abs(Layers.shell_fraction(p0, p1, O, lo, hi) - brute))
            end
        end
        @test maxerr < 3e-5          # ≤ 4 boundary crossings, each ≤ 1/(2K)
    end

    @testset "cylinder_exponent: both sides" begin
        layer = Layers.CylinderLayer(O, 1.0, Layers.LayerSide(0.02, 0.2), Layers.LayerSide(0.02, 0.3))
        rin = 0.02 / (0.2 * 0.9)
        rout = 0.02 / (0.3 * 1.15)
        # radial piece r = 0.5 → 1.5 of duration 2 ms: 0.2 in the inner shell, 0.3 in the outer shell (of 1.0)
        @test Layers.cylinder_exponent(layer, O, P(0.5, 0.), P(1.5, 0.), 2.0) ≈ 2.0 * (0.2 * rin + 0.3 * rout)
        # inner side off
        layer_out = Layers.CylinderLayer(O, 1.0, Layers.LayerSide(0.0, 0.0), Layers.LayerSide(0.02, 0.3))
        @test Layers.cylinder_exponent(layer_out, O, P(0.5, 0.), P(0.9, 0.), 1.0) == 0.0
    end
end
```

Add `"layer_cylinders",` after `"layer",` in `all_tests` in `test/runtests.jl`.

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders"])'`
Expected: FAIL with `UndefVarError: centroid not defined` (or `CylinderLayer`).

- [ ] **Step 3: Implement**

In `src/geometries/internal/layers.jl`:
- add `import LinearAlgebra: norm` after `import StaticArrays: SVector`;
- update the module docstring's "Implemented:" line to "planar walls (all profiles) and cylinders (step profile)";
- append after `side_exponent`:

```julia
"""
    CylinderLayer(centre, radius, inner, outer)

Near-surface layers of one cylinder of radius `radius` (um) with its axis at `centre` (local 2D coordinates of its
group, um): `inner` inside the cylinder (0 ≤ d = R − r ≤ h), `outer` outside it (0 ≤ d = r − R ≤ h).
ρ is per unit surface area (D3.1): ΔR2(d) = (ρ/h)·g(d/h)/κ with κ = 1 ∓ c̄/R (− inside, + outside).
"""
struct CylinderLayer
    centre :: SVector{2, Float64}
    radius :: Float64
    inner :: LayerSide
    outer :: LayerSide
end

"""
    centroid(side)

Profile centroid c̄ = h∫u g(u) du (um). Only the step profile (c̄ = h/2) is supported on cylinders so far (D3.2).
"""
function centroid(side::LayerSide)
    side.shape == :step && return side.h / 2
    error("Layer profile $(side.shape) is not supported on cylinders yet (step only).")
end

"""
    kappa(side, radius, inside)

Curvature normalisation κ = ∫g(u)(1 ∓ u·h/R) du = 1 ∓ c̄/R (D3.1): − inside the cylinder, + outside.
"""
kappa(side::LayerSide, radius::Float64, inside::Bool) = 1 + (inside ? -1 : 1) * centroid(side) / radius

"""
    cylinder_surface_rate(side, radius, inside)

ΔR2 at the cylinder surface on one side, (ρ/h)·g(0)/κ (1/ms); 0 without a layer.
"""
cylinder_surface_rate(side::LayerSide, radius::Float64, inside::Bool) =
    iszero(side.rho) ? 0.0 : surface_rate(side) / kappa(side, radius, inside)

"""
    disk_length(a, b, q, X)

Length of {s ∈ [0, 1] : a s² + b s + q ≤ X} for a ≥ 0. The squared distance from a cylinder axis along a straight
piece is this convex quadratic, so the set is one interval; its ends are the roots, computed with the
cancellation-free form t = −(b + sign(b)√disc)/2, roots t/a and (q − X)/t.
"""
function disk_length(a::Float64, b::Float64, q::Float64, X::Float64)
    if iszero(a)
        return q <= X ? 1.0 : 0.0
    end
    disc = b * b - 4 * a * (q - X)
    disc <= 0 && return 0.0
    t = -(b + copysign(sqrt(disc), b)) / 2
    (s1, s2) = minmax(t / a, (q - X) / t)
    return max(0.0, min(1.0, s2) - max(0.0, s1))
end

"""
    shell_fraction(p0, p1, centre, r_lo, r_hi)

Fraction of the straight piece p(s) = p0 + s(p1 − p0), s ∈ [0, 1] (2D, um), with r_lo ≤ |p(s) − centre| ≤ r_hi.
Exact: |p(s) − centre|² = a s² + b s + q, so the fraction is disk_length(r_hi²) − disk_length(r_lo²).
"""
function shell_fraction(p0::SVector{2, Float64}, p1::SVector{2, Float64}, centre::SVector{2, Float64}, r_lo::Float64, r_hi::Float64)
    u = p0 - centre
    Δ = p1 - p0
    a = Δ[1] * Δ[1] + Δ[2] * Δ[2]
    b = 2 * (u[1] * Δ[1] + u[2] * Δ[2])
    q = u[1] * u[1] + u[2] * u[2]
    if iszero(a)
        return (r_lo^2 <= q <= r_hi^2) ? 1.0 : 0.0
    end
    inner = r_lo > 0 ? disk_length(a, b, q, r_lo^2) : 0.0
    return max(0.0, disk_length(a, b, q, r_hi^2) - inner)
end

"""
    cylinder_exponent(layer, centre, p0, p1, dt)

∫ΔR2 dt of both layers of one cylinder (axis at `centre`, local 2D) over the straight piece p0 → p1 of duration `dt`
(ms). Step profile: ΔR2(0)·dt·(fraction of the piece inside the shell), exact.
"""
function cylinder_exponent(layer::CylinderLayer, centre::SVector{2, Float64}, p0::SVector{2, Float64}, p1::SVector{2, Float64}, dt::Float64)
    R = layer.radius
    total = 0.0
    if !iszero(layer.inner.rho)
        total += cylinder_surface_rate(layer.inner, R, true) * dt * shell_fraction(p0, p1, centre, R - layer.inner.h, R)
    end
    if !iszero(layer.outer.rho)
        total += cylinder_surface_rate(layer.outer, R, false) * dt * shell_fraction(p0, p1, centre, R, R + layer.outer.h)
    end
    return total
end
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders"])'`
Expected: PASS (all tests in "shell maths").

- [ ] **Step 5: Commit**

```bash
git add src/geometries/internal/layers.jl test/test_layer_cylinders.jl test/runtests.jl
git commit -m "Cylinder layer maths: exact shell fraction from line–circle roots, κ normalisation (P3.2.2)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Cylinder layer fields, resolution and validation (P3.2.4)

**Files:**
- Modify: `src/geometries/user/obstructions/obstructions.jl` (Cylinder fields)
- Modify: `src/geometries/user/fix.jl` (`obstruction_shifts`, `cylinder_layers`, `check_outer_reach`, cylinder `fix_type`)
- Modify: `src/geometries/internal/internal.jl` (import `CylinderLayer`)
- Test: `test/test_layer_cylinders.jl`, `test/test_cli.jl`

**Interfaces:**
- Consumes: `Internal.CylinderLayer`, `Internal.LayerSide` (Task 1).
- Produces:
  - user fields `layer_rho`, `layer_h`, `layer_rho_inner`, `layer_rho_outer`, `layer_h_inner`, `layer_h_outer`, `layer_shape`, `layer_shape_inner`, `layer_shape_outer` on `Cylinders`;
  - the fixed group's `layer::Union{Nothing, Vector{CylinderLayer}}`;
  - `obstruction_shifts(user_obstructions)::Vector`.

- [ ] **Step 1: Write the failing tests**

Append to `test/test_layer_cylinders.jl`:

```julia
@testset "test_layer_cylinders.jl: cylinder layer parameters" begin
    geom(c) = mr.Simulation([]; geometry=c, verbose=false).geometry

    @testset "no layer gives nothing" begin
        @test geom(mr.Cylinders(radius=1.))[1].layer === nothing
        @test geom(mr.Cylinders(radius=1., layer_h=0.2))[1].layer === nothing       # ρ = 0
    end

    @testset "same layer on both sides" begin
        l = geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2))[1].layer[1]
        @test l.centre == SVector(0., 0.)
        @test l.radius == 1.0
        @test (l.inner.rho, l.inner.h, l.inner.shape) == (0.02, 0.2, :step)
        @test (l.outer.rho, l.outer.h, l.outer.shape) == (0.02, 0.2, :step)
    end

    @testset "side-specific values override both-sides values" begin
        l = geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2, layer_rho_outer=0.05, layer_h_inner=0.4))[1].layer[1]
        @test (l.inner.rho, l.inner.h) == (0.02, 0.4)
        @test (l.outer.rho, l.outer.h) == (0.05, 0.2)
        l = geom(mr.Cylinders(radius=1., layer_rho_inner=0.02, layer_h_inner=0.2))[1].layer[1]
        @test l.outer.rho == 0.0
    end

    @testset "per-cylinder values and centres" begin
        layers = geom(mr.Cylinders(radius=[1., 0.5], position=[[0., 0.], [5., 0.]], layer_rho=[0.02, 0.], layer_h=0.2))[1].layer
        @test [l.centre for l in layers] == [SVector(0., 0.), SVector(5., 0.)]
        @test [l.radius for l in layers] == [1.0, 0.5]
        @test layers[2].inner.rho == 0.0
    end

    @testset "centres are wrapped into the unit cell as the obstructions are" begin
        l = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], position=[2.0, 0.0], layer_rho=0.02, layer_h=0.2))[1].layer[1]
        @test l.centre ≈ SVector(-0.5, 0.0)
    end

    @testset "h_inner = R is allowed (uniform inner rate)" begin
        l = geom(mr.Cylinders(radius=1., layer_rho_inner=0.05, layer_h_inner=1.0))[1].layer[1]
        @test l.inner.h == 1.0
    end

    @testset "validation" begin
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=-0.1, layer_h=0.2))
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.1, layer_h=-0.2))
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.1))                              # h = 0
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.1, layer_h=0.1, layer_h_inner=1.5))  # h_inner > R
        @test_throws ErrorException geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.1, layer_h=0.6))  # h_outer > gap 0.5
        @test_throws ErrorException geom(mr.Cylinders(radius=1., position=[[0., 0.], [2.4, 0.]], layer_rho=0.1, layer_h=0.3))  # gap 0.4
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.1, layer_h=0.2, layer_shape="linear"))   # step only
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.1, layer_h=0.2, layer_shape="bogus"))
        # h_outer = gap exactly is allowed; a side with ρ = 0 is not checked against the gap
        @test geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.1, layer_h=0.5))[1].layer !== nothing
        @test geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho_inner=0.1, layer_h_inner=0.2, layer_h_outer=0.9))[1].layer !== nothing
    end

    @testset "warning when spins can stick to the surface" begin
        @test_logs (:warn, r"not applied while a spin is stuck") geom(
            mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2, density=1., dwell_time=1.))
    end

    @testset "JSON round trip keeps the layer fields" begin
        c = mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2, layer_h_outer=0.3)
        io = IOBuffer()
        mr.write_geometry(io, c)
        back = mr.read_geometry_json(String(take!(io)))
        @test back.layer_rho.value == 0.02
        @test back.layer_h.value == 0.2
        @test back.layer_h_outer.value == 0.3
        @test isnothing(back.layer_rho_inner.value)
        @test back.layer_shape.value == "step"
    end
end
```

In `test/test_cli.jl`, after the testset "create walls with an exponential layer", add:

```julia
        @testset "create cylinders with a layer" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create cylinders 1 test.json --radius 1 --repeats 2.5,2.5 --layer_rho 0.02 --layer_h 0.2 --layer_h_outer 0.3")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["layer_rho"] == 0.02
                @test result["layer_h_outer"] == 0.3
            end
        end
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders", "cli"])'`
Expected: FAIL. The constructor errors on the unknown keyword `layer_rho` for `Cylinders`, and the CLI reports an unknown option `--layer_rho`.

- [ ] **Step 3: Implement**

(a) In `src/geometries/user/obstructions/obstructions.jl`, extend the `:Cylinder` field list after `:lorentz_radius`:

```julia
            Field{Float64}(:layer_rho, "Near-surface layer: integrated relaxivity ρ (um/ms) per unit surface area, inside and outside the cylinder unless a side-specific value is set. Zero means no layer.", 0.),
            Field{Float64}(:layer_h, "Near-surface layer: thickness h (um) inside and outside the cylinder, unless a side-specific value is set.", 0.),
            Field{Float64}(:layer_rho_inner, "Near-surface layer: ρ (um/ms) inside the cylinder. Overrides `layer_rho`."),
            Field{Float64}(:layer_rho_outer, "Near-surface layer: ρ (um/ms) outside the cylinder. Overrides `layer_rho`."),
            Field{Float64}(:layer_h_inner, "Near-surface layer: h (um) inside the cylinder (at most the radius). Overrides `layer_h`."),
            Field{Float64}(:layer_h_outer, "Near-surface layer: h (um) outside the cylinder (at most the surface gap to the nearest cylinder). Overrides `layer_h`."),
            Field{String}(:layer_shape, "Near-surface layer: profile g inside and outside, unless a side-specific value is set. Only step is supported on cylinders so far.", "step"),
            Field{String}(:layer_shape_inner, "Near-surface layer: profile inside the cylinder. Overrides `layer_shape`."),
            Field{String}(:layer_shape_outer, "Near-surface layer: profile outside the cylinder. Overrides `layer_shape`."),
```

(b) In `src/geometries/internal/internal.jl`, change the Layers import line to:

```julia
import .Layers: LayerSide, WallLayer, CylinderLayer, layer_exponent, has_layer, max_layer_rate
```

(c) In `src/geometries/user/fix.jl`:

Refactor the shift computation out of `apply_properties`. Replace its first block with:

```julia
    # apply shifts
    if apply_shift && hasproperty(user_obstructions, :position)
        internal_obstructions = Internal.Shift.(internal_obstructions, obstruction_shifts(user_obstructions))
    end
```

and add above `apply_properties`:

```julia
"""
    obstruction_shifts(user_obstructions)

Position of every obstruction as applied by [`apply_properties`](@ref): for a repeating group wrapped into the unit
cell [−repeats/2, repeats/2).
"""
function obstruction_shifts(user_obstructions::ObstructionGroup)
    shifts = isglobal(user_obstructions.position) ? fill(user_obstructions.position.value, length(user_obstructions)) : user_obstructions.position.value
    if hasproperty(user_obstructions, :repeats) && ~isnothing(user_obstructions.repeats.value)
        shifts = [mod.(s .+ user_obstructions.repeats.value/2, user_obstructions.repeats.value) .- user_obstructions.repeats.value/2 for s in shifts]
    end
    return shifts
end
```

Replace the cylinder `fix_type` with:

```julia
function fix_type(cylinders::Cylinders, index::Int, original_index::Int; kwargs...)
    if isglobal(cylinders.radius)
        base_obstructions = fill(Internal.Cylinder(cylinders.radius.value), length(cylinders))
    else
        base_obstructions = [Internal.Cylinder(r) for r in cylinders.radius.value]
    end
    layer = cylinder_layers(cylinders; density=kwargs[:density])
    apply_properties(cylinders, base_obstructions, index, original_index; surface="surface", volume="inside", layer=layer, kwargs...)
end

"""
    cylinder_layers(cylinders; density=0.)

Resolves the near-surface layer parameters of every cylinder into `Internal.CylinderLayer` objects (D3.1–D3.4):
`inner` is the layer inside the cylinder, `outer` the one outside. Side-specific fields (`layer_rho_inner`, ...)
override the both-sides fields. Only the step profile is supported on cylinders so far.

Returns `nothing` if no cylinder has a layer, so that simulations without a layer run the unmodified code.
"""
function cylinder_layers(cylinders::Cylinders; density=0.)
    function resolve(quantity::String, side::String, i::Int)
        specific = getproperty(cylinders, Symbol("layer_" * quantity * "_" * side))[i]
        isnothing(specific) || return Float64(specific)
        both = getproperty(cylinders, Symbol("layer_" * quantity))[i]
        return isnothing(both) ? 0.0 : Float64(both)
    end
    function resolve_shape(side::String, i::Int)
        specific = getproperty(cylinders, Symbol("layer_shape_" * side))[i]
        shape = isnothing(specific) ? cylinders.layer_shape[i] : specific
        shape = isnothing(shape) ? "step" : shape
        shape == "step" || error("Cylinder $i, $side side: only layer_shape = \"step\" is supported on cylinders so far, got \"$shape\".")
        return :step
    end
    radii = value_as_vector(cylinders.radius)
    centres = obstruction_shifts(cylinders)
    layers = Internal.CylinderLayer[]
    for i in 1:length(cylinders)
        R = Float64(isglobal(cylinders.radius) ? cylinders.radius.value : radii[i])
        (inner, outer) = map(("inner", "outer")) do side
            rho = resolve("rho", side, i)
            h = resolve("h", side, i)
            rho < 0 && error("Cylinder $i, $side side: layer_rho must be >= 0, got $rho.")
            h < 0 && error("Cylinder $i, $side side: layer_h must be >= 0, got $h.")
            rho > 0 && iszero(h) && error("Cylinder $i, $side side: layer_h must be > 0 when layer_rho > 0.")
            Internal.LayerSide(rho, h, resolve_shape(side, i))
        end
        inner.rho > 0 && inner.h > R && error("Cylinder $i: layer_h_inner = $(inner.h) um exceeds the radius ($R um).")
        push!(layers, Internal.CylinderLayer(SVector{2, Float64}(centres[i]), R, inner, outer))
    end
    if all(l -> iszero(l.inner.rho) && iszero(l.outer.rho), layers)
        return nothing
    end
    check_outer_reach(layers, cylinders.repeats.value)
    stuck_possible = density > 0 || any(d -> !isnothing(d) && d > 0, value_as_vector(cylinders.density))
    stuck_possible && @warn "The near-surface layer is not applied while a spin is stuck to a cylinder (surface density > 0)."
    return layers
end

"""
    check_outer_reach(layers, repeats)

Errors if an outer layer reaches inside another cylinder or an image of one: layer_h_outer must not exceed the surface
gap |c_i − c_j| − R_i − R_j to any neighbour (the wall rule h ≤ spacing). Centres are wrapped into the unit cell,
so the nearest images are within one cell.
"""
function check_outer_reach(layers::Vector{Internal.CylinderLayer}, repeats)
    shifts = isnothing(repeats) ? [SVector(0., 0.)] : [SVector(i * repeats[1], j * repeats[2]) for i in -1:1 for j in -1:1]
    for (i, li) in enumerate(layers)
        iszero(li.outer.rho) && continue
        for (j, lj) in enumerate(layers), shift in shifts
            (i == j && iszero(shift)) && continue
            v = lj.centre + shift - li.centre
            gap = sqrt(v[1]^2 + v[2]^2) - li.radius - lj.radius
            li.outer.h > gap && error("Cylinder $i: layer_h_outer = $(li.outer.h) um exceeds the surface gap ($gap um) to cylinder $j, so the outer layer would reach inside that cylinder.")
        end
    end
end
```

If `value_as_vector` or `SVector` is not in scope in `fix.jl`, check the imports at the top of the file. Walls use `value_as_vector` and `SVector` is used in the sphere code, so both should be there; add the import if a test reports `UndefVarError`.

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders", "cli", "layer", "collisions"])'`
Expected: PASS. The walls layer tests and the collision tests still pass after the `obstruction_shifts` refactor.

- [ ] **Step 5: Commit**

```bash
git add src/geometries/user/obstructions/obstructions.jl src/geometries/user/fix.jl src/geometries/internal/internal.jl test/test_layer_cylinders.jl test/test_cli.jl
git commit -m "Cylinders: layer fields (inner/outer), resolution, validation h_inner ≤ R and h_outer ≤ gap (P3.2.4)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Layer exponent on cylinder groups, images, overlap sum and timestep (P3.2.1, P3.2.3, P3.2.6)

**Files:**
- Modify: `src/geometries/internal/layers.jl`
- Test: `test/test_layer_cylinders.jl`

**Interfaces:**
- Consumes: `CylinderLayer`, `cylinder_exponent`, `cylinder_surface_rate` (Task 1); `group.layer::Vector{CylinderLayer}` (Task 2).
- Produces:
  - `image_centres(g, centre, lo, hi, reach)` (iterator of `SVector{2,Float64}`);
  - `layer_exponent(::Vector{CylinderLayer}, g, start, dest, dt)`;
  - `max_layer_rate(g::FixedObstructionGroup)`, which for cylinders is an upper bound on the summed rate.

- [ ] **Step 1: Write the failing tests**

Append to `test/test_layer_cylinders.jl`:

```julia
@testset "test_layer_cylinders.jl: layer exponent on cylinder groups" begin
    Layers = mr.Geometries.Internal.Layers
    geom(c) = mr.Simulation([]; geometry=c, verbose=false).geometry
    V3(x, y, z=0.) = SVector{3, Float64}(x, y, z)
    rin(ρ, h, R) = ρ / (h * (1 - h / (2R)))           # D3.1, written independently of the package
    rout(ρ, h, R) = ρ / (h * (1 + h / (2R)))

    "Brute-force summed rate at (x, y) over cylinders (centres, radii), images within ±2 cells if repeating."
    function brute_rate(x, y, centres, radii, ρ, h_in, h_out; repeats=nothing)
        cells = isnothing(repeats) ? [(0., 0.)] : [(i * repeats[1], j * repeats[2]) for i in -2:2 for j in -2:2]
        total = 0.0
        for (c, R) in zip(centres, radii), (ox, oy) in cells
            r = hypot(x - c[1] - ox, y - c[2] - oy)
            (R - h_in <= r <= R) && (total += rin(ρ, h_in, R))
            (R <= r <= R + h_out) && (total += rout(ρ, h_out, R))
        end
        return total
    end
    function brute_exponent(p0, p1, dt, args...; kw...)
        K = 20_000
        return dt * sum(brute_rate((p0 + (k - 0.5) / K * (p1 - p0))[1:2]..., args...; kw...) for k in 1:K) / K
    end
    function max_error(g, centres, radii, ρ, h_in, h_out, box; repeats=nothing, n=200)
        Random.seed!(77)
        err = 0.0
        for _ in 1:n
            p0 = V3((box .* (2 .* rand(2) .- 1))..., rand())
            p1 = V3((box .* (2 .* rand(2) .- 1))..., rand())
            err = max(err, abs(Layers.layer_exponent(g, p0, p1, 1.0) - brute_exponent(p0, p1, 1.0, centres, radii, ρ, h_in, h_out; repeats=repeats)))
        end
        return err
    end

    @testset "single cylinder: stationary points (rate) and a radial piece" begin
        g = geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h_inner=0.2, layer_h_outer=0.3))
        @test Layers.layer_exponent(g, V3(0.9, 0.), V3(0.9, 0.), 1.0) ≈ rin(0.02, 0.2, 1.0)
        @test Layers.layer_exponent(g, V3(1.2, 0.), V3(1.2, 0.), 1.0) ≈ rout(0.02, 0.3, 1.0)
        @test Layers.layer_exponent(g, V3(0.5, 0.), V3(0.5, 0.), 1.0) == 0.0
        @test Layers.layer_exponent(g, V3(1.5, 0.), V3(1.5, 0.), 1.0) == 0.0
        # axial motion does not change r
        @test Layers.layer_exponent(g, V3(0.9, 0., 0.), V3(0.9, 0., 5.), 2.0) ≈ 2 * rin(0.02, 0.2, 1.0)
    end

    @testset "Review Focus 3: rotated axis (cylinder along x)" begin
        g = geom(mr.Cylinders(radius=1., rotation=:x, layer_rho_inner=0.02, layer_h_inner=0.2))
        @test Layers.layer_exponent(g, V3(5., 0.9, 0.), V3(5., 0.9, 0.), 1.0) ≈ rin(0.02, 0.2, 1.0)
        @test Layers.layer_exponent(g, V3(0., 0., 0.9), V3(0., 0., 0.9), 1.0) ≈ rin(0.02, 0.2, 1.0)
        @test Layers.layer_exponent(g, V3(0.9, 0., 0.), V3(0.9, 0., 0.), 1.0) == 0.0     # on the axis
        @test Layers.layer_exponent(g, V3(-3., 0.9, 0.), V3(3., 0.9, 0.), 1.0) ≈ rin(0.02, 0.2, 1.0)  # along the axis
    end

    @testset "two overlapping outer layers (non-repeating) vs brute force" begin
        centres, radii = [(0., 0.), (2.2, 0.)], [1.0, 0.8]          # gap 0.4, outer h 0.3 → overlap 0.2 wide
        g = geom(mr.Cylinders(radius=radii, position=[[0., 0.], [2.2, 0.]], layer_rho=0.03, layer_h_inner=0.2, layer_h_outer=0.3))
        @test Layers.layer_exponent(g, V3(1.15, 0.), V3(1.15, 0.), 1.0) ≈ rout(0.03, 0.3, 1.0) + rout(0.03, 0.3, 0.8)
        @test max_error(g, centres, radii, 0.03, 0.2, 0.3, (2.5, 1.5)) < 1e-4
    end

    @testset "four overlapping cylinders (non-repeating) vs brute force" begin
        centres = [(1.2, 1.2), (-1.2, 1.2), (1.2, -1.2), (-1.2, -1.2)]
        g = geom(mr.Cylinders(radius=1., position=[[c...] for c in centres], layer_rho=0.03, layer_h=0.3))
        @test max_error(g, centres, fill(1.0, 4), 0.03, 0.3, 0.3, (2.5, 2.5)) < 1e-4
    end

    @testset "square array (images) vs brute force" begin
        g = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.03, layer_h_inner=0.2, layer_h_outer=0.3))
        @test max_error(g, [(0., 0.)], [1.0], 0.03, 0.2, 0.3, (4.0, 4.0); repeats=(2.5, 2.5)) < 1e-4
    end

    @testset "Review Focus 4: two cylinders per cell, rectangular array vs brute force" begin
        centres = [(0., 0.), (1.5, 1.25)]
        g = geom(mr.Cylinders(radius=[0.6, 0.5], position=[[c...] for c in centres], repeats=[3.0, 2.5],
            layer_rho=0.03, layer_h_inner=0.2, layer_h_outer=0.25))
        @test max_error(g, centres, [0.6, 0.5], 0.03, 0.2, 0.25, (4.0, 4.0); repeats=(3.0, 2.5)) < 1e-4
    end

    @testset "timestep: max_layer_rate is the summed-rate bound" begin
        # separated (gap 0.5 > 2h = 0.4): the largest single rate, inside
        g = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2))
        @test Layers.max_layer_rate(g) ≈ rin(0.02, 0.2, 1.0)
        # overlapping outer layers (2h = 0.6 > gap 0.5): own + 4 nearest neighbours (diagonals at 3.54 > 2.6 do not reach)
        g = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.1, layer_h=0.3))
        @test Layers.max_layer_rate(g) ≈ 5 * rout(0.1, 0.3, 1.0)
        s = mr.Simulation([]; geometry=mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.1, layer_h=0.3), diffusivity=3., verbose=false)
        @test s.timestep.max_timestep ≈ 0.005 / (5 * rout(0.1, 0.3, 1.0))
        # single cylinder
        @test Layers.max_layer_rate(geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h_inner=0.2, layer_h_outer=0.1))) ≈ rin(0.02, 0.2, 1.0)
        @test Layers.max_layer_rate(geom(mr.Cylinders(radius=1.))) == 0.0
    end
end
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders"])'`
Expected: FAIL with `MethodError: no method matching layer_exponent(::Vector{CylinderLayer}, ...)` and `MethodError` in `max_layer_rate`.

- [ ] **Step 3: Implement**

In `src/geometries/internal/layers.jl`, append after `layer_exponent(layers::Vector{WallLayer}, ...)`:

```julia
"""
    image_centres(g, centre, lo, hi, reach)

Axis positions (local 2D) of the images of a cylinder at `centre` in group `g` that lie within `reach` of the box
[lo, hi] in each coordinate: lattice images for a repeating group, `centre` itself otherwise.
"""
function image_centres(g::FixedObstructionGroup, centre::SVector{2, Float64}, lo::SVector{2, Float64}, hi::SVector{2, Float64}, reach::Float64)
    rep = repeating(g)
    sx = rep ? g.repeats[1] : 0.0
    sy = rep ? g.repeats[2] : 0.0
    kx = rep ? (ceil(Int, (lo[1] - reach - centre[1]) / sx):floor(Int, (hi[1] + reach - centre[1]) / sx)) : (0:0)
    ky = rep ? (ceil(Int, (lo[2] - reach - centre[2]) / sy):floor(Int, (hi[2] + reach - centre[2]) / sy)) : (0:0)
    return (centre + SVector(i * sx, j * sy) for i in kx for j in ky)
end

function layer_exponent(layers::Vector{CylinderLayer}, g::FixedObstructionGroup, start::SVector{3, Float64}, dest::SVector{3, Float64}, dt::Float64)
    p0 = SVector{2, Float64}(rotate_from_global(g, start))
    p1 = SVector{2, Float64}(rotate_from_global(g, dest))
    lo = min.(p0, p1)
    hi = max.(p0, p1)
    total = 0.0
    for layer in layers
        for c in image_centres(g, layer.centre, lo, hi, layer.radius + layer.outer.h)
            total += cylinder_exponent(layer, c, p0, p1, dt)
        end
    end
    return total
end
```

Replace the three `max_layer_rate` methods at the end of the module with:

```julia
max_layer_rate(::Nothing) = 0.0
max_layer_rate(layers::Vector{WallLayer}) = maximum((surface_rate(s) for l in layers for s in (l.positive, l.negative)); init=0.0)
max_layer_rate(g::FixedObstructionGroup) = max_layer_rate(g.layer, g)
max_layer_rate(::Nothing, ::FixedObstructionGroup) = 0.0
max_layer_rate(layers::Vector{WallLayer}, ::FixedObstructionGroup) = max_layer_rate(layers)

"""
Upper bound on the summed ΔR2 for cylinders: for each cylinder, the larger of its inner surface rate and its outer
surface rate plus the outer surface rates of every other cylinder (or image) whose outer shell intersects its outer
shell, |c_i − c_j| < R_i + h_i + R_j + h_j. Exact when no outer shells overlap. Inner shells overlap nothing
(validation: h_outer ≤ surface gap).
"""
function max_layer_rate(layers::Vector{CylinderLayer}, g::FixedObstructionGroup)
    best = 0.0
    for (i, li) in enumerate(layers)
        best = max(best, cylinder_surface_rate(li.inner, li.radius, true))
        own = cylinder_surface_rate(li.outer, li.radius, false)
        iszero(own) && continue
        total = own
        for (j, lj) in enumerate(layers)
            reach = li.radius + li.outer.h + lj.radius + lj.outer.h
            for c in image_centres(g, lj.centre, li.centre, li.centre, reach)
                (j == i && c == li.centre) && continue
                norm(c - li.centre) < reach && (total += cylinder_surface_rate(lj.outer, lj.radius, false))
            end
        end
        best = max(best, total)
    end
    return best
end

max_layer_rate(geometry::FixedGeometry) = maximum((max_layer_rate(g) for g in geometry); init=0.0)
```

Update the `max_layer_rate` docstring above the methods: "Largest near-surface rate (1/ms) in the geometry; for cylinders an upper bound on the summed rate where outer layers overlap (P3.2.6)". Also update item 6 of the `Simulation(timestep=...)` docstring in `src/timesteps.jl` to "`layer` / ΔR2_max, where ΔR2_max is the largest near-surface layer rate in the geometry (for cylinders with overlapping outer layers, an upper bound on their summed rate)".

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders", "layer"])'`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/geometries/internal/layers.jl src/timesteps.jl test/test_layer_cylinders.jl
git commit -m "Cylinder layer exponent: images within reach, overlap sum, rotated axes; summed-rate timestep bound (P3.2.1, P3.2.3, P3.2.6)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Layer on cylinders during the simulation (package tests, P3.2.5 and exact cases)

**Files:**
- Test: `test/test_layer_cylinders.jl`
- (No source change expected: `apply_layer!` already dispatches on `layer_exponent`. If a test fails, use superpowers:systematic-debugging.)

**Interfaces:**
- Consumes: everything from Tasks 1–3, `build_sequence` (test helper in `test/runtests.jl`).
- Produces: nothing new. This task gates the research tasks.

- [ ] **Step 1: Write the tests**

Append to `test/test_layer_cylinders.jl`:

```julia
@testset "test_layer_cylinders.jl: layer applied during the simulation" begin
    per_spin(snap) = [s.orientations[1].transverse for s in snap.spins]
    inside_mask(snap; R=1.) = [hypot(s.position[1], s.position[2]) < R for s in snap.spins]
    rin(ρ, h, R) = ρ / (h * (1 - h / (2R)))
    rout(ρ, h, R) = ρ / (h * (1 + h / (2R)))
    box = mr.BoundingBox([-1.5, -1.5, -1.], [1.5, 1.5, 1.])

    @testset "radial step with a reflection, inside and outside" begin
        cyl = mr.Cylinders(radius=1., layer_rho=0.02, layer_h_inner=0.2, layer_h_outer=0.3)
        seq = build_sequence([1.2, :readout])
        sim = mr.Simulation(seq; geometry=cyl, diffusivity=3., verbose=false)
        part = first(p for p in mr.parts([seq], 0., mr.TimeStep(1.2, Inf)) if p.duration > 0)
        # inside: 0.5 → 1.3 reflects at 1 and ends at 0.7; path 0.8 um in 1.2 ms; in [0.8, 1]: 0.2 + 0.2 um → 0.6 ms
        spin = mr.Spin(position=[0.5, 0., 0.], transverse=1., longitudinal=0.)
        mr.Evolve.draw_step!(spin, sim, part, [3.], [1.3, 0., 0.])
        @test spin.position[1] ≈ 0.7 atol=1e-12
        @test spin.orientations[1].transverse ≈ exp(-rin(0.02, 0.2, 1.) * 0.6) rtol=1e-12
        # outside: 1.5 → 0.7 reflects at 1 and ends at 1.3; in [1, 1.3]: 0.3 + 0.3 um → 0.9 ms
        spin = mr.Spin(position=[1.5, 0., 0.], transverse=1., longitudinal=0.)
        mr.Evolve.draw_step!(spin, sim, part, [3.], [0.7, 0., 0.])
        @test spin.position[1] ≈ 1.3 atol=1e-12
        @test spin.orientations[1].transverse ≈ exp(-rout(0.02, 0.3, 1.) * 0.9) rtol=1e-12
    end

    @testset "ρ = 0: bit-identical to no layer (P3.2.5)" begin
        function run(c)
            sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=c, diffusivity=3., R2=1/80, verbose=false)
            Random.seed!(7)
            snap = mr.Snapshot(2000, sim, 500; transverse=1., longitudinal=0.)
            return mr.readout(snap, sim, [50.]; return_snapshot=true)[1]
        end
        a = run(mr.Cylinders(radius=1., repeats=[2.5, 2.5]))
        b = run(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_h=0.3))
        @test mr.position.(a) == mr.position.(b)
        @test per_spin(a) == per_spin(b)
    end

    @testset "h_inner = R: every inside spin decays at exactly 2ρ/R; outside spins untouched" begin
        cyl = mr.Cylinders(radius=1., layer_rho_inner=0.05, layer_h_inner=1.0)       # κ = 1/2 → ΔR2 = 0.1 /ms
        sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=cyl, diffusivity=3., verbose=false)
        Random.seed!(1)
        snap = mr.Snapshot(2000, sim, box; transverse=1., longitudinal=0.)
        mask = inside_mask(snap)
        res = mr.readout(snap, sim, [10., 50.]; return_snapshot=true)
        for (t, s) in zip((10., 50.), res)
            m = per_spin(s)
            @test all(isapprox.(m[mask], exp(-0.1 * t); rtol=1e-10))
            @test all(m[.!mask] .== 1.0)
        end
    end

    @testset "bulk R2 and the layer multiply exactly" begin
        cyl = mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2)
        function run(R2)
            sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=cyl, diffusivity=3., R2=R2, verbose=false)
            Random.seed!(3)
            snap = mr.Snapshot(2000, sim, 500; transverse=1., longitudinal=0.)
            return mr.readout(snap, sim, [50.]; return_snapshot=true)[1]
        end
        a = run(0.)
        b = run(1/80)
        @test mr.position.(a) == mr.position.(b)
        @test all(isapprox.(per_spin(b), per_spin(a) .* exp(-50 / 80); rtol=1e-10))
        @test minimum(per_spin(a)) < 1.0
    end

    @testset "side-specific: each side acts only in its own compartment (P3.4.7)" begin
        function final(c)
            sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=c, diffusivity=3., verbose=false)
            Random.seed!(5)
            snap = mr.Snapshot(2000, sim, box; transverse=1., longitudinal=0.)
            mask = inside_mask(snap)
            return (per_spin(mr.readout(snap, sim, [50.]; return_snapshot=true)[1]), mask)
        end
        (m, mask) = final(mr.Cylinders(radius=1., layer_rho_inner=0.02, layer_h_inner=0.2))
        @test all(m[.!mask] .== 1.0)
        @test minimum(m[mask]) < 1.0
        (m, mask) = final(mr.Cylinders(radius=1., layer_rho_outer=0.02, layer_h_outer=0.2))
        @test all(m[mask] .== 1.0)
        @test minimum(m[.!mask]) < 1.0
    end

    @testset "Review Focus 5: permeable cylinder with a layer runs, M⊥ finite and ≤ 1" begin
        cyl = mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2, permeability=1.)
        sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=cyl, diffusivity=3., verbose=false)
        Random.seed!(9)
        snap = mr.Snapshot(1000, sim, box; transverse=1., longitudinal=0.)
        m = per_spin(mr.readout(snap, sim, [20.]; return_snapshot=true)[1])
        @test all(isfinite.(m))
        @test all(0 .< m .<= 1.0)
        @test minimum(m) < 1.0
    end

    @testset "finite RF pulse with a cylinder layer is rejected" begin
        seq = mr.SequenceParts.SequenceWaveform(
            (([], []), ([], []), ([], [])),
            [(0., 9., [mr.SequenceParts.ConstantPulse(0.25 / 9, 0., 0.)])],
            [], [10.], 10.,
        )
        sim = mr.Simulation(seq; geometry=mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2), verbose=false)
        @test_throws Exception mr.readout(100, sim, [10.])
    end
end
```

- [ ] **Step 2: Run them**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders"])'`
Expected: PASS. These tests exercise code from Tasks 1–3 end to end, so they are expected to pass on first run. If one fails, treat it as a bug in Tasks 1–3 (superpowers:systematic-debugging), not as a test to adjust.

Ledger ruling: these are integration tests written after the units they exercise. Their RED is the corresponding unit tests of Tasks 1–3. Two tests check that they can fail:
- multiply `cylinder_surface_rate` by 1.01 temporarily and confirm "radial step" and "h_inner = R" fail;
- then revert.

- [ ] **Step 3: Run the full suite without plots**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["no-plots"])'`
Expected: PASS, all test sets.

- [ ] **Step 4: Commit**

```bash
git add test/test_layer_cylinders.jl
git commit -m "Cylinder layer in the simulation: reflection step, V0b, uniform h = R, additivity, side-specific, permeable, RF rejected

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Phase 3 research harness

**Files:**
- Create: `research/phase3/harness.jl`, `research/phase3/test/runtests.jl`, `research/phase3/test/test_harness.jl`

**Interfaces:**
- Consumes:
  - `research/baseline/common.jl`: `REPO_ROOT`, `mr`, `Simulation`, `Snapshot`, `readout`, `Cylinders`, `write_json`, `write_csv`, `provenance`;
  - the package from Tasks 1–4.
- Produces:
  - constants `R_REF = 1.0`, `S_REF = 2.5`, `D_REF = 3.0`, `TAU_REF = 1e-2`, `TIMES_REF = 0:5:50`, `SEEDS_REF = 1:10`, `NSPINS_REF = 100_000`;
  - `phase3_outdir(name)`;
  - `sv_inside(R)`, `sv_outside(R, sx, sy)`, `shell_area_inside(R, h)`, `shell_area_outside(R, h, sx, sy)`;
  - `rate_inside(ρ, h, R)`, `rate_outside(ρ, h, R)`;
  - `radial_distance(pos; sx=Inf, sy=Inf)`;
  - `brute_rate(x, y; R, s, h_in, h_out, rate_in, rate_out)`;
  - `roundoff_bound(times, tau, nspins)`;
  - `expected_neglog_conventional(θ, D, sv, t)`, `expected_neglog_layer(ρ, sv, t)`;
  - `run_cylinders(geometry; region, R, sx, sy, bounding_box, R2_bulk, timestep, times, seeds, nspins, D)`, which returns the named tuple `(S, mean, sigma, sem, neglog, neglog_mean, neglog_sem, times, tau, runtime, nregion, masks, finals, sim)`.

- [ ] **Step 1: Write the failing tests**

`research/phase3/test/runtests.jl`:

```julia
using Test
include(joinpath(@__DIR__, "..", "harness.jl"))

@testset "Phase 3 research code" begin
    for f in sort(filter(f -> startswith(f, "test_") && endswith(f, ".jl"), readdir(@__DIR__)))
        include(joinpath(@__DIR__, f))
    end
end
```

`research/phase3/test/test_harness.jl`:

```julia
@testset "harness: geometry formulas" begin
    @test sv_inside(1.0) == 2.0
    @test sv_outside(1.0, 2.5, 2.5) ≈ 2π / (6.25 - π)
    @test shell_area_inside(1.0, 1.0) == 1.0
    @test shell_area_inside(1.0, 0.2) ≈ 0.36
    @test shell_area_outside(1.0, 0.2, 2.5, 2.5) ≈ π * 0.44 / (6.25 - π)
    @test rate_inside(0.02, 0.2, 1.0) ≈ 0.02 / 0.18
    @test rate_outside(0.02, 0.2, 1.0) ≈ 0.02 / 0.22
    @test radial_distance([3.0, 0.2, 7.0]; sx=2.5, sy=2.5) ≈ hypot(0.5, 0.2)
    @test radial_distance([3.0, 4.0, 0.0]) == 5.0
    @test roundoff_bound([0.0, 50.0], 1e-2, 100_000) ≈ (5000 + 100_000) * eps()
    @test expected_neglog_layer(0.02, 2.0, 50.0) ≈ 2.0
    @test expected_neglog_conventional(0.01, 3.0, 2.0, 50.0) ≈ 0.01 * sqrt(3 / π) * 100
end

@testset "harness: brute_rate" begin
    kw = (R=1.0, s=2.5, h_in=0.2, h_out=0.3, rate_in=1.0, rate_out=10.0)
    @test brute_rate(0.9, 0.0; kw...) == 1.0
    @test brute_rate(1.1, 0.0; kw...) == 10.0
    @test brute_rate(1.25, 0.0; kw...) == 20.0        # in the outer shells of two neighbours
    @test brute_rate(0.5, 0.0; kw...) == 0.0
    @test brute_rate(1.25 + 2.5 * 7, 0.0; kw...) == 20.0
end

@testset "harness: run_cylinders" begin
    R = 1.0
    bb = mr.BoundingBox([-R, -R, -R], [R, R, R])
    r = run_cylinders(Cylinders(radius=R); region=:inside, bounding_box=bb, times=[0.0, 1.0], seeds=1:2, nspins=4000)
    @test size(r.S) == (2, 2)
    @test all(r.S .== 1.0)
    @test all(r.neglog .== 0.0)
    @test all(abs.(r.nregion ./ 4000 .- π / 4) .< 5 * sqrt(π / 4 * (1 - π / 4) / 4000))
    # uniform inner layer h = R: −ln M⊥ = 0.1·t exactly for every inside spin
    u = run_cylinders(Cylinders(radius=R, layer_rho_inner=0.05, layer_h_inner=R); region=:inside, bounding_box=bb,
        times=[0.0, 5.0], seeds=1:2, nspins=2000)
    @test all(isapprox.(u.neglog[:, 2], 0.5; rtol=1e-9))
    o = run_cylinders(Cylinders(radius=R, repeats=[S_REF, S_REF]); region=:outside, sx=S_REF, sy=S_REF,
        times=[0.0, 1.0], seeds=1:1, nspins=4000)
    @test abs(o.nregion[1] / 4000 - (1 - π / S_REF^2)) < 5 * sqrt(0.5 * 0.5 / 4000)
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`
Expected: FAIL with `could not open file .../research/phase3/harness.jl`.

- [ ] **Step 3: Implement `research/phase3/harness.jl`**

```julia
# Phase 3 research harness: near-surface layer on cylinders (single and arrays), no RF.
# Settings, formulas and decisions: research/phase3_record.md (D3.1–D3.6). Baselines (MCMR θ) are references, not truth.

include(joinpath(@__DIR__, "..", "baseline", "common.jl"))
using Statistics, Random

const R_REF = 1.0                       # cylinder radius (um)
const S_REF = 2.5                       # square-array centre spacing (um): surface gap 0.5 um
const D_REF = 3.0                       # um²/ms
const TAU_REF = 1e-2                    # ms
const TIMES_REF = collect(0.0:5.0:50.0)
const SEEDS_REF = 1:10
const NSPINS_REF = 100_000

phase3_outdir(name::AbstractString) = mkpath(joinpath(REPO_ROOT, "research", "results", "phase3", name))

"Surface-to-volume ratio inside a cylinder: 2/R."
sv_inside(R) = 2 / R
"Surface-to-volume ratio outside, one cylinder per sx × sy cell: 2πR/(sx·sy − πR²)."
sv_outside(R, sx, sy) = 2π * R / (sx * sy - π * R^2)
"Area fraction of the inside within h of the surface: (R² − (R − h)²)/R²."
shell_area_inside(R, h) = (R^2 - (R - h)^2) / R^2
"Area fraction of the outside within h of the surface (shells not overlapping): π((R + h)² − R²)/(sx·sy − πR²)."
shell_area_outside(R, h, sx, sy) = π * ((R + h)^2 - R^2) / (sx * sy - π * R^2)
"Step layer rate inside (D3.1): ρ/(h(1 − h/2R)); written independently of the package."
rate_inside(ρ, h, R) = ρ / (h * (1 - h / (2R)))
"Step layer rate outside (D3.1): ρ/(h(1 + h/2R))."
rate_outside(ρ, h, R) = ρ / (h * (1 + h / (2R)))

"Distance from the nearest cylinder axis (one cylinder per cell at the origin; sx = sy = Inf for a single cylinder)."
function radial_distance(pos; sx=Inf, sy=Inf)
    x = isfinite(sx) ? mod(pos[1] + sx / 2, sx) - sx / 2 : pos[1]
    y = isfinite(sy) ? mod(pos[2] + sy / 2, sy) - sy / 2 : pos[2]
    return hypot(x, y)
end

"""
    brute_rate(x, y; R, s, h_in, h_out, rate_in, rate_out)

Summed step-layer ΔR2 at (x, y) for a square array of spacing `s` (one cylinder per cell at the origin; `s = Inf` for
a single cylinder), from the 5 × 5 nearest images. Independent of the package's `Layers` code.
"""
function brute_rate(x, y; R, s, h_in, h_out, rate_in, rate_out)
    cells = isfinite(s) ? (-2:2) : (0:0)
    cx = isfinite(s) ? round(x / s) * s : 0.0
    cy = isfinite(s) ? round(y / s) * s : 0.0
    total = 0.0
    for i in cells, j in cells
        r = hypot(x - cx - i * (isfinite(s) ? s : 0.0), y - cy - j * (isfinite(s) ? s : 0.0))
        (R - h_in <= r <= R) && (total += rate_in)
        (R <= r <= R + h_out) && (total += rate_out)
    end
    return total
end

"Floating-point bound for an ensemble mean after maximum(times)/τ multiplications per spin and a sum over nspins."
roundoff_bound(times, tau, nspins) = (maximum(times) / tau + nspins) * eps()

"Conventional model (MCMR θ): E[−ln M⊥] = θ√(D/π)·(S/V)·t, if the wall-collision rate holds on curved surfaces (P3.1.6)."
expected_neglog_conventional(θ, D, sv, t) = θ * sqrt(D / π) * sv * t
"Layer, ρ per unit surface area (D3.1): E[−ln M⊥] = ρ·(S/V)·t (mean-exponent identity, P3.4.3)."
expected_neglog_layer(ρ, sv, t) = ρ * sv * t

"""
    run_cylinders(geometry; region, R=R_REF, sx=Inf, sy=Inf, bounding_box=500, R2_bulk=0.0, timestep=TAU_REF,
                  times=TIMES_REF, seeds=SEEDS_REF, nspins=NSPINS_REF, D=D_REF)

Free decay (no RF) for every seed (`Random.seed!(seed)` before each snapshot, so seed k gives the same initial
positions in every run). `region` (:inside or :outside) selects spins by their distance from the nearest axis at t = 0;
surfaces are impermeable, so the set does not change. `timestep=nothing` uses the simulator default.
Returns `S[i, j]` (mean M⊥ of seed i's region spins at times[j]) and `neglog[i, j]` (mean −ln M⊥ over them), with
mean, sigma and SEM over seeds; `nregion[i]`, `masks[i]` and `finals[i]` (final Snapshot) per seed.
"""
function run_cylinders(geometry; region::Symbol, R=R_REF, sx=Inf, sy=Inf, bounding_box=500, R2_bulk=0.0,
                       timestep=TAU_REF, times=TIMES_REF, seeds=SEEDS_REF, nspins=NSPINS_REF, D=D_REF)
    region in (:inside, :outside) || error("region must be :inside or :outside, got $region")
    kw = isnothing(timestep) ? (;) : (; timestep=timestep)
    sim = Simulation([mr.SequenceParts.empty_sequence()]; geometry=geometry, diffusivity=D, R2=R2_bulk, verbose=false, kw...)
    times = collect(Float64, times)
    S = zeros(length(seeds), length(times))
    L = zeros(length(seeds), length(times))
    nregion = zeros(Int, length(seeds))
    masks = BitVector[]
    finals = Snapshot[]
    runtime = 0.0
    for (i, seed) in enumerate(seeds)
        Random.seed!(seed)
        snap = Snapshot(nspins, sim, bounding_box; transverse=1.0, longitudinal=0.0)
        r0 = [radial_distance(s.position; sx=sx, sy=sy) for s in snap.spins]
        mask = BitVector(region == :inside ? r0 .< R : r0 .> R)
        t0 = time()
        res = readout(snap, sim, times; return_snapshot=true)
        runtime += time() - t0
        for (j, r) in enumerate(res)
            m = [r.spins[k].orientations[1].transverse for k in findall(mask)]
            S[i, j] = mean(m)
            L[i, j] = mean(-log.(m))
        end
        nregion[i] = count(mask)
        push!(masks, mask)
        push!(finals, res[end])
    end
    n = length(seeds)
    sd(X) = n > 1 ? vec(std(X, dims=1)) : zeros(size(X, 2))
    return (S=S, mean=vec(mean(S, dims=1)), sigma=sd(S), sem=sd(S) ./ sqrt(n), neglog=L,
        neglog_mean=vec(mean(L, dims=1)), neglog_sem=sd(L) ./ sqrt(n), times=times,
        tau=sim.timestep.max_timestep, runtime=runtime, nregion=nregion, masks=masks, finals=finals, sim=sim)
end
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add research/phase3/harness.jl research/phase3/test/
git commit -m "Phase 3 harness: cylinder runs by region, S/V and shell formulas, brute-force rate, expected identities

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: C-B1 and C-B2 (P3.1.1, P3.1.2)

**Files:**
- Create: `research/phase3/p3_1_1_null_bulk.jl` → `research/results/phase3/p3_1_1_null_bulk/`

**Interfaces:**
- Consumes: `run_cylinders`, `roundoff_bound`, `phase3_outdir` (Task 5).
- Produces: `summary.json` with per-configuration `tau_default`, `max_dev_per_spin`, `max_dev_ensemble`, `pass`.

- [ ] **Step 1: Write the script**

```julia
# P3.1.1 C-B1 (null run) and P3.1.2 C-B2 (bulk R2 only) on cylinders, old code (no layer).
# Configurations: single cylinder R = 1 (inside), square array s = 2.5, R = 1 (inside, outside). Default timestep
# (recorded: C-B1 asks for τ_max), seeds 1–10, 100000 spins per seed, readouts 0:5:50.
# C-B1: every spin keeps M⊥ = 1 exactly. C-B2 (R2 = 1/80): positions identical to C-B1 (same seeds); per spin
# |M⊥/e^(−t/80) − 1| ≤ 1e-10; ensemble ≤ roundoff_bound.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_1_1_null_bulk.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_1_1_null_bulk")
const R2B = 1 / 80

configs = [
    ("single_inside", Cylinders(radius=R_REF), :inside, mr.BoundingBox([-R_REF, -R_REF, -R_REF], [R_REF, R_REF, R_REF]), Inf),
    ("array_inside", Cylinders(radius=R_REF, repeats=[S_REF, S_REF]), :inside, 500, S_REF),
    ("array_outside", Cylinders(radius=R_REF, repeats=[S_REF, S_REF]), :outside, 500, S_REF),
]
results = Dict{String, Any}()
allpass = true
for (name, geometry, region, bb, s) in configs
    kw = (region=region, bounding_box=bb, sx=s, sy=s, timestep=nothing)
    null = run_cylinders(geometry; kw...)
    bulk = run_cylinders(geometry; R2_bulk=R2B, kw...)
    null_exact = all(all(sp.orientations[1].transverse == 1.0 for sp in f.spins) for f in null.finals)
    same_pos = all(mr.position.(a) == mr.position.(b) for (a, b) in zip(null.finals, bulk.finals))
    target = exp(-R2B * TIMES_REF[end])
    spin_dev = maximum(maximum(abs(sp.orientations[1].transverse / target - 1) for (k, sp) in enumerate(f.spins) if m[k])
                       for (f, m) in zip(bulk.finals, bulk.masks))
    ens_dev = maximum(abs.(bulk.mean ./ exp.(-R2B .* TIMES_REF) .- 1))
    bound = roundoff_bound(TIMES_REF, bulk.tau, NSPINS_REF)
    pass = null_exact && all(null.S .== 1.0) && same_pos && spin_dev <= 1e-10 && ens_dev <= bound
    global allpass &= pass
    results[name] = Dict("tau_default" => null.tau, "null_exact" => null_exact, "positions_identical" => same_pos,
        "max_dev_per_spin" => spin_dev, "max_dev_ensemble" => ens_dev, "bound" => bound,
        "nregion_mean" => mean(null.nregion), "runtime_s" => null.runtime + bulk.runtime, "pass" => pass)
    @printf("%-14s τ = %.4g ms: C-B1 exact %s; C-B2 positions identical %s, per spin %.1e, ensemble %.1e (bound %.1e) → %s\n",
        name, null.tau, null_exact, same_pos, spin_dev, ens_dev, bound, pass ? "PASS" : "FAIL")
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "results" => results, "pass" => allpass,
    "params" => Dict("R" => R_REF, "s" => S_REF, "D" => D_REF, "R2_bulk" => R2B, "seeds" => collect(SEEDS_REF), "nspins" => NSPINS_REF)))
println("saved to $OUT")
```

- [ ] **Step 2: Run**

Run: `julia --project=research/baseline -t 8 research/phase3/p3_1_1_null_bulk.jl`
Expected: three PASS lines. τ_default is 0.01 ms (tortuosity 0.03·R²/D with size scale R = 1). Per-spin deviation ≈ 1e-13; ensemble ≤ 2.3e-11.

- [ ] **Step 3: Commit**

```bash
git add research/phase3/p3_1_1_null_bulk.jl research/results/phase3/p3_1_1_null_bulk
git commit -m "P3.1.1/P3.1.2: cylinder C-B1 null and C-B2 bulk baselines

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: C-B3 density and shell occupancy, without and with the layer (P3.1.3, P3.4.5)

**Files:**
- Create: `research/phase3/p3_1_3_density.jl` → `research/results/phase3/p3_1_3_density/` (`summary.json`, `occupancy.csv`, `p3_1_3_density.png`)

**Interfaces:**
- Consumes: `run_cylinders`, `radial_distance`, `shell_area_inside`, `shell_area_outside` (Task 5).
- Produces: record numbers (χ²/dof, z_χ, occupancy z per h, positions-identical flags).

- [ ] **Step 1: Write the script**

```julia
# P3.1.3 C-B3: density across the cylinder (inside and outside) and the fraction within h of the surface; P3.4.5: the
# same with a step layer on (positions must equal the no-layer run). Default timestep, seeds 1–10, 1e5 spins, t = 50 ms.
#
# Inside (single cylinder R = 1): r²/R² is uniform on [0, 1] for a uniform density → 100-bin χ² of r²/R²,
#   z_χ = (χ² − dof)/√(2 dof), dof = 99, |z_χ| ≤ 3; occupancy within h: p = (R² − (R − h)²)/R², z = (n − Np)/√(Np(1 − p)), |z| ≤ 2.
# Outside (square array s = 2.5): for R ≤ r ≤ s/2 whole annuli lie in the cell, so P(r ∈ [r_k, r_k+1] | r ≤ s/2) =
#   (r_k+1² − r_k²)/((s/2)² − R²) → 50-bin χ² on spins with r ≤ s/2; occupancy p = π((R + h)² − R²)/(s² − πR²).
# Layer on: ρ = 0.02 um/ms, h = 0.2 um both sides.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_1_3_density.jl

include(joinpath(@__DIR__, "harness.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_1_3_density")
const HS_IN = [0.05, 0.1, 0.2, 0.5, 1.0]
const HS_OUT = [0.05, 0.1, 0.2]

zchi(counts, p) = (N = sum(counts); e = N .* p; χ2 = sum((counts .- e) .^ 2 ./ e); dof = length(p) - 1; (χ2 / dof, (χ2 - dof) / sqrt(2dof)))
binom_z(n, N, p) = (n - N * p) / sqrt(N * p * (1 - p))

function radii_at_end(r; s)
    vcat([[radial_distance(sp.position; sx=s, sy=s) for (k, sp) in enumerate(f.spins) if m[k]] for (f, m) in zip(r.finals, r.masks)]...)
end

bb = mr.BoundingBox([-R_REF, -R_REF, -R_REF], [R_REF, R_REF, R_REF])
rows = NamedTuple[]
summary = Dict{String, Any}()
allpass = true
radii_store = Dict{String, Vector{Float64}}()
for (label, layer_kw) in (("no_layer", (;)), ("layer_on", (layer_rho=0.02, layer_h=0.2)))
    ins = run_cylinders(Cylinders(radius=R_REF; layer_kw...); region=:inside, bounding_box=bb, timestep=nothing, times=[0.0, 50.0])
    out = run_cylinders(Cylinders(radius=R_REF, repeats=[S_REF, S_REF]; layer_kw...); region=:outside, sx=S_REF, sy=S_REF,
        timestep=nothing, times=[0.0, 50.0])
    rin = radii_at_end(ins; s=Inf)
    rout = radii_at_end(out; s=S_REF)
    radii_store[label * "_in"] = rin
    radii_store[label * "_out"] = rout
    # inside histogram of r²/R²
    u = clamp.((rin ./ R_REF) .^ 2, 0, prevfloat(1.0))
    cin = [count(x -> (k - 1) / 100 <= x < k / 100, u) for k in 1:100]
    (chi_in, z_in) = zchi(cin, fill(1 / 100, 100))
    # outside histogram of r on [R, s/2]
    edges = sqrt.(range(R_REF^2, (S_REF / 2)^2, length=51))
    ro = filter(r -> r <= S_REF / 2, rout)
    cout = [count(r -> edges[k] <= r < edges[k + 1], ro) for k in 1:50]
    (chi_out, z_out) = zchi(cout, fill(1 / 50, 50))           # equal-area bins
    occ_ok = true
    for h in HS_IN
        p = shell_area_inside(R_REF, h)
        z = binom_z(count(r -> r >= R_REF - h, rin), length(rin), p)
        occ_ok &= abs(z) <= 2
        push!(rows, (case=label, region="inside", h=h, p=p, fraction=count(r -> r >= R_REF - h, rin) / length(rin), z=z))
    end
    for h in HS_OUT
        p = shell_area_outside(R_REF, h, S_REF, S_REF)
        n = count(r -> r <= R_REF + h, rout)
        z = binom_z(n, length(rout), p)
        occ_ok &= abs(z) <= 2
        push!(rows, (case=label, region="outside", h=h, p=p, fraction=n / length(rout), z=z))
    end
    pass = abs(z_in) <= 3 && abs(z_out) <= 3 && occ_ok
    summary[label] = Dict("inside_chi2_per_dof" => chi_in, "inside_z_chi2" => z_in, "outside_chi2_per_dof" => chi_out,
        "outside_z_chi2" => z_out, "occupancy_all_within_2sigma" => occ_ok, "n_inside" => length(rin), "n_outside" => length(rout),
        "tau" => ins.tau, "pass" => pass)
    global allpass &= pass
    @printf("%-8s: inside χ²/dof %.3f (z %+.2f), outside χ²/dof %.3f (z %+.2f), occupancy ok %s → %s\n",
        label, chi_in, z_in, chi_out, z_out, occ_ok, pass ? "PASS" : "FAIL")
    summary[label * "_runs"] = (ins, out)
end
# P3.4.5: positions with the layer on equal the no-layer run
(ni, no) = summary["no_layer_runs"]
(li, lo) = summary["layer_on_runs"]
same = all(mr.position.(a) == mr.position.(b) for (a, b) in zip(vcat(ni.finals, no.finals), vcat(li.finals, lo.finals)))
@printf("P3.4.5 positions identical with the layer on: %s\n", same)
delete!(summary, "no_layer_runs"); delete!(summary, "layer_on_runs")
summary["layer_on_positions_identical"] = same
write_csv(joinpath(OUT, "occupancy.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "results" => summary, "pass" => allpass && same))

fig = Figure(size=(1000, 380))
ax1 = Axis(fig[1, 1], xlabel="r²/R²", ylabel="density / uniform", title="inside, t = 50 ms (no layer)")
u = (radii_store["no_layer_in"] ./ R_REF) .^ 2
hist!(ax1, u, bins=50, normalization=:pdf)
hlines!(ax1, [1.0], color=:black)
ax2 = Axis(fig[1, 2], xlabel="r (um)", ylabel="density per unit area / uniform", title="outside, r ≤ s/2 (no layer)")
ro = filter(r -> r <= S_REF / 2, radii_store["no_layer_out"])
edges = range(R_REF, S_REF / 2, length=31)
c = [count(r -> edges[k] <= r < edges[k + 1], ro) for k in 1:30]
area = [edges[k + 1]^2 - edges[k]^2 for k in 1:30]
scatter!(ax2, (edges[1:end-1] .+ edges[2:end]) ./ 2, (c ./ area) ./ (length(ro) / ((S_REF / 2)^2 - R_REF^2)))
hlines!(ax2, [1.0], color=:black)
save(joinpath(OUT, "p3_1_3_density.png"), fig)
println("saved to $OUT")
```

- [ ] **Step 2: Run**

Run: `julia --project=research/baseline -t 8 research/phase3/p3_1_3_density.jl`
Expected: "no_layer … PASS", "layer_on … PASS", "P3.4.5 positions identical with the layer on: true".

If a χ² or occupancy z fails, check it is not a binning artefact: r exactly at the edges, or the r ≤ s/2 cut. Then ledger it as a finding and record it. Do not loosen the criteria. A real non-uniform density near a curved surface is a result.

- [ ] **Step 3: Commit**

```bash
git add research/phase3/p3_1_3_density.jl research/results/phase3/p3_1_3_density
git commit -m "P3.1.3/P3.4.5: cylinder density and shell occupancy, without and with the layer

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: C-B4 conventional θ reference curves (P3.1.4)

**Files:**
- Create: `research/phase3/p3_1_4_theta_reference.jl` → `research/results/phase3/p3_1_4_theta_reference/` (`reference_curves.csv`, `summary.json`)

**Interfaces:**
- Consumes: `run_cylinders` (Task 5).
- Produces: `reference_curves.csv` with columns `config, R, s, region, theta, tau, t, mean, sem`. These are the conventional reference curves for later comparisons.

- [ ] **Step 1: Write the script**

```julia
# P3.1.4 C-B4: conventional surface relaxation (MCMR θ) reference curves on cylinders. A reference, not the truth.
# Configurations: single cylinder inside R = 0.5, 1, 2; square array s = 2.5, R = 1 inside and outside.
# θ = 0.003, 0.01, 0.03; D = 3; readouts 0:5:50; seeds 1–10, NSPINS (default 20000) spins per seed.
# Timesteps τ_a = min(1e-2, 0.03R²/D) and τ_b = τ_a/4. Convergence: |S_a(50) − S_b(50)|/√(SEM_a² + SEM_b²) ≤ 2
# (one test per curve; the 11 readouts are correlated, so they are reported, not each graded).
# Bulk factorisation (θ = 0.01, single R = 1, τ_a): same seeds with R2 = 1/80 → per spin ratio e^(−t/80) within 1e-10.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_1_4_theta_reference.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_1_4_theta_reference")
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const THETAS = [0.003, 0.01, 0.03]

configs = [
    ("single_R0.5", 0.5, Inf, :inside), ("single_R1", 1.0, Inf, :inside), ("single_R2", 2.0, Inf, :inside),
    ("array_inside", R_REF, S_REF, :inside), ("array_outside", R_REF, S_REF, :outside),
]
geom(R, s, θ) = isfinite(s) ? Cylinders(radius=R, repeats=[s, s], surface_relaxation=θ) : Cylinders(radius=R, surface_relaxation=θ)
bbox(R, s) = isfinite(s) ? 500 : mr.BoundingBox([-R, -R, -R], [R, R, R])

rows = NamedTuple[]
conv = NamedTuple[]
for (name, R, s, region) in configs, θ in THETAS
    τa = min(1e-2, 0.03 * R^2 / D_REF)
    for τ in (τa, τa / 4)
        r = run_cylinders(geom(R, s, θ); region=region, R=R, sx=s, sy=s, bounding_box=bbox(R, s), timestep=τ, nspins=NSPINS)
        append!(rows, [(config=name, R=R, s=s, region=String(region), theta=θ, tau=τ, t=t, mean=m, sem=e)
                       for (t, m, e) in zip(r.times, r.mean, r.sem)])
    end
    a = filter(x -> x.config == name && x.theta == θ && x.tau == τa && x.t == 50.0, rows)[1]
    b = filter(x -> x.config == name && x.theta == θ && x.tau == τa / 4 && x.t == 50.0, rows)[1]
    z = (a.mean - b.mean) / sqrt(a.sem^2 + b.sem^2)
    push!(conv, (config=name, theta=θ, tau_a=τa, S50_a=a.mean, S50_b=b.mean, z=z, pass=abs(z) <= 2))
    @printf("%-14s θ = %.3f: S(50) %.5f (τ %.2e) vs %.5f (τ/4), z = %+.2f\n", name, θ, a.mean, τa, b.mean, z)
end
# bulk factorisation
g = geom(1.0, Inf, 0.01)
kw = (region=:inside, bounding_box=bbox(1.0, Inf), timestep=1e-2, nspins=NSPINS, seeds=1:3)
a = run_cylinders(g; kw...)
b = run_cylinders(g; R2_bulk=1 / 80, kw...)
fac = maximum(maximum(abs(sb.orientations[1].transverse / (sa.orientations[1].transverse * exp(-50 / 80)) - 1)
                      for (sa, sb) in zip(fa.spins, fb.spins)) for (fa, fb) in zip(a.finals, b.finals))
@printf("bulk factorisation: max per-spin deviation %.1e\n", fac)
write_csv(joinpath(OUT, "reference_curves.csv"), rows)
write_csv(joinpath(OUT, "convergence.csv"), conv)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_converged" => count(c -> c.pass, conv),
    "n_curves" => length(conv), "bulk_factorisation_max_dev" => fac, "nspins" => NSPINS))
println("saved to $OUT")
```

- [ ] **Step 2: Run**

Run: `julia --project=research/baseline -t 8 research/phase3/p3_1_4_theta_reference.jl`
Expected:
- 15 lines with |z| ≤ 2 for most curves;
- bulk factorisation ≤ 1e-10.

A curve with |z| > 2 is not a code failure. It is a timestep dependence of the conventional model on a curved surface. Rerun that curve with τ_a/16. If τ_a/4 and τ_a/16 agree, ledger "Ruling: reference uses τ_a/4 for <config>", and record the τ dependence.

- [ ] **Step 3: Commit**

```bash
git add research/phase3/p3_1_4_theta_reference.jl research/results/phase3/p3_1_4_theta_reference
git commit -m "P3.1.4: cylinder C-B4 conventional θ reference curves with τ convergence and bulk factorisation

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 9: Conventional mapping on curved surfaces (P3.1.6)

**Files:**
- Create: `research/phase3/p3_1_6_mapping.jl` → `research/results/phase3/p3_1_6_mapping/` (`mapping.csv`, `summary.json`, `p3_1_6_mapping.png`)

**Interfaces:**
- Consumes: `run_cylinders`, `expected_neglog_conventional`, `sv_inside`, `sv_outside` (Task 5).
- Produces: `mapping.csv` with columns `config, R, D, tau, sqrt_Dtau_over_R, expected, measured, sem, z, ratio`. A later plan uses the measured ratio as the mapping if it differs from 1.

- [ ] **Step 1: Write the script**

```julia
# P3.1.6: does the conventional model's relaxivity ρ = θ√(D/π) hold on curved surfaces? The per-collision factor is
# e^(−θ√τ), so −ln M⊥ = θ√τ·(number of collisions) for every path, and the test is the collision rate:
#     E[−ln M⊥(t)] = θ√(D/π)·(S/V)·t        (exact for walls, P2.5.0)
# On a cylinder the finite step √(2Dτ) can bias the collision rate by O(√(Dτ)/R). Measured at t = 50 ms for single
# cylinders (inside) R = 0.5, 1, 2 and the square array (inside, outside); D = 1, 3, 6; τ = 1e-2 and 2.5e-3 ms;
# θ = 0.01; seeds 1–10, NSPINS (default 20000). z = (mean − expected)/SEM; |z| ≤ 3 passes, otherwise the measured
# ratio is recorded as the mapping (tracker P3.1.6).
# Run: julia --project=research/baseline -t 8 research/phase3/p3_1_6_mapping.jl

include(joinpath(@__DIR__, "harness.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_1_6_mapping")
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const θ = 0.01
const T = 50.0

configs = [("single_R0.5", 0.5, Inf, :inside), ("single_R1", 1.0, Inf, :inside), ("single_R2", 2.0, Inf, :inside),
           ("array_inside", R_REF, S_REF, :inside), ("array_outside", R_REF, S_REF, :outside)]
rows = NamedTuple[]
for (name, R, s, region) in configs, D in (1.0, 3.0, 6.0), τ in (1e-2, 2.5e-3)
    (name != "single_R1" && name != "array_outside" && D != 3.0) && continue     # D scan on two configurations
    g = isfinite(s) ? Cylinders(radius=R, repeats=[s, s], surface_relaxation=θ) : Cylinders(radius=R, surface_relaxation=θ)
    bb = isfinite(s) ? 500 : mr.BoundingBox([-R, -R, -R], [R, R, R])
    r = run_cylinders(g; region=region, R=R, sx=s, sy=s, bounding_box=bb, timestep=τ, D=D, times=[0.0, T], nspins=NSPINS)
    sv = region == :inside ? sv_inside(R) : sv_outside(R, s, s)
    expct = expected_neglog_conventional(θ, D, sv, T)
    m, e = r.neglog_mean[end], r.neglog_sem[end]
    z = (m - expct) / e
    push!(rows, (config=name, R=R, D=D, tau=τ, sqrt_Dtau_over_R=sqrt(D * τ) / R, expected=expct, measured=m, sem=e, z=z, ratio=m / expct))
    @printf("%-14s D = %.0f τ = %.1e: E[−ln M] %.5f vs %.5f (ratio %.4f, z = %+.2f)\n", name, D, τ, m, expct, m / expct, z)
end
write_csv(joinpath(OUT, "mapping.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_pass" => count(r -> abs(r.z) <= 3, rows),
    "n" => length(rows), "max_abs_z" => maximum(abs(r.z) for r in rows), "theta" => θ, "nspins" => NSPINS))
fig = Figure(size=(600, 420))
ax = Axis(fig[1, 1], xlabel="√(Dτ)/R", ylabel="measured / θ√(D/π)(S/V)t", title="conventional mapping on cylinders")
for (name, mk) in zip(unique(r.config for r in rows), (:circle, :rect, :diamond, :utriangle, :dtriangle))
    q = filter(r -> r.config == name, rows)
    errorbars!(ax, [r.sqrt_Dtau_over_R for r in q], [r.ratio for r in q], [r.sem / r.expected for r in q])
    scatter!(ax, [r.sqrt_Dtau_over_R for r in q], [r.ratio for r in q], marker=mk, label=name)
end
hlines!(ax, [1.0], color=:black)
axislegend(ax, position=:lb, labelsize=9)
save(joinpath(OUT, "p3_1_6_mapping.png"), fig)
println("saved to $OUT")
```

- [ ] **Step 2: Run**

Run: `julia --project=research/baseline -t 8 research/phase3/p3_1_6_mapping.jl`
Expected: 18 lines. The ratio is close to 1 where √(Dτ)/R is small. Report any trend with √(Dτ)/R as found; do not assume its sign. |z| > 3 is a measured correction, not a failure. Ledger "Ruling: mapping on curved surfaces = measured ratio(√(Dτ)/R)" if any |z| > 3.

- [ ] **Step 3: Commit**

```bash
git add research/phase3/p3_1_6_mapping.jl research/results/phase3/p3_1_6_mapping
git commit -m "P3.1.6: conventional mapping θ√(D/π)(S/V) tested on curved surfaces vs √(Dτ)/R

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 10: C-B5 runtime vs τ (P3.1.5)

**Files:**
- Create: `research/phase3/p3_1_5_runtime.jl` → `research/results/phase3/p3_1_5_runtime/` (`runtime.csv`, `summary.json`)

**Interfaces:**
- Consumes: `run_cylinders` (Task 5).
- Produces: `runtime.csv` with columns `config, layer, tau, us_per_spin_step`. Used later by P3.6.4.

- [ ] **Step 1: Write the script**

```julia
# P3.1.5 C-B5: runtime per spin-step on cylinders vs τ, without the layer (old code) and with the step layer
# (ρ = 0.02, h = 0.2 both sides). 1 seed, 10000 spins, 1000 steps per run. Single runs on a shared machine: ratios carry
# a load uncertainty of order 10–20 %.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_1_5_runtime.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_1_5_runtime")
const N = 10_000
rows = NamedTuple[]
for (name, s, region) in (("single", Inf, :inside), ("array", S_REF, :outside)), layer in (false, true), τ in (0.1, 1e-2, 1e-3, 1e-4)
    kw = layer ? (layer_rho=0.02, layer_h=0.2) : (;)
    g = isfinite(s) ? Cylinders(radius=R_REF, repeats=[s, s]; kw...) : Cylinders(radius=R_REF; kw...)
    bb = isfinite(s) ? 500 : mr.BoundingBox([-R_REF, -R_REF, -R_REF], [R_REF, R_REF, R_REF])
    run_cylinders(g; region=region, sx=s, sy=s, bounding_box=bb, timestep=τ, times=[0.0, 10τ], seeds=1:1, nspins=100)   # warm-up
    r = run_cylinders(g; region=region, sx=s, sy=s, bounding_box=bb, timestep=τ, times=[0.0, 1000τ], seeds=1:1, nspins=N)
    us = r.runtime / (N * 1000) * 1e6
    push!(rows, (config=name, layer=layer, tau=τ, us_per_spin_step=us))
    @printf("%-6s layer %-5s τ = %.0e: %.4f µs per spin-step\n", name, layer, τ, us)
end
write_csv(joinpath(OUT, "runtime.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rows" => rows))
println("saved to $OUT")
```

- [ ] **Step 2: Run**

Run: `julia --project=research/baseline -t 8 research/phase3/p3_1_5_runtime.jl`
Expected: 16 lines of µs per spin-step. The layer ratio is reported, with no pass criterion. The τ = 0.1 run on a single R = 1 cylinder is far above the tortuosity limit and is only a cost point.

- [ ] **Step 3: Commit**

```bash
git add research/phase3/p3_1_5_runtime.jl research/results/phase3/p3_1_5_runtime
git commit -m "P3.1.5: cylinder C-B5 runtime vs τ, with and without the layer

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 11: Step-layer verification on cylinders (P3.4.1–P3.4.4, P3.2.5, P3.4.7)

**Files:**
- Create: `research/phase3/p3_4_1_trace.jl` → `research/results/phase3/p3_4_1_trace/`
- Create: `research/phase3/p3_4_exact.jl` → `research/results/phase3/p3_4_exact/`
- Create: `research/phase3/p3_4_3_identity.jl` → `research/results/phase3/p3_4_3_identity/`

**Interfaces:**
- Consumes: `run_cylinders`, `brute_rate`, `rate_inside`, `rate_outside`, `roundoff_bound`, `expected_neglog_layer`, `sv_inside`, `sv_outside` (Task 5).
- Produces: three `summary.json` files with `pass` flags for the record.

- [ ] **Step 1: Write the trace script `p3_4_1_trace.jl`**

```julia
# P3.4.1: single-spin trace on a square array with overlapping outer layers. Each step's M⊥ change is predicted from
# the simulator's own reflected path, the layer exponent of every straight piece by brute-force sampling of the summed
# rate (brute_rate, 10⁴ points per piece, independent of the package's Layers code):
#     M⊥,pred = M⊥(previous) · exp(−R2_bulk τ − Σ_j dt_j · mean_k brute_rate(p_j(s_k))),  dt_j = τ ℓ_j/Σℓ.
# Square array s = 2.5, R = 1, step layer h = 0.3 both sides (outer shells overlap: gap 0.5 < 2h), ρ = 0.03,
# R2_bulk = 1/80, D = 3, τ = 0.01, 100 spins × 500 steps from uniform positions in a cell, Random.seed!(20261024).
# Pass: max |M⊥,sim/M⊥,pred − 1| ≤ 1e-6, ≥ 1 step with a reflection inside a shell, ≥ 1 step through an overlap.
# Sampling error per boundary crossing ≤ rate·τ/(2·10⁴) ≈ 6e-8.
# Run: julia --project=research/baseline -t 1 research/phase3/p3_4_1_trace.jl

include(joinpath(@__DIR__, "harness.jl"))
using StaticArrays, LinearAlgebra
const OUT = phase3_outdir("p3_4_1_trace")
const TAU, NSTEPS, NTRACE, H, RHO, R2B, TOL, K = 0.01, 500, 100, 0.3, 0.03, 1 / 80, 1e-6, 10_000
const RIN, ROUT = rate_inside(RHO, H, R_REF), rate_outside(RHO, H, R_REF)
rate(x, y) = brute_rate(x, y; R=R_REF, s=S_REF, h_in=H, h_out=H, rate_in=RIN, rate_out=ROUT)

seq = mr.SequenceParts.SequenceWaveform((([], []), ([], []), ([], [])), [], [], [NSTEPS * TAU], NSTEPS * TAU)
sim = Simulation(seq; geometry=Cylinders(radius=R_REF, repeats=[S_REF, S_REF], layer_rho=RHO, layer_h=H), diffusivity=D_REF, R2=R2B, verbose=false)
part = first(p for p in mr.parts([seq], 0., mr.TimeStep(TAU, Inf)) if p.duration > 0)
@assert part.duration ≈ TAU

Random.seed!(20261024)
max_err, n_reflect_shell, n_overlap = 0.0, 0, 0
rows = NamedTuple[]
for i in 1:NTRACE
    spin = Spin(position=[(rand() - 0.5) * S_REF, (rand() - 0.5) * S_REF, rand()], transverse=1.0, longitudinal=0.0)
    for step in 1:NSTEPS
        m0 = spin.orientations[1].transverse
        proposed = spin.position .+ randn(SVector{3, Float64}) .* sqrt(2 * D_REF * TAU)
        path = mr.Evolve.draw_step!(spin, sim, part, [3.], proposed)
        lens = [norm(path[k + 1] - path[k]) for k in 1:length(path) - 1]
        total = sum(lens)
        expo, overlap = 0.0, false
        for k in eachindex(lens)
            iszero(lens[k]) && continue
            rates = [rate((path[k] + (j - 0.5) / K * (path[k + 1] - path[k]))[1:2]...) for j in 1:K]
            overlap |= any(r -> r > ROUT * 1.5 && r < RIN, rates) || any(r -> r ≈ 2ROUT, rates)
            expo += TAU * lens[k] / total * sum(rates) / K
        end
        pred = m0 * exp(-expo) * exp(-R2B * TAU)
        err = abs(spin.orientations[1].transverse / pred - 1)
        global max_err = max(max_err, err)
        (length(path) > 2 && expo > 0) && (global n_reflect_shell += 1)
        overlap && (global n_overlap += 1)
        i == 1 && push!(rows, (step=step, n_pieces=length(lens), layer_exponent=expo, transverse=spin.orientations[1].transverse, rel_error=err))
    end
end
pass = max_err <= TOL && n_reflect_shell > 0 && n_overlap > 0
@printf("max error %.2e (tol %.0e); steps with a reflection in a shell %d; through an overlap %d → %s\n",
    max_err, TOL, n_reflect_shell, n_overlap, pass ? "PASS" : "FAIL")
write_csv(joinpath(OUT, "trace_spin1.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "max_rel_error" => max_err, "tolerance" => TOL,
    "n_steps_reflection_in_shell" => n_reflect_shell, "n_steps_overlap" => n_overlap, "pass" => pass,
    "params" => Dict("tau" => TAU, "nsteps" => NSTEPS, "nspins" => NTRACE, "h" => H, "rho" => RHO, "R" => R_REF, "s" => S_REF,
        "rate_in" => RIN, "rate_out" => ROUT, "samples_per_piece" => K)))
println("saved to $OUT")
```

The overlap detector counts a sample whose summed rate is about 2·ROUT, i.e. inside two outer shells. RIN/ROUT ≈ 1.35, so `r ≈ 2ROUT` is unambiguous. The first condition is redundant but harmless. If an executor simplifies it to `any(r -> r ≈ 2ROUT, rates)`, that needs no ruling.

- [ ] **Step 2: Write the exact-checks script `p3_4_exact.jl`**

```julia
# P3.2.5 V0b, P3.4.2 uniform inner layer h = R, P3.4.4 additivity, P3.4.7 side-specific. Exact in any diffusion regime.
# Seeds 1–10, 100000 spins, readouts 0:5:50, τ = 1e-2 (V0b: default τ).
# - V0b: ρ = 0, h = 0.3 vs no layer (square array, R2 = 1/80): positions and every M⊥ identical (==).
# - Uniform: single cylinder, inner step h = R = 1, ρ = 0.05 → κ = 1/2, ΔR2 = 2ρ/R = 0.1 /ms everywhere inside, so every
#   inside spin has M⊥ = e^(−(0.1 + R2)t): per spin ≤ 1e-10, ensemble ≤ roundoff_bound; bulk off and on (1/80).
# - Additivity: square array ρ = 0.02, h = 0.2, R2 0 vs 1/80: positions identical, per spin ratio e^(−t/80) ≤ 1e-10.
# - Side-specific: inner-only layer leaves every outside spin at M⊥ = 1 exactly, outer-only leaves every inside spin.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_exact.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_4_exact")
transverse(f, m) = [sp.orientations[1].transverse for (k, sp) in enumerate(f.spins) if m[k]]
bb = mr.BoundingBox([-R_REF, -R_REF, -R_REF], [R_REF, R_REF, R_REF])
res = Dict{String, Any}()
report(name, pass, d) = (res[name] = merge(d, Dict("pass" => pass)); @printf("%-22s %s\n", name, pass ? "PASS" : "FAIL"))

# V0b (both regions)
for region in (:inside, :outside)
    kw = (region=region, sx=S_REF, sy=S_REF, R2_bulk=1 / 80, timestep=nothing)
    a = run_cylinders(Cylinders(radius=R_REF, repeats=[S_REF, S_REF]); kw...)
    b = run_cylinders(Cylinders(radius=R_REF, repeats=[S_REF, S_REF], layer_h=0.3); kw...)
    same = all(mr.position.(x) == mr.position.(y) && [s.orientations[1].transverse for s in x.spins] == [s.orientations[1].transverse for s in y.spins]
               for (x, y) in zip(a.finals, b.finals))
    report("V0b_$(region)", same, Dict("identical" => same))
end
# uniform inner layer
for R2 in (0.0, 1 / 80)
    r = run_cylinders(Cylinders(radius=R_REF, layer_rho_inner=0.05, layer_h_inner=R_REF); region=:inside, bounding_box=bb, R2_bulk=R2)
    target = exp.(-(0.1 + R2) .* TIMES_REF)
    spin_dev = maximum(maximum(abs.(transverse(f, m) ./ target[end] .- 1)) for (f, m) in zip(r.finals, r.masks))
    ens_dev = maximum(abs.(r.mean ./ target .- 1))
    bound = roundoff_bound(TIMES_REF, r.tau, NSPINS_REF)
    report("uniform_R2_$(R2)", spin_dev <= 1e-10 && ens_dev <= bound,
        Dict("max_dev_per_spin" => spin_dev, "max_dev_ensemble" => ens_dev, "bound" => bound))
end
# additivity
g = Cylinders(radius=R_REF, repeats=[S_REF, S_REF], layer_rho=0.02, layer_h=0.2)
for region in (:inside, :outside)
    a = run_cylinders(g; region=region, sx=S_REF, sy=S_REF)
    b = run_cylinders(g; region=region, sx=S_REF, sy=S_REF, R2_bulk=1 / 80)
    same = all(mr.position.(x) == mr.position.(y) for (x, y) in zip(a.finals, b.finals))
    dev = maximum(maximum(abs.(transverse(fb, m) ./ (transverse(fa, m) .* exp(-50 / 80)) .- 1)) for (fa, fb, m) in zip(a.finals, b.finals, a.masks))
    report("additivity_$(region)", same && dev <= 1e-10, Dict("positions_identical" => same, "max_dev_per_spin" => dev,
        "S50_layer_only" => a.mean[end]))
end
# side-specific
for (label, kw, untouched) in (("inner_only", (layer_rho_inner=0.02, layer_h_inner=0.2), :outside),
                               ("outer_only", (layer_rho_outer=0.02, layer_h_outer=0.2), :inside))
    g = Cylinders(radius=R_REF, repeats=[S_REF, S_REF]; kw...)
    u = run_cylinders(g; region=untouched, sx=S_REF, sy=S_REF)
    acted = run_cylinders(g; region=(untouched == :inside ? :outside : :inside), sx=S_REF, sy=S_REF)
    exact = all(all(transverse(f, m) .== 1.0) for (f, m) in zip(u.finals, u.masks))
    report("sides_$(label)", exact && acted.mean[end] < 1, Dict("untouched_exact" => exact, "S50_acted" => acted.mean[end]))
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "results" => res, "pass" => all(v["pass"] for v in values(res))))
println("saved to $OUT")
```

- [ ] **Step 3: Write the identity script `p3_4_3_identity.jl`**

```julia
# P3.4.3: mean-exponent identity on cylinders. With ρ per unit surface area (D3.1) and a uniform stationary density,
#     E[−ln M⊥(t)] = t·⟨ΔR2⟩_compartment = ρ·(S/V)·t
# for any h and τ, including overlapping outer layers (the summed rate integrates to ρ per unit area of every surface
# whose layer reaches the compartment; validation keeps outer layers out of other cylinders). Tests the κ normalisation,
# the shell fractions, the image loop and the overlap sum end to end.
# Configurations (ρ = 0.02, step): single inside R = 1, h = 0.2 and R = 0.5, h = 0.2; square array s = 2.5 inside and
# outside, h = 0.2 (separated) and h_outer = 0.4 (overlap, gap 0.5 < 0.8); rectangular sx = 2.4, sy = 3.0, h_outer = 0.3
# (overlap in x only: gaps 0.4 and 1.0 vs 2h = 0.6). τ = 1e-2 and 4e-2 ms, t = 50 ms, seeds 1–10, 100000 spins.
# z = (mean −ln M − ρ(S/V)t)/SEM, |z| ≤ 3.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_3_identity.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_4_3_identity")
const RHO, T = 0.02, 50.0
configs = [
    ("single_R1", Cylinders(radius=1.0, layer_rho=RHO, layer_h=0.2), 1.0, Inf, Inf, :inside),
    ("single_R0.5", Cylinders(radius=0.5, layer_rho=RHO, layer_h=0.2), 0.5, Inf, Inf, :inside),
    ("square_sep_in", Cylinders(radius=1.0, repeats=[2.5, 2.5], layer_rho=RHO, layer_h=0.2), 1.0, 2.5, 2.5, :inside),
    ("square_sep_out", Cylinders(radius=1.0, repeats=[2.5, 2.5], layer_rho=RHO, layer_h=0.2), 1.0, 2.5, 2.5, :outside),
    ("square_overlap_out", Cylinders(radius=1.0, repeats=[2.5, 2.5], layer_rho=RHO, layer_h_inner=0.2, layer_h_outer=0.4), 1.0, 2.5, 2.5, :outside),
    ("rect_overlap_x_out", Cylinders(radius=1.0, repeats=[2.4, 3.0], layer_rho=RHO, layer_h_inner=0.2, layer_h_outer=0.3), 1.0, 2.4, 3.0, :outside),
]
rows = NamedTuple[]
for (name, g, R, sx, sy, region) in configs, τ in (1e-2, 4e-2)
    bb = isfinite(sx) ? 500 : mr.BoundingBox([-R, -R, -R], [R, R, R])
    r = run_cylinders(g; region=region, R=R, sx=sx, sy=sy, bounding_box=bb, timestep=τ, times=[0.0, T])
    sv = region == :inside ? sv_inside(R) : sv_outside(R, sx, sy)
    e = expected_neglog_layer(RHO, sv, T)
    z = (r.neglog_mean[end] - e) / r.neglog_sem[end]
    push!(rows, (config=name, tau=τ, expected=e, measured=r.neglog_mean[end], sem=r.neglog_sem[end], z=z, S50=r.mean[end], pass=abs(z) <= 3))
    @printf("%-20s τ = %.0e: E[−ln M] %.6f vs ρ(S/V)t %.6f, z = %+.2f\n", name, τ, r.neglog_mean[end], e, z)
end
write_csv(joinpath(OUT, "identity.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_pass" => count(r -> r.pass, rows), "n" => length(rows),
    "pass" => all(r -> r.pass, rows)))
println("saved to $OUT")
```

- [ ] **Step 4: Run all three**

Run:
```
julia --project=research/baseline -t 1 research/phase3/p3_4_1_trace.jl
julia --project=research/baseline -t 8 research/phase3/p3_4_exact.jl
julia --project=research/baseline -t 8 research/phase3/p3_4_3_identity.jl
```
Expected:
- **Trace:** PASS, with max error ≈ 1e-7 or less.
- **Exact checks:** every line PASS.
- **Identity:** 12 lines with |z| ≤ 3. All lines share seeds, so their deviations are correlated.

A failed identity row points at κ, the image loop or the overlap sum. Debug with superpowers:systematic-debugging, starting from the brute-force test of Task 3 with that row's geometry.

- [ ] **Step 5: Commit**

```bash
git add research/phase3/p3_4_1_trace.jl research/phase3/p3_4_exact.jl research/phase3/p3_4_3_identity.jl research/results/phase3/p3_4_1_trace research/results/phase3/p3_4_exact research/results/phase3/p3_4_3_identity
git commit -m "P3.4.1–P3.4.4, P3.2.5, P3.4.7: step layer on cylinders — trace, exact cases, mean-exponent identity

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 12: Record and tracker

**Files:**
- Modify: `research/phase3_record.md`, `research/project-progress-tracker.md`

**Interfaces:**
- Consumes: every `summary.json` and CSV from Tasks 6–11.
- Produces: record sections and tracker ticks.

- [ ] **Step 1: Write the record sections**

Each section of `research/phase3_record.md` follows the record standard: **Design** (choice, alternatives, why), **Method** (formulas, estimators), **Result** (table copied from the JSON/CSV; numbers only from the files), **Conclusion** (tied to the formulas). Sections, in this order, after D3.1–D3.6:
1. **Implementation (P3.2.1–P3.2.4, P3.2.6).**
   - r(s)² = a s² + b s + q, and fraction = L(r_hi²) − L(r_lo²) with L from the stable roots;
   - κ and the surface rates;
   - images within reach R + h_out;
   - validation h_inner ≤ R, h_outer ≤ gap;
   - the summed-rate bound for τ (exact without overlap; upper bound with overlap);
   - the Review Focus items and the tests that pin them.
2. **C-B1/C-B2 (P3.1.1, P3.1.2)** from `p3_1_1_null_bulk`, including τ_default.
3. **C-B3 and P3.4.5** from `p3_1_3_density`:
   - χ² formulas, binomial z, the equal-area bins of r² and r;
   - figure `p3_1_3_density.png`.
4. **C-B4 (P3.1.4)** from `p3_1_4_theta_reference`: a convergence table, and the curves file named as the reference set.
5. **P3.1.6 conventional mapping on curved surfaces:**
   - the collision-rate formula;
   - the table, and figure `p3_1_6_mapping.png`;
   - the conclusion: either "θ√(D/π)(S/V) holds within |z| ≤ 3" or "measured ratio vs √(Dτ)/R, used as the mapping".
6. **C-B5 (P3.1.5)** from `p3_1_5_runtime`.
7. **Verification of the step layer on cylinders:**
   - the trace (P3.4.1), with the sampling bound rate·τ/(2K);
   - the exact cases (P3.2.5, P3.4.2, P3.4.4, P3.4.7), with B = (t/τ + N)ε;
   - the identity (P3.4.3), with E[−ln M] = ρ(S/V)t and why D3.1 makes it exact, including overlap.

Add rows to the Summary table at the top of `phase3_record.md` for each item, with outcome and key number.

- [ ] **Step 2: Update the tracker**

In `research/project-progress-tracker.md`:
- tick (☑) P3.1.1–P3.1.6, P3.2.1–P3.2.6, P3.4.1–P3.4.5 and P3.4.7. Fill the "Output" or pass column with the outcome, as for Phase 2 rows;
- add "(step)" where the item will be repeated for other profiles;
- update the Phase 3 row of the overview table: tasks done count, and status "in progress".

- [ ] **Step 3: Run all test suites**

Run:
```
julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["no-plots"])'
julia --project=research/baseline -t 8 research/phase3/test/runtests.jl
julia --project=research/baseline -t 8 research/phase2/test/runtests.jl
```
Expected: all PASS. The Phase 2 research tests confirm that the walls path is unchanged.

- [ ] **Step 4: Commit**

```bash
git add research/phase3_record.md research/project-progress-tracker.md
git commit -m "Phase 3 record: cylinder baselines, conventional mapping on curved surfaces, step-layer verification; tracker

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## After the last task

Run the final whole-branch review with a fresh subagent reviewer on the most capable model (user rule; see the executing-plans skill). Range: from the commit before Task 1 to HEAD. Give the reviewer this plan, `research/phase3_record.md` (spec D3.1–D3.6) and the Review Focus list. Grade the findings, fix Critical/Important ones with RED→GREEN tests, and ledger minors.
