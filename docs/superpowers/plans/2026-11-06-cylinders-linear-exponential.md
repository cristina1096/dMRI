# Linear and exponential profiles on cylinders (Phase 3.3, 3.4 per profile) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Support the linear and exponential profiles on cylinders: inner and outer side, arrays, and overlapping layers summed. Integrate them along straight pieces by adaptive quadrature between the exact shell roots (D3.3), and verify them as the step profile was verified: closed form, brute force, trace, identity, continuum.

**Architecture:**
- **Exponent of a piece.** The exponent is the integral of g(d(s)/h) over the part of the piece inside the support, multiplied by (ρ/h)/κ.
  - The support is cut into sub-intervals at the exact shell roots (reusing the step code's quadratic roots) and at the point of closest approach s* = −b/2a.
  - On each sub-interval r(s) is smooth and monotone, so 8-point Gauss–Legendre with bisection to a fixed tolerance converges fast.
  - The step profile keeps its exact closed form, unchanged.
- **Exponential support (D3.5).** Inside, the support reaches the axis (support R). Outside, a new `layer_cutoff` field sets it; the default is the smallest surface gap to a neighbouring cylinder.
- **Research code.** The radial continuum gains exact ring averages for linear and exponential (removing the deferred thin-layer issue). The research scripts get a `SHAPE` switch.

**Tech Stack:** Julia 1.12, MCMRSimulator.jl (package `src/`, tests `test/`), research environment `research/baseline` (CairoMakie, SpecialFunctions).

**Spec:**
- **Binding:**
  - `research/phase3_record.md`, design decisions D3.1 (κ = 1 ∓ c̄/R), D3.2 (step first, then linear and exponential), D3.3 (adaptive quadrature between exact roots; linear closed form as its test), D3.5 (exponential support);
  - `research/project-progress-tracker.md` Phase 3: P3.3.1–P3.3.4, P3.4.1/3/6/8 for each profile, P3.6.3.
- **Background:** the Phase 2 exponential conventions (`src/geometries/user/fix.jl`, `wall_layers`: cutoff validation, an inert ρ = 0 side skips validation) and Review C1 (fixed-ρ differences need small τ).

## Global Constraints

- **Branch.** Work only on branch `cc/near-surface-r2`. Never modify or commit on other branches. Do not push; the user pushes.
- **Commits.** Every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Step path unchanged.** The step profile's results must stay **bit-identical** to before this plan (positions and M⊥). There is a test for this in Task 2.
- **No new package dependency.** Gauss–Legendre nodes are computed with LinearAlgebra (Golub–Welsch), which the package already depends on.
- **κ and support.**
  - κ = 1 ∓ c̄/R with c̄ = h/2 (step), h/3 (linear), h(1 − (c + 1)e^(−c))/(1 − e^(−c)) (exponential, c = support/h).
  - Support: h for step and linear. For the exponential, R inside and `layer_cutoff` outside.
- **Validation.**
  - Inner support ≤ R. Step and linear need h_inner ≤ R; the inner exponential support is R by definition.
  - Outer support ≤ the smallest surface gap to every other cylinder or image (existing rule, now on the support).
  - The outer exponential's default cutoff is that gap. A cylinder with no neighbour and an outer exponential layer needs `layer_cutoff`.
- **Quadrature tolerance.** `LAYER_QUAD_TOL[] = 1e-10`, absolute on ∫₀¹ g ds (dimensionless), held in a `Ref` so the tolerance study can vary it. P3.6.3 fixes it with the study's numbers.
- **Model framing.**
  - Conventional model = reference, not truth.
  - No fast-diffusion assumptions.
  - Record signs and orderings as found.
  - Records: Design / Method / Result / Conclusion with formulas; figures where useful.
- **Monte Carlo.** h is large enough that τ_a = 1e-2 (user, 2026-10-25: avoid thin layers). Fixed-ρ differences are checked at two τ.
- **Commands.**
  - Package tests: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders", "layer"])'`; full suite `test_args=["no-plots"]`.
  - Research tests: `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`.
  - Long scripts run in the background, in chunks of ≤ 2 h.
- **Plan end.** A fresh subagent reviewer (most capable model) reviews the whole branch (user rule).

## Review Focus

1. **A piece through the axis** (the inner exponential support reaches r = 0, and r(s) has a kink at s*). The quadrature must split at s* and stay exact against brute force. Test in Task 1.
2. **Pieces tangent to a shell edge, or grazing the support edge** (the linear g is continuous there, the exponential jumps at the truncation). No NaN, and agreement with the closed form or brute force. Test in Task 1.
3. **Outer exponential on a single non-repeating cylinder without `layer_cutoff`.** A clear error, not an infinite support. Test in Task 2.
4. **Bit-identical step results** after the refactor (`disk_interval`, `cylinder_side_exponent`). Test in Task 2.
5. **Image reach and the timestep bound must use the support, not h**, for the exponential, whose support can be larger than h. Test in Task 2 (an outer exponential with h = 0.05 but cutoff 0.5 must still be summed at 0.4 from the surface).

---

## File Structure

- **Modify** `src/geometries/internal/layers.jl`:
  - `centroid` for all profiles; `gauss_legendre`, `GL_X`, `GL_W`, `LAYER_QUAD_TOL`, `quad_adaptive`;
  - `disk_interval`, with `disk_length` built on it;
  - `shell_profile_integral`, `cylinder_side_exponent`, `cylinder_exponent`;
  - support-based reach in `layer_exponent` and `max_layer_rate`.
- **Modify** `src/geometries/user/obstructions/obstructions.jl`: the Cylinder `layer_cutoff` field and the shape description.
- **Modify** `src/geometries/user/fix.jl`: `cylinder_layers` (all shapes, supports, default cutoff), `surface_gaps`, and `check_outer_reach` on the support.
- **Modify** `test/test_layer_cylinders.jl` and `test/test_cli.jl`.
- **Modify** `research/phase3/continuum_radial.jl` (exact ring averages, exponential) and `research/phase3/test/test_continuum_radial.jl`.
- **Modify** `research/phase3/harness.jl` (`g_profile`, `side_prefactor`, `brute_rate_profile`) and `research/phase3/test/test_harness.jl`.
- **Create the scripts:**
  - `research/phase3/p3_3_quadrature.jl`;
  - `research/phase3/p3_4_6d_profiles.jl`;
  - `research/phase3/p3_4_8d_profiles_inside.jl`.
- **Modify the scripts:**
  - `research/phase3/p3_4_1_trace.jl` (`SHAPE`);
  - `research/phase3/p3_4_3_identity.jl` (`SHAPE`).
- **Modify** `research/phase3_record.md` and `research/project-progress-tracker.md`.

---

### Task 1: Quadrature core and profile integral along a piece (P3.3.1–P3.3.3)

**Files:**
- Modify: `src/geometries/internal/layers.jl`
- Test: `test/test_layer_cylinders.jl`

**Interfaces:**
- Produces:
  - `centroid(side)` for :step / :linear / :exponential;
  - `quad_adaptive(f, a, b, tol)::Float64`; `LAYER_QUAD_TOL::Base.RefValue{Float64}`;
  - `disk_interval(a, b, q, X)::Tuple{Float64, Float64}`;
  - `shell_profile_integral(side, R, inside, centre, p0, p1)::Float64`, equal to ∫₀¹ g(d(s)/h, c) ds over the support.

- [ ] **Step 1: Write the failing tests** (append to `test/test_layer_cylinders.jl`)

```julia
@testset "test_layer_cylinders.jl: non-step profiles along a piece (P3.3)" begin
    Layers = mr.Geometries.Internal.Layers
    P(x, y) = SVector{2, Float64}(x, y)
    O = P(0., 0.)
    expside(h, L) = Layers.LayerSide(0.02, h, :exponential, L)

    @testset "centroid and κ for all profiles (D3.1)" begin
        @test Layers.centroid(Layers.LayerSide(0.02, 0.3, :linear)) ≈ 0.1
        c = 2.0                                           # support 0.4, h = 0.2
        @test Layers.centroid(expside(0.2, 0.4)) ≈ 0.2 * (1 - (c + 1) * exp(-c)) / (1 - exp(-c))
        n = 200_000                                       # κ by midpoint quadrature, exponential outside, R = 1
        κq = sum(exp(-(k - 0.5) / n * c) / (1 - exp(-c)) * (1 + (k - 0.5) / n * 0.4 / 1.0) for k in 1:n) / n * c
        @test Layers.kappa(expside(0.2, 0.4), 1.0, false) ≈ κq rtol=1e-8
    end

    @testset "quad_adaptive on polynomials, a kink and a jump-free exponential" begin
        @test Layers.quad_adaptive(s -> s^7, 0.0, 1.0, 1e-12) ≈ 1 / 8 rtol=1e-12
        @test Layers.quad_adaptive(s -> exp(-3s), 0.2, 0.9, 1e-12) ≈ (exp(-0.6) - exp(-2.7)) / 3 rtol=1e-12
        @test Layers.quad_adaptive(s -> sqrt(abs(s - 0.3)), 0.0, 0.3, 1e-10) ≈ 2 / 3 * 0.3^1.5 rtol=1e-8
    end

    @testset "disk_interval and disk_length (refactor keeps disk_length)" begin
        @test Layers.disk_interval(1.0, 0.0, 0.0, 0.25) == (0.0, 0.5)
        @test Layers.disk_length(1.0, 0.0, 0.0, 0.25) ≈ 0.5
        (lo, hi) = Layers.disk_interval(1.0, 0.0, 2.0, 1.0)
        @test hi <= lo
    end

    "∫₀¹ √(a s² + b s + q) ds in closed form (independent of the package; a > 0)."
    function sqrtquad(a, b, q, s0, s1)
        F(s) = (Q = a * s^2 + b * s + q; (2a * s + b) * sqrt(Q) / (4a) + (4a * q - b^2) / (8a^1.5) * log(abs(2sqrt(a) * sqrt(Q) + 2a * s + b)))
        return F(s1) - F(s0)
    end

    @testset "linear profile vs closed form (P3.3.2)" begin
        # inside, h = 0.3, support [0.7, 1]: g = 2(1 − (R − r)/h) is linear in r, so ∫g ds = ∫(2 − 2(R − r)/h) ds over the
        # shell parts, with ∫ r ds from the √(quadratic) closed form on each part.
        side = Layers.LayerSide(0.02, 0.3, :linear)
        R, h = 1.0, 0.3
        for (p0, p1) in ((P(-2., 0.5), P(2., 0.5)), (P(-0.9, -0.2), P(0.95, 0.4)), (P(0.75, 0.), P(0.95, 0.3)))
            Δ = p1 - p0; a = Δ ⋅ Δ; b = 2(p0 ⋅ Δ); q = p0 ⋅ p0
            (h0, h1) = Layers.disk_interval(a, b, q, R^2)
            (l0, l1) = Layers.disk_interval(a, b, q, (R - h)^2)
            parts = l1 > l0 ? [(h0, l0), (l1, h1)] : [(h0, h1)]
            closed = sum((x1 - x0) * (2 - 2R / h) + 2 / h * sqrtquad(a, b, q, x0, x1) for (x0, x1) in parts if x1 > x0; init=0.0)
            @test Layers.shell_profile_integral(side, R, true, O, p0, p1) ≈ closed atol=1e-10
        end
    end

    @testset "wall regression: linear d(s) reproduces side_exponent (P3.3.3)" begin
        # along a radial piece r(s) is linear in s, so d(s) is exactly linear, as on a wall
        R = 2.0
        for shape in (:linear, :exponential)
            side = shape == :linear ? Layers.LayerSide(0.02, 0.3, :linear) : Layers.LayerSide(0.02, 0.3, :exponential, 2.0)
            d0, d1 = 0.05, 0.8
            p0, p1 = P(R - d0, 0.), P(R - d1, 0.)          # radial, inside: d = R − r, exactly linear in s
            wall = Layers.side_exponent(side, d0, d1, 1.0) / (side.rho / side.h)
            @test Layers.shell_profile_integral(side, R, true, O, p0, p1) ≈ wall rtol=1e-9
        end
    end

    @testset "brute force, inside and outside, linear and exponential (P3.3.1)" begin
        Random.seed!(11)
        K = 200_000
        g(shape, u, c) = Layers.profile_g(shape, u, c)
        maxerr = 0.0
        for _ in 1:300, (shape, inside, L, h) in ((:linear, true, 0.3, 0.3), (:linear, false, 0.4, 0.4),
                                                    (:exponential, true, 1.0, 0.2), (:exponential, false, 0.5, 0.15))
            side = shape == :linear ? Layers.LayerSide(0.02, h, :linear) : Layers.LayerSide(0.02, h, :exponential, L)
            p0 = P((3 .* rand(2) .- 1.5)...); p1 = P((3 .* rand(2) .- 1.5)...)
            dist(r) = inside ? 1.0 - r : r - 1.0
            brute = sum(k -> (r = norm(p0 + (k - 0.5) / K * (p1 - p0)); d = dist(r); 0 <= d <= L ? g(shape, d / h, L / h) : 0.0), 1:K) / K
            maxerr = max(maxerr, abs(Layers.shell_profile_integral(side, 1.0, inside, O, p0, p1) - brute))
        end
        @test maxerr < 2e-4            # exponential truncation jump: ≤ g(c)·(crossings)/(2K) + smooth part O(K⁻²)
    end

    @testset "Review Focus 1: piece through the axis (inner exponential reaches r = 0)" begin
        side = Layers.LayerSide(0.02, 0.2, :exponential, 1.0)
        # radial through the centre, from (−1, 0) to (1, 0): ds = dx/2, d = 1 − |x|, so
        # ∫₀¹ g(d/h) ds = ∫₀¹ g(d/h) dd = h ∫₀^c g(u) du = h (g normalised over its support)
        expected = 0.2
        @test Layers.shell_profile_integral(side, 1.0, true, O, P(-1., 0.), P(1., 0.)) ≈ expected rtol=1e-9
        @test isfinite(Layers.shell_profile_integral(side, 1.0, true, O, P(-1., 1e-14), P(1., 1e-14)))
    end

    @testset "Review Focus 2: tangent and grazing pieces" begin
        lin = Layers.LayerSide(0.02, 0.2, :linear)
        @test Layers.shell_profile_integral(lin, 1.0, true, O, P(-2., 1.0), P(2., 1.0)) == 0.0          # tangent to r = R
        v = Layers.shell_profile_integral(lin, 1.0, true, O, P(-2., 0.8), P(2., 0.8))                    # tangent to r = R − h
        @test isfinite(v) && v > 0
        ex = Layers.LayerSide(0.02, 0.1, :exponential, 0.3)
        @test isfinite(Layers.shell_profile_integral(ex, 1.0, false, O, P(-2., 1.3), P(2., 1.3)))        # tangent to the support edge
    end

    @testset "stationary points" begin
        lin = Layers.LayerSide(0.02, 0.2, :linear)
        @test Layers.shell_profile_integral(lin, 1.0, true, O, P(0.9, 0.), P(0.9, 0.)) ≈ Layers.profile_g(:linear, 0.5)
        @test Layers.shell_profile_integral(lin, 1.0, true, O, P(0.5, 0.), P(0.5, 0.)) == 0.0
    end
end
```

(`test/runtests.jl` already imports `LinearAlgebra` and `Random`; `⋅` and `norm` come from LinearAlgebra.)

The expected value in Review Focus 1 is derived in the test's comment: ds = dx/2 and d = 1 − |x|, so ∫₀¹ g(d/h) ds = h ∫₀^c g(u) du = h = 0.2, because g is normalised over its support.

- [ ] **Step 2: Run to verify they fail**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders"])'`
Expected: FAIL with `UndefVarError: quad_adaptive` / `disk_interval` / `shell_profile_integral`, and a centroid error for :linear.

- [ ] **Step 3: Implement** (in `layers.jl`)

Imports: `import LinearAlgebra: norm, eigen, SymTridiagonal`.

Replace `centroid`:

```julia
"""
    centroid(side)

Profile centroid c̄ = h∫u g(u) du (um): step h/2, linear h/3, exponential with support L = c·h
h(1 − (c + 1)e^(−c))/(1 − e^(−c)).
"""
function centroid(side::LayerSide)
    side.shape == :step && return side.h / 2
    side.shape == :linear && return side.h / 3
    c = side.cutoff / side.h
    return side.h * (1 - (c + 1) * exp(-c)) / -expm1(-c)
end
```

Add after `kappa`:

```julia
"n-point Gauss–Legendre nodes and weights on [−1, 1] (Golub–Welsch)."
function gauss_legendre(n::Int)
    β = [k / sqrt(4k^2 - 1) for k in 1:n - 1]
    F = eigen(SymTridiagonal(zeros(n), β))
    return (SVector{n}(F.values), SVector{n}(2 .* F.vectors[1, :] .^ 2))
end
const (GL_X, GL_W) = gauss_legendre(8)

"Absolute tolerance of `quad_adaptive` on ∫₀¹ g ds along a piece (dimensionless); fixed by the P3.6.3 study."
const LAYER_QUAD_TOL = Ref(1e-10)

"8-point Gauss–Legendre estimate of ∫_a^b f."
function gl8(f, a::Float64, b::Float64)
    m = (a + b) / 2
    w = (b - a) / 2
    s = 0.0
    for k in 1:8
        s += GL_W[k] * f(m + w * GL_X[k])
    end
    return s * w
end

"""
    quad_adaptive(f, a, b, tol)

∫_a^b f by 8-point Gauss–Legendre with bisection: accept the two-halves estimate when it agrees with the
whole-interval estimate within tol·(b − a); depth ≤ 30.
"""
function quad_adaptive(f, a::Float64, b::Float64, tol::Float64, depth::Int=0)
    m = (a + b) / 2
    whole = gl8(f, a, b)
    halves = gl8(f, a, m) + gl8(f, m, b)
    (abs(whole - halves) <= tol * (b - a) || depth >= 30) && return halves
    return quad_adaptive(f, a, m, tol, depth + 1) + quad_adaptive(f, m, b, tol, depth + 1)
end
```

Replace `disk_length` by:

```julia
"""
    disk_interval(a, b, q, X)

{s ∈ [0, 1] : a s² + b s + q ≤ X} = [lo, hi] for a > 0 (empty when hi ≤ lo). Roots by the cancellation-free form
t = −(b + sign(b)√disc)/2, roots t/a and (q − X)/t.
"""
function disk_interval(a::Float64, b::Float64, q::Float64, X::Float64)
    disc = b * b - 4 * a * (q - X)
    disc <= 0 && return (0.0, 0.0)
    t = -(b + copysign(sqrt(disc), b)) / 2
    (s1, s2) = minmax(t / a, (q - X) / t)
    return (max(0.0, s1), min(1.0, s2))
end

"""
    disk_length(a, b, q, X)

Length of {s ∈ [0, 1] : a s² + b s + q ≤ X} for a ≥ 0 (see [`disk_interval`](@ref)).
"""
function disk_length(a::Float64, b::Float64, q::Float64, X::Float64)
    iszero(a) && return q <= X ? 1.0 : 0.0
    (lo, hi) = disk_interval(a, b, q, X)
    return max(0.0, hi - lo)
end
```

(This is the same arithmetic as before, so step results are unchanged.)

Add after `shell_fraction`:

```julia
"""
    shell_profile_integral(side, R, inside, centre, p0, p1)

∫₀¹ g(d(s)/h) ds over the part of the straight piece p0 → p1 (local 2D) inside the side's support, with
d = R − r inside (support R − L ≤ r ≤ R) and d = r − R outside (R ≤ r ≤ R + L), L = side.cutoff. The support is split
at the exact shell roots and at the closest approach s* = −b/2a, so r(s) = √(a s² + b s + q) is smooth and monotone on
each part, which is integrated by `quad_adaptive` to `LAYER_QUAD_TOL[]` (D3.3).
"""
function shell_profile_integral(side::LayerSide, R::Float64, inside::Bool, centre::SVector{2, Float64},
                                p0::SVector{2, Float64}, p1::SVector{2, Float64})
    h = side.h
    c = side.cutoff / h
    umax = side.shape == :exponential ? c : 1.0
    u = p0 - centre
    Δ = p1 - p0
    a = Δ[1] * Δ[1] + Δ[2] * Δ[2]
    b = 2 * (u[1] * Δ[1] + u[2] * Δ[2])
    q = u[1] * u[1] + u[2] * u[2]
    r_lo = inside ? max(R - side.cutoff, 0.0) : R
    r_hi = inside ? R : R + side.cutoff
    gd(r) = profile_g(side.shape, clamp((inside ? R - r : r - R) / h, 0.0, umax), c)
    if iszero(a)
        r = sqrt(q)
        return (r_lo <= r <= r_hi) ? gd(r) : 0.0
    end
    (h0, h1) = disk_interval(a, b, q, r_hi^2)
    h1 <= h0 && return 0.0
    (l0, l1) = r_lo > 0 ? disk_interval(a, b, q, r_lo^2) : (0.0, 0.0)
    f(s) = gd(sqrt(max(a * s * s + b * s + q, 0.0)))
    sstar = -b / (2a)
    tol = LAYER_QUAD_TOL[]
    total = 0.0
    parts = l1 > l0 ? ((h0, l0), (l1, h1)) : ((h0, h1), (0.0, 0.0))
    for (x0, x1) in parts
        x1 > x0 || continue
        if x0 < sstar < x1
            total += quad_adaptive(f, x0, sstar, tol) + quad_adaptive(f, sstar, x1, tol)
        else
            total += quad_adaptive(f, x0, x1, tol)
        end
    end
    return total
end
```

- [ ] **Step 4: Run to verify they pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders"])'`
Expected: PASS. All earlier `layer_cylinders` tests still pass, since step is unaffected.

- [ ] **Step 5: Commit**

```bash
git add src/geometries/internal/layers.jl test/test_layer_cylinders.jl
git commit -m "Cylinders: adaptive Gauss–Legendre quadrature of g(d(s)/h) between exact shell roots; linear closed-form and wall checks (P3.3.1–P3.3.3)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Linear and exponential in the cylinder layer: fields, support, validation, simulation (D3.5)

**Files:**
- Modify: `src/geometries/internal/layers.jl`, `src/geometries/user/fix.jl`, `src/geometries/user/obstructions/obstructions.jl`
- Test: `test/test_layer_cylinders.jl`, `test/test_cli.jl`

**Interfaces:**
- Consumes: `shell_profile_integral`, `centroid`, `kappa` (Task 1).
- Produces:
  - `cylinder_side_exponent(side, R, inside, centre, p0, p1, dt)`; the reach R + outer.cutoff in `layer_exponent`;
  - support-based `max_layer_rate`;
  - the user field `layer_cutoff` on `Cylinders`; `surface_gaps(centres, radii, repeats)::Vector{Float64}`.

- [ ] **Step 1: Write the failing tests** (append to `test/test_layer_cylinders.jl`)

```julia
@testset "test_layer_cylinders.jl: linear and exponential on cylinders (D3.5)" begin
    Layers = mr.Geometries.Internal.Layers
    geom(c) = mr.Simulation([]; geometry=c, verbose=false).geometry
    V3(x, y, z=0.) = SVector{3, Float64}(x, y, z)
    per_spin(snap) = [s.orientations[1].transverse for s in snap.spins]

    @testset "resolution: shapes and supports" begin
        l = geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.3, layer_shape="linear"))[1].layer[1]
        @test (l.inner.shape, l.inner.cutoff, l.outer.shape, l.outer.cutoff) == (:linear, 0.3, :linear, 0.3)
        l = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2, layer_shape="exponential"))[1].layer[1]
        @test l.inner.cutoff == 1.0                      # inside: support to the axis (D3.5)
        @test l.outer.cutoff ≈ 0.5                       # outside: default = smallest surface gap
        l = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=0.3))[1].layer[1]
        @test l.outer.cutoff == 0.3
        @test l.inner.cutoff == 1.0                      # layer_cutoff does not change the inner support
    end

    @testset "validation" begin
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2, layer_shape="exponential"))   # Review Focus 3
        @test geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=2.))[1].layer !== nothing
        @test_throws ErrorException geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2,
            layer_shape="exponential", layer_cutoff=0.6))                                                                       # > gap
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=1.5, layer_shape="linear"))          # h_inner > R
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=-1.))
        # inert outer side (ρ = 0) needs no cutoff
        @test geom(mr.Cylinders(radius=1., layer_rho_inner=0.02, layer_h_inner=0.2, layer_shape="exponential"))[1].layer !== nothing
        @test_throws ErrorException geom(mr.Cylinders(radius=1., layer_rho=0.02, layer_h=0.2, layer_shape="bogus"))
    end

    @testset "Review Focus 4: step results bit-identical" begin
        function run(c)
            sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=c, diffusivity=3., verbose=false)
            Random.seed!(21)
            snap = mr.Snapshot(1000, sim, 500; transverse=1., longitudinal=0.)
            return mr.readout(snap, sim, [20.]; return_snapshot=true)[1]
        end
        a = run(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2))
        b = run(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2, layer_shape="step"))
        @test per_spin(a) == per_spin(b)
        # reference values recorded before this plan (Task 2 Step 2 writes them; see the brief)
        @test sum(per_spin(a)) == STEP_REFERENCE_SUM
    end

    @testset "radial step with a reflection: linear and exponential, hand-computed" begin
        seq = build_sequence([1.2, :readout])
        part(sim) = first(p for p in mr.parts([seq], 0., mr.TimeStep(1.2, Inf)) if p.duration > 0)
        n = 200_000
        κq(gfun, umax, h, R) = sum(gfun((k - 0.5) / n * umax) * (1 - (k - 0.5) / n * umax * h / R) for k in 1:n) / n * umax
        # linear inside h = 0.2: radial path 0.5 → 1 → 0.7 (1.5 ms/um); in the shell 0.4 um → 0.6 ms; mean g over a radial
        # traverse of the whole shell is 1, so E = (ρ/(hκ))·0.6 with κ = 1 − h/(3R)
        sim = mr.Simulation(seq; geometry=mr.Cylinders(radius=1., layer_rho_inner=0.02, layer_h_inner=0.2, layer_shape="linear"),
            diffusivity=3., verbose=false)
        spin = mr.Spin(position=[0.5, 0., 0.], transverse=1., longitudinal=0.)
        mr.Evolve.draw_step!(spin, sim, part(sim), [3.], [1.3, 0., 0.])
        κ = κq(u -> 2(1 - u), 1.0, 0.2, 1.0)
        @test spin.orientations[1].transverse ≈ exp(-0.02 / (0.2 * κ) * 0.6) rtol=1e-9
        # exponential inside h = 0.2 (support R = 1, c = 5): d runs 0.5 → 0 then 0 → 0.3, dt = 1.5 dd, so
        # ∫g dt = 1.5 h [(1 − e^(−2.5)) + (1 − e^(−1.5))]/(1 − e^(−5))
        sim = mr.Simulation(seq; geometry=mr.Cylinders(radius=1., layer_rho_inner=0.02, layer_h_inner=0.2, layer_shape="exponential"),
            diffusivity=3., verbose=false)
        spin = mr.Spin(position=[0.5, 0., 0.], transverse=1., longitudinal=0.)
        mr.Evolve.draw_step!(spin, sim, part(sim), [3.], [1.3, 0., 0.])
        κ = κq(u -> exp(-u) / (1 - exp(-5)), 5.0, 0.2, 1.0)
        I = 1.5 * 0.2 * ((1 - exp(-2.5)) + (1 - exp(-1.5))) / (1 - exp(-5))
        @test spin.orientations[1].transverse ≈ exp(-0.02 / (0.2 * κ) * I) rtol=1e-8
    end

    @testset "Review Focus 5: reach and timestep use the support" begin
        g = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho_outer=0.02, layer_h_outer=0.05, layer_shape="exponential"))
        # point on the diagonal at r = 1.4 (0.4 from this surface: beyond h = 0.05, inside the support 0.5); its distances to
        # the neighbouring surfaces are ≥ 0.80, outside their supports, so only this cylinder contributes
        c = 0.5 / 0.05
        κ = Layers.kappa(Layers.LayerSide(0.02, 0.05, :exponential, 0.5), 1.0, false)
        p = V3(1.4 / sqrt(2), 1.4 / sqrt(2))
        @test Layers.layer_exponent(g, p, p, 1.0) ≈ 0.02 / (0.05 * κ) * exp(-0.4 / 0.05) / (1 - exp(-c)) rtol=1e-10
        # overlap bound counts the 4 neighbours whose outer support (R + 0.5) intersects: 2(R + 0.5) = 3 > 2.5
        @test Layers.max_layer_rate(g) ≈ 5 * Layers.cylinder_surface_rate(Layers.LayerSide(0.02, 0.05, :exponential, 0.5), 1.0, false)
    end

    @testset "layer exponent vs brute force, linear and exponential, square array with overlap" begin
        for (shape, kw) in ((:linear, (layer_h_inner=0.3, layer_h_outer=0.4)), (:exponential, (layer_h_inner=0.2, layer_h_outer=0.15)))
            g = geom(mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.03, layer_shape=String(shape); kw...))
            l = g[1].layer[1]
            pre_in = l.inner.rho / l.inner.h / Layers.kappa(l.inner, 1.0, true)
            pre_out = l.outer.rho / l.outer.h / Layers.kappa(l.outer, 1.0, false)
            function rate(x, y)
                t = 0.0
                for i in -2:2, j in -2:2
                    r = hypot(x - i * 2.5, y - j * 2.5)
                    d = 1.0 - r
                    0 <= d <= l.inner.cutoff && (t += pre_in * Layers.profile_g(shape, d / l.inner.h, l.inner.cutoff / l.inner.h))
                    d = r - 1.0
                    0 <= d <= l.outer.cutoff && (t += pre_out * Layers.profile_g(shape, d / l.outer.h, l.outer.cutoff / l.outer.h))
                end
                return t
            end
            Random.seed!(5)
            K = 20_000
            err = 0.0
            for _ in 1:100
                p0 = V3((2.5 .* (2 .* rand(2) .- 1))..., 0.); p1 = V3((2.5 .* (2 .* rand(2) .- 1))..., 0.)
                brute = sum(k -> rate((p0 + (k - 0.5) / K * (p1 - p0))[1:2]...), 1:K) / K
                err = max(err, abs(Layers.layer_exponent(g, p0, p1, 1.0) - brute))
            end
            @test err < 2e-5                     # truncation jumps ≤ 0.007/(2K) per crossing, smooth parts O(K⁻²)
        end
    end

    @testset "simulation: side-specific and additivity hold for linear and exponential" begin
        for shape in ("linear", "exponential")
            c = mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2, layer_shape=shape)
            function run(R2)
                sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=c, diffusivity=3., R2=R2, verbose=false)
                Random.seed!(3)
                snap = mr.Snapshot(1000, sim, 500; transverse=1., longitudinal=0.)
                return mr.readout(snap, sim, [20.]; return_snapshot=true)[1]
            end
            a = run(0.); b = run(1 / 80)
            @test mr.position.(a) == mr.position.(b)
            @test all(isapprox.(per_spin(b), per_spin(a) .* exp(-20 / 80); rtol=1e-10))
            @test minimum(per_spin(a)) < 1.0
        end
    end

    @testset "JSON round trip keeps layer_cutoff and shape" begin
        c = mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2, layer_shape="exponential", layer_cutoff=0.4)
        io = IOBuffer()
        mr.write_geometry(io, c)
        back = mr.read_geometry_json(String(take!(io)))
        @test back.layer_cutoff.value == 0.4
        @test back.layer_shape.value == "exponential"
    end
end
```

`STEP_REFERENCE_SUM`: before changing any source in this task, compute it on the current code and insert the literal value:

```bash
julia --project -e 'import MCMRSimulator as mr; import Random; sim = mr.Simulation([mr.SequenceParts.empty_sequence()]; geometry=mr.Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho=0.02, layer_h=0.2), diffusivity=3., verbose=false); Random.seed!(21); snap = mr.Snapshot(1000, sim, 500; transverse=1., longitudinal=0.); r = mr.readout(snap, sim, [20.]; return_snapshot=true)[1]; println(repr(sum(s.orientations[1].transverse for s in r.spins)))'
```

Replace `STEP_REFERENCE_SUM` with the printed literal (e.g. `998.1234567890123`).

CLI test in `test/test_cli.jl` (after "create cylinders with a layer"):

```julia
        @testset "create cylinders with an exponential layer" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create cylinders 1 test.json --radius 1 --repeats 2.5,2.5 --layer_rho 0.02 --layer_h 0.2 --layer_shape exponential --layer_cutoff 0.4")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["layer_shape"] == "exponential"
                @test result["layer_cutoff"] == 0.4
            end
        end
```

- [ ] **Step 2: Run to verify they fail**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders", "cli"])'`
Expected: FAIL:
- the layer_shape "linear" error from `cylinder_layers` ("only step");
- unknown keyword `layer_cutoff`;
- the bit-identical test passes (a guard, not RED).

- [ ] **Step 3: Implement**

(a) `obstructions.jl`, Cylinder fields: change the `layer_shape` description to "… One of step, linear, exponential." Add after `layer_shape_outer`:

```julia
            Field{Float64}(:layer_cutoff, "Near-surface layer: support (um) of the exponential profile outside the cylinder, measured from the surface; g is truncated there and renormalised. Default: the smallest surface gap to a neighbouring cylinder. Inside, the exponential always reaches the axis (support = radius). Ignored by step and linear."),
```

(b) `layers.jl`: replace `cylinder_exponent` by

```julia
"""
    cylinder_side_exponent(side, R, inside, centre, p0, p1, dt)

∫ΔR2 dt of one side of a cylinder over the straight piece p0 → p1 (local 2D, axis at `centre`) of duration `dt` (ms):
(ρ/h)/κ · dt · ∫₀¹ g(d(s)/h) ds. Step: exact shell fraction (g = 1); linear and exponential: `shell_profile_integral`.
"""
function cylinder_side_exponent(side::LayerSide, R::Float64, inside::Bool, centre::SVector{2, Float64},
                                p0::SVector{2, Float64}, p1::SVector{2, Float64}, dt::Float64)
    iszero(side.rho) && return 0.0
    if side.shape == :step
        (r_lo, r_hi) = inside ? (R - side.h, R) : (R, R + side.h)
        return cylinder_surface_rate(side, R, inside) * dt * shell_fraction(p0, p1, centre, r_lo, r_hi)
    end
    return side.rho / side.h / kappa(side, R, inside) * dt * shell_profile_integral(side, R, inside, centre, p0, p1)
end

"∫ΔR2 dt of both layers of one cylinder (axis at `centre`) over the straight piece p0 → p1 of duration `dt` (ms)."
cylinder_exponent(layer::CylinderLayer, centre::SVector{2, Float64}, p0::SVector{2, Float64}, p1::SVector{2, Float64}, dt::Float64) =
    cylinder_side_exponent(layer.inner, layer.radius, true, centre, p0, p1, dt) +
    cylinder_side_exponent(layer.outer, layer.radius, false, centre, p0, p1, dt)
```

In `layer_exponent(::Vector{CylinderLayer}, …)` use the reach `layer.radius + layer.outer.cutoff`. In `max_layer_rate(::Vector{CylinderLayer}, g)` replace both `outer.h` by `outer.cutoff`, in `reach`. Update the docstrings ("step" → "step, linear, exponential"; "R + h_out" → "R + outer support").

(c) `fix.jl`, `cylinder_layers`:
- `resolve_shape` accepts any of `Internal.Layers.SHAPES`, with the walls' error message;
- compute `surface_gaps` before the sides;
- build the sides as follows. Keep `check_outer_reach`, but compare `li.outer.cutoff` (support) with the gap, and say "support" in its message.

```julia
"""
    surface_gaps(centres, radii, repeats)

Smallest surface gap |c_i − c_j| − R_i − R_j from each cylinder to any other cylinder or image (centres wrapped into
the cell, images within ±1 cell); Inf without a neighbour.
"""
function surface_gaps(centres, radii, repeats)
    shifts = isnothing(repeats) ? [SVector(0., 0.)] : [SVector(i * repeats[1], j * repeats[2]) for i in -1:1 for j in -1:1]
    n = length(centres)
    gaps = fill(Inf, n)
    for i in 1:n, j in 1:n, shift in shifts
        (i == j && iszero(shift)) && continue
        v = SVector{2, Float64}(centres[j]) + shift - SVector{2, Float64}(centres[i])
        gaps[i] = min(gaps[i], sqrt(v[1]^2 + v[2]^2) - radii[i] - radii[j])
    end
    return gaps
end
```

Inside the loop over cylinders (after `R`), with `gaps = surface_gaps(centres, [radius_of(i) for i in 1:n], cylinders.repeats.value)` computed once:

```julia
        (inner, outer) = map(("inner", "outer")) do side
            rho = resolve("rho", side, i)
            h = resolve("h", side, i)
            rho < 0 && error("Cylinder $i, $side side: layer_rho must be >= 0, got $rho.")
            h < 0 && error("Cylinder $i, $side side: layer_h must be >= 0, got $h.")
            rho > 0 && iszero(h) && error("Cylinder $i, $side side: layer_h must be > 0 when layer_rho > 0.")
            shape = resolve_shape(side, i)
            shape == :exponential || return Internal.LayerSide(rho, h, shape)
            iszero(rho) && return Internal.LayerSide(rho, h, shape, h)                  # inert side: no support needed
            side == "inner" && return Internal.LayerSide(rho, h, shape, R)              # D3.5: inside, support to the axis
            given = cylinders.layer_cutoff[i]
            L = isnothing(given) ? gaps[i] : Float64(given)
            isfinite(L) || error("Cylinder $i, outer side: an exponential layer on a cylinder without neighbours needs `layer_cutoff`.")
            L > 0 || error("Cylinder $i, outer side: layer_cutoff must be > 0, got $L.")
            Internal.LayerSide(rho, h, shape, L)
        end
        inner.rho > 0 && inner.shape != :exponential && inner.h > R && error("Cylinder $i: layer_h_inner = $(inner.h) um exceeds the radius ($R um).")
```

Docstring: "Linear and exponential are supported (D3.2); the exponential's support is R inside and `layer_cutoff` (default: the smallest surface gap) outside (D3.5)."

- [ ] **Step 4: Run to verify they pass**

Run: `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders", "layer", "cli", "collisions"])'`
Expected: PASS.

Then run the full suite in the background: `test_args=["no-plots"]` → PASS.

- [ ] **Step 5: Commit**

```bash
git add src test
git commit -m "Cylinders: linear and exponential layers (inner/outer), exponential support (D3.5), layer_cutoff, support-based reach and timestep

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Radial continuum for linear and exponential (exact ring averages)

**Files:**
- Modify: `research/phase3/continuum_radial.jl`, `research/phase3/test/test_continuum_radial.jl`

**Interfaces:**
- Produces: `radial_signal(...; layer=(shape, ρ, h))` with shape ∈ {"step", "linear", "exponential"}; the inner exponential's support is R (c = R/h). Also `radial_g(shape, u, c=Inf)` and `radial_kappa(shape, h, R)`.

- [ ] **Step 1: Write the failing tests** (append to `test_continuum_radial.jl`)

```julia
@testset "radial continuum: linear and exponential (exact ring averages)" begin
    R, D, ρ = 1.0, 3.0, 0.02
    @testset "∫ΔR2 dA = 2πRρ for every profile, including h below the grid spacing" begin
        for (shape, h) in (("linear", 1e-4), ("linear", 0.3), ("exponential", 1e-4), ("exponential", 0.2), ("exponential", 3.0))
            (M, A) = radial_operator(R=R, D=D, layer=(shape, ρ, h), N=1000)
            (M0, _) = radial_operator(R=R, D=D, N=1000)
            V = (M.dv .- M0.dv)            # diagonal difference = V (the stiffness part is identical)
            @test sum(A .* V) ≈ 2π * R * ρ rtol=1e-10
        end
    end
    @testset "κ for the exponential with support R" begin
        c = R / 0.2
        @test radial_kappa("exponential", 0.2, R) ≈ 1 - 0.2 * (1 - (c + 1) * exp(-c)) / (1 - exp(-c)) / R rtol=1e-8
    end
    @testset "initial slope −2ρ/R and Jensen bound" begin
        for (shape, h) in (("linear", 0.3), ("exponential", 0.2))
            δ = 1e-6
            @test (radial_signal(R=R, D=D, times=[δ], layer=(shape, ρ, h))[1] - 1) / δ ≈ -2ρ / R rtol=1e-3
            @test radial_signal(R=R, D=D, times=[50.0], layer=(shape, ρ, h))[1] >= exp(-2ρ / R * 50) * (1 - 1e-9)
        end
    end
    @testset "grid convergence" begin
        for (shape, h) in (("linear", 0.3), ("exponential", 0.2))
            a = radial_signal(R=R, D=D, times=[50.0], layer=(shape, ρ, h), N=1000)[1]
            b = radial_signal(R=R, D=D, times=[50.0], layer=(shape, ρ, h), N=2000)[1]
            @test abs(a - b) < 2e-8
        end
    end
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`
Expected: FAIL. The exponential gives "radial_g: profile exponential is not supported". The linear thin-h conservation fails (0.996).

- [ ] **Step 3: Implement** (in `continuum_radial.jl`)

```julia
"Normalised profile g(u) on its support (written independently of the package); c = support/h (exponential only)."
function radial_g(shape::AbstractString, u::Real, c::Real=Inf)
    shape == "step" && return 0 <= u <= 1 ? 1.0 : 0.0
    shape == "linear" && return 0 <= u <= 1 ? 2 * (1 - u) : 0.0
    shape == "exponential" && return 0 <= u <= c ? exp(-u) / -expm1(-c) : 0.0
    error("radial_g: profile $shape is not supported")
end

"Support of the inner layer (um): h for step and linear, R for the exponential (D3.5)."
radial_support(shape, h, R) = shape == "exponential" ? R : h

"κ = ∫g(u)(1 − u h/R) du over the support (inside, D3.1), by midpoint quadrature."
function radial_kappa(shape, h, R; n=100_000)
    umax = radial_support(shape, h, R) / h
    return sum(radial_g(shape, (k - 0.5) / n * umax, umax) * (1 - (k - 0.5) / n * umax * h / R) for k in 1:n) / n * umax
end

"""
    ring_integral(shape, h, R, lo, hi)

∫_lo^hi g((R − r)/h) r dr exactly (antiderivatives in r), clipped to the support [R − L, R].
step: r²/2; linear: 2[(1 − R/h) r²/2 + r³/(3h)]; exponential: h e^((r − R)/h)(r − h)/(1 − e^(−c)).
"""
function ring_integral(shape, h, R, lo, hi)
    L = radial_support(shape, h, R)
    a, b = max(lo, R - L), min(hi, R)
    b <= a && return 0.0
    if shape == "step"
        return (b^2 - a^2) / 2
    elseif shape == "linear"
        F(r) = 2 * ((1 - R / h) * r^2 / 2 + r^3 / (3h))
        return F(b) - F(a)
    else
        c = L / h
        G(r) = h * exp((r - R) / h) * (r - h)
        return (G(b) - G(a)) / -expm1(-c)
    end
end
```

In `radial_operator`, replace the whole `if !isnothing(layer)` block by:

```julia
    if !isnothing(layer)
        (shape, ρ, h) = layer
        κ = radial_kappa(shape, h, R)
        for i in 1:N
            V[i] = ρ / (h * κ) * 2π * ring_integral(shape, h, R, redge[i], redge[i + 1]) / A[i]
        end
    end
```

(The step result is unchanged: 2π(b² − a²)/2 = π(b² − a²). Remove the now unused `subsamples` keyword.)

- [ ] **Step 4: Run to verify they pass**

Run: `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`
Expected: PASS. All earlier continuum tests still pass, including the exact step eigenvalue.

- [ ] **Step 5: Commit**

```bash
git add research/phase3/continuum_radial.jl research/phase3/test/test_continuum_radial.jl
git commit -m "Radial continuum: exact ring averages for step, linear and exponential (support R inside)

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Harness brute-force rate for any profile

**Files:**
- Modify: `research/phase3/harness.jl`, `research/phase3/test/test_harness.jl`

**Interfaces:**
- Produces:
  - `g_profile(shape, u, c)`;
  - `side_prefactor(ρ, shape, h, L, R, inside)` = (ρ/h)/κ, with κ by quadrature;
  - `brute_rate_profile(x, y; R, s, inner, outer)`, where `inner`/`outer` are `(shape=..., h=..., L=..., pre=...)` or `nothing`.

- [ ] **Step 1: Write the failing tests** (append to `test_harness.jl`)

```julia
@testset "harness: brute_rate_profile" begin
    pin = side_prefactor(0.02, "linear", 0.3, 0.3, 1.0, true)
    @test pin ≈ 0.02 / (0.3 * (1 - 0.1)) rtol=1e-8
    inner = (shape="linear", h=0.3, L=0.3, pre=pin)
    @test brute_rate_profile(0.85, 0.0; R=1.0, s=Inf, inner=inner, outer=nothing) ≈ pin * 2 * (1 - 0.15 / 0.3)
    pout = side_prefactor(0.02, "exponential", 0.1, 0.5, 1.0, false)
    outer = (shape="exponential", h=0.1, L=0.5, pre=pout)
    v = brute_rate_profile(1.25, 0.0; R=1.0, s=2.5, inner=nothing, outer=outer)       # 0.25 from both neighbours
    @test v ≈ 2 * pout * exp(-2.5) / (1 - exp(-5)) rtol=1e-12
    # the step case agrees with the old brute_rate
    pst = side_prefactor(0.02, "step", 0.2, 0.2, 1.0, true)
    @test pst ≈ rate_inside(0.02, 0.2, 1.0) rtol=1e-8
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`
Expected: FAIL with `UndefVarError: side_prefactor`.

- [ ] **Step 3: Implement** (append to `harness.jl`)

```julia
"Normalised profile g(u), c = support/h (independent of the package)."
g_profile(shape, u, c) = shape == "step" ? (0 <= u <= 1 ? 1.0 : 0.0) :
                         shape == "linear" ? (0 <= u <= 1 ? 2 * (1 - u) : 0.0) :
                         (0 <= u <= c ? exp(-u) / -expm1(-c) : 0.0)

"(ρ/h)/κ with κ = ∫g(u)(1 ∓ u h/R) du over the support L by midpoint quadrature (D3.1; − inside, + outside)."
function side_prefactor(ρ, shape, h, L, R, inside; n=200_000)
    c = L / h
    κ = sum(g_profile(shape, (k - 0.5) / n * c, c) * (1 + (inside ? -1 : 1) * (k - 0.5) / n * c * h / R) for k in 1:n) / n * c
    return ρ / (h * κ)
end

"""
    brute_rate_profile(x, y; R, s, inner, outer)

Summed ΔR2 at (x, y) for a square array of spacing `s` (one cylinder per cell at the origin; `s = Inf` for a single
cylinder) from the 5 × 5 nearest images; `inner`/`outer` = (shape, h, L, pre) or `nothing`. Independent of the package.
"""
function brute_rate_profile(x, y; R, s, inner, outer)
    cells = isfinite(s) ? (-2:2) : (0:0)
    sp = isfinite(s) ? s : 0.0
    cx = isfinite(s) ? round(x / s) * s : 0.0
    cy = isfinite(s) ? round(y / s) * s : 0.0
    total = 0.0
    for i in cells, j in cells
        r = hypot(x - cx - i * sp, y - cy - j * sp)
        if !isnothing(inner) && 0 <= R - r <= inner.L
            total += inner.pre * g_profile(inner.shape, (R - r) / inner.h, inner.L / inner.h)
        end
        if !isnothing(outer) && 0 <= r - R <= outer.L
            total += outer.pre * g_profile(outer.shape, (r - R) / outer.h, outer.L / outer.h)
        end
    end
    return total
end
```

- [ ] **Step 4: Run to verify it passes**

Run: `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add research/phase3/harness.jl research/phase3/test/test_harness.jl
git commit -m "Phase 3 harness: brute-force summed rate and κ prefactor for any profile

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Quadrature tolerance and cost (P3.3.4, P3.6.3)

**Files:**
- Create: `research/phase3/p3_3_quadrature.jl` → `research/results/phase3/p3_3_quadrature/` (`tolerance.csv`, `runtime.csv`, `summary.json`)

- [ ] **Step 1: Write the script**

```julia
# P3.3.4 / P3.6.3: error of the per-piece profile integral vs the quadrature tolerance, and runtime per spin-step.
# (1) Tolerance: 20000 random pieces (endpoints uniform in [−1.6, 1.6]², R = 1) for linear (h = 0.3) and exponential
#     (h = 0.2, support R inside; h = 0.15, support 0.5 outside); I(tol) = shell_profile_integral at LAYER_QUAD_TOL = tol
#     against tol = 1e-14; reported: max |error| and mean number of g evaluations (counted by a wrapper).
#     The exponent error of a piece is (ρ/h)/κ · dt · |error| ≤ ΔR2(0) τ |error|/g(0), so |error| ≤ 1e-6 keeps the
#     per-step exponent error ≤ 1e-6 ΔR2(0) τ, far below the timestep error ≈ 1e-3 relative in R (P2.6).
# (2) Runtime: MCMR µs per spin-step, square array s = 2.5, R = 1, ρ = 0.02, h = 0.2 both sides, step / linear /
#     exponential, τ = 1e-2, at the adopted tolerance and at 1e-6; 10000 spins × 1000 steps, 1 seed.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_3_quadrature.jl

include(joinpath(@__DIR__, "harness.jl"))
using StaticArrays
const OUT = phase3_outdir("p3_3_quadrature")
Layers = mr.Geometries.Internal.Layers
P(x, y) = SVector{2, Float64}(x, y)
O = P(0.0, 0.0)
Random.seed!(2026)
pieces = [(P((3.2 .* rand(2) .- 1.6)...), P((3.2 .* rand(2) .- 1.6)...)) for _ in 1:20_000]
cases = [("linear_in", Layers.LayerSide(0.02, 0.3, :linear), true), ("linear_out", Layers.LayerSide(0.02, 0.3, :linear), false),
         ("exp_in", Layers.LayerSide(0.02, 0.2, :exponential, 1.0), true), ("exp_out", Layers.LayerSide(0.02, 0.15, :exponential, 0.5), false)]
rows = NamedTuple[]
default_tol = Layers.LAYER_QUAD_TOL[]
for (name, side, inside) in cases
    Layers.LAYER_QUAD_TOL[] = 1e-14
    ref = [Layers.shell_profile_integral(side, 1.0, inside, O, p0, p1) for (p0, p1) in pieces]
    for tol in (1e-4, 1e-6, 1e-8, 1e-10, 1e-12)
        Layers.LAYER_QUAD_TOL[] = tol
        t0 = time_ns()
        v = [Layers.shell_profile_integral(side, 1.0, inside, O, p0, p1) for (p0, p1) in pieces]
        ns = (time_ns() - t0) / length(pieces)
        err = maximum(abs.(v .- ref))
        push!(rows, (case=name, tol=tol, max_abs_error=err, ns_per_piece=ns))
        @printf("%-10s tol %.0e: max |error| %.1e, %.0f ns per piece\n", name, tol, err, ns)
    end
end
Layers.LAYER_QUAD_TOL[] = default_tol
write_csv(joinpath(OUT, "tolerance.csv"), rows)

rt = NamedTuple[]
for shape in ("step", "linear", "exponential"), tol in (default_tol, 1e-6)
    shape == "step" && tol != default_tol && continue
    Layers.LAYER_QUAD_TOL[] = tol
    g = Cylinders(radius=R_REF, repeats=[S_REF, S_REF], layer_rho=0.02, layer_h=0.2, layer_shape=shape)
    run_cylinders(g; region=:outside, sx=S_REF, sy=S_REF, timestep=1e-2, times=[0.0, 0.1], seeds=1:1, nspins=100)
    r = run_cylinders(g; region=:outside, sx=S_REF, sy=S_REF, timestep=1e-2, times=[0.0, 10.0], seeds=1:1, nspins=10_000)
    us = r.runtime / (10_000 * 1000) * 1e6
    push!(rt, (shape=shape, tol=tol, us_per_spin_step=us))
    @printf("%-11s tol %.0e: %.4f µs per spin-step\n", shape, tol, us)
end
Layers.LAYER_QUAD_TOL[] = default_tol
write_csv(joinpath(OUT, "runtime.csv"), rt)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "default_tol" => default_tol, "tolerance" => rows, "runtime" => rt))
println("saved to $OUT")
```

- [ ] **Step 2: Run** (≈ 10 min)

Run: `julia --project=research/baseline -t 8 research/phase3/p3_3_quadrature.jl`
Expected:
- the max error falls with tol and is ≤ tol at every tolerance (the per-interval tolerance is conservative);
- the runtime of linear and exponential is a modest multiple of step.

If the max error at the default tolerance exceeds 1e-8, rule a smaller default (ledger) and update `LAYER_QUAD_TOL`. If runtime at 1e-10 is more than 2× that at 1e-6 and the error at 1e-6 is ≤ 1e-6, the default may be relaxed to 1e-8 (ledger ruling with the numbers).

- [ ] **Step 3: Commit**

```bash
git add research/phase3/p3_3_quadrature.jl research/results/phase3/p3_3_quadrature
git commit -m "P3.3.4/P3.6.3: quadrature error vs tolerance and runtime per profile on cylinders

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Trace and identity for linear and exponential (P3.4.1, P3.4.3)

**Files:**
- Modify: `research/phase3/p3_4_1_trace.jl`, `research/phase3/p3_4_3_identity.jl` (a `SHAPE` switch; step results and paths unchanged)

- [ ] **Step 1: Add the SHAPE switch to the trace**

In `p3_4_1_trace.jl`:
- `const SHAPE = get(ENV, "SHAPE", "step")`;
- `OUT = phase3_outdir(SHAPE == "step" ? "p3_4_1_trace" : "p3_4_1_trace_" * SHAPE)`;
- geometry `Cylinders(radius=R_REF, repeats=[S_REF, S_REF], layer_rho=RHO, layer_h=H, layer_shape=SHAPE)`;
- the predicted rate:

```julia
const LIN = SHAPE == "exponential" ? R_REF : H                     # inner support
const LOUT = SHAPE == "exponential" ? S_REF - 2R_REF : H          # outer support (default cutoff = gap)
const INNER = (shape=SHAPE, h=H, L=LIN, pre=side_prefactor(RHO, SHAPE, H, LIN, R_REF, true))
const OUTER = (shape=SHAPE, h=H, L=LOUT, pre=side_prefactor(RHO, SHAPE, H, LOUT, R_REF, false))
rate(x, y) = brute_rate_profile(x, y; R=R_REF, s=S_REF, inner=INNER, outer=OUTER)
```

The overlap detector becomes "≥ 2 outer supports contain the sample". Count them directly:

```julia
nout(x, y) = count(((i, j),) -> (r = hypot(x - round(x / S_REF) * S_REF - i * S_REF, y - round(y / S_REF) * S_REF - j * S_REF); 0 <= r - R_REF <= LOUT), [(i, j) for i in -2:2 for j in -2:2])
```

Then `overlap |= any(k -> nout(sample_k...) >= 2, …)`. Evaluate it on every 100th sample to keep the cost down. Write `samples_per_piece`, `shape`, `L_in`, `L_out` to the summary. For the exponential the truncation jump at the support edge (r = R + 0.5, the neighbour's surface, where spins cannot be) does not occur inside the domain, so the sampling bound stays ≈ rate·τ/(2K) per crossing of a linear kink.

- [ ] **Step 2: Add the SHAPE switch to the identity**

In `p3_4_3_identity.jl`:
- `SHAPE` env, `OUT` suffix;
- every `Cylinders(...)` gets `layer_shape=SHAPE`;
- `expected_neglog_layer(RHO, sv, T)` is unchanged, since ρ is per unit surface area for every profile (D3.1).

The configurations stay as they are. Two notes for the exponential:
- inside, h is the decay length and the support is R;
- outside, the default cutoff is the gap: 0.5 in the square array, 0.4 in the rectangular one.

The single-cylinder rows are inside only, so no `layer_cutoff` is needed.

- [ ] **Step 3: Run** (background; trace ≈ 5 min per shape at 1 thread; identity ≈ 40 min per shape)

```
SHAPE=linear      julia --project=research/baseline -t 1 research/phase3/p3_4_1_trace.jl
SHAPE=exponential julia --project=research/baseline -t 1 research/phase3/p3_4_1_trace.jl
SHAPE=linear      julia --project=research/baseline -t 8 research/phase3/p3_4_3_identity.jl
SHAPE=exponential julia --project=research/baseline -t 8 research/phase3/p3_4_3_identity.jl
```

Expected:
- **trace:** PASS, max error ≤ 1e-6, reflections in shells and overlap steps both > 0;
- **identity:** 12/12, |z| ≤ 3 for each shape.

- [ ] **Step 4: Commit**

```bash
git add research/phase3/p3_4_1_trace.jl research/phase3/p3_4_3_identity.jl research/results/phase3/p3_4_1_trace_* research/results/phase3/p3_4_3_identity_*
git commit -m "P3.4.1, P3.4.3 for linear and exponential on cylinders: single-spin trace and mean-exponent identity

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Monte Carlo vs continuum for linear and exponential, inside (P3.4.6 per profile)

**Files:**
- Create: `research/phase3/p3_4_6d_profiles.jl` → `research/results/phase3/p3_4_6d_profiles/`

- [ ] **Step 1: Write the script**

```julia
# P3.4.6 for linear and exponential: Monte Carlo vs the radial continuum inside a single cylinder R = 1, D = 3, ρ = 0.02.
# Models: conventional (θ = ρ√(π/D)), linear h = 0.3, 0.6 (centroids 0.1, 0.2), exponential h = 0.1, 0.2 (support R;
# centroids ≈ 0.1, 0.2). No thin layers (user request 2026-10-25): ΔR2(0) ≤ 0.23 /ms, so τ_a = 1e-2 and τ_b = 2.5e-3.
# All models share τ (identical trajectories, checked), ΔS = S_layer − S_conv paired per seed.
# Criterion (tracker P3.4.6): at τ_b, max over readouts |z| ≤ 3 for S and for ΔS vs the continuum (3, not 2: 11
# correlated readouts). Seeds 1–10, NSPINS (default 20000) in the box [−R, R]³, readouts 0:5:50.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_6d_profiles.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
const OUT = phase3_outdir("p3_4_6d_profiles")
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const RHO = 0.02
const θ = RHO * sqrt(π / D_REF)
const TAUS = [1e-2, 2.5e-3]
bb = mr.BoundingBox([-R_REF, -R_REF, -R_REF], [R_REF, R_REF, R_REF])
models = [("conventional", Cylinders(radius=R_REF, surface_relaxation=θ), nothing),
          ("linear_h0.3", Cylinders(radius=R_REF, layer_rho_inner=RHO, layer_h_inner=0.3, layer_shape="linear"), ("linear", 0.3)),
          ("linear_h0.6", Cylinders(radius=R_REF, layer_rho_inner=RHO, layer_h_inner=0.6, layer_shape="linear"), ("linear", 0.6)),
          ("exponential_h0.1", Cylinders(radius=R_REF, layer_rho_inner=RHO, layer_h_inner=0.1, layer_shape="exponential"), ("exponential", 0.1)),
          ("exponential_h0.2", Cylinders(radius=R_REF, layer_rho_inner=RHO, layer_h_inner=0.2, layer_shape="exponential"), ("exponential", 0.2))]
cont = Dict(n => (isnothing(p) ? radial_signal(R=R_REF, D=D_REF, times=TIMES_REF, surface=RHO) :
                                 radial_signal(R=R_REF, D=D_REF, times=TIMES_REF, layer=(p[1], RHO, p[2]))) for (n, _, p) in models)
runs = Dict((n, τ) => run_cylinders(g; region=:inside, bounding_box=bb, timestep=τ, nspins=NSPINS) for τ in TAUS for (n, g, _) in models)
identical = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs[("conventional", τ)].finals, runs[(n, τ)].finals))
                for τ in TAUS for (n, _, _) in models)
@printf("identical trajectories across models (per τ): %s\n", identical)
rows = NamedTuple[]
for τ in TAUS, (n, _, p) in models
    r = runs[(n, τ)]
    zS = (r.mean .- cont[n]) ./ r.sem
    d = r.S .- runs[("conventional", τ)].S
    dm = vec(mean(d, dims=1)); ds = vec(std(d, dims=1)) ./ sqrt(size(d, 1))
    dc = cont[n] .- cont["conventional"]
    zD = isnothing(p) ? zeros(length(dm)) : (dm .- dc) ./ ds
    push!(rows, (model=n, tau=τ, S50_mc=r.mean[end], S50_sem=r.sem[end], S50_cont=cont[n][end], maxabs_zS=maximum(abs.(zS[2:end])),
        dS50_mc=dm[end], dS50_sem=ds[end], dS50_cont=dc[end], maxabs_zdS=maximum(abs.(zD[2:end]))))
    @printf("%-18s τ = %.1e: S(50) %.6f ± %.1e vs %.6f (max|z| %.2f); ΔS(50) %+.2e ± %.1e vs %+.2e (max|z| %.2f)\n",
        n, τ, r.mean[end], r.sem[end], cont[n][end], maximum(abs.(zS[2:end])), dm[end], ds[end], dc[end], maximum(abs.(zD[2:end])))
end
final = filter(r -> r.tau == TAUS[end], rows)
pass = identical && all(r -> r.maxabs_zS <= 3 && r.maxabs_zdS <= 3, final)
write_csv(joinpath(OUT, "comparison.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "pass" => pass, "identical_trajectories" => identical,
    "taus" => TAUS, "nspins" => NSPINS, "criterion" => "tau_b, max|z| <= 3 for S and paired dS"))
println(pass ? "PASS" : "FAIL")
```

- [ ] **Step 2: Run** (background, ≈ 1 h)

Expected: identical trajectories; at τ_b every model passes. If a profile fails at τ_b, add a τ_a/16 pass for it and for the conventional model, following P3.4.6 (ledger ruling); the criterion is unchanged.

- [ ] **Step 3: Commit**

```bash
git add research/phase3/p3_4_6d_profiles.jl research/results/phase3/p3_4_6d_profiles
git commit -m "P3.4.6 for linear and exponential: Monte Carlo vs radial continuum inside, paired ΔS at two τ

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Attenuation vs h and the centroid comparison inside, all profiles (P3.4.8 per profile, inside)

**Files:**
- Create: `research/phase3/p3_4_8d_profiles_inside.jl` → `research/results/phase3/p3_4_8d_profiles_inside/` (`sweep.csv`, `centroid.csv`, `summary.json`, `p3_4_8d_profiles_inside.png`)

- [ ] **Step 1: Write the script**

```julia
# P3.4.8 inside for all profiles, from the radial continuum (validated for step in P3.4.6, linear/exponential in Task 7).
# Fixed ρ = 0.02, D = 3, t = 50 ms, R = 0.5, 1, 2.
# (1) h sweep: step and linear h/R = 0.005 … 1; exponential h/R = 0.005 … 5 (support R). S(50), J(50) = ln S + ρ(2/R)·50,
#     ΔS(50) vs the conventional model; monotonicity: sign of neighbouring steps (resolved if > 10 × grid error).
# (2) Same centroid (Phase 2 P2.5.6 found profiles nearly coincide at the same centroid on walls): for c̄/R = 0.02, 0.05,
#     0.1, 0.2, choose h per profile with that centroid (step h = 2c̄, linear h = 3c̄, exponential by bisection with
#     support R), and report ΔS(50) per profile and the spread max|a − b|/max(|a|, |b|) over profile pairs; also the same
#     comparison at the same h.
# Run: julia --project=research/baseline -t 4 research/phase3/p3_4_8d_profiles_inside.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_4_8d_profiles_inside")
const RHO = 0.02
const T = 50.0
HR = Dict("step" => [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0],
          "linear" => [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0],
          "exponential" => [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 5.0])
centroid_of(shape, h, R) = shape == "step" ? h / 2 : shape == "linear" ? h / 3 :
                           (c = R / h; h * (1 - (c + 1) * exp(-c)) / -expm1(-c))
function h_for_centroid(shape, cbar, R)
    shape == "step" && return 2cbar
    shape == "linear" && return 3cbar
    lo, hi = 1e-6, 100R                     # exponential centroid increases monotonically with h, from 0 to R/2
    for _ in 1:200
        mid = sqrt(lo * hi)
        centroid_of("exponential", mid, R) < cbar ? (lo = mid) : (hi = mid)
    end
    return sqrt(lo * hi)
end
S50(shape, h, R) = radial_signal(R=R, D=D_REF, times=[T], layer=(shape, RHO, h))[1]

rows = NamedTuple[]
steps = NamedTuple[]
for R in (0.5, 1.0, 2.0)
    Sc = radial_signal(R=R, D=D_REF, times=[T], surface=RHO)[1]
    E = RHO * (2 / R) * T
    for shape in ("step", "linear", "exponential")
        prev = nothing
        for hr in HR[shape]
            h = hr * R
            S = S50(shape, h, R)
            ge = abs(S - radial_signal(R=R, D=D_REF, times=[T], layer=(shape, RHO, h), N=4000)[1])
            row = (R=R, shape=shape, h=h, h_over_R=hr, centroid_over_R=centroid_of(shape, h, R) / R, S50=S, J50=log(S) + E,
                   dS50=S - Sc, grid_err=ge)
            push!(rows, row)
            if !isnothing(prev)
                d = S - prev.S50
                push!(steps, (R=R, shape=shape, h_from=prev.h, h_to=h, dS=d, resolved=abs(d) > 10 * max(ge, prev.grid_err), sign=sign(d)))
            end
            prev = row
        end
    end
end
cen = NamedTuple[]
for R in (0.5, 1.0, 2.0), cb in (0.02, 0.05, 0.1, 0.2)
    Sc = radial_signal(R=R, D=D_REF, times=[T], surface=RHO)[1]
    vals = Dict(sh => S50(sh, h_for_centroid(sh, cb * R, R), R) - Sc for sh in ("step", "linear", "exponential"))
    hvals = Dict(sh => S50(sh, 2cb * R, R) - Sc for sh in ("step", "linear", "exponential"))     # same h = 2c̄ (step's)
    spread(v) = maximum(abs(v[a] - v[b]) / max(abs(v[a]), abs(v[b])) for a in keys(v) for b in keys(v) if a < b)
    push!(cen, (R=R, centroid_over_R=cb, dS_step=vals["step"], dS_linear=vals["linear"], dS_exponential=vals["exponential"],
        spread_same_centroid=spread(vals), spread_same_h=spread(hvals)))
    @printf("R = %.1f c̄/R = %.2f: ΔS step %+.3e linear %+.3e exponential %+.3e; spread same c̄ %.3f, same h %.3f\n",
        R, cb, vals["step"], vals["linear"], vals["exponential"], spread(vals), spread(hvals))
end
signs = unique(s.sign for s in steps if s.resolved)
@printf("h steps resolved: %d of %d; signs: %s\n", count(s -> s.resolved, steps), length(steps), signs)
write_csv(joinpath(OUT, "sweep.csv"), rows)
write_csv(joinpath(OUT, "steps.csv"), steps)
write_csv(joinpath(OUT, "centroid.csv"), cen)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rho" => RHO, "D" => D_REF, "t" => T,
    "n_steps" => length(steps), "n_resolved" => count(s -> s.resolved, steps), "resolved_signs" => signs,
    "max_spread_same_centroid" => maximum(c.spread_same_centroid for c in cen), "max_spread_same_h" => maximum(c.spread_same_h for c in cen)))

fig = Figure(size=(1100, 400))
ax1 = Axis(fig[1, 1], xlabel="h/R", ylabel="ΔS(50) = S − S_conv", xscale=log10, title="R = 1, fixed ρ = 0.02: same h")
ax2 = Axis(fig[1, 2], xlabel="c̄/R", ylabel="ΔS(50) = S − S_conv", xscale=log10, title="R = 1: same centroid")
for (sh, col) in zip(("step", "linear", "exponential"), (:steelblue, :darkorange, :seagreen))
    q = filter(r -> r.R == 1.0 && r.shape == sh, rows)
    scatterlines!(ax1, [r.h_over_R for r in q], [r.dS50 for r in q], color=col, label=sh)
    scatterlines!(ax2, [r.centroid_over_R for r in q], [r.dS50 for r in q], color=col, label=sh)
end
axislegend(ax1, position=:lb, labelsize=9)
save(joinpath(OUT, "p3_4_8d_profiles_inside.png"), fig)
println("saved to $OUT")
```

- [ ] **Step 2: Run** (≈ 5–10 min, continuum only)

Expected:
- every resolved h step has one sign (recorded as found);
- the centroid spread is reported. Phase 2 walls found ≤ 11 %; whether the same holds on cylinders is a result.

- [ ] **Step 3: Commit**

```bash
git add research/phase3/p3_4_8d_profiles_inside.jl research/results/phase3/p3_4_8d_profiles_inside
git commit -m "P3.4.8 inside for all profiles (continuum): attenuation vs h and same-centroid comparison at fixed ρ

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8b: Attenuation vs h outside for linear and exponential (P3.4.8 per profile, outside; user request 2026-11-06)

**Files:**
- Create: `research/phase3/p3_4_8e_profiles_outside.jl` → `research/results/phase3/p3_4_8e_profiles_outside_<shape>/` (`sweep.csv`, `steps.csv`, `summary.json`)

**Interfaces:**
- Consumes: `run_cylinders`, `sv_outside` (harness); the package's linear and exponential layers (Task 2).
- Produces: per shape, S(50), J(50), E[−ln M] and the paired ΔS vs the conventional model per h, at τ_a and τ_b. Also the paired neighbour steps.

- [ ] **Step 1: Write the script**

```julia
# P3.4.8 outside for linear and exponential (user request 2026-11-06): square array s = 3.0, R = 1 (gap 1.0), outer layer
# only, ρ = 0.02, D = 3, plus the conventional model at the same τ (paired: identical trajectories, checked).
#   linear: h = 0.1, 0.2, 0.3, 0.5, 0.7, 1.0 (support h; outer layers overlap for h > 0.5);
#   exponential: h = 0.05, 0.1, 0.2, 0.3, 0.5, 1.0 (decay length; support = default cutoff = gap 1.0, so outer supports
#   always reach the neighbouring surfaces and overlap, as the exponential on walls).
# Linear starts at h = 0.1: its surface rate 2ρ/(hκ) at h = 0.05 (0.79 /ms) would force τ < 1e-2 (user: avoid thin layers).
# τ_a = 1e-2 and τ_b = 2.5e-3 (Review C1: fixed-ρ differences checked at two τ). Seeds 1–10, NSPINS (default 20000), readouts
# 0:5:50 — the same seeds, spins and τ as the step run at s = 3.0 (p3_4_8_h_sweep_outside.jl), so all three profiles can be
# compared at the same h and the same centroid (each paired with its own conventional run).
# Monotonicity: neighbour step resolved when |paired mean| > 3 paired SEM at τ_b; sign recorded, not assumed.
# Run: SHAPE=linear|exponential julia --project=research/baseline -t 8 research/phase3/p3_4_8e_profiles_outside.jl

include(joinpath(@__DIR__, "harness.jl"))
const SHAPE = get(ENV, "SHAPE", "linear")
const OUT = phase3_outdir("p3_4_8e_profiles_outside_" * SHAPE)
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const RHO = 0.02
const S = 3.0
const HS = SHAPE == "linear" ? [0.1, 0.2, 0.3, 0.5, 0.7, 1.0] : [0.05, 0.1, 0.2, 0.3, 0.5, 1.0]
const θ = RHO * sqrt(π / D_REF)
const TAUS = [1e-2, 2.5e-3]
centroid(h) = SHAPE == "linear" ? h / 3 : (c = (S - 2R_REF) / h; h * (1 - (c + 1) * exp(-c)) / -expm1(-c))
paired(a, b) = (d = a.S[:, end] .- b.S[:, end]; (mean(d), std(d) / sqrt(length(d))))
models = vcat([("conventional", Cylinders(radius=R_REF, repeats=[S, S], surface_relaxation=θ))],
              [("$(SHAPE)_h$(h)", Cylinders(radius=R_REF, repeats=[S, S], layer_rho_outer=RHO, layer_h_outer=h, layer_shape=SHAPE)) for h in HS])
runs = Dict((n, τ) => run_cylinders(g; region=:outside, sx=S, sy=S, timestep=τ, nspins=NSPINS) for τ in TAUS for (n, g) in models)
identical = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs[("conventional", τ)].finals, runs[(n, τ)].finals))
                for τ in TAUS for (n, _) in models)
@printf("identical trajectories across models (per τ): %s\n", identical)
E50 = RHO * sv_outside(R_REF, S, S) * 50
rows = NamedTuple[]
for τ in TAUS, (k, (n, _)) in enumerate(models)
    r = runs[(n, τ)]
    (dm, ds) = paired(r, runs[("conventional", τ)])
    h = k == 1 ? 0.0 : HS[k - 1]
    push!(rows, (model=n, h=h, centroid=k == 1 ? 0.0 : centroid(h), tau=τ, S50=r.mean[end], sem=r.sem[end], J50=log(r.mean[end]) + E50,
        neglog50=r.neglog_mean[end], neglog_sem=r.neglog_sem[end], neglog_expected=E50, dS_vs_conv=dm, dS_vs_conv_sem=ds))
    @printf("%-18s τ = %.1e: S(50) %.6f ± %.1e, E[−ln M] %.5f ± %.1e (exp. %.5f), ΔS vs conv %+.2e ± %.1e\n",
        n, τ, r.mean[end], r.sem[end], r.neglog_mean[end], r.neglog_sem[end], E50, dm, ds)
end
names = [n for (n, _) in models]
steps = NamedTuple[]
for k in 1:length(names) - 1
    (ma, sa) = paired(runs[(names[k + 1], TAUS[1])], runs[(names[k], TAUS[1])])
    (mb, sb) = paired(runs[(names[k + 1], TAUS[2])], runs[(names[k], TAUS[2])])
    push!(steps, (from=names[k], to=names[k + 1], d_tau_a=ma, sem_tau_a=sa, d_tau_b=mb, sem_tau_b=sb, z_b=mb / sb,
        resolved=abs(mb / sb) > 3, sign=sign(mb), z_tau_shift=(ma - mb) / sqrt(sa^2 + sb^2)))
    @printf("%-18s → %-18s: τ_a %+.2e ± %.1e, τ_b %+.2e ± %.1e (z %+.1f), τ shift z %+.2f\n",
        names[k], names[k + 1], ma, sa, mb, sb, mb / sb, (ma - mb) / sqrt(sa^2 + sb^2))
end
layer_steps = steps[2:end]
signs = unique(s.sign for s in layer_steps if s.resolved)
write_csv(joinpath(OUT, "sweep.csv"), rows)
write_csv(joinpath(OUT, "steps.csv"), steps)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "shape" => SHAPE, "identical_trajectories" => identical,
    "taus" => TAUS, "nspins" => NSPINS, "n_layer_steps" => length(layer_steps), "n_resolved" => count(s -> s.resolved, layer_steps),
    "resolved_signs" => signs, "max_abs_tau_shift_z" => maximum(abs(s.z_tau_shift) for s in steps)))
println("saved to $OUT")
```

- [ ] **Step 2: Run** (background, two runs of ≈ 1.5 h)

```
SHAPE=linear      julia --project=research/baseline -t 8 research/phase3/p3_4_8e_profiles_outside.jl
SHAPE=exponential julia --project=research/baseline -t 8 research/phase3/p3_4_8e_profiles_outside.jl
```

Expected:
- identical trajectories;
- E[−ln M] = ρ(S/V)t within noise for every h, a re-check of P3.4.3 for these profiles outside;
- the neighbour steps' signs recorded.

If a step's τ shift has |z| > 3, report τ_b with the C1 caveat; a τ/16 pass is optional (ledger ruling with cost).

- [ ] **Step 3: Commit**

```bash
git add research/phase3/p3_4_8e_profiles_outside.jl research/results/phase3/p3_4_8e_profiles_outside_*
git commit -m "P3.4.8 outside for linear and exponential: attenuation vs h in the square array s = 3.0, paired Monte Carlo at two τ

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 9: Record and tracker

- [ ] **Step 1: Record sections** in `research/phase3_record.md` (Design / Method / Result / Conclusion; numbers only from the result files):

1. **"Linear and exponential profiles on cylinders (P3.3, D3.3, D3.5)."**
   - The integral (ρ/h)/κ·dt·∫g(d(s)/h)ds, d(s) = R ∓ √(as² + bs + q).
   - The sub-intervals: shell roots and s*.
   - Gauss–Legendre 8 with bisection, LAYER_QUAD_TOL.
   - The checks: linear closed form ∫√(quadratic), the wall regression, brute force.
   - Supports and validation; bit-identical step.
2. **"P3.3.4/P3.6.3: quadrature tolerance and cost"**: the tables, and the default tolerance as fixed.
3. **"P3.4.1/P3.4.3 for linear and exponential"**: the trace and identity tables.
4. **"P3.4.6 for linear and exponential"**: the comparison table.
5. **"P3.4.8 inside for all profiles"**: the h sweep, monotonicity and the same-centroid comparison, with the figure.
   - **Outside (Task 8b):** paired ΔS vs conventional and the neighbour steps for linear and exponential in s = 3.0, compared with the step profile's s = 3.0 run at the same h and the same centroid (each profile paired with its own conventional run; same seeds, spins and τ).
   - Conclusion: does the Phase 2 centroid collapse hold on cylinders? With which spread?
6. **Summary rows.**
7. **The radial-continuum minor** fixed in Task 3 (thin linear layers now conserve ∫V dA exactly).

- [ ] **Step 2: Tracker.** Tick P3.3.1–P3.3.4 and P3.6.3 if the criteria are met (tick only when met, otherwise ◐ with the reason). Add notes "(linear, exponential)" to P3.4.1, P3.4.3, P3.4.6 and P3.4.8 (inside and outside), with outcomes. Update the Phase 3 done count.

- [ ] **Step 3: Suites**

```
julia --project=research/baseline -t 8 research/phase3/test/runtests.jl
julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["no-plots"])'
```
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add research/phase3_record.md research/project-progress-tracker.md
git commit -m "Phase 3 record: linear and exponential on cylinders — quadrature, tolerance, trace, identity, continuum, h and centroid comparison; tracker

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## After the last task

Final whole-branch review by a fresh subagent on the most capable model (user rule). Give it this plan, the spec pointers, the Review Focus list and the ledger rulings. Grade the findings; fix Critical/Important ones with RED→GREEN tests; ledger the minors.
