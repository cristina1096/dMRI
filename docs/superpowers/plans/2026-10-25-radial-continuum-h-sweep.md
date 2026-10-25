# Radial continuum benchmark and attenuation vs h on cylinders (P3.4.6, P3.4.8) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:**
- Build an exact radial continuum benchmark for a single cylinder (inside), so that S(t) for the step layer and for the conventional model can be computed without Monte Carlo and without a timestep. Validate it against exact eigenvalues and against Monte Carlo.
- Record attenuation against h at fixed ρ on cylinders:
  - inside, from the continuum;
  - outside, in a square array, from paired Monte Carlo.

**Architecture:**
- **Continuum.** A finite-volume discretisation in r of ∂m/∂t = D(1/r)∂_r(r∂_r m) − ΔR₂(r)m on rings, solved through a symmetric eigenproblem. It is research code only (`research/phase3/continuum_radial.jl`) and independent of the package's layer code; κ is computed by quadrature.
- **Monte Carlo.** Uses the existing harness `run_cylinders`. All models and all h of one comparison share one timestep, so their trajectories are identical and their differences are paired per seed.

**Tech Stack:** Julia 1.12, research environment `research/baseline` (MCMRSimulator, SpecialFunctions, CairoMakie, LinearAlgebra).

**Spec:**
- `research/project-progress-tracker.md`, Phase 3: P3.4.6 and P3.4.8.
- `research/phase3_record.md`: D3.1 (ρ per unit surface area, κ = 1 − c̄/R inside), C-B3/P3.1.6 (θ = ρ√(π/D) holds on cylinders), and the lesson of Review C1 in `research/phase2_record.md`: at fixed ρ, profile differences are as small as the timestep error, so τ must be set by the error in the difference.

## Global Constraints

- **Branch.** Work only on branch `cc/near-surface-r2`. Never modify or commit on other branches. Do not push; the user pushes.
- **Commits.** Every commit message ends with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- **Package source.** No change to package source (`src/`). If a package bug is found, stop that task, ledger it and fix it test-first in `test/`.
- **Profiles.** Step profile only for the Monte Carlo (cylinders support only step). The continuum also accepts `"linear"`, for later use.
- **Fixed ρ.** ρ = 0.02 µm/ms throughout, so that ρ·(S/V)·t is the same for every h. The conventional model uses θ = ρ√(π/D) (P3.1.6).
- **Model framing.**
  - The conventional model is a reference, not the truth.
  - No fast-diffusion or Brownstein–Tarr assumptions. Bessel-function eigenvalues are used only as exact tests of the solver's lowest mode, never as a model of the signal.
  - Monotonicity: record the sign; do not assume it.
- **Records.** Design / Method / Result / Conclusion with formulas; figures where useful.
- **Research tests.** `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`
- **Research scripts.** `julia --project=research/baseline -t 8 research/phase3/<script>.jl`. Each writes `research/results/phase3/<name>/summary.json` with `provenance()`. Long runs go in the background, in chunks of ≤ 2 h.
- **Plan end.** A fresh subagent reviewer (most capable model) reviews the whole branch (user rule).

## Review Focus

1. **Thin layers below the grid spacing** (h < dr). The step's exact ring-overlap average must still integrate to 2πRρ and give a finite S that satisfies the Jensen bound S ≥ e^(−ρ(S/V)t). Test in Task 1 (Jensen bound, thin h).
2. **h = R** (the shell reaches the axis, r_a = 0). The uniform rate must be exact. Test in Task 1.
3. **Strong relaxation** (ρ = 0.5, so q ≈ 5 ms⁻¹). The exact matching must still bracket the root and the solver must agree. Test in Task 1.
4. **Pairing assumption.** The conventional model and the layer, run with the same seeds and τ, must have identical trajectories, otherwise paired SEMs are invalid. The scripts assert identical final positions in Tasks 2 and 4.
5. **Common τ across h.** Every h of a sweep must use the same τ (set by the thinnest layer), or pairing breaks. Asserted in Tasks 2 and 4.

---

## File Structure

- Create `research/phase3/continuum_radial.jl`:
  - `radial_g`, `radial_kappa`;
  - `radial_operator`, `radial_rates`, `radial_signal`;
  - `exact_rate_conventional`, `exact_rate_step`.
- Create `research/phase3/test/test_continuum_radial.jl`. Modify `research/phase3/test/runtests.jl` to include `continuum_radial.jl`.
- Create the scripts:
  - `research/phase3/p3_4_6_continuum_check.jl` → `research/results/phase3/p3_4_6_continuum_check/`;
  - `research/phase3/p3_4_8_h_sweep_inside.jl` → `research/results/phase3/p3_4_8_h_sweep_inside/`;
  - `research/phase3/p3_4_8_h_sweep_outside.jl` → `research/results/phase3/p3_4_8_h_sweep_outside/`.
- Modify `research/phase3_record.md` and `research/project-progress-tracker.md`.

---

### Task 1: Radial continuum solver with exact tests (P3.4.6, solver part)

**Files:**
- Create: `research/phase3/continuum_radial.jl`, `research/phase3/test/test_continuum_radial.jl`
- Modify: `research/phase3/test/runtests.jl`

**Interfaces:**
- Produces:
  - `radial_signal(; R, D, times, layer=nothing, surface=0.0, R2_bulk=0.0, N=2000)::Vector{Float64}`, with `layer = (shape::String, ρ, h)` and `surface = ρ` (conventional);
  - `radial_rates(; R, D, layer, surface, N)::Vector{Float64}` (ascending);
  - `radial_kappa(shape, h, R)`;
  - `exact_rate_conventional(R, D, ρ)`, `exact_rate_step(R, D, ρ, h)`.

- [ ] **Step 1: Write the failing tests**

`research/phase3/test/runtests.jl`: after the `include` of `harness.jl`, add

```julia
include(joinpath(@__DIR__, "..", "continuum_radial.jl"))
```

`research/phase3/test/test_continuum_radial.jl`:

```julia
@testset "radial continuum: exact cases" begin
    R, D, ρ = 1.0, 3.0, 0.02
    ts = [0.0, 10.0, 50.0]

    @testset "free diffusion and bulk only" begin
        @test all(isapprox.(radial_signal(R=R, D=D, times=ts), 1.0; atol=1e-12))
        @test all(isapprox.(radial_signal(R=R, D=D, times=ts, R2_bulk=1 / 80), exp.(-ts ./ 80); rtol=1e-12))
    end

    @testset "κ by quadrature (D3.1)" begin
        @test radial_kappa("step", 0.2, 1.0) ≈ 0.9 rtol=1e-8
        @test radial_kappa("step", 1.0, 1.0) ≈ 0.5 rtol=1e-8
        @test radial_kappa("linear", 0.3, 1.0) ≈ 1 - 0.1 rtol=1e-8      # c̄ = h/3
    end

    @testset "Review Focus 2: h = R gives the exact uniform rate 2ρ/R" begin
        @test all(isapprox.(radial_signal(R=R, D=D, times=ts, layer=("step", ρ, R)), exp.(-2ρ / R .* ts); rtol=1e-10))
    end

    @testset "initial slope −ρ·(S/V) = −2ρ/R for every model (uniform start)" begin
        δ = 1e-6
        for kw in ((layer=("step", ρ, 0.2),), (layer=("linear", ρ, 0.3),), (surface=ρ,))
            S = radial_signal(; R=R, D=D, times=[δ], kw...)[1]
            @test (S - 1) / δ ≈ -2ρ / R rtol=1e-3
        end
    end

    @testset "lowest rate: conventional (Robin) vs D k² with D k J₁(kR) = ρ J₀(kR)" begin
        for (Rr, Dd, ρρ) in ((1.0, 3.0, 0.02), (0.5, 1.0, 0.1), (2.0, 3.0, 0.5))
            @test radial_rates(R=Rr, D=Dd, surface=ρρ, N=4000)[1] ≈ exact_rate_conventional(Rr, Dd, ρρ) rtol=1e-3
        end
    end

    @testset "lowest rate: step layer vs exact Bessel J/I/K matching (Review Focus 3: strong ρ)" begin
        for (Rr, Dd, ρρ, h) in ((1.0, 3.0, 0.02, 0.2), (1.0, 3.0, 0.5, 0.1), (0.5, 1.0, 0.05, 0.25))
            @test radial_rates(R=Rr, D=Dd, layer=("step", ρρ, h), N=2000)[1] ≈ exact_rate_step(Rr, Dd, ρρ, h) rtol=1e-5
        end
    end

    @testset "Review Focus 1: Jensen bounds, including h below the grid spacing" begin
        E = 2ρ / R * 50
        for h in (1e-4, 0.01, 0.2, 0.7)
            S = radial_signal(R=R, D=D, times=[50.0], layer=("step", ρ, h), N=1000)[1]     # dr = 1e-3
            @test isfinite(S)
            @test S >= exp(-E) * (1 - 1e-9)                       # J ≥ 0 (Jensen, a theorem; the ordering vs the conventional model is a result, not tested)
        end
    end

    @testset "grid convergence" begin
        a = radial_signal(R=R, D=D, times=[50.0], layer=("step", ρ, 0.2), N=1000)[1]
        b = radial_signal(R=R, D=D, times=[50.0], layer=("step", ρ, 0.2), N=4000)[1]
        @test abs(a / b - 1) < 1e-6
        a = radial_signal(R=R, D=D, times=[50.0], surface=ρ, N=1000)[1]
        b = radial_signal(R=R, D=D, times=[50.0], surface=ρ, N=4000)[1]
        @test abs(a / b - 1) < 1e-4
    end
end
```

- [ ] **Step 2: Run to verify they fail**

Run: `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`
Expected: FAIL with "could not open file … continuum_radial.jl".

- [ ] **Step 3: Implement `research/phase3/continuum_radial.jl`**

```julia
# Radial continuum benchmark (P3.4.6): exact signal inside a single cylinder for the near-surface layer and for the
# conventional surface relaxation, without Monte Carlo or timestep. Independent of the package's layer code.
#
#     ∂m/∂t = D (1/r) ∂_r(r ∂_r m) − [R2_bulk + ΔR2(r)] m,   0 ≤ r ≤ R,   uniform start m = 1,
#     layer:        ∂_r m = 0 at r = R,  ΔR2(r) = (ρ/h)·g((R − r)/h)/κ,  κ = ∫g(u)(1 − u h/R) du   (D3.1; κ by quadrature),
#     conventional: D ∂_r m = −ρ m at r = R (Robin), ρ = θ√(D/π) (P3.1.6).
# Finite volumes: N rings of width dr = R/N, areas A_i = π(r_{i+1/2}² − r_{i−1/2}²), face conductances 2πD r_{i+1/2}/dr,
# the Robin sink 2πRρ on the last ring; ΔR2 averaged over each ring (step: exact area overlap; other profiles:
# r-weighted sub-sampling). A ṁ = −(K + A V) m with y = √A m becomes ẏ = −M y, M = A^(−1/2)(K + A V)A^(−1/2) symmetric
# tridiagonal = Σ λ_k u_k u_kᵀ, so S(t) = Σ_i A_i m_i / Σ A_i = Σ_k (u_kᵀ√A)² e^(−λ_k t) / Σ A_i.

using LinearAlgebra
import SpecialFunctions: besselj0, besselj1, besseli, besselk

"Normalised profile g(u), ∫g = 1 on [0, 1] (written independently of the package)."
function radial_g(shape::AbstractString, u::Real)
    shape == "step" && return 0 <= u <= 1 ? 1.0 : 0.0
    shape == "linear" && return 0 <= u <= 1 ? 2 * (1 - u) : 0.0
    error("radial_g: profile $shape is not supported")
end

"κ = ∫g(u)(1 − u h/R) du inside a cylinder (D3.1), by midpoint quadrature."
radial_kappa(shape, h, R; n=100_000) = sum(radial_g(shape, (k - 0.5) / n) * (1 - (k - 0.5) / n * h / R) for k in 1:n) / n

"Symmetric operator M (1/ms) and ring areas A for the radial problem (see the file header)."
function radial_operator(; R, D, layer=nothing, surface=0.0, N=2000, subsamples=64)
    dr = R / N
    redge = collect(0:N) .* dr
    A = [π * (redge[i + 1]^2 - redge[i]^2) for i in 1:N]
    V = zeros(N)
    if !isnothing(layer)
        (shape, ρ, h) = layer
        κ = radial_kappa(shape, h, R)
        ra = max(R - h, 0.0)
        for i in 1:N
            if shape == "step"
                lo, hi = max(redge[i], ra), redge[i + 1]
                V[i] = hi > lo ? ρ / (h * κ) * π * (hi^2 - lo^2) / A[i] : 0.0
            else
                s = 0.0
                w = 0.0
                for j in 1:subsamples
                    r = redge[i] + (j - 0.5) / subsamples * dr
                    s += r * radial_g(shape, (R - r) / h)
                    w += r
                end
                V[i] = ρ / (h * κ) * s / w
            end
        end
    end
    c = [2π * D * redge[i + 1] / dr for i in 1:N - 1]
    K = zeros(N)
    for i in 1:N - 1
        K[i] += c[i]
        K[i + 1] += c[i]
    end
    K[N] += 2π * R * surface
    M = SymTridiagonal((K .+ A .* V) ./ A, -c ./ sqrt.(A[1:N - 1] .* A[2:N]))
    return (M, A)
end

"Decay rates λ_k (1/ms, ascending) of the discretised radial operator."
radial_rates(; kw...) = eigvals(radial_operator(; kw...)[1])

"""
    radial_signal(; R, D, times, layer=nothing, surface=0.0, R2_bulk=0.0, N=2000)

S(t)/S(0) inside a cylinder of radius R for a uniform start. `layer = (shape, ρ, h)` (ρ per unit surface area, D3.1),
`surface = ρ` for the conventional model.
"""
function radial_signal(; R, D, times, layer=nothing, surface=0.0, R2_bulk=0.0, N=2000)
    (M, A) = radial_operator(; R=R, D=D, layer=layer, surface=surface, N=N)
    F = eigen(M)
    c2 = (F.vectors' * sqrt.(A)) .^ 2 ./ sum(A)
    return [sum(c2 .* exp.(-F.values .* t)) * exp(-R2_bulk * t) for t in times]
end

"Exact lowest rate with surface relaxivity ρ (test only): D k², D k J₁(kR) = ρ J₀(kR), k ∈ (0, j₀₁/R), by bisection."
function exact_rate_conventional(R, D, ρ)
    f(k) = D * k * besselj1(k * R) - ρ * besselj0(k * R)
    lo, hi = 0.0, 2.404825557695773 / R
    for _ in 1:200
        mid = (lo + hi) / 2
        f(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    return D * ((lo + hi) / 2)^2
end

"""
    exact_rate_step(R, D, ρ, h)

Exact lowest rate for the step layer inside (test only), 0 < h < R. Core r < r_a = R − h: m = J₀(kr), λ = Dk². Shell
(rate q = ρ/(h(1 − h/2R))): m = I₀(βr) + (I₁(βR)/K₁(βR))K₀(βr), β = √((q − λ)/D), which satisfies ∂_r m(R) = 0. The
lowest λ lies in (0, q) (λ₁ ≤ mean rate 2ρ/R < q) and makes m'/m continuous at r_a; found by bisection.
"""
function exact_rate_step(R, D, ρ, h)
    q = ρ / (h * (1 - h / (2R)))
    ra = R - h
    function F(λ)
        k = sqrt(λ / D)
        β = sqrt((q - λ) / D)
        Lin = -k * besselj1(k * ra) / besselj0(k * ra)
        ratio = besseli(1, β * R) / besselk(1, β * R)
        Lsh = β * (besseli(1, β * ra) - ratio * besselk(1, β * ra)) / (besseli(0, β * ra) + ratio * besselk(0, β * ra))
        return Lin - Lsh
    end
    lo, hi = 0.0, q * (1 - 1e-12)
    @assert F(lo) > 0 && F(hi) < 0 "exact_rate_step: root not bracketed"
    @assert sqrt(q / D) * ra < 2.404825557695773 "exact_rate_step: core J₀ would change sign"
    for _ in 1:200
        mid = (lo + hi) / 2
        F(mid) > 0 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `julia --project=research/baseline -t 8 research/phase3/test/runtests.jl`
Expected: PASS (23 harness tests + the new continuum tests).

Two failure cases to recognise:
- **The step-rate test fails with a constant relative offset.** Check κ: the continuum uses quadrature, the exact rate uses 1 − h/2R. A mismatch above 1e-8 means the quadrature is wrong.
- **The conventional test is off by about dr/R.** That is the first-order Robin cell error. Raise N in the test only if the error falls as 1/N, and ledger it.

- [ ] **Step 5: Commit**

```bash
git add research/phase3/continuum_radial.jl research/phase3/test/
git commit -m "P3.4.6: radial continuum benchmark for a single cylinder, tested against exact Bessel eigenvalues

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Monte Carlo vs continuum, inside a single cylinder (P3.4.6, validation)

**Files:**
- Create: `research/phase3/p3_4_6_continuum_check.jl` → `research/results/phase3/p3_4_6_continuum_check/` (`comparison.csv`, `summary.json`, `p3_4_6_continuum_check.png`)

**Interfaces:**
- Consumes: `radial_signal` (Task 1); `run_cylinders`, `rate_inside`, `phase3_outdir` (harness).
- Produces: the validation verdict, and per-seed Monte Carlo S(t) for the step layer h ∈ {0.02, 0.05, 0.2, 0.5, 1.0} and for the conventional model. Task 3 cites these as the Monte Carlo points on the continuum curves.

- [ ] **Step 1: Write the script**

```julia
# P3.4.6 validation: Monte Carlo vs the radial continuum, inside a single cylinder R = 1, D = 3, ρ = 0.02 (fixed ρ),
# step layer (inner side only) h = 0.02, 0.05, 0.2, 0.5, 1.0 and the conventional model θ = ρ√(π/D).
# All models share one timestep per pass, τ_a = min(1e-2, 0.005/ΔR2(0) of the thinnest layer), then τ_b = τ_a/4, so the
# trajectories are identical across models (asserted) and ΔS = S_layer − S_conv is paired per seed (lesson of Review C1:
# at fixed ρ the model difference is as small as the timestep error, so it is checked at two τ).
# Seeds 1–10, NSPINS (default 20000) spins per seed in the box [−R, R]³ (≈ π/4 inside), readouts 0:5:50.
# Criteria (tracker P3.4.6): at τ_b, max over readouts |S_MC − S_cont|/SEM ≤ 3 for every model, and
# |ΔS_MC − ΔS_cont|/SEM_paired ≤ 3 for every layer (3, not 2: 11 correlated readouts).
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_6_continuum_check.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_4_6_continuum_check")
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const RHO = 0.02
const HS = [0.02, 0.05, 0.2, 0.5, 1.0]
const θ = RHO * sqrt(π / D_REF)
const τa = min(1e-2, 0.005 / rate_inside(RHO, minimum(HS), R_REF))
const TAUS = [τa, τa / 4]
bb = mr.BoundingBox([-R_REF, -R_REF, -R_REF], [R_REF, R_REF, R_REF])
models = vcat([("conventional", Cylinders(radius=R_REF, surface_relaxation=θ), nothing)],
              [("step_h$(h)", Cylinders(radius=R_REF, layer_rho_inner=RHO, layer_h_inner=h), h) for h in HS])
cont = Dict(name => (isnothing(h) ? radial_signal(R=R_REF, D=D_REF, times=TIMES_REF, surface=RHO) :
                                    radial_signal(R=R_REF, D=D_REF, times=TIMES_REF, layer=("step", RHO, h)))
            for (name, _, h) in models)

rows = NamedTuple[]
runs = Dict{Tuple{String, Float64}, Any}()
for τ in TAUS, (name, g, h) in models
    runs[(name, τ)] = run_cylinders(g; region=:inside, bounding_box=bb, timestep=τ, nspins=NSPINS)
    @assert runs[(name, τ)].tau == τ
end
identical = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs[("conventional", τ)].finals, runs[(name, τ)].finals))
                for τ in TAUS for (name, _, _) in models)
@printf("identical trajectories across models (per τ): %s\n", identical)
allpass = identical
for τ in TAUS, (name, _, h) in models
    r = runs[(name, τ)]
    zS = (r.mean .- cont[name]) ./ max.(r.sem, 1e-15)
    if isnothing(h)
        dmean = dsem = dcont = zeros(length(TIMES_REF)); zD = zeros(length(TIMES_REF))
    else
        d = r.S .- runs[("conventional", τ)].S
        dmean = vec(mean(d, dims=1)); dsem = vec(std(d, dims=1)) ./ sqrt(size(d, 1))
        dcont = cont[name] .- cont["conventional"]
        zD = (dmean .- dcont) ./ max.(dsem, 1e-15)
    end
    pass = maximum(abs.(zS[2:end])) <= 3 && maximum(abs.(zD[2:end])) <= 3
    τ == TAUS[end] && (global allpass &= pass)
    push!(rows, (model=name, h=isnothing(h) ? 0.0 : h, tau=τ, S50_mc=r.mean[end], S50_sem=r.sem[end], S50_cont=cont[name][end],
        maxabs_zS=maximum(abs.(zS[2:end])), dS50_mc=dmean[end], dS50_sem=dsem[end], dS50_cont=dcont[end],
        maxabs_zdS=maximum(abs.(zD[2:end])), pass=pass))
    @printf("%-14s τ = %.2e: S(50) MC %.6f ± %.1e vs %.6f (max|z| %.2f); ΔS(50) MC %+.2e ± %.1e vs %+.2e (max|z| %.2f)\n",
        name, τ, r.mean[end], r.sem[end], cont[name][end], maximum(abs.(zS[2:end])), dmean[end], dsem[end], dcont[end], maximum(abs.(zD[2:end])))
end
write_csv(joinpath(OUT, "comparison.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "pass" => allpass, "identical_trajectories" => identical,
    "taus" => TAUS, "nspins" => NSPINS, "rho" => RHO, "theta" => θ))
fig = Figure(size=(700, 420))
ax = Axis(fig[1, 1], xlabel="t (ms)", ylabel="ΔS = S_layer − S_conv", title="inside R = 1, ρ = 0.02: continuum (lines), MC τ_b (points)")
for (name, _, h) in models[2:end]
    lines!(ax, TIMES_REF, cont[name] .- cont["conventional"], label="h = $h")
    d = runs[(name, TAUS[end])].S .- runs[("conventional", TAUS[end])].S
    errorbars!(ax, TIMES_REF, vec(mean(d, dims=1)), vec(std(d, dims=1)) ./ sqrt(size(d, 1)))
    scatter!(ax, TIMES_REF, vec(mean(d, dims=1)), markersize=6)
end
axislegend(ax, position=:lb, labelsize=9)
save(joinpath(OUT, "p3_4_6_continuum_check.png"), fig)
println(allpass ? "PASS" : "FAIL")
println("saved to $OUT")
```

- [ ] **Step 2: Run (background, ≈ 1 h)**

Run: `julia --project=research/baseline -t 8 research/phase3/p3_4_6_continuum_check.jl`
Expected:
- "identical trajectories … true";
- at τ_b, every model has max|z| ≤ 3 for S and ΔS → "PASS";
- at τ_a, the thinnest layers may show a ΔS offset (Review C1). It is recorded, not graded.

If τ_b fails for a thin layer while thicker ones pass, rerun that layer and the conventional model at τ_a/16 with the same script logic (ledger a ruling). Do not change the criterion.

- [ ] **Step 3: Commit**

```bash
git add research/phase3/p3_4_6_continuum_check.jl research/results/phase3/p3_4_6_continuum_check
git commit -m "P3.4.6: Monte Carlo vs radial continuum inside a cylinder, paired ΔS at two τ

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Attenuation vs h inside a cylinder, from the continuum (P3.4.8 inside)

**Files:**
- Create: `research/phase3/p3_4_8_h_sweep_inside.jl` → `research/results/phase3/p3_4_8_h_sweep_inside/` (`sweep.csv`, `summary.json`, `p3_4_8_h_sweep_inside.png`)

**Interfaces:**
- Consumes: `radial_signal` (Task 1); `phase3_outdir`, `write_csv`, `write_json` (harness and common).
- Produces: `sweep.csv` with columns `R, h, h_over_R, S50, J50, dS50_vs_conv, curvature, grid_err`.

- [ ] **Step 1: Write the script**

```julia
# P3.4.8 (inside): attenuation vs h at fixed ρ inside single cylinders, from the radial continuum (validated in P3.4.6).
# R = 0.5, 1, 2; D = 3; ρ = 0.02 (so E[−ln M] = ρ(2/R)t is the same for every h); step; h/R = 0.005 … 1.
# Quantities at readouts 0:5:50:
#   S(50); Jensen gap J(50) = ln S(50) + ρ(2/R)·50 ≥ 0 (Phase 2 record, M3: S = e^(−E[X]) e^J);
#   ΔS(50) = S − S_conv (conventional: Robin, same ρ); curvature c of ln S = a + bt + ct² (least squares);
#   grid error = |S_N=2000 − S_N=4000|.
# Monotonicity: the sign of S(50; h_{k+1}) − S(50; h_k) is recorded for every neighbouring pair; a step counts as resolved
# when |difference| > 10 × grid error. h → 0 approaches the conventional model; h = R is uniform (J = 0).
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_8_h_sweep_inside.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_4_8_h_sweep_inside")
const RHO = 0.02
const HR = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]
curv(ts, S) = (X = hcat(ones(length(ts)), ts, ts .^ 2); (X \ log.(S))[3])

rows = NamedTuple[]
steps = NamedTuple[]
conv = Dict{Float64, Float64}()
for R in (0.5, 1.0, 2.0)
    E50 = RHO * (2 / R) * 50
    Sc = radial_signal(R=R, D=D_REF, times=TIMES_REF, surface=RHO)
    conv[R] = Sc[end]
    prev = nothing
    for hr in HR
        h = hr * R
        S = radial_signal(R=R, D=D_REF, times=TIMES_REF, layer=("step", RHO, h))
        S4 = radial_signal(R=R, D=D_REF, times=[50.0], layer=("step", RHO, h), N=4000)[1]
        row = (R=R, h=h, h_over_R=hr, S50=S[end], J50=log(S[end]) + E50, dS50_vs_conv=S[end] - Sc[end],
               curvature=curv(TIMES_REF, S), grid_err=abs(S[end] - S4))
        push!(rows, row)
        if !isnothing(prev)
            d = row.S50 - prev.S50
            push!(steps, (R=R, h_from=prev.h, h_to=h, dS50=d, resolved=abs(d) > 10 * max(row.grid_err, prev.grid_err), sign=sign(d)))
        end
        prev = row
        @printf("R = %.1f h/R = %.3f: S(50) %.6f  J %.3e  ΔS vs conv %+.2e  curvature %+.1e  grid %.1e\n",
            R, hr, row.S50, row.J50, row.dS50_vs_conv, row.curvature, row.grid_err)
    end
    @printf("R = %.1f conventional: S(50) %.6f  J %.3e\n", R, Sc[end], log(Sc[end]) + E50)
end
write_csv(joinpath(OUT, "sweep.csv"), rows)
write_csv(joinpath(OUT, "steps.csv"), steps)
signs = unique(s.sign for s in steps if s.resolved)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rho" => RHO, "D" => D_REF,
    "conventional_S50" => conv, "n_steps" => length(steps), "n_resolved" => count(s -> s.resolved, steps),
    "resolved_signs" => signs, "monotone" => length(signs) == 1))
@printf("neighbour steps resolved: %d of %d; signs of resolved steps: %s\n", count(s -> s.resolved, steps), length(steps), signs)

fig = Figure(size=(1100, 400))
ax1 = Axis(fig[1, 1], xlabel="h/R", ylabel="J(50) = ln S + ρ(S/V)t", xscale=log10, title="Jensen gap (fixed ρ = 0.02)")
ax2 = Axis(fig[1, 2], xlabel="h/R", ylabel="ΔS(50) = S − S_conv", xscale=log10, title="difference from the conventional model")
for (R, col) in zip((0.5, 1.0, 2.0), (:steelblue, :darkorange, :seagreen))
    q = filter(r -> r.R == R, rows)
    scatterlines!(ax1, [r.h_over_R for r in q], [r.J50 for r in q], color=col, label="R = $R")
    hlines!(ax1, [log(conv[R]) + RHO * (2 / R) * 50], color=col, linestyle=:dash)
    scatterlines!(ax2, [r.h_over_R for r in q], [r.dS50_vs_conv for r in q], color=col, label="R = $R")
end
hlines!(ax2, [0.0], color=:black)
axislegend(ax1, labelsize=9)
save(joinpath(OUT, "p3_4_8_h_sweep_inside.png"), fig)
println("saved to $OUT")
```

- [ ] **Step 2: Run (seconds)**

Run: `julia --project=research/baseline -t 8 research/phase3/p3_4_8_h_sweep_inside.jl`
Expected:
- 30 rows;
- J(50) ≥ 0 everywhere, falling to ≈ 0 at h/R = 1;
- conventional J above all layers;
- the resolved signs are reported. Record them as found (from Phase 2 one expects S(50) to fall with h at fixed ρ, but do not assert it).

- [ ] **Step 3: Commit**

```bash
git add research/phase3/p3_4_8_h_sweep_inside.jl research/results/phase3/p3_4_8_h_sweep_inside
git commit -m "P3.4.8 inside: attenuation vs h at fixed ρ in single cylinders (continuum), Jensen gap and ΔS vs conventional

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Attenuation vs h outside, square array, paired Monte Carlo (P3.4.8 outside)

**Files:**
- Create: `research/phase3/p3_4_8_h_sweep_outside.jl` → `research/results/phase3/p3_4_8_h_sweep_outside/` (`sweep.csv`, `steps.csv`, `summary.json`, `p3_4_8_h_sweep_outside.png`)

**Interfaces:**
- Consumes: `run_cylinders`, `rate_outside`, `sv_outside` (harness).
- Produces: `sweep.csv` (S(t) per h, τ) and `steps.csv` (paired neighbour differences and their convergence between τ).

- [ ] **Step 1: Write the script**

```julia
# P3.4.8 (outside): attenuation vs h at fixed ρ in the square array s = 2.5, R = 1 (gap 0.5), outer layer only,
# h = 0.02, 0.05, 0.1, 0.2 (separated), 0.3, 0.4, 0.5 (outer layers overlap: 2h > gap), plus the conventional model.
# No continuum exists outside (2D), so Monte Carlo with pairing: all h and the conventional model share one τ (set by
# the thinnest layer), hence identical trajectories (asserted) and per-seed paired differences between neighbouring h.
# Two passes τ_a = min(1e-2, 0.005/ΔR2(0) of h = 0.02) and τ_b = τ_a/4; the difference between the paired estimates
# of the two passes measures the timestep bias of the h dependence (Review C1).
# Seeds 1–10, NSPINS (default 20000) spins per seed in the 500 µm box (≈ 50 % outside), readouts 0:5:50.
# Monotonicity: a neighbour step is resolved when |paired mean| > 3 paired SEM at τ_b; its sign is recorded, not assumed.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_8_h_sweep_outside.jl

include(joinpath(@__DIR__, "harness.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_4_8_h_sweep_outside")
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const RHO = 0.02
const HS = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5]
const θ = RHO * sqrt(π / D_REF)
const τa = min(1e-2, 0.005 / rate_outside(RHO, minimum(HS), R_REF))
const TAUS = [τa, τa / 4]
models = vcat([("conventional", Cylinders(radius=R_REF, repeats=[S_REF, S_REF], surface_relaxation=θ))],
              [("step_h$(h)", Cylinders(radius=R_REF, repeats=[S_REF, S_REF], layer_rho_outer=RHO, layer_h_outer=h)) for h in HS])
runs = Dict{Tuple{String, Float64}, Any}()
for τ in TAUS, (name, g) in models
    runs[(name, τ)] = run_cylinders(g; region=:outside, sx=S_REF, sy=S_REF, timestep=τ, nspins=NSPINS)
    @assert runs[(name, τ)].tau == τ
end
identical = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs[("conventional", τ)].finals, runs[(name, τ)].finals))
                for τ in TAUS for (name, _) in models)
@printf("identical trajectories across models (per τ): %s\n", identical)

E50 = RHO * sv_outside(R_REF, S_REF, S_REF) * 50
rows = NamedTuple[]
for τ in TAUS, (name, _) in models
    r = runs[(name, τ)]
    push!(rows, (model=name, tau=τ, S50=r.mean[end], sem=r.sem[end], J50=log(r.mean[end]) + E50,
        neglog50=r.neglog_mean[end], neglog_expected=E50))
    @printf("%-14s τ = %.2e: S(50) %.6f ± %.1e, J %.3e, E[−ln M] %.5f (expected %.5f)\n", name, τ, r.mean[end], r.sem[end],
        log(r.mean[end]) + E50, r.neglog_mean[end], E50)
end
paired(a, b) = (d = a.S[:, end] .- b.S[:, end]; (mean(d), std(d) / sqrt(length(d))))
steps = NamedTuple[]
names = [n for (n, _) in models]
for k in 1:length(names) - 1
    (ma, sa) = paired(runs[(names[k + 1], TAUS[1])], runs[(names[k], TAUS[1])])
    (mb, sb) = paired(runs[(names[k + 1], TAUS[2])], runs[(names[k], TAUS[2])])
    push!(steps, (from=names[k], to=names[k + 1], d_tau_a=ma, sem_tau_a=sa, d_tau_b=mb, sem_tau_b=sb, z_b=mb / sb,
        resolved=abs(mb / sb) > 3, sign=sign(mb), z_tau_shift=(ma - mb) / sqrt(sa^2 + sb^2)))
    @printf("%-14s → %-14s: ΔS(50) τ_a %+.2e ± %.1e, τ_b %+.2e ± %.1e (z %+.1f), τ shift z %+.2f\n",
        names[k], names[k + 1], ma, sa, mb, sb, mb / sb, (ma - mb) / sqrt(sa^2 + sb^2))
end
signs = unique(s.sign for s in steps[2:end] if s.resolved)      # steps[1] is conventional → thinnest layer
write_csv(joinpath(OUT, "sweep.csv"), rows)
write_csv(joinpath(OUT, "steps.csv"), steps)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "identical_trajectories" => identical,
    "taus" => TAUS, "nspins" => NSPINS, "rho" => RHO, "n_layer_steps_resolved" => count(s -> s.resolved, steps[2:end]),
    "resolved_signs_layer_steps" => signs, "monotone" => length(signs) <= 1))
fig = Figure(size=(650, 420))
ax = Axis(fig[1, 1], xlabel="h (µm)", ylabel="S(50) − S_conv(50) (paired)", title="outside, square array s = 2.5, R = 1, ρ = 0.02")
for (τ, col) in zip(TAUS, (:gray, :black))
    d = [paired(runs[(n, τ)], runs[("conventional", τ)]) for n in names[2:end]]
    errorbars!(ax, HS, first.(d), last.(d), color=col)
    scatterlines!(ax, HS, first.(d), color=col, label="τ = $(round(τ, sigdigits=3)) ms")
end
vlines!(ax, [0.25], color=:red, linestyle=:dash, label="overlap starts (2h = gap)")
hlines!(ax, [0.0], color=:black, linewidth=0.5)
axislegend(ax, labelsize=9)
save(joinpath(OUT, "p3_4_8_h_sweep_outside.png"), fig)
println("saved to $OUT")
```

- [ ] **Step 2: Run (background, ≈ 2 h)**

Run: `julia --project=research/baseline -t 8 research/phase3/p3_4_8_h_sweep_outside.jl`
Expected:
- "identical trajectories … true";
- E[−ln M] equal to ρ(S/V)t within noise for every h, a free re-check of P3.4.3;
- paired neighbour steps with recorded signs.

If |z_tau_shift| > 3 for a step, the h dependence is timestep-biased at τ_a (as in Review C1). Then run a τ_a/16 pass for the two h of that step and the conventional model (ledger a ruling), and report τ_b, or τ_a/16 if they differ, as the value.

- [ ] **Step 3: Commit**

```bash
git add research/phase3/p3_4_8_h_sweep_outside.jl research/results/phase3/p3_4_8_h_sweep_outside
git commit -m "P3.4.8 outside: attenuation vs h at fixed ρ in the square array, paired Monte Carlo at two τ

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Record and tracker

**Files:**
- Modify: `research/phase3_record.md`, `research/project-progress-tracker.md`

- [ ] **Step 1: Record sections** (Design / Method / Result / Conclusion, numbers only from the result files)

1. **"P3.4.6: radial continuum benchmark"**:
   - the equation and boundary conditions;
   - the finite-volume scheme (ring areas, conductances, Robin sink), the symmetric eigenproblem and the S(t) formula;
   - the exact tests with their tolerances (the Bessel eigen-equations written out);
   - the Monte Carlo comparison table (S(50) and ΔS(50) at τ_a and τ_b, max |z|) and the figure;
   - the conclusion: whether the continuum stands in for Monte Carlo inside single cylinders, and the size of the τ_a bias in ΔS (lesson of C1).
2. **"P3.4.8: attenuation vs h at fixed ρ"**:
   - why fixed ρ: E[−ln M] = ρ(S/V)t is the same at every h, so only J changes (Phase 2 M3);
   - inside: a continuum table for R = 1 (all h/R), a summary for R = 0.5 and 2, and the figure;
   - outside: paired Monte Carlo table, τ shift and figure;
   - the monotonicity outcome with the signs as found;
   - the conclusion relating the h → 0 limit to the conventional model and h = R to the uniform case.
3. **Summary table rows** for P3.4.6 and P3.4.8.

- [ ] **Step 2: Tracker.** Tick P3.4.6 and P3.4.8 with outcomes, marking P3.4.8 "(step)". Update the Phase 3 done count.

- [ ] **Step 3: Run all suites**

```
julia --project=research/baseline -t 8 research/phase3/test/runtests.jl
julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["layer_cylinders", "layer"])'
```
Expected: PASS. The package is unchanged, so the targeted package suites suffice.

- [ ] **Step 4: Commit**

```bash
git add research/phase3_record.md research/project-progress-tracker.md
git commit -m "Phase 3 record: radial continuum benchmark (P3.4.6) and attenuation vs h at fixed ρ (P3.4.8); tracker

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

## After the last task

Final whole-branch review by a fresh subagent on the most capable model (user rule), over this plan's range. Give it the plan, the spec pointers above, this Review Focus list and the ledger rulings. Grade the findings; fix Critical/Important ones with RED→GREEN tests; ledger the minors.
