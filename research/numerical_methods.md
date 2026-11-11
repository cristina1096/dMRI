# Numerical methods reference: near-surface R₂ layer project

This is a single reference for every numerical method used in the project (Phases 2–3). For each method it gives:
- what it computes;
- the formula or algorithm;
- its accuracy, and how that was established;
- where it is implemented.

Results and decisions live in the records: [phase2_record.md](phase2_record.md) and [phase3_record.md](phase3_record.md). User documentation is in [docs/src/near_surface_layer.md](../docs/src/near_surface_layer.md).

**The model throughout:**

  R₂,total(x) = R₂,bulk + ΔR₂(d(x)),  ΔR₂(d) = (ρ/h)·g(d/h)/κ,

- d is the distance from the nearest surface, measured into the compartment.
- g is a normalised profile (∫g = 1 over its support).
- κ = 1 on walls and 1 ∓ c̄/R on cylinders.
- The conventional surface-relaxation model of MCMR (per-collision loss θ√τ) is a reference, not the truth.

**Contents**
1. Monte Carlo simulation (package)
2. Path integral of the layer rate along a straight step (package)
3. Periodic geometry: images, overlap and the timestep bound (package)
4. Continuum reference solvers (research)
5. Exact analytic references
6. Independent verification (brute force)
7. Statistics
8. Root finding, fitting and threshold location
9. Rules derived from the numerics

---

## 1. Monte Carlo simulation (package)

**What.** Transverse magnetisation of many spins diffusing in a geometry. This is MCMRSimulator's existing engine.

**Algorithm.**
- **Random walk.** Each spin takes Gaussian steps of RMS length √(2Dτ) per dimension and reflects specularly off surfaces. Every straight piece between reflections is a "piece".
- **Relaxation per piece.** The transverse magnetisation is multiplied by e^(−R₂,bulk·dt) and by the layer factor e^(−∫ΔR₂ dt), with the integral computed exactly or to quadrature tolerance along the piece (§2).
- **Conventional model.** Each collision multiplies by 1 − θ√τ (MCMR's `surface_relaxation`).
- **Signal.** S(t) is the mean of M⊥ over the spins of a compartment.
- **Restriction.** Instantaneous RF pulses only; a finite RF pulse with a layer raises an error.

**Where.**
- Layer factor applied per piece: [src/evolve.jl:494-498](../src/evolve.jl#L494); finite-RF check at [:492](../src/evolve.jl#L492).
- Research runners:
  - walls: [research/phase2/harness.jl](phase2/harness.jl);
  - cylinders (`run_cylinders`): [research/phase3/harness.jl:74](phase3/harness.jl#L74);
  - paired comparison runners: [p3_7_mc.jl](phase3/p3_7_mc.jl) (inside), [p3_7b_mc.jl](phase3/p3_7b_mc.jl) (outside, rectangular arrays).

**Timestep.** τ = min of several limits ([src/timesteps.jl](../src/timesteps.jl)). The layer limit is τ ≤ c/ΔR₂,max with c = 0.005 (§3); for comparison studies τ is set explicitly (§9).

---

## 2. Path integral of the layer rate along a straight step (package)

Along a straight piece p(s) = p₀ + s(p₁ − p₀), s ∈ [0, 1], of duration dt, the exponent is ∫ΔR₂(d(s)) dt.

### 2.1 Walls: closed form for every profile

**Why it is exact.** The distance to a wall is linear along the piece, d(s) = d₀ + s(d₁ − d₀).

**Method.**
- **Clip.** The piece is clipped to the support 0 ≤ d ≤ L (`segment_overlap`, the fraction inside).
- **Integrate exactly:**
  - step: ΔR₂(0)·dt·fraction;
  - linear: g is linear in s, so its mean over the part inside equals its value at the midpoint;
  - exponential: the mean of e^(−d/h) for d uniform on [d_a, d_b] is e^(−d_near/h)·(1 − e^(−x))/x, with x = |d_b − d_a|/h. It is evaluated with `expm1`, switching to the series 1 − x/2 for x < 1e-8, and from the end nearer the wall so it never computes 0·∞.

**Where.** `segment_overlap`, `segment_fraction`, `side_exponent`: [src/geometries/internal/layers.jl:80-125](../src/geometries/internal/layers.jl#L80).

**Accuracy.** Exact to roundoff: single-spin trace 2.4e-12 (linear) and 7.1e-13 (exponential), P2.2.9 / P2.4.

### 2.2 Cylinders, step profile: exact line–circle intersection

**Setup.** With u = p₀ − c and Δ = p₁ − p₀, the squared distance to the axis is quadratic in s:

  |p(s) − c|² = a s² + b s + q,  a = |Δ|², b = 2 u·Δ, q = |u|².

**Method.**
- The set {s : a s² + b s + q ≤ X} is an interval. Its ends are the roots, computed in the cancellation-free form t = −(b + sign(b)√disc)/2, s = t/a and (q − X)/t (`disk_interval`).
- The fraction inside the shell r_lo ≤ r ≤ r_hi is the length for X = r_hi² minus the length for X = r_lo² (`shell_fraction`).
- The exponent is the surface rate × dt × that fraction.
- Purely axial motion (a = 0) keeps r fixed.

**Where.** `disk_interval`, `disk_length`, `shell_fraction`, `cylinder_side_exponent`: [layers.jl:215-258, 303-321](../src/geometries/internal/layers.jl#L215).

**Accuracy.** Exact; checked against brute-force sampling (P3.2, package tests) and by the single-spin trace (max 1.0e-7, set by the trace's own sampling; P3.4.1).

### 2.3 Cylinders, linear and exponential: adaptive Gauss–Legendre quadrature

**What.** ∫₀¹ g(d(s)/h) ds over the part of the piece inside the support, where d(s) = |r(s) − R| and r(s) = √(a s² + b s + q) is not polynomial in s.

**Algorithm (decision D3.3).**
1. Split the piece at the exact shell roots (r = R and r = R ± L) from `disk_interval`, and at the closest approach to the axis, s* = −b/(2a). On each part r(s) is smooth and monotone.
2. Integrate each part with **8-point Gauss–Legendre**. The nodes and weights come from the Golub–Welsch method: the eigenvalues of the symmetric tridiagonal Jacobi matrix with off-diagonal β_k = k/√(4k² − 1) are the nodes, and the weights are w_k = 2·v_{1k}² (`gauss_legendre`, `gl8`).
3. **Adaptive bisection** (`quad_adaptive`). The two-halves estimate is accepted when it agrees with the whole-interval estimate to tol·(b − a); otherwise the interval is bisected, up to depth 30. The tolerance is `LAYER_QUAD_TOL = 1e-10`.
4. Prefactor (ρ/h)/κ·dt.

**Where.**
- Quadrature: `gauss_legendre`, `gl8`, `quad_adaptive` at [layers.jl:163-213](../src/geometries/internal/layers.jl#L163).
- Splitting and integration: `shell_profile_integral` at [:261-301](../src/geometries/internal/layers.jl#L261).
- The constants are type-stable and the path allocates nothing (package test "quadrature path does not allocate").

**Accuracy (P3.3, P3.6.3).**
- Linear against its closed form ∫√(quadratic): ≤ 9e-11.
- Exponential: ≤ 3.2e-11 at tol = 1e-10.
- Single-spin traces: 6.7e-12 (linear), 2.4e-8 (exponential).

**Cost.** ×1.2 (linear) and ×2.5 (exponential) against the step.

### 2.4 Normalisation: centroid and κ

**Profiles** (`profile_g`):
- step: g = 1 on [0, 1];
- linear: g = 2(1 − u) on [0, 1];
- exponential: g = e^(−u)/(1 − e^(−c)) on [0, c], with c = L/h, truncated and renormalised.

**Centroid** c̄ = h∫u g du, in closed form:
- step h/2;
- linear h/3;
- exponential h(1 − (c + 1)e^(−c))/(1 − e^(−c)).

**Curvature factor.** κ = ∫g(u)(1 ∓ u h/R) du = 1 ∓ c̄/R, with − inside and + outside. It makes ρ the relaxation per unit surface area (D3.1).

**Where.** `profile_g` [layers.jl:46](../src/geometries/internal/layers.jl#L46); `centroid`, `kappa` [:148-160](../src/geometries/internal/layers.jl#L148).

**Checks.**
- κ against quadrature: max deviation 1.6e-9 (D3.1 check).
- In the research code κ is also computed independently: by 200 000-point midpoint quadrature (`side_prefactor`, [phase3/harness.jl:115](phase3/harness.jl#L115)), and from exact ring integrals (§4.2).

---

## 3. Periodic geometry: images, overlap and the timestep bound (package)

**Method of images.**
- **Walls:** every image wall c + k·spacing whose support can reach the piece is summed ([layers.jl:329-352](../src/geometries/internal/layers.jl#L329)).
- **Cylinders:** `image_centres` lists the lattice images within reach R + L of the piece's bounding box ([layers.jl:360-383](../src/geometries/internal/layers.jl#L360)).
- Overlapping layers **add**: the rate is the sum over every surface within its support.

**Validation of supports** ([src/geometries/user/fix.jl](../src/geometries/user/fix.jl)).
- **Walls:** step and linear need h ≤ the smallest wall spacing; for the exponential the cutoff must be ≤ the spacing.
- **Cylinders:** step and linear inside need h ≤ R. Outside, the support must not exceed the surface gap to any neighbour (`surface_gaps`, `check_outer_reach`), so a layer never reaches inside another cylinder.
- **Exponential defaults:** the cutoff is the smallest gap outside and R inside (D3.5).

**Summed peak-rate bound `max_layer_rate`** ([layers.jl:397-460](../src/geometries/internal/layers.jl#L397)). An upper bound on ΔR₂ summed over overlapping layers. It is exact when nothing overlaps.
- **Walls:** the largest single surface rate, or across each gap the sum of the two facing surface rates when L₊ + L₋ > gap·(1 + 1e-12).
- **Cylinders:** for each cylinder, its surface rate plus those of every neighbour whose outer shell intersects its own. This is conservative: up to 5× the true peak in an s = 2.5 array.

**Timestep.** τ ≤ c/ΔR₂,max with c = 0.005 (`DEFAULT_LAYER_SCALING`, [src/timesteps.jl:13](../src/timesteps.jl#L13)); adjustable with `timestep=(layer=…)`. This keeps the error of a single model's decay rate at ≈ 0.15 % (P2.6).

---

## 4. Continuum reference solvers (research)

**Equation.** All three solve the deterministic diffusion–relaxation equation for a uniform start, with no Monte Carlo and no timestep:

  ∂m/∂t = D∇²m − [R₂,bulk + ΔR₂(x)]m,  S(t) = ⟨m⟩/⟨m₀⟩.

**Common scheme.**
- Finite volumes give a semi-discrete system A ṁ = −(K + AV)m (A cell volumes, K conductances, V cell-averaged ΔR₂).
- The substitution y = √A·m makes the operator symmetric: M = A^(−½)(K + AV)A^(−½).
- M is tridiagonal, and its eigendecomposition gives the signal as a spectral sum:

  S(t) = Σ_k (u_kᵀ√A)² e^(−λ_k t) / ΣA · e^(−R₂,bulk t).

### 4.1 Walls, 1D (Phase 2)

- **Domain.** x ∈ [0, w] between reflecting walls, with N = 3000 cells and a Neumann Laplacian.
- **Layer.** ΔR₂ comes from both faces, averaged over each cell by 16-point midpoint sub-sampling.
- **Conventional model.** A sink ρ/dx in the two wall cells (first order in dx).
- **Where.** [phase2/continuum.jl](phase2/continuum.jl) (`continuum_operator`, `continuum_signal`), Review C1.

### 4.2 Single cylinder, inside: radial (Phase 3)

- **Grid.** N = 2000 rings of width dr = R/N, with ring areas A_i = π(r_{i+½}² − r_{i−½}²) and face conductances 2πD r_{i+½}/dr.
- **Layer.** Averaged over each ring with **exact ring integrals** ∫g((R − r)/h) r dr, using antiderivatives for step, linear and exponential. κ is taken from the same exact integral, after a midpoint rule failed for thin exponential layers.
- **Conventional model.** A **face-centred Robin** condition on r = R, with conductance 2πRρ/(1 + ρ dr/(2D)). This is second order; an earlier ring-centre sink was first order, with an error ≈ 5e-7.
- **Accuracy.** Grid error |S(N = 2000) − S(N = 4000)| ≤ 5e-8 (P3.4.6–P3.4.8); ≤ 9e-7 at h/R = 0.01.
- **Where.** [phase3/continuum_radial.jl](phase3/continuum_radial.jl) (`radial_operator`, `radial_signal`, `radial_rates`, `ring_integral`, `radial_kappa`).

### 4.3 Outside a cylinder: annulus (Phase 3.7b, approximate)

- **Domain.** A Krogh-type cell R ≤ r ≤ R_out with equal area (π R_out² = s², or s_x·s_y), reflecting at R_out.
- **Scheme.** The same as §4.2, with outer ring integrals and the Robin face at r = R.
- **Truncated layers.** A layer cut at R_out is **renormalised** by 1/kept, so that ∫ΔR₂ dA = 2πRρ (user decision after validation).
- **Roundoff floor.** Eigenvalue error ≈ ε‖M‖ ≈ 1e-8 /ms at N = 2000, so S is exact only to ≈ 1e-8·t.
- **Validity.** Validated against paired Monte Carlo only at wide spacing and weak relaxation (s = 4, ρR/D ≤ 0.1). At s = 2.5, or at strong ρ, the square cell differs from the circle by as much as the effect being measured.
- **Where.** [phase3/continuum_annulus.jl](phase3/continuum_annulus.jl). Validation: [p3_7b_validate_analysis.jl](phase3/p3_7b_validate_analysis.jl).

---

## 5. Exact analytic references

| Reference | Formula | Method | Where |
|---|---|---|---|
| Conventional decay inside a cylinder | lowest D k² with D k J₁(kR) = ρ J₀(kR) | bisection on k ∈ (0, j₀₁/R) | [continuum_radial.jl:99](phase3/continuum_radial.jl#L99) |
| Step layer inside, exact rate | core r < R − h: J₀(kr), λ = Dk²; shell: I₀(βr) + (I₁(βR)/K₁(βR))K₀(βr), β = √((q − λ)/D), q = ρ/(h(1 − h/2R)); λ makes m′/m continuous at R − h | bisection on λ ∈ (0, q) | [continuum_radial.jl:116](phase3/continuum_radial.jl#L116) |
| Annulus, conventional and ceiling | u(r) = J₀(kr)Y₁(kR_out) − Y₀(kr)J₁(kR_out) satisfies u′(R_out) = 0. Robin: D u′(R) = ρu(R); absorbing (ceiling): u(R) = 0. Written without dividing by Y₁, which has zeros in the bracket. | scan for the first sign change, then bisection (`first_root`) | [continuum_annulus.jl](phase3/continuum_annulus.jl) |
| Decay ceilings (fastest conventional decay) | walls D·π²/w²; cylinder D·j₀₁²/R² (j₀₁ = 2.4048); annulus Dirichlet root | closed form or root | [analysis.jl:64](phase3/analysis.jl#L64) |
| Mean-exponent identity | E[−ln M⊥] = ρ(S/V)t for any profile and h; conventional θ√(D/π)(S/V)t | closed form | [phase3/harness.jl:60-62](phase3/harness.jl#L60) |
| Uniform-layer limit | step with h = R inside (or filling the annulus): S = e^(−ρt/(hκ)) | closed form | package tests; test_continuum_annulus.jl |
| Mapping θ ↔ ρ | θ = ρ√(π/D) | closed form; MC ratio 0.99976–1.00013 for √(Dτ)/R ≤ 0.35 | [p3_1_6_mapping.jl](phase3/p3_1_6_mapping.jl) |

---

## 6. Independent verification (brute force)

Each check is written without the package's layer code, except where noted.

| Check | Method | Accuracy target | Where |
|---|---|---|---|
| Summed rate at a point | direct sum over ±3 image cells, inclusive support, rectangular (s_x, s_y) | — | [phase3/harness.jl:128](phase3/harness.jl#L128) `brute_rate_profile` |
| Path exponent, midpoint rule | 20 000-point midpoint average of the brute-force rate | ≤ 1e-4 (step), 2e-5 (profiles) | [test/test_layer_cylinders.jl](../test/test_layer_cylinders.jl) |
| Path exponent, exact | split at every line–circle crossing of r = R and r = R + L of every image, then 64-point Gauss–Legendre per smooth piece; stationary points exponent = rate·dt | ≤ 1e-10 (segments), ≤ 1e-12 (rates) | test_layer_cylinders.jl, P3.5.3 testset (uses `Layers.kappa`, `profile_g` and `gauss_legendre`; κ is pinned independently elsewhere) |
| Single-spin trace | the MC exponent of every step against the exact or sampled exponent | ≤ 1e-6 | [p3_4_1_trace.jl](phase3/p3_4_1_trace.jl), phase2/p2_2_9_layer_trace*.jl |
| κ by quadrature | 200 000-point midpoint quadrature of ∫g(1 ∓ uh/R)du | 1.6e-9 against closed form | [harness.jl:115](phase3/harness.jl#L115) |
| Continuum grid convergence | S at N against 2N | ≤ 1e-6 | continuum tests, map analyses |
| Realism of summed rates (2D) | 801² grid plus a **ray test**: from every surface point along its outward normal to the own Voronoi cell, the summed rate never exceeds its surface value. t is sampled log-spaced from 1e-9 µm and angles are refined to ±1e-6 rad near facing directions; the sampled excess is checked to reach the analytic e^(−gap/h). | excess ≤ 1e-9 | [p3_5_realism.jl](phase3/p3_5_realism.jl) |
| Coverage and overlap | 600² midpoint grid of the unit cell: fraction covered by ≥ 1 and ≥ 2 supports, farthest distance from every surface | — | [p3_5_configs.jl](phase3/p3_5_configs.jl) |
| Roundoff bound | (T/τ + N)·ε for null and bulk runs | per spin 5e-14 | [harness.jl:57](phase3/harness.jl#L57) |

---

## 7. Statistics

| Quantity | Method | Where |
|---|---|---|
| Mean and SEM | over independent seeds (typically 10 × 2·10⁴ spins); SEM = std/√n | `run_cylinders` [harness.jl:74](phase3/harness.jl#L74) |
| Paired differences | every model of a condition runs with the same seeds and τ, so the trajectories are identical (checked on final positions); ΔS = S_layer − S_conv is formed per seed, which cancels the common noise | [p3_7_mc.jl](phase3/p3_7_mc.jl), [p3_7b_mc.jl](phase3/p3_7b_mc.jl) |
| Independent conditions | seeds 1000·k + (1:10) per condition (3.7b on); in 3.7a, shared seeds made conditions correlated, which a fresh-seed check exposed | [p3_7b_run.jl](phase3/p3_7b_run.jl), [p3_7_seedcheck.jl](phase3/p3_7_seedcheck.jl) |
| Summed compartments | S = f·S_in + (1 − f)·S_out, f = πR²/(s_x s_y); SEMs in quadrature (independent runs); identical readout times asserted | `sum_compartments` in [p3_7b_mc.jl](phase3/p3_7b_mc.jl) |
| z-scores | z = (MC − reference)/SEM; ratio ΔS_mc/ΔS_ref with σ = SEM_ΔS/\|ΔS_ref\| | analysis scripts |
| Pooling | inverse-variance and equal-weight means (stated with the method) | [p3_7_seedcheck.jl](phase3/p3_7_seedcheck.jl), record P3.7a |
| Error propagation (first order) | σ(T₂) = TE·SEM_S/(S ln²S); σ(ΔT₂/T₂) from the paired SEM_ΔS; σ_Δq = (1 + q_c)√(σ_δ₀.₅² + σ_δ₃²) | [p3_7_shared.jl](phase3/p3_7_shared.jl), [p3_7b_ddep_analysis.jl](phase3/p3_7b_ddep_analysis.jl) |
| Two-point T₂ noise | T₂ = −TE/ln S(TE); σ_T₂/T₂ = √(1 + r²)/(SNR·r·\|ln r\|), r = S(TE). At TE* (r = e⁻¹) the 2σ level is 2.9 % at SNR 200. Two-measurement band 2√(n₁² + n₂²). | [analysis.jl:55-62](phase3/analysis.jl#L55) |
| Jackknife | leave-one-seed-out standard error of the fitted h* | [phase2/fit.jl](phase2/fit.jl) `jackknife_h` |
| Monotonicity steps | a step is resolved when it exceeds 2√(σ₁² + σ₂²), or unpaired σ (stated as conservative) | [phase2/fit.jl:17](phase2/fit.jl#L17), [p3_5_overlap_analysis.jl](phase3/p3_5_overlap_analysis.jl) |

---

## 8. Root finding, fitting and threshold location

| Method | What | Where |
|---|---|---|
| Bisection (200 steps) | TE* where S(TE*) = e⁻¹ | `te_star` [analysis.jl:67](phase3/analysis.jl#L67) |
| Log-scale bisection | the h giving a target centroid (exponential); errors if the target is unreachable (c̄ ≥ L/2) | `h_for_cbar` [p3_7_shared.jl](phase3/p3_7_shared.jl), `h_for_cbar_out` [p3_7b_mc.jl](phase3/p3_7b_mc.jl) |
| Minimax fit plus golden section | best conventional ρ_c = argmin max_t \|S_layer − S_conv(ρ_c)\|: a 61-point log grid over e^(−4)…e^8·ρ, then 60 golden-section iterations. Flags: at the grid edge or on a plateau reaching it, or layer rate above the ceiling. Flagged fits report no ρ_c or residual; the edge residual is an upper bound. | `best_fit_conventional` [analysis.jl:87](phase3/analysis.jl#L87) |
| Weighted least squares in log–log | power law y = a xᵖ with the slope's standard error from the residuals | `fit_power` [analysis.jl:10](phase3/analysis.jl#L10) |
| Threshold crossing | first x where y > level, interpolated log–log between bracketing points; status crossed / below_range (value is a lower bound) / above_range. One convention everywhere. | `threshold_status` [analysis.jl:38](phase3/analysis.jl#L38) |
| χ² scan plus refinements (Phase 2) | h*(θ): χ² over a 10 001-point h grid with the model interpolated per readout; parabola vertex; Gauss–Newton step | [phase2/fit.jl](phase2/fit.jl) (`fit_h`, `parabola_vertex`, `gauss_newton_step`) |
| Timestep convergence (Phase 2) | τ_conv: the largest τ with relative rate error ≤ tol, interpolated linearly in log τ; τ_conv = c_h·h²/D fitted as a mean ratio | [phase2/timestep.jl](phase2/timestep.jl) |
| Quadratic fit of ln S | curvature c of ln S = a + bt + ct²; 0 for a single exponential | `curvature` [phase2/fit.jl:11](phase2/fit.jl#L11) |
| Interpolation-bias check | threshold from 3 decade points against 49 log points on a smooth test curve: ≤ 7.2 % | [p3_7b_interp_bias.jl](phase3/p3_7b_interp_bias.jl) |

---

## 9. Rules derived from the numerics

- **Default layer timestep.** ΔR₂,max·τ ≤ 0.005 keeps the error of a single model's decay rate at ≈ 0.15 % (P2.6). It is the package default.
- **Comparison timestep rule** (P3.6, user decision option (b)):
  - for differences between models at the same ρ, the error is controlled by √(2Dτ)/c̄, not by ΔR₂τ;
  - the rule is τ = min(0.135·c̄_min², 0.005/ΔR₂,max) at D = 3 (√(2Dτ) ≤ 0.9 c̄), for ≲ 10 % error in ΔS. It was measured at c̄ = 0.1 and confirmed at 0.05; below that it is extrapolated.
  - Where: [p3_6_1_tau_scan.jl](phase3/p3_6_1_tau_scan.jl), [p3_6_2_analysis.jl](phase3/p3_6_2_analysis.jl).
- **Annulus use rule** (P3.7b, pre-registered):
  - a model passes if |S_mc − S_ann| ≤ 1 %·S_ann + 2·SEM at every readout;
  - the ΔS ratio must lie in [0.8, 1.2] within 2σ;
  - it is applied per spacing, and at s = 4 per map condition.

---

*Compiled 2026-11-11 from the code and the Phase 2–3 records. Line numbers refer to commit f447f5f or later; if code moves, search for the function names given.*
