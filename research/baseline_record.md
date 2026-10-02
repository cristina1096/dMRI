# Phase 1 · Baseline Record (P1.3, P1.4)

*Naming:* Rep1–Rep5 are the Phase 1 reproduction targets P1.3.1–P1.3.5 (formerly labelled R1–R5; renamed 2026-10-01 to avoid confusion with the relaxation rates R₁, R₂). Script and result folder names (`r1_…`, `r2_…`) keep the old prefix. "Walls (w = …)" is the parallel-plane geometry `Walls(repeats=w)`, previously called "slab".

Baseline results of the **unmodified** MCMRSimulator. Each run lists every setting needed to regenerate it, so that any later signal change can be attributed to the R₂(d) modification rather than to the setup.

- Code map and trace: [phase1_code_log.md](phase1_code_log.md)
- Scripts: [research/baseline/](baseline/)
- Raw outputs: [research/results/baseline/](results/baseline/). Each run folder holds CSV data, `summary.json` with provenance, `run.log` and a figure.

## Summary

| ID | Result | Reference | Criterion | Outcome |
|---|---|---|---|---|
| P1.3.1 | **Rep1** surface relaxation, walls (w = 2 µm) | exact solution (Brownstein–Tarr lowest mode); Eq. 17 recorded | within 2 SEM of the exact rate | **PASS**: +1.7 SEM vs exact (relative +2.4 × 10⁻⁴). Eq. 17 (first-order approximation) is 0.11% above the exact rate: −5.9 SEM, resolved by the data |
| P1.3.1 | **Rep1** surface relaxation, single cylinder | exact solution (Brownstein–Tarr lowest mode); Eq. 17 recorded | within 2 SEM of the exact rate | **PASS**: −0.5 SEM vs exact (relative −4.6 × 10⁻⁵). Eq. 17 is 0.08% above the exact rate: −10 SEM, resolved by the data |
| P1.3.2 | **Rep2** compartment T₂ (spin echo, MT = 5e-3) | Fig. 8 (notebook `se_grid_plot`) | 5% (read off figure) | **PASS**: every TE compared within 5% of the digitised notebook figure (max 4.9%, mean \|Δ\| ≤ 2.4%), with and without myelin susceptibility |
| P1.3.3 | **Rep3** timestep plateau, ADC vs τ | Supp. Fig. S1 | plateau onset near τ ≈ 0.01 ms | **PASS**: onset (2σ) 6.3 × 10⁻³ – 2.5 × 10⁻² ms in all 4 packings |
| P1.3.4 | **Rep4** correction schemes, walls (w = 1 µm) | Supp. Fig. S2C | log scheme flat to τ ≈ 10⁻³ ms; record θ_relax | **PASS**: flat within 2σ to 3.3 × 10⁻³ ms (toy) and 1.25 × 10⁻² ms (MCMR) at θ = 10; MCMR matches the toy log scheme at every τ |
| P1.3.5 | Rep5 (optional) | Fig. 5 | — | not run |
| B0 | Existing test suite (collisions, evolve, permeability, transfer, known_sequences) | package tests | all pass | **PASS**: 442 / 442 |
| B1 | Wall reference config, null run (no relaxation) | exact: M⊥ = 1 | \|M⊥/M⊥0 − 1\| ≤ 1e-12 | **PASS**: deviation exactly 0 at every readout, both τ |
| B2 | Wall reference config, bulk R2 only | exact: exp(−R2·t) | roundoff bound (n_steps + N)·ε; positions == R2 = 0 | **PASS** for R2 = 1/80 – 0.1 ms⁻¹: 2.6e-12 ensemble, 3e-14 per spin; trajectories unchanged by R2 |
| B3 | Wall reference config, density profile and layer occupancy | uniform, fraction = λ/w | \|z\| ≤ 2 (occupancy), \|z_χ²\| ≤ 3 (profile) | **PASS** at τ = 0.04, 1e-2, 1e-3, 1e-4 ms: no density artefact near the walls, largest occupancy deviation 1.4% (z = 2.2), no systematic sign |
| B4 | Baseline θ_relax reference curves S_θ(t) on walls (both faces), θ = 0, 0.001–0.05 | reference (not a formula) | τ-converged: τ = 1e-3 vs 1e-2 within 2 SEM; bulk factorisation ≤ roundoff bound | **Reference set** at τ = 1e-2 ms; bulk **PASS** (2.5e-12); default τ = 0.04 ms shows a small systematic bias (≤ +0.055% in S, up to 2.7 SEM at θ = 0.05), so it is not used as reference |
| B5 (slim) | Runtime vs τ, walls, θ = 0 and 0.01 | — (cost reference) | — | **Done**: per-step cost flat at ≈ 0.016 µs (8 threads); 1e5 spins × 50 ms costs 81 s / 13 min / 2.1 h per seed at τ = 1e-3 / 1e-4 / 1e-5 ms |

## Baseline test index (what, why, result)

All tests run on unmodified `src/`. P1.x tests reproduce the paper (Phase 1). B-tests are the wall-geometry baselines for Phase 2 (reference configuration: `Walls(repeats=2)`, D = 3, no RF, free decay, readouts 0–50 ms).

| ID | What is tested | Why (what it protects later) | Result | Section |
|---|---|---|---|---|
| P1.2.6 | One cylinder, 100 spins, spin echo; one spin's position and M logged every step | Confirms the reading of the code (step draw, reflection, `θ_relax` application) before editing it | PASS (max rel. error 6e-16) | [trace](#run-p126-single-spin-trace) |
| Rep1 walls | Surface relaxation `θ_relax` on walls, w = 2 µm, vs Eq. 17 | Paper reproduction. Phase 2 uses B4 (reference curves in the wall reference config) instead | PASS vs exact (+1.7 SEM); Eq. 17 off by 0.11% (−5.9 SEM) | [Rep1](#run-rep1-surface-relaxation-p131) |
| Rep1 cylinder | Same, one cylinder a = 1 µm | Curved-geometry reference for Phase 3 | PASS vs exact (−0.5 SEM); Eq. 17 off by 0.08% (−10 SEM). Decided 2026-10-02 | [Rep1](#run-rep1-surface-relaxation-p131) |
| Rep2 | Compartment T₂ in myelinated white matter | Confirms the full simulator reproduces a published result with realistic settings | PASS (≤ 4.9% vs digitised Fig. 8) | [Rep2](#run-rep2-compartment-t-in-myelinated-white-matter-p132) |
| Rep3 | Extra-axonal ADC vs τ, 4 cylinder packings | Where diffusion/collision results stop depending on τ; the default τ sits at the knee for tight packings | PASS (onset 6e-3–2.5e-2 ms) | [Rep3](#run-rep3-timestep-plateau-of-extra-axonal-adc-p133) |
| Rep4 | Linear / logarithmic / log+return correction schemes, walls (w = 1 µm) | Identifies the `θ_relax` scheme used (logarithmic), i.e. what the B4 reference curves are made of | PASS (log scheme flat to 1.25e-2 ms at θ = 10) | [Rep4](#run-rep4-surface-relaxation-correction-schemes-walls-w--1-µm-p134) |
| Rep5 | Diffusion vs Mitra / van Gelderen (optional) | — | not run | — |
| **B0** | Existing package tests for collisions, evolve, permeability, transfer, known_sequences | A test that fails after the R₂(d) change must be one you broke, not one already failing | **PASS** 442 / 442 | [B0](#run-b0-existing-test-suite) |
| **B1** | Wall reference config with nothing relaxing: M⊥(t)/M⊥(0) and phase; τ_max chosen; runtime | Target for V0b (λ = 0 must reproduce this) and for P2.2.8 (λ = 0 must not slow down) | **PASS**: exactly 1 and 0°; τ_max = 0.04 ms (tortuosity); 0.013 µs per spin-step | [B1](#run-b1-wall-reference-configuration-null-run) |
| **B2** | Bulk R2 alone (1/80, 0.025, 0.05, 0.1 ms⁻¹) in the reference config: signal, phase, trajectories, runtime | Fixes how the old code applies R2_bulk, the first term of R2_total = R2_bulk + ΔR2(d). Phase 2 must keep it intact (P2.2.5, P2.2.7, P2.2.8, V0b, additivity test) | **PASS**: S = exp(−R2·t) to 2.6e-12 (per spin 3e-14); positions bit-identical to R2 = 0; runtime +0–6% | [B2](#run-b2-bulk-r2-in-the-wall-reference-configuration) |
| **B3** | Density across the gap and fraction of spins within λ of each face, relaxation off, 4 timesteps | The layer relaxes a spin only while it is within h of a face, so the time spins spend there depends on the density near the wall. If the old code distorts that density, every Phase 2 result inherits it. Baseline for V9 / P2.3.3 (the Phase 2 gate) | **PASS** at all τ; no near-wall artefact even at τ = 0.04 ms (step 0.49 µm > λ) | [B3](#run-b3-density-and-layer-occupancy-between-walls) |
| **B4** | Baseline reference curves S_θ(t), θ = 0 and 0.001–0.05, walls, both faces; timestep convergence; bulk factorisation | The curves the R₂(d) model is fitted to (h\*(θ) at fixed ΔR₂(0)); θ = 0 is the no-surface-relaxation reference. Also a non-zero noise floor for P2.1.4 | **Reference set** (τ = 1e-2 ms, converged: τ = 1e-3 within 2 SEM). S(50) spans 0.95 → 0.088. Bulk factorises to 2.5e-12. Default τ = 0.04 ms slightly under-attenuates (≤ 0.055%) | [B4](#run-b4-baseline-reference-curves-s_θt) |
| **B5** (slim) | Runtime of the old code vs τ = 1e-1 … 1e-5 ms, θ = 0 and 0.01 | Cost reference for the timestep study (P2.6) and feasibility of small τ required by τ ≤ c_h·h²/D (proposal Eq. 11). Signal vs τ is covered by B4 | **Done**: 0.0155–0.017 µs per spin-step for τ ≤ 1e-2 (flat); surface relaxation adds ≤ 2%. One Phase 2 seed (1e5 spins, 50 ms): 81 s at 1e-3, 13 min at 1e-4, 2.1 h at 1e-5 | [B5](#run-b5-slim-runtime-vs-timestep) |
| B6 | Route A: layer built from permeable boundaries + inside R2 | Independent check on the step profile (V10). Walls have no inside volume, so this probably needs meshes | **deferred to future work** (decision 2026-09-30) | — |

**Two findings about the paper's figure notebooks** (they affect how references are read):
1. **Supp. Fig. S1 x-axis mislabelled.** `turtuosity.ipynb` computes ADC at τ = 10^(−3:0.2:1) (cell 17) but plots it against τ = 10^(−4:0.25:1) (cell 18). The saved figure is therefore shifted: its first point is τ = 10⁻³ ms, not 10⁻⁴. Rep3 is compared at the true τ.
2. **Supp. Fig. S2C x-axis.** In the saved image, all three rows use the x-axis and dashed line of the T = 0.1 ms row (dashed line at τ = 10⁻³ ms = 0.01·T). The tracker's "flat to τ ≈ 10⁻³ ms" therefore refers to T = 0.1 ms, θ_relax = 10.

## Common environment (all runs)

| Item | Value |
|---|---|
| Julia | 1.12.7 (juliaup, aarch64-apple-darwin) |
| MCMRSimulator | 1.1.0, this repository, path source `../..`; `src/` identical to commit `d9c5d0e` (`provenance.src_modified = false`) |
| Git | branch `cc/near-surface-r2`, HEAD `c73c723` at run time (only `research/` differs from `d9c5d0e`) |
| MRIBuilder | 0.4.1 (`https://git.fmrib.ox.ac.uk/ndcn0236/mribuilder.jl.git#v0.4.1`) |
| Other packages | pinned in [research/baseline/Manifest.toml](baseline/Manifest.toml) (CairoMakie 0.15.15, JSON 1.10, SpecialFunctions 2.9, Bessels 0.2.8, …) |
| Paper reference environment | notebooks `mcmr_paper_figures` @ `c7b86f9`: Julia 1.11.7, MCMRSimulator v1.0.0, MRIBuilder v0.3.0 |
| Hardware | Apple M1 Pro, 10 cores (Julia reports 8 `CPU_THREADS`), 16 GB RAM, macOS (Darwin 25.5.0 arm64) |
| Threads | `julia -t 8` (P1.2.6 trace: `-t 1`) |
| Run date | 2026-09-29 |
| Reproducibility | `Random.seed!(seed)` before each `readout` / `Snapshot`. Spin positions and per-spin RNG streams come from the global RNG, so results do not depend on thread count |

How to regenerate everything (≈ 25 min on the machine above):
```sh
julia --project=research/baseline -e 'using Pkg; Pkg.instantiate()'
julia --project=research/baseline -t 1 research/baseline/p1_2_6_single_spin_trace.jl
julia --project=research/baseline -t 8 research/baseline/r1_surface_relaxation.jl
MYELIN=false julia --project=research/baseline -t 8 research/baseline/r2_compartment_t2.jl
MYELIN=true  julia --project=research/baseline -t 8 research/baseline/r2_compartment_t2.jl
julia --project=research/baseline -t 8 research/baseline/r3_timestep_plateau.jl
julia --project=research/baseline -t 8 research/baseline/r3b_geometry_realisation.jl
julia --project=research/baseline -t 8 research/baseline/r4_correction_schemes.jl
# Wall baselines for Phase 2 (≈ 30 min):
julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["collisions","evolve","permeability","transfer","known_sequences"])'   # B0
julia --project=research/baseline -t 8 research/baseline/b1_walls_null.jl
julia --project=research/baseline -t 8 research/baseline/b2_walls_bulk_r2.jl
julia --project=research/baseline -t 8 research/baseline/b3_walls_density.jl
julia --project=research/baseline -t 8 research/baseline/b4_theta_reference.jl   # ≈ 45 min
julia --project=research/baseline -t 8 research/baseline/b5_runtime_vs_tau.jl    # ≈ 3 min
# Rep2 figure comparison (needs a clone of mcmr_paper_figures):
python3 research/baseline/r2_digitise_notebook_fig8.py <mcmr_paper_figures>/Figure_8_9/gradient_spin_echo.ipynb
```
`NSPINS`, `NSEEDS` (and `NSPINS_TOY`, `NSPINS_MCMR`, `NGEOM`) can be overridden by environment variables. The values below are the defaults.

---

## Run Rep1: surface relaxation (P1.3.1)

Script [r1_surface_relaxation.jl](baseline/r1_surface_relaxation.jl) → [results/baseline/r1_surface_relaxation/](results/baseline/r1_surface_relaxation/)

| Field | Walls (w = 2 µm) | Single cylinder |
|---|---|---|
| Geometry | `Walls(repeats=2)`: walls every w = 2 µm, both faces active | `Cylinders(radius=1)`: one cylinder, not repeating |
| Exact S/V | 2/w = **1.0 µm⁻¹** | 2/a = **2.0 µm⁻¹** (inside) |
| θ_relax (`surface_relaxation`) | 0.01 µm⁻¹ ms^−1/2 | 0.01 |
| ρ = θ√(D/π) | 0.0097721 µm/ms | 0.0097721 |
| Physics | D = 3 µm²/ms; R1 = R2 = 0; off-resonance 0; permeability 0; no sticking/MT | same |
| Sequence | `GradientEcho(TE=200, TR=1000, Siemens_Prisma, excitation=(phase=90,))`, instant 90°, no gradients before TE | same |
| Readouts | t = 10, 20, …, 100 ms | same |
| τ_max / binding constraint | **0.04 ms** / tortuosity 0.03·(2 µm)²/3 (surface-relaxation bound 0.01/θ² = 100 ms) | **0.01 ms** / tortuosity 0.03·(1 µm)²/3 |
| N_spins | 10 000 per seed (all count) | 10 000 initial in the 2×2×2 µm box; inside subset ≈ 7 860 per seed |
| Seeds | 1–10 | 1–10 |
| Signal S(100 ms)/S(0) | 0.3767 | 0.1419 |
| **R₂ fit** (slope of −log S vs t, with intercept) | **0.0097638 ± 0.0000044 ms⁻¹** (σ, 10 seeds; SEM 1.4e-6) | **0.0195273 ± 0.0000053 ms⁻¹** (SEM 1.7e-6) |
| Apparent T₂ | 102.4 ms | 51.2 ms |
| Eq. 17 R₂ = θ√(D/π)·S/V | 0.0097721 | 0.0195441 |
| Exact Brownstein–Tarr lowest mode | 0.0097614 (ξ tan ξ = ρa/D, a = w/2) | 0.0195282 (ξJ₁/J₀ = ρa/D) |
| Deviation vs Eq. 17 | −1.87σ → **PASS** | −3.15σ → **2σ not met** |
| Deviation vs Brownstein–Tarr | +0.53σ | −0.17σ → PASS |
| Eq. 17 / BT − 1 | 1.09 × 10⁻³ | 8.1 × 10⁻⁴ |
| Runtime (10 seeds) | 14 s | 61 s |

**Criterion decision (2026-10-02).** The original tracker criterion was "within 2σ of Eq. 17", with σ the spread of one seed's fitted rate. Two problems: the tested value is the mean of 10 seeds, whose uncertainty is the SEM = σ/√10 (so 2σ ≈ 6 SEM, lenient); and Eq. 17 is a first-order approximation, measurably above the exact rate at this precision. Rep1 is therefore judged **against the exact Brownstein–Tarr rate within 2 SEM**: walls +1.7 SEM, cylinder −0.5 SEM → **PASS** for both. Against Eq. 17 both geometries are resolved as different (walls −5.9 SEM, cylinder −10 SEM), which is the approximation error of Eq. 17, not a simulator error. The σ-based numbers in the table above are kept for the record.

**Interpretation.** The simulator reproduces surface relaxation to 0.03% of the exact solution. The cylinder "failure" against Eq. 17 is a real, resolvable difference between Eq. 17 (the first-order fast-diffusion approximation, κ = ρa/D = 3.3 × 10⁻³) and the exact eigenvalue. It is not a simulator error. (Phase 2 does not compare against Eq. 17 or Brownstein–Tarr: the R₂(d) model is fitted to the baseline reference curves of B4 instead.)

---

## Run Rep2: compartment T₂ in myelinated white matter (P1.3.2)

Script [r2_compartment_t2.jl](baseline/r2_compartment_t2.jl) → [results/baseline/r2_compartment_t2_myelin_false/](results/baseline/r2_compartment_t2_myelin_false/) and [..._myelin_true/](results/baseline/r2_compartment_t2_myelin_true/). Figure comparison: [r2_digitise_notebook_fig8.py](baseline/r2_digitise_notebook_fig8.py) → `notebook_fig8_comparison.csv` in each folder.

| Field | Value |
|---|---|
| Geometry | `Random.seed!(123)`; `random_positions_radii((60,60), 0.7, 2; mean=1, variance=0.2, min_radius=0.3)` → **651 annuli** (same count as the notebook), outer-radius area fraction 0.7004; `Annuli(inner=0.7r, outer=r, rotation=:y, repeats=(60,60), myelin=MYELIN)` |
| Exact S/V (per unit length) | intra 2.358 µm⁻¹ · myelin 5.502 µm⁻¹ (inner + outer) · extra 3.859 µm⁻¹ |
| Volume fractions | intra 0.343 · myelin 0.357 · extra 0.300 |
| θ_relax | MT = 5 × 10⁻³ on inner and outer surfaces |
| Physics | D = 3; R1 = 1/4000 kHz; R2 (global, CSF) = 1/2000 kHz; susceptibility χ_I = χ_A = −0.1 ppm (defaults) when `myelin=true`; B0 = 3 T |
| Sequences | 40 × `SpinEcho(TE=5:5:200, TR=1000, Siemens_Prisma, excitation=(phase=90,))`: instant 90°/180°, one readout at TE each |
| Subsets | `inside=0` extra, `inside=1` myelin, `inside=2` intra |
| τ_max / binding | **0.01 ms** / tortuosity with user `size_scale=1` µm (as notebook) |
| N_spins / seeds | 1 000 per seed (as notebook) ≈ 305 extra / 347 myelin / 348 intra; seeds 1–10 |
| Runtime per seed | myelin=false ≈ 6.3 s; myelin=true ≈ 38 s |

Apparent T₂ = −TE / log(S(TE)/N) (notebook definition), mean ± σ over 10 seeds:

| Compartment | T₂ at TE = 5 ms | TE = 50 | TE = 100 | TE = 200 | Mono-exp fit (all TE) | Tracker "≈" | vs digitised notebook figure (mean / max \|Δ\|) |
|---|---|---|---|---|---|---|---|
| intra (myelin=false) | 83.3 ± 2.2 | 86.7 ± 2.0 | 90.2 ± 1.9 | 96.3 ± 2.1 | 97.1 ± 2.1 | 80 | −1.9% / 2.9% |
| extra (myelin=false) | 52.7 ± 1.9 | 54.3 ± 1.0 | 55.4 ± 0.7 | 56.7 ± 0.4 | 57.0 ± 0.3 | 60 | −0.4% / 2.8% |
| myelin (myelin=false) | 36.8 ± 0.6 | 40.3 ± 0.7 | 43.3 ± 0.7 | 48.4 ± 0.9 | 49.3 ± 0.9 | 45 | +2.1% / 4.6% |
| intra (myelin=true) | 83.3 ± 2.2 | 86.7 ± 2.0 | 90.3 ± 2.0 | 96.4 ± 2.2 | 97.2 ± 2.2 | 80 | −2.4% / 3.3% |
| extra (myelin=true) | 53.5 ± 1.7 | 54.2 ± 0.7 | 54.9 ± 0.9 | 56.2 ± 0.9 | 56.4 ± 0.9 | 60 | +1.4% / 4.9% |
| myelin (myelin=true) | 36.7 ± 0.5 | 39.9 ± 0.6 | 42.9 ± 0.7 | 47.7 ± 0.9 | 48.5 ± 0.9 | 45 | +1.4% / 2.7% |

**Interpretation.**
- The apparent T₂ rises with TE in every compartment, because heterogeneous axon sizes make each compartment's decay multi-exponential. The tracker's single values (80 / 60 / 45) are therefore only approximate readings.
- The correct reference is the notebook's own curve. Digitised pixel by pixel (±1 px ≈ ±2.6% in T₂), every compared TE agrees within 5%. → **PASS**.
- Myelin susceptibility changes the spin-echo T₂ by < 1 ms, as the notebook states.

---

## Run Rep3: timestep plateau of extra-axonal ADC (P1.3.3)

Scripts [r3_timestep_plateau.jl](baseline/r3_timestep_plateau.jl) and [r3b_geometry_realisation.jl](baseline/r3b_geometry_realisation.jl) → [results/baseline/r3_timestep_plateau/](results/baseline/r3_timestep_plateau/)

| Field | Value |
|---|---|
| Geometries | ordered_high `Cylinders(radius=1, repeats=[2.1,2.1])` (density 0.712); ordered_low `repeats=[3,3]` (0.349); random_high / random_low: positions from `random_positions_radii([50,50], φ, 2; mean=1, variance=0.01)` with `Random.seed!(1)`, cylinders `radius=1` (as notebook), `repeats=[50,50]` |
| Exact S/V (extra-axonal) | high: 2π/(2.1² − π) = **4.95 µm⁻¹**; low: 2π/(9 − π) = **1.07 µm⁻¹** (random: same by construction) |
| Physics | D = 3; no relaxation, no sequence (`Simulation([])`) |
| Observable | ADC = var(Δx ∪ Δy) / (2 · 10 ms), extra-axonal spins only |
| τ | fixed, 10^(−3:0.2:1) ms (21 values) |
| Simulator default τ_max | **0.01 ms** (tortuosity, size_scale 1 µm) for all four |
| N_spins / seeds | 10 000 initial (extra ≈ 2 900 high / 6 500 low); seeds 1–10 per τ |
| Runtime | 63–81 s per configuration (21 τ × 10 seeds) |

| Configuration | Plateau ADC (τ ≤ 1.6e-3) | Onset (2σ) | Onset (1%) | ADC at default τ = 0.01 | σ range |
|---|---|---|---|---|---|
| ordered_high | 1.249 ± 0.013 | 6.3e-3 ms | 2.5e-3 ms | 1.160 (−7.1%) | 0.004–0.029 |
| ordered_low | 2.220 ± 0.013 | 2.5e-2 ms | 4.0e-3 ms | 2.199 (−1.0%) | 0.007–0.037 |
| random_high | 0.788 ± 0.009 | 6.3e-3 ms | 4.0e-3 ms | 0.756 (−4.2%) | 0.004–0.018 |
| random_low | 2.092 ± 0.013 | 1.6e-2 ms | 4.0e-3 ms | 2.052 (−1.9%) | 0.012–0.032 |

**Comparison with the notebook** (true τ; see finding 1 above): the notebook's printed random_low values (one run, unseeded geometry) are 2–3% below ours at small τ. Across 10 geometry realisations the plateau is **2.039 ± 0.030** (`geometry_realisation.json`), and the notebook value 2.021 lies at −0.6σ. The offset is therefore geometry-realisation variance, not a simulator difference. Visually, the ordered_high (≈1.26) and ordered_low (≈2.2) plateaus match the notebook figure.

**Interpretation.** The plateau onset is at 0.006–0.025 ms → **PASS** (near 0.01 ms). At the simulator's default τ = 0.01 ms, the dense packings are already 4–7% below their plateau. Record this for Phase 2: the default tortuosity constant sits at the knee for tight packings.

---

## Run Rep4: surface-relaxation correction schemes, walls (w = 1 µm) (P1.3.4)

Script [r4_correction_schemes.jl](baseline/r4_correction_schemes.jl) → [results/baseline/r4_correction_schemes/](results/baseline/r4_correction_schemes/)

| Field | Value |
|---|---|
| Geometry | walls every 1 µm (toy: compartments [k, k+1); MCMR: `Walls(repeats=1)`) → **S/V = 2 µm⁻¹** |
| Physics | D = 0.5 µm²/ms (toy step `randn·√τ`); no bulk relaxation |
| **θ_relax used** | θ = 1/T: **0.1** (T = 10 ms), **1** (T = 1 ms), **10** (T = 0.1 ms) |
| Schemes | linear 1 − θ√τ; logarithmic e^(−θ√τ) (= MCMRSimulator, evolve.jl:557); log + return e^(−θ√τ)·I₀(θ√τ) |
| τ | T / N, N ∈ {1, 2, 3, 5, 8, 10, 20, 30, 50, 80, 100, 300, 1000, 3000} |
| Observable | attenuation 1 − S(T)/S(0) |
| Sequence (MCMR) | `GradientEcho(TE=T, TR=10T, excitation=(phase=90,))`, readout at T |
| MCMR default τ_max | T = 10: 0.06 ms (tortuosity, size_scale 1); T = 1: 0.01 ms (surface relaxation 0.01/θ²); T = 0.1: 1e-4 ms (surface relaxation) |
| N_spins / seeds | toy 100 000; MCMR 10 000; seeds 1–10 (toy seed = 1000·seed + j) |
| σ (MCMR) | 0.0005–0.0017 (T = 10), 0.0010–0.0027 (T = 1), 0.0021–0.0053 (T = 0.1) |

Attenuation (mean over 10 seeds) at selected N:

| T (θ) | Scheme | N = 1 | N = 8 | N = 50 | N = 100 | N = 3000 |
|---|---|---|---|---|---|---|
| 10 ms (0.1) | linear | 0.5224 | 0.5560 | 0.5520 | 0.5505 | 0.5462 |
| | log | 0.4729 | 0.5358 | 0.5440 | 0.5448 | 0.5452 |
| | log+return | 0.4510 | 0.5261 | 0.5400 | 0.5420 | 0.5447 |
| | **MCMR** | 0.4738 | 0.5353 | 0.5436 | 0.5446 | 0.5458 |
| 1 ms (1) | linear | 0.6318 | 0.5626 | 0.5296 | 0.5232 | 0.5094 |
| | log | 0.4355 | 0.4998 | 0.5063 | 0.5070 | 0.5065 |
| | log+return | 0.3766 | 0.4727 | 0.4953 | 0.4993 | 0.5051 |
| | **MCMR** | 0.4376 | 0.4998 | 0.5063 | 0.5073 | 0.5077 |
| 0.1 ms (10) | linear | 0.7956 | 0.4101 | 0.3564 | 0.3473 | 0.3286 |
| | log | 0.2412 | 0.3127 | 0.3241 | 0.3254 | 0.3250 |
| | log+return | 0.1925 | 0.2787 | 0.3097 | 0.3152 | 0.3231 |
| | **MCMR** | 0.2419 | 0.3135 | 0.3231 | 0.3259 | 0.3233 |

Flatness of the log scheme (largest τ with all smaller τ within tolerance of the smallest-τ value):

| T (θ) | toy log, 2σ | toy log, 5% | MCMR, 2σ | MCMR, 5% | toy linear, 5% | toy log+return, 5% |
|---|---|---|---|---|---|---|
| 10 ms (0.1) | 0.125 ms | 3.3 ms | 0.125 ms | 3.3 ms | 10 ms | 2.0 ms |
| 1 ms (1) | 0.05 ms | 0.33 ms | 0.05 ms | 0.33 ms | 0.02 ms | 0.05 ms |
| **0.1 ms (10)** | **3.3e-3 ms** | 1.25e-2 ms | **1.25e-2 ms** | 1.25e-2 ms | 3.3e-4 ms | 2.0e-3 ms |

**Interpretation.**
- The toy reproduces the notebook's S2C values: e.g. at N = 1 for T = 0.1 the figure reads ≈0.80 / 0.245 / 0.19, against 0.796 / 0.241 / 0.193 here.
- The unmodified simulator agrees with the toy log scheme within σ at every τ, confirming that evolve.jl:557 implements the logarithmic scheme.
- The log scheme is the flattest of the three and stays flat to beyond τ = 10⁻³ ms at θ = 10 → **PASS**.

---

## Run P1.2.6: single-spin trace

Script [p1_2_6_single_spin_trace.jl](baseline/p1_2_6_single_spin_trace.jl) → [results/baseline/p1_2_6_single_spin_trace/](results/baseline/p1_2_6_single_spin_trace/)

| Field | Value |
|---|---|
| Geometry | single cylinder R = 1 µm (S/V = 2 µm⁻¹), θ_relax = 0.1 |
| Physics | D = 3; R2 = 1/80; R1 = 1/1000 kHz |
| Sequences | SpinEcho TE = 5 ms (a) instant pulses, (b) hard pulses 0.5 / 1 ms |
| τ | fixed 0.01 ms (`TimeStep(0.01, Inf)`) |
| N_spins / seed | 100 / `Random.seed!(20260929)` |
| Outcome | no-RF step prediction max relative error 6.2e-16 (a) / 5.1e-16 (b); 112 / 121 collisions; **PASS** |

---

# Wall-geometry baselines for Phase 2 (B-series)

Run 2026-09-30 on HEAD `677a708`, `src/` unmodified (`provenance.src_modified = false`), same environment as above.

**Reference configuration** (Phase 2 tracker) used by B1 and B3:

| Field | Value |
|---|---|
| Geometry | `Walls(repeats=2)`: planes x = 2k µm, **normal along x**, both faces present, no surface parameters |
| Physics | D = 3 µm²/ms; R1 = R2 = 0; θ_relax = 0; permeability 0; off-resonance 0; no MT / sticking |
| Sequence | none: `mr.SequenceParts.empty_sequence()` (no RF, no gradients). Spins start with `transverse=1, longitudinal=0, phase=0` via `Snapshot` keywords, i.e. free decay from t = 0 |
| Spins | uniform in a ±500 µm box (`Snapshot(N, sim, 500)`); `Random.seed!(seed)` before each snapshot, so the same seed gives the same initial positions in every configuration |
| Default τ_max | **0.04 ms**, set by tortuosity 0.03·w²/D with size_scale = w = 2 µm; step length √(2Dτ) = 0.49 µm |

## Run B0: existing test suite

Output → [results/baseline/b0_test_suite/](results/baseline/b0_test_suite/) (`run.log`, `summary.md`)

| Field | Value |
|---|---|
| Command | `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["collisions","evolve","permeability","transfer","known_sequences"])'` |
| Why these suites | They cover the code the R₂(d) change will touch or depend on: wall collisions and reflection, the evolve loop, permeable surfaces (Route A), surface relaxation / MT (`transfer`, including `θ_relax` on walls) and known-signal sequences |
| Result | **442 / 442 pass**, 0 fail / error / broken. Test time 4 min 36 s (collisions 55 s, evolve 6 s, permeability 25 s, transfer 55 s, known_sequences 136 s); 9 min 23 s wall clock including precompilation |

Not run: `meshes`, `offresonance`, `radio_frequency`, `hierarchical_mri`, `various`, `subsets`, `swc`, `plots`, `cli`. They do not touch wall relaxation, and `plots` needs a display.

## Run B1: wall reference configuration, null run

Script [b1_walls_null.jl](baseline/b1_walls_null.jl) → [results/baseline/b1_walls_null/](results/baseline/b1_walls_null/)

**What / why.** The reference configuration with nothing that relaxes. The unmodified simulator must return M⊥(t)/M⊥(0) = 1 and phase 0 exactly. This is the reference V0b (P2.3.1) compares against at λ = 0. The runtime is the reference for P2.2.8 (the λ = 0 short-circuit must add no measurable cost) and the starting point of the cost study (P2.6.8).

| Field | default τ | τ = 1e-3 ms |
|---|---|---|
| τ_max / binding | 0.04 ms / tortuosity | 0.001 ms / user-set |
| Readouts | 0, 10, 20, 30, 40, 50 ms | same |
| N_spins / seeds | 100 000 / 1–10 | 10 000 / 1–3 (runtime scaling only) |
| max \|M⊥/M⊥0 − 1\| | **0** | **0** |
| max \|phase\| | **0°** | **0°** |
| Runtime per seed (excl. first, which includes compilation) | 1.62 s | 8.10 s |
| µs per spin per ms of sequence | 0.324 | 16.2 |
| µs per spin-step | 0.0130 | 0.0162 |
| Outcome | **PASS** | **PASS** |

**Interpretation.** With every relaxation parameter zero, the old code leaves M⊥ and its phase untouched, bit for bit. Any non-zero deviation at λ = 0 after the modification is therefore introduced by the new code. Cost is about 0.013–0.016 µs per spin-step on 8 threads. The small-τ run has 10× fewer spins, so it uses the threads less efficiently and its per-step cost is slightly higher. Estimate for later sweeps: a run costs about N_spins × (T/τ) × 0.015 µs, e.g. 10⁵ spins × 50 ms at τ = 1e-5 ms ≈ 12 min per seed.

## Run B2: bulk R2 in the wall reference configuration

Script [b2_walls_bulk_r2.jl](baseline/b2_walls_bulk_r2.jl) → [results/baseline/b2_walls_bulk_r2/](results/baseline/b2_walls_bulk_r2/)

**What / why.** The proposed model is R2_total = R2_bulk + ΔR2(d) (proposal Eq. 5). This run fixes how the unmodified simulator applies R2_bulk on its own, so Phase 2 can show the new code leaves that term intact. The old code multiplies M⊥ by exp(−R2·Δt) once per free sub-segment ([relax.jl:60-63](../src/relax.jl#L60-L63)), and bulk R2 does not enter τ_max.

| Field | Value |
|---|---|
| Configuration | reference config (B1) with global `R2` = R2_bulk |
| R2_bulk | 1/80 (T₂ = 80 ms), 0.025, 0.05, 0.1 ms⁻¹ |
| τ_max | 0.04 ms for every R2_bulk (unchanged from B1) |
| N_spins / seeds | 100 000 / 1–10; readouts 0, 10, …, 50 ms |
| Criteria | ensemble \|S/exp(−R2·t) − 1\| ≤ (n_steps + N)·ε = 2.2 × 10⁻¹¹; per spin ≤ 10⁻¹³; phase ≤ 10⁻⁹°; spin positions at 50 ms `==` the R2 = 0 run (seed 1) |

| R2_bulk (ms⁻¹) | Ensemble max rel. deviation | Per-spin max rel. deviation | Phase | Positions identical to R2 = 0 | Runtime vs R2 = 0 | Outcome |
|---|---|---|---|---|---|---|
| 0.0125 | 2.5 × 10⁻¹² | 2.1 × 10⁻¹⁴ | 0 | yes | 1.00 | **PASS** |
| 0.025 | 1.9 × 10⁻¹² | 2.0 × 10⁻¹⁴ | 0 | yes | 1.03 | **PASS** |
| 0.05 | 2.6 × 10⁻¹² | 2.6 × 10⁻¹⁴ | 0 | yes | 1.06 | **PASS** |
| 0.1 | 2.6 × 10⁻¹² | 3.0 × 10⁻¹⁴ | 0 | yes | 1.04 | **PASS** |

**Criterion change (recorded).** The first run used an ensemble tolerance of 10⁻¹² and reported FAIL at 2.5 × 10⁻¹². Diagnosis (seed 1, R2 = 1/80, t = 50 ms): every spin matches exp(−R2·t) to ≤ 2 × 10⁻¹⁴ (the 54 distinct per-spin values come from different numbers of sub-segments per step, i.e. repeated multiplication); the exact mean of the per-spin values deviates by 1.5 × 10⁻¹⁴; MCMR's own ensemble summation (`SpinOrientationSum`) returns 1.8 × 10⁻¹². The 10⁻¹² tolerance was below the floating-point roundoff of that summation. The criterion was replaced by the roundoff bound (n_steps + N)·ε plus the per-spin check, and the run repeated.

**Interpretation.**
- Bulk R2 is applied exactly (to roundoff) and does not change the phase.
- It does not change the random trajectories: with the same seed, positions are bit-identical to R2 = 0. In Phase 2, bulk-on and bulk-off runs with the same seed can therefore be compared spin by spin, and the additivity check S(bulk + layer) = exp(−R2_bulk·t)·S(layer) can be tested to roundoff, not to 2σ.
- Runtime differences (0–6%) are run-to-run timing noise; no criterion is attached to them.

## Run B3: density and layer occupancy between walls

Script [b3_walls_density.jl](baseline/b3_walls_density.jl) → [results/baseline/b3_walls_density/](results/baseline/b3_walls_density/) (`occupancy.csv`, `histogram.csv`, `histogram_chi2.csv`, figure `b3_walls_density.png`)

**What / why.** The layer relaxes a spin only while it is within λ (= h) of a face, so the relaxation it produces depends on how much time spins spend near the wall, and therefore on the density there. If the old code's reflection piles spins up at the wall or depletes them, every layer result inherits that distortion, whatever the new code does. Reflection errors show up within about one step length of a wall. At the default τ that step length (0.49 µm) is longer than every λ tested. This run checks, on the old code and with relaxation off, that (a) the density profile across the gap is flat and (b) the fraction within λ of each face is λ/w. It is the baseline for V9 (P2.3.2) and the occupancy check (P2.3.3, the Phase 2 gate).

| Field | Value |
|---|---|
| Observable | u = x mod w ∈ [0, 2) µm; distance to lower face u, to upper face w − u |
| λ | 0.05, 0.1, 0.2, 0.4 µm (f = 2.5, 5, 10, 20%) |
| Histogram | 100 bins of 0.02 µm, pooled over seeds (10⁴ spins expected per bin, Poisson σ = 1%) |
| N_spins / seeds | 100 000 per seed / 1–10 (10⁶ pooled) |
| Pooled occupancy resolution (1σ, relative) | 0.62% (λ = 0.05), 0.44% (0.1), 0.30% (0.2), 0.20% (0.4) |
| Criteria | occupancy \|z\| ≤ 2 with z = (n − Np)/√(Np(1−p)); profile χ² vs uniform \|z_χ²\| ≤ 3 and the bin touching each wall within 3σ |

| τ (ms) | Step √(2Dτ) (µm) | Readouts (ms) | Profile χ²/dof range | Wall-bin z range | Occupancy outside 2σ | Worst occupancy z | Largest rel. deviation (t > 0) | Runtime (10 seeds) | Outcome |
|---|---|---|---|---|---|---|---|---|---|
| 0.04 (default) | 0.49 | 0–50 | 0.82–1.27 | −1.92 … +2.61 | 2 / 48 | +2.20 | +1.37% (λ = 0.05, lower, t = 30) | 19 s | **PASS** |
| 1e-2 | 0.245 | 0–50 | 0.78–1.10 | −1.33 … +0.61 | 1 / 48 | −2.23 | −1.39% (λ = 0.05, lower, t = 10) | 58 s | **PASS** |
| 1e-3 | 0.078 | 0–50 | 0.89–1.21 | −1.34 … +1.16 | 0 / 48 | −1.85 | −1.09% (λ = 0.05, upper, t = 50) | 484 s | **PASS** |
| 1e-4 | 0.025 | 0, 1, 2, 5 | 0.75–1.16 | −1.40 … +1.16 | 0 / 32 | +1.84 | +1.15% (λ = 0.05, lower, t = 5) | 469 s | **PASS** |

The τ = 1e-4 ms run stops at 5 ms to keep the runtime reasonable. Spins cross the gap in w²/2D ≈ 0.7 ms, so 5 ms is about 7 mixing times. That is enough to show whether a near-wall artefact builds up.

**Interpretation.**
- The density profile is flat at every τ and readout, including the bins touching the walls. The largest bin deviations are within the 2σ band (±2%).
- Occupancy equals λ/w at every λ, face, readout and τ. The 3 of 176 comparisons outside 2σ are fewer than the ≈ 8 expected by chance, and they have no consistent sign. The seed-to-seed scatter of the occupancy fraction matches the binomial prediction (ratio 0.8–1.1), so the spins behave as independent uniform samples.
- **Therefore:** the unmodified reflection code produces no density artefact near walls, even when a single step (0.49 µm) is longer than the layer. The detection limit is about 1.2% (2σ) in occupancy at λ = 0.05 µm, and about 2% per 0.02 µm bin in the profile. If P2.3.2 or P2.3.3 fails after the modification, the cause is in the new code.

## Run B4: baseline reference curves S_θ(t)

Script [b4_theta_reference.jl](baseline/b4_theta_reference.jl) → [results/baseline/b4_theta_reference/](results/baseline/b4_theta_reference/) (`reference_curves.csv`, `timestep_check.csv`, `bulk_check.csv`, figure `b4_theta_reference.png`)

**What / why.** The R₂(d) model (R2_total = R2_bulk + ΔR2(d), ΔR2(d) = (ρ/h)·g(d/h)) is fitted to the baseline θ_relax model: for each θ, Phase 2 finds the layer thickness h\*(θ) (fixed profile g, fixed surface excess rate ΔR2(0)) whose signal best matches S_θ(t), and checks that attenuation increases with h. The baseline is a **reference, not the true answer**, so this run records reference curves with Monte Carlo σ and compares them with no analytic formula. θ = 0 is the no-surface-relaxation reference (identical to B1). The run also checks that the reference does not depend on the timestep and that bulk R2 multiplies it exactly.

| Field | Value |
|---|---|
| Geometry | `Walls(repeats=2, surface_relaxation=θ)`: θ acts on **both faces** of every gap (the R₂(d) layer will also be on both faces) |
| Physics | D = 3; R1 = 0; R2_bulk = 0 (reference) and 1/80 ms⁻¹ (factorisation check); no permeability / off-resonance / MT |
| Sequence | none (no RF); spins start transverse = 1; readouts 0, 5, …, 50 ms |
| θ | 0, 0.001, 0.0015, 0.002, 0.003, 0.005, 0.007, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05 (range set by a 1-seed pilot so S(50) covers ≈ 0.95 → 0.09) |
| τ | default 0.04 ms and 1e-2 ms for every θ (10 seeds); 1e-3 ms for θ = 0.001, 0.005, 0.02, 0.05 (5 seeds). Surface-relaxation τ bound 0.01/θ² ≥ 4 ms, never binding |
| N_spins / seeds | 100 000 per seed; seed k gives the same initial positions in every run |
| Reference curve | **τ = 1e-2 ms**, 10 seeds, all θ |
| Runtime | ≈ 45 min total (default 17 s, 1e-2 57 s, 1e-3 ≈ 250 s per θ) |

**Reference curves** (τ = 1e-2 ms, mean ± σ over 10 seeds; all 11 readouts in `reference_curves.csv`):

| θ_relax | S(10 ms) | S(25 ms) | S(50 ms) |
|---|---|---|---|
| 0 | 1.00000 ± 0.00000 | 1.00000 ± 0.00000 | 1.00000 ± 0.00000 |
| 0.001 | 0.99028 ± 0.00000 | 0.97587 ± 0.00001 | 0.95232 ± 0.00001 |
| 0.0015 | 0.98545 ± 0.00001 | 0.96402 ± 0.00001 | 0.92934 ± 0.00001 |
| 0.002 | 0.98065 ± 0.00001 | 0.95232 ± 0.00001 | 0.90692 ± 0.00002 |
| 0.003 | 0.97112 ± 0.00001 | 0.92935 ± 0.00002 | 0.86369 ± 0.00002 |
| 0.005 | 0.95234 ± 0.00002 | 0.88507 ± 0.00003 | 0.78335 ± 0.00004 |
| 0.007 | 0.93393 ± 0.00003 | 0.84291 ± 0.00004 | 0.71051 ± 0.00005 |
| 0.01 | 0.90700 ± 0.00004 | 0.78344 ± 0.00005 | 0.61379 ± 0.00006 |
| 0.015 | 0.86386 ± 0.00006 | 0.69358 ± 0.00006 | 0.48106 ± 0.00007 |
| 0.02 | 0.82282 ± 0.00008 | 0.61410 ± 0.00007 | 0.37713 ± 0.00007 |
| 0.03 | 0.74660 ± 0.00011 | 0.48160 ± 0.00009 | 0.23196 ± 0.00006 |
| 0.04 | 0.67758 ± 0.00013 | 0.37789 ± 0.00009 | 0.14281 ± 0.00005 |
| 0.05 | 0.61507 ± 0.00015 | 0.29666 ± 0.00009 | 0.08802 ± 0.00004 |

**Timestep check** (difference vs τ = 1e-2, in units of the combined SEM; max over readouts):

| Comparison | θ | max \|z\| | S(50) relative difference | Pattern |
|---|---|---|---|---|
| τ = 1e-3 vs 1e-2 | 0.001, 0.005, 0.02, 0.05 | 1.80, 1.79, 1.75, 1.66 | −2.8e-6, −1.3e-5, −3.0e-5, +2.1e-5 | no consistent sign → **converged** |
| default 0.04 vs 1e-2 | 0.001 → 0.05 | 1.62 → 2.74 (exceeds 2 from θ = 0.02) | +3.7e-6 → **+5.5e-4** | always positive, grows with θ → **small systematic bias** |

Caveat: runs at different τ share initial positions (same seed), so they are positively correlated and the combined SEM somewhat overestimates the noise of the difference. The z values are therefore conservative, which makes the default-τ trend more certain, not less.

**Bulk factorisation:** S(θ, R2_bulk = 1/80) / [S(θ, 0)·e^(−t/80)] − 1 ≤ 2.5 × 10⁻¹² for every θ and readout (bound 2.2 × 10⁻¹¹) → **PASS**. Adding R2_bulk in the old code multiplies the surface-relaxation signal exactly.

**Interpretation.**
- The reference set is ready: S_θ(t) at τ = 1e-2 ms, with σ ≤ 1.5 × 10⁻⁴ at every point. It covers S(50) from 1 (θ = 0) down to 0.088 (θ = 0.05), so any model attenuation in that range can be mapped to an equivalent θ by interpolation, and any θ in the grid can be used as a target for fitting h.
- At the default τ = 0.04 ms the baseline slightly under-attenuates (S too high by up to 0.055% at θ = 0.05). This is the baseline's own timestep dependence (the per-collision loss e^(−θ√τ) is not exactly τ-independent). It is below 3 SEM, but it has a consistent sign, so the default-τ curves are **not** used as the reference.
- For Phase 2 fitting, compare model runs against the τ = 1e-2 reference and use the same seeds (1–10) and spin count for pairing.
- Signal noise is small here (σ/S ≈ 10⁻⁵–5 × 10⁻⁴), so h\* will be determined precisely. Differences in curve *shape* between the two models (whether one h matches all readouts) will be resolvable.

## Run B5 (slim): runtime vs timestep

Script [b5_runtime_vs_tau.jl](baseline/b5_runtime_vs_tau.jl) → [results/baseline/b5_runtime_vs_tau/](results/baseline/b5_runtime_vs_tau/) (`runtime.csv`)

**What / why.** The R₂(d) timestep constraint τ ≤ c_h·h²/D (proposal Eq. 11) can require very small τ for thin layers. This measures the per-step cost of the unmodified code so that (a) the feasible τ range of the Phase 2 timestep study can be planned and (b) the new code's cost can later be quoted relative to the old one at the same τ. Runtime only; signal vs τ is covered by B4. (The full B5 originally planned, with signal over τ = 1e-5 – 1e-1 ms, was reduced to this after B4 showed the reference converged at τ = 1e-2 ms.)

| Field | Value |
|---|---|
| Configuration | `Walls(repeats=2, surface_relaxation=θ)`, θ ∈ {0, 0.01}; D = 3; no RF; fixed user τ |
| Workload | 10 000 spins × 20 000 steps per spin (simulated time T = 20 000·τ); one readout at T |
| Timing | median of 3 timed runs after one untimed warm-up; `julia -t 8` |

| τ (ms) | µs per spin-step, θ = 0 | µs per spin-step, θ = 0.01 | µs per spin per ms | One seed, 1e5 spins × 50 ms | 10 seeds |
|---|---|---|---|---|---|
| 1e-1 | 0.0265 | 0.0206 | 0.2–0.27 | ≈ 1 s | ≈ 15 s |
| 1e-2 | 0.0169 | 0.0172 | 1.7 | ≈ 9 s | ≈ 1.5 min |
| 1e-3 | 0.0163 | 0.0163 | 16 | ≈ 81 s | ≈ 14 min |
| 1e-4 | 0.0155 | 0.0157 | 155 | ≈ 13 min | ≈ 2.2 h |
| 1e-5 | 0.0153 | 0.0155 | 1 530 | ≈ 2.1 h | ≈ 21 h |

**Interpretation.**
- For τ ≤ 1e-2 ms the cost per spin-step is flat at ≈ 0.016 µs, so runtime is simply N_spins × (T/τ) × 0.016 µs. Surface relaxation (θ = 0.01) adds ≤ 2%.
- At τ = 0.1 ms the per-step cost is higher (0.02–0.027 µs): each step (√(6Dτ) ≈ 1.3 µm) often hits a wall, and collision handling dominates. The θ = 0 vs 0.01 difference at that τ is within run-to-run timing variation (min–max 4.8–5.4 s).
- **Feasibility for Phase 2 (old-code cost; the new code will add to it):** at the reference size (1e5 spins, 50 ms, 10 seeds), τ = 1e-3 ms takes ≈ 14 min and τ = 1e-4 ms ≈ 2 h, both routine. τ = 1e-5 ms takes ≈ 21 h, feasible only with fewer spins, fewer seeds or a shorter sequence. The tracker's V1 sweep down to τ = 1e-7 ms (≈ 9 days per seed at full size) is not feasible at full size.
- With Eq. 11, the τ needed depends on c_h (still to be determined in P2.6). For example, h = 0.1 µm gives h²/D = 3.3 × 10⁻³ ms; c_h = 0.03 → τ = 1e-4 ms (13 min per seed), c_h = 0.003 → τ = 1e-5 ms (2 h per seed).

---

## P1.4.2 Noise floor σ per baseline configuration (≥ 10 seeds)

| Configuration | Observable | σ across seeds | Relative σ | Spins per seed |
|---|---|---|---|---|
| Rep1 walls | R₂ fit | 4.4 × 10⁻⁶ ms⁻¹ | 0.045% | 10 000 |
| Rep1 walls | S(100 ms) | see `summary.json` `sigma_signal` | — | 10 000 |
| Rep1 cylinder | R₂ fit | 5.3 × 10⁻⁶ ms⁻¹ | 0.027% | ≈ 7 860 inside |
| Rep2 intra | T₂(TE = 50) | 2.0 ms | 2.3% | ≈ 348 |
| Rep2 extra | T₂(TE = 50) | 1.0 ms | 1.9% | ≈ 305 |
| Rep2 myelin | T₂(TE = 50) | 0.7 ms | 1.6% | ≈ 347 |
| Rep3 ordered_high | ADC plateau | 0.013 µm²/ms | 1.0% | ≈ 2 880 extra |
| Rep3 ordered_low | ADC plateau | 0.013 | 0.6% | ≈ 6 500 |
| Rep3 random_high | ADC plateau | 0.009 | 1.1% | ≈ 3 030 |
| Rep3 random_low | ADC plateau | 0.013 (geometry-to-geometry: 0.030) | 0.6% (1.5%) | ≈ 6 620 |
| Rep4 MCMR, T = 10 / 1 / 0.1 | attenuation | ≤ 0.0017 / ≤ 0.0027 / ≤ 0.0053 | ≤ 0.3% / 0.5% / 1.6% | 10 000 |

Per-τ and per-TE σ are stored in each `summary.json` (`*_sigma` arrays).

## P1.4.3 Raw outputs

All raw outputs are in `research/results/baseline/<run>/`. They are regenerable from the scripts with the seeds listed above. Each `summary.json` carries the provenance block (date, Julia and package versions, git branch and commit, `src_modified`, CPU, memory, threads, script name).

## Open items (Phase 1 exit)

- **P1.4.4 tag `baseline-v1`:** created 2026-09-30 on commit `06eccf3`.
- **Rep1 criterion:** decided 2026-10-02: judged against the exact rate within 2 SEM; both geometries pass (see Run Rep1).
- **P1.1.1 version pin:** decided 2026-10-02 to keep MCMRSimulator 1.1.0 (this repo, the code Phase 2 modifies) instead of the paper's v1.0.0. Rep1–Rep4 agree with the v1.0.0 notebook outputs, so no behavioural difference was found for these features.
- **Rep5** (optional) not run.
