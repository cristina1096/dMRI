# Phase 1 · Baseline Record (P1.3, P1.4)

Baseline results of the **unmodified** MCMRSimulator. Each run lists every setting needed to regenerate it, so that any later signal change can be attributed to the R₂(d) modification rather than to the setup.

- Code map and trace: [phase1_code_log.md](phase1_code_log.md)
- Scripts: [research/baseline/](baseline/)
- Raw outputs: [research/results/baseline/](results/baseline/). Each run folder holds CSV data, `summary.json` with provenance, `run.log` and a figure.

## Summary

| ID | Result | Reference | Criterion | Outcome |
|---|---|---|---|---|
| P1.3.1 | **R1** surface relaxation, slab | Eq. 17 | within 2σ | **PASS**: −1.87σ (vs exact Brownstein–Tarr: +0.53σ) |
| P1.3.1 | **R1** surface relaxation, single cylinder | Eq. 17 | within 2σ | **Not met against Eq. 17 (−3.15σ). Cause identified:** Eq. 17 is the fast-diffusion limit and exceeds the exact Brownstein–Tarr rate by 0.081%. At 7.9 × 10⁴ spins/config, σ = 0.027% resolves that gap. Against the exact rate: **−0.17σ, PASS** |
| P1.3.2 | **R2** compartment T₂ (spin echo, MT = 5e-3) | Fig. 8 (notebook `se_grid_plot`) | 5% (read off figure) | **PASS**: every TE compared within 5% of the digitised notebook figure (max 4.9%, mean \|Δ\| ≤ 2.4%), with and without myelin susceptibility |
| P1.3.3 | **R3** timestep plateau, ADC vs τ | Supp. Fig. S1 | plateau onset near τ ≈ 0.01 ms | **PASS**: onset (2σ) 6.3 × 10⁻³ – 2.5 × 10⁻² ms in all 4 packings |
| P1.3.4 | **R4** correction schemes, 1D slab | Supp. Fig. S2C | log scheme flat to τ ≈ 10⁻³ ms; record θ_relax | **PASS**: flat within 2σ to 3.3 × 10⁻³ ms (toy) and 1.25 × 10⁻² ms (MCMR) at θ = 10; MCMR matches the toy log scheme at every τ |
| P1.3.5 | R5 (optional) | Fig. 5 | — | not run |

**Two findings about the paper's figure notebooks** (they affect how references are read):
1. **Supp. Fig. S1 x-axis mislabelled.** `turtuosity.ipynb` computes ADC at τ = 10^(−3:0.2:1) (cell 17) but plots it against τ = 10^(−4:0.25:1) (cell 18). The saved figure is therefore shifted: its first point is τ = 10⁻³ ms, not 10⁻⁴. R3 is compared at the true τ.
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
# R2 figure comparison (needs a clone of mcmr_paper_figures):
python3 research/baseline/r2_digitise_notebook_fig8.py <mcmr_paper_figures>/Figure_8_9/gradient_spin_echo.ipynb
```
`NSPINS`, `NSEEDS` (and `NSPINS_TOY`, `NSPINS_MCMR`, `NGEOM`) can be overridden by environment variables. The values below are the defaults.

---

## Run R1: surface relaxation (P1.3.1)

Script [r1_surface_relaxation.jl](baseline/r1_surface_relaxation.jl) → [results/baseline/r1_surface_relaxation/](results/baseline/r1_surface_relaxation/)

| Field | Slab | Single cylinder |
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

**Interpretation.** The simulator reproduces surface relaxation to 0.03% of the exact solution. The cylinder "failure" against Eq. 17 is a real, resolvable difference between Eq. 17 (the first-order fast-diffusion approximation, κ = ρa/D = 3.3 × 10⁻³) and the exact eigenvalue. It is not a simulator error. For Phase 2's V4 (λ → 0 limit), compare against the exact Brownstein–Tarr rate, or choose κ small enough that the Eq. 17 bias is below σ.

---

## Run R2: compartment T₂ in myelinated white matter (P1.3.2)

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

## Run R3: timestep plateau of extra-axonal ADC (P1.3.3)

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

## Run R4: surface-relaxation correction schemes, 1D slab (P1.3.4)

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

## P1.4.2 Noise floor σ per baseline configuration (≥ 10 seeds)

| Configuration | Observable | σ across seeds | Relative σ | Spins per seed |
|---|---|---|---|---|
| R1 slab | R₂ fit | 4.4 × 10⁻⁶ ms⁻¹ | 0.045% | 10 000 |
| R1 slab | S(100 ms) | see `summary.json` `sigma_signal` | — | 10 000 |
| R1 cylinder | R₂ fit | 5.3 × 10⁻⁶ ms⁻¹ | 0.027% | ≈ 7 860 inside |
| R2 intra | T₂(TE = 50) | 2.0 ms | 2.3% | ≈ 348 |
| R2 extra | T₂(TE = 50) | 1.0 ms | 1.9% | ≈ 305 |
| R2 myelin | T₂(TE = 50) | 0.7 ms | 1.6% | ≈ 347 |
| R3 ordered_high | ADC plateau | 0.013 µm²/ms | 1.0% | ≈ 2 880 extra |
| R3 ordered_low | ADC plateau | 0.013 | 0.6% | ≈ 6 500 |
| R3 random_high | ADC plateau | 0.009 | 1.1% | ≈ 3 030 |
| R3 random_low | ADC plateau | 0.013 (geometry-to-geometry: 0.030) | 0.6% (1.5%) | ≈ 6 620 |
| R4 MCMR, T = 10 / 1 / 0.1 | attenuation | ≤ 0.0017 / ≤ 0.0027 / ≤ 0.0053 | ≤ 0.3% / 0.5% / 1.6% | 10 000 |

Per-τ and per-TE σ are stored in each `summary.json` (`*_sigma` arrays).

## P1.4.3 Raw outputs

All raw outputs are in `research/results/baseline/<run>/`. They are regenerable from the scripts with the seeds listed above. Each `summary.json` carries the provenance block (date, Julia and package versions, git branch and commit, `src_modified`, CPU, memory, threads, script name).

## Open items (Phase 1 exit)

- **P1.4.4 tag `baseline-v1`:** not created yet. It needs a decision on the R1 cylinder criterion (Eq. 17 vs exact Brownstein–Tarr) and a commit of `research/`.
- **P1.1.1 version pin:** the baseline was taken on MCMRSimulator 1.1.0 (this repo), not the paper's v1.0.0. R1–R4 agreement with the v1.0.0 notebook outputs suggests no behavioural difference for these features.
- **R5** (optional) not run.
