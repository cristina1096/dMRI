# Phase 2 · Record (R₂(d) on walls, step profile, no RF)

Results of the modified simulator (branch `cc/near-surface-r2`, layer implemented by plan `docs/superpowers/plans/2026-10-01-near-surface-layer-walls.md`).
Baselines are in [baseline_record.md](baseline_record.md). The baseline θ_relax model is a **reference**, not the true answer.

- Code: [research/phase2/](phase2/) (harness, expectations, scripts); tests: `julia --project=research/baseline -t 8 research/phase2/test/runtests.jl`
- Raw outputs: [research/results/phase2/](results/phase2/), each with `summary.json` (provenance: git commit, versions, hardware)

**Reference configuration:** `Walls(repeats=2)` (w = 2 µm), D = 3 µm²/ms, no RF, spins start transverse in a ±500 µm box, readouts 0, 5, …, 50 ms, τ = 1e-2 ms, seeds 1–10, 10⁵ spins per seed, step-profile layer on both faces, ΔR₂(0) = 0.1 ms⁻¹ (ρ = ΔR₂(0)·h), R₂_bulk = 0 unless stated. Uncertainty of a mean = SEM over seeds.

## Summary

| ID | Test | Criterion | Outcome |
|---|---|---|---|
| P2.2.9 | Single-spin trace with the layer on | every step: \|M⊥_sim/M⊥_pred − 1\| ≤ 1e-6; ≥ 1 reflection inside the layer | **PASS**: max error 5.0e-08; 4848 steps with a reflection inside the layer |
| P2.1.3 | Harness reproduces B1, B2 and one B4 curve (same seeds) | exact / CSV precision | **PASS**: B1 = 1 exactly, B2 = e^(−t/80) to 1e-11, B4 θ = 0.01 curve to 1e-8 (`research/phase2/test/`) |
| P2.1.5 | Exact expectations and B4 interpolation | unit-tested before any run | **Done**: `research/phase2/expected.jl` (21 tests) |
| P2.1.4 | Noise floor and spin count | neighbouring h in the sweep differ by > 10 SEM | σ(50) = 4.8e-05 (h = 0.1), 1.6e-05 (h = 0.5); worst neighbour separation 4287 SEM → **100000** spins per seed (no increase needed) |
| P2.3.1 | V0b: ρ = 0 vs no layer (bulk off / on), default τ | bit-identical positions and M⊥ | **PASS** / **PASS** |
| P2.3.4(a) | Exact end point h = w/2 (bulk off / on) | per spin ≤ 1e-10; ensemble ≤ roundoff bound | **PASS** / **PASS**: per spin ≤ 1.4e-13, ensemble ≤ 2.6e-12 |
| P2.3.11 | Additivity R₂_bulk + ΔR₂ (h = 0.2 µm) | same positions; per spin ≤ 1e-10; ensemble ≤ roundoff bound | **PASS**: per spin 6.6e-14, ensemble 7.5e-14 |
| P2.3.2 | V9 density with the layer on (τ = 0.04, 1e-2) | positions == no-layer run; histogram \|z_χ²\| ≤ 3 | **PASS** |
| P2.3.3 | Occupancy within h of either face (**gate**) | fraction = 2h/w, \|z\| ≤ 2 for h = 0.05–0.5 | **PASS** |
| P2.3.4(b) | Attenuation vs h (characterisation), τ = 1e-2 | recorded with SEM | S(50) from 1 (h = 0) to 0.00674 (h = 1.0); 100000 spins × 10 seeds |
| P2.3.5 | Monotonicity in h (proposal §3.2 VI) | every step > 2 combined SEM | **PASS**: 13 of 13 steps |
| P2.3.6 | Decay shape (curvature of ln S) | characterisation | curvature 0 within noise at every h (\|c\| ≤ 5e-08 ms⁻², ≤ 1.3 SEM) |

## P2.2.9: single-spin trace with the layer on

Script [p2_2_9_layer_trace.jl](phase2/p2_2_9_layer_trace.jl) → [results/phase2/p2_2_9_layer_trace/](results/phase2/p2_2_9_layer_trace/) (`trace_spin1.csv`: every step of spin 1)

**What / why.** The P1.2.6 trace repeated with the layer on. Each step's change in M⊥ is predicted from the simulator's own reflected path, with the layer time on every straight piece computed by brute-force sampling (10⁴ points per piece), independently of the `Layers` code. Agreement shows that distance, crossings, reflection handling, segmentation and the M⊥ update (P2.2.1–P2.2.7) work together.

| Field | Value |
|---|---|
| Configuration | walls every 2 µm; step layer on both faces, h = 0.5 µm, ΔR₂(0) = 0.1 ms⁻¹; R₂_bulk = 1/80 ms⁻¹; D = 3 µm²/ms; τ = 0.01 ms |
| Spins / steps | 100 spins × 500 steps (50 000 steps); `Random.seed!(20261003)` |
| Max relative error per step | **5.0e-08** (tolerance 1e-6; the sampling resolution is about 1e-7 per step, so this is sampling error) |
| Steps with a reflection inside the layer | 4848 |
| Runtime | ≈ 15 s (1 thread) |
| Outcome | **PASS** |

## P2.1.4: noise floor and spin count

Script [p2_1_4_noise_floor.jl](phase2/p2_1_4_noise_floor.jl) → [results/phase2/p2_1_4_noise_floor/](results/phase2/p2_1_4_noise_floor/) (`pilot_pairs.csv`)

| h (µm) | S(50) | σ(50) (one seed) | SEM(50) (10 seeds) | max σ over readouts | runtime (10 seeds) |
|---|---|---|---|---|---|
| 0.1 | 0.60686 | 4.78e-05 | 1.51e-05 | 6.58e-05 | 73 s |
| 0.5 | 0.08240 | 1.56e-05 | 4.93e-06 | 9.00e-05 | 79 s |

Pilot sweep (3 seeds, 10⁵ spins): the closest pair of neighbouring h values (0.02 → 0.03 µm) is separated by **4287 SEM**, far above the required 10, so the sweep keeps **100000 spins per seed**. The per-spin variance is small because every spin spends a similar fraction of 50 ms in the layer.

H_GRID coverage of the B4 range (1 − S(50) from 0.048 to 0.912): pilot 1 − S(50) = 0.0488 at h = 0.01 and 0.9176 at h = 0.5 → **covered**; h = 1.0 gives S(50) = 0.00674 (= e⁻⁵ = 0.00674, the exact end point).

## P2.3.1, P2.3.4(a), P2.3.11: exact checks at full size

Script [p2_3_exact_checks.jl](phase2/p2_3_exact_checks.jl) → [results/phase2/p2_3_exact_checks/](results/phase2/p2_3_exact_checks/)

These hold in any diffusion regime, so they are pass/fail to floating-point precision. Same seeds (1–10) and 100000 spins per seed throughout.

| Check | Configuration | Max deviation per spin | Max deviation ensemble | Bound | Outcome |
|---|---|---|---|---|---|
| V0b, bulk off | ρ = 0, h = 0.5 vs no layer; default τ = 0.04 ms | identical (==) | identical | — | **PASS** |
| V0b, bulk on (1/80) | same | identical (==) | identical | — | **PASS** |
| h = w/2, bulk off | h = 1.0, ΔR₂(0) = 0.1, τ = 1e-2 | 8.5e-14 | 2.6e-12 | 2.3e-11 | **PASS** |
| h = w/2, bulk on | same + R₂_bulk = 1/80 | 1.4e-13 | 1.9e-12 | 2.3e-11 | **PASS** |
| additivity | h = 0.2, ΔR₂(0) = 0.1, R₂_bulk 0 vs 1/80; positions identical: True | 6.6e-14 | 7.5e-14 | 2.3e-11 | **PASS** |

## P2.3.2, P2.3.3: density and occupancy with the layer on

Script [p2_3_2_density_occupancy.jl](phase2/p2_3_2_density_occupancy.jl) → [results/phase2/p2_3_2_density_occupancy/](results/phase2/p2_3_2_density_occupancy/)

The layer changes |Mxy| only, so final positions with the layer on (h = 0.5 µm, ΔR₂(0) = 0.1 ms⁻¹) are compared with the no-layer run (same seeds); the B3 statistics are then recomputed on the layer-on positions at t = 50 ms (10⁶ spins pooled).

| τ (ms) | Positions identical to no layer | χ²/dof (z) | Occupancy z for h = 0.05 / 0.1 / 0.2 / 0.4 / 0.5 | Outcome |
|---|---|---|---|---|
| 0.04 (default) | True | 0.948 (-0.36) | +0.08 / +0.09 / -1.03 / -0.90 / -0.58 | **PASS** |
| 1e-2 | True | 0.949 (-0.36) | -0.01 / +0.54 / -0.37 / -1.52 / -1.11 | **PASS** |

At the default τ the histogram statistic equals B3's value at t = 50 ms (χ²/dof 0.948), as expected: same seeds, same positions.

## P2.3.4(b), P2.3.5, P2.3.6: h sweep

Script [p2_3_4_h_sweep.jl](phase2/p2_3_4_h_sweep.jl) → [results/phase2/p2_3_4_h_sweep/](results/phase2/p2_3_4_h_sweep/) (`sweep.csv` per seed, figure `p2_3_4_h_sweep.png`). Step profile, both faces, ΔR₂(0) = 0.1 ms⁻¹, τ = 1e-2 ms (provisional until P2.6), 100000 spins × 10 seeds.

| h (µm) | S(25) | S(50) ± SEM | 1 − S(50) | curvature c ± SEM (ms⁻²) |
|---|---|---|---|---|
| 0.0 | 1.00000 | 1.00000 ± 0.0e+00 | 0.00000 | -0.0e+00 ± 0.0e+00 |
| 0.005 | 0.98758 | 0.97531 ± 1.9e-06 | 0.02469 | 1.7e-09 ± 1.4e-09 |
| 0.01 | 0.97531 | 0.95123 ± 3.9e-06 | 0.04877 | 3.7e-09 ± 2.8e-09 |
| 0.02 | 0.95124 | 0.90486 ± 6.1e-06 | 0.09514 | 6.5e-09 ± 6.1e-09 |
| 0.03 | 0.92777 | 0.86076 ± 7.7e-06 | 0.13924 | 7.9e-09 ± 8.6e-09 |
| 0.05 | 0.88256 | 0.77893 ± 9.8e-06 | 0.22107 | 1.1e-08 ± 1.4e-08 |
| 0.07 | 0.83957 | 0.70489 ± 1.2e-05 | 0.29511 | 1.7e-08 ± 1.9e-08 |
| 0.1 | 0.77900 | 0.60686 ± 1.5e-05 | 0.39314 | 2.7e-08 ± 2.5e-08 |
| 0.15 | 0.68764 | 0.47286 ± 1.5e-05 | 0.52714 | 2.7e-08 ± 3.7e-08 |
| 0.2 | 0.60701 | 0.36847 ± 1.6e-05 | 0.63153 | 3.0e-08 ± 4.8e-08 |
| 0.3 | 0.47300 | 0.22374 ± 1.2e-05 | 0.77626 | 4.5e-08 ± 6.7e-08 |
| 0.5 | 0.28705 | 0.08240 ± 4.9e-06 | 0.91760 | 4.5e-08 ± 7.1e-08 |
| 0.7 | 0.17401 | 0.03028 ± 1.6e-06 | 0.96972 | 3.3e-08 ± 5.5e-08 |
| 1.0 | 0.08208 | 0.00674 ± 2.9e-19 | 0.99326 | -4.5e-16 ± 0.0e+00 |

**Interpretation.**
- Attenuation rises with h at every step (13 of 13, each by thousands of SEM) → **P2.3.5 PASS** (proposal §3.2 VI: fixed profile and ΔR₂(0), larger h gives more attenuation).
- ln S(t) is a straight line at every h: the fitted curvature is zero within 1.3 SEM everywhere, so in this configuration (w = 2 µm, D = 3, up to 50 ms) the decay of the whole ensemble is single-exponential at this precision. h = 0 and h = w/2 = 1.0 are exactly single-exponential by construction.
- An empirical observation (not an assumption used anywhere): S(50) is close to exp(−ΔR₂(0)·(2h/w)·50 ms), i.e. ΔR₂(0) weighted by the measured layer occupancy 2h/w (P2.3.3):

| h (µm) | occupancy 2h/w | S(50) simulated | exp(−ΔR₂(0)·(2h/w)·50) |
|---|---|---|---|
| 0.005 | 0.0050 | 0.97531 | 0.97531 |
| 0.01 | 0.0100 | 0.95123 | 0.95123 |
| 0.02 | 0.0200 | 0.90486 | 0.90484 |
| 0.03 | 0.0300 | 0.86076 | 0.86071 |
| 0.05 | 0.0500 | 0.77893 | 0.77880 |
| 0.07 | 0.0700 | 0.70489 | 0.70469 |
| 0.1 | 0.1000 | 0.60686 | 0.60653 |
| 0.15 | 0.1500 | 0.47286 | 0.47237 |
| 0.2 | 0.2000 | 0.36847 | 0.36788 |
| 0.3 | 0.3000 | 0.22374 | 0.22313 |
| 0.5 | 0.5000 | 0.08240 | 0.08208 |

It agrees to 4–5 digits for thin layers and drifts slightly for thicker ones (0.08240 vs 0.08209 at h = 0.5). This is a description of this geometry and diffusivity only. The differences from the baseline are what V4 measures.
