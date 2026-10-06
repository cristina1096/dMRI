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
| P2.3.7 | V11: endpoint rule vs exact overlap (toy, h = 0.1) | bias quantified (sign not assumed); toy exact = MCMR within 2 SEM | bias at τ = 1e-3 / 1e-2 / 4e-2: -0.1 / -1.4 / +5.2 SEM; cross-check **PASS** (+1.8 SEM) |
| P2.3.12 | Side-specific layer: exact mirror check + characterisation | positive layer on −x walls == negative layer on +x walls (identical) | **PASS** (bit-identical); paired statistical difference reported (max 3.2 SEM, 1.3 SEM at 50 ms); S(50): both 0.36847, one side 0.60771, asymmetric 0.36834 |
| P2.3.9 | V4: h*(θ) fit to the B4 reference (τ = 1e-2, ΔR₂(0) = 0.1) | h* found for every θ, with uncertainty and residuals (no pass/fail) | **12 of 12 θ fitted**; h* from 0.00977 to 0.48681 µm; at h* the model matches the reference curve within 4e-05 relative (≤ 0.4 combined SEM) at every readout |
| P2.6.1 | V1: timestep plateau (toy, h = 0.1, ΔR₂(0) = 0.1, τ = 1e-5–0.1 ms) + MCMR at 4 τ | τ_conv at 1 % / 0.1 %; MCMR = toy within 2 SEM | 1 % never exceeded (τ_conv ≥ 0.1 ms); τ_conv(0.1 %) = 0.086 ms; bias at 0.1 ms −1.1e-3 (less attenuation); MCMR vs toy **PASS** (\|z\| ≤ 1.3) |
| P2.6.2, P2.6.3 | V2 (h²) and V3 (1/D) scaling of τ_conv; decision point (rate) | τ_conv·D/h² and τ_conv·D constant within ±30 % | **V2 FAIL, V3 FAIL** (trends opposite to h²/D; 1 % not reached at ΔR₂(0) = 0.1). τ_conv set by ΔR₂(0): ΔR₂(0)·τ_conv ≈ 0.05–0.06 (1 %), ≈ 0.003 (0.1 %, ΔR₂(0) ≥ 0.5); bias at 0.1 ms ∝ ΔR₂(0) |
| P2.6.5 | Layer constraint in the simulator's timestep | constraint active; binding term reported | **Done** (user decision): option `layer`, τ ≤ c/ΔR₂(0)_max, default c = 0.005; verbose message names it when it binds; suite 2894/2894 |
| P2.6.8, P2.6.10–P2.6.12 | Cost with the layer; consequences for τ = 1e-2 results | cost tabulated; earlier results judged | layer overhead ×1.2–1.4 vs B5 (×2.1 at h = 0.01); full run ≤ 0.11 h up to ΔR₂(0) = 2, 0.63 h at 10; all earlier τ = 1e-2 results (ΔR₂(0) = 0.1) biased ≤ 6.4e-4 in R → **no reruns needed** |

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

## P2.3.7: V11, endpoint rule vs exact overlap

Script [p2_3_7_endpoint.jl](phase2/p2_3_7_endpoint.jl) → [results/phase2/p2_3_7_endpoint/](results/phase2/p2_3_7_endpoint/). 1D toy walk identical to MCMR's x-motion between walls ([toy.jl](phase2/toy.jl)); h = 0.1 µm, ΔR₂(0) = 0.1 ms⁻¹, T = 50 ms, seeds 1–10. The **exact** rule integrates the layer time along every reflected straight piece (as MCMR does); the **endpoint** rule charges a whole step to the layer if the step's end point is inside it.

| τ (ms) | spins per seed | S(50) exact ± SEM | S(50) endpoint ± SEM | bias (endpoint − exact) | bias in SEM |
|---|---|---|---|---|---|
| 0.001 | 20000 | 0.60679 ± 2.6e-05 | 0.60678 ± 2.5e-05 | -4.9e-06 | -0.1 |
| 0.01 | 100000 | 0.60690 ± 1.4e-05 | 0.60687 ± 1.1e-05 | -2.5e-05 | -1.4 |
| 0.04 | 100000 | 0.60707 ± 1.9e-05 | 0.60721 ± 1.8e-05 | +1.4e-04 | +5.2 |

Cross-check: toy exact vs MCMR (P2.3.4(b), h = 0.1, τ = 1e-2): 0.60690 vs 0.60686 (+1.8 SEM) → **PASS** (independent random streams, so the check is statistical).

**Interpretation.**
- The endpoint rule's bias is small here and **changes sign with τ**: indistinguishable from zero at τ = 1e-3 ms, slightly negative (more attenuation) at 1e-2, and positive (less attenuation, +1.4 × 10⁻⁴) at 4e-2 ms. The tracker's expectation "endpoint under-counts" is therefore not confirmed as a general rule for this configuration; the bias depends on τ and is small because the end points of steps sample the layer occupancy without bias in a uniform ensemble.
- The exact rule itself also moves with τ (S(50) = 0.60679 → 0.60690 → 0.60707 for τ = 1e-3 → 1e-2 → 4e-2 ms, i.e. +4.6 × 10⁻⁴ relative over the range). That is the straight-line timestep effect of the proposal §3.2 V, to be quantified in V1 (P2.6). At τ = 1e-2 it is +1.8 × 10⁻⁴ relative to τ = 1e-3, larger than the 10⁻⁵ SEM, so results at τ = 1e-2 are not yet timestep-converged at this precision.

## P2.3.12: side-specific layer

Script [p2_3_12_sides.jl](phase2/p2_3_12_sides.jl) → [results/phase2/p2_3_12_sides/](results/phase2/p2_3_12_sides/) (first version kept in `first_pass/`). h = 0.2 µm, ΔR₂(0) = 0.1 ms⁻¹ unless stated; asymmetric: positive h = 0.2 / ΔR₂(0) = 0.1, negative h = 0.4 / ΔR₂(0) = 0.05 (same ρ = 0.02 µm/ms on both sides).

| Case | S(25) ± SEM | S(50) ± SEM |
|---|---|---|
| both sides | 0.60701 ± 2.0e-05 | 0.36847 ± 1.6e-05 |
| positive only | 0.77959 ± 2.4e-05 | 0.60771 ± 2.8e-05 |
| negative only | 0.77948 ± 2.4e-05 | 0.60765 ± 2.4e-05 |
| asymmetric | 0.60691 ± 1.9e-05 | 0.36834 ± 1.5e-05 |

**Exact mirror check → PASS.** A positive-side layer on walls whose normal points along −x is geometrically the same as a negative-side layer on walls along +x. With the same seeds the two runs give **bit-identical** signals at every readout (and identical positions), so the code treats the two sides exactly symmetrically.

**Statistical mirror comparison (reported, not graded).** Positive-only minus negative-only on the same +x walls, paired per seed, z at t = 5…50 ms: -0.8 / -0.4 / +0.3 / +2.1 / +2.7 / +3.2 / +2.2 / +2.3 / +2.0 / +1.3. The difference is at most 1.1 × 10⁻⁴ (relative 2 × 10⁻⁴), changes sign between early and late readouts, and is 1.3 SEM at 50 ms. Because the exact check rules out the code, this is a property of the finite spin sample carried through strongly correlated readouts (the same spins at every time).

**Criterion change (recorded).** The plan graded the statistical comparison as "|z| ≤ 2 at every readout", with unpaired SEMs. The first run gave max |z| = 3.15 → FAIL (kept in `first_pass/`). That criterion was mis-specified: the two runs share their spins, and 11 correlated readouts were each required to stay within 2σ. The pass criterion is now the exact flipped-normal identity, which tests the code directly.

**Characterisation.** Layers on both sides attenuate more than a layer on one side (0.368 vs 0.608 at 50 ms). The asymmetric layer, with the same ρ on both sides but different h and ΔR₂(0), gives almost the same S(50) as the symmetric one (0.36834 vs 0.36847, a 1.3 × 10⁻⁴ difference that is resolved at this precision). No baseline exists for asymmetric surfaces, because MCMR's θ_relax is the same on both sides.

## P2.3.9: V4, h*(θ) fitted to the baseline reference

Scripts [p2_3_9_fit_h.jl](phase2/p2_3_9_fit_h.jl) (start and parabola refinement) and [p2_3_9b_gauss_newton.jl](phase2/p2_3_9b_gauss_newton.jl) (final step) → [results/phase2/p2_3_9_fit_h/](results/phase2/p2_3_9_fit_h/) (`h_star_final.csv`, figure `p2_3_9_fit_h.png`; earlier stages in `first_pass/` and `parabola_pass/`). Step profile, both faces, ΔR₂(0) = 0.1 ms⁻¹, τ = 1e-2 ms, seeds 1–10, 10⁵ spins. The baseline is a reference, not the truth: residuals show where the two models differ.

**Method (three stages, recorded as rulings).**
1. Start: h from the P2.3.4(b) sweep, S(t; h) interpolated linearly in h. The residuals of direct runs at this h were 8–660 SEM. That was interpolation error between grid points, because the SEM (about 10⁻⁵) is far below the interpolation accuracy (`first_pass/`).
2. Parabola: χ² from direct runs at h·(1 − 5%, 1, 1 + 5%), vertex as the new h. The residuals still had one sign at every θ (up to 23 SEM), so this was not yet the minimum at this precision (`parabola_pass/`). The jackknife σ(h*) below comes from this stage.
3. One Gauss–Newton step (slope dS/dh from the sweep), then a final direct run. Residuals are now mixed in sign and below 0.4 SEM everywhere.

| θ_relax | h* (µm) | σ_jack (µm) | h*/θ | χ²/dof | max \|residual\| (SEM) | max \|S_model/S_ref − 1\| | S(50) reference | S(50) model |
|---|---|---|---|---|---|---|---|---|
| 0.001 | 0.00977 | 4e-06 | 9.772 | 0.05 | 0.4 | 1.9e-06 | 0.95232 | 0.95232 |
| 0.0015 | 0.01466 | 4e-06 | 9.772 | 0.04 | 0.4 | 2.3e-06 | 0.92934 | 0.92934 |
| 0.002 | 0.01954 | 8e-06 | 9.772 | 0.03 | 0.3 | 2.6e-06 | 0.90692 | 0.90692 |
| 0.003 | 0.02932 | 5e-06 | 9.773 | 0.02 | 0.3 | 2.5e-06 | 0.86369 | 0.86370 |
| 0.005 | 0.04887 | 2e-05 | 9.773 | 0.02 | 0.2 | 2.5e-06 | 0.78335 | 0.78335 |
| 0.007 | 0.06841 | 3e-05 | 9.773 | 0.03 | 0.4 | 5.3e-06 | 0.71051 | 0.71051 |
| 0.01 | 0.09772 | 5e-05 | 9.772 | 0.02 | 0.3 | 7.7e-06 | 0.61379 | 0.61380 |
| 0.015 | 0.14655 | 4e-05 | 9.770 | 0.01 | 0.1 | 6.7e-06 | 0.48106 | 0.48106 |
| 0.02 | 0.19534 | 3e-05 | 9.767 | 0.03 | 0.2 | 1.2e-05 | 0.37713 | 0.37713 |
| 0.03 | 0.29277 | 6e-05 | 9.759 | 0.03 | 0.2 | 1.7e-05 | 0.23196 | 0.23196 |
| 0.04 | 0.38994 | 3e-04 | 9.748 | 0.04 | 0.3 | 2.9e-05 | 0.14281 | 0.14281 |
| 0.05 | 0.48681 | 2e-04 | 9.736 | 0.05 | 0.3 | 4.0e-05 | 0.08802 | 0.08802 |

**Timestep check** (h* fixed at the parabola value, τ = 1e-2 → 1e-3 ms, S(50) shift): θ = 0.002: +0.38 SEM (+3.6e-06); θ = 0.01: -0.99 SEM (-3.8e-05); θ = 0.05: -1.00 SEM (-1.2e-04). At these h the fitted model moves by at most about 1 SEM when τ is reduced tenfold.

**Bulk on** (θ = 0.01, R₂_bulk = 1/80, at the parabola h*): the residuals equal the bulk-off residuals at the same h (max 15.79 in both), as expected from the exact factorisation (P2.3.11).

**Caveat.** Model and reference share seeds 1–10 (same initial positions), so √(SEM₁² + SEM₂²) overestimates the noise of their difference. That is why χ²/dof is far below 1. Residuals in SEM are therefore conservative.

**Interpretation.**
- For every θ in the B4 grid, a single layer thickness h* reproduces the **whole** baseline curve S_θ(t), 0–50 ms, to within a few 10⁻⁵ relative. In this configuration (walls 2 µm apart, D = 3 µm²/ms, up to 50 ms) the step-profile layer and the baseline θ_relax model cannot be told apart from the signal once h is fitted. Any difference in curve shape is below about 10⁻⁵.
- h* is almost exactly proportional to θ: h*/θ = 9.74–9.77 µm per unit θ over the whole range, drifting by 0.3% from the smallest to the largest θ. Equivalently, the fitted integrated relaxivity ρ* = ΔR₂(0)·h* = 0.974–0.977·θ in these units. This is an empirical result of the fit, reported as found.
- Whether the shape of the decay ever distinguishes the two models (narrower gaps, thicker layers, other profiles, curved geometry) is the subject of the shape comparison (P2.5) and Phase 3.

## Gate for section 2.4 (2026-10-03)

P2.2.9 single-spin trace: **PASS**. P2.3.3 occupancy: **PASS**. P2.3.4(a) exact end point h = w/2: **PASS**. → Section 2.4 (other profiles) may start. Characterisation results above are at τ = 1e-2 ms. The V4 timestep check moved the fitted model by at most about 1 SEM, while the V11 toy shows a τ-shift of the exact rule of order 10⁻⁴ relative. The timestep study (P2.6) settles convergence.

## P2.6.1: V1, timestep plateau

Toy 1D walk (exact segment rule, step profile, both faces, w = 2 µm, D = 3 µm²/ms), h = 0.1 µm, ΔR₂(0) = 0.1 ms⁻¹, T = 10 ms, 100000 spins × 10 seeds, τ = 1e-5 … 0.1 ms (3 per decade). Observable R = −ln S(T_eff)/T_eff; reference = τ = 1e-5. Script `research/phase2/p2_6_scan.jl`; output `research/results/phase2/p2_6_scan_V1_h0.1/`.

| τ (ms) | R (ms⁻¹) | SEM | R/R_ref − 1 |
|---|---|---|---|
| 1e-5 | 9.98899e-03 | 1.1e-06 | 0 |
| 2.15e-5 | 9.98741e-03 | 8.7e-07 | −1.6e-04 |
| 4.64e-5 | 9.98803e-03 | 1.5e-06 | −9.6e-05 |
| 1e-4 | 9.99133e-03 | 1.7e-06 | +2.3e-04 |
| 2.15e-4 | 9.99159e-03 | 1.3e-06 | +2.6e-04 |
| 4.64e-4 | 9.99037e-03 | 1.2e-06 | +1.4e-04 |
| 1e-3 | 9.99025e-03 | 1.5e-06 | +1.3e-04 |
| 2.15e-3 | 9.98823e-03 | 1.1e-06 | −7.6e-06 |
| 4.64e-3 | 9.98849e-03 | 1.3e-06 | +5.0e-05 |
| 1e-2 | 9.98642e-03 | 1.6e-06 | −2.6e-04 |
| 2.15e-2 | 9.98557e-03 | 1.6e-06 | −3.4e-04 |
| 4.64e-2 | 9.98424e-03 | 1.7e-06 | −4.8e-04 |
| 0.1 | 9.97774e-03 | 1.5e-06 | −1.13e-03 |

- τ_conv at 1 %: **≥ 0.1 ms** (no τ in the range fails; lower bound). τ_conv at 0.1 %: **0.086 ms**.
- Bias at the largest τ: −1.13e-3 relative, i.e. **less attenuation** at large τ, the sign the proposal (§3.2 V) expects.
- Points from 1e-5 to 1e-3 scatter by up to ±2.6e-4 (about 2 SEM): this is the precision floor of the scan, and it is below the 0.1 % tolerance.

**MCMR confirmation** (`p2_6_1_mcmr_check.jl`, same configuration, 100000 spins × 10 seeds, independent random stream):

| τ (ms) | R MCMR | R toy | z |
|---|---|---|---|
| 0.1 | 9.97614e-03 ± 2.0e-06 | 9.97774e-03 ± 1.5e-06 | −0.65 |
| 1e-2 | 9.98929e-03 ± 1.5e-06 | 9.98642e-03 ± 1.6e-06 | +1.31 |
| 1e-3 | 9.99003e-03 ± 1.1e-06 | 9.99025e-03 ± 1.5e-06 | −0.12 |
| 1e-4 | 9.99071e-03 ± 1.6e-06 | 9.99133e-03 ± 1.7e-06 | −0.27 |

→ **PASS**: the toy represents MCMR's τ dependence at this precision, so the dense scans below use the toy.

Sizing (ledger rulings): 100000 spins per seed (relative SEM of R ≈ 1.6e-4 < tol/3 at 0.1 %); reference τ = 1e-5 for V1 only and 1e-4 for the other scans, since R is flat from 1e-5 to 1e-3 (V1, and a probe at h = 0.01).

## P2.6.2, P2.6.3: scaling of τ_conv

Same toy, 100000 spins × 10 seeds, T = 10 ms, τ = 1e-4 … 0.1 ms, reference τ = 1e-4 (V1: 1e-5). Besides the planned V2/V3 scans, the probes during sizing showed no τ dependence for a thin layer but a growing one with the rate, so the tracker's decision point (τ_conv vs ΔR₂(0)) was measured too: ΔR₂(0) = 0.5, 1, 2 at h = 0.1, and h = 0.01 at ΔR₂(0) = 10 (same ρ = ΔR₂(0)·h as h = 0.1 at 1) to separate ΔR₂(0) from ρ. Script `research/phase2/p2_6_2_scaling.jl`; output `research/results/phase2/p2_6_2_scaling/`.

| scan | h (µm) | ΔR₂(0) (ms⁻¹) | ρ (µm/ms) | D | bias at τ = 0.1 ms | τ_conv 1 % (ms) | τ_conv 0.1 % (ms) |
|---|---|---|---|---|---|---|---|
| V2_h0.01 | 0.01 | 0.1 | 0.001 | 3 | −4.6e-04 | ≥ 0.1 | ≥ 0.1 |
| V2_h0.05 | 0.05 | 0.1 | 0.005 | 3 | −1.08e-03 | ≥ 0.1 | 0.086 |
| V1_h0.1 | 0.1 | 0.1 | 0.01 | 3 | −1.13e-03 | ≥ 0.1 | 0.086 |
| V2_h0.2 | 0.2 | 0.1 | 0.02 | 3 | −1.56e-03 | ≥ 0.1 | 0.057 |
| V3_D1 | 0.1 | 0.1 | 0.01 | 1 | −1.48e-03 | ≥ 0.1 | 0.043 |
| V3_D6 | 0.1 | 0.1 | 0.01 | 6 | −1.04e-03 | ≥ 0.1 | 0.089 |
| V2_rate0.5 | 0.1 | 0.5 | 0.05 | 3 | −6.80e-03 | ≥ 0.1 | 6.3e-03 |
| V2_rate1 | 0.1 | 1 | 0.1 | 3 | −1.34e-02 | 0.063 | 3.0e-03 |
| V2_rate2 | 0.1 | 2 | 0.2 | 3 | −2.60e-02 | 0.025 | 1.5e-03 |
| V2_h0.01_rate10 | 0.01 | 10 | 0.1 | 3 | −3.11e-02 | 0.013 | 5.6e-04 |

Every bias is negative: large τ gives less attenuation.

- **V2 (τ_conv ∝ h²/D): FAIL.** At 1 % nothing is testable (no scan at ΔR₂(0) = 0.1 reaches 1 %). At 0.1 %, τ_conv·D/h² = 103, 26, 4.3 for h = 0.05, 0.1, 0.2 (h = 0.01: lower bound): τ_conv *falls* slightly as h grows instead of rising as h². A thin layer (h = 0.01 µm, step length up to 25× h) shows no τ effect beyond 5e-4.
- **V3 (τ_conv ∝ 1/D): FAIL.** At 0.1 %, τ_conv = 0.043, 0.086, 0.089 ms for D = 1, 3, 6: τ_conv rises with D instead of falling.
- **Decision point (τ_conv vs 1/ΔR₂(0)).** The bias at τ = 0.1 ms is proportional to ΔR₂(0) (−1.1e-3, −6.8e-3, −1.34e-2, −2.6e-2 for ΔR₂(0) = 0.1, 0.5, 1, 2; figure, right panel). ΔR₂(0)·τ_conv is 0.063 and 0.051 at 1 % (ΔR₂(0) = 1, 2; within 11 %) and 0.0032, 0.0030, 0.0031 at 0.1 % (ΔR₂(0) = 0.5, 1, 2; within 4 %). At ΔR₂(0) = 0.1 the 0.1 % value (0.0086) is set by the small h- and D-dependent part of the bias (≈ 1e-3 at 0.1 ms), not by the rate.
- **ΔR₂(0), not ρ.** At the same ρ = 0.1 µm/ms, h = 0.01/ΔR₂(0) = 10 has a larger bias (−3.1e-2) than h = 0.1/ΔR₂(0) = 1 (−1.3e-2). ΔR₂(0)·τ_conv there is 0.13 (1 %) and 0.0056 (0.1 %), about 2× the h = 0.1 values: ΔR₂(0)·τ is the leading control, with a factor-2 dependence on h left over.
- Empirical summary: the per-step layer exponent ΔR₂(0)·τ controls the timestep error, roughly ΔR₂(0)·τ ≲ 0.05 for 1 % and ≲ 0.003 for 0.1 % (h = 0.01–0.1 µm, D = 3). The spec's form τ ≤ c_h·h²/D (Eq. 11) is not supported in this geometry. The form of the simulator constraint (P2.6.5) is left to the user.

![scaling](results/phase2/p2_6_2_scaling/scaling.png)

## P2.6.5: layer timestep constraint

Decision (user, 2026-10-06) after P2.6.2/P2.6.3: the constraint follows the measured control, ΔR₂(0)·τ, not the spec's c_h·h²/D.

- `src/timesteps.jl`: new option `layer` in `TimeStep`, giving τ ≤ `layer` / ΔR₂(0)_max, where ΔR₂(0)_max = max ρ/h over all layer sides with ρ > 0 (`Internal.max_layer_rate`; 0 without a layer, so the option is Inf and existing timesteps are unchanged). Override with `Simulation(...; timestep=(layer=...,))`.
- Default c = **0.005**: ΔR₂(0)·τ ≈ 0.05 gives 1 % error in R (P2.6.2), divided by 10. At ΔR₂(0)·τ = 0.005 the measured error is 1.3e-3 (ΔR₂(0) = 0.5, τ = 1e-2) and 1.3e-3 (ΔR₂(0) = 1, τ ≈ 4.6e-3 grid point), i.e. about 0.1 %.
- When it binds: for w = 2 µm, D = 3 the tortuosity default is 0.04 ms, so the layer constraint binds for ΔR₂(0) > 0.125 ms⁻¹. At the Phase 2 reference ΔR₂(0) = 0.1 it does not bind.
- Tests: `test/test_layer.jl` "layer timestep constraint" (max rate with side-specific values, ρ = 0 ignored, binding value, unchanged without a layer, default, verbose message).

## P2.6.8, P2.6.10–P2.6.12: cost and consequences

Script `research/phase2/p2_6_8_cost.jl`; output `research/results/phase2/p2_6_8_cost/`. Runtimes are single runs (1 seed, 10000 spins) on the same machine; B5 was measured on an earlier day, so ratios carry a load uncertainty of order 10–20 %.

**Runtime per spin-step with the layer** (h = 0.1, ΔR₂(0) = 0.1) vs B5 (no layer):

| τ (ms) | with layer (µs) | B5 (µs) | ratio |
|---|---|---|---|
| 0.1 | 0.0319 | 0.0265 | 1.21 |
| 1e-2 | 0.0236 | 0.0169 | 1.40 |
| 1e-3 | 0.0226 | 0.0163 | 1.39 |
| 1e-4 | 0.0195 | 0.0155 | 1.26 |

**Runtime vs h** (τ = 1e-2, ΔR₂(0) = 0.1, constraint not binding): 0.0353 / 0.0246 / 0.0207 / 0.0201 / 0.0215 µs per spin-step for h = 0.01 / 0.05 / 0.1 / 0.2 / 0.5 (×2.1 … ×1.2 vs B5). A full reference run (1e5 spins × 50 ms × 10 seeds) takes ≈ 0.03–0.05 h at every h. The h = 0.01 value is the first run after the τ loop and may include load noise.

**Runtime vs ΔR₂(0) at the new default timestep** (h = 0.1, τ = min(0.005/ΔR₂(0), 0.04)):

| ΔR₂(0) (ms⁻¹) | τ (ms) | µs per spin-step | full run (h) | σ penalty at fixed budget √(0.04/τ) |
|---|---|---|---|---|
| 0.1 | 0.04 | 0.0219 | 0.01 | 1.0 |
| 0.5 | 0.01 | 0.0191 | 0.03 | 2.0 |
| 1 | 5e-3 | 0.0200 | 0.06 | 2.8 |
| 2 | 2.5e-3 | 0.0204 | 0.11 | 4.0 |
| 10 | 5e-4 | 0.0227 | 0.63 | 8.9 |

- Noise cost (P2.6.11): at a fixed compute budget the number of spins scales with τ, so σ grows as √(τ_default/τ) (last column). The cost of the constraint grows linearly with ΔR₂(0) and does not depend on h.
- Accuracy–cost (P2.6.12): `cost.png` (toy error vs MCMR runtime per spin per ms, V1 and ΔR₂(0) = 1, 2). Recommended setting: the default (c = 0.005, ≈ 0.1 % in R); `timestep=(layer=0.05,)` for ≈ 1 % at a tenth of the cost.

![cost](results/phase2/p2_6_8_cost/cost.png)

**Consequences for earlier results.** All Phase 2 characterisation runs (P2.3.4(b) sweep, P2.3.9 V4) used τ = 1e-2 ms at ΔR₂(0) = 0.1, i.e. ΔR₂(0)·τ = 1e-3, below the default constraint (0.005). Measured bias of R at τ = 1e-2 in the toy scans at ΔR₂(0) = 0.1: −2.6e-4 to −6.4e-4 (h = 0.01–0.2, D = 1–6). For V4 this corresponds to a relative shift in h* of the same order (≈ 2.6e-4 at h = 0.1, i.e. 2.5e-5 µm against σ_jack = 5e-5 µm), within about 1 σ, consistent with the V4 timestep check (≤ 1 SEM). → **No reruns needed.** The linear-profile plan may keep τ = 1e-2 at ΔR₂(0) = 0.1 (its shape dependence is P2.6.4).
