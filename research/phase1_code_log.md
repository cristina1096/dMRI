# Phase 1 · Task 1.2 — Code Study Log

Code map of the **unmodified** MCMRSimulator source on branch `cc/near-surface-r2`. Every file, function and line number below refers to commit `d9c5d0e` ("Iintial commit"); `src/` has not been modified.

- Package version in `Project.toml`: **1.1.0**. The paper's figure notebooks pin **v1.0.0** with MRIBuilder v0.3.0 and Julia 1.11.7 (see the §0 note).
- Environment used here: Julia 1.12.7, MRIBuilder 0.4.1 (git tag `v0.4.1`), `research/baseline/Project.toml` + `Manifest.toml`.

Status key: ☑ done · ◐ partial · ☐ not started

| ID | Task | Status | Output |
|---|---|---|---|
| P1.2.1 | Code structure | ☑ | §1 |
| P1.2.2 | Simulation workflow, τ_max | ☑ | §2 |
| P1.2.3 | Geometry representation | ☑ | §3 |
| P1.2.4 | Spin trajectories | ☑ | §4 |
| P1.2.5 | Magnetisation update (with / without RF), θ_relax lines | ☑ | §5 |
| P1.2.6 | Single-spin trace | ☑ | §6, `research/results/baseline/p1_2_6_single_spin_trace/` |
| P1.2.7 | Insertion points for R₂(d) | ☑ | §7 |

---

## §0 Version note (affects P1.1.1)

| Item | Paper notebooks (`mcmr_paper_figures`, commit `c7b86f9`) | This work |
|---|---|---|
| Julia | 1.11.7 | 1.12.7 |
| MCMRSimulator | v1.0.0 (`repo-rev = "v1.0.0"`) | 1.1.0 (this repo, `src/` unmodified) |
| MRIBuilder | v0.3.0 | v0.4.1 (root `Project.toml` compat requires 0.4.1) |

MCMRSimulator 1.1.0 is not the paper's version. The tracker asks for v1.0, but this repo's code (including local SWC / overlapping-sphere additions in `src/geometries/user/obstructions/obstructions.jl:47-56`) is what Phase 2 will modify, so the baseline is taken on this code. The R1–R4 results in `baseline_record.md` show whether the differences matter.

---

## §1 Code structure (P1.2.1)

**Entry point:** [src/MCMRSimulator.jl](../src/MCMRSimulator.jl). Lines 15–29 `include` the submodules in dependency order:

| Order | File | Module | Role |
|---|---|---|---|
| 1 | `constants.jl` | `Constants` | `gyromagnetic_ratio` |
| 2 | `scanners.jl` | `Scanners` | B0, max gradient, slew rate |
| 3 | `methods.jl` | `Methods` | `get_time`, `get_rotation` |
| 4 | `properties.jl` | `Properties` | `GlobalProperties(R1, R2, off_resonance)`, `stick_probability` |
| 5 | `geometries/geometries.jl` | `Geometries` (`.User`, `.Internal`) | user geometry → internal `FixedGeometry`; collision detection; susceptibility |
| 6 | `spins.jl` | `Spins` | `Spin`, `SpinOrientation`, `Snapshot`, per-spin RNG |
| 7 | `timesteps.jl` | `TimeSteps` | `TimeStep` (τ_max computation) |
| 8 | `pulseq/pulseq.jl` | `Pulseq` | reads `.seq` files |
| 9 | `sequence_parts.jl` | `SequenceParts` | cuts sequences into timesteps (`parts`) |
| 10 | `simulations.jl` | `Simulations` | `Simulation` object |
| 11 | `relax.jl` | `Relax` | Bloch update for one segment (`relax!`) |
| 12 | `subsets.jl` | `Subsets` | `Subset` (inside/outside/bound signal selection) |
| 13 | `evolve.jl` | `Evolve` | main loop: `readout`, `evolve`, `draw_step!` |
| 14 | `plot.jl`, `cli/cli.jl` | `Plot`, `CLI` | plotting (Makie extension), `mcmr` CLI |

**Geometry split:**
- `src/geometries/user/`: user-facing `ObstructionGroup` types (`Walls`, `Cylinders`, `Annuli`, `Spheres`, `Mesh`, `BendyCylinder`), `Field` definitions, and `fix` (user → internal).
- `src/geometries/internal/`: `FixedObstructionGroup`, the obstruction primitives (`walls.jl`, `rounds.jl`, `triangles.jl`, `shifts.jl`), hit grids, intersections, reflections, property look-ups and susceptibility fields.

**How MRIBuilder sequences get in.** MRIBuilder is a weak dependency. [ext/MRIBuilderToMCMRSimulator/](../ext/MRIBuilderToMCMRSimulator/) defines `gradient_waveform`, `get_pulses`, `get_instants` and `readout_times` for `MRIBuilder.Sequence`. [sequence_parts.jl:205-231](../src/sequence_parts.jl#L205-L231) `SequenceWaveform(sequence)` calls these to build gradient waveforms, finite RF pulses (`(t0, t1, Vector{ConstantPulse})`), instant events (`PulseEvent`, `GradientEvent`) and ADC sample times. Pulseq files take the same route via `src/pulseq/`.

**Dependency direction:** `Simulation` holds the geometry, properties and `TimeStep` → `Evolve` drives the loop → it calls `SequenceParts.parts` for timing, `Geometries.Internal.detect_intersection` for collisions, and `Relax.relax!` for magnetisation.

---

## §2 Simulation workflow (P1.2.2)

### 2.1 Setup: `Simulation(...)`, [simulations.jl:96-138](../src/simulations.jl#L96-L138)

| Line | What happens |
|---|---|
| 114 | `fix_susceptibility(geometry)` → susceptibility sources (myelin / iron) |
| 115 | `fix(geometry; permeability, density=surface_density, dwell_time, relaxation=surface_relaxation)` → `FixedGeometry` (global collision defaults filled in where not set per obstruction) |
| 117–120 | `inside_geometry` = groups with non-zero volume R1/R2/off-resonance; `prepare_isinside!` |
| 122 | `GlobalProperties(R1, R2, off_resonance)` |
| 134 | `timestep isa Number ? TimeStep(timestep, Inf) : TimeStep(; diffusivity, geometry, timestep...)` |

Line 126 builds the same `TimeStep` into a tuple (trailing comma) and never uses it. It is dead code with no effect.

### 2.2 τ_max from the constraint list (Supp. S1): [timesteps.jl:29-77](../src/timesteps.jl#L29-L77)

`TimeStep(; diffusivity, geometry, size_scale, tortuosity=3e-2, gradient=1e-4, permeability=0.5, surface_relaxation=0.01, transfer=0.01, dwell_time=0.1)`. The static bound is the **minimum** of (lines 39–45):

| # | Term | Code | Physical form |
|---|---|---|---|
| 1 | tortuosity | `tortuosity * size_scale^2 / D` | 0.03 · l²/D; `size_scale` = min radius / wall spacing / repeat ([size_scales.jl:25-54](../src/geometries/user/size_scales.jl#L25-L54)) |
| 2 | sticking / MT | `max_timestep_sticking` ([properties.jl:137-145](../src/geometries/internal/properties.jl#L137-L145)) | keeps `log(1 − p_stick)` ≤ scaling |
| 3 | permeability | `0.5 · κ_max⁻²` ([timesteps.jl:94](../src/timesteps.jl#L94)) | |
| 4 | surface relaxation | `0.01 · θ_max⁻²` ([timesteps.jl:101](../src/timesteps.jl#L101)) | per-collision loss θ√τ ≤ 0.1 |
| 5 | dwell time | `0.1 · min dwell_time` | |

`idx = argmin(options)` (line 46) names the **binding constraint** in the verbose log (lines 47–75). The gradient term is applied per interval: `(ts::TimeStep)(G) = min(max_timestep, (gradient/D / G²)^(1/3))` ([timesteps.jl:84-87](../src/timesteps.jl#L84-L87)). **Finite RF does not limit τ_max.** RF sub-division happens inside `relax!` (§5.3).

### 2.3 Cutting the sequence into steps: `parts(...)`, [sequence_parts.jl:246-360](../src/sequence_parts.jl#L246-L360)

1. Control times (lines 254–305) are readout times, instant events, gradient slew changes (`compress_timeseries`) and finite-RF start/end times.
2. Lines 322–332 split each interval `[t0, t1]` into `ceil((t1 − t0)/timestep(max|G|))` equal pieces.
3. `build_blocks` (lines 380–395) → one `MultSequencePart` per piece. Its fields are `duration`, one `SequencePart` per sequence (`EmptyPart` / `ConstantPart` / `LinearPart`, or `PulsePart` if RF is on, from `get_block` lines 397–431), and the instants at the **end** of the piece.

### 2.4 Main loop: [evolve.jl](../src/evolve.jl)

```
readout(nspins, sim, times)            evolve.jl:355-387  (splits into runs of ≤1e7 magnetisations)
 └ run_readout!(snapshot, sim, acc)    evolve.jl:396-405
    └ for part in parts(sequences, t0, sim.timestep)       ← sequence_parts.jl:246
        process_sequence_step!(spins, sim, part, B0s, acc)  evolve.jl:462-465
          ├ draw_step!(spins, …)       evolve.jl:475-479   Threads.@threads over spins
          │   └ draw_step!(spin, …)    evolve.jl:482-595   (§4, §5)
          └ apply_instants!(spins, part.instants, acc)      evolve.jl:610-660
              ├ PulseEvent   → rotation (L646-660)
              ├ GradientEvent→ phase += x·q (L630-636)
              └ IndexedReadout → readout!(accumulator) (L638-644)   ← signal summed here
```

`evolve(snapshot, sim, t)` ([evolve.jl:415-460](../src/evolve.jl#L415-L460)) is `readout_internal(…, return_snapshot=true)`. The output is the per-spin state at `t`.

**Readout ordering:** at a control time, the `IndexedReadout` is `pushfirst!`-ed before other instants (sequence_parts.jl:353). A readout at the same time as a pulse therefore sees the state **before** the pulse. A readout at t = 0 before a t = 0 excitation returns M⊥ = 0.

---

## §3 Geometry representation (P1.2.3)

### 3.1 User → internal
- User types: `ObstructionGroup` with `FieldValue`s ([obstruction_types.jl:42-63](../src/geometries/user/obstructions/obstruction_types.jl#L42-L63), [obstructions.jl](../src/geometries/user/obstructions/obstructions.jl)).
- `fix_type` per type ([fix.jl:49-125](../src/geometries/user/fix.jl#L49-L125)):

| User type | Internal primitive | Dimensionality | Notes |
|---|---|---|---|
| `Walls` | `Internal.Wall()` ([walls.jl:12](../src/geometries/internal/obstructions/walls.jl#L12)) | 1 (coordinate along the normal) | `has_inside = false` |
| `Cylinders` | `Round{2}(radius)` ([rounds.jl:12-17](../src/geometries/internal/obstructions/rounds.jl#L12-L17)) | 2 | inside if r² < R² |
| `Spheres` | `Round{3}` / `OverlappingSphere` | 3 | |
| `Annuli` | **two** cylinder groups: inner (index i) + outer (index i+1) ([fix.jl:93-113](../src/geometries/user/fix.jl#L93-L113)) | 2 | inner-volume R1/R2 stored as (inner − outer) so the sum is correct |
| `Mesh` / `BendyCylinder` | `IndexTriangle` ([triangles.jl](../src/geometries/internal/obstructions/triangles.jl)) | 3 | split into connected components |

- Position offsets wrap each primitive in `Shift` ([shifts.jl:17](../src/geometries/internal/obstructions/shifts.jl#L17)). Repeat positions are folded into [−r/2, r/2) ([fix.jl:141-147](../src/geometries/user/fix.jl#L141-L147)).
- `FixedObstructionGroup` ([fixed_obstruction_groups.jl:46-89](../src/geometries/internal/fixed_obstruction_groups.jl#L46-L89)) stores `repeats`, `rotation` (3×N: local N-D ↔ global 3D), `hit_grid`, `volume::(R1,R2,off_resonance)`, `surface::(R1,R2,off_resonance,permeability,surface_density,dwell_time,surface_relaxation)` and `size_scale`. `FixedGeometry = NTuple{N, FixedObstructionGroup}` (line 99).

### 3.2 Surface vs volume parameters (paper Table 1)
- [fields.jl:163-172](../src/geometries/user/obstructions/fields.jl#L163-L172) `property_fields`:
  - `R1`, `R2` and `off_resonance` are `per_volume` **and** `per_surface`.
  - `dwell_time`, `density`, `permeability` and `relaxation` are `per_surface` only.
- [fix.jl:158-177](../src/geometries/user/fix.jl#L158-L177) `apply_properties` builds the `volume` and `surface` NamedTuples. `nothing` (not set per object) falls back to the `Simulation` keyword (`fix_array`, lines 149–156). `relaxation` is renamed `surface_relaxation` and `density` is renamed `surface_density`.
- **Look-up during simulation:**
  - Volume and stuck-surface MRI properties: `MRIProperties(full_geometry, inside_geometry, glob, position, stuck_to)` ([properties.jl:38-57](../src/geometries/internal/properties.jl#L38-L57)) = global + stuck-surface value + Σ over groups whose `isinside` contains the position.
  - Collision properties: `surface_relaxation(geometry, intersection)` etc. ([properties.jl:94-125](../src/geometries/internal/properties.jl#L94-L125)). Mesh gaps return `gap_*` constants (permeability ∞, relaxation 0).

### 3.3 Collision grid and repeating geometry
- `HitGrid` ([hit_grids.jl:71-130](../src/geometries/internal/hit_grids.jl#L71-L130)) precomputes, per voxel, which obstruction bounding boxes (and which repeat shifts) overlap it. The default resolution gives about 1 voxel per obstruction ([size_scales.jl:65-76](../src/geometries/user/size_scales.jl#L65-L76)).
- `detect_intersection_grid` ([hit_grids.jl:252-281](../src/geometries/internal/hit_grids.jl#L252-L281)) walks voxels along the ray (`ray_grid_intersections`) and stops once the best hit lies before the voxel exit.
- Repeating: `detect_intersection_repeating` ([fixed_obstruction_groups.jl:219-241](../src/geometries/internal/fixed_obstruction_groups.jl#L219-L241)) takes positions modulo `repeats` and walks the repeat cells along the ray. `isinside` folds the position into the base cell (lines 132–147).
- Global dispatch: `detect_intersection(geometries::FixedGeometry, start, dest, previous)` ([fixed_obstruction_groups.jl:244-253](../src/geometries/internal/fixed_obstruction_groups.jl#L244-L253)) = the closest hit over all groups, returned as `Intersection(distance ∈ [0,1], normal, inside, geometry_index, obstruction_index, hit_gap)`.

---

## §4 Spin trajectories (P1.2.4): `draw_step!(spin, …)`, [evolve.jl:482-595](../src/evolve.jl#L482-L595)

| Lines | Step |
|---|---|
| 499–510 | Inside `@spin_rng` (per-spin RNG, [spins.jl:204-213](../src/spins.jl#L204-L213)): **Gaussian step** `displacement = randn(SVector{3}) · sqrt(2Dτ)` (L502–504). `Reflection(ratio_displaced = |Δ|/sqrt(2Dτ))` remembers the step length for later re-scaling. |
| 511 | Up to 10⁶ bounces per step, then error (L588–590). |
| 512–528 | If stuck (bound pool): exponential dwell time, `relax!` with `new_pos = inside::Bool`, then release along `direction(reflection, …)`. |
| 530–535 | `detect_intersection(sim.geometry, current_pos, new_pos, phit)`: **segment–surface intersection**, not endpoint testing. `phit` = previously hit obstruction, which is handled specially to avoid re-detecting the same hit (e.g. walls.jl:40-46 ignores a wall within 1 nm of the start). |
| 537–544 | `use_distance = prevfloat(collision.distance)`: the spin stops **just before** the surface. `next_fraction_timestep` is the time fraction at contact. |
| 547 | `relax!(spin, collision_pos, …, fraction_timestep, next_fraction_timestep)`: magnetisation advanced **to the collision time** (§5). |
| 549–553 | No intersection → accept `new_pos`, done. |
| 555–561 | **Surface relaxation θ_relax** (§5.4). |
| 563–565 | **Permeability** (§2.10): `x = sqrt(τ)·κ`; `P(pass) = 1 − e^(−x) I0(x)`. κ = ∞ gives `P = 1` and needs no random draw (perfectly permeable: L565 `isone(permeability_prob)`). |
| 566–570 | `Reflection(collision, direction, …, passes_through)`. Reflected direction `d − 2(n·d)n/|n|²` ([reflections.jl:49-61](../src/geometries/internal/reflections.jl#L49-L61)). If permeable: direction unchanged, `inside` flipped. |
| 576–582 | Sticking with `stick_probability` ([properties.jl:60-68](../src/properties.jl#L60-L68)); otherwise the **remaining** displacement is `direction(reflection, remaining_time, D)`, whose length is `ratio·sqrt(2D(t_moved + t_left)) − distance_moved` ([reflections.jl:99-103](../src/geometries/internal/reflections.jl#L99-L103)). |
| 584–585 | Loop back: re-test the remaining segment from the collision point. |

Spin initialisation: `Snapshot(nspins, bounding_box, geometry)` ([spins.jl:340-361](../src/spins.jl#L340-L361)) places spins uniformly in the box using the **global** RNG. Each spin gets its own `FixedXoshiro` seeded from the global RNG ([spins.jl:48-53](../src/spins.jl#L48-L53)). `Random.seed!(s)` before `readout`/`Snapshot` therefore reproduces a run regardless of thread count.

---

## §5 Magnetisation update (P1.2.5): [relax.jl](../src/relax.jl)

### 5.1 Spin-level entry: `relax!(spin, new_pos, sim, parts, t1, t2, B0s)`, [relax.jl:31-49](../src/relax.jl#L31-L49)

Called from `draw_step!` at **evolve.jl:547** (each free sub-segment up to a collision or the end of the step) and **evolve.jl:515** (stuck spins).
- L33–34: `mean_pos = (spin.position + new_pos)/2`, then `props = MRIProperties(...)`. **R1, R2 and off-resonance are evaluated once at the segment midpoint and held constant along the segment.**
- L36: susceptibility off-resonance at a **random point** along the segment ([simulations.jl:206-212](../src/simulations.jl#L206-L212)).
- L41–44: precomputes `R1_att = exp(−R1·Δt)`, `R2_att = exp(−R2·Δt)` with Δt = (t2 − t1)·duration.
- L46–48: dispatches on the sequence-part type, once per sequence.

### 5.2 Case 1: **no RF** (`NoPulsePart` = `EmptyPart` / `ConstantPart` / `LinearPart`), [relax.jl:53-63](../src/relax.jl#L53-L63)

```
relax!(orient, old_pos, new_pos, part::NoPulsePart, props, duration, t1, t2, off_res, pre_comp)   L53
    phase += (off_res + props.off_resonance + grad_off_resonance(part, old,new, t1,t2)) · Δt · 360   L54-56
    relax_single_step!(orient, props, pre_comp)                                                        L57
        longitudinal = 1 − (1 − Mz)·R1_att                                                             L61
        transverse  *= R2_att                                                                          L62   ← bulk/volume R2 applied here
```
The gradient phase uses the exact path average for constant and linear gradients (L154–172).

### 5.3 Case 2: **finite RF on** (`PulsePart`), [relax.jl:71-152](../src/relax.jl#L71-L152)

```
relax!(…, pulse::PulsePart, …, ::NamedTuple)                     L71-86
    relax_time = 1/max(R1,R2);  nsplit = ceil(Δt_block / (relax_time/10))   L72-78   (extra sub-steps if relaxation is fast)
    map (t1,t2) onto the RF block grid (first/last partial blocks)      L80-85
relax!(…, pulse::PulsePart, …, split::Val)                        L88-119
    for each ConstantPulse block overlapping [t1,t2]: apply_pulse!(…)   L97-117
apply_pulse!(…, ::Val{1})                                         L121-142   (Strang splitting)
    relax_single_step!(orient, props, dt/2)       L136   (relaxation, timestep form L65-68)
    rotate by RotationVec(B1 cos φ, B1 sin φ, Δω)  L130-139
    relax_single_step!(orient, props, dt/2)       L141
apply_pulse!(…, ::Val{N})                         L144-152   (N equal sub-chunks of the block)
```
RF chunking splits **only the magnetisation update**. The spin still takes one straight step per timestep, and `props` is still evaluated at the segment midpoint.

### 5.4 Surface relaxation θ_relax: **exact lines**

[evolve.jl:555-561](../src/evolve.jl#L555-L561), inside `draw_step!` after the segment update (L547) and before the permeability decision (L563):
```julia
relaxation = surface_relaxation(simulation.geometry, collision)     # L555  (properties.jl:94-125)
if ~iszero(relaxation)
    collision_attenuation = exp(-sqrt(timestep) * relaxation)       # L557  logarithmic scheme, τ = full step duration
    for orientation in spin.orientations
        orientation.transverse *= collision_attenuation            # L559  M⊥ only; Mz and phase untouched
    end
end
```
- It is applied at **every** collision, **including** ones where the spin then passes through a permeable surface. Mesh gaps have relaxation 0.
- `timestep` here is `parts.duration`, the full step length even when the collision happens part-way through the step.
- Timestep constraint: `0.01 · θ_max⁻²` ([timesteps.jl:101](../src/timesteps.jl#L101)).
- Instant pulses do not relax; they rotate only ([evolve.jl:646-660](../src/evolve.jl#L646-L660)).

---

## §6 Single-spin trace (P1.2.6)

Script: [research/baseline/p1_2_6_single_spin_trace.jl](baseline/p1_2_6_single_spin_trace.jl). Output: `research/results/baseline/p1_2_6_single_spin_trace/` (`trace.log`, `trace_*.csv`, `summary.json`).

**Setup:**

| Variable | Value |
|---|---|
| Geometry | single cylinder, R = 1 µm, `surface_relaxation` θ = 0.1 |
| D | 3 µm²/ms |
| Global R2, R1 | 1/80, 1/1000 kHz |
| Timestep | fixed τ = 0.01 ms |
| Spins | 100 spins inside the cylinder, `Random.seed!(20260929)`; spin #1 traced |
| Sequences | (a) SpinEcho TE = 5 ms, instant 90°/180°; (b) the same with finite hard pulses (0.5 ms / 1 ms) |

**Method:**
- The script walks `SequenceParts.parts(...)` itself, calling `Evolve.draw_step!` and `Evolve.apply_instants!` as the simulator does.
- For each step it replays the same proposed Gaussian displacement through `draw_step!(…, test_new_pos)` on a copy to get the collision points.
- For every no-RF step it predicts M⊥(after) = M⊥(before) · e^(−R2·Δt) · e^(−θ√τ)^n_hits.

**Result:**

| | instant pulses | finite hard pulses |
|---|---|---|
| sequence parts (steps) | 501 | 526 |
| steps with ≥1 collision / total collisions | 111 / 112 | 120 / 121 |
| no-RF steps checked against prediction | 500 | 375 |
| max relative error | **6.2 × 10⁻¹⁶** | **5.1 × 10⁻¹⁶** |
| spin stayed inside cylinder (r ≤ R) | yes | yes |
| collision point radius | r = 1.000000 (stopped at `prevfloat`) | same |
| 90° instant pulse (phase 90°) | M⊥ 0 → 1.000000 (expected 1.000000) | n/a (finite) |
| 180° instant pulse | M⊥ 0.576301 → 0.576229, Mz −0.002497 unchanged (matches direct rotation) | n/a |

**Pass: the trace matches the reading of the code.**
- The per-collision factor e^(−θ√τ) (= 0.99005 for θ = 0.1, τ = 0.01) is applied exactly once per detected intersection.
- Bulk R2 is applied as e^(−R2·Δt) over each sub-segment.
- The replayed path (start → collision point(s) → end) ends exactly where the real step ended.

---

## §7 Insertion points for R₂(d) (P1.2.7)

These are the minimum places a distance-dependent layer term touches. No code is changed in Phase 1.

| # | Purpose | File · function (lines) | Why here |
|---|---|---|---|
| 1a | Parameter definition (`rho`, `lambda`, `shape` as **surface** fields) | [fields.jl](../src/geometries/user/obstructions/fields.jl#L163-L172) `property_fields` | Same mechanism as `relaxation`: `per_surface=true` gives per-object and per-surface values (`inner_surface` / `outer_surface` for annuli) for free |
| 1b | Carry into the internal representation | [fix.jl](../src/geometries/user/fix.jl#L168-L177) `apply_properties` (surface symbol list L169–170); [fixed_obstruction_groups.jl](../src/geometries/internal/fixed_obstruction_groups.jl#L52) NamedTuple key list in the `S` type parameter; global defaults via [simulations.jl](../src/simulations.jl#L96-L115) keywords → `fix(...)` (L115, fix.jl L25) | `A = rho/lambda` can be derived here, never user-settable |
| 1c | Accessors and gap values | [properties.jl](../src/geometries/internal/properties.jl#L94-L130) (loop over collision properties + `gap_*` constants) | |
| 1d | Construction-time validation (P2.1.2) | [fix.jl](../src/geometries/user/fix.jl#L25-L34) `fix(::ObstructionGroup, …)` | Already checks required fields here |
| 2a | Distance evaluation (point → nearest layered face) | new method per primitive next to `isinside` / `detect_intersection`: [walls.jl](../src/geometries/internal/obstructions/walls.jl#L18-L38) (d = \|x\| in the 1-D local coordinate), [rounds.jl](../src/geometries/internal/obstructions/rounds.jl#L21-L53) (d = \|\|r\| − R\|) | Local coordinates are already available |
| 2b | Map a global position → local primitive (rotation, repeat folding, grid look-up) | reuse [fixed_obstruction_groups.jl](../src/geometries/internal/fixed_obstruction_groups.jl#L128-L147) `isinside(g, pos, …)` pattern and [hit_grids.jl](../src/geometries/internal/hit_grids.jl#L183-L241) `get_objects` / `isinside(grid, …)` | Handles `Shift`, `repeats` and the grid for free |
| 2c | Crossing detection of the offset surface d = λ | analogue of `detect_intersection` with the plane shifted to ±λ ([walls.jl:24-38](../src/geometries/internal/obstructions/walls.jl#L24-L38)) / radius R ± λ ([rounds.jl:28-53](../src/geometries/internal/obstructions/rounds.jl#L28-L53)); grid walk via [hit_grids.jl:252-313](../src/geometries/internal/hit_grids.jl#L252-L313) with bounding boxes enlarged by λ (`HitGrid(...; extend=λ)`, [hit_grids.jl:94-101](../src/geometries/internal/hit_grids.jl#L94-L101), already supported) | Segment–surface intersection, not endpoint testing (P2.2.2) |
| 3a | Relaxation application (no RF) | [evolve.jl:530-547](../src/evolve.jl#L530-L547) `draw_step!`: split the free segment `current_pos → collision_pos` at layer crossings **before** `relax!` (L547); the integral ∫ΔR₂(d(t))dt goes into M⊥ like `relax_single_step!` ([relax.jl:60-63](../src/relax.jl#L60-L63)) | `relax!` currently uses one midpoint `props` per segment (relax.jl:33-34), which cannot resolve a thin layer |
| 3b | Relaxation application (finite RF) | [relax.jl:121-142](../src/relax.jl#L121-L142) `apply_pulse!` (relaxation half-steps L136/L141 with `props.R2`) | Phase 3 (V13); out of scope for walls / no-RF |
| 3c | Keep θ_relax path intact for the λ → 0 check (V4) | [evolve.jl:555-561](../src/evolve.jl#L555-L561) | Baseline comparator |
| 4 | Timestep constraint τ ≤ c_λ·λ²/D (P2.6.5) | [timesteps.jl:39-45](../src/timesteps.jl#L39-L45) `options` tuple (+ verbose message L47–75) and a `max_timestep_layer` helper next to [timesteps.jl:101](../src/timesteps.jl#L101) | Same pattern as the existing constraints |
| 5 | λ = 0 short-circuit (P2.2.8) | branch at the top of the new code in `draw_step!` / geometry type parameter | Keeps baseline runtime |

**Status (2026-10-02):** 1a, 1b, 1d, 2a, 2c, 3a and 5 are implemented for walls (step profile, per-side ρ and h) in `src/geometries/internal/layers.jl`, `src/geometries/user/fix.jl` (`wall_layers`) and `src/evolve.jl` (`apply_layer!`). 2c uses the exact linear overlap (`segment_overlap`) instead of the hit grid. 3b and 4 are not yet implemented.

---

## §8 Observations and possible issues (recorded, not changed)

1. **Per-spin RNG reconstruction drops one state word.** [spins.jl:54](../src/spins.jl#L54): `Random.Xoshiro(rng.s0, rng.s2, rng.s2, rng.s3)` uses `s2` twice and never `s1`. Streams stay deterministic and pass the statistical checks in R1–R4, but the state entropy is reduced. Worth reporting upstream; do not change in Phase 1.
2. **Dead code:** [simulations.jl:126](../src/simulations.jl#L126) builds an unused tuple.
3. **`SequenceWaveform.TR` for MRIBuilder sequences** is the last sample or gradient time, not the sequence's `TR` (probe: SpinEcho TE = 20, TR = 100 gives `TR = 20`). `override_repetition_time` only exists for Pulseq ([sequence_parts.jl:200-201](../src/sequence_parts.jl#L200-L201)). This only matters for multi-TR runs.
4. **Midpoint property evaluation:** volume R2 inside an obstruction is decided by `isinside(midpoint)` of each free sub-segment (relax.jl:33-34). This is exact for both impermeable and permeable surfaces, because every surface crossing ends a sub-segment (evolve.jl:537-547). A near-surface *layer* is not a surface, so its boundary would fall inside a sub-segment and midpoint evaluation would be wrong; hence insertion point 3a.
5. **Surface relaxation uses the full step τ** even for later collisions in the same step. This is by design: the per-step hit-count correction.
