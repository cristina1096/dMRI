# Baseline test log

## Set up gitHub

- Git repo: https://github.com/cristina1096/dMRI
- Create a new branch: `codex/near-surface-r2`

## Set up Julia and run the unmodified code
- Julia: 1.12.7
- Code commit tested: d9c5d0e
- Command: julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["evolve"])'
- Result: PASS — 113 of 113 tests passed in 16.9 seconds.
- Dependency issues: None observed during this test.

- Command: julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["radio_frequency"])'
- Result: PASS — 69 of 69 tests passed in 8.4 seconds.
- Dependency issues: None observed during this test.

## Learn the code path

### Read `src/MCMRSimulator.jl`.
- `include(...)` loads the package's source files in order.
- `import .Module: name` brings a name from a submodule into `MCMRSimulator`.
- `export` lets users access selected names directly after `using MCMRSimulator`.
- `Walls`, `Cylinders`, and `Annuli` are user-facing geometry constructors.
- `BoundingBox` defines a region, including where initial spins may be placed.
- Plotting names ending in `!` add to an existing plot; plotting requires Makie.

### Read `src/simulations.jl`
- `struct Simulation` stores the settings used by a simulation; defining it does not move spins.
- The constructor beginning at line 96 prepares the user's inputs. The internal constructor beginning at line 70 creates the object with `new(...)`.
- The constructor's `R2` becomes the global baseline rate in `GlobalProperties`, stored as `sim.properties`. Geometry may also contribute local R2.
- `fix(geometry; ...)` converts user geometry into the internal form stored as `sim.geometry`.
- The timestep setting becomes a `TimeStep` controller stored as `sim.timestep`.
- `flatten` records whether the sequence was supplied directly or in a vector; it affects the shape of returned results.
- A later call to `evolve(spins, sim, time)` uses these stored settings to run the spins.

### Read the main flow in `src/evolve.jl`
- `evolve(spins, simulation, new_time)` returns a Snapshot at the requested target time. `new_time` is not one timestep; `TR` can select a particular sequence repetition.
- `run_readout!` processes the sequence as many short parts, limited by the timestep controller and sequence events.
- For each part, `process_sequence_step!` calls `draw_step!` and then applies instantaneous events.
- `draw_step!` proposes a random movement, detects any physical wall collision, and calls `relax!` to update magnetisation over the travelled segment. At a collision it then handles surface attenuation, permeability, and the remaining movement.
- To revisit later: the proposed near-surface layer boundary at distance λ is different from a physical wall collision.

### Read the magnetisation update paths in `src/relax.jl`
- The spin-level `relax!` calculates local MRI properties at the movement segment's midpoint, then updates one `SpinOrientation` per MRI sequence.
- With no active RF (`NoPulsePart`), it updates phase and applies R1/R2 relaxation over the segment.
- With finite RF (`PulsePart`), it divides the pulse into waveform pieces. `apply_pulse!` applies half the relaxation, an RF rotation, then half the relaxation.
- Instantaneous RF events are handled separately in `src/evolve.jl`.
- Relevant to my project: the current local R2 is treated as constant over a movement segment. A distance-dependent near-surface R2 will require evaluating how it changes along that segment.

### Read `src/timesteps.jl` and the `parts(...)` function in `src/sequence_parts.jl`.
- `TimeStep(...)` calculates a maximum outer timestep from limits for tortuosity, sticking/transfer, permeability, collision-based surface relaxation, and dwell time.
- `parts(...)` first divides the sequence at readouts, instantaneous events, gradient changes, and finite-RF start/end times.
- It then divides intervals that exceed the timestep limit; gradient strength can make this limit shorter.
- Each resulting part has an actual duration (`t1 - t0`). `draw_step!` uses that duration for one proposed random movement.
- Subdivisions inside finite-RF `relax!` refine magnetisation updates; they do not add random movements.
- Relevant to my project: resolving a thin near-surface layer may require a new limit on the outer timestep.

### Read `walls.jl`, `rounds.jl`, and the wall-related functions in `fixed_obstruction_groups.jl`.
- An internal `Wall` uses one local coordinate, but represents an infinite plane in 3-D. `rotation` selects its perpendicular direction, and `Shift` places it at the requested position.
- `detect_intersection` in `walls.jl` checks whether a proposed straight movement crosses the physical wall. Its result gives the collision as a fraction of that movement.
- The geometry-level code converts 3-D positions to the wall's local coordinate and returns the earliest physical collision.
- `previous_intersection` identifies a recently hit surface. The wall detector uses it to avoid immediately counting the same collision again.
- The current detector finds physical wall crossings; it does not find crossings of my proposed near-surface boundary at distance λ.
- I am deferring `rounds.jl` and repeating cylinder geometry until the cylinder stage.

### Simulation flow summary

1. `Simulation(...)` prepares and stores the sequences, geometry, MRI properties, and timestep controller.
2. `SequenceParts.parts(...)` divides the sequence into short time intervals.
3. For each interval, `draw_step!` proposes spin movement and handles physical collisions; it calls `relax!` along travelled segments.
4. `relax!` updates each spin's magnetisation and phase using its position and the active sequence part.
5. At a readout, spin magnetisations are combined into the signal; `evolve(...)` instead returns a Snapshot at the requested time.