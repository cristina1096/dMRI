# P1.2.6 — Single-spin trace through the unmodified simulator.
#
# Goal: confirm that the code-reading in research/phase1_code_log.md (§2–§5) matches
# what the simulator actually does, by stepping one spin through every sequence part
# and predicting each step's magnetisation change independently.
#
# Setup (variables used):
#   geometry          : one non-repeating cylinder, radius R = 1 µm, along z
#   surface_relaxation: θ = 0.1 (per sqrt(ms)) on the cylinder surface
#   diffusivity       : D = 3 µm²/ms
#   R2 (global)       : 1/80 kHz  (T2 = 80 ms);  R1 = 1/1000 kHz
#   timestep          : τ = 0.01 ms fixed (TimeStep(τ, Inf))
#   sequences         : (a) spin echo, instant 90°/180° pulses, TE = 5 ms
#                       (b) spin echo, finite hard pulses (0.5 ms / 1 ms), TE = 5 ms
#   spins             : 100 spins inside the cylinder; spin #1 traced
#   seed              : Random.seed!(20260929)
#
# For every no-RF step the trace predicts
#     M⊥(after) = M⊥(before) · exp(-R2·dt) · exp(-θ·sqrt(τ))^n_hits
# where n_hits comes from replaying the same proposed displacement through
# `draw_step!(…, test_new_pos)` (relax.jl L60–63 and evolve.jl L555–561).
# Instant pulses are checked against a direct rotation (evolve.jl L646–660).
#
# Run:  julia --project=research/baseline -t 1 research/baseline/p1_2_6_single_spin_trace.jl

include(joinpath(@__DIR__, "common.jl"))
using MRIBuilder
import StaticArrays: SVector
import Rotations

const OUT = outdir("p1_2_6_single_spin_trace")
const SEED = 20260929
const R = 1.0
const θ = 0.1
const D = 3.0
const R2g = 1 / 80
const R1g = 1 / 1000
const τ = 0.01
const TE = 5.0

Ev = mr.Evolve
SP = mr.SequenceParts

"Proposed Gaussian step exactly as draw_step! draws it (evolve.jl L499–504), without consuming the spin's RNG."
function proposed_position(spin, timestep)
    local p
    saved = copy(Random.TaskLocalRNG())
    copy!(Random.TaskLocalRNG(), spin.rng)
    p = spin.position + randn(SVector{3, Float64}) .* sqrt(2 * D * timestep)
    copy!(Random.TaskLocalRNG(), saved)
    return p
end

function run_trace(label, seq, io)
    tee(io, "\n==== Sequence: $label ====")
    sim = Simulation(seq; diffusivity=D, R2=R2g, R1=R1g,
        geometry=Cylinders(radius=R, surface_relaxation=θ), timestep=τ, verbose=false)
    tee(io, "sim.timestep = ", sim.timestep)
    Random.seed!(SEED)
    # 100 spins uniformly inside the cylinder (x,y); z arbitrary
    spins = []
    while length(spins) < 100
        p = (rand(3) .- 0.5) .* 2R
        if p[1]^2 + p[2]^2 < R^2
            push!(spins, Spin(position=p))
        end
    end
    snap = Snapshot([s for s in spins])
    spin = snap.spins[1]
    parts = SP.parts(sim.sequences, 0.0, sim.timestep)
    B0s = [3.0]
    tee(io, "number of sequence parts = ", length(parts))

    rows = NamedTuple[]
    n_checked = 0
    max_err = 0.0
    n_hit_steps = 0
    total_hits = 0
    first_hits_logged = 0
    for (k, part) in enumerate(parts)
        p = part.parts[1]
        before_T = spin.orientations[1].transverse
        before_L = spin.orientations[1].longitudinal
        before_pos = spin.position

        # replay the proposed step on a copy to find collision points
        nhits = 0
        path = nothing
        if part.duration > 0
            prop = proposed_position(spin, part.duration)
            copy_spin = deepcopy(spin)
            path = Ev.draw_step!(copy_spin, sim, part, B0s, prop)
            nhits = length(path) - 2
        end

        # the real update (all spins, as the simulator does it)
        Ev.draw_step!(snap.spins, sim, part, B0s)
        if !isnothing(path)
            @assert spin.position ≈ path[end] "replayed path must end where the real step ended"
        end
        after_T = spin.orientations[1].transverse

        kind = typeof(p).name.name
        pred = NaN
        if part.duration > 0 && p isa SP.NoPulsePart
            pred = before_T * exp(-R2g * part.duration) * exp(-θ * sqrt(τ))^nhits
            err = abs(after_T - pred) / max(abs(pred), 1e-300)
            max_err = max(max_err, err)
            n_checked += 1
        end
        if nhits > 0
            n_hit_steps += 1
            total_hits += nhits
        end

        # instantaneous events (pulses) — apply the same way apply_instants! does, check rotation
        for inst in part.instants.instants[1]
            if inst isa SP.PulseEvent
                o = spin.orientations[1]
                rot = Rotations.RotationVec(deg2rad(inst.flip_angle) * cosd(inst.phase), deg2rad(inst.flip_angle) * sind(inst.phase), 0.0)
                expected = mr.SpinOrientation(rot * mr.orientation(o))
                Ev.apply_instants!(snap.spins, 1, inst, nothing)
                o2 = spin.orientations[1]
                tee(io, @sprintf("  t=%.4f ms: instant pulse flip=%.0f° phase=%.0f° → L %.6f→%.6f (exp %.6f), T %.6f→%.6f (exp %.6f)",
                    sum(q.duration for q in parts[1:k]), inst.flip_angle, inst.phase,
                    o.longitudinal, o2.longitudinal, expected.longitudinal, before_T, o2.transverse, expected.transverse))
            end
        end

        t_end = sum(q.duration for q in parts[1:k])
        push!(rows, (step=k, t_end_ms=t_end, dt_ms=part.duration, part_type=string(kind),
            x=spin.position[1], y=spin.position[2], z=spin.position[3], r=hypot(spin.position[1], spin.position[2]),
            n_hits=nhits, longitudinal=spin.orientations[1].longitudinal,
            transverse=spin.orientations[1].transverse, phase=spin.orientations[1].phase,
            predicted_transverse=pred))

        if nhits > 0 && first_hits_logged < 3
            first_hits_logged += 1
            tee(io, @sprintf("  step %d (t=%.3f ms, %s): %d collision(s). path:", k, t_end, kind, nhits))
            for q in path
                tee(io, @sprintf("      (%.5f, %.5f, %.5f)  r=%.6f", q[1], q[2], q[3], hypot(q[1], q[2])))
            end
            tee(io, @sprintf("      M⊥ %.8f → %.8f ; predicted %.8f  (factor exp(-R2 dt)=%.8f, per-hit exp(-θ√τ)=%.8f)",
                before_T, after_T, pred, exp(-R2g * part.duration), exp(-θ * sqrt(τ))))
        end
    end
    inside = all(r -> r.r <= R + 1e-9, rows)
    tee(io, @sprintf("steps with collision: %d, total collisions: %d", n_hit_steps, total_hits))
    tee(io, @sprintf("no-RF steps checked: %d, max relative error vs prediction: %.3e", n_checked, max_err))
    tee(io, "spin stayed inside cylinder for all steps: ", inside)
    tee(io, @sprintf("final: L=%.6f  T=%.6f  phase=%.3f°", spin.orientations[1].longitudinal, spin.orientations[1].transverse, spin.orientations[1].phase))
    fname = replace(lowercase(label), r"[^a-z0-9]+" => "_")
    write_csv(joinpath(OUT, "trace_$(fname).csv"), rows)
    return (label=label, nparts=length(parts), n_hit_steps=n_hit_steps, total_hits=total_hits,
            n_checked=n_checked, max_rel_err=max_err, stayed_inside=inside)
end

open(joinpath(OUT, "trace.log"), "w") do io
    prov = provenance()
    tee(io, "provenance: ", prov)
    seq_instant = SpinEcho(TE=TE, scanner=MRIBuilder.Siemens_Prisma, excitation=(phase=90.,))
    seq_finite = SpinEcho(TE=TE, scanner=MRIBuilder.Siemens_Prisma,
        excitation=(shape=:hard, duration=0.5, phase=90.), refocus=(shape=:hard, duration=1.0))
    res = [run_trace("instant pulses", seq_instant, io), run_trace("finite hard pulses", seq_finite, io)]
    write_json(joinpath(OUT, "summary.json"), Dict("provenance" => prov, "results" => res,
        "params" => Dict("R" => R, "theta" => θ, "D" => D, "R2" => R2g, "R1" => R1g, "tau" => τ, "TE" => TE, "seed" => SEED)))
end
