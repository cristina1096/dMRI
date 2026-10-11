# P2.5.6 prerequisite — idealised spin echo (instantaneous pulses, MRIBuilder SpinEcho as in Rep2) vs free decay at
# TE = 50 ms, for the conventional model and the exponential layer at T2_REF (w = 2, h = 0.2, D = 3, ρ = 1/70,
# R2_bulk = 1/2000), seeds 1–10, 50000 spins. No off-resonance: refocusing should not change |M⊥|, so S_SE(TE) and
# S_free(TE) agree up to sampling. Pass: |S_SE − S_free| ≤ 3 combined SEM for both models.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_5_t2_spin_echo_check.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "sensitivity.jl"))
using MRIBuilder

const OUT = phase2_outdir("p2_5_t2_spin_echo_check")
const TE = 50.0
const NSPINS = 50_000
const C = T2_REF

function se_signal(geometry)
    seq = SpinEcho(TE=TE, TR=1000.0, scanner=MRIBuilder.Siemens_Prisma, excitation=(phase=90.0,))
    sim = Simulation(seq; geometry=geometry, diffusivity=C.D, R2=1 / 2000, verbose=false, timestep=1e-2)
    vals = Float64[]
    for seed in 1:10
        Random.seed!(seed)
        snap = Snapshot(NSPINS, sim, 500)                 # longitudinal start; the excitation pulse tips into M⊥
        push!(vals, per_spin(readout(snap, sim)))        # one sequence, one readout at TE (as Rep2: per_spin)
    end
    return vals
end

rows = NamedTuple[]
for (label, geom) in (("conventional", baseline_walls(theta=theta_for(C.rho; D=C.D), w=C.w)),
                      ("exponential", rho_walls(rho=C.rho, h=C.h, shape="exponential", w=C.w)))
    se = se_signal(geom)
    fd = run_walls(geom; D=C.D, R2_bulk=1 / 2000, timestep=1e-2, times=[0.0, TE], nspins=NSPINS)
    a, b = mean(se), fd.mean[end]
    σ = hypot(std(se) / sqrt(length(se)), fd.sem[end])
    z = (a - b) / σ
    push!(rows, (model=label, S_spin_echo=a, S_free=b, sigma=σ, z=z, pass=abs(z) <= 3))
    @printf("%-12s S_SE(50) = %.6f, S_free(50) = %.6f, z = %+.2f → %s\n", label, a, b, z, abs(z) <= 3 ? "PASS" : "FAIL")
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rows" => rows, "pass" => all(r.pass for r in rows)))
println("saved to $OUT")
