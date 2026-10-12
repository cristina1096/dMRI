# P2.5 review check (C1) — is the fixed-ρ profile effect ΔS = S_layer − S_conv biased by the timestep?
#
# The independent review compared P2.5 with an exact 1D continuum solution of the same model on walls and found the
# Monte Carlo biased positive for thin layers at the adopted τ (= min(1e-2, 0.005/ΔR2(0))). Here: the same paired
# comparison at τ = 1e-2, 2.5e-3 and 6.25e-4 ms (factor 4 steps), w = 2, D = 3, ρ = 0.01, R2_bulk = 0, readouts
# 0:10:100 ms, seeds 1–10, 50000 spins. Configurations: step and linear at h = 0.05, and step, linear, exponential at
# h = 0.2 (the P2.5 reference point). Per τ the conventional run (θ = ρ√(π/D)) is shared by all layers.
# Output per configuration and τ: ΔS(t) mean and paired SEM; ΔS at 50 and 90 ms. Characterisation (no pass/fail):
# the result is the τ trend, compared with the continuum value in the record.
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_5_tau_check.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "sensitivity.jl"))

const OUT = phase2_outdir("p2_5_tau_check")
const NSPINS = parse(Int, get(ENV, "NSPINS", "50000"))
const TIMES = collect(0.0:10.0:100.0)
const TAUS = [1e-2, 2.5e-3, 6.25e-4]
const C = (w=2.0, D=3.0, rho=0.01)
const CONFIGS = [("step", 0.05), ("linear", 0.05), ("step", 0.2), ("linear", 0.2), ("exponential", 0.2)]

rows = NamedTuple[]
for τ in TAUS
    base = run_walls(baseline_walls(theta=theta_for(C.rho; D=C.D), w=C.w); D=C.D, timestep=τ, times=TIMES, nspins=NSPINS)
    for (shape, h) in CONFIGS
        lay = run_walls(rho_walls(rho=C.rho, h=h, shape=shape, w=C.w); D=C.D, timestep=τ, times=TIMES, nspins=NSPINS)
        d = paired_difference(lay.S, base.S)
        i50, i90 = findfirst(==(50.0), TIMES), findfirst(==(90.0), TIMES)
        push!(rows, (shape=shape, h_um=h, tau_ms=τ, dS_50=d.mean[i50], dS_50_sem=d.sem[i50], dS_90=d.mean[i90], dS_90_sem=d.sem[i90],
            runtime_s=lay.runtime))
        @printf("τ = %.2e  %-11s h = %.2f: ΔS(50) = %+.2e ± %.0e, ΔS(90) = %+.2e ± %.0e\n", τ, shape, h, d.mean[i50], d.sem[i50], d.mean[i90], d.sem[i90])
        flush(stdout)
    end
end
write_csv(joinpath(OUT, "tau_check.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "taus" => TAUS, "times" => TIMES, "nspins" => NSPINS,
    "config" => Dict(pairs(C)), "rows" => rows))
println("saved to $OUT")
