# P2.5.9 spot check — Monte Carlo vs continuum for overlapping layers: linear h/w = 0.75 and 1.0, exponential h/w = 0.5
# and 1.0, at ρ = 0.1 and 1.0 (w = 2, D = 3, R2_bulk = 0). τ = (0.005 / max_x ΔR2_summed(x)) / 4, capped at 1e-2, so
# the summed (overlapping) maximum sets the step (review M1). Readouts: 10 points up to the time where the conventional
# continuum signal is 0.05. Seeds 1–10, 50000 spins. Pass: |z| ≤ 3 at every readout (z = (MC − continuum)/SEM).
#
# Run: julia --project=research/baseline -t 8 research/phase2/p2_5_overlap_mc_check.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum.jl"))
include(joinpath(@__DIR__, "overlap.jl"))

const OUT = phase2_outdir("p2_5_overlap_mc_check")
const W, D = 2.0, 3.0
rows = NamedTuple[]
for ρ in (0.1, 1.0), (shape, hw) in (("linear", 0.75), ("linear", 1.0), ("exponential", 0.5), ("exponential", 1.0))
    h = hw * W
    rmax = maximum(summed_rate_closed(shape, ρ, h, W, x) for x in range(0.0, W, length=2001))
    τ = min(0.005 / rmax / 4, 1e-2)
    tend = log(20) / continuum_rates(w=W, D=D, surface=ρ, N=1500)[1]
    T = collect(range(0.0, tend, length=10))
    cont = continuum_signal(w=W, D=D, times=T, layer=(shape, ρ, h))
    r = run_walls(rho_walls(rho=ρ, h=h, shape=shape, w=W); D=D, timestep=τ, times=T, nspins=50_000)
    z = (r.mean[2:end] .- cont[2:end]) ./ max.(r.sem[2:end], 1e-9)   # floor: the uniform case has SEM ≈ 1e-16 (exact)
    push!(rows, (rho=ρ, shape=shape, h_over_w=hw, tau_ms=τ, max_abs_z=maximum(abs.(z)), max_abs_diff=maximum(abs.(r.mean .- cont)), pass=maximum(abs.(z)) <= 3))
    @printf("ρ=%.1f %-11s h/w=%.2f τ=%.1e: max |z| = %.2f, max |MC − continuum| = %.1e → %s\n", ρ, shape, hw, τ, rows[end].max_abs_z, rows[end].max_abs_diff, rows[end].pass ? "PASS" : "FAIL")
end
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "rows" => rows, "pass" => all(r.pass for r in rows)))
println("saved to $OUT")
