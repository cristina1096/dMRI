# P3.4.3: mean-exponent identity on cylinders. With ρ per unit surface area (D3.1) and a uniform stationary density,
#     E[−ln M⊥(t)] = t·⟨ΔR2⟩_compartment = ρ·(S/V)·t
# for any h and τ, including overlapping outer layers (the summed rate integrates to ρ per unit area of every surface
# whose layer reaches the compartment; validation keeps outer layers out of other cylinders). Tests the κ normalisation,
# the shell fractions, the image loop and the overlap sum end to end.
# Configurations (ρ = 0.02, step): single inside R = 1, h = 0.2 and R = 0.5, h = 0.2; square array s = 2.5 inside and
# outside, h = 0.2 (separated) and h_outer = 0.4 (overlap, gap 0.5 < 0.8); rectangular sx = 2.4, sy = 3.0, h_outer = 0.3
# (overlap in x only: gaps 0.4 and 1.0 vs 2h = 0.6). τ = 1e-2 and 4e-2 ms, t = 50 ms, seeds 1–10, 100000 spins.
# z = (mean −ln M − ρ(S/V)t)/SEM, |z| ≤ 3.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_3_identity.jl

include(joinpath(@__DIR__, "harness.jl"))
const OUT = phase3_outdir("p3_4_3_identity")
const RHO, T = 0.02, 50.0
configs = [
    ("single_R1", Cylinders(radius=1.0, layer_rho=RHO, layer_h=0.2), 1.0, Inf, Inf, :inside),
    ("single_R0.5", Cylinders(radius=0.5, layer_rho=RHO, layer_h=0.2), 0.5, Inf, Inf, :inside),
    ("square_sep_in", Cylinders(radius=1.0, repeats=[2.5, 2.5], layer_rho=RHO, layer_h=0.2), 1.0, 2.5, 2.5, :inside),
    ("square_sep_out", Cylinders(radius=1.0, repeats=[2.5, 2.5], layer_rho=RHO, layer_h=0.2), 1.0, 2.5, 2.5, :outside),
    ("square_overlap_out", Cylinders(radius=1.0, repeats=[2.5, 2.5], layer_rho=RHO, layer_h_inner=0.2, layer_h_outer=0.4), 1.0, 2.5, 2.5, :outside),
    ("rect_overlap_x_out", Cylinders(radius=1.0, repeats=[2.4, 3.0], layer_rho=RHO, layer_h_inner=0.2, layer_h_outer=0.3), 1.0, 2.4, 3.0, :outside),
]
rows = NamedTuple[]
for (name, g, R, sx, sy, region) in configs, τ in (1e-2, 4e-2)
    bb = isfinite(sx) ? 500 : mr.BoundingBox([-R, -R, -R], [R, R, R])
    r = run_cylinders(g; region=region, R=R, sx=sx, sy=sy, bounding_box=bb, timestep=τ, times=[0.0, T])
    sv = region == :inside ? sv_inside(R) : sv_outside(R, sx, sy)
    e = expected_neglog_layer(RHO, sv, T)
    z = (r.neglog_mean[end] - e) / r.neglog_sem[end]
    push!(rows, (config=name, tau=τ, expected=e, measured=r.neglog_mean[end], sem=r.neglog_sem[end], z=z, S50=r.mean[end], pass=abs(z) <= 3))
    @printf("%-20s τ = %.0e: E[−ln M] %.6f vs ρ(S/V)t %.6f, z = %+.2f\n", name, τ, r.neglog_mean[end], e, z)
end
write_csv(joinpath(OUT, "identity.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "n_pass" => count(r -> r.pass, rows), "n" => length(rows),
    "pass" => all(r -> r.pass, rows)))
println("saved to $OUT")
