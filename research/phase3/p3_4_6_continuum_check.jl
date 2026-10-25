# P3.4.6 validation: Monte Carlo vs the radial continuum, inside a single cylinder R = 1, D = 3, ρ = 0.02 (fixed ρ),
# step layer (inner side only) h = 0.05, 0.1, 0.2, 0.5, 1.0 and the conventional model θ = ρ√(π/D). (No h < 0.05:
# thin layers force τ ∝ h at fixed ρ and are covered by the continuum, user request 2026-10-25.)
# All models share one timestep per pass, τ_a = min(1e-2, 0.005/ΔR2(0) of the thinnest layer), then τ_b = τ_a/4, so the
# trajectories are identical across models (asserted) and ΔS = S_layer − S_conv is paired per seed (lesson of Review C1:
# at fixed ρ the model difference is as small as the timestep error, so it is checked at two τ).
# Seeds 1–10, NSPINS (default 20000) spins per seed in the box [−R, R]³ (≈ π/4 inside), readouts 0:5:50.
# Criteria (tracker P3.4.6): at τ_b, max over readouts |S_MC − S_cont|/SEM ≤ 3 for every model, and
# |ΔS_MC − ΔS_cont|/SEM_paired ≤ 3 for every layer (3, not 2: 11 correlated readouts).
# Run: julia --project=research/baseline -t 8 research/phase3/p3_4_6_continuum_check.jl

include(joinpath(@__DIR__, "harness.jl"))
include(joinpath(@__DIR__, "continuum_radial.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_4_6_continuum_check")
const NSPINS = parse(Int, get(ENV, "NSPINS", "20000"))
const RHO = 0.02
const HS = [0.05, 0.1, 0.2, 0.5, 1.0]
const θ = RHO * sqrt(π / D_REF)
const τa = min(1e-2, 0.005 / rate_inside(RHO, minimum(HS), R_REF))       # = 1e-2 (ΔR2(0) = 0.41 /ms at h = 0.05)
const TAUS = [τa, τa / 4]
bb = mr.BoundingBox([-R_REF, -R_REF, -R_REF], [R_REF, R_REF, R_REF])
models = vcat([("conventional", Cylinders(radius=R_REF, surface_relaxation=θ), nothing)],
              [("step_h$(h)", Cylinders(radius=R_REF, layer_rho_inner=RHO, layer_h_inner=h), h) for h in HS])
cont = Dict(name => (isnothing(h) ? radial_signal(R=R_REF, D=D_REF, times=TIMES_REF, surface=RHO) :
                                    radial_signal(R=R_REF, D=D_REF, times=TIMES_REF, layer=("step", RHO, h)))
            for (name, _, h) in models)

rows = NamedTuple[]
runs = Dict{Tuple{String, Float64}, Any}()
for τ in TAUS, (name, g, h) in models
    runs[(name, τ)] = run_cylinders(g; region=:inside, bounding_box=bb, timestep=τ, nspins=NSPINS)
    @assert runs[(name, τ)].tau == τ
end
identical = all(all(mr.position.(a) == mr.position.(b) for (a, b) in zip(runs[("conventional", τ)].finals, runs[(name, τ)].finals))
                for τ in TAUS for (name, _, _) in models)
@printf("identical trajectories across models (per τ): %s\n", identical)
allpass = identical
for τ in TAUS, (name, _, h) in models
    r = runs[(name, τ)]
    zS = (r.mean .- cont[name]) ./ max.(r.sem, 1e-15)
    if isnothing(h)
        dmean = dsem = dcont = zeros(length(TIMES_REF)); zD = zeros(length(TIMES_REF))
    else
        d = r.S .- runs[("conventional", τ)].S
        dmean = vec(mean(d, dims=1)); dsem = vec(std(d, dims=1)) ./ sqrt(size(d, 1))
        dcont = cont[name] .- cont["conventional"]
        zD = (dmean .- dcont) ./ max.(dsem, 1e-15)
    end
    pass = maximum(abs.(zS[2:end])) <= 3 && maximum(abs.(zD[2:end])) <= 3
    τ == TAUS[end] && (global allpass &= pass)
    push!(rows, (model=name, h=isnothing(h) ? 0.0 : h, tau=τ, S50_mc=r.mean[end], S50_sem=r.sem[end], S50_cont=cont[name][end],
        maxabs_zS=maximum(abs.(zS[2:end])), dS50_mc=dmean[end], dS50_sem=dsem[end], dS50_cont=dcont[end],
        maxabs_zdS=maximum(abs.(zD[2:end])), pass=pass))
    @printf("%-14s τ = %.2e: S(50) MC %.6f ± %.1e vs %.6f (max|z| %.2f); ΔS(50) MC %+.2e ± %.1e vs %+.2e (max|z| %.2f)\n",
        name, τ, r.mean[end], r.sem[end], cont[name][end], maximum(abs.(zS[2:end])), dmean[end], dsem[end], dcont[end], maximum(abs.(zD[2:end])))
end
write_csv(joinpath(OUT, "comparison.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "pass" => allpass, "identical_trajectories" => identical,
    "taus" => TAUS, "nspins" => NSPINS, "rho" => RHO, "theta" => θ))
fig = Figure(size=(700, 420))
ax = Axis(fig[1, 1], xlabel="t (ms)", ylabel="ΔS = S_layer − S_conv", title="inside R = 1, ρ = 0.02: continuum (lines), MC τ_b (points)")
for (name, _, h) in models[2:end]
    lines!(ax, TIMES_REF, cont[name] .- cont["conventional"], label="h = $h")
    d = runs[(name, TAUS[end])].S .- runs[("conventional", TAUS[end])].S
    errorbars!(ax, TIMES_REF, vec(mean(d, dims=1)), vec(std(d, dims=1)) ./ sqrt(size(d, 1)))
    scatter!(ax, TIMES_REF, vec(mean(d, dims=1)), markersize=6)
end
axislegend(ax, position=:lb, labelsize=9)
save(joinpath(OUT, "p3_4_6_continuum_check.png"), fig)
println(allpass ? "PASS" : "FAIL")
println("saved to $OUT")
