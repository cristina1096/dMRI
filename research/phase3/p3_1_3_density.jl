# P3.1.3 C-B3: density across the cylinder (inside and outside) and the fraction within h of the surface; P3.4.5: the
# same with a step layer on (positions must equal the no-layer run). Default timestep, seeds 1–10, 1e5 spins, t = 50 ms.
#
# Inside (single cylinder R = 1): r²/R² is uniform on [0, 1] for a uniform density → 100-bin χ² of r²/R²,
#   z_χ = (χ² − dof)/√(2 dof), dof = 99, |z_χ| ≤ 3; occupancy within h: p = (R² − (R − h)²)/R², z = (n − Np)/√(Np(1 − p)), |z| ≤ 2.
# Outside (square array s = 2.5): for R ≤ r ≤ s/2 whole annuli lie in the cell, so P(r ∈ [r_k, r_k+1] | r ≤ s/2) =
#   (r_k+1² − r_k²)/((s/2)² − R²) → 50-bin χ² on spins with r ≤ s/2; occupancy p = π((R + h)² − R²)/(s² − πR²).
# Layer on: ρ = 0.02 um/ms, h = 0.2 um both sides.
# Run: julia --project=research/baseline -t 8 research/phase3/p3_1_3_density.jl

include(joinpath(@__DIR__, "harness.jl"))
using CairoMakie
const OUT = phase3_outdir("p3_1_3_density")
const HS_IN = [0.05, 0.1, 0.2, 0.5, 1.0]
const HS_OUT = [0.05, 0.1, 0.2]

zchi(counts, p) = (N = sum(counts); e = N .* p; χ2 = sum((counts .- e) .^ 2 ./ e); dof = length(p) - 1; (χ2 / dof, (χ2 - dof) / sqrt(2dof)))
binom_z(n, N, p) = (n - N * p) / sqrt(N * p * (1 - p))

function radii_at_end(r; s)
    vcat([[radial_distance(sp.position; sx=s, sy=s) for (k, sp) in enumerate(f.spins) if m[k]] for (f, m) in zip(r.finals, r.masks)]...)
end

bb = mr.BoundingBox([-R_REF, -R_REF, -R_REF], [R_REF, R_REF, R_REF])
rows = NamedTuple[]
summary = Dict{String, Any}()
allpass = true
radii_store = Dict{String, Vector{Float64}}()
for (label, layer_kw) in (("no_layer", (;)), ("layer_on", (layer_rho=0.02, layer_h=0.2)))
    ins = run_cylinders(Cylinders(radius=R_REF; layer_kw...); region=:inside, bounding_box=bb, timestep=nothing, times=[0.0, 50.0])
    out = run_cylinders(Cylinders(radius=R_REF, repeats=[S_REF, S_REF]; layer_kw...); region=:outside, sx=S_REF, sy=S_REF,
        timestep=nothing, times=[0.0, 50.0])
    rin = radii_at_end(ins; s=Inf)
    rout = radii_at_end(out; s=S_REF)
    radii_store[label * "_in"] = rin
    radii_store[label * "_out"] = rout
    # inside histogram of r²/R²
    u = clamp.((rin ./ R_REF) .^ 2, 0, prevfloat(1.0))
    cin = [count(x -> (k - 1) / 100 <= x < k / 100, u) for k in 1:100]
    (chi_in, z_in) = zchi(cin, fill(1 / 100, 100))
    # outside histogram of r on [R, s/2]
    edges = sqrt.(range(R_REF^2, (S_REF / 2)^2, length=51))
    ro = filter(r -> r <= S_REF / 2, rout)
    cout = [count(r -> edges[k] <= r < edges[k + 1], ro) for k in 1:50]
    (chi_out, z_out) = zchi(cout, fill(1 / 50, 50))           # equal-area bins
    occ_ok = true
    for h in HS_IN
        p = shell_area_inside(R_REF, h)
        n = count(r -> r >= R_REF - h, rin)
        # h = R: the shell is the whole disk (p = 1), so the fraction must be exactly 1 and z is undefined
        z = p == 1 ? 0.0 : binom_z(n, length(rin), p)
        occ_ok &= p == 1 ? n == length(rin) : abs(z) <= 2
        push!(rows, (case=label, region="inside", h=h, p=p, fraction=count(r -> r >= R_REF - h, rin) / length(rin), z=z))
    end
    for h in HS_OUT
        p = shell_area_outside(R_REF, h, S_REF, S_REF)
        n = count(r -> r <= R_REF + h, rout)
        z = binom_z(n, length(rout), p)
        occ_ok &= abs(z) <= 2
        push!(rows, (case=label, region="outside", h=h, p=p, fraction=n / length(rout), z=z))
    end
    pass = abs(z_in) <= 3 && abs(z_out) <= 3 && occ_ok
    summary[label] = Dict("inside_chi2_per_dof" => chi_in, "inside_z_chi2" => z_in, "outside_chi2_per_dof" => chi_out,
        "outside_z_chi2" => z_out, "occupancy_all_within_2sigma" => occ_ok, "n_inside" => length(rin), "n_outside" => length(rout),
        "tau" => ins.tau, "pass" => pass)
    global allpass &= pass
    @printf("%-8s: inside χ²/dof %.3f (z %+.2f), outside χ²/dof %.3f (z %+.2f), occupancy ok %s → %s\n",
        label, chi_in, z_in, chi_out, z_out, occ_ok, pass ? "PASS" : "FAIL")
    summary[label * "_runs"] = (ins, out)
end
# P3.4.5: positions with the layer on equal the no-layer run
(ni, no) = summary["no_layer_runs"]
(li, lo) = summary["layer_on_runs"]
same = all(mr.position.(a) == mr.position.(b) for (a, b) in zip(vcat(ni.finals, no.finals), vcat(li.finals, lo.finals)))
@printf("P3.4.5 positions identical with the layer on: %s\n", same)
delete!(summary, "no_layer_runs"); delete!(summary, "layer_on_runs")
summary["layer_on_positions_identical"] = same
write_csv(joinpath(OUT, "occupancy.csv"), rows)
write_json(joinpath(OUT, "summary.json"), Dict("provenance" => provenance(), "results" => summary, "pass" => allpass && same))

fig = Figure(size=(1000, 380))
ax1 = Axis(fig[1, 1], xlabel="r²/R²", ylabel="density / uniform", title="inside, t = 50 ms (no layer)")
u = (radii_store["no_layer_in"] ./ R_REF) .^ 2
hist!(ax1, u, bins=50, normalization=:pdf)
hlines!(ax1, [1.0], color=:black)
ax2 = Axis(fig[1, 2], xlabel="r (um)", ylabel="density per unit area / uniform", title="outside, r ≤ s/2 (no layer)")
ro = filter(r -> r <= S_REF / 2, radii_store["no_layer_out"])
edges = range(R_REF, S_REF / 2, length=31)
c = [count(r -> edges[k] <= r < edges[k + 1], ro) for k in 1:30]
area = [edges[k + 1]^2 - edges[k]^2 for k in 1:30]
scatter!(ax2, (edges[1:end-1] .+ edges[2:end]) ./ 2, (c ./ area) ./ (length(ro) / ((S_REF / 2)^2 - R_REF^2)))
hlines!(ax2, [1.0], color=:black)
save(joinpath(OUT, "p3_1_3_density.png"), fig)
println("saved to $OUT")
