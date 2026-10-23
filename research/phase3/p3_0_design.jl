# Phase 3, design decision D3.1 — what ρ means on a curved surface, and its consequence for the profile.
#
# Per unit length of a cylinder of radius R, a layer ΔR2(d) inside (d = R − r) or outside (d = r − R) contains
#     ∫ ΔR2 dA = 2πR ∫ ΔR2(d) (1 ∓ d/R) dd        (− inside, + outside).
# Choosing ρ per unit SURFACE AREA (D3.1 option b) means ρ = ∫ ΔR2(d)(1 ∓ d/R) dd, so the profile is
#     ΔR2(d) = (ρ/h) g(d/h) / κ,   κ = ∫ g(u)(1 ∓ u h/R) du = 1 ∓ c̄/R,
# with c̄ = h∫u g(u) du the profile centroid (step h/2, linear h/3, exponential truncated). Then the volume mean of
# ΔR2 over a compartment of volume V and surface A is ρ·A/V, the same as the conventional relaxivity ρ times S/V.
# This script (1) checks κ = 1 ∓ c̄/R by quadrature, (2) plots κ vs h/R, (3) plots ΔR2 inside a cylinder for both
# definitions of ρ.
#
# Run: julia --project=research/baseline research/phase3/p3_0_design.jl

include(joinpath(@__DIR__, "..", "baseline", "common.jl"))
using CairoMakie

const OUT = mkpath(joinpath(REPO_ROOT, "research", "results", "phase3", "p3_0_design"))

g(shape, u, c) = shape == "step" ? (0 <= u <= 1 ? 1.0 : 0.0) :
                 shape == "linear" ? (0 <= u <= 1 ? 2 * (1 - u) : 0.0) :
                 (0 <= u <= c ? exp(-u) / (1 - exp(-c)) : 0.0)
"Centroid c̄ (µm) of the profile with support c·h (exponential) — step h/2, linear h/3."
cbar(shape, h, c) = shape == "step" ? h / 2 : shape == "linear" ? h / 3 : h * (1 - (c + 1) * exp(-c)) / (1 - exp(-c))
"κ by midpoint quadrature; sgn = −1 inside, +1 outside; support c (in units of h) for the exponential."
function κ_quad(shape, h, R, sgn; c=R / h, n=200_000)
    umax = shape == "exponential" ? c : 1.0
    du = umax / n
    return sum(g(shape, (k - 0.5) * du, c) * (1 + sgn * (k - 0.5) * du * h / R) for k in 1:n) * du
end

R = 1.0
rows = NamedTuple[]
maxdev = 0.0
for shape in ("step", "linear", "exponential"), hR in (0.05, 0.1, 0.25, 0.5, 1.0), (side, sgn) in (("inside", -1), ("outside", +1))
    h = hR * R
    c = side == "inside" ? R / h : 2.0 / h          # exponential support: inside to the centre (R); outside an example cutoff of 2 µm
    shape == "step" && side == "inside" && hR > 1 && continue
    kq = κ_quad(shape, h, R, sgn; c=c)
    kc = 1 + sgn * cbar(shape, h, c) / R
    global maxdev = max(maxdev, abs(kq - kc))
    push!(rows, (shape=shape, side=side, h_over_R=hR, kappa_quadrature=kq, kappa_closed=kc))
    @printf("%-11s %-7s h/R=%.2f: κ = %.6f (quadrature), 1 ∓ c̄/R = %.6f\n", shape, side, hR, kq, kc)
end
@printf("max |κ_quadrature − (1 ∓ c̄/R)| = %.1e\n", maxdev)
write_csv(joinpath(OUT, "kappa.csv"), rows)

fig = Figure(size=(1250, 380))
ax = Axis(fig[1, 1], xlabel="h/R", ylabel="κ = ∫ g(u)(1 ∓ u h/R) du", title="curvature normalisation (inside: −, outside: +)")
colors = Dict("step" => :steelblue, "linear" => :darkorange, "exponential" => :seagreen)
for shape in ("step", "linear", "exponential"), (side, ls) in (("inside", :solid), ("outside", :dash))
    q = [r for r in rows if r.shape == shape && r.side == side]
    scatterlines!(ax, [r.h_over_R for r in q], [r.kappa_closed for r in q], color=colors[shape], linestyle=ls, label="$shape, $side")
end
hlines!(ax, [1.0], color=:black, linewidth=0.5)
Legend(fig[1, 2], ax, labelsize=8, rowgap=0)
ax2 = Axis(fig[1, 3], xlabel="r (µm), cylinder R = 1", ylabel="ΔR₂(r) / (ρ/h)", title="linear layer inside, h = 0.5 µm: two definitions of ρ")
h = 0.5
rs = range(0.0, R, length=400)
lines!(ax2, rs, [g("linear", (R - r) / h, Inf) for r in rs], color=:gray, label="(a) ρ along the normal")
κ = 1 - cbar("linear", h, Inf) / R
lines!(ax2, rs, [g("linear", (R - r) / h, Inf) / κ for r in rs], color=:darkorange, label="(b) ρ per surface area (÷κ = $(round(κ, digits=3)))")
axislegend(ax2, position=:lt, labelsize=8)
save(joinpath(OUT, "p3_0_kappa.png"), fig)
println("saved to $OUT")
