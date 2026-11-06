# P3.3.2 accuracy (Group 1 review): the linear profile's quadrature against its independent closed form, rather than
# against itself at a tighter tolerance. Inside, g = 2(1 − (R − r)/h), so ∫g ds = ∫(2 − 2R/h) ds + (2/h)∫√(as² + bs + q) ds
# over the shell parts, with ∫√(as² + bs + q) ds = (2as + b)√Q/(4a) + (4aq − b²)/(8a^(3/2))·ln|2√a√Q + 2as + b|.
# 20000 random pieces (endpoints uniform in [−1.6, 1.6]², R = 1) for h = 0.3 and h = R (pieces through the axis), at the
# default tolerance. Reported: max |quadrature − closed form|.
# Run: julia --project=research/baseline research/phase3/p3_3b_linear_closed_form.jl

include(joinpath(@__DIR__, "harness.jl"))
using StaticArrays
const OUT = phase3_outdir("p3_3_quadrature")
Layers = mr.Geometries.Internal.Layers
P(x, y) = SVector{2, Float64}(x, y)
function sqrtquad(a, b, q, s0, s1)
    F(s) = (Q = max(a * s^2 + b * s + q, 0.0); (2a * s + b) * sqrt(Q) / (4a) + (4a * q - b^2) / (8a^1.5) * log(abs(2sqrt(a) * sqrt(Q) + 2a * s + b)))
    return F(s1) - F(s0)
end
function closed(p0, p1, R, h)
    Δ = p1 - p0; a = Δ[1]^2 + Δ[2]^2; b = 2(p0[1] * Δ[1] + p0[2] * Δ[2]); q = p0[1]^2 + p0[2]^2
    (h0, h1) = Layers.disk_interval(a, b, q, R^2)
    (l0, l1) = R - h > 0 ? Layers.disk_interval(a, b, q, (R - h)^2) : (0.0, 0.0)
    parts = l1 > l0 ? [(h0, l0), (l1, h1)] : [(h0, h1)]
    return sum((x1 - x0) * (2 - 2R / h) + 2 / h * sqrtquad(a, b, q, x0, x1) for (x0, x1) in parts if x1 > x0; init=0.0)
end
Random.seed!(7)
rows = NamedTuple[]
for h in (0.3, 1.0)
    side = Layers.LayerSide(0.02, h, :linear)
    err = 0.0
    for _ in 1:20_000
        p0 = P((3.2 .* rand(2) .- 1.6)...); p1 = P((3.2 .* rand(2) .- 1.6)...)
        err = max(err, abs(Layers.shell_profile_integral(side, 1.0, true, P(0.0, 0.0), p0, p1) - closed(p0, p1, 1.0, h)))
    end
    push!(rows, (h=h, max_abs_error_vs_closed_form=err))
    @printf("linear inside, h = %.1f: max |quadrature − closed form| = %.1e\n", h, err)
end
write_csv(joinpath(OUT, "linear_closed_form.csv"), rows)
