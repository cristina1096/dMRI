# Annulus continuum (P3.7b, support model outside cylinders): a Krogh-type cell around one cylinder of radius R in a
# square array of spacing s, R ≤ r ≤ R_out with π R_out² = s² (equal area), reflecting at R_out. An APPROXIMATION of
# the outside compartment (no neighbour layers, no square cell corners); validated against Monte Carlo in P3.7b Task 3.
# Independent of the package's layer code.
#
#     ∂m/∂t = D (1/r) ∂_r(r ∂_r m) − [R2_bulk + ΔR2(r)] m,   R ≤ r ≤ R_out,   uniform start m = 1,   ∂_r m = 0 at R_out,
#     layer:        ∂_r m = 0 at r = R,  ΔR2(r) = (ρ/h)·g((r − R)/h)/κ_out on 0 ≤ r − R ≤ L,  κ_out = 1 + c̄/R
#                   (L = h for step/linear, the cutoff for the exponential; c̄ the centroid with that support, D3.1/D3.5),
#     conventional: D ∂_r m = ρ m at r = R (Robin; the sink is on the inner face of the fluid region).
# If R + L > R_out the layer is cut at R_out; by default (`renormalise=true`, user decision 2026-11-09 after the P3.7b
# validation) the kept part is scaled by 1/kept so that ∫ΔR2 dA = 2πRρ, the mean-exponent identity E[−ln M] = ρ(S/V)t
# that holds exactly in the true cell. Without it, a 1.5 % loss (exponential at s = 2.5) swamped the profile effect.
# `annulus_layer_kept` reports the kept fraction; `renormalise=false` gives the cut, unscaled layer. Finite volumes as continuum_radial.jl: N rings of width dr = (R_out − R)/N, exact ring integrals of g, the
# Robin sink on the inner face as 2πRρ/(1 + ρ dr/2D) (second order); y = √A m gives a symmetric tridiagonal M and
# S(t) = Σ_k (u_kᵀ√A)² e^(−λ_k t) / Σ A_i.
# Roundoff floor: eigenvalues carry an absolute error ≈ ε‖M‖ ≈ 1e-8 /ms at N = 2000 (‖M‖ ∝ D/dr²), so S(t) is exact
# only to ≈ 1e-8·t (e.g. 1e-6 at 100 ms); differences between models share most of it.
# Exact references (Bessel): with u(r) = J₀(kr)Y₁(kR_out) − Y₀(kr)J₁(kR_out) (∂_r u = 0 at R_out), the Robin rate is
# the lowest D k² with D u′(R) = ρ u(R); the ceiling (perfectly absorbing cylinder) the lowest D k² with u(R) = 0.

using LinearAlgebra
import SpecialFunctions: besselj0, besselj1, bessely0, bessely1

"Equal-area (Krogh) outer radius of a square-array cell of spacing s: π R_out² = s²."
krogh_radius(s) = s / sqrt(π)

"Centroid h∫u g du (µm) of an outer profile with support L: step h/2, linear h/3, exponential h(1 − (c+1)e^(−c))/(1 − e^(−c)), c = L/h."
function outer_centroid(shape, h, L)
    shape == "step" && return h / 2
    shape == "linear" && return h / 3
    shape == "exponential" || error("outer_centroid: profile $shape is not supported")
    c = L / h
    return h * (1 - (c + 1) * exp(-c)) / -expm1(-c)
end

"""
    outer_ring_integral(shape, h, R, L, lo, hi)

∫_lo^hi g((r − R)/h) r dr exactly, clipped to the support [R, R + L]. Antiderivatives in r: step r²/2;
linear 2[(1 + R/h) r²/2 − r³/(3h)]; exponential −h e^(−(r − R)/h)(r + h)/(1 − e^(−c)), c = L/h.
"""
function outer_ring_integral(shape, h, R, L, lo, hi)
    a, b = max(lo, R), min(hi, R + L)
    b <= a && return 0.0
    if shape == "step"
        return (b^2 - a^2) / 2
    elseif shape == "linear"
        F(r) = 2 * ((1 + R / h) * r^2 / 2 - r^3 / (3h))
        return F(b) - F(a)
    elseif shape == "exponential"
        G(r) = -h * exp(-(r - R) / h) * (r + h)
        return (G(b) - G(a)) / -expm1(-L / h)
    end
    error("outer_ring_integral: profile $shape is not supported")
end

"Fraction of ∫ΔR2 dA of the outer layer that lies inside the annulus (1 when R + L ≤ R_out)."
annulus_layer_kept(shape, h, L, R, R_out) =
    outer_ring_integral(shape, h, R, L, R, min(R + L, R_out)) / outer_ring_integral(shape, h, R, L, R, R + L)

"Symmetric operator M (1/ms) and ring areas A; `layer = (shape, ρ, h, L)`, `surface = ρ` for the conventional model."
function annulus_operator(; R, R_out, D, layer=nothing, surface=0.0, N=2000, renormalise=true)
    R_out > R || error("annulus_operator: R_out = $R_out must exceed R = $R")
    dr = (R_out - R) / N
    redge = R .+ collect(0:N) .* dr
    A = [π * (redge[i + 1]^2 - redge[i]^2) for i in 1:N]
    V = zeros(N)
    if !isnothing(layer)
        (shape, ρ, h, L) = layer
        κ = 1 + outer_centroid(shape, h, L) / R
        scale = renormalise ? 1 / annulus_layer_kept(shape, h, L, R, R_out) : 1.0
        for i in 1:N
            V[i] = scale * ρ / (h * κ) * 2π * outer_ring_integral(shape, h, R, L, redge[i], redge[i + 1]) / A[i]
        end
    end
    c = [2π * D * redge[i + 1] / dr for i in 1:N - 1]
    K = zeros(N)
    for i in 1:N - 1
        K[i] += c[i]
        K[i + 1] += c[i]
    end
    K[1] += 2π * R * surface / (1 + surface * dr / (2D))        # Robin sink on the inner face r = R
    M = SymTridiagonal((K .+ A .* V) ./ A, -c ./ sqrt.(A[1:N - 1] .* A[2:N]))
    return (M, A)
end

"Decay rates (1/ms, ascending) of the discretised annulus operator."
annulus_rates(; R2_bulk=0.0, kw...) = eigvals(annulus_operator(; kw...)[1]) .+ R2_bulk

"S(t) closure from one eigendecomposition of the annulus operator."
function annulus_signal_fun(; R, R_out, D, layer=nothing, surface=0.0, bulk=0.0, N=2000, renormalise=true)
    (M, A) = annulus_operator(; R=R, R_out=R_out, D=D, layer=layer, surface=surface, N=N, renormalise=renormalise)
    F = eigen(M)
    c2 = (F.vectors' * sqrt.(A)) .^ 2 ./ sum(A)
    return t -> sum(c2 .* exp.(-F.values .* t)) * exp(-bulk * t)
end

"""
    annulus_signal(; R, R_out, D, times, layer=nothing, surface=0.0, R2_bulk=0.0, N=2000)

S(t)/S(0) in the annulus for a uniform start. `layer = (shape, ρ, h, L)` (ρ per unit surface area, D3.1), `surface = ρ`
for the conventional model. `renormalise=true` (default) scales a layer cut at R_out by 1/kept so that ∫ΔR2 dA = 2πRρ.
"""
annulus_signal(; R, R_out, D, times, layer=nothing, surface=0.0, R2_bulk=0.0, N=2000, renormalise=true) =
    annulus_signal_fun(R=R, R_out=R_out, D=D, layer=layer, surface=surface, bulk=R2_bulk, N=N, renormalise=renormalise).(times)

"First root of f on (k0, kmax) by a scan in steps dk then bisection (f continuous, no poles)."
function first_root(f, k0, dk, kmax)
    a, fa = k0, f(k0)
    while a < kmax
        b = a + dk
        fb = f(b)
        if sign(fb) != sign(fa)
            for _ in 1:200
                m = (a + b) / 2
                fm = f(m)
                if sign(fm) == sign(fa)
                    a, fa = m, fm
                else
                    b = m
                end
            end
            return (a + b) / 2
        end
        a, fa = b, fb
    end
    error("first_root: no sign change below $kmax")
end
annulus_u(k, r, R_out) = besselj0(k * r) * bessely1(k * R_out) - bessely0(k * r) * besselj1(k * R_out)
annulus_du(k, r, R_out) = -k * (besselj1(k * r) * bessely1(k * R_out) - bessely1(k * r) * besselj1(k * R_out))

"Exact lowest rate (1/ms) of the annulus with Robin ρ at R (test reference): D k², D u′(R) = ρ u(R)."
function exact_rate_annulus_robin(R, R_out, D, ρ)
    k = first_root(k -> D * annulus_du(k, R, R_out) - ρ * annulus_u(k, R, R_out), 1e-6, 1e-3 * π / (R_out - R), 1e3 / (R_out - R))
    return D * k^2
end

"Conventional ceiling outside (perfectly absorbing cylinder, reflecting R_out): lowest D k² with u(R) = 0."
function annulus_ceiling(R, R_out, D)
    k = first_root(k -> annulus_u(k, R, R_out), 1e-6, 1e-3 * π / (R_out - R), 1e3 / (R_out - R))
    return D * k^2
end
