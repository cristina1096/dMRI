# Radial continuum benchmark (P3.4.6): exact signal inside a single cylinder for the near-surface layer and for the
# conventional surface relaxation, without Monte Carlo or timestep. Independent of the package's layer code.
#
#     ∂m/∂t = D (1/r) ∂_r(r ∂_r m) − [R2_bulk + ΔR2(r)] m,   0 ≤ r ≤ R,   uniform start m = 1,
#     layer:        ∂_r m = 0 at r = R,  ΔR2(r) = (ρ/h)·g((R − r)/h)/κ,  κ = ∫g(u)(1 − u h/R) du   (D3.1; κ by quadrature),
#     conventional: D ∂_r m = −ρ m at r = R (Robin), ρ = θ√(D/π) (P3.1.6).
# Finite volumes: N rings of width dr = R/N, areas A_i = π(r_{i+1/2}² − r_{i−1/2}²), face conductances 2πD r_{i+1/2}/dr,
# the Robin sink on the outer face, 2πRρ/(1 + ρ dr/2D) (face value from the last ring, second order); ΔR2 averaged over each ring (step: exact area overlap; other profiles:
# r-weighted sub-sampling). A ṁ = −(K + A V) m with y = √A m becomes ẏ = −M y, M = A^(−1/2)(K + A V)A^(−1/2) symmetric
# tridiagonal = Σ λ_k u_k u_kᵀ, so S(t) = Σ_i A_i m_i / Σ A_i = Σ_k (u_kᵀ√A)² e^(−λ_k t) / Σ A_i.

using LinearAlgebra
import SpecialFunctions: besselj0, besselj1, besseli, besselk

"Normalised profile g(u) on its support (written independently of the package); c = support/h (exponential only)."
function radial_g(shape::AbstractString, u::Real, c::Real=Inf)
    shape == "step" && return 0 <= u <= 1 ? 1.0 : 0.0
    shape == "linear" && return 0 <= u <= 1 ? 2 * (1 - u) : 0.0
    shape == "exponential" && return 0 <= u <= c ? exp(-u) / -expm1(-c) : 0.0
    error("radial_g: profile $shape is not supported")
end

"Support of the inner layer (um): h for step and linear, R for the exponential (D3.5)."
radial_support(shape, h, R) = shape == "exponential" ? R : h

"""
    radial_kappa(shape, h, R)

κ = ∫g(u)(1 − u h/R) du over the support (inside, D3.1), from the exact ring integral of the whole support:
∫_{R−L}^R g((R − r)/h) r dr = h R κ. (A midpoint rule over u ∈ [0, c] failed for thin exponential layers, c = R/h ≫ 1.)
"""
radial_kappa(shape, h, R) = ring_integral(shape, h, R, R - radial_support(shape, h, R), R) / (h * R)

"""
    ring_integral(shape, h, R, lo, hi)

∫_lo^hi g((R − r)/h) r dr exactly (antiderivatives in r), clipped to the support [R − L, R].
step: r²/2; linear: 2[(1 − R/h) r²/2 + r³/(3h)]; exponential: h e^((r − R)/h)(r − h)/(1 − e^(−c)).
"""
function ring_integral(shape, h, R, lo, hi)
    L = radial_support(shape, h, R)
    a, b = max(lo, R - L), min(hi, R)
    b <= a && return 0.0
    if shape == "step"
        return (b^2 - a^2) / 2
    elseif shape == "linear"
        F(r) = 2 * ((1 - R / h) * r^2 / 2 + r^3 / (3h))
        return F(b) - F(a)
    else
        c = L / h
        G(r) = h * exp((r - R) / h) * (r - h)
        return (G(b) - G(a)) / -expm1(-c)
    end
end

"Symmetric operator M (1/ms) and ring areas A for the radial problem (see the file header)."
function radial_operator(; R, D, layer=nothing, surface=0.0, N=2000)
    dr = R / N
    redge = collect(0:N) .* dr
    A = [π * (redge[i + 1]^2 - redge[i]^2) for i in 1:N]
    V = zeros(N)
    if !isnothing(layer)
        (shape, ρ, h) = layer
        κ = radial_kappa(shape, h, R)
        for i in 1:N
            V[i] = ρ / (h * κ) * 2π * ring_integral(shape, h, R, redge[i], redge[i + 1]) / A[i]
        end
    end
    c = [2π * D * redge[i + 1] / dr for i in 1:N - 1]
    K = zeros(N)
    for i in 1:N - 1
        K[i] += c[i]
        K[i + 1] += c[i]
    end
    # Robin condition on the face r = R (final review I3): centre → face conductance 2πDR/(dr/2) in series with the
    # surface sink 2πRρ, i.e. 2πRρ/(1 + ρ dr/(2D)); second order in dr (the sink on the ring centre was first order)
    K[N] += 2π * R * surface / (1 + surface * dr / (2D))
    M = SymTridiagonal((K .+ A .* V) ./ A, -c ./ sqrt.(A[1:N - 1] .* A[2:N]))
    return (M, A)
end

"Decay rates λ_k (1/ms, ascending) of the discretised radial operator."
radial_rates(; kw...) = eigvals(radial_operator(; kw...)[1])

"""
    radial_signal(; R, D, times, layer=nothing, surface=0.0, R2_bulk=0.0, N=2000)

S(t)/S(0) inside a cylinder of radius R for a uniform start. `layer = (shape, ρ, h)` (ρ per unit surface area, D3.1),
`surface = ρ` for the conventional model.
"""
function radial_signal(; R, D, times, layer=nothing, surface=0.0, R2_bulk=0.0, N=2000)
    (M, A) = radial_operator(; R=R, D=D, layer=layer, surface=surface, N=N)
    F = eigen(M)
    c2 = (F.vectors' * sqrt.(A)) .^ 2 ./ sum(A)
    return [sum(c2 .* exp.(-F.values .* t)) * exp(-R2_bulk * t) for t in times]
end

"Exact lowest rate with surface relaxivity ρ (test only): D k², D k J₁(kR) = ρ J₀(kR), k ∈ (0, j₀₁/R), by bisection."
function exact_rate_conventional(R, D, ρ)
    f(k) = D * k * besselj1(k * R) - ρ * besselj0(k * R)
    lo, hi = 0.0, 2.404825557695773 / R
    for _ in 1:200
        mid = (lo + hi) / 2
        f(mid) < 0 ? (lo = mid) : (hi = mid)
    end
    return D * ((lo + hi) / 2)^2
end

"""
    exact_rate_step(R, D, ρ, h)

Exact lowest rate for the step layer inside (test only), 0 < h < R. Core r < r_a = R − h: m = J₀(kr), λ = Dk². Shell
(rate q = ρ/(h(1 − h/2R))): m = I₀(βr) + (I₁(βR)/K₁(βR))K₀(βr), β = √((q − λ)/D), which satisfies ∂_r m(R) = 0. The
lowest λ lies in (0, q) (λ₁ ≤ mean rate 2ρ/R < q) and makes m'/m continuous at r_a; found by bisection.
"""
function exact_rate_step(R, D, ρ, h)
    q = ρ / (h * (1 - h / (2R)))
    ra = R - h
    function F(λ)
        k = sqrt(λ / D)
        β = sqrt((q - λ) / D)
        Lin = -k * besselj1(k * ra) / besselj0(k * ra)
        ratio = besseli(1, β * R) / besselk(1, β * R)
        Lsh = β * (besseli(1, β * ra) - ratio * besselk(1, β * ra)) / (besseli(0, β * ra) + ratio * besselk(0, β * ra))
        return Lin - Lsh
    end
    lo, hi = 0.0, q * (1 - 1e-12)
    @assert F(lo) > 0 && F(hi) < 0 "exact_rate_step: root not bracketed"
    @assert sqrt(q / D) * ra < 2.404825557695773 "exact_rate_step: core J₀ would change sign"
    for _ in 1:200
        mid = (lo + hi) / 2
        F(mid) > 0 ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end
