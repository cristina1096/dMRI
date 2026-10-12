# Continuum benchmark (review C1): exact signal of the same models on walls, without Monte Carlo or timestep.
#
# 1D diffusion on [0, w] with reflecting walls (the x-motion between Walls(repeats=w); y, z do not matter), uniform start,
# M⊥ obeying ∂m/∂t = D ∂²m/∂x² − [R2_bulk + ΔR2(x)] m. The layer contributes ΔR2(x) = (ρ/h)[g(x/h) + g((w − x)/h)]
# (both faces; exponential truncated at w and renormalised, as in the package). The conventional model contributes a
# surface sink with relaxivity ρ on both walls (MCMR's θ with ρ = θ√(D/π), P2.5.0). Discretisation: N finite-volume cells,
# Neumann Laplacian, ΔR2 averaged over each cell by midpoint sub-sampling, the sink as ρ/dx in the two wall cells;
# S(t) = Σ_k ⟨u0, φ_k⟩² e^(−λ_k t) e^(−R2_bulk t) from the eigen-decomposition. Independent of the package's layer code.

using LinearAlgebra

"Normalised profile g(u) (∫ = 1 over its support), written independently of `Layers.profile_g`."
function continuum_g(shape::AbstractString, u::Real, c::Real)
    shape == "step" && return 0 <= u <= 1 ? 1.0 : 0.0
    shape == "linear" && return 0 <= u <= 1 ? 2 * (1 - u) : 0.0
    shape == "exponential" && return 0 <= u <= c ? exp(-u) / (1 - exp(-c)) : 0.0
    error("Unknown shape $shape")
end

function continuum_operator(; w, D, layer=nothing, surface=0.0, N=3000, subsamples=16)
    dx = w / N
    V = zeros(N)
    if !isnothing(layer)
        shape, ρ, h = layer
        for i in 1:N, j in 1:subsamples
            x = (i - 1 + (j - 0.5) / subsamples) * dx
            V[i] += ρ / h * (continuum_g(shape, x / h, w / h) + continuum_g(shape, (w - x) / h, w / h)) / subsamples
        end
    end
    V[1] += surface / dx
    V[N] += surface / dx
    dg = fill(2D / dx^2, N); dg[1] = dg[N] = D / dx^2
    return SymTridiagonal(dg .+ V, fill(-D / dx^2, N - 1))
end

"Decay rates λ_k (1/ms, ascending) of the discretised operator."
continuum_rates(; kw...) = eigvals(continuum_operator(; kw...))

"""
    continuum_signal(; w, D, times, layer=nothing, surface=0.0, R2_bulk=0.0, N=3000)

Signal S(t)/S(0) for a uniform start. `layer = (shape, ρ, h)` for the near-surface layer, `surface = ρ` for the
conventional surface relaxation (both may be combined; neither gives free diffusion).
"""
function continuum_signal(; w, D, times, layer=nothing, surface=0.0, R2_bulk=0.0, N=3000)
    F = eigen(continuum_operator(; w=w, D=D, layer=layer, surface=surface, N=N))
    c2 = (F.vectors' * fill(1 / sqrt(N), N)) .^ 2
    return [sum(c2 .* exp.(-F.values .* t)) * exp(-R2_bulk * t) for t in times]
end
