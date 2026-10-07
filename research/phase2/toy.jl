# 1D toy random walk between reflecting walls at 0 and w, with a layer of thickness h on both faces
# (step profile, or linear for P2.6.4), for the V11 comparison (P2.3.7). The x-motion of MCMR's 3D walk between Walls(repeats=w) is exactly
# this 1D walk (step ~ N(0, 2Dτ), mirror reflection), so the toy's :exact rule must agree with MCMR.

const _Layers = mr.Geometries.Internal.Layers

"""
    toy_layer_signal(; h, rate, tau, T=50.0, nspins, seed, method, shape="step")

Per-spin final M⊥ after time T; `rate` is ΔR2(0). `method = :exact`: on every reflected straight piece, add
∫ΔR2 dt over the piece for each face: for the step profile rate·dt_piece·(time fraction within h) using
`segment_fraction` (unchanged from P2.3.7); linear or exponential (cutoff = w) use MCMR's `side_exponent` with
ρ = rho_for(rate, h; shape). `rho_for` assumes cutoff = W_REF, so the toy's w must equal W_REF for the exponential.
`:endpoint`: add τ·ΔR2(end point of the step), summed over both faces. Spin i uses its own RNG Xoshiro(seed·10⁷ + i),
so results do not depend on the thread count.
"""
function toy_layer_signal(; h, rate, tau, T=50.0, nspins, seed, method, w=W_REF, D=D_REF, shape="step")
    method in (:exact, :endpoint) || error("method must be :exact or :endpoint, got $method")
    shape in ("step", "linear", "exponential") || error("shape must be step, linear or exponential, got $shape")
    side = _Layers.LayerSide(rho_for(rate, h; shape=shape), h, Symbol(shape), shape == "exponential" ? w : h)
    nsteps = round(Int, T / tau)
    σ = sqrt(2 * D * tau)
    out = ones(nspins)
    iszero(h) && return out
    Threads.@threads for i in 1:nspins
        rng = Random.Xoshiro(seed * 10_000_000 + i)
        x = rand(rng) * w
        expo = 0.0
        pts = Float64[]
        for _ in 1:nsteps
            step = randn(rng) * σ
            empty!(pts); push!(pts, x)
            while true
                if x + step < 0
                    step = -(step + x); x = 0.0; push!(pts, x)
                elseif x + step > w
                    step = -(step - (w - x)); x = w; push!(pts, x)
                else
                    x += step; push!(pts, x); break
                end
            end
            if method == :exact
                total = sum(abs(pts[k + 1] - pts[k]) for k in 1:length(pts) - 1)
                total > 0 || continue
                for k in 1:length(pts) - 1
                    dt = tau * abs(pts[k + 1] - pts[k]) / total
                    if shape == "step"
                        expo += rate * dt * (_Layers.segment_fraction(pts[k], pts[k + 1], 0.0, h) +
                                             _Layers.segment_fraction(w - pts[k], w - pts[k + 1], 0.0, h))
                    else
                        expo += _Layers.side_exponent(side, pts[k], pts[k + 1], dt) +
                                _Layers.side_exponent(side, w - pts[k], w - pts[k + 1], dt)
                    end
                end
            elseif shape == "step"
                expo += rate * tau * ((x <= h) + (w - x <= h))
            else
                expo += tau * side.rho / h * (_Layers.profile_g(side.shape, x / h, side.cutoff / h) +
                                              _Layers.profile_g(side.shape, (w - x) / h, side.cutoff / h))
            end
        end
        out[i] = exp(-expo)
    end
    return out
end
