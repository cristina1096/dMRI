# 1D toy random walk between reflecting walls at 0 and w, with a step-profile layer of thickness h on both
# faces, for the V11 comparison (P2.3.7). The x-motion of MCMR's 3D walk between Walls(repeats=w) is exactly
# this 1D walk (step ~ N(0, 2Dτ), mirror reflection), so the toy's :exact rule must agree with MCMR.

const _Layers = mr.Geometries.Internal.Layers

"""
    toy_layer_signal(; h, rate, tau, T=50.0, nspins, seed, method)

Per-spin final M⊥ after time T. `method = :exact`: on every reflected straight piece, add
rate·dt_piece·(time fraction within h of each face) using `segment_fraction`; `:endpoint`: add
rate·τ·(number of layers containing the end point of the step). Spin i uses its own RNG Xoshiro(seed·10⁷ + i),
so results do not depend on the thread count.
"""
function toy_layer_signal(; h, rate, tau, T=50.0, nspins, seed, method, w=W_REF, D=D_REF)
    method in (:exact, :endpoint) || error("method must be :exact or :endpoint, got $method")
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
                    expo += rate * dt * (_Layers.segment_fraction(pts[k], pts[k + 1], 0.0, h) +
                                         _Layers.segment_fraction(w - pts[k], w - pts[k + 1], 0.0, h))
                end
            else
                expo += rate * tau * ((x <= h) + (w - x <= h))
            end
        end
        out[i] = exp(-expo)
    end
    return out
end
