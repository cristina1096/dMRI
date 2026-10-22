# Overlapping near-surface layers (P2.5.9): the summed rate ΔR2(x) across the gap from both faces, from the package and in
# closed form; realism check (never higher mid-gap than at the wall); overlap measures.

using StaticArrays: SVector
const _LY = mr.Geometries.Internal.Layers

"ΔR2(x) as the simulator applies it: layer exponent of a stationary spin at x for dt = 1 (both faces, all wall images)."
function summed_rate_package(shape, ρ, h, w, x)
    g = Simulation([]; geometry=Walls(repeats=w, layer_rho=ρ, layer_h=h, layer_shape=shape), verbose=false).geometry
    p = SVector{3, Float64}(x, 0.3, -1.2)
    return _LY.layer_exponent(g, p, p, 1.0)
end

"Closed-form ΔR2(x) = (ρ/h)[g(x/h) + g((w − x)/h)]; exponential truncated at w and renormalised."
function summed_rate_closed(shape, ρ, h, w, x)
    c = w / h
    g(u) = shape == "step" ? (0 <= u <= 1 ? 1.0 : 0.0) :
           shape == "linear" ? (0 <= u <= 1 ? 2 * (1 - u) : 0.0) :
           shape == "exponential" ? (0 <= u <= c ? exp(-u) / (1 - exp(-c)) : 0.0) : error("Unknown shape $shape")
    return ρ / h * (g(x / h) + g((w - x) / h))
end

"""
    overlap_metrics(shape, h, w)

overlap_width: length of the region where both faces contribute (linear/step: max(0, 2h − w); exponential: w).
mid_over_wall: ΔR2(w/2) / ΔR2(0). opposite_share: fraction of one face's integral ∫(ρ/h)g(d/h)dd lying beyond the
middle of the gap (d > w/2).
"""
function overlap_metrics(shape, h, w)
    width = shape == "exponential" ? w : max(0.0, 2h - w)
    r(x) = summed_rate_closed(shape, 1.0, h, w, x)
    if shape == "exponential"                                    # one face, truncated at w: share beyond w/2 in closed form
        c = w / h
        share = (exp(-c / 2) - exp(-c)) / (1 - exp(-c))
    else
        share = shape == "linear" ? (h <= w / 2 ? 0.0 : (1 - w / (2h))^2) : max(0.0, 1 - w / (2h))
    end
    # mid value just left of w/2: a step support edge exactly at w/2 (h = w/2) would otherwise be counted from both faces
    return (overlap_width=width, mid_over_wall=r(w / 2 * (1 - 1e-12)) / r(0.0), opposite_share=share)
end

"ΔR2 non-increasing from the wall to the middle of the gap (n points), relative tolerance 1e-12."
function is_realistic(shape, ρ, h, w; n=401)
    # interior samples: at exactly x = 0 or x = w/2 a support edge can be counted from both sides (measure zero)
    v = [summed_rate_closed(shape, ρ, h, w, x) for x in range(1e-9, w / 2 - 1e-9, length=n)]
    return all(v[k + 1] <= v[k] * (1 + 1e-12) + 1e-15 for k in 1:n - 1)
end
