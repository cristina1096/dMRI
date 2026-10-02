"""
Near-surface transverse relaxation layer (proposal Eqs. 5–10):

    R2_total(d) = R2_bulk + ΔR2(d),    ΔR2(d) = (ρ/h)·g(d/h)

Implemented: planar walls, step profile g(u) = 1 for 0 ≤ u ≤ 1, no finite RF pulses.
Each side of a wall has its own (ρ, h).
"""
module Layers

"""
    LayerSide(rho, h)

Near-surface layer on one side of a surface: integrated relaxivity `rho` (um/ms) and length scale `h` (um).
A side with `rho == 0` has no layer.
"""
struct LayerSide
    rho :: Float64
    h :: Float64
end

"""
    surface_rate(side)

Excess transverse relaxation rate at the surface, ΔR2(0) = (ρ/h)·g(0), in 1/ms (g(0) = 1 for the step profile).
"""
surface_rate(side::LayerSide) = iszero(side.rho) ? 0.0 : side.rho / side.h

"""
    WallLayer(position, positive, negative)

Layers on the two sides of one wall at local coordinate `position` (um).
`positive` is the side where the local coordinate is larger than `position`.
"""
struct WallLayer
    position :: Float64
    positive :: LayerSide
    negative :: LayerSide
end

"""
    segment_overlap(a0, a1, lo, hi)

A coordinate changes linearly from `a0` (s = 0) to `a1` (s = 1) along a straight segment.
Returns `(s_enter, s_exit)`, the part of s ∈ [0, 1] where lo ≤ a(s) ≤ hi.
There is no overlap when `s_exit <= s_enter`.
"""
function segment_overlap(a0::Float64, a1::Float64, lo::Float64, hi::Float64)
    if a0 == a1
        return (lo <= a0 <= hi) ? (0.0, 1.0) : (0.0, 0.0)
    end
    s_lo = (lo - a0) / (a1 - a0)
    s_hi = (hi - a0) / (a1 - a0)
    return (max(0.0, min(s_lo, s_hi)), min(1.0, max(s_lo, s_hi)))
end

"""
    segment_fraction(a0, a1, lo, hi)

Fraction of the segment (see [`segment_overlap`](@ref)) spent with lo ≤ a(s) ≤ hi.
"""
function segment_fraction(a0::Float64, a1::Float64, lo::Float64, hi::Float64)
    (s_enter, s_exit) = segment_overlap(a0, a1, lo, hi)
    return max(0.0, s_exit - s_enter)
end

"""
    side_exponent(side, d0, d1, dt)

∫ΔR2(d(t)) dt over a straight segment of duration `dt` (ms).
`d0` and `d1` are the distances from the surface at the start and end, positive into this side.
For the step profile this is ΔR2(0)·(time spent with 0 ≤ d ≤ h).
"""
function side_exponent(side::LayerSide, d0::Float64, d1::Float64, dt::Float64)
    iszero(side.rho) && return 0.0
    return surface_rate(side) * dt * segment_fraction(d0, d1, 0.0, side.h)
end

end
