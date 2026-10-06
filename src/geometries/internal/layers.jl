"""
Near-surface transverse relaxation layer (proposal Eqs. 5–10):

    R2_total(d) = R2_bulk + ΔR2(d),    ΔR2(d) = (ρ/h)·g(d/h)

Implemented: planar walls, step profile g(u) = 1 for 0 ≤ u ≤ 1, no finite RF pulses.
Each side of a wall has its own (ρ, h).
"""
module Layers

import StaticArrays: SVector
import ..FixedObstructionGroups: FixedObstructionGroup, FixedGeometry, repeating, rotate_from_global

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

"""
    layer_exponent(group/geometry, start, dest, dt)

∫ΔR2(d(t)) dt (dimensionless) for a spin moving in a straight line from `start` to `dest` (global coordinates, um)
in `dt` ms. For a geometry, the contributions of all groups and all walls are added (proposal Eq. 10).
Returns 0 when no layer is defined.
"""
layer_exponent(g::FixedObstructionGroup, start::SVector{3, Float64}, dest::SVector{3, Float64}, dt::Float64) = layer_exponent(g.layer, g, start, dest, dt)

layer_exponent(::Nothing, ::FixedObstructionGroup, ::SVector{3, Float64}, ::SVector{3, Float64}, ::Float64) = 0.0

function layer_exponent(layers::Vector{WallLayer}, g::FixedObstructionGroup, start::SVector{3, Float64}, dest::SVector{3, Float64}, dt::Float64)
    x0 = rotate_from_global(g, start)[1]
    x1 = rotate_from_global(g, dest)[1]
    lo = min(x0, x1)
    hi = max(x0, x1)
    spacing = repeating(g) ? g.repeats[1] : 0.0
    total = 0.0
    for layer in layers
        if repeating(g)
            # wall images c = position + k·spacing whose layers can reach [lo, hi]
            reach = max(layer.positive.h, layer.negative.h)
            kmin = ceil(Int, (lo - reach - layer.position) / spacing)
            kmax = floor(Int, (hi + reach - layer.position) / spacing)
        else
            kmin = kmax = 0
        end
        for k in kmin:kmax
            c = layer.position + k * spacing
            total += side_exponent(layer.positive, x0 - c, x1 - c, dt)
            total += side_exponent(layer.negative, c - x0, c - x1, dt)
        end
    end
    return total
end

function layer_exponent(geometry::FixedGeometry, start::SVector{3, Float64}, dest::SVector{3, Float64}, dt::Float64)
    return sum(g -> layer_exponent(g, start, dest, dt), geometry; init=0.0)
end

"""
    has_layer(geometry)

Whether any group in the [`FixedGeometry`](@ref) has a near-surface layer.
"""
has_layer(geometry::FixedGeometry) = any(g -> !isnothing(g.layer), geometry)


"""
    max_layer_rate(geometry)

Largest near-surface rate ΔR2(0) = ρ/h (1/ms) over all layer sides with ρ > 0 in the geometry; 0 if there is no layer.
"""
max_layer_rate(::Nothing) = 0.0
max_layer_rate(layers::Vector{WallLayer}) = maximum((surface_rate(s) for l in layers for s in (l.positive, l.negative)); init=0.0)
max_layer_rate(geometry::FixedGeometry) = maximum((max_layer_rate(g.layer) for g in geometry); init=0.0)

end
