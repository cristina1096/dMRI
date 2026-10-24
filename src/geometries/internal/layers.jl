"""
Near-surface transverse relaxation layer (proposal Eqs. 5–10):

    R2_total(d) = R2_bulk + ΔR2(d),    ΔR2(d) = (ρ/h)·g(d/h)

Implemented: planar walls (all profiles below) and cylinders (step profile only); profiles step g(u) = 1 (0 ≤ u ≤ 1), linear g(u) = 2(1 − u) (0 ≤ u ≤ 1) and
exponential g(u) = e^(−u)/(1 − e^(−c)) (0 ≤ u ≤ c = cutoff/h, truncated and renormalised); no finite RF pulses.
Each side of a wall has its own (ρ, h).
"""
module Layers

import StaticArrays: SVector
import LinearAlgebra: norm
import ..FixedObstructionGroups: FixedObstructionGroup, FixedGeometry, repeating, rotate_from_global

"Profiles implemented in `profile_g`."
const SHAPES = (:step, :linear, :exponential)

"""
    LayerSide(rho, h[, shape[, cutoff]])

Near-surface layer on one side of a surface: integrated relaxivity `rho` (um/ms), length scale `h` (um), profile
`shape` (one of `SHAPES`; default `:step`) and `cutoff` (um), the distance from the surface beyond which g = 0.
For `:step` and `:linear` the cutoff is `h`. `:exponential` needs an explicit cutoff. A side with `rho == 0` has no layer.
"""
struct LayerSide
    rho :: Float64
    h :: Float64
    shape :: Symbol
    cutoff :: Float64
end
LayerSide(rho, h) = LayerSide(rho, h, :step)
function LayerSide(rho, h, shape::Symbol)
    shape == :exponential && error("The exponential profile needs a cutoff: LayerSide(rho, h, :exponential, cutoff).")
    return LayerSide(rho, h, shape, h)
end

"""
    profile_g(shape, u, c=Inf)

Normalised profile g(u), ∫ g(u) du = 1 over its support (proposal §3.2 III). Zero outside the support.
`:step`: g = 1 on [0, 1]. `:linear`: g = 2(1 − u) on [0, 1]. `:exponential`: g = e^(−u)/(1 − e^(−c)) on [0, c],
with c = cutoff/h (c = Inf: untruncated). `c` is ignored by step and linear.
"""
function profile_g(shape::Symbol, u::Real, c::Real=Inf)
    shape == :exponential && return (0 <= u <= c) ? exp(-u) / -expm1(-c) : 0.0
    (0 <= u <= 1) || return 0.0
    shape == :step && return 1.0
    shape == :linear && return 2 * (1 - u)
    error("Unknown layer profile $shape")
end

"""
    surface_rate(side)

Excess transverse relaxation rate at the surface, ΔR2(0) = (ρ/h)·g(0), in 1/ms.
"""
surface_rate(side::LayerSide) = iszero(side.rho) ? 0.0 : side.rho / side.h * profile_g(side.shape, 0.0, side.cutoff / side.h)

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

∫ΔR2(d(t)) dt over a straight segment of duration `dt` (ms), where `d0` and `d1` are the distances from the
surface at the start and end, positive into this side. d is linear along the segment, so the integral is exact
in closed form over the part of the segment inside the support 0 ≤ d ≤ cutoff:
- `:step`: ΔR2(0)·(time inside);
- `:linear`: g is linear along the segment, so its average over the part inside is its value at the midpoint;
- `:exponential`: the mean of e^(−d/h) for d uniform on [d_a, d_b] is e^(−d_near/h)·(1 − e^(−x))/x with
  d_near = min(d_a, d_b) and x = |d_b − d_a|/h, evaluated with `expm1` (series 1 − x/2 for x < 1e-8).
"""
function side_exponent(side::LayerSide, d0::Float64, d1::Float64, dt::Float64)
    iszero(side.rho) && return 0.0
    (s_enter, s_exit) = segment_overlap(d0, d1, 0.0, side.cutoff)
    fraction = max(0.0, s_exit - s_enter)
    iszero(fraction) && return 0.0
    side.shape == :step && return surface_rate(side) * dt * fraction
    if side.shape == :exponential
        # evaluate from the end nearer the wall, so the factor is direction-independent and never 0·Inf
        d_near = min(d0 + (d1 - d0) * s_enter, d0 + (d1 - d0) * s_exit)
        x = abs(d1 - d0) * fraction / side.h
        mean_e = exp(-d_near / side.h) * (x < 1e-8 ? 1 - x / 2 : -expm1(-x) / x)
        return side.rho / side.h / -expm1(-side.cutoff / side.h) * mean_e * dt * fraction
    end
    d_mid = d0 + (d1 - d0) * (s_enter + s_exit) / 2
    return side.rho / side.h * profile_g(side.shape, d_mid / side.h) * dt * fraction
end

"""
    CylinderLayer(centre, radius, inner, outer)

Near-surface layers of one cylinder of radius `radius` (um) with its axis at `centre` (local 2D coordinates of its
group, um): `inner` inside the cylinder (0 ≤ d = R − r ≤ h), `outer` outside it (0 ≤ d = r − R ≤ h).
ρ is per unit surface area (D3.1): ΔR2(d) = (ρ/h)·g(d/h)/κ with κ = 1 ∓ c̄/R (− inside, + outside).
"""
struct CylinderLayer
    centre :: SVector{2, Float64}
    radius :: Float64
    inner :: LayerSide
    outer :: LayerSide
end

"""
    centroid(side)

Profile centroid c̄ = h∫u g(u) du (um). Only the step profile (c̄ = h/2) is supported on cylinders so far (D3.2).
"""
function centroid(side::LayerSide)
    side.shape == :step && return side.h / 2
    error("Layer profile $(side.shape) is not supported on cylinders yet (step only).")
end

"""
    kappa(side, radius, inside)

Curvature normalisation κ = ∫g(u)(1 ∓ u·h/R) du = 1 ∓ c̄/R (D3.1): − inside the cylinder, + outside.
"""
kappa(side::LayerSide, radius::Float64, inside::Bool) = 1 + (inside ? -1 : 1) * centroid(side) / radius

"""
    cylinder_surface_rate(side, radius, inside)

ΔR2 at the cylinder surface on one side, (ρ/h)·g(0)/κ (1/ms); 0 without a layer.
"""
cylinder_surface_rate(side::LayerSide, radius::Float64, inside::Bool) =
    iszero(side.rho) ? 0.0 : surface_rate(side) / kappa(side, radius, inside)

"""
    disk_length(a, b, q, X)

Length of {s ∈ [0, 1] : a s² + b s + q ≤ X} for a ≥ 0. The squared distance from a cylinder axis along a straight
piece is this convex quadratic, so the set is one interval; its ends are the roots, computed with the
cancellation-free form t = −(b + sign(b)√disc)/2, roots t/a and (q − X)/t.
"""
function disk_length(a::Float64, b::Float64, q::Float64, X::Float64)
    if iszero(a)
        return q <= X ? 1.0 : 0.0
    end
    disc = b * b - 4 * a * (q - X)
    disc <= 0 && return 0.0
    t = -(b + copysign(sqrt(disc), b)) / 2
    (s1, s2) = minmax(t / a, (q - X) / t)
    return max(0.0, min(1.0, s2) - max(0.0, s1))
end

"""
    shell_fraction(p0, p1, centre, r_lo, r_hi)

Fraction of the straight piece p(s) = p0 + s(p1 − p0), s ∈ [0, 1] (2D, um), with r_lo ≤ |p(s) − centre| ≤ r_hi.
Exact: |p(s) − centre|² = a s² + b s + q, so the fraction is disk_length(r_hi²) − disk_length(r_lo²).
"""
function shell_fraction(p0::SVector{2, Float64}, p1::SVector{2, Float64}, centre::SVector{2, Float64}, r_lo::Float64, r_hi::Float64)
    u = p0 - centre
    Δ = p1 - p0
    a = Δ[1] * Δ[1] + Δ[2] * Δ[2]
    b = 2 * (u[1] * Δ[1] + u[2] * Δ[2])
    q = u[1] * u[1] + u[2] * u[2]
    if iszero(a)
        return (r_lo^2 <= q <= r_hi^2) ? 1.0 : 0.0
    end
    inner = r_lo > 0 ? disk_length(a, b, q, r_lo^2) : 0.0
    return max(0.0, disk_length(a, b, q, r_hi^2) - inner)
end

"""
    cylinder_exponent(layer, centre, p0, p1, dt)

∫ΔR2 dt of both layers of one cylinder (axis at `centre`, local 2D) over the straight piece p0 → p1 of duration `dt`
(ms). Step profile: ΔR2(0)·dt·(fraction of the piece inside the shell), exact.
"""
function cylinder_exponent(layer::CylinderLayer, centre::SVector{2, Float64}, p0::SVector{2, Float64}, p1::SVector{2, Float64}, dt::Float64)
    R = layer.radius
    total = 0.0
    if !iszero(layer.inner.rho)
        total += cylinder_surface_rate(layer.inner, R, true) * dt * shell_fraction(p0, p1, centre, R - layer.inner.h, R)
    end
    if !iszero(layer.outer.rho)
        total += cylinder_surface_rate(layer.outer, R, false) * dt * shell_fraction(p0, p1, centre, R, R + layer.outer.h)
    end
    return total
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
            reach = max(layer.positive.cutoff, layer.negative.cutoff)
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

"""
    image_centres(g, centre, lo, hi, reach)

Axis positions (local 2D) of the images of a cylinder at `centre` in group `g` that lie within `reach` of the box
[lo, hi] in each coordinate: lattice images for a repeating group, `centre` itself otherwise.
"""
function image_centres(g::FixedObstructionGroup, centre::SVector{2, Float64}, lo::SVector{2, Float64}, hi::SVector{2, Float64}, reach::Float64)
    rep = repeating(g)
    sx = rep ? g.repeats[1] : 0.0
    sy = rep ? g.repeats[2] : 0.0
    kx = rep ? (ceil(Int, (lo[1] - reach - centre[1]) / sx):floor(Int, (hi[1] + reach - centre[1]) / sx)) : (0:0)
    ky = rep ? (ceil(Int, (lo[2] - reach - centre[2]) / sy):floor(Int, (hi[2] + reach - centre[2]) / sy)) : (0:0)
    return (centre + SVector(i * sx, j * sy) for i in kx for j in ky)
end

function layer_exponent(layers::Vector{CylinderLayer}, g::FixedObstructionGroup, start::SVector{3, Float64}, dest::SVector{3, Float64}, dt::Float64)
    p0 = SVector{2, Float64}(rotate_from_global(g, start))
    p1 = SVector{2, Float64}(rotate_from_global(g, dest))
    lo = min.(p0, p1)
    hi = max.(p0, p1)
    total = 0.0
    for layer in layers
        for c in image_centres(g, layer.centre, lo, hi, layer.radius + layer.outer.h)
            total += cylinder_exponent(layer, c, p0, p1, dt)
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

Largest near-surface rate (1/ms) in the geometry: ΔR2(0) = (ρ/h)·g(0) over all wall layer sides with ρ > 0; for cylinders an upper bound on the summed rate where outer layers overlap (P3.2.6). 0 if there is no layer.
"""
max_layer_rate(::Nothing) = 0.0
max_layer_rate(layers::Vector{WallLayer}) = maximum((surface_rate(s) for l in layers for s in (l.positive, l.negative)); init=0.0)
max_layer_rate(g::FixedObstructionGroup) = max_layer_rate(g.layer, g)
max_layer_rate(::Nothing, ::FixedObstructionGroup) = 0.0
max_layer_rate(layers::Vector{WallLayer}, ::FixedObstructionGroup) = max_layer_rate(layers)

"""
Upper bound on the summed ΔR2 for cylinders: for each cylinder, the larger of its inner surface rate and its outer
surface rate plus the outer surface rates of every other cylinder (or image) whose outer shell intersects its outer
shell, |c_i − c_j| < R_i + h_i + R_j + h_j. Exact when no outer shells overlap. Inner shells overlap nothing
(validation: h_outer ≤ surface gap).
"""
function max_layer_rate(layers::Vector{CylinderLayer}, g::FixedObstructionGroup)
    best = 0.0
    for (i, li) in enumerate(layers)
        best = max(best, cylinder_surface_rate(li.inner, li.radius, true))
        own = cylinder_surface_rate(li.outer, li.radius, false)
        iszero(own) && continue
        total = own
        for (j, lj) in enumerate(layers)
            reach = li.radius + li.outer.h + lj.radius + lj.outer.h
            for c in image_centres(g, lj.centre, li.centre, li.centre, reach)
                (j == i && c == li.centre) && continue
                norm(c - li.centre) < reach && (total += cylinder_surface_rate(lj.outer, lj.radius, false))
            end
        end
        best = max(best, total)
    end
    return best
end

max_layer_rate(geometry::FixedGeometry) = maximum((max_layer_rate(g) for g in geometry); init=0.0)

end
