# [Near-surface R2 layer](@id near_surface_layer)

Besides the per-collision [`MCMRSimulator.surface_relaxation`](@ref), [`Walls`](@ref) and [`Cylinders`](@ref) can carry a
**near-surface layer**: an extra transverse relaxation rate that depends on the distance ``d`` of a spin from the surface,

```math
R_{2,\mathrm{total}}(\mathbf{x}) = R_{2,\mathrm{bulk}} + \Delta R_2(d(\mathbf{x})), \qquad
\Delta R_2(d) = \frac{\rho}{h}\,\frac{g(d/h)}{\kappa},
```

where ``d`` is measured from the surface into the compartment, ``g`` is a normalised profile (``\int g(u)\,du = 1`` over
its support), ``h`` its length scale and ``\rho`` the integrated relaxivity (µm/ms).
The layer relaxes spins continuously while they are near the surface, instead of at the moment of a collision.

## Normalisation

- **Walls.** ``\kappa = 1``, so ``\rho = \int \Delta R_2\,\mathrm{d}d`` is the rate integrated along the normal.
- **Cylinders.** ``\rho`` is the relaxivity **per unit surface area**, ``\rho = (1/A)\int \Delta R_2\,\mathrm{d}V``. The curvature factor is
  ```math
  \kappa = \int g(u)\,(1 \mp u\,h/R)\,du = 1 \mp \bar c/R,
  ```
  with ``-`` inside and ``+`` outside the cylinder of radius ``R``. Here ``\bar c = h\int u\,g(u)\,du`` is the profile's centroid: ``h/2`` for step, ``h/3`` for linear, and ``h\,(1-(c+1)e^{-c})/(1-e^{-c})`` for an exponential with support ``L = c\,h``.
- **Mean exponent.** With this normalisation, for any profile and thickness, the mean exponent of a compartment of surface-to-volume ratio ``S/V`` grows as ``E[-\ln M_\perp] = \rho\,(S/V)\,t``. The conventional surface model obeys the same law when ``\theta = \rho\sqrt{\pi/D}``.

## Profiles

| `layer_shape` | ``g(u)`` | support |
|---|---|---|
| `"step"` (default) | 1 | ``0 \le d \le h`` |
| `"linear"` | ``2(1-u)`` | ``0 \le d \le h`` |
| `"exponential"` | ``e^{-u}/(1-e^{-L/h})`` | ``0 \le d \le L`` |

Exponential support ``L``:
- **walls:** `layer_cutoff`, by default the smallest spacing between walls (repeating or not);
- **cylinders, inside:** the radius, so the profile reaches the axis;
- **cylinders, outside:** `layer_cutoff`, by default the smallest surface gap to a neighbouring cylinder.

A single wall or cylinder without neighbours needs an explicit `layer_cutoff` for an exponential layer.

## Keywords

**Walls** (positive side: where the local coordinate is larger than the wall position; negative side: smaller):

| keyword | meaning |
|---|---|
| `layer_rho`, `layer_h`, `layer_shape` | ``\rho`` (µm/ms), ``h`` (µm) and the profile on both sides; `layer_rho = 0` (default) means no layer |
| `layer_rho_positive`, `layer_rho_negative`, `layer_h_positive`, `layer_h_negative`, `layer_shape_positive`, `layer_shape_negative` | side-specific values, overriding the both-sides ones |
| `layer_cutoff` | support of the exponential profile (µm) |

**Cylinders** (inner = inside the cylinder, outer = outside):

| keyword | meaning |
|---|---|
| `layer_rho`, `layer_h`, `layer_shape` | ``\rho`` (µm/ms, per unit surface area), ``h`` (µm) and the profile inside and outside |
| `layer_rho_inner`, `layer_rho_outer`, `layer_h_inner`, `layer_h_outer`, `layer_shape_inner`, `layer_shape_outer` | side-specific values, overriding the both-sides ones |
| `layer_cutoff` | outside support of the exponential profile (µm) |

**Limits on the parameters:**
- walls: for step and linear ``h`` (the support) must not exceed the spacing between walls; for the exponential ``h`` is a decay length and only `layer_cutoff` is limited by the spacing;
- cylinders, inside: step and linear need ``h \le R`` (for the exponential ``h`` is a decay length and may exceed ``R``);
- cylinders, outside: the support (``h`` for step and linear, `layer_cutoff` for the exponential) must not exceed the surface gap to any neighbouring cylinder, so that no layer reaches inside another cylinder.

From the command line the same keywords are flags, e.g. `mcmr geometry create cylinders --layer_rho 0.02 --layer_h 0.2 ...`
(see `mcmr geometry create cylinders --help`).

```julia
using MCMRSimulator
# square array of cylinders with an exponential layer outside only
geometry = Cylinders(radius=1., repeats=[2.5, 2.5], layer_rho_outer=0.02, layer_h_outer=0.2, layer_shape="exponential")
simulation = Simulation([]; geometry=geometry, diffusivity=3., R2=1/2000)
```

## Behaviour

- **Overlap.** Where the layers of several surfaces overlap (narrow gaps, thick layers), their rates **add**.
- **Integration along a step.**
  - Walls: every profile is integrated in closed form (the distance changes linearly along a straight step).
  - Cylinders: the step profile is integrated exactly from line–circle intersections; linear and exponential use adaptive Gauss–Legendre quadrature between the exact shell boundaries, to an internal tolerance of 1e-10.
  - Axial motion along a cylinder does not change ``d``.
- **Timestep.** The default maximum timestep keeps ``\Delta R_{2,\max}\,\tau \le 0.005`` (set by `Simulation(...; timestep=(layer=0.005,))`). This bounds the error of a single model's decay rate to about 0.15 %.
  - ``\Delta R_{2,\max}`` is an upper bound on the summed rate where layers overlap: for walls, the sum of the two facing surface rates across a gap whose supports overlap; for cylinders, the surface rate plus those of every neighbour whose outer support reaches (conservative, up to several times the true peak). Without overlap it is the largest single surface rate. Note that with the default exponential cutoff (the spacing) the facing exponential layers of neighbouring walls always overlap, so for such setups the default timestep is set by the summed bound (e.g. half the single-layer value for equal layers on both sides).
  - For **comparing** models or profiles at the same ``\rho`` (differences of a few 1e-3 in the signal), set the timestep explicitly, ``\tau = \min(0.135\,\bar c_{\min}^2,\ 0.005/\Delta R_{2,\max})`` at ``D = 3`` (in general ``\sqrt{2D\tau} \le 0.9\,\bar c_{\min}``), which gives ≲ 10 % error in the difference. An explicit `timestep=τ` replaces every automatic limit, so the rate limit must be included by hand.
- **Not supported:**
  - finite RF pulses together with a layer (an error is raised); use instantaneous pulses;
  - spins stuck to a surface (surface density > 0) do not feel the layer (a warning is given).

## Caveats

- **Step layers that overlap double the rate** in the overlap region, which is not a realistic near-surface profile.
- **A truncated exponential on curved arrays has an interior jump.**
  - The support of a neighbouring cylinder (a circle of radius ``R+L``) cuts into the gap away from the line of centres. There the neighbour's tail switches on with a jump of ``e^{-L/h}`` relative to the surface value.
  - Near the surface the jump is negligible for ``h \ll L``, up to ``e^{-1}`` for ``h = L``, and tends to 1 for ``h \gg L``. At the edges of the overlap lens it reaches twice the local rate.
  - Walls with the default cutoff (the spacing) have no such jump, because the cut falls on the opposite wall; a shorter user `layer_cutoff` creates the same jump on walls.
- **Linear layers remain realistic** when they overlap.
