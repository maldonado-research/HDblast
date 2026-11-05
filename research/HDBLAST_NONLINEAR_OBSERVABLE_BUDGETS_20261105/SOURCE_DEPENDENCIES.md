# Retained analytic inputs and the missing numerical joins

The source copies in [INPUT_PROVENANCE.json](INPUT_PROVENANCE.json) are
authenticated retained premises. Their programs and nonlinear PDE campaigns
are not rerun by this checkpoint's algebra checkers.

| Retained input | What it supplies conditionally | What is still needed here |
| --- | --- | --- |
| [Exact quasilinear family](inputs/EXACT_QUASILINEAR_FAMILY.md), §§1–5 | Exact finite-regularity evolution, state/hyperbolic/subsidiary assumptions, and a common sufficiently small parameter interval | Evaluated existence/state constants, a numerical parameter radius, and the particular reflected past-domain application |
| [Variation completion](inputs/FINITE_VARIATION_COMPLETION.md), §§2–4 | Genuine uniform first `B13` and second `B12` variation bounds, and the final `C^2` parameter curve into `X10` | Numerical values of their constants and the component/norm comparison bounds needed for full shell fields and observer transport |
| [Registered family binding](inputs/REGISTERED_EXACT_FAMILY_BINDING.md), §§2–4 | One selected affine scalar seed and weighted exact constraint curve | Quantitative inverse/state/flow margins and an admitted extension of the stated future slab to the required reflected past tube |
| [Intrinsic observer and signs](inputs/INTRINSIC_OBSERVER_AND_SIGNS.md), §§1, 4–6 | Anchored intrinsic geodesic observers and family-wise qualitative finite-parameter consequences | A numerical common noncrossing proper-time tube, observer/Jacobi derivative bounds, and uniform parameter bounds on the present past-central interval |
| [Nonlinear bridge](inputs/NONLINEAR_BRIDGE.md), §§2–5 | Precise selected observable, conditional past-slab application, coarse-cell endpoint guard and selection-dependent oddness | Evaluated error budgets and all common-domain hypotheses; the guard does not transfer unchanged to the newer fine linear bracket |

The source audit distinguishes those actual uniform analytic variation bounds
from quantities evaluated only at parameter zero. A first tangent, a complete
second Taylor coefficient at zero, a numerical first-field gradient budget,
or a restricted difference between preparations does not by itself bound the
chosen nonlinear family's derivatives uniformly over a parameter interval.

The constructive dependency chain for the notes is

```
specified exact preparation + admitted reflected evolution and parameter tube
 -> uniform full shell metric/scalar parameter and spatial derivative budgets
 -> augmented proper-time observer/Jacobi parameter-variation budgets
 -> full cut metric G',G'' and full scalar A',A'' bounds
 -> positive area/denominator and normalized moment remainder bounds
 -> sufficient central moment-zero guard, if its numerical inputs pass.
```

An odd-residual route additionally needs an effective continuity modulus of
the appropriate full second parameter derivative, or a directly validated
projected residual budget. Conditional `C^2` regularity supplies qualitative
continuity, not a numerical Hölder/Lipschitz constant or a chosen power law.

The final notes deliberately preserve this boundary. The previous complete
second-order construction and conditional exact nonlinear existence results
already exist under their stated premises. This checkpoint does not replace
them with a new theorem, evaluate their constants, or convert coordinate
coefficients into physical energy, radiation, or observations.
