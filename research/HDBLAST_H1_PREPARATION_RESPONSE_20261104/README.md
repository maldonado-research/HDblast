# Integrated geometric control of shell-volume response

4 November 2026. This checkpoint proves a direct conditional bound from an
integrated norm of the initial geometric preparation to the two proper shell
observables. It uses the retained coupled five-mean equations and all their
physical boundary rows. Earlier papers and source capsules remain unchanged.

For two admitted families with the same background and first perturbation, let
`D` be their difference in second scalar-level area derivatives, divided by the
background area, on the specified physical initial Cauchy slice. Continue the
selected scalar-level component and use the fixed background scalar labels:

```
I=[1,5/4], ell=10^-6,
D_r=(-ell*rho*p)*partial_scalar_label D,
R_H1=sqrt(integral_I (D^2+D_r^2) dr).
```

This is the norm of the **complete added preparation**: for `base+lambda*h`,
use `|lambda|R_H1(h)`. The slice, level component, background calibration and
parameter normalization are part of the data. The norm is not a physical energy.

For every retained admissible preparation with the closed bound `R_H1<=2e31`,
the second-parameter-derivative coefficients at the physical shell and fixed
proper time `tau_star=rho_b*ell*(5/4)` satisfy:

| Complete observable | Safe outward enclosure |
|---|---|
| Relative volume `v2=V2/V0` | `[-7.2400e16,-8.6400e15]` |
| Proper-time expansion `Theta2` | `[-1.0113e22,-3.3540e21]` |

Each admitted exact family separately therefore has positive volume and positive
expansion below the background values at that time for sufficiently small
nonzero perturbations of either sign. The bounds supply no common or effective
nonlinear amplitude interval and no signed proper-time history. These large
numbers are model-normalized second coefficients, not finite fields or SI data.

## What changed

A positive tilted wave energy, including the shell boundary value as a state,
controls the scalar shell trace directly from `R_H1`. Its boundary flux treats
the actual mixed Dirichlet and nonsymmetric Robin rows; no symmetric or positive
Robin matrix is assumed. The resulting trace estimate is `|w(x,0)|<=R_H1/300`.
The complete proper-time Green formulas then give the displayed observable bounds.
This energy is a mathematical estimate, not a new Noether energy or positivity law.

The earlier pointwise derivative condition cannot be recovered from an H1 budget:
an explicit smooth, compact, exactly zero-moment sequence in the actual admissible
image has bounded `R_H1` but unbounded derivative supremum and old initial-supremum
quantity `N`. Fixed rescaling places that sequence within the new class. Thus the
new route covers preparations outside the previous pointwise criterion. No claim
that either specified norm ball contains the other is needed. Large derivatives
alone do not imply a large or unstable observable response.

The earlier opposite-sign preparation is outside this new class as well. Its
initial-level area magnitude, together with the exact one-dimensional H1 embedding,
implies `R_H1>1.7914e33`. This initial bulk label is distinct from the future proper-time
target with the same numerical value of T.

First-tangent quadratic Noether tests cannot select this new budget: their input
is shared while the second preparation varies. That is a restricted, standard
nonselection observation. A physical quantity sensitive to the second preparation,
or another physical preparation rule, remains necessary.

## Assumptions and reproduction

The continuous operator norm caps are explicitly retained **conditional
assertions** from the earlier analytic preflight. Authenticating their source
bytes does not reprove them. Background bounds, the constraint inverse, exact-family
development, regularity, mean projection, finite propagation and mixed linear
uniqueness are also retained. The new contribution is the direct energy argument,
initial norm estimate and quantitative transport under those premises.

Read the [derivation](DERIVATION.md), [primary-source comparison](REFERENCES.md),
[supporting preparation-selection note](support/PREPARATION_SELECTION.md), and
[next scientific targets](NEXT_STEP.md). The working manuscript is available as
[PDF](manuscript/H1_PREPARATION_RESPONSE.pdf) and [editable TeX](manuscript/H1_PREPARATION_RESPONSE.tex)
for author and specialist review.

With Python 3.12.14, choose a fresh output directory outside the checkout:

```sh
python -I -S -B research/HDBLAST_H1_PREPARATION_RESPONSE_20261104/checks/run_checkpoint.py --repo-root . --output-dir /tmp/hdblast-h1-preparation-response-fresh
```

The replay authenticates the new inventory and all 29 inherited closed source
capsules before scientific imports and afterward. It reproduces new rational
reports in normal and optimized Python and counts actual subprocess controls
separately from in-process checks. It does not repeat the inherited proof campaigns
or independently solve the Einstein–scalar PDE. New arithmetic uses no numerical
background, evolved-field or bump oracle queries.

The method follows established wave-energy and boundary-estimate ideas; originality
and significance require external specialist assessment. No physical rule selecting
the class or its radius, matter/radiation generation, physical scale, observational
test or higher-dimensional Big Bang origin is established. No new Zenodo edition
or journal submission accompanies this GitHub checkpoint.
