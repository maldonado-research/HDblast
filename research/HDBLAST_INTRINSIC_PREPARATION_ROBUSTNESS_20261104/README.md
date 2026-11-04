# Intrinsic preparation data and robust shell-volume coefficients

4 November 2026. The preceding checkpoint showed that different admissible second
preparations of the same first packet can give opposite shell responses. Here we
define a geometrical initial-data class and prove that the original negative
volume and expansion coefficients persist throughout that class, under the
retained analytic assumptions.

The data are areas of scalar-field level surfaces **inside the specified initial
Cauchy slice**, continuing the selected annular three-spheres. The actual first
fields may be nonspherical. For two families sharing those first fields, let
`D` be the difference of their second area derivatives divided by the background
area, at the same fixed scalar label. Its derivative uses the declared background
scalar-label parametrization, not a freely chosen perturbed radial coordinate.

On the causal initial segment `1<=r<=5/4`, define

```
R_A=max(sup|D|,sup|D_r|),
D_r=(-ell*rho*p)*partial_phi_label D, ell=10^-6.
```

For the complete added preparation, the closed class `R_A<=2e33` has the following
safe outward second-parameter-derivative bounds at the physical shell and fixed
proper time `tau_star=rho_b*ell*(5/4)`:

| Shell observable | Whole-class enclosure |
|---|---|
| Relative volume `v2=V2/V0` | `[-7.6900e16,-4.1600e15]` |
| Proper-time expansion `Theta2` | `[-1.1070e22,-2.3980e21]` |

For each admitted exact family separately, sufficiently small nonzero
perturbations of either sign therefore give positive volume and positive
expansion below the background values at that time. The class supplies no common
numerical perturbation-amplitude interval and no signed finite-time history.
Its large numbers are model-normalized second coefficients, not finite fields
or SI measurements.

## What is new in this checkpoint

The scalar-level area difference is unchanged by coordinate transformations
transporting the chosen initial slice, scalar and metric. This does not make it
independent of the choice of physical slice, level component, background labels
or parameter normalization. Those choices are declared parts of the preparation.

In the retained compact metric-mean family, the area profile reconstructs the
initial mean through a two-component linear constraint system. A continuous
Gronwall estimate, all five conditioned initial channels, and the retained coupled
mean resolvent give `N<=R_A/2000`. Thus the closed class has `N<=1e30`, strictly
inside the previously established sufficient sign guard. This is a bound on the
complete increment: for `base+lambda*h`, use `|lambda|R_A(h)`.

The opposite-sign example from the earlier checkpoint lies outside this class.
At the initial scalar label corresponding to `r=5/4`, its area-profile magnitude
alone exceeds `3.7621e33`. This initial bulk label is distinct from the future
proper-time target with the same numerical value of T. The two results are
compatible: unrestricted preparation freedom can change the signs, while a
specified smaller class has stable coefficients.

The geometric construction uses established scalar–metric and hypersurface
perturbation machinery. The contribution is the constrained reconstruction and
quantitative class estimate in this model. No physical mechanism selecting the
class, and no physical justification of its radius, is supplied.

Read the [derivation and assumptions](DERIVATION.md), [primary-source comparison](REFERENCES.md),
and [next targets](NEXT_STEP.md). The working paper is available as
[PDF](manuscript/INTRINSIC_PREPARATION_ROBUSTNESS.pdf) and
[editable TeX](manuscript/INTRINSIC_PREPARATION_ROBUSTNESS.tex) for author and
specialist review. Earlier papers and capsules remain unchanged.

## Reproduce and interpret

Use Python 3.12.14 with the pinned source closures:

```sh
python -I -S -B research/HDBLAST_INTRINSIC_PREPARATION_ROBUSTNESS_20261104/checks/run_checkpoint.py --repo-root . --output-dir /tmp/hdblast-intrinsic-preparation-fresh
```

Choose a fresh output directory outside the checkout. The root replay
authenticates the complete new and inherited inventories before project imports
and afterward. New canonical calculations run in normal and optimized Python.
Actual rejected subprocesses, successful controls and in-process mutants are
counted separately. All new arithmetic is standard-library exact rational
calculation; there are no new PDE solves, background/evolved-field oracle calls
or bump queries. Background enclosures and analytic theorems remain shared inputs.

The statement is conditional on the retained background, constraint inverse,
exact-family development, regularity, mixed linear uniqueness and anchored
observer premises. The level-area functional also requires the continuing regular
scalar-level component on the declared physical initial slice. Whole-packet
injectivity and causal-segment reconstruction are distinct: data known only on
the causal segment do not determine exterior preparations or higher-order fields.

Internal AI-assisted checking and exact arithmetic do not establish external
peer review, formal proof-assistant verification or originality. The higher-
dimensional Big Bang origin remains unverified; matter/radiation generation,
physical scale and observational tests remain open. No Zenodo edition or journal
submission accompanies this GitHub checkpoint.
