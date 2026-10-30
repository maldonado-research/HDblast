# Original bulk-source derivative and an evolved C2 error bound

This checkpoint bounds the complete time derivative of the original actual
quadratic bulk source. It then closes a conditional second-derivative error
estimate for the already published `W1_quadratic` reconstruction. The actual
source, target, numerical candidate, initial data and physical shell law are
unchanged. The separately archived source-averaged candidate is not used here.

## The missing original-source input

The new actual first-gauge estimate bounds all 15 mixed derivative slots
through order four for its three amplitudes. Every component's scaled
fourth-order mixed derivative is strictly below `0.0499` on
`0<=x<=3,r>=0,x+r<=8`. The proof retains the exceptional initial Picard
iterates, nonzero differentiated initial forcing, time-independent physical
boundary laws, and ordinary derivatives rather than normalized Taylor
coefficients. The earlier master fourth-derivative bound is retained.

These bounds supply ordinary third derivatives of all eight first-field
amplitudes needed by the complete geometric-source chain rule. The physical
radial conversion includes the derivative of `1/(ell rho)`. All explicit
temporal background factors and all 23 metric/scalar forcing normalizations
are differentiated. The result, in the fixed positive angular H norm, is

```
sup ||F_x||_H <905405545954284 <9.055e14,
x=t/ell, ell=1e-6.
```

This is a continuous bound on the original scaled-time source derivative
over the retained first-field triangle. It is not the response quantity
`X(KF+LG)`, a finite-difference stencil, or a derivative of independently
selected interval values. It does not assign a physical force or sign.

See the [actual gauge proof](theory/ACTUAL_GAUGE_FOURTH_DERIVATIVES.md),
[complete source derivation](theory/ORIGINAL_BULK_FX_NOTE.md),
[source report](reference/FX.json), and
[independent source review](review/INDEPENDENT_BULK_FX_REVIEW.md).

## Consequence for the unchanged reconstruction

Let `U_target` be the same original formal second-coefficient target, and
`e=U_target-W1_quadratic`. On

```
D+={0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x},
```

the existing C1 comparison and the new differentiated estimate give

| Error quantity | Uniform positive H-norm strict upper bound |
| --- | ---: |
| Each first scaled time/radial derivative | 3.541e25 |
| Each second derivative `e_xx,e_xr,e_rr` | 6.222e33 |

The existing value bound is also retained. The proof first establishes the
candidate's physical-side second derivatives using the actual source
derivative, then differentiates the error equation with the same actual
coefficient and boundary inverse. The actual affine shell source cancels
from the error law. The error boundary source and its time derivative are
zero, and `M_x=0`.

The intermediate bulk correction has a generally nonzero initial
acceleration; this is retained. For this original candidate the final
forcing defect itself has `r1(0,r)=0`, which is also retained. Neither
identity is imported from the source-averaged reconstruction.

The result concerns continuous derivatives on the physical side, with
one-sided shell and domain-edge traces. It asserts no C2 regularity across
the artificial reflected extension. See the
[error derivation](theory/ORIGINAL_TARGET_C2_NOTE.md),
[exact error report](reference/C2.json),
[independent numerical review](review/INDEPENDENT_C2_GATE_REVIEW.md), and
[regularity review](review/INDEPENDENT_C2_REGULARITY_REVIEW.md).

## What the size of the bound tells us

These are large absolute coordinate error bounds. No nonzero physical
comparison scale or relative accuracy has been established. In the exact
upper-bound ledger, the initial free-wave term accounts for more than
`0.999999999` of the differentiated numerator. The contribution of the new
`F_x` cap to the second-error upper bound is only about `4.37e11`.

Those fractions describe the proof majorant, not a decomposition or lower
bound of the unknown true error. They identify the initial-data estimates
as the next useful place to improve this certificate. Two inner nominal
second-derivative jumps and coarse interval dependencies remain. A changed
cutoff or initial continuation would define a new candidate and would need
its own comparison estimates.

## Reproduction and scope

From a complete checkout, using only Python's standard library:

```
python -I -S -B research/HDBLAST_ORIGINAL_BULK_SOURCE_DERIVATIVE_20261030/checks/run_checkpoint.py --repo-root . --output-dir /tmp/hdblast-original-source-fresh
```

Use a fresh output directory. The runner recomputes eight producer/review
reports in normal and optimized Python and compares them byte for byte
with their references. It also rejects three corrupted gauge formulas in
both modes. The source-role guard separately rejects four incorrect input
types in both modes. All mathematical acceptance decisions use explicit
checks and exact rational arithmetic. The independent source review uses
ordinary binomial derivative calculus instead of the producer's normalized
Taylor arithmetic and imports no producer code.

The results retain the registered actual background, master/gauge equations,
physical boundary laws, exact harmonic identity on a neighbourhood, source
map identities, compatible initial/corner data and upstream certificate
premises. Internal reproducibility and mathematical review are not outside
specialist approval or a novelty assessment.

This closes the previously missing original `F_x` input and the conditional
coordinate C2 error gate. It does not establish small covariant Einstein
residuals, independently verify every physical initial/subsidiary constraint,
construct a finite-amplitude nonlinear family, give a cosmological prediction
or establish a higher-dimensional Big Bang origin. Frozen manuscripts,
earlier candidates and public checkpoints are preserved. No Zenodo version
or journal submission accompanies this checkpoint.
