# Tighter initial-data and evolved coordinate error certificates

The conditional continuous second-derivative error upper bounds for the
existing `W1_quadratic` reconstruction improve from `6.222e33` to
`1.187e28`, more than **524,000 times tighter**. This changes the estimates;
the reconstruction and its actual error are unchanged. These remain large
absolute coordinate bounds, without demonstrated relative physical accuracy.

## The source of the improvement

The earlier order-five paired comparison was dominated by interval
dependency in one initial radial cell, `[65/64,66/64]`. Ten exact dyadic
leaves now cover that cell. All other 47 cells are retained. The three
actual-minus-nominal telescoping channels are bounded on every leaf,
packed in the same positive 23-amplitude angular `H` norm, maximized over
the complete cover and summed by the triangle inequality.

The same convolution and inner-jump allowances are retained, including
the zero nominal curvature jump at `65/64`. The same initial displacement,
outer quadratic continuation, rounded-node velocity, smoothing widths,
original target, first-order fields, bulk source and physical shell law
are used. No new candidate or cutoff is introduced.

| Uniform error quantity | Previous strict upper | New strict upper |
| --- | ---: | ---: |
| Initial displacement curvature `e0_rr` | 6.219e33 | 1.185e28 |
| Each evolved scaled second derivative `e_xx,e_xr,e_rr` | 6.222e33 | 1.187e28 |
| Each evolved scaled first derivative | 3.541e25 | 3.541e25, retained |

Here `e=U_target-W1_quadratic`, `x=t/ell`, `ell=1e-6`, and

```
D+={0<=x<=3/2, max(0,x-5/4)<=r<=7/4-x}.
```

The initial comparison applies on the retained `0<=r<=7/4` footprint.
The evolved estimates include physical-side one-sided shell and domain-edge
traces. The original complete `F_x` input and all absolute candidate and
forcing bounds are retained. Replacing only the initial curvature-error
bound in the authenticated differentiated comparison theorem gives the
new estimates. A second implementation checks the exact affine difference
from the preceding C2 theorem. No rounded-bound subtraction is used.

See the [initial refinement derivation](theory/REFINED_INITIAL_CURVATURE_NOTE.md),
[transfer theorem](theory/REFINED_INITIAL_C2_TRANSFER.md),
[initial report](reference/INITIAL.json), [evolved report](reference/C2.json),
and [independent initial audit](review/INDEPENDENT_INITIAL_REFINEMENT_REVIEW.md).

## Actual mismatch at the two inner joins

A separate calculation evaluates the actual initial curvature and the
smoothed nominal curvature at each exact join. For the positive even kernel
of width `delta=2^-30`,

```
u_s,rr(rc)=(kappa_minus+kappa_plus)/2+remainder,
abs(remainder_i)<=35delta/256*(M3_minus_i+M3_plus_i).
```

The exact one-sided nominal branches supply the traces and local third
derivative caps. Signed actual order-five jets and an independent ODE
recurrence supply actual curvature. The resulting error intervals exclude
zero as follows:

| Initial join | Sign-excluding components | Scalar mean error interval | Full H-error strict lower |
| --- | --- | --- | ---: |
| 65/64 | J0,J2,J4,w0,w2,w4 | (0.6114,0.6115) | 0.6942 |
| 66/64 | w0,w2,w4 | (-5.195e20,-5.193e20) | 5.896e20 |

Other component signs are inconclusive. At the first join both nominal
curvature traces are exactly zero; its nonzero error comes from clipping
the actual tail. At the second join curvature averaging contributes to
the scalar error. A nominal jump alone would not establish either result:
the actual-minus-smoothed intervals establish these pointwise lower bounds.
The corresponding zero-width limits also exclude zero. Narrowing the
kernel alone therefore does not remove these coordinate mismatches.

These lower bounds concern the fixed reconstruction's initial data. They
do not imply a physical obstruction, force or measured signal. They also
explain why an upper-bound refinement cannot remove every actual error.
See the [join derivation](theory/INNER_JOIN_INITIAL_CURVATURE_ERROR.md),
[complete signed report](reference/JOINS.json), and
[independent ODE review](review/INDEPENDENT_INNER_JOIN_REVIEW.md).

## Reproduction and remaining work

From a complete checkout, using Python's standard library:

```
python -I -S -B research/HDBLAST_INITIAL_CURVATURE_REFINEMENT_20261030/checks/run_checkpoint.py --repo-root . --output-dir /tmp/hdblast-initial-curvature-fresh
```

Use a fresh output directory outside the capsule. The runner authenticates
the capsule and retained inputs, reproduces seven reports in normal and
optimized Python, and compares their bytes with the references. It rejects
six corrupted join formulas and eight corrupted C2 input files in each mode.
Independent internal implementations recompute the join ODE jets, all ten
refined leaves, complete-domain aggregation and C2 transfer identities.
They share authenticated upstream mathematical premises; this is not
external specialist certification or a formal proof-assistant result.

The retained nominal third-derivative convolution allowance and the refined
cutoff comparison remain useful targets for further improvement. A changed
cutoff or initial reconstruction would define a distinct candidate and
require new comparison estimates and point tables. Relevant established
methods and a bounded current literature check are described in the
[methods comparison](review/FOCUSED_METHODS_COMPARISON.md).

All results retain the registered background, master/gauge and initial-data
identities, source maps, harmonic-neighbourhood, compatible initial/corner
and physical boundary assumptions of the preceding checkpoints. They do
not establish small covariant Einstein residuals, independently complete
every physical/subsidiary constraint, construct a finite-amplitude nonlinear
family, predict a cosmological observable or establish a higher-dimensional
origin of the Big Bang. Frozen manuscripts, earlier candidates and public
checkpoints are preserved. No Zenodo version or journal submission
accompanies this checkpoint.
