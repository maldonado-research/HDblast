# Causal quadrature of the first actual-bulk correction

This checkpoint numerically encloses `C=W1-W0` at the shell point
`(x,r)=(9/8,0)`. Here `W0` is the earlier boundary-compatible reconstruction,
and `W1` is the first actual-bulk update defined in
[the preceding checkpoint](../HDBLAST_ACTUAL_BULK_UPDATE_20261030/).
Both are auxiliary, integral-defined candidates with the same smoothed nominal
initial data. The conditional target `U` retains the original initial data.

The computation integrates interval ranges over entire panels. Every panel
weight is an exact positive rational number. Actual variable coefficient
functions are assembled over each panel; the actual constant shell matrix is
retained. The actual affine sources `F,G` are enclosed using their authenticated
uniform component bounds. Their point values and signs are not newly evaluated.

The result does not change the inherited whole-cone error certificate
`X(U-W1)<3.541e25`. It adds local information about the same candidate. A sign
of `W1-W0` is not a sign of `U-W1` or a physical force.

## Numerical result at the selected point

For the scalar mean, exact outward enclosures prove

```
1.596e21 < C_w:0     < 6.733e21,
8.355e22 < (C_w:0)_x < 1.058e23,
4.147e17 < (C_w:0)_r < 1.750e18.
```

The old scalar-mean forcing defect is negative at this endpoint, while this
accumulated correction is positive. A one-point residual sign does not
determine the sign of the causal response.

The calculation determines correction signs for 16 value components, 16 time
derivatives and 17 normal derivatives. It tightens 53 of the 69 best target
point boxes; all 69 remain inside the preceding boxes. Scalar-mean target
widths shrink by factors exceeding 1.125 for value/normal derivative and 1.330
for time derivative. These are comparisons of interval widths, not measured
improvements in unknown true errors.

The complete 23-component correction at this point has H-norm upper bounds
`7.646e21` for value, `1.202e23` for time derivative and `1.061e21` for normal
derivative. The first two upper certificates are over 330 and 178 times
smaller, respectively, than the prior local analytic caps. Remaining initial
data uncertainty limits the improvement in the target boxes.

## Construction

With `La Y=A Y_x+B Y_r+D Y`, the complete physical source is

```
r0=(La-L0)Bs+La(W0-Bs)+F,
C=K r0+L(M PN C).
```

The first identity avoids subtracting two large free-wave source terms. The
interior evaluator computes `W0-Bs` directly, including the reflected frozen
source and the actual shell lift. Displacement primitives use whole-panel
nominal ranges, positive quadrature and explicit smoothing error bounds.
Velocity primitives use exact convolution of piecewise linear hinges.
Every finite shell-history series includes a factorial tail.

Put `a=9/8`, `alpha=s+y`, `theta=(y-s)/(s+y)`. The active physical triangle
has `1<alpha<a`, `-1<theta<1` and measure
`ds dy=(alpha/2) d alpha d theta`. For its weighted moments,

```
Qj=integral_Omega (a-s-y)^j/j! r0(s,y) ds dy,
PN C(a,0)=PN sum_j>=0 (-M)^j Qj.
```

This folds the nested Robin history into one family of weighted integrals.
For `p=integral_0^a r0(s,a-s) ds`, the shell identities give

```
PN C_x=PN p-M PN C,      PN C_r=M PN C,
PD C=PD C_x=0,          PD C_r=PD p.
```

The full physical coupling is formed before output projection or reflection.
No coefficient matrix is commuted through the parity operator. The inherited
zero wedge `s+y<=1` is a support premise separate from magnitude bounds.

## Scope and verification

The mathematical and source-projection assumptions are inherited from the
registered reduced system. This remains a finite update; the forcing in `W1`
is evaluated at `W0`. It does not solve the self-consistent second-order bulk
problem or establish a full Einstein/constraint residual, a quantitative C2
error bound, an exact nonlinear family, an invariant observable, or a
higher-dimensional origin of the Big Bang.

The frozen public manuscript and private editorial history are untouched.
No journal submission or new Zenodo version is part of this checkpoint.

## Reproduction

From the repository root, using only the Python standard library:

```sh
python -I -S -B research/HDBLAST_CAUSAL_BULK_QUADRATURE_20261030/checks/check_bulk_quadrature.py --repo-root . > /tmp/hdblast-quadrature.json
python -I -S -B research/HDBLAST_CAUSAL_BULK_QUADRATURE_20261030/checks/check_point_transfer.py --repo-root . --quadrature-report /tmp/hdblast-quadrature.json > /tmp/hdblast-transfer.json
python -I -S -B research/HDBLAST_CAUSAL_BULK_QUADRATURE_20261030/checks/check_release_summary.py --repo-root .
```

[QUADRATURE.json](reference/QUADRATURE.json) contains the exact rational
enclosures and integrated source moments; [TRANSFER.json](reference/TRANSFER.json)
contains the resulting candidate and target point boxes.
[SUMMARY.json](reference/SUMMARY.json) checks the readable statements above.
The [derivation](theory/QUADRATURE_NOTE.md),
[target transfer](theory/POINT_TRANSFER_NOTE.md), and
[independent formula review](theory/FOLDED_GEOMETRY_REVIEW.md) state the
assumptions and error bounds.

The published workflow reproduces the fine and coarse quadratures, target
transfer, inherited-point regressions, release statements and independent
geometry controls in normal and optimized execution modes, comparing their
reports byte for byte. The fine quadrature uses 544 physical source panels,
34 edge panels and degree-three Robin history plus its infinite tail.
The independent manufactured fixtures comprise 27 exact groups and seven
wrong-formula controls. Point-overlap and interval-inclusion regressions are
diagnostics; the enclosure proofs use whole-panel ranges and explicit tails.
