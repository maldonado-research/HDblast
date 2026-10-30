# Quantitative harmonic-constraint values for the finite bulk update

This checkpoint assembles the full linear harmonic-covector map and bounds
`C2[W1]=DC[W1]+QC` for the same finite update `W1` evaluated in
[the causal quadrature checkpoint](../HDBLAST_CAUSAL_BULK_QUADRATURE_20261030/).
Here `C2` denotes a second parameter variation of the harmonic constraint;
it does **not** denote a bound on second spacetime derivatives.

At the selected shell point `(x,r)=(9/8,0)`, the direct interval calculation
gives

```
||C2[W1]||_positive sphere L2 < 1.319e24,
sup_n |C2[W1](n)|_positive background frame < 7.801e24.
```

These are absolute bounds in the declared model normalization. They do not
establish relative smallness, an observational tolerance, or an invariant
physical quantity. All eight raw modal intervals contain zero, so no nonzero
constraint sign is proved. The normal component is exactly zero at the shell
by the actual boundary/source identity. The remaining four tangential
components are enclosed; their exact vanishing has not been established.

## What was calculated

The 23 metric/scalar amplitudes map to eight covector amplitudes: three
normal scalars, three time scalars and two scalar-gradient coefficients on
the unit S3. The scalar matter amplitudes and both transverse trace-free
metric channels are annihilated by the linear metric-gauge operator.
The calculation retains the complete registered input and the correct
different Gram weights for metric tensors and covectors.

The direct point calculation uses W1's enclosed first jets and genuine
uniform caps for the same actual, full metric-lowered quadratic gauge source
QC. QC is retained in physical covector units; it is not multiplied by the
wave-forcing scale `ell^2 rho^2`, set to zero, or assigned a point value.
This direct bound does not require assuming that a separate formal target
has zero harmonic constraint.

The shell identity is independently checked using

```
GN=-2 tr(QR)-2 QC_y,
G_K=-ell rho GN,  G_N=2ell rho QR00,  G_J=2ell rho QRJ.
```

Together with the actual homogeneous Robin rows and Dirichlet traces, these
give `DC[W1]_y+QC_y=0`. This verifies the imposed normal harmonic boundary
condition. Dirichlet normal derivatives `L_r,v_r` remain in the tangential
constraint rows; a zero boundary value does not set them to zero.

## Uniform conditional transfer and angular interpretation

For comparison, the existing global first-derivative error certificate can
also be transferred to harmonic values. Conditional on the same fixed QC
and on the formal target having `DC[U]+QC=0`,

```
C2[W1]=DC[W1-U].
```

The resulting positive sphere-L2 bounds are `<5.100e29` throughout the
retained cone and `<4.975e29` at the selected shell point. These are general
error-transfer bounds. Additional direct point information produces the
smaller bound above. The global certificate `X(U-W1)<3.541e25` is unchanged.

The derivative-only scaled operator has the sharp bound `sqrt(7/2)` for
independently bounded time and radial first jets. The full irreducible angular
hull of these covector modes has evaluation kernel
`diag(35,35,(34/3)I_3)`, so a sphere-L2 bound converts to an every-angular-point
positive frame bound by `sqrt(35)`. The sphere-L2 norm itself must not be
called an angular supremum. Both norms use a chosen positive background
frame comparison; neither is a Lorentz-invariant scalar observable.

## Reproduce and inspect

The [direct report](reference/DIRECT.json) contains all eight physical
constraint intervals, source caps, exact boundary refinements and norm sums.
The [transfer report](reference/TRANSFER.json) records actual background
intervals, complete matrices, weights and conditional uniform bounds.
The [summary](reference/SUMMARY.json) independently checks the displayed
numbers and the angular conversion.

From the repository root, using only the Python standard library:

```sh
python -I -S -B research/HDBLAST_HARMONIC_CONSTRAINT_VALUES_20261030/checks/check_direct_constraint.py --repo-root .
python -I -S -B research/HDBLAST_HARMONIC_CONSTRAINT_VALUES_20261030/checks/check_constraint_transfer.py --repo-root .
python -I -S -B research/HDBLAST_HARMONIC_CONSTRAINT_VALUES_20261030/checks/check_constraint_summary.py
```

The workflow reproduces every reference report in normal and optimized
execution modes and compares them byte for byte. Independent exact frame
contractions, angular-polynomial fixtures and deliberately wrong formulas
test signs, frame connections, the half-trace factor, radial orientation,
mixed-component weights, source convention and retained normal derivatives.

The [formula derivation](theory/CONSTRAINT_FORMULA_NOTE.md),
[direct point note](theory/DIRECT_CONSTRAINT_NOTE.md),
[uniform transfer](theory/CONSTRAINT_TRANSFER_NOTE.md), and
[focused literature comparison](literature/LITERATURE_SCOPE.md) state the
assumptions and distinctions from established work. Reviews are internal
mathematical checks, not external specialist or journal review.

## Remaining scientific gap

This quantifies harmonic-gauge constraint **values** for one auxiliary finite
update. It does not prove exact tangential constraint preservation, bound
constraint derivatives, verify Hamiltonian/momentum or full natural Einstein
residuals, provide a quantitative second-derivative error certificate, solve
the complete self-consistent second-order problem, or construct an exact
nonlinear family. Complete differentiated bulk forcing and higher initial
error bounds remain priorities. No evidence establishing a higher-dimensional
origin of the Big Bang is supplied.

The frozen manuscript, existing work and websites are preserved. No new
Zenodo version or journal submission is included.
