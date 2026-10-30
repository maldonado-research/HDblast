# Independent numerical and normalization audit

The numerical transfer code and reports in `checks/check_constraint_transfer.py`
were reviewed independently of their producer. Normal/optimized reports are
byte-identical. All actual source and coefficient assumptions remain inherited
conditional premises.

## Complete first-order map and norm transfer

The three normal/time and two gradient covector modes agree with the direct
five-frame contraction in `check_constraint_formulas.py`. The output norm has
weights `S_l,S_l,lambda*S_l`, while the metric mixed entries in the input norm
have their factor-two symmetric-tensor multiplicity. Scalar-field `w` and both
TT slots are absent from every constraint map.

The derivative bound `sqrt(7/2)` follows from the exact joint weighted Gram
identity. The shell bound correctly retains the Dirichlet normal derivatives
`L_r,v_r`; it uses `Rx PN` and `Rr PD`, with operator bounds
`sqrt(5/4)` and `sqrt(1/2)`. A zero Dirichlet value/time trace does not set
its normal derivative to zero.

I independently reconstructed all six weighted Frobenius block squares from
the reported interval matrices and positive input/output weights. Each lies
below the recorded block cap squared. I recomputed the complete scaled and
physical-frame bounds from the source report, rather than relying on its
readable thresholds. The resulting upper bounds are approximately

```
global scaled: 6.624862913317297e25,
global frame:  5.099128555799713e29,
shell scaled:  6.462826265751241e25,
shell frame:   4.974283673104305e29.
```

The physical-frame result divides the scaled result by `ell*rho`, not by
`ell` alone. Certified rho lower endpoints are positive, approximately
129.92147 globally and 129.92476 at the shell. Actual background ranges
and hyperbolic enclosures cover the declared containing rectangle and point.
The exact symbolic normal cancellation uses the shared primitive Robin
coefficients; it is not inferred by canceling independent interval entries.

The error transfer assumes the same quadratic gauge source for U and W1 and
the retained formal subsidiary completion `DC[U]+QC=0`. It gives a covector
value upper certificate. It does not itself establish that formal completion,
source derivatives, quantitative C2 control, a full Einstein residual, or a
physical observable. The absolute sizes do not establish relative smallness.

## Source normalization for a direct point audit

The authenticated geometric source report states that `QC` is the
actual-metric-lowered quadratic gauge Hessian. Its `QC/y`, `QC/0` and `QC/a`
component caps therefore already have the physical parallel-frame covector
normalization. The source inventory explicitly states that `QC_0` is lower,
not raised. The metric/scalar wave-forcing factors
`-2ell^2rho^2` and `ell^2rho^2` apply to QEred and Qphi, respectively; they must
not be applied to QC in a direct constraint sum.

A direct sum must first convert `ell*rho*DC[W1]` to `DC[W1]`, then add these
physical QC component intervals. The analytic shell normal cancellation can
be proved from the exact actual source-packing identity
`G_K+G_N-3G_J=2ell*rho*QC_y`. A finite interval sum merely containing zero
does not prove an exact cancellation. Tangential intervals containing zero
also supply no nonzero lower bound or equality theorem.

The all-angular-point conversion factor `sqrt(35)` follows from the full
irreducible-hull addition kernel. It converts the normalized angular L2
covector norm to the positive parallel-frame value norm. It does not create
a Lorentz-invariant observable.

## Completed direct-point audit

I reviewed the final direct checker in `checks/check_direct_constraint.py`
and independently recomputed its eight complete physical `DC[W1]+QC` interval
sums and their raw/refined weighted norm squares. Normal/optimized reports
are byte-identical. All eight raw modal intervals include zero. The direct
checker retains nontrivial Dirichlet normal derivatives and rejects an
explicit omission of `L_r`; its analytic normal refinement uses authenticated
normal-gauge trace joins at all three scalar degrees and the retained exact
boundary law, not cancellation of independently widened intervals.

The raw and boundary-refined positive angular L2 norm upper values are
approximately `1.3186232124e24` and `1.3185609304e24`, respectively. Their
squared endpoint sums are both strictly below `(1.319e24)^2`. There are
three exact normal refinements, while the five tangential modal intervals
are left unchanged. The physical QC source enters once with its correct
lower-covector normalization; it is not multiplied by a unit-principal wave
forcing factor.

This local direct upper certificate is available without assuming that a
separate formal target U has zero constraint. It shares the authenticated
actual background, first-field/source-projection and candidate-enclosure
premises. Containment of zero does not prove the four tangential physical
components vanish or produce a nonzero lower bound. The exact normal
cancellation is a boundary component identity, not full physical constraint
verification. Comparing this direct point cap to the transfer cap at the
same point compares upper certificates; it does not measure the unknown
true defect or show that it is small relative to a physical scale.

The direct point report SHA256 reviewed here is
`08c2849ea32071821a1b4cd5550c20ca3031e2ed538148e22d9a59f17379c872`.
No source-unit, radial orientation, missing normal derivative, norm or claim
inconsistency was found in either final numerical route.
