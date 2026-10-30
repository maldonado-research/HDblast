# Actual initial-curvature mismatch at the two inner cutoff joins

30 October 2026. Conditional coordinate calculation for the **unchanged**
`W1_quadratic` (`W1_dagger`) reconstruction of PR100/PR102. The actual target,
actual first field, original bulk and shell sources, nominal initial profile,
positive Gram norm, and smoothing widths are preserved. This calculation does
not change the finite update or re-establish a smaller evolved upper bound.

## What the signed calculation establishes

For `e=U-W1_quadratic`, at the initial time `x=0`:

| Initial radial point | Actual-minus-candidate scalar mean curvature `e_rr,w:0` | Other sign-excluding rows | Complete positive `H` norm lower bound |
| --- | --- | --- | --- |
| `65/64` | `0.6114 < e_rr,w:0 < 0.6115` | `J:0,J:2,J:4,w:2,w:4` | `||e_rr||_H > 0.6942` |
| `66/64=33/32` | `-5.195e20 < e_rr,w:0 < -5.193e20` | `w:2,w:4` | `||e_rr||_H > 5.896e20` |

These are bounds on the **actual target-minus-fixed-candidate error**, not just
nominal curvature jumps or nonzero upper allowances. All 23 rows are retained in
the JSON report. At the first join six rows have determined signs; at the second
join three do. The other row intervals include zero and establish no nonzero
error there. Exact-zero rows and overlapping intervals are distinguished in the
machine-readable values.

At `65/64`, every one-sided nominal curvature trace is exactly zero. The
nonzero actual error there diagnoses the discarded analytic packet tail,
including the six signed rows above. At `66/64`, the signed nominal curvature
jump is nonzero in the scalar rows; averaging its two sides contributes to the
remaining actual curvature mismatch. A nonzero nominal jump alone would not
prove the error statements in the table: the actual target is separately
bounded and subtracted.

The fixed-width and zero-width-limit error intervals both exclude zero in the
listed rows. Consequently merely taking a narrower even normalized kernel,
while retaining the same nominal initial profile, does not make these pointwise
initial curvatures converge to the actual target curvatures. This is an
approximation property of this reconstruction. It is not a model obstruction,
a physical force, or a statement that an exact spacetime fails to exist.

## Proof of the local convolution enclosure

The candidate displacement is the convolution of its retained nominal profile
with the even, nonnegative normalized kernel

```
k_delta(s)=35/(32 delta) (1-(s/delta)^2)^3, |s|<=delta,
delta=2^-30,
integral_0^delta k_delta(s) ds=1/2,
integral_0^delta s k_delta(s) ds=35 delta/256.
```

The two joins are separated by `1/64`, and each kernel footprint contains just
its one join and lies strictly inside the unmodified outer-extension region.
The endpoint identities of the packet cutoff polynomial through derivative
order four give matching displacement and first derivative in the registered
initial formulas. Thus the nominal displacement is `C1`; its distributional
second derivative is its piecewise classical second derivative. There is no
Dirac delta term in the second derivative.

Let `kappa_minus,kappa_plus` be the exact nominal second-derivative traces at a
join `rc`, and let `M3_minus,M3_plus` bound the absolute nominal third derivative
on the two closed side intervals of width `delta`, using the forced analytic
branches on each side. Then

```
u_s,rr(rc)=(kappa_minus+kappa_plus)/2+R,
|R_i|<=35delta/256*(M3_minus_i+M3_plus_i).
```

This follows by applying the fundamental theorem of calculus on each analytic
side and integrating its distance-to-trace error against the positive kernel.
The local intervals include the respective one-sided traces. They do not
assume that the third derivative agrees across the join.

The checker evaluates branches `zero/cutoff` at `65/64` and `cutoff/plain` at
`66/64` exactly, without floating-point epsilon stencils. Nominal order-six
Taylor jets provide each side's local third derivative; ordinary scaled radial
curvature is `2 ell^2 c[2]` and nominal third derivative is `6 ell^3 c[3]`, with
`ell=10^-6`. The factors are ordinary-derivative factorials, not normalized
Taylor coefficients.

## Actual-point evaluation and cancellation

The authenticated retained order-five background representation bounds each
background derivative separately. No unknown remainder function is
differentiated. The actual analytic packet derivatives use their signed bump
polynomials through order five. This signed point evaluation is useful at
`65/64`: the retained broad tail cap in the paired packet route has symmetric
derivative bounds and therefore cannot by itself determine these signs.
Evaluating the same actual analytic function with a sharper signed enclosure
does not change the target.

The direct actual curvature interval is intersected with the retained paired
actual interval. To preserve cancellation where it matters, the error is also
computed as

```
2ell^2 * paired.nominal.c[2]
+ sum_channels(2ell^2 * paired.difference[channel].c[2])
- (kappa_minus+kappa_plus)/2 + [-Rcap,+Rcap].
```

This interval is intersected with the signed direct-actual-minus-average
interval. Both enclose the same actual error. The resulting fixed-width error
and the analogous zero-width limit are recorded separately for every row.

For each error interval `E_i`, define `m_i=0` if it includes zero and otherwise
`m_i=min(abs(E_i.lower),abs(E_i.upper))`. The retained diagonal positive Gram
weights then give the rigorous pointwise lower bound

```
||e_rr(rc)||_H^2 >= sum_i weight_i*m_i^2.
```

An exact rational dyadic square-root lower bound is computed with integer
arithmetic. Correlations among interval coordinates do not invalidate this
inequality; each component lower bound holds simultaneously for the actual
vector.

## Reproduction and limits

Run with the original authenticated input repository and adjacent pinned core:

```
python -I -S -B checks/check_inner_join_error.py --repo-root /path/to/HDblast
python -I -S -B -O checks/check_inner_join_error.py --repo-root /path/to/HDblast
```

The normal and optimized JSON reports are identical. The checker uses explicit
runtime checks rather than Python `assert` statements. Its six corrupted
formula modes reject wrong half mass, half first moment, curvature factorial,
one-sided branch choice, cutoff derivative matching, and Gram weight. The
complete retained source closure is authenticated before scientific imports and
on final readback.

These initial-point lower bounds do not replace the PR102 conditional global
second-error upper bound, do not explain the full size of that conservative
majorant, and do not establish relative physical accuracy. A separate continuous
refinement can reduce interval overestimation without altering this candidate.
A different cutoff or nominal initial profile would define a new candidate and
would require a new error and transfer ledger. No nonlinear family,
cosmological observable, experimental evidence, journal acceptance, or Big
Bang-origin confirmation is claimed.
