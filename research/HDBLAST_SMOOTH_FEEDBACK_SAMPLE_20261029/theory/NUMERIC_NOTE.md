# One complete point evaluation of the tighter smooth nominal response

This calculation evaluates `V_s=B0_s+B1_s` in all 23 modal slots at
`(x,r)=(1/1000,1499/1000)`. It uses exactly the tighter integral-defined
candidate from `HDBLAST_SMOOTH_CANDIDATES_20261026`: the expression-defined
nominal displacement `u0` convolved with the normalized compact C2 kernel of
width `delta_u=2^-30`, and the P1 velocity convolved with the same kernel shape
of width `delta_v=1/2048`; its radial profile is `u_s'`. The candidate has only
one constant-matrix feedback, and omits the affine bulk and shell sources as
in the prior construction. The prior full C1 error certificate and its source,
background, exact-formal-completion and subsidiary assumptions are retained.
This does not construct an exact nonlinear spacetime or compute a C2 error.

At this point `a=r+x=3/2`, `b=r-x=749/500`, and `0<x<r`. There is no shell
lift and no reflected branch in the local first-feedback integral. The entire
interval `[b-delta_u,a+delta_u]` is inside a single analytic branch of the
nominal initial data and inside the retained footprint. The local formulas in
the companion derivation reduce first feedback to endpoint values and first
radial derivatives of `u_s`, the single local integral `Iu=int_b^a u_s`, and
exact moments of the smoothed P1 velocity. This avoids subtracting two large
global initial-data primitives.

## Interval construction

The executable authenticates the published initial producer, initial report,
response profiles and response report before use. The initial producer's
Context authenticates all seven upstream roots before scientific imports and
again afterward. It evaluates outward nominal interval boxes at both endpoints
and throughout each local integration panel.

For each component i, let M2_i be its inherited absolute bound on `u0,rr` over
the complete footprint. Since the original nominal displacement is C1 across
its finitely many branch joins, its derivative is Lipschitz; affine extension
past the last point preserves that bound. The exact normalized-kernel moments
give

```
|u_s-u0|_i <= M2_i delta_u^2/18,
|u_s'-u0'|_i <= 35 M2_i delta_u/128.
```

These remainders enlarge the raw endpoint boxes. To enclose Iu, the script
partitions `[b,a]` into 32 panels, evaluates `u0` on each entire panel using
interval arithmetic, multiplies by its exact positive width, sums the boxes,
and adds `(a-b) M2_i delta_u^2/18`. Thus it uses interval integration, not an
uncertified point quadrature. There is no unaccounted quadrature remainder.
The velocity and its first three primitives are exact rational compact-kernel
hinge sums. The source matrices and all signs are the frozen 23-mode ones.

The exact output contains separate B0, B1 and total value/x-derivative/radial-
derivative boxes for every mode. Both TT candidate slots remain exactly zero;
this does not say their actual target is zero. The complete numerical table
is a point evaluation and does not supply a continuous candidate enclosure on
the full domain.

## Conditional comparison with the exact formal response

The program independently recomputes the previously published C1 transfer
ledger from the same pinned reports, including all constant-matrix block
Frobenius caps in the positive Gram norm. Its exact rational upper bound is
identical to the old tighter-candidate checker and is below `E=2.184e28`.
This verifies the ledger arithmetic; it does not re-prove the underlying
actual-field error certificate or the analytic smoothing theorem.

The compatible zero-wedge estimate gives `|U-V_s|_H <= (x+r-1) E = E/2` at
this point; derivative errors are at most E. An individual mode is bounded by
`H_i^-1/2` times these norm errors, using a rigorous upward rational square
root. For the primary scalar mean w:0 the Gram weight is one. Exact arithmetic
then gives the readable strict bounds

```
4.17e28 < U_w:0 < 6.36e28,
-5.288e31 < partial_r U_w:0 < -5.283e31.
```

The x-derivative interval contains zero. The script emits the stronger rational
interval endpoints. The value sign and radial decrease are conditional signs
of one formal modal coefficient in the chosen normalization, not signs of a
physical force, a total scalar observable or an Einstein residual. Such signs
may already be implied by earlier coarser candidates; no originality or new
physical discovery is asserted. The broad all23 error table does not supersede
earlier sharper specialized TT certificates.

## Reproduce

```
python check_complete_feedback.py --repo-root /path/to/HDblast
python -O check_complete_feedback.py --repo-root /path/to/HDblast
```

Both modes use explicit failure checks and have byte-identical output. The
optional `--panels` setting changes the width of the interval integral but not
the definition of the candidate; the default 32 panels passed every sign and
scope check. This is a one-point interval calculation conditional on retained
mathematics. Feedback truncation is included in the inherited C1 error; a
quantitative C2 error, full pointwise Einstein residual and nonlinear solution
remain unestablished.
