# A conservative complete-candidate enclosure at one unreflected point

This is an independent coarse enclosure for the tighter integral-defined smooth
candidate from `HDBLAST_SMOOTH_CANDIDATES_20261026`. It combines a previously
published validated free-wave point calculation with elementary bounds for its
entire first feedback. The result is conditional on the inherited formal
full23 response, background, first-field and C1 error certificates. It is not
a quantitative C2 error, Einstein-equation residual, nonlinear existence theorem
or physical scalar observable.

## The point and the candidate

Use `(x*,r*)=(1/1000,1499/1000)`. The tighter candidate is
`V_s=B0_s+B1_s`, with

```
u_s=k_(2^-30)*u0,
v_s=k_(1/2048)*vP1,
B1_s=K(A B0_s,x+B B0_s,r+D B0_s)+L(M PN B0_s|shell).
```

Here `u0` is the exact expression-defined nominal initial displacement, `vP1`
is the old continuous piecewise-linear initial x-velocity, and both have the
same parity and endpoint continuation as in the published construction. The
candidate omits the affine bulk/shell sources, whose effects are in the
inherited error ledger.

For every point `(s,z)` of the backward cone of `(x*,r*)`, its initial free-wave
arguments `z-s` and `z+s` lie in
`[r*-x*,r*+x*]=[749/500,3/2]`. This interval, enlarged by either convolution
width, lies inside the retained initial footprint, away from its cutoff and
outer endpoint. The cone never reaches `r=0`; no source reflection or shell
lift contributes here. These assertions are checked with exact fractions.

## Coarse feedback bound

Let `U0,U1` be the inherited positive H-norm caps on `u0,u0'`, and let `Vcap`
be an outward bound on the largest H norm of the exact P1 velocity node vectors.
Convexity bounds the norm of the P1 interpolant by its nodes. The smoothing
kernels are positive and normalized, so the same caps bound `u_s,u_s',v_s`
throughout this local footprint. The reflected free-wave formula, used here
only in its unreflected branch, gives

```
|B0_s|_H <= U0+x* Vcap,
|B0_s,x|_H, |B0_s,r|_H <= U1+Vcap.
Fcap=(nA+nB)(U1+Vcap)+nD(U0+x* Vcap).
```

The constants `nA,nB,nD` are the frozen comparison-norm operator bounds. The
checker independently rebuilds the 23 Gram weights, verifies exact off-block
zeros, and checks each operator's weighted block Frobenius bound against those
constants. The full coupled source is bounded before any reflection; this cone
has none.

For a source with norm at most Fcap, d'Alembert's integral and its endpoint
derivatives yield

```
|B1_s(x*,r*)|_H <= (x*)^2 Fcap/2,
|B1_s,x(x*,r*)|_H, |B1_s,r(x*,r*)|_H <= x* Fcap.
```

Exact arithmetic gives approximately

```
Fcap = 2.6754727599e29,
B1 value cap = 1.3377363799e23,
B1 derivative cap = 2.6754727599e26.
```

These are symmetric upper bounds, not signed quadrature evaluations of B1.
They bound the complete first feedback of this specified candidate. The checker
also directly integrates a nonzero affine synthetic source over the cone and
tests its independent derivative formulas, including the half factor and radial
sign.

## Transfer to one actual formal coefficient

The pinned published interval checker encloses the `w:0` free-wave value and
its first derivatives. `w:0` is the registered scalar degree-zero coefficient,
with unit H Gram weight, so any H-norm error bounds that component directly.
Adding the symmetric B1 caps encloses the complete candidate's selected value
and first derivatives at this point.

Recomputing the tighter candidate's inherited C1 transfer gives
`X(U-V_s)<=E<2.184e28`. Both fields vanish on the zero wedge `x+r<=1`.
At fixed x*, the segment from `r0=1-x*` to r* lies in D+ and has length 1/2.
Integrating the H-bounded radial error derivative therefore gives local value
error at most E/2, while the derivative errors are at most E. Combining them
with the candidate intervals proves, strictly,

```
U_w:0(x*,r*) > 4.17e28,
partial_r U_w:0(x*,r*) < -5.283e31.
```

The full exact rational intervals are the primary computed result. Approximate
displays for the actual formal coefficient are

```
U_w:0       in [4.174089534e28, 6.357771581e28],
U_w:0,r     in [-5.287815946e31,-5.283394376e31].
```

The x-derivative interval includes zero. The displayed scalar sign is a
registered second parameter coefficient, not a gauge-invariant field value,
finite-amplitude scalar field, force or cosmological signal. Previous work
already certified other evolved scalar-mean signs on different whole cells;
this point result must not be described as the first such sign or a physical
discovery.

## Reproduction and interpretation

Run `check_coarse_feedback.py --repo-root /path/to/HDblast`. An optional
`--smooth-capsule` lets local audit use the exact published smoother files from
a separate checkout; its B0 checker must match the immutable SHA256. Default
location is the published `research/HDBLAST_SMOOTH_CANDIDATES_20261026` folder.
The checker runs that authenticated B0 interval program, which authenticates
its upstream interval inputs, and imports no feedback project code. It
reconstructs the error transfer independently from the pinned earlier reports.
Normal and optimized runs must agree byte for byte.

This completes a conservative one-component candidate enclosure and a
conditional transfer to the formal response. It does not numerically evaluate
the full response, its complete inhomogeneous boundary defect, differentiated
bulk sources, constraints, pointwise Einstein residual or nonlinear remainder.
