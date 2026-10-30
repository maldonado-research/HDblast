# Matching the outer second derivative before convolution

Let f be the fixed nominal initial displacement on `[0,R]`, with its existing
odd/even component reflection. The old continuation beyond R is affine,
`f(R)+f'(R)h`. Define the new quadratic continuation by adding
`kappa h^2/2`, where `kappa=f''(R-)` is the exact left analytic value.
Intervals enclose this fixed coefficient; no arbitrary midpoint is chosen.
Let `k_delta(s)=35/(32 delta)(1-(s/delta)^2)^3` on `|s|<=delta`.

The difference of the two smoothed profiles is

```
Delta u_s(r)=(kappa/2) integral k_delta(s)(r-s-R)_+^2 ds.
```

It is zero for `r<=R-delta`. At the endpoint, the half-kernel moments give

```
Delta u_s(R)=kappa delta^2/36,
Delta u_s'(R)=35 kappa delta/256,
Delta u_s''(R)=kappa/2.
```

Positivity and monotonicity give the same value/first-derivative expressions
as component magnitude bounds throughout `[0,R]`. They also give the H-norm
bounds after multiplying by the norm of kappa. Reflection preserves this
comparison norm and does not move any coupling matrix through parity.

On the inner half of the endpoint convolution, the existing one-sided
third-derivative bound M3 gives

```
u_s,affine''(R)=kappa/2 + remainder,
u_s,quadratic''(R)=kappa + remainder,
|remainder_i| <= (35 delta/256) M3_i.
```

The paired actual-minus-nominal endpoint derivative error can therefore be
combined with the remainder before subtracting large interval quantities.
For the scalar mean, the old error is strictly larger than `1.6245e33` in
magnitude and the new error is strictly smaller than `1.820e27`. Exact
cross-multiplication proves the improvement exceeds 800,000 for this single
initial endpoint derivative. The claim is not based on comparing two upper
error certificates.

The new continuation is C2 at R, while its third derivative can jump there.
Its outside second derivative is constant kappa. Since every `|kappa_i|`
lies below the old component M2 cap, its first derivative remains Lipschitz
with the same cap. The old even-kernel bounds `M2 delta^2/18` for displacement
and `35 M2 delta/128` for its derivative therefore remain valid. The separate
transfer calculation rechecks every global majorant needed to apply the same
readable C1 error certificate to the new finite candidate.

For the uniform initial second derivative, only the R curvature jump is
removed from the previous jump-aware error bound. The interior cutoff jumps
remain. The nominal function and its paired actual error inside the footprint
are unchanged before convolution. Reusing those continuous paired caps and
the existing piecewise M3 gives the new conservative initial H bound
`6.219e33`. P1 velocity and its error bound are unchanged.

Define the new reduced finite candidates consistently, using the same fixed
actual coefficients and sources:

```
P=(I-L M PN)^(-1),
W0_quadratic=P[Bs_quadratic+K L0 Bs_quadratic+LG],
W1_quadratic=P[Bs_quadratic+K(La W0_quadratic+F)+LG].
```

The old candidates are not redefined. Since all differences start in the
outer initial strip and the operators are causal, they vanish whenever
`x+r<R-delta`. The selected point `9/8` and retained shell times through
`5/4` precede this strip's arrival. All inherited point enclosures there
describe exactly the same fields for the new candidate, not merely nearby
ones. The full actual affine shell boundary law is retained by construction.

This repairs an amplitude-profile derivative artifact at the initial outer
endpoint. It is not a statement about a spacetime curvature invariant,
uniform evolved C2 accuracy, physical constraints, nonlinear realization or
the origin of the Big Bang.
