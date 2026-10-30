# A real initial-curvature defect caused by affine continuation

This audit concerns the existing smoothed nominal displacement used by the
integral-defined candidates W0 and W1. The retained initial footprint ends at
`R=7/4`. Beyond R, that nominal displacement is continued by its affine
tangent before convolution with the even compact kernel

```
k_delta(s)=35/(32 delta) (1-(s/delta)^2)^3,  |s|<=delta,
delta=2^-30.
```

The continuation is C1, so it supports the previously published C1 smoothing
estimates. It need not be C2: its exterior second derivative is zero, whereas
its interior second derivative at R need not vanish. Every fixed convolution
is sufficiently smooth for the qualitative reconstruction. The issue is
accuracy of its initial second derivative and its limiting behavior as the
kernel width tends to zero, not existence of that derivative.

## Exact half-curvature identity

Let `f` be one nominal component, let `g=f''` on the interior, and assume
`|g'(r)|<=M3` on `[R-delta,R]`. Here that interval lies entirely in the plain
analytic packet branch, well beyond the cutoff join `66/64`. The affine
tangent continuation has exterior second derivative zero. There are no
Dirac contributions in the distributional second derivative because f and
f' match continuously at R. Consequently convolution differentiates as

```
(k_delta*f)''(R)=integral_0^delta k_delta(s) g(R-s) ds.
```

The kernel is positive and even, with exact moments

```
integral_0^delta k_delta(s) ds = 1/2,
integral_0^delta s k_delta(s) ds = 35 delta/256.
```

It follows that

```
(k_delta*f)''(R)=g(R-)/2 + h,
|h|<=35 M3 delta/256.                                      (1)
```

In particular the limit is half the interior curvature, not the full
interior curvature. A piecewise third-derivative bound alone cannot discard
the jump: the C1 function `r_+^2/2` has zero third derivative on either open
side, while convolution of its second derivative at the join is exactly 1/2.
The same issue can occur at interior cutoff joins. For separated joins, a
valid uniform smoothing bound includes a half-jump charge, rather than
assuming that nominal second derivatives are globally Lipschitz.

## Independent actual initial jet

The numeric checker authenticates the existing full23 initial producer and
all seven complete scientific roots before scientific imports. It takes the
actual background value enclosure at R and computes its normalized order-five
jet using the autonomous ODE recurrence. This independently checks the
paired Taylor-difference implementation: it does not differentiate a state
box or treat an unknown Taylor remainder as a polynomial function.

The actual compact packet derivatives through order five come from the
already authenticated degree-six packet-derivative table. Their exact
polynomial recurrence and continuous envelope are checked independently.
The actual exponential jet is enclosed by the retained exact-rational
exponential procedure. Substitution into the registered physical initial
prescription supplies the complete 23-component displacement second jet.

The jet variable is the ordinary physical conformal coordinate z. Its c2
coefficient is normalized by 2!, so the displayed scaled radial derivative
is `U0_rr=2 ell^2 c2`. It is not `U0_yy`; physical normal second derivatives
require the extra warp/connection conversion, including the first-derivative
term. The nominal curvature and its piecewise third-radial-derivative caps
come from the already authenticated initial report in the same convention.

The scalar mean at the endpoint satisfies, with outward intervals,

```
actual U0_rr(R)  approximately [-3.249082574e33,-3.249079963e33],
nominal f_rr(R-) approximately -3.249080947e33.
```

Equation (1), using the registered width and the existing third-derivative
cap, proves the strict fixed-width error enclosure

```
-1.6246e33 < (actual U0_rr - smoothed affine nominal U0_rr)_w:0
           < -1.6245e33.
```

The same strict outer interval also encloses the nonzero zero-width limit
`actual U0_rr(R)-f_rr(R-)/2`. The full positive H norm of the coordinate
initial-curvature error is therefore strictly greater than `1.6245e33`.
This is a genuine error lower bound from enclosed actual data, rather than
a floor in an upper-error budget. The complete report gives every modal
interval, the actual background jets and the exact squared H lower bound.

W1 has the same initial displacement as its smoothed free wave, because its
bulk and boundary correction has zero Cauchy data. Thus the initial finding
applies to the finite update without needing a differentiated bulk source.
It does not supply a spacetime C2 upper bound.

## A precisely defined repair, with separate candidate accounting

A quadratic continuation matching the interior curvature kappa adds
`(kappa/2)*(r-R)_+^2` to the affine extension. At the endpoint, its smoothed
changes relative to the old affine continuation are exactly

```
displacement:       kappa delta^2/36,
first derivative:  35 kappa delta/256,
second derivative: kappa/2.
```

This shows why a substantial curvature correction can coexist with a small
C1 change at a narrow smoothing width. It is a change of candidate, so its
initial/free-wave and propagated C1 charges must be certified separately.
The published affine candidate and its certificates cannot be silently
relabelled. Even with curvature matching, actual-minus-nominal approximation
error remains; matching does not make the physical initial prescription
exact. This audit records the algebraic repair identity and the unchanged
candidate's defect. A separate checker defines and bounds a new continuation.

The defect is at the outer initial endpoint of the chosen footprint. It does
not imply an interior nonlinear obstruction, invalidate earlier C1 bounds,
or refute the physical hypothesis. An interior causal region separated from
that endpoint can remain unaffected by the continuation change. The result
does not establish an Einstein residual, a nonlinear family, an invariant
observable or a higher-dimensional origin of the Big Bang.

## Verification

`check_endpoint_curvature.py` uses exact-rational or outward interval
arithmetic throughout. Its normal/optimized reports are byte-identical. It
includes 35 exact algebra groups and three explicit wrong-formula controls
for the endpoint half-mass, jump omission, and half-first-moment factor.
The packet checks include its derivative recurrence through six and
continuous shifted-polynomial positivity; the cutoff checks distinguish
matching derivatives through four from a possible fifth-derivative jump.
