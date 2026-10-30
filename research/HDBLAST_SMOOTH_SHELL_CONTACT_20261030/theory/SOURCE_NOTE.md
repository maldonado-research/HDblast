# Actual affine shell-source enclosure for the active shell sample

The inherited analytic source bounds cover `0<=x<=3, r=0`, hence the new
sample `x=9/8, r=0`. They concern the actual formal affine source `G`, not
the approximate response, a full boundary residual, a force, or an observable.
The frozen factorization is `G=Bc`, with constant rational 23-by-8 matrix `B`.
This check authenticates the two existing source reports before using their
numbers; it does not replay their upstream first-field evolution certificates.

Write positive magnitude caps with overbars. The registered profiles are

```
beta = a*h + sigma1*v/3
alpha = a*h*h/4 - sigma1*h*v/6 - sigma2*v*v/6
pi = p*h*h/4 + sigma2*h*v/2 + sigma3*v*v/2
g0 = 4*sigma2*v*v/3 + 4*sigma1*h*v - 10*a*h*h
     + 8*a*h*tau - 4*a*(4*X*X + 12*L*L - 64*L*d)
c = ell*rho * (-g0, -8*a*V*V, 16*a*d*d,
               2*(-alpha+2*beta*X), pi, -2*beta*V,
               2*(alpha-2*beta*L), -4*beta*d).
```

The triangle inequality, with each factor replaced by its inherited magnitude
cap and `rho` by its interval upper endpoint, gives eight nonnegative rational
profile caps `C`. There is exactly one factor `ell` here; these are source
values, not time derivatives. The separate physical-to-scaled-time factor used
when bounding `G_x` does not enter this value calculation.

For each modal coefficient and for the fixed positive comparison norm,

```
|G_i| <= sum_j |B_ij| C_j,
||G||_H^2 = c^T (B^T H B) c <= sum_ij |(B^T H B)_ij| C_i C_j.
```

The checker rebuilds all 64 entries of `B^T H B` from the signed matrix `B`
and the registered section-space weights. Taking absolute values of the
completed Gram entries is valid for all unresolved signs and correlations;
this is an upper bound, not a statement that they are independently realized.
No signs of the actual source at the sample are claimed.

The exact output establishes the readable bounds

```
||G||_H < 487122,
max_i |G_i| < 1508337,
|G_w:0| < 94279,
max(|G_T:+|, |G_T:-|) < 4.933e-21.
```

The two computed TT magnitude upper bounds are equal; that equality alone
does not determine actual source signs or saturation.
All five homogeneous Dirichlet source rows vanish exactly. The machine report
contains the exact rational cap and symmetric enclosure for every row.

The readable alternative obtained by integrating the prior uniform
`||G_x||_H < 58807538` from the zero collar `x=1` would be `58807538/8` at
this point, which is larger than `487122`. The source value enclosure is
therefore used without a derivative-based improvement.

The first-feedback candidate uses the frozen matrix `M0`, while the actual
formal boundary law uses an interval-enclosed coefficient `M`. Its defect is
`PN partial_r V - M PN V - G = -M0 PN B1 - (M-M0)PN V - G`.
The present source intervals can be subtracted after independently enclosing
both matrix terms. They do not by themselves certify that boundary defect to
vanish. Replacing actual `M` by `M0` without a coefficient-error term would
omit a documented uncertainty.

The reports pinned by the checker are:

- `COMMON_HARMONIC_SHELL_NORMAL.json`: `559f7587feca262cc68990f2ee79a3aeea22df6955e5aafb698f3326fdacd94a`.
- `FULL23_SHELL_NORMAL.json`: `0b4a7fd639deec8bfbe65f7b2427e6b084f18195ef954d0add39471528b6f6fe`.

Run `check_affine_shell_source.py --repo-root /path/to/HDblast` under ordinary
Python and `python -O`; their JSON outputs must match byte for byte.
