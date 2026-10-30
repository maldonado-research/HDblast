# Higher initial-error bounds and an affine-extension curvature defect

The candidate analyzed here is the fixed W1 used in the earlier bulk-update,
causal-quadrature and harmonic-value checkpoints. Its initial displacement is
the convolution `u_s=k_du*u_nom`, with `du=2^-30`, and its velocity is
`v_s=k_dv*P1v`, with `dv=1/2048`. Neither data prescription is changed here.
The actual target has the registered physical harmonic-completed initial
data. This calculation bounds actual-minus-candidate **initial-data** second
radial derivatives and velocity gradients. It is not a full spacetime C2
error certificate.

All derivatives in the numerical statements below are with respect to the
scaled initial coordinate r. The H norm is the same positive full23 angular
comparison norm as the preceding response certificates, including the scalar
and both transverse tensor slots.

## Order-five paired algebra

The retained initial producer uses four-input jets to bound initial values,
their first radial derivative and velocity. The scalar initial displacement
contains three packet derivatives. Therefore its second radial derivative
requires packet jets through order five. Differentiating a previously reported
error magnitude would not supply these jets.

The local helper extends the same telescoping paired algebra to order five.
No scientific input is replaced. The authenticated bump table already contains
derivative polynomials through order six. If z=ell*r and a background function
has its Taylor polynomial through degree four and a uniform normalized fifth
derivative box F5, then the normalized nth jet has the Taylor remainder

```
binom(5,n) F5 z^(5-n),  0<=n<=5.
```

For n<=4 this is Taylor's theorem applied separately to the nth derivative;
for n=5 it is the normalized fifth-derivative enclosure itself. This does
not differentiate an unknown remainder coefficient as if it were constant.

The exponential-versus-rational-kernel difference formula is extended to the
same order and checked against the actual bump-derivative oracle. The three
paired channels remain the inherited shell-Taylor coefficient difference,
rational packet/cutoff and exponential approximation, and validated fifth-ODE
derivative remainder/outward arithmetic. Ordinary radial derivatives are
`u_rr=2 ell^2 jet[2]` and `v_r=ell jet[1]`.

The complete closures are authenticated before project imports and afterward.
The checker covers [1,7/4] by 48 closed cells, split at both cutoff joins;
actual and smoothed nominal data vanish through r=1. Channel maxima are taken
componentwise before the positive Gram norm. Consequently the resulting
continuous higher-initial-error bounds cover the full [0,7/4] footprint.

## Jump-aware displacement smoothing

The cutoff polynomial matches packet derivatives only through order four.
Since scalar displacement uses three packet derivatives, the nominal full23
displacement is guaranteed C1 across the cutoff joins, rather than globally
C2. Its piecewise second derivative may jump. The exact nominal displacement
is extended beyond R=7/4 by its affine tangent, so its outer second derivative
is zero even when its interior curvature is nonzero.

Let f denote one nominal displacement component including this affine
extension. Since f and f' are continuous, the smoothed second derivative is
`(k_du*f)''=k_du*f''`, with f'' interpreted piecewise/almost everywhere.
There are no Dirac masses from jumps in f'. The inherited piecewise third
derivative cap M3 controls variation between joins, but cannot by itself
bound the convolution difference across a jump in f''.

Write J_i for the absolute curvature jumps at `65/64`, `66/64` and R. Their
distances exceed `2du`, so a kernel neighborhood crosses at most one join.
For a nonnegative even normalized kernel, the probability of crossing any
single join from a fixed base point is at most one half. The compact kernel
has `E|s|=35du/128`. Therefore, componentwise,

```
|f''(r)-(k_du*f)''(r)|
  <= M3*(35du/128) + (1/2)*max_i J_i.
```

At a join either consistent one-sided representative can be chosen for the
intermediate f'' comparison. The actual initial profile and the convolution
have classical derivatives; the final error bound is independent of this
measure-zero convention. One-sided jump intervals are computed directly from
the same nominal expression, with the appropriate analytic branch forced at
the join. No sampled difference or derivative of a magnitude cap is used.

Adding this convolution error to the paired actual-minus-nominal curvature
error gives the continuous full23 bound

```
sup_[0,7/4] ||u_actual,rr-u_s,rr||_H < 8.064e33.
```

This is an upper certificate, not a claim of small C2 error.

## Velocity-gradient error, including rounded P1 nodes

The stored rational P1 nodes are not assumed to be exact samples of the
analytic nominal velocity. On a cell of width h, the derivative difference
from the exact analytic secant is at most `h*M2v/2`. Nodal radii add
`(radius_left+radius_right)/h`. The authenticated cell interpolation cap is
at least each endpoint radius, so the checker safely uses

```
Dv_cell = h*M2v_cell/2 + 2*P1_cell_error_cap/h.
```

For the proof only, extend the analytic nominal velocity beyond R by its own
tangent line. This auxiliary function has a globally Lipschitz derivative
with the inherited piecewise M2v cap, since the nominal velocity is C1 at its
joins. The candidate still uses the original P1 velocity, extended with its
final P1 slope. Beyond R the difference between these two slopes is bounded
by the same final-cell secant and nodal-radius estimate. No analytic nominal
formula is evaluated outside its certified footprint.

Positivity of convolution then bounds the derivative discrepancy by
`max_cell Dv_cell + M2v*(35dv/128)`. Adding the paired actual-minus-nominal
velocity-gradient channels gives

```
sup_[0,7/4] ||v_actual,r-v_s,r||_H < 9.672e24.
```

## A genuine endpoint lower bound

At R, half of the even kernel sees the affine exterior, whose curvature is
zero. The interior half sees the nominal one-sided curvature N2(R). No other
join is within du. Thus

```
u_s,rr(R) = (1/2) N2(R) + remainder,
|remainder| <= M3*(35du/256).
```

If D2 is the signed paired actual-minus-nominal curvature, the actual error
is enclosed without subtracting two large independent fields:

```
e0_rr(R) = (1/2)N2(R) + D2(R) - remainder.
```

For the scalar mean, exact interval arithmetic proves

```
-1.624543e33 < (e0_rr)_w:0(R) < -1.624538e33,
||(e0_rr)(R)||_H > 1.6245e33.
```

Nine endpoint components are certified nonzero: the N, J and w coefficients
at degrees zero, two and four. The signed scalar-mean error is negative; it
is an initial-curvature mismatch, not a physical force.

With the nominal function and its affine tangent extension fixed, shrinking
the kernel width alone gives the limit
`u_s,delta,rr(R)->N2(R)/2`. The actual-minus-half-nominal scalar curvature
interval is strictly separated from zero by more than `1.6245e33`. Hence
this family of numerical smoothing widths cannot make the uniform initial
second-derivative error vanish. This is a theorem about the specified
numerical extension and mollifier, not an obstruction to the physical model,
an exact nonlinear spacetime family, or the Big Bang hypothesis.

Bulk updates with zero Cauchy data retain the same initial mismatch. A
different exterior extension would define a distinct numerical candidate and
requires its own propagation/error accounting. No extension or frozen
candidate is silently replaced in this calculation.

## Remaining scope

These inputs can enter the differentiated-error gate, but they do not close
it: a complete differentiated bulk-forcing/defect bound is still required.
The lower bound at R is a true initial-data error lower bound; the other
reported majorants remain upper bounds. Neither is a full Einstein,
Hamiltonian/momentum, or quantitative spacetime C2 residual certificate.
The model's background/first-field premises remain conditional. No physical
origin, finite-amplitude nonlinear solution or observational result follows.
