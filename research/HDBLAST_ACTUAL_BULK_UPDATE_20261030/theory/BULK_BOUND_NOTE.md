# Conditional majorant for a boundary-preserving actual-bulk update

This calculation bounds a new **integral-defined** candidate. It is not a
quadrature evaluation of the actual bulk correction or affine bulk source.
The fixed actual coefficient functions, actual Robin matrix and actual affine
sources are the same throughout the construction. Their earlier conditional
enclosures and the retained formal-response premises remain in force.

## Construction and operator split

Use the registered norm

```
X(V)=max(sup ||V||_H/3, sup ||V_x||_H, sup ||V_r||_H)
```

on the retained cone. Let `Bs` be the smoothed nominal free wave, `Ba` the
actual-initial-data free wave, `La Y=A Y_x+B Y_r+D Y`, `Bdry Y=L(M PN Y)`,
and `P=(I-Bdry)^(-1)`. `Bdry` is an operator, distinct from the bulk matrix
`B`. No bulk or boundary matrices are commuted through each other or through
the physical-source reflection.

The previous reconstruction is

```
W0=P[Bs+K(L0 Bs)+L G].
```

Define the actual-bulk update map and iterates by

```
Psi(Y)=P[Bs+K(La Y+F)+L G],
W(n+1)=Psi(Wn).
```

Thus `W1` retains the actual affine boundary law and includes the actual bulk
coefficients and affine bulk forcing, acting on `W0` once. Its forcing is
`La W0+F`, rather than `La W1+F`. It is therefore not yet a self-consistent
full bulk solution. The construction includes `F` as the same actual function;
it does not replace `F` by zero or by a midpoint.

## Recomputed norm caps

The inherited zero-wedge estimates give a whole-cone field-value factor `3/4`
and a shell-trace factor `1/4`. The checker verifies every raw actual-minus-
frozen coefficient interval and each weighted block Frobenius bound, then sets

```
k=(3/2)(na+nb+(3/4)nd),
b=nm/4,
q=k+b,
beta=k/(1-b).
```

Here `na,nb,nd,nm` are actual operator upper bounds, including coefficient
uncertainty. The resulting numbers are approximately

```
k     =0.000321529951748768,
b     =0.0000684876118606234,
q     =0.000390017563609391,
beta  =0.000321551974075561.
```

Since `b<1`, the boundary inverse has `X` operator norm at most `1/(1-b)`.
Hence `Psi` has Lipschitz bound `beta<1`. This inversion is a causal Volterra
or Neumann construction; it requires no inverse of the matrix `M` and no
spectral positivity claim.

For a separate conservative calculation, retain the older full cap
`qR=3409/8000000`, assign its unused slack to the bulk, and use

```
kR=qR-b >= k,
betaR=kR/(1-b) <1.
```

The exact and retained splits are computed separately. The checker does not
mix the smaller `beta` with the larger denominator `1-qR` without a proof.

## Residual and initial-data contributions

The existing full residual after the wave propagator is

```
r0=(La-L0)Bs+La(W0-Bs)+F,
X(K r0) <= R.
```

The checker reconstructs the initial smoothing and velocity interpolation,
refined continuous initial cover, whole-cone first feedback, global Robin
correction, coefficient mismatch and lower-order feedback residual. It
reproduces the earlier value

```
R =1.667200340838006...e26
```

as an `X` response cap. This number is not a pointwise bound on `r0` or `F`.
The affine source estimate remains in its earlier units: `SOURCE` bounds
`X(KF+LG)`, and `X(KF)<=SOURCE+sup||G||_H`. The separately inherited shell
source cap is retained. These units are checked explicitly.

The first correction satisfies

```
W1-W0=P K r0,
D0 := R/(1-b) >= X(W1-W0).
```

Numerically, `D0` is approximately `1.667314531228471e26`. The actual target
`U` has different initial data from all `Wn`. Its free-wave initial mismatch
has recomputed bound

```
X(Ba-Bs)<=e_initial=3.532863599878737...e25.
```

This term must not be dropped merely because the bulk update is contractive.
The mismatch in initial data persists through all iterations.

## One update and the subsequent majorant sequence

The error identity gives the recurrence

```
X(U-W(n+1)) <= e_initial/(1-b)+beta X(U-Wn).
```

Equivalently, let `Us` be the reduced integral problem's fixed point with
smoothed nominal initial data and the same actual coefficients and sources.
Contraction bounds give

```
X(Us-Wn) <= beta^n D0/(1-beta),
X(U-Us)  <= e_initial/(1-q),
En := e_initial/(1-q)+beta^n D0/(1-beta)
    = (e_initial+beta^n R)/(1-q).
```

These statements concern the reduced integral problem under the inherited
premises. They do not prove that `Us` is an exact nonlinear spacetime or that
all Einstein subsidiary constraints have been solved.

The same formulas with `qR,betaR` give conservative majorants. For `n=0`, both
versions reproduce the previous direct certificate exactly. For `n=1`,

```
E1_actual_split  =3.539605023613945...e25,
E1_retained_split=3.540335165342622...e25 <3.541e25.
```

The previous retained bound was approximately `2.021348047762733e26`.
Its ratio to the new readable one-update bound exceeds `5.7`. This is an
upper certificate for a different integral-defined candidate, not a measured
decrease in the true error or a numerical integration result.

| Actual-bulk update count | Retained-split error majorant, approximately |
| --- | ---: |
| 0 | `2.021348047763e26` |
| 1 | `3.540335165343e25` |
| 2 | `3.534371816786e25` |
| 3 | `3.534369683924e25` |
| 4 | `3.534369683161e25` |

The limiting initial-data majorant is approximately `3.534369683159973e25`.
It is the remaining term of this estimate, **not a lower bound on the unknown
true error**. Additional iterations have diminishing returns in this bound.

## What norm information alone does not improve

There is an exact arithmetic identity for both splits:

```
D0+E1=E0.
```

If `W1` is enclosed only by widening a known `W0` sample by the norm bound
`D0`, and then the target is enclosed by widening again by `E1`, the total
norm radius is exactly the old one. Therefore this norm-only calculation does
not by itself give narrower target point intervals. Resolving the actual bulk
correction by validated quadrature or another sharper enclosure is required
for new point information. The machine report states this limitation rather
than presenting the old target intervals as a new numerical evolution.
The separate [local causal estimate](LOCAL_NOTE.md) provides such additional
information and gives narrower target point intervals. It does not alter
the global identity or quadrature-evaluate the correction.

## Reproduction and limits

Run `check_bulk_update_bound.py --repo-root PATH` normally and with `-O`;
their JSON reports must be byte-identical. Optional `--initial-capsule`,
`--robin-capsule` and `--source-report` arguments select separately located
pinned inputs. The checker uses standard-library JSON and exact fractions;
it imports no scientific interval evaluator. It independently checks the
continuous refinement cover and norm aggregation, raw coefficient boxes,
smoothing jump norm, initial mismatch, residual units and old certificate,
then verifies the update recurrence and norm-enclosure identity in exact
rational arithmetic.

The analytic operator estimates and prior first-field/background/source
certificates remain inherited conditional assumptions. Qualitative regularity
needed to read the integral identities as classical wave equations must be
verified separately; this norm calculation is not a quantitative C2 estimate.
No pointwise Einstein/constraint residual, nonlinear family, invariant physical
observable or higher-dimensional Big Bang evidence is established here.
