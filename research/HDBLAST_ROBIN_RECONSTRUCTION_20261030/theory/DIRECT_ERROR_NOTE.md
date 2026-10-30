# Direct C1 residual estimate for the Robin reconstruction

This estimate is for the integral-defined boundary-compatible reconstruction
`W`, rather than for its finite numerical Neumann sum. It compares `W` with
the retained formal full23 response `U`, conditional on the inherited actual
initial-data, coefficient, affine-source, zero-wedge and integral-equation
premises. It does not establish a nonlinear spacetime or a numerical C2 error.

The new direct certificate is strictly

```
X(U-W) < 1.735e27,
```

where `X(V)=max(sup|V|_H/3,sup|V_x|_H,sup|V_r|_H)` on the full retained cone.
This upper certificate is more than twelve times smaller than the earlier
tighter-candidate certificate. The comparison concerns proved upper bounds;
it does not measure the unknown approximation error or prove a twelvefold
decrease in that unknown error.

## Exact residual identity

Let `B0` be the derivative-consistent smoothed free wave and let

```
L0 V = A0 V_x + B0mat V_r + D0 V,
La V = A V_x + B V_r + D V,
T V  = K(La V) + L(M PN V).
```

Here `B0mat` denotes the frozen bulk matrix, distinct from the free wave `B0`.
The actual coefficient functions are enclosed by the authenticated actual
boxes. The actual Robin matrix `M` need not equal the frozen `M0`.
The exact formal target satisfies

```
U = Bactual + K(La U + F) + L(M PN U + G).
```

The new reconstruction, with trace solved by the Volterra equation, satisfies

```
W = B0 + K(L0 B0) + L(M PN W + G).
```

Subtracting these integral identities, with `e=U-W`, gives

```
e = (Bactual-B0)
    +K[(La-L0)B0 + La(W-B0) + F] + T e.
```

The affine boundary source `G` cancels because the target and reconstruction
use the same actual source and the same actual Robin coefficient. No boundary
defect term remains. The bulk forcing remains frozen at `L0 B0`, so its
remaining lack of self-consistency is included explicitly in the residual.
No bulk equation has been declared solved by enforcing the boundary law.

## Direct initial-data error

The initial report gives actual-minus-expression-defined-nominal caps
`delta_u,delta_ur,delta_v`, positive norm caps `nom_u,nom_ur,nom_v`, and P1
velocity interpolation error `eps_v`. Let the proved smoothing errors be

```
a_u = M2 delta_u_smooth^2/18,
b_u = 35 M2 delta_u_smooth/128,
a_v = (35 delta_v_smooth/256) max_knots ||slope_jump||_H,
delta_u_smooth=2^-30, delta_v_smooth=1/2048.
```

The same affine tangent/last-slope continuations and parity extension as the
published smooth reconstruction are used. In particular, the global C1,1
nominal displacement extension justifies the convolution remainder even near
the end of the initial footprint. Disjoint velocity smoothing bands make the
single-hinge error bound applicable globally.

The actual initial-data differences from this smoothed free wave are bounded
by

```
Delta_u  = delta_u  + a_u,
Delta_ur = delta_ur + b_u,
Delta_v  = delta_v  + eps_v + a_v.
e_initial = max((Delta_u+(3/2)Delta_v)/3, Delta_ur+Delta_v).
```

The corresponding free-wave size is

```
B0cap = max((nom_u+a_u+(3/2)(nom_v+eps_v+a_v))/3,
            nom_ur+b_u+nom_v+eps_v+a_v).
```

These formulas use only one radial derivative of displacement. `eps_v`
remains because the new velocity is the smoothed P1 velocity. The older
upstream feedback evaluator used independent P1 displacement and radial
derivative tables. The published PR94 smooth candidate and the present `W`
already share the derivative-consistent smooth nominal bulk source. Their
earlier transferred certificate nevertheless compared through that older P1
candidate and charged its large interpolation differences. The present
direct residual identity avoids this detour. It does not assume those
interpolation errors were zero, change the PR94 nominal bulk source, or show
that the boundary correction alone reduced the unknown true error by twelve.

## Whole-cone feedback and bulk residual

For any compatible field on the retained cone, the inherited zero-wedge lemma
gives `sup|V|_H <= 3 X(V)/4` and shell trace `sup|V|_H <= X(V)/4`.
Write `na0,nb0,nd0,nm0` for the frozen matrix norm caps and
`da,db,dd,dm` for actual-minus-frozen coefficient norm caps. The checker
verifies all block Gram-weighted Frobenius inequalities and every actual
coefficient interval difference before using these numbers. Set

```
na=na0+da, nb=nb0+db, nd=nd0+dd, nm=nm0+dm,
q0=(3/2)(na0+nb0+(3/4)nd0)+nm0/4,
qa=(3/2)(na+nb+(3/4)nd)+nm/4.
```

The frozen first-feedback candidate is `V=B0+B1`, where
`B1=K(L0 B0)+L(M0 PN B0)`. Consequently the whole-cone bound is

```
X(B1) <= q0 B0cap.
```

This is a whole-cone X bound, not merely the earlier shell trace bound.
The independently established Robin correction certificate gives
`X(W-V)<=Ccap`, where `Ccap<6.344e24`. Thus

```
X(W-B0) <= q0 B0cap+Ccap =: feedbackcap.
```

It follows that the lower-order bulk residual, before adding the affine `F`,
has uniform H-norm bound

```
Rcoeff=(da+db+(3/4)dd) B0cap,
Rfeedback=(na+nb+(3/4)nd) feedbackcap,
Rlower=Rcoeff+Rfeedback.
```

The affine-source convention matters. The inherited number
`SOURCE=52162547172519/2` bounds `X(K F+L G)`, not the pointwise norm of `F`.
The separately authenticated actual shell cap gives `sup||G||_H<487122`.
By the triangle inequality and `X(LG)<=sup||G||_H`,

```
X(K F) <= X(K F+L G)+X(L G) <= SOURCE+487122.
```

Therefore the inhomogeneous part of the exact error identity has X norm at
most

```
r = e_initial + (3/2)Rlower + SOURCE + 487122.
```

The affine response number is already in X units and is not multiplied by
`3/2` as though it were a pointwise bulk-source cap. Its full-H scope and the
zero-wedge integral-operator estimates are inherited from the authenticated
`HDBLAST_TOTAL_RESPONSE_PREFLIGHT_20261023/analytic/LOCAL_TOTAL_REVIEW.md`.

## Inverting the error operator

The arithmetic gives `qa < 3409/8000000 < 1`; the latter is the inherited
full23 contraction certificate. Both direct estimates are computed:

```
X(e) <= r/(1-qa) <= r/(1-3409/8000000) < 1.735e27.
```

The conservative second estimate alone suffices for the readable statement.
The prior bound transferred to `B0+B1` exceeds this new conservative upper
certificate by a factor greater than twelve. All smoothing, initial velocity
interpolation, coefficient uncertainties, remaining bulk feedback, boundary
correction and affine bulk source have been included. No separate two-term
Neumann tail should be added: the direct residual and contraction inversion
already account for all remaining coupled iterations.

The sign of `Rlower` is not determined; its number is an upper bound. The
actual reduced equations' first-derivative integral error is controlled, but
individual second derivatives, pointwise Einstein/constraint residuals,
finite-amplitude remainders and invariant physical observables are not.

## Reproduction and trust boundary

`checks/check_robin_direct_error.py --repo-root PATH` uses Python's standard
library and JSON arithmetic only. It imports no point-evaluation code. The
optional `--transfer-report` and `--source-report` arguments select separately
located authenticated report copies. The default transfer report is
`../reference/TRANSFER.json`; the default source report is the published PR94
`reference/AFFINE_SOURCE.json` under the supplied repository root.

Normal and optimized output are byte-identical. The checker authenticates the
earlier initial/profile/response reports, the Robin correction report and the
source/review inputs. It recomputes the complete matrix norm/error checks,
smoothing jump norm, initial-data differences, boundary-majorant arithmetic,
whole-cone feedback, residual and contraction division. The analytic meaning
and earlier full-field enclosures remain stated inherited premises; hashes
authenticate bytes rather than replace those proofs. The calculation is not
external specialist review, a new physical detection or a novelty claim.
