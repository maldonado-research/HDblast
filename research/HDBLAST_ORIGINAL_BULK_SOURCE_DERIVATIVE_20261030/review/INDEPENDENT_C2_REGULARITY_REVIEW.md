# Independent analytic review of the original-target C2 error gate

30 October 2026. **Assessment: PASS under the stated inherited physical-side
regularity, source and compatibility premises.** This is an internal analytic
review. Separate independent arithmetic reconstruction is required for the
printed rational caps. It does not establish external originality or a
physical interpretation.

## Candidate and target remain the original ones

Write `Bs` for the PR100 quadratic-continued smoothed-data free wave,
`L0=A0 partial_x+B0 partial_r+D0`, and
`La=A partial_x+B partial_r+D`. Let `P=(I-L M PN)^(-1)` be the boundary-only
inverse. The candidate is exactly the original-source finite update

```
W0=P[Bs+K L0 Bs+LG],
W1=P[Bs+K(La W0+F)+LG].
```

Neither its bulk source F nor its affine shell source G is averaged or
replaced. In particular, this is not the PR101 auxiliary candidate
`Z_epsilon`. The original actual second-coefficient target uses the same
F,G,M and actual coefficients.

Set `Dfield=W0-Bs`, `C=W1-W0`. The exact cancellation identities are

```
r0=La Dfield+(La-L0)Bs+F,
C=P K r0,
r1=La C.
```

The residual convention is actual RHS minus principal wave for the error
equation. No absolute W0/Bs second-jet subtraction is needed.

## Two differentiated prescribed constructions

`Dfield` has zero Cauchy displacement/velocity, but its time derivative has
initial velocity `(L0 Bs)(0)`. Subtracting the free-wave shell law first,
then differentiating gives boundary source `M PN Bs_x+G_x`. Thus

```
Dfield_x=P[Bfree[0,(L0 Bs)(0)]
           +K(L0 Bs)_x+L(M PN Bs_x+G_x)].
```

The actual M is constant in time. Dropping `M Bs_x` would omit part of the
affine boundary residual of `W0-Bs`, despite Bs satisfying its own homogeneous
free-wave reflection law. The free-wave X bound for zero displacement and
initial velocity f is `sup||f||_H` on the length-3/2 time interval.

The next correction has

```
C_x=P[Bfree[0,r0(0)]+K r0_x].
```

Here `r0(0)` is generally nonzero in the active initial packet and contains
the actual F(0). It must be charged. The gate differentiates the complete
cancellation-preserving expression for r0 before applying norm bounds:

```
r0_x=La Dfield_x+A_x Dfield_x+B_x Dfield_r+D_x Dfield
      +(La-L0)Bs_x+A_x Bs_x+B_x Bs_r+D_x Bs+F_x.
```

The separate `D Dfield_x` and `D Bs_x` lower-order terms are already part
of the two displayed operators. In the norm ledger their bounds use the
known first jets of Dfield/Bs, rather than paying for their second derivatives.
The complete actual coefficient commutator and the actual original F_x are
retained. This gives an explicit bound for `X(C_x)` and hence for

```
r1_x=La C_x+A_x C_x+B_x C_r+D_x C.
```

These prescribed-forcing constructions establish physical-side C2 regularity
of W0 and W1 before using a differentiated target-error equation. The initial
smoothed data are C2, G has its retained continuous time derivative, and
the now-proved complete actual F_x exists continuously. The physical source
is continuous in the radial coordinate. The original open zero collar makes
every required initial/shell corner compatible. Odd reflected output forcing
can jump on an auxiliary negative-depth extension; no C2 extension across
the doubled shell is used.

## Original target-error equation and initial acceleration

Let `e=U_target-W1`. Both fields satisfy the same actual affine shell law,
so the error has homogeneous physical D/Robin rows. Its principal equation is

```
e_xx-e_rr=La e+r1,
PN e_r=M PN e,       PD e=0.
```

Although `C_x` has the nonzero initial velocity r0(0), C itself has zero
initial value and first jets. Since La is a first-order operator,
`r1(0)=La C(0)=0` exactly. This differs from the filtered candidate's initial
target-error forcing and does not set the actual F(0) to zero.

For `v=e_x`, differentiating the original error equation gives

```
v_xx-v_rr=La v+A_x e_x+B_x e_r+D_x e+r1_x,
PN v_r=M PN v,       PD v=0,
v(0)=e_velocity,
v_x(0,r)=e0_rr+A(0,r)e_velocity+B(0,r)e0_r+D(0,r)e0.
```

The last formula retains the full homogeneous lower-order initial
acceleration. No unfiltered/filtered source cancellation is interchanged.
All initial error profiles retain the flat initial shell collar.

The retained C1 error contraction can also justify existence of this
derivative, not merely bound a derivative that is assumed to exist.
Differentiate its Picard iterates in time. Every iterate has the same error
Cauchy value/velocity, so from the first updated iterate onward the
differentiated fields have the common Cauchy data displayed above; the free
zeroth iterate is the finite exception. Their differences satisfy

```
delta v_n <= q delta v_(n-1)+T comm delta e_(n-1),
```

where q is the same full coupled bulk-plus-boundary contraction and
`comm=||A_x||+||B_x||+(3/4)||D_x||`. The lower differences decay geometrically;
the new differences are bounded by summable tails of the form `n q^n`.
Thus the differentiated iterates converge in C1, and the limit is e_x.
The principal equation then recovers e_rr continuously. This supplies
classical physical-side C2 error jets with one-sided shell/domain traces.

## Norm constants and scope

All fields vanish in the same zero wedge, so `||e||_H<=(3/4)X(e)` and the
same value factor applies to the differentiated fixed-point operator.
With T=3/2, retained C1 error cap E, initial bounds `ev,evr,eurr,eur,eu`,
and `F1>=sup||r1_x||_H`, define

```
w0=eurr+||A||ev+||B||eur+||D||eu,
J=max((ev+T w0)/3,evr+w0),
Z=[J+T(F1+comm E)]/(1-q).
```

Then `||e_xx||_H,||e_xr||_H<=Z`. With
`c=||A||+||B||+(3/4)||D||` and the unchanged original defect cap F0,

```
||e_rr||_H<=Z+c E+F0.
```

The gate's boundary-only inverses for Dfield and C use `1/(1-b)`; the final
error uses the full feedback inverse `1/(1-q)`. Confusing those inverses
would change the claimed theorem. The individual initial error caps for the
quadratic candidate were re-established in PR100; they are not the absolute
candidate jets from the filtered construction.

The numerical result, once its arithmetic and F_x input are independently
checked, is an absolute scaled-coordinate error bound for the original
candidate. A large majorant dominated by the initial-curvature allowance is
not a lower bound on the actual error or proof that its true error has that
distribution. Physical relative accuracy, covariant equation residuals,
constraint propagation, nonlinear realization and cosmological origin
require their separate conversions and evidence.
