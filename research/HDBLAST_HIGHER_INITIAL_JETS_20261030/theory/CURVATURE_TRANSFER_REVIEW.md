# Curvature-matched endpoint extension: independent comparison proof

Let R be the right endpoint of the physical initial footprint, and let the
fixed nominal displacement have endpoint data u_R, u'_R and left curvature
`kappa=u''_nom(R-)`. The old exterior extension is affine. Replace it, only
for h>0, by

```
u_ext_dagger(R+h)=u_R+u'_R h+(1/2)kappa h^2.
```

The interior unsmoothed nominal function is unchanged. Exact curvature
matching requires this fixed nominal endpoint curvature, not a midpoint of
an unrelated curvature interval. A verified bound `kappa_H>=||kappa||_H`
can be used to enclose the resulting fixed candidate change. The actual
physical target and its initial prescription are not changed.

## Half-kernel charges on the physical initial interval

Use the same positive, even compact kernel as before,

```
k_delta(s)=35/(32delta) [1-(s/delta)^2]^3,    |s|<=delta.
```

Its mass is one, second moment is `delta^2/9`, and first absolute moment is
`35delta/128`. The exterior difference is `(1/2)kappa(z-R)_+^2`, so for
physical r<=R,

```
Delta u_s(r)=(1/2)kappa integral (r-s-R)_+^2 k_delta(s) ds,
Delta u'_s(r)=kappa integral (r-s-R)_+ k_delta(s) ds.
```

Both scalar weights are nonnegative. Their largest values on r<=R occur at
r=R, where only the negative half of the even kernel contributes. Therefore

```
||Delta u_s||_H <= kappa_H delta^2/36,
||Delta u'_s||_H <= kappa_H 35delta/256.
```

The displacement factor has two separate halves: the coefficient one half
in the quadratic extension, and the half-kernel second moment. The derivative
factor uses one half of the absolute first moment.

The exact checker also represents these weights by polynomials. With
q=(R-r)/delta in [0,1], write

```
F_n(q)=35/(32n!) integral_-1^-q (-q-t)^n(1-t^2)^3 dt.
```

The normalized value and derivative weights are F_2 and F_1. The checker
verifies `F_2'=-F_1`, `F_1'=-F_0`, the endpoint constants, and that their
Bernstein coefficients lie between zero and the claimed maxima. Convexity
of the Bernstein basis proves the bounds continuously over every q in
[0,1], independently of sampled values. For q>=1 all weights vanish.

The velocity prescription is unchanged. On the retained cone `x+r<=R`,
both free-wave initial arguments lie in `[-R,R]`. The reflected free-wave
formula has two half-weighted displacement terms; reflection preserves H. Hence,
in the registered norm
`X=max(sup||value||_H/3,sup||partial_x||_H,sup||partial_r||_H)`,

```
X(Delta Bs) <= C,
C=kappa_H max(delta^2/108,35delta/256).
```

This charges the change to the auxiliary smoothed initial data. It does
not silently identify the modified candidate with its predecessor.

## Noncommuting boundary and bulk comparison

Let `P=(I-L M PN)^(-1)`, and let b<1 bound the boundary operator in X.
Let k0 bound K L0 and k bound K La. The actual coefficient functions, M,
F and G are fixed and identical for both candidates. The frozen operator
L0 is also unchanged. Then the two definitions give exactly

```
Delta W0=P(I+K L0)Delta Bs,
Delta W1=P[Delta Bs+K La Delta W0].
```

The common affine G cancels in the first comparison; the same F and G
cancel in the second. No operator is commuted through P. The Neumann
bound `||P||<=1/(1-b)` and triangle inequalities imply

```
X(Delta W0) <= C(1+k0)/(1-b),
X(Delta W1) <= C[1+k(1+k0)/(1-b)]/(1-b).
```

If U is the same retained target and `X(U-W1)<=E_old`, the new comparison
certificate is `X(U-W1_dagger)<=E_old+X(Delta W1)`. This is a change of
candidate with an explicit cost, not an improvement in the unknown true
error or a replacement of physical initial data.

The manufactured operator test uses noncommuting rational bulk and boundary
matrices and a nilpotent boundary inverse. It checks the exact affine-source
cancellations and both norm bounds. Controls reject moving P through either
bulk operator, mismatching the fixed affine sources, and understating either
half-kernel charge.

## Causal support and preserved earlier point calculations

The smoothed displacement difference and its derivative vanish for
r<=R-delta. Their reflected negative-argument extensions have the same
support threshold in |r|, with the usual opposite derivative parity.
At any physical output (x,r) with `x+r<R-delta`, both free-wave initial
arguments obey `|r+x|,|r-x|<=r+x<R-delta`. Thus `Delta Bs=0` there.

Local coefficient operators cannot spread support out of this zero region.
For the reflected retarded bulk integral, every physical source point
(s,y) in its past domain satisfies `s+|y|<=x+r`; this follows from
`|y-r|<=x-s` and the triangle inequality. It remains valid when the complete
coupled output is reflected. A boundary lift at (x,r) uses only shell times
at most `x-r<=x+r`. Therefore K, the lift L, and every term of the boundary
Neumann inverse preserve the zero region. Uniform convergence of the
inverse preserves exact zero. Induction proves that every fixed finite
bulk update of the candidate difference vanishes when `x+r<R-delta`.

For `R=7/4` and `delta=2^-30`, the earlier event has `x+r=9/8<R-delta`.
Its field values and first jets are therefore unchanged exactly. This is
a support argument, not an inference that the numerical change is merely
too small to notice. The retained cone's shell horizon is only x<=5/4, so
the new outer-endpoint difference cannot reach that shell horizon either;
the generic boundary-inverse norm charges above remain safe conservative
bounds. A finite filtration test checks the same causal operator property
without replacing this analytic continuum argument.

These identities do not compute the actual endpoint curvature norm, supply
F_x, close the quantitative C2 gate, certify a full Einstein residual or
establish nonlinear or physical Big Bang realization. Those remain separate
scientific inputs and obligations.
