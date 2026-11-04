# An independent elementary bound for the preparation slope

This component uses the same single preparation direction as the main
whole-cell calculation, with a different analytic bounding method. It uses no
quadrature, scientific producer import, or new oracle call. The background and
five-mean contraction premises are the same retained premises, so it is not
an independent PDE observation. The full preparation theorem and physical
interpretation belong to the parent derivation.

Set `ell=10^-6`, `T=5/4`, `Delta=1/4`, and fix

```
Q(r)=exp(4-1/[(r-1)(2-r)]) for 1<r<2, zero otherwise,
chi=ell rho_b^3 Q,
h=rho^-2 chi_y=-(rho_b/rho)^3 Q_r.
```

The zero moment `integral rho^2h dy=0` is exact. In the conditioned mean
order `(K,L,N,J,w)`, the homogeneous difference initial data are

```
W0=(h,0,0,h,s), W0_x=(0,V_L,0,0,0), W0_r=(h_r,0,0,h_r,s_r),
s=-p zeta,
zeta=3(h_y-ah)/(2p^2)+9chi/(2rho^4p^2),
V_L=ell rho(-h_y+ah).
```

Write `R=(rho_b/rho)^3`, `b=ell rho a`,
`c=-3rho_b^3/(2ell rho^4p)>0`, and `d=ell rho p_y/p`. Direct
differentiation, using `rho_r=-ell a rho^2`, gives

```
h=-R Q', h_r=-R(Q''+3bQ'), V_L=-R(Q''+4bQ'),
s=c(Q''+4bQ'+3ell^2Q),
s_r=c[Q'''+4bQ''+3ell^2Q'
      +(4b+d)(Q''+4bQ'+3ell^2Q)
      -4ell^2rho^2(a^2+a_y)Q'].
```

The scalar's initial coordinate velocity is zero. Its free even-reflected
boundary value is exactly `s(x)`. The complete mean couplings remain in the
homogeneous problem; none is dropped from the error bound.

For the backward-closed triangle `x+r<=T`, only the initial segment
`0<=r<=T` is needed. The data vanish on `r<=1`. Let `t=r-1` and
`D=t(1-t)`. On `1<r<=T`, `0<D<=3/16`, and

```
Q'=Q(1-2t)/D^2,
Q''=Q P2(D)/D^4, P2(D)=1-6D+6D^2,
Q'''=Q(1-2t)P3(D)/D^6, P3(D)=1-10D+24D^2-12D^3.
```

These identities follow from `D'=1-2t`, `(D')^2=1-4D`, `D''=-2`.
In particular

```
P2=(1-4D)(1-2D)-2D^2,
P3=P2+D^2P2'-4DP2.
```

The producer checks those polynomial identities with exact rational
polynomial arithmetic. It checks `P2(3/16)=11/128` and
`P3(3/16)=-113/1024`. Their derivative upper bounds are respectively
`-6+12(3/16)<0` and `-10+48(3/16)<0`, where dropping the negative quadratic
term in `P3'` is an upper bound. Thus `0<P2<=1` and `|P3|<=1` on the full
interval. Therefore `Q,Q',Q''` and `s` are positive on the open history.

Put `w=1/D>=16/3`. The maximum of `e^-w w^j` on this interval occurs at
`w=16/3` for `j=2,4`, and at `w=6` for `j=6`; this follows by
differentiating, whose sign is that of `j-w`. Rational Taylor estimates give

```
15/4 < exp(4/3) <4,
exp(2)>6^6/6400.
```

For the first lower bound use the exponential sum through degree five.
For the upper bound use that sum plus its degree-six term divided by
`1-(4/3)/7`; successive tail term ratios are at most `(4/3)/7`.
For the last bound use the exponential sum through degree six. The producer
checks all these rational inequalities exactly. They imply

```
Q<=1/3, |Q'|<=8, |Q''|<=216, |Q'''|<=6400
```

continuously on the entire causal initial segment, including the flat
endpoint. Indeed `exp(-4/3)<4/15`, and
`(4/15)(16/3)^2<8`, `(4/15)(16/3)^4<216`; the third derivative is at most
`exp(-2)6^6`. There is no unresolved tail or sampling argument.

The two pinned published JSON sources imply, in exact interval arithmetic,

```
129<=rho,rho_b<=131,
3/25000<=-p<=1/7000,
111/1000<=a<=14/125,
|p_y|<=1/4000, |a_y|<=1/10000,
1/1200<=-F'_b<=1/1000.
```

The annulus static source explicitly covers every point of the packet
annulus. Its first jet entries equal first derivatives; no factorial changes
them. The producer verifies byte hashes, the continuous domain, and each
inclusion before using the broad constants. The shell `rho_b,F'_b` inputs
are exact enclosures of the selected background, not replacement background
values.

Insert these bounds in the profile formulas to obtain rational upper bounds
`H,Hr,V,S,Sr` for `|h|,|h_r|,|V_L|,|s|,|s_r|`. With Gram matrix
`diag(1,2,1,3,1)`, safe initial norm caps are

```
D0=2H+S, Dr=2Hr+Sr, Dv=2V,
N=max((D0+T Dv)/3,Dr+Dv),
E=qN/(1-q), q=1201/3200000.
```

These use only the triangle inequality and `sqrt(2)<=2`. The retained
homogeneous resolvent and exact zero wedge give

```
|m(x)-s(x)|<=(x-1)_+ E,
|integral k(m-s)|<=Delta^3E/6,
|integral k_T(m-s)|<=Delta^2E/2.
```

All mean couplings and homogeneous Robin feedback are included in `q`.
The actual first forcing, proper clock/measure correction, and every `Q1`
term cancel only when forming this difference. They remain in the published
absolute base intervals.

The exact kernels satisfy `0<=k<=T-x`, `0<k_T<=1`, and on this target
interval both comparison ratios are at least `1-eta`, `eta=10^-11`:

```
k/(T-x)>=1-eta, k_T>=1-eta.
```

For the derivative, with `v=ell T,b=ell x`, the useful exact identity is
`k_T=cosh(b)^2cosh(v-2b)/cosh(v)^3`. This is at most one because
`0<=b<=v` and `|v-2b|<=v`. Its lower bound is
`(1+v^2)^-3>=1-3v^2>1-eta`. The analogous lower bound for `k/(T-x)` is
`(1+v^2)^-2>=1-2v^2`; the ratio is extended continuously at the endpoint.

Let `c_min=3*129^3/(2ell*131^4*(1/7000))`. Since `s>=c_min Q''>0`,
integration by parts and flatness give

```
integral_1^T (T-x)s dx>=c_min Q(T)>c_min/4,
integral_1^T s dx>=c_min Q'(T)>32c_min/9,
s(T)>1408c_min/81.
```

Here the exact endpoint values are
`Q(T)=exp(-4/3)`, `Q'(T)=(128/9)exp(-4/3)`, and
`Q''(T)=(5632/81)exp(-4/3)`. The strict coupled sign tests are therefore

```
1408c_min/81-Delta E>0,
(1-eta)c_min/4-Delta^3E/6>0,
(1-eta)32c_min/9-Delta^2E/2>0.
```

All pass by exact rational arithmetic. Multiplication by the retained
negative `F'_b/2`, with the positive proper-time factors, proves
`L_M>0`, `L_v<0`, and `L_Theta<0` for the complete homogeneous response.

The producer also applies the registered deterministic rule: choose
`lambda=-10^n`, where `10^n` is the least nonnegative-integer power of ten
strictly above both ratios of the absolute negative base lower bound to
the absolute negative slope upper bound. It selects `lambda=-10^21`.
Propagating the old absolute base intervals and the new slope intervals
proves positive `v2_lambda` and `Theta2_lambda`. The exact JSON records the
full intervals. This algebraic postprocessing does not change the direction,
target, proof method or sign tests.

The exact-family conclusion is conditional and local in the original packet
parameter for this fixed finite `lambda`. There is no uniform or effective
parameter radius and no contracting total-expansion claim.
