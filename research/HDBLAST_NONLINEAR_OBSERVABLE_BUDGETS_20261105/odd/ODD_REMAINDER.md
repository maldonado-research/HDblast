# Sharp odd C2 remainder and an area-compatible error interface

Independent mathematical working note, 5 November 2026 (UTC).

The two supplied notes define a selected equivariant actual-area moment and
prove, conditionally on the common observer/clock domain, that
`M:(-a0,a0)->C(I)` is `C2`, `M(-e)=-M(e)`, `M(0)=0`, and `M'(0)=S`.
This note uses precisely those hypotheses. All parameter derivatives hold the
actual proper time and transported angular labels fixed. No coordinate-time
substitute or new nonlinear existence assertion is used.

The completion, affine parameter, quarter-turn equivariance and common admitted
observer domain are the same selected ones as in the supplied bridge. In its
notation the retained coarse interval is
`I=[-rho_b*6145ell/4096,-rho_b*6143ell/4096]`. Its endpoint signs used below
are not reassigned to a newer or smaller linear-cell bracket. The full physical
scalar first tangent on these fixed actual-clock cuts and labels is
`A'_0(tau,Omega)=S(tau/rho_b)Y(Omega)`; for brevity `S` below denotes this
proper-time function.

## 1. Signed identity and sharp anchored modulus

Write `h=|e|>0`, and use the supremum norm on `C(I)`. Differentiating oddness
twice gives `M''(-e)=-M''(e)`, in particular `M''(0)=0`. Integral Taylor, with
the oriented integral also valid at negative `e`, gives

```
M(e)/e-S = e integral_0^1 (1-t) M''(te) dt
          = (1/h) integral_0^h (h-s) M''(s) ds.             (1)
```

The second expression follows because the two minus signs cancel when
`e=-h`. Thus the normalized error is even; averaging the positive and
negative measurements does not cancel it.

Suppose a nonnegative nondecreasing integrable envelope `omega(s)` has been validated with
`||M''(u)-M''(0)|| <= omega(s)` whenever `|u|<=s<=a`. A usual choice is the
nondecreasing anchored modulus
`omega(s)=sup_{|u|<=s} ||M''(u)-M''(0)||`. Then

```
||M(e)/e-S|| <= E_omega(h),
E_omega(h) = (1/h) integral_0^h (h-s) omega(s) ds
           <= h omega(h)/2.                              (2)
```

Continuity makes the actual anchored modulus tend to zero. Using that modulus,
or any certified envelope tending to zero, proves `M(e)=eS+o(e^2)` uniformly
on `I`. An arbitrary loose majorant need not tend to zero. No particular power
is supplied by `C2`. The expression
`E_omega(h)=h integral_0^1(1-t)omega(ht)dt` is nondecreasing.

The integral bound is sharp as an assertion from an anchored modulus alone.
For any continuous nondecreasing `omega` with `omega(0)=0`, define a scalar
odd function by

```
M(e)=eS+sgn(e) integral_0^|e| (|e|-s) omega(s) ds.
```

Its second derivative is `sgn(e)omega(|e|)`, extended by zero at zero. It is
`C2`, has the stated anchoring and symmetry, and attains (2)'s integral
bound. No derivative of `omega` is needed.

For the anchored Hölder envelope `omega(s)=L s^alpha`, `0<alpha<=1`,

```
E_omega(h)= L h^(alpha+1)/((alpha+1)(alpha+2)).             (3)
```

Indeed the integral is
`h^(alpha+2)[1/(alpha+1)-1/(alpha+2)]`. Equality is attained by

```
M(e)=eS + L sgn(e)|e|^(alpha+2)/((alpha+1)(alpha+2)).
```

This is `C2` for every `alpha>0`; for `alpha<1` it is not `O(e^3)` after
subtracting `eS`. For `alpha=1`, a Lipschitz bound on `M''` gives the cubic
unnormalized bound `L|e|^3/6` without requiring a third derivative.

At the endpoint `alpha=0`, the integral of a constant budget `L` gives
`Lh/2`. A positive constant is not a continuity modulus tending to zero.
Also `sgn(e)|e|^2` is not `C2`, so it cannot establish equality in the
selected class. Nevertheless `1/2` is the optimal universal coefficient
for the coarse budget `||M''||<=L`: take the continuous odd second
derivative `L clip(e/delta,-1,1)`, and integrate twice with zero constant
and linear parts. For `h>=delta`, its normalized error is exactly

```
L[h/2-delta/2+delta^2/(6h)] -> Lh/2 as delta->0.
```

This endpoint distinction prevents an incorrect sharpness claim.

## 2. A full two-point modulus contains more information

If instead the certified input is a nonnegative nondecreasing pairwise
modulus `sigma`, with

```
||M''(u)-M''(v)|| <= sigma(|u-v|),  u,v in [-a,a],
```

oddness gives two valid anchored envelopes:

```
||M''(s)|| <= sigma(s),
||M''(s)|| = ||M''(s)-M''(-s)||/2 <= sigma(2s)/2.
```

Consequently (2) can use
`eta(s)=min{sigma(s),sigma(2s)/2}`. The modulus must be certified up to
separation `2a`; no value beyond its certified domain is implied.

For a full pairwise Hölder estimate `sigma(r)=L r^alpha`, `0<alpha<=1`,
this gives the sharp bound

```
||M(e)/e-S|| <= 2^(alpha-1)L h^(alpha+1)
                         /((alpha+1)(alpha+2)).           (4)
```

To verify sharpness, the scalar function `g(e)=sgn(e)|e|^alpha` has
pairwise Hölder seminorm `2^(1-alpha)`. For equal signs use
`|x^alpha-y^alpha|<=|x-y|^alpha`; for opposite signs concavity gives
`x^alpha+y^alpha<=2^(1-alpha)(x+y)^alpha`, with equality at `x=y`.
Thus choosing `M''=2^(alpha-1)L g` attains (4). This improvement must not
be silently applied to an anchored estimate with constant `L`.

## 3. What the quarter-turn actually cancels

Let `T f=f o R`, where
`R(n1,n2,n3,n4)=(-n2,n1,n3,n4)`. Then `T^4=1`, `TY=-Y`,
`T(Y^2)=Y^2`, and round measure is invariant. Use the full observer scalar
`A_e=Phi_e^obs-phi_b` and the full area density `j_e`, so

```
A_-e=T A_e,  j_-e=T j_e,
N_e=c integral Y A_e j_e dOmega,
D_e=c integral Y^2 j_e dOmega.
```

A change of variables gives `N_-e=-N_e`, `D_-e=D_e`; hence `M=N/D`
is odd. It does not assert `A_-e(Ω)=-A_e(Ω)` or
`j_-e(Ω)=j_e(Ω)` at a fixed angular label.

For arbitrary angular fields the genuine cyclic character projectors are

```
P_-=(1-T+T^2-T^3)/4,   P_+=(1+T+T^2+T^3)/4.
```

Since `T` has order four, `(1-T)/2` is not a projection on arbitrary
fields. Equivariance applied twice does imply `T^2 A_e=A_e` and
`T^2 j_e=j_e`, and products and their parameter derivatives inherit
this property. On these fields the projectors reduce to half sums and
differences. With `F_e=A_e j_e`, exactly

```
N_e=c integral Y P_-F_e dOmega,
D_e=c integral Y^2 P_+j_e dOmega,
P_-F''_s=(F''_s-F''_-s)/2,    P_-F''_0=0.                 (5)
```

Only the indicated angular sector vanishes pointwise at zero. The full
`F''_0`, `A''_0`, and `j'_0` can all be nonzero. In particular the
usable integrated estimate is

```
|N''_s| <= (c/2) integral |Y| |F''_s-F''_-s| dOmega.       (6)
```

It improves an absolute budget if that difference is certified smaller.
Parity by itself supplies neither its numerical size nor a power law.
Likewise a central difference of an odd `M` is
`[M(h)-M(-h)]/(2h)=M(h)/h`; it retains the entire error in (1).
An odd integrand integrated symmetrically vanishes, but the Taylor kernel
and the sign of the parameter transform together in the normalized
remainder. Such a vanished symmetric integral cannot replace (1).

## 4. Exact normalized numerator/denominator interface

Normalize out the common positive background factor using
`Q=integral Y^2 dOmega=pi^2/3`,
`n=N/(cQ)`, `d=D/(cQ)`. Then `n` is odd, `d` is even,
`n'(0)=S`, `d(0)=1`, `d'(0)=0`. Let a certified `m(a)>0` satisfy
`d(e,tau)>=m(a)` on the common domain. Pointwise `j_e>=m(a)` is one
sufficient way to obtain it. Exactly,

```
M(e)/e-S = [n(h)/h-S-S(d(h)-1)]/d(h)
 = {h^-1 integral_0^h(h-s)n''(s)ds
       -S integral_0^h(h-s)d''(s)ds}/d(h).                (7)
```

Thus envelopes `||n''(s)||<=eta_n(s)` and
`||d''(s)||<=beta_d(s)` give the robust bound

```
||M(e)/e-S|| <= [h^-1 integral_0^h(h-s)eta_n(s)ds
     +||S|| integral_0^h(h-s)beta_d(s)ds]/m(a).           (8)
```

When taking separate absolute values this may overcount real
numerator/denominator cancellation. The denominator term is generally
`O(h^2)` if `d''` is bounded; parity does not remove it. One may also
use `d(h)>=1-integral_0^h(h-s)beta_d(s)ds` when positive.

## 5. A residual interface preserves the exact cancellation

Define

```
B_e=A_e-e S Y,       H_e=B_e j_e,
q_e=(1/Q) integral Y H_e dOmega = n_e-e S d_e.
```

The retained first field gives `B_0=B'_0=0` pointwise. Equivariance also
holds for `B,H`; hence `q` is odd, `q_0=q'_0=q''_0=0`. Exactly,

```
M(e)/e-S = q_e/(e d_e),
H''_e=A''_e j_e+2(A'_e-SY)j'_e+(A_e-eSY)j''_e.           (9)
```

If a certified envelope satisfies

```
eta_H(s) >= sup_tau (1/Q) integral |Y| |P_-H''_s| dOmega,
```

then the sharp triangular-kernel interface is

```
||M(e)/e-S|| <= [1/(h m(a))]
                   integral_0^h(h-s)eta_H(s)ds.           (10)
```

Equivalently use the half-difference
`P_-H''_s=(H''_s-H''_-s)/2`; a certified pairwise field modulus can
therefore improve this input. In contrast to independent estimates in
(8), (10) is exactly zero whenever `A_e=eSY`, for any admissible area
density. Thus area variation alone does not produce an error for that
exact scalar profile.

A concrete conditional raw-field input for (10) is especially economical.
On `|e|<=a` suppose the common-domain uniform bounds

```
E2 >= ||A''_e||_infinity,
J1 >= ||j'_e||_infinity,    J2 >= ||j''_e||_infinity,
w_A2(s) >= sup_|u|<=s ||A''_u-A''_0||_infinity
```

have been certified. From `j_0=1`, `|j_e|<=1+sJ1`, and integrating
`B''=A''` gives `||B'_e||<=sE2`, `||B_e||<=s^2 E2/2`.
Subtracting `H''_0=A''_0` in (9) yields

```
||H''_e-H''_0||_infinity
 <= (1+sJ1)w_A2(s)+3sE2J1+(s^2/2)E2J2,   |e|<=s.      (11)
```

The coefficient three bounds one contribution from
`A''_0(j_e-1)` and two from `B'_e j'_e`. A slightly finer coefficient is
`(||A''_0||+2E2)J1`. As `P_-H''_0=0`, the right side of (11), multiplied
by `r=(integral|Y|)/Q=6/pi`, is a valid `eta_H(s)`. Integrating gives

```
||M(e)/e-S|| <= (r/m(a))[
 h^-1 integral_0^h(h-s)(1+sJ1)w_A2(s)ds
 + E2J1 h^2/2 + E2J2 h^3/24 ].                          (12)
```

This requires a modulus of the full pulled-back scalar's second
parameter derivative, not a third derivative; a bound on `j''` suffices
here. If `w_A2(s)<=L_A s^alpha`, the integral contribution in brackets
is explicitly

```
L_A h^(alpha+1)/((alpha+1)(alpha+2))
 + L_A J1 h^(alpha+2)/((alpha+2)(alpha+3)).               (13)
```

The rest of (12) remains present. Bounds of the form (11)–(13) are
stronger small-parameter information than a fixed absolute `K`, provided
their missing inputs are actually computed.

For example `m(a)=1-aJ1` is available when that value is positive, from
`j_0=1` and the same full-domain derivative bound. It is still a conditional
margin until the bound and domain are validated.

For a check that the area coupling is real, take the admissible algebraic
test fields `A_e=eSY+C e^2/2` and `j_e=1+b eY`, with `|be|<1`.
They obey the required equivariance; `A''_0=C` and `j'_0=bY` are nonzero.
Using `integral Y=integral Y^3=0` gives

```
d_e=1,       M(e)/e-S = Cb e^2/2.
```

It is the integrated coupling of an even scalar second variation to an
odd area first variation. Discarding those pointwise fields because
their separate projected expansion-point terms vanish would miss it.
This is a test of the observable algebra, not a claimed PDE solution.

## 6. What follows, and what is still missing

The retained endpoint facts are `S(t_L)>63` and `S(t_R)<-53`. Any
validated uniform bound from (2), (4), (8), (10), or (12) that is below
`53` preserves opposite signs of `M(e)/e`, for either sign of `e`.
Continuity then gives a moment zero in the open actual proper-time
interval. For a monotone anchored modulus, a sufficient all-parameter
guard is

```
(1/a) integral_0^a(a-s)omega(s)ds < 53,
```

which improves the coarse sufficient guard `a omega(a)<106` whenever
the interior modulus is smaller than its endpoint value. Negative
parameters reverse the signs of `M` itself and preserve the zero.

None of the two source notes evaluates the modulus, the raw-field bounds
used here, the positive denominator margin, or the common admitted
observer/parameter radius. This calculation therefore supplies a precise
error interface and conditional guards, without a numerical constant or
effective amplitude window. It proves no full-field zero, uniqueness,
new nonlinear evolution result, or physical origin claim. No PDE solver or
field oracle was used.
