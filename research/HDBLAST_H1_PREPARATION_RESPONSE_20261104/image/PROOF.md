# H1 control does not place the actual preparation image in the declared C1 class

4 November 2026. Analytic image component of the accompanying H1 preparation-response checkpoint. This is an analytic construction in the retained admissible compact
mean-preparation image, not a construction of arbitrary scalar functions that
might fail its constraints.

The conclusion has a limited target: every positive H1 budget still permits
unbounded `R_A=max(sup|D|,sup|D_r|)`, and permits unbounded initial supremum
quantity `N` used by the existing five-mean resolvent argument. Thus an H1
budget alone cannot imply the already declared sufficient `R_A<=2e33` class
condition by that route. This does **not** rule out direct H1-to-observable
bounds, a different negative-sign certificate, or energy/trace estimates that
use Green-function smoothing. No new sign-reversing family is asserted.

## Retained objects and actual image

Keep the selected background, first fields, parameter normalization, physical
initial Cauchy section and background scalar-label parametrization of the
4 November intrinsic-class theorem. Write

```
I=[1,T], T=5/4, ell=10^-6, dy/dr=-ell rho,
b=ell rho a, g=ell rho p^2/(3a)>0.
```

Here `rho>0`, `a>0`, `p!=0` on the retained compact annulus. `D` is the
relative difference in the actual second area variation of the connected
initial scalar-level component, with its label held fixed. Its derivative is
the fixed background-calibrated derivative `D_r=(-ell rho p)partial_v D`.

For an admissible compact mean direction, the exact image relation is

```
z_r=3bz-h,
h_r=3ell^2 z+(-b+g)h-(2/3)gD,
z=chi/(ell rho^3), chi_y=rho^2 h.                              (1)
```

Conversely, choose any smooth `z` with support strictly inside `(1,T)` and
define

```
chi=ell rho^3 z, h=3bz-z_r,
D=(3/(2g))[3ell^2 z+(-b+g)h-h_r].                             (2)
```

The identity `rho_r=-ell a rho^2` gives
`chi_r=-ell rho^3 h`, hence `chi_y=rho^2 h` exactly. Both `chi` and `h`
are compact and flat before either packet endpoint; therefore

```
integral_(yL)^(yU) rho^2 h dy=chi(yU)-chi(yL)=0.               (3)
```

This is exact, with no moment quadrature or numerical cancellation. The
retained local construction of `zeta` then supplies the homogeneous spatial
metric datum

```
H_h=h gamma0+Lie_(zeta partial_y)gamma0,
zeta=3(h_y-a h)/(2p^2)+9chi/(2rho^4 p^2).
```

Its constraint is the same exact homogeneous constraint as before. Added
canonical momentum and both scalar canonical sectors are zero; they are not
varied by this sequence. The scalar appearing in the conditioned adapter is
the same pure second scalar pullback `s=-p zeta`, not a new raw canonical
scalar perturbation. The initial collar, first packet, first sources and
second TT pair stay fixed. Thus (2) is in the actual retained image with all
those canonical data fixed.

As in the retained theorem, each fixed direction is completed by its own
conditional exact family. The sequence gives no common nonlinear family
neighborhood as its high derivatives grow.

## Eliminate h and isolate the leading derivative

Put

```
A=3/(2g), c=-2b-g, d=3(ell^2-b^2+bg-b_r).
```

Substitution into (2) yields exactly

```
D=A(z_rr+c z_r+d z),                                          (4)
D_r=A z_rrr+(A_r+Ac)z_rr+(A_r c+A c_r+Ad)z_r
                  +(A_r d+A d_r)z.                           (5)
```

For coefficient bounds in these formulas it suffices that `b` is C2 and
`g` is C1 and bounded away from zero on a compact interior interval. The
retained C4 background supplies these derivatives: `a,p` use one spatial
background derivative, and `d_r` uses at most two derivatives of `b`.
The estimates do not require higher background derivative bounds or any new
background query.

Literal smooth admissibility also uses the smooth selected background of
the retained preparation construction. Its fixed static equations have
smooth coefficients on this nonsingular annulus; the usual ODE bootstrap
gives smoothness there. This is distinct from the retained lower C4 topology
for parameter dependence. Bare C4 data alone should not be described as C∞;
if that inherited smooth-background setting were discarded, the formulas
would first give a finite-regularity sequence and smooth admissibility would
need a separate argument. No such weakening is made here.

## One explicit smooth compact bump

Take

```
psi(u)=u^3 exp(1-1/(1-u^2)) for |u|<1, and zero otherwise.
```

This is C∞ and flat at `u=+-1`. To see this without numerical evaluation,
let `w=(1-u^2)^-1>=1` and `E=exp(1-w)`. Every derivative is `E` times a
polynomial in `u,w`; `E w^k` tends to zero as `w` tends to infinity, for
every fixed integer `k`. This proves the flat extension inductively. At the
center the exact jets are

```
psi(0)=psi'(0)=psi''(0)=0, psi'''(0)=6.                        (6)
```

There are explicit conservative continuous bounds, not sampled bounds:

```
sup |psi^(j)| <= L_j, (L_0,L_1,L_2,L_3)=(1,11,1302,550638).
```

Indeed `d(EP)/du=E(P_u+2u w^2 P_w-2u w^2 P)`, and the maximum of
`exp(1-w)w^k` for `w>=1` is at most `k^k` for `k>=1`; for `k=0` it is
at most 1. Expanding the first three derivatives and summing absolute
monomial bounds gives the displayed integers. Thus, with
`B_j=2L_j`, `||psi^(j)||L2(-1,1)<=B_j`. The factor 2 safely bounds the
square root of the support length. The small rational verifier checks these
polynomials, constants and center jets; it evaluates no exponential.

## Concentration inside the admissible image

For the proposed sequence let `r0=9/8` and

```
z_n(r)=n^(-5/2) psi(n(r-r0)).
```

All integers `n>=9` have support strictly inside `(1,T)`. It is convenient
to take the exact subsequence `n=m^2`, with integers `m>=4`:

```
z_m(r)=m^-5 psi(m^2(r-r0)).                                   (7)
```

Its support lies in the fixed interior interval `[17/16,19/16]`. Define
`h_m,chi_m,D_m` by (2). Their exact support, flatness and zero moment follow
from (3). They are actual admissible directions, not independently chosen
profiles.

For `j=0,1,2,3`, change of variables gives the exact scaling

```
||z_n^(j)||2=n^(j-3)||psi^(j)||2,
||z_m^(j)||2=m^(2j-6)||psi^(j)||2.                            (8)
```

In particular the L2 norm of the third derivative is bounded uniformly even
though its pointwise size grows. Let `F_0=Ad,F_1=Ac,F_2=A`, and on the fixed
interior interval define finite constants

```
M_j=sup |F_j|, J_j=sup |(F_j)_r|, j=0,1,2,
C=sum_(j=0)^2 [(M_j+J_j) B_j+M_j B_(j+1)].                    (9)
```

`C` is finite and positive since `A>0`. By (4), (5), (8), and `m>=1`,

```
||D_m||H1 <= ||D_m||2+||(D_m)_r||2 <= C,                     (10)
||D||H1^2 := integral_I (|D|^2+|D_r|^2)dr.
```

For example `||D_m||2=O(m^-2)`, whereas its derivative's L2 norm stays
bounded. No numerical value of a background derivative is needed for (9):
its finite value defines the scaling below.

At `r0`, (6) annihilates every lower-order term in (4), (5), so there is
no remainder or lower-order asymptotic estimate:

```
z_m=z_m,r=z_m,rr=0, h_m=h_m,r=D_m=0 at r0,
(D_m)_r(r0)=6A(r0)m=9m/g(r0).                                (11)
```

Therefore `R_A(h_m)>=9m/g(r0)` tends to infinity. Also
`||D_m||infinity=O(m^-1)` tends to zero; this is consistent with the
one-dimensional H1 embedding into C0. The failure concerns derivative
suprema, not the profile amplitude.

## Every positive H1 budget and the existing initial N

Fix **any** `B>0`. Multiplying the entire homogeneous datum by
`alpha=B/(2C)>0` preserves support, exact zero moment, fixed canonical
sectors and the same retained completion argument. Its actual geometric
profile is `alpha D_m`, with

```
||alpha D_m||H1 <= B/2 < B,
R_A(alpha h_m) >= 9 alpha m/g(r0) -> infinity.                (12)
```

For any finite cutoff `K`, an integer
`m>max(3,K g(r0)/(9alpha))` violates `R_A<=K`. Thus no positive H1 budget,
however small, implies the declared `R_A<=2e33` condition on this actual
image. In particular there is no finite function of the H1 budget alone
that bounds `R_A` over it.

The same sequence obstructs feeding an H1 bound into the **existing**
initial-supremum `N` estimate. In the conditioned mean data,

```
s=(p/(2a))(h-2D/3),
s_r=(p/(2a))_r(h-2D/3)+(p/(2a))(h_r-2D_r/3).
```

At the center (11) gives exactly

```
(s_m)_r(r0)=-3p(r0)m/(a(r0)g(r0)).                            (13)
```

The Gram norm of `W0_r=(h_r,0,0,h_r,s_r)` there is `|s_r|`, so the genuine
initial derivative supremum, and hence every sufficient `N` of the retained
form `max((D0+T Dv)/3,Dr+Dv)`, diverge even after the fixed scaling in (12).
The nonzero `p` bound makes the factor in (13) strictly positive in absolute
value. This is a failure of a particular supremum-norm bridge, not a proof
that the observable differences diverge.

## H2 supplies a different, extra bridge

For any `f` in H1 on an interval of length `L`, select a point whose value
is no larger than its RMS value and integrate its derivative. This gives

```
||f||infinity <= L^(-1/2)||f||2+L^(1/2)||f_r||2.
```

Applying this to `D` and `D_r`, for `L=1/4`, gives the generic bound

```
R_A <= (sqrt(17)/2)||D||H2.
```

The actual retained admissible image is flat at first contact `r=1`:
`D(1)=D_r(1)=0`. Therefore a sharper anchored estimate is available whenever
an H2 norm is controlled:

```
D_r(r)=integral_1^r D_rr(s)ds,
D(r)=integral_1^r (r-s)D_rr(s)ds,
||D_r||infinity <= (1/2)||D_rr||2,
||D||infinity <= (1/(8sqrt(3)))||D_rr||2,
R_A <= (1/2)||D_rr||2 <= (1/2)||D||H2.                        (14)
```

Consequently `||D||H2<=4e33` would be one sufficient extra assumption for
the already declared class cutoff. It is a background-calibrated geometric
coefficient regularity assumption. It is not automatically a physical
energy, a preparation mechanism, or a justified selection law.

## Relation to the accompanying direct H1 result

An H1 bound can control initial data in weaker norms even while (13) grows.
The accompanying checkpoint addresses this alternative with a separate direct
H1 energy and boundary-history estimate for the smoothed observables. That
positive estimate is not proved or certified by this image component; its
own proof, analytic premises and exact gates must be assessed separately.
The older supremum resolvent involving `W_r` alone does not establish the
alternative. This construction supplies no broader no-go theorem, no new
observable sign result of its own, and no uniform finite-epsilon family
conclusion.

`verify_algebra.py` is a standard-library, exact-rational check of the
operator/derivative coefficients, bump polynomials and center jets,
continuous bump constants, localization exponents, support and H2 constants.
Its canonical certificate is `references/NORMAL.json` (byte-identical to `references/OPTIMIZED.json`). These checks make no actual
background query, numerical bump query, oracle call or PDE solve.
