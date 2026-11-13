# Global radial bounds for the actual affine initial source

6 November 2026. This calculation bounds the actual canonical affine source
identified in [the prior source derivation](inputs/PRIOR_ACTUAL_SOURCE.md).
The model, branch and first seed are fixed by
[the registered binding](inputs/REGISTERED_EXACT_FAMILY_BINDING.md),
[the canonical normalization](inputs/CANONICAL_CONSTRAINTS.md) and
[the compact initial jet](inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md).
The immutable public input copies are identified by [SOURCE_PINS.json](SOURCE_PINS.json).

The smooth actual background and its all-point collar bounds remain analytic
premises. The derivative bounds below follow from these premises by a global
Taylor-majorant argument. The finite rational calculations are recomputed by
[derive_source_jets.py](derive_source_jets.py), following
[REPRODUCE.md](REPRODUCE.md); no saved numerical result is an input to the
majorant calculation.

## Fixed background and exact canonical source

Set `ell=10^-6`, `dy=-rho dz` and `eta=phi0-1`. The retained collar includes
the complete correction annulus `ell/2<z<5ell/2` and obeys

```
129 <= rho <= 131, 111/1000 <= a <= 112/1000,
-14/100000 <= p <= -12/100000,
-9/100000 <= eta <= -8/100000.
```

These are the boxes in [the actual background certificate](inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json).
All background derivatives are derived from the coupled equations

```
rho_y=a rho, a_y=-rho^(-2)-p²/3,
p_y=V'(eta)-4ap, eta_y=p,
V=-2/27+(14/9)eta²+(50/27)eta³-eta⁴/6-(4/9)eta⁵-(2/27)eta⁶.
```

For `ell<z<2ell` put

```
Q=exp(4-ell²/[(z-ell)(2ell-z)]),
alpha=(rho/rho_b)^(3/2) p/p_b, B_amp=alpha/rho²,
u=B_amp Q, w_velocity=exp(-z)u,
v=-3(u_y+2au)/p,
vt=-3((w_velocity)_y+2a w_velocity)/p,
t_m=-3w_velocity/(2rho), chi=vt/rho.
```

The bump is zero elsewhere. The actual first canonical datum uses metric
components `h_yy=-4uY`, `h_ab=2rho²uY kappa_ab`, scalar `phi1=vY`, ordinary
momentum `pi1/rho0=t_m Y partial_y tensor partial_y`, and scalar momentum
`varpi1/rho0=chi Y`. Here `rho0` denotes the fixed background volume density,
distinct from the warp factor `rho`; `Y=n1²-n2²`. Scalar and scalar-density
second data are prescribed to vanish in this affine curve.

The exact affine Hessian functions are

```
F=-12u u_yy-12u_y²-144a u u_y
  -(192a²+96a_y+168/rho²)u²
  -(v_y+4up)²-V''v²-(6w_velocity²+vt²)/rho²,
J=(18u²-v²)/rho²,
R=8u_y t_m+2chi(v_y+4up),
L=-8u t_m+2chi v.
```

They satisfy `D²H=F Y²+J|DY|²` and
`D²(M/rho0)=R Y²dy+L Y dY`. The required second source is their negative.
This uses the canonical density prescription of the cited source derivation.
Every source term is smooth, flat and supported in `K={ell<=z<=2ell}`.

For the round probability sphere define

```
H0=1,
H2=(n1²+n2²)/2-1/4,
H4=Y²-(n1²+n2²)/2+1/12.
```

The eigenvalues of `-Delta` are `lambda=0,8,24`, their squared norms are
`1,1/48,19/720`, and

```
Y²=1/6+H2+H4, |DY|²=4/3+4H2-4H4.
```

Thus the unnormalized mode coefficients of the negative Hessian are

| Degree | Hamiltonian `f_l` | Radial momentum `r_l` | Gradient momentum `g_l` |
| --- | --- | --- | --- |
| 0 | `-F/6-4J/3` | `-R/6` | absent |
| 2 | `-F-4J` | `-R` | `-L/2` |
| 4 | `-F+4J` | `-R` | `-L/2` |

For a unit probability harmonic, each coefficient is multiplied by
`||H_l||_prob<=1`. For momentum, the orthonormal radial components of the
source for a unit harmonic are

```
f_M=(r_l,sqrt(lambda) g_l/rho), lambda=8,24,
f_M=r0, lambda=0.
```

All vector norms below use these components.

## Global bump derivative bounds

Write `x=(z-ell)/ell` and `b_Q=1/[x(1-x)]` on `0<x<1`. Then

```
Q=exp(4-b_Q), b_Q>=4,
b_Q=1/x+1/(1-x),
|d_x^k b_Q|/k! <= x^(-k-1)+(1-x)^(-k-1) <= b_Q^(k+1).
```

The exponential chain rule gives the global bound

```
|d_x^n Q|/n! <= exp(4-b_Q) P_n(b_Q),
P_0=1,
P_n(B)=(1/n) sum_{k=1}^n k B^(k+1) P_(n-k)(B).
```

Indeed `P_n(B)` is the coefficient of `t^n` in
`exp(sum_{k>=1} B^(k+1)t^k)`. Its coefficients are nonnegative. For each
monomial `B^m`, the maximum of `B^m exp(4-B)` over `B>=4` is at `B=4` when
`m<=4`, and at `B=m` when `m>=4`. Since `e>2`, respective rational upper
bounds are

```
4^m,                  m<=4,
m^m/2^(m-4),          m>=4.
```

Termwise maximization gives exact rational bounds through order sixteen.
No bump sample is used. Each derivative also tends to zero at both support
endpoints, because `exp(-b_Q)` dominates every finite polynomial in `b_Q`.
Consequently the zero extension is smooth and the global bounds hold on the
closed support.

## Simultaneous coupled background jets

Fix an arbitrary center `y0` in the collar and set `t=(y-y0)/ell`. A
nonnegative Taylor coefficient vector for a quantity `h` bounds
`|d_t^n h(y0)|/n!`. Products use ordinary Cauchy products. Differentiation
replaces coefficient `h_n` by `(n+1)h_(n+1)`. The triangle inequality makes
these operations valid for signed quantities.

Introduce `b=1/rho`, `c=1/p`, `H=rho^(3/2)` and `E=exp(-z)`. Alongside the
original background equations use

```
b_y=-ab, c_y=-V'c²+4ac,
(B_amp)_y=(V'c-9a/2)B_amp,
H_y=(3a/2)H, E_y=bE.
```

The equation for `B_amp` follows by differentiating the actual alpha
formula. At every center the following zeroth coefficients are valid:

```
rho:131, a:112/1000, p:14/100000, eta:9/100000,
b:1/129, c:100000/12, B_amp:6/(5*129²), H:1500, E:1.
```

The reconstruction bound `4/5<alpha<6/5` follows from the same retained
boxes and the fixed normalization at the seed anchor. For example, its
squared endpoint comparisons are bounded by
`(129/131)^3(12/14)^2>(4/5)^2` and
`(131/129)^3(14/12)^2<(6/5)^2`. Also `131³<1500²` proves the displayed `H`
bound. The interval for `p` excludes zero, so the signed reciprocal `c`
and its differential equation are legitimate.

At order `n`, insert the coefficient bounds already known through order
`n` into these polynomial ODEs. Replace signs and polynomial coefficients
by their absolute values, multiply the right side by `ell`, and divide its
coefficient of order `n` by `n+1`. This computes the next coefficient for
every variable simultaneously. Induction proves the background Taylor
majorants through order sixteen. This procedure preserves the coupled
background equations rather than selecting unrelated derivative intervals.

The helper performs all arithmetic with rational numbers. After an
operation it may round a nonnegative coefficient upward to fourteen
significant decimal digits by exact integer division. The rounding check
ensures that this never decreases the coefficient. No floating-point
arithmetic enters these majorants.

Because `x_t=-b`, the zero-constant increment of `x` has coefficient bounds
`h_n=b_(n-1)/n` for `n>=1`. The coefficient of order `n` in

```
sum_{k=0}^{16} (sup |d_x^k Q|/k!) h(t)^k
```

bounds the normalized `t` derivative of `Q` for `n<=16`. This is the finite
Taylor chain rule; it requires no convergence radius. Source envelopes now
follow from the exact formulas above by products and `d_y=ell^-1 d_t`.

There is an explicit derivative margin. A source uses at most two
proper-radial derivatives of the bump, so source jets through order twelve
require bump and seed jets only through order fourteen. Twelve later
applications of `D_s=q d_y` require source jets only through twelve. For the
known product `u vt/rho`, derivative thirteen uses at most derivative
fourteen of `w_velocity`. Thus stored order `N=16` covers both the source
H12 calculation and the known quadratic profiles through thirteen. A
truncated series derivative can create an unknown top coefficient, but no
such coefficient enters any claimed bound.

## Physical measure and transformed source H12

Use the explicitly distinguished weight member

```
kappa=ell/100=10^-8,
r(z)=(z-ell/2)(5ell/2-z)/ell²=1-(x-1/2)²,
q=kappa²r², wloc=exp(-1/r), ds/dy=1/q.
```

On `K`, `3/4<=r<=1`. Since `r0<=1` and `|r_x(y0)|<=1`, the Taylor
coefficients of `r` are bounded by `1+h+h²`. For the reciprocal, the
coefficient recursion is

```
I0=4/3, I_n=(4/3) sum_{k=1}^n r_k I_(n-k).
```

The exponential coefficient recursion then bounds `exp(1/r)`, starting
from zeroth coefficient `4`. The inequality `exp(4/3)<4` is verified by a
rational Taylor sum and a geometric upper bound for the remaining tail.

The original radial measure is `rho³dy`, while `dy=q ds`. The unitary
multiplier for the weighted potential is therefore
`U=sqrt(q rho³) wloc`. In particular the transformed scalar and momentum
sources are

```
G_H=sqrt(q rho³) wloc^(-1) f_l,
G_M=sqrt(q rho³) wloc^(-1) (r_l,sqrt(lambda) g_l/rho).
```

Their common factor is `kappa r rho^(3/2) exp(1/r)`. The volume factor
`rho^(3/2)` is included in every source jet. For momentum the further
factor `1/rho` in the tangential component is also differentiated; the
rational bound `sqrt(lambda)<=5` applies to `lambda=8,24`.

Given a Taylor envelope `A`, one envelope for `D_s A` is obtained by a
normalized series derivative followed by multiplication by `q/ell`.
Iteration gives pointwise bounds through twelve for both transformed
sources. Only source coefficients through order twelve enter the resulting
zeroth coefficients. All source derivatives vanish outside `K`.

Since `|K|_dy<=131 ell` and `q>=kappa²(3/4)²` on `K`,

```
|K|_ds <= 131 ell/[kappa²(3/4)²] = 20960000000000/9,
sqrt(|K|_ds)<10^7,
||D_s^j G||_L2(ds) <= 10^7 sup_K |D_s^j G|.
```

For a vector, the sum of component bounds controls its Euclidean norm.
The rational recurrence yields

| Channel | `sum_{j=0}^{12} ||D_s^j G_H||_2` | `sum_{j=0}^{12} ||D_s^j G_M||_2` |
| --- | --- | --- |
| 0 | `<10^24` | `<10^17` |
| 2 | `<10^25` | `<10^17` |
| 4 | `<10^25` | `<10^17` |

These hold for the public coefficients and hence for unit probability
harmonics. The derivative l1 sum bounds the usual square-sum Sobolev norm.
Applying the triangle inequality to the three channels gives the convenient
aggregate ceilings

```
||G_H||_H12(ds) < 3*10^25,
||G_M||_H12(ds) < 3*10^17.
```

These are probability-sphere conventions. The transformation already
includes the full `rho³dy` measure; no additional warp-volume factor is
inserted.

## Whole-annulus momentum coefficient derivatives

For coefficient bounds use the full correction annulus. Here
`-1/2<x<3/2`, so `|r_x(y0)|<=2` and the geometric Taylor envelope is
`r_ann=1+2h+h²`. No reciprocal of `r` occurs in the following whole-annulus
estimates.

Define the distinct logarithmic derivatives

```
b0=D_s log(sqrt(q rho³) wloc)
  =kappa²[(1+r)r_y+(3/2)ar²],
alpha_phys=D_s log(kappa r wloc rho^(-3/2))
  =kappa²[(1+r)r_y-(3/2)ar²].
```

The two expressions have a common valid absolute Taylor envelope, despite
their different signs. For the matrix in normalized nonconstant momentum
components,

```
C=[[-2b0,0],
   [sqrt(2)qk,-sqrt(2)(b0+qa)],
   [sqrt(12)qa,-sqrt(12)qk/3],
   [0,2q sqrt(d_lambda)/rho]],
k=sqrt(lambda)/rho, d_lambda=2(lambda-3)/3,
```

use `sqrt(2)<=2`, `sqrt(12)<=4`, `sqrt(lambda)<=5` and
`sqrt(d_lambda)<=4`. A sum of absolute entry envelopes controls the operator
norm of each derivative. Applying `D_s` to these envelopes gives, uniformly
on the entire annulus,

```
||D_s^j C||op < 10^-10, 0<=j<=14.
```

The zeroth entry-sum envelope is `1.2403343896643*10^-11`; all subsequent
ones are smaller. The same induction bounds the derivatives `D_s^j` of
`q/rho`, `qa`, `q_y`, `b0` and `alpha_phys` by one through order fourteen.
These bounds use at most fifteen background/geometric coefficients and are
therefore covered by `N=16`.

## Known quadratic initial-conversion profiles

For the two known radial products

```
A_quad=u w_velocity/rho,
B_quad=u vt/rho,
```

the same series products give exact proper-y derivative caps through
thirteen. They imply

```
max_{0<=j<=13} ||d_y^j A_quad||infty <= 10^73,
max_{0<=j<=13} ||d_y^j B_quad||infty <= 10^84.
```

These products are distinct from the background amplitude coefficient
`B_amp`. Their derivative margin is established above; no higher bump
derivative is inferred from a truncated coefficient. The physical ADM
conversion and its angular and covariant derivative estimates are treated
separately.

All bounds in this note are finite source and coefficient estimates,
conditional on the retained smooth branch. The weighted inverse and
zero-extension conclusions of [the initial theorem](inputs/EXACT_INITIAL_THEOREM.md)
are separate analytic inputs. This note supplies neither an evolved-field
bound nor a finite-amplitude interval or nonlinear remainder.
