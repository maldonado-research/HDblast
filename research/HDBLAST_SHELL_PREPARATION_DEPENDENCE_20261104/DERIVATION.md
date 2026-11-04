# Conditional preparation-dependence theorem and an elementary sign certificate

4 November 2026. Research checkpoint with two separately checked sufficient methods.
The [published 1 November prepared-family capsule](../HDBLAST_PREPARED_SHELL_VOLUME_20261101/) is preserved unchanged. The theorem
below is conditional on precisely the retained background, weighted constraint
inverse, exact-family development, linear uniqueness, and anchored observable
premises used by that published result. The new calculation is a homogeneous
five-mean difference calculation. It is not an independent observation of the
Einstein–scalar PDE and does not select a physical preparation mechanism.

## 1. Fixed objects and hypotheses

Fix the exact selected background, first packet and parameter normalization

```
x_e=x0+e v_seed+e^2 x2/2+o(e^2),
g_e=g0+e q+e^2 k/2+o(e^2).
```

Retain the same first-field solution, its complete quadratic constraint and
evolution sources, initial shell anchor, annulus, and prescribed development
slab. Write `p=phi0_y`, `a=rho_y/rho`, `dy=-rho dz`, `r=z/ell`, and
`ell=10^-6`. The exact packet interval is `1<=r<=2`; the target is
`T=5/4`, `tau_star=rho_b ell T`, where `rho_b` is the exact fixed-background
shell radius, not a rounded value selected from an enclosure.

For each fixed finite preparation parameter, assume the same weighted IFT
isomorphism and nonlinear mapping assumptions as the published preparation,
and the same conditional C2 parameter dependence in the retained lower C4
spacetime topology. Retain a common intrinsic normal tube through the target
after family-dependent shrinkage in `e`. The mean linear mixed problem retains
all five metric/scalar channels, the physical Robin rows and normal-gauge row.
Its uniqueness and coefficient bound on

```
D_T={0<=x<=T, 0<=r<=T-x}
```

are the published ones, including `q=1201/3200000<1` in the Gram norm
`H=diag(1,2,1,3,1)` and
`X(U)=max(sup|U|_H/3,sup|U_x|_H,sup|U_r|_H)`.

These are hypotheses inherited from the fixed-preparation theorem; the new
finite algebra and rational bounds do not prove those analytic hypotheses.

## 2. Exact admissible preparation freedom

Let `h` be smooth on the radial line, supported in the closed packet interval,
flat at both endpoints, and satisfy

```
integral_(yL)^(yU) rho^2 h dy=0,
yL=y(2ell), yU=y(ell).
```

For every fixed finite real `lambda`, define

```
f_lambda=f+lambda h,
zeta=-(rho^4 p^2)^-1 integral_(yL)^y rho^4 C[h] ds,
xi_lambda=xi+lambda zeta,
C[h]=-3(h_yy+3a h_y)/2-Vbulk h.
```

Leave the canonical momentum mean `diag(0,-QM/(6a),-QM/(6a),-QM/(6a))`
and both scalar canonical sectors unchanged. The preparation difference is

```
B_lambda-B_base=lambda(H_h,0,0,0),
H_h=h gamma0+Lie_(zeta partial_y)gamma0.
```

The retained background identity gives the exact flux formula

```
rho^4 C[h]=-(3/2)[rho^4(h_y-a h)]_y-(9/2)rho^2 h.
```

Indeed, expanding its right side and using
`a_y+4a^2=3/rho^2-2Vbulk/3` gives the displayed definition of `C`.
Consequently the zero moment makes the primitive vanish on both exterior
sides. Flatness and the nonzero `p` bound imply smooth compact extension.
The metric constraint difference is exactly

```
P_g H_h=C[h]+rho^-4(rho^4 p^2 zeta)_y=0.
```

The momentum and scalar canonical variations are zero, so the whole second
constraint equation remains `P B_lambda=-P0 Q` with the identical full first
constraint Hessian `Q`. This construction adds a homogeneous second datum;
it neither changes nor removes a quadratic source.

Insert `e^2 B_lambda/2` in the same weighted IFT trial. Its derivative in the
correction variable at the background is unchanged. For fixed weights,

```
L U_lambda''(0)=-w^-2(Q+P B_lambda)
```

is independent of `lambda`. Thus injectivity gives the same full correction
second derivative for these fixed weights. Haar averaging and the retained
equivariance argument give zero second monopole for that correction, so the
exact completed family's second mean is exactly `B_lambda`.

The initial shell collar, first tangent and second TT initial pair persist.
The first field and all its quadratic sources therefore persist by the
retained linear uniqueness. The support and slab arguments apply to each
chosen finite `lambda` after qualitative shrinkage in `e`.

No uniform radius as `|lambda|` grows is asserted. In particular, `lambda`
is fixed independently of `e`; an unbounded choice `lambda(e)` would not be
covered. The full finite-`e` family and its higher derivatives need not be
affine in `lambda`.

## 3. The freedom is physical and is not the conditioning gauge

The class can be parameterized without numerical moment cancellation:

```
h=rho^-2 chi_y,
chi(y)=integral_(yL)^y rho^2 h ds.
```

The zero moment is equivalent to a smooth compact flat `chi`. The local
formula following from the flux identity is

```
zeta=3(h_y-a h)/(2p^2)+9chi/(2rho^4p^2),
s=-p zeta.
```

A nonzero `h` does not yield an anchored pure second diffeomorphism of the
background. On the time-symmetric initial section a normal generator changes
the induced spatial metric and scalar by multiples of `K0` and `Pi0`, which
both vanish. Haar averaging any proposed tangential generator leaves a radial
one, since there is no invariant angular one-form. Its scalar variation is
`p X^y`. The raw scalar difference is zero and `p` is nonzero, so this radial
generator vanishes. It could therefore account for the metric difference
only if `H_h=0`. Its radial and angular equations would be

```
h+2zeta_y=0, h+2a zeta=0.
```

They imply `zeta_y=a zeta`, and compact support implies `zeta=0`, then `h=0`.
This proves nontriviality in the anchored pure second gauge quotient with
the common first representative fixed. It does not state that every physical
direction is detected by every observable.

For each individual preparation, the published pure second harmonic vector
with initial `Z^y=-xi_lambda`, `Z^t=0`, `Z^y_t=0`, `Z^t_t=a xi_lambda`
changes only its representative and transports the shell clock/anchor.
Its scalar pullback is essential. Subtracting the two conditioned adapters
gives, in `(K,L,N,J,w)` order,

```
W(0,r)=(h,0,0,h,s),
W_x(0,r)=(0,ell rho(-h_y+a h),0,0,0),
W_r(0,r)=(h_r,0,0,h_r,s_r).
```

The lapse-square term in the absolute adapter cancels. The mixed metric
velocity is retained; the scalar coordinate velocity is zero. The gauge
conditioning alone cannot change any anchored proper observable.

## 4. Exact affine observable relations

The actual forced first-field source, operator, and shell rows are identical
for the two preparations. Retained linear uniqueness therefore gives

```
U_lambda=U_base+lambda W,
W_xx-W_rr=A0 W_x+B0 W_r+D0 W,
L|shell=0, PN W_r=M0 PN W.
```

This is the full coupled homogeneous mean problem, including metric–scalar,
connection and Robin feedback. Set `m(x)=W_w(x,0)`. The normalized proper
scalar mean's first clock and measure corrections are identical, and the
complete first-only quadratic Gauss term `Q1` is identical. They cancel
exactly when subtracting the two responses. Thus

```
M2_lambda(x)=M2_base(x)+lambda m(x),
v2_lambda(x)=v2_base(x)+lambda nu(x),
Theta2_lambda(x)=Theta2_base(x)+lambda theta(x),
nu_tautau+4H nu_tau-2nu/A^2=(F'_b/2)m,
nu(0)=nu_tau(0)=0, theta=nu_tau.
```

The cancellation is a difference identity only. Every absolute response
continues to contain the actual forcing, clock correction, measure correction,
and every `Q1` term. No additional factor `1/2` multiplies the preparation
slopes: these observables are already second parameter derivatives.

At the target the three continuous linear slopes are

```
L_M[h]=m(T),
L_v[h]=(rho_b ell)^2(F'_b/2) integral_1^T k(T,x)m(x)dx,
L_Theta[h]=rho_b ell(F'_b/2) integral_1^T k_T(T,x)m(x)dx,
k(T,x)=sinh(2ell(T-x))/(2ell)*cosh(ell x)^2/cosh(ell T)^2.
```

The exact constant `F'_b` is negative on the authenticated background.
For `v=ell T`, `b=ell x`,

```
k_T=cosh(b)^2 cosh(v-2b)/cosh(v)^3.
```

Since `0<=b<=v` and `|v-2b|<=v`, it follows that `0<k_T<=1`.
Also `0<=k<=T-x`. On the retained target interval,

```
1-eta <= k/(T-x) <=1,
1-eta <= k_T <=1, eta=10^-11,
```

where the first ratio is extended continuously at `x=T`. For example the
lower ratios are at least `(1+v^2)^-2` and `(1+v^2)^-3`, respectively, and
`3v^2<eta`. These are properties of the exact kernels, not numerical samples.

## 5. A rigorous homogeneous error bound

Let `B` be the free reflected wave of the displayed initial data and let
`D0,Dr,Dv` bound their three Gram norms. Only the initial segment `0<=r<=T`
is needed for the backward-closed cone; the exact data vanish on `r<=1`.
Then

```
N=max((D0+T Dv)/3,Dr+Dv) >= X(B),
E=qN/(1-q) >= X(W-B).
```

All coupled feedback is bounded through the resolvent `1/(1-q)`. There is
no source-history term in this homogeneous difference estimate. Both `W`
and `B` vanish in `x+r<=1`, giving the stronger boundary estimate

```
|m(x)-s(x)|<=(x-1)_+ E.
```

Consequently, for `Delta=1/4`,

```
|L_M-s(T)| <= Delta E,
|integral k(m-s)| <= Delta^3 E/6,
|integral k_T(m-s)| <= Delta^2 E/2.
```

The last bound uses the exact `k_T<=1`; an extra `eta` on this term is
unnecessary. A free-history sign alone is insufficient, but a strict lower
bound after these errors certifies the full coupled slope.

## 6. One explicit direction and its derivative checks

Choose the single fixed direction

```
Q(r)=exp(4-1/[(r-1)(2-r)]) for 1<r<2, zero otherwise,
chi=ell rho_b^3 Q(r), h=-(rho_b/rho)^3 Q_r(r).
```

This direction has exact zero moment, independent of quadrature. Write

```
R=(rho_b/rho)^3, b=ell rho a,
c=-3rho_b^3/(2ell rho^4 p)>0, d=ell rho p_y/p.
```

Then the exact local profiles are

```
h=-R Q',
h_r=-R(Q''+3bQ'),
V_L=-R(Q''+4bQ'),
s=c(Q''+4bQ'+3ell^2 Q),
s_r=c[Q'''+4bQ''+3ell^2Q'
      +(4b+d)(Q''+4bQ'+3ell^2Q)
      -4ell^2rho^2(a^2+a_y)Q'].
```

Primes here are dimensionless `r` derivatives. The identities use
`rho_r=-ell a rho^2`, `R_r=3bR`, `c_r=c(4b+d)`, and
`b_r=-ell^2rho^2(a^2+a_y)`. Independent symbolic differentiation confirms
the formulas. They imply the initial Gram norm bounds

```
|W0|_H=sqrt(4h^2+s^2),
|W0_r|_H=sqrt(4h_r^2+s_r^2),
|W0_x|_H=sqrt(2)|V_L|.
```

For `t=r-1`, `D=t(1-t)`, and `1<r<=T`, the exact bump derivatives are

```
Q'=Q(1-2t)/D^2,
Q''=Q(1-6D+6D^2)/D^4,
Q'''=Q(1-2t)(1-10D+24D^2-12D^3)/D^6.
```

Here `0<D<=3/16`. The second numerator is decreasing and at least
`11/128>0`. The third numerator decreases from `1` to `-113/1024`, hence
its absolute value is at most `1`. Thus `Q,Q',Q''` and `s` are strictly
positive throughout this open causal history. This does not require the
actual `m(x)` to have a sign arbitrarily close to first contact.

Put `w=1/D>=16/3`. For `j=2,4`, `e^-w w^j` is decreasing on this interval.
For `j=6` its maximum occurs at `w=6`. Rational Taylor bounds prove

```
1/4 < exp(-4/3) < 4/15,
exp(-2)6^6 <6400.
```

They give the conservative, continuous bounds

```
0<=Q<=1/3, |Q'|<=8, |Q''|<=216, |Q'''|<=6400
```

on the whole causal initial segment, including the flat endpoint by smooth
extension. No mesh sampling or unresolved tail is used. In addition,

```
Q(T)=exp(-4/3), Q'(T)=(128/9)exp(-4/3),
Q''(T)=(5632/81)exp(-4/3).
```

If `c_min` is any uniform lower bound for `c`, positivity and integration
by parts give

```
integral_1^T (T-x)s(x)dx >= c_min Q(T)>c_min/4,
integral_1^T s(x)dx >= c_min Q'(T)>32c_min/9,
s(T)>1408c_min/81.
```

These inequalities avoid cancellation in the free moments.

## 7. Independent coarse certificate

[theory/produce.py](theory/produce.py) reads only the two byte-pinned published JSON inputs
named in its report. It checks that their exact rational intervals imply

```
129<=rho,rho_b<=131,
3/25000<=-p<=1/7000,
111/1000<=a<=14/125,
|p_y|<=1/4000, |a_y|<=1/10000,
1/1200<=-F'_b<=1/1000.
```

The derivative bounds are global annulus bounds; their first-order Taylor
jet entries are the actual first derivatives. Using these caps in the exact
profile formulas gives rational upper bounds for `D0,Dr,Dv`. The script uses
the conservative triangle inequalities `sqrt(4h^2+s^2)<=2|h|+|s|` and
`sqrt(2)|V_L|<=2|V_L|`.

The exact lower sign gates are

```
m(T) > 1408c_min/81 - Delta E,
integral k m > (1-eta)c_min/4 - Delta^3 E/6,
integral k_T m > (1-eta)32c_min/9 - Delta^2 E/2.
```

All three are positive. Safe shortened outward slope enclosures from the
exact JSON report are

| Slope | Outward enclosure |
|---|---|
| `L_M` | `[1.269e9,2.198e10]` |
| `L_v` | `[-0.005883,-0.0001282]` |
| `L_Theta` | `[-359.5,-14.21]` |

Thus the coupled volume and expansion response functionals are nonzero on
this same admissible physical direction and both are negative. These broad
bounds are an alternative analytic sufficient method to the tighter new
whole-cell direction report. They reuse the same background and coefficient
premises, so their agreement is not an independent PDE observation.

## 8. Deterministic simultaneous reversal and qualitative exact-family meaning

Suppose the base enclosures are `[b_v^-,b_v^+]` and `[b_t^-,b_t^+]`, both
strictly negative, and certified slope enclosures are `[l_v^-,l_v^+]` and
`[l_t^-,l_t^+]`, also strictly negative. Define the fixed algebraic choice

```
L = least 10^n, n in nonnegative integers, strictly greater than
    max((-b_v^-)/(-l_v^+),(-b_t^-)/(-l_t^+)),
lambda=-L.
```

Then

```
v2_lambda >= b_v^- + L(-l_v^+) >0,
Theta2_lambda >= b_t^- + L(-l_t^+) >0.
```

This rule is postprocessing of certified intervals for one previously fixed
direction. It changes no shape, target, input, error method, or sign test.
The affine relations are exact for every finite chosen `lambda`.

Using the independent coarse exact intervals selects `L=10^21` and proves

```
v2_(-10^21) >7e16,
Theta2_(-10^21) >5e21.
```

The complete exact rational bounds are in [theory/references/NORMAL.json](theory/references/NORMAL.json); the
quoted lower bounds are shortened outward. The published base family's
coefficients remain negative. Therefore their signs are not determined by
the retained first packet, constraints, shell collar and TT data alone.

For each of these two fixed finite preparations, the first averaged volume
and expansion variations vanish. Hence C2 parameter dependence gives

```
V_e/V0=1+(e^2/2)v2+o(e^2),
Theta_e=Theta0+(e^2/2)Theta2+o(e^2).
```

Since `V0>0` and `Theta0>0` at this same proper-clock target, the new positive
coefficients imply, for all sufficiently small nonzero `e` of either sign,

```
V_e>V0>0, Theta_e>Theta0>0.
```

The negative base coefficients instead imply smaller positive volume and
slower positive expansion. Each assertion has a qualitative, uncomputed,
family-dependent parameter interval; a common interval for these two fixed
families can be obtained by taking the minimum. Neither assertion supplies
an effective physical amplitude, a contracting branch, a signed finite-time
history, or a finite-amplitude cosmological mechanism.

Separate nonzero slopes allow arbitrary individual real second coefficients
as `lambda` varies. This single line supplies simultaneous reversal because
its two slopes have the same sign. It does not establish rank two or arbitrary
independent pairs `(v2,Theta2)`. Many directions supported outside the target
causal region are invisible there.

The predictive conclusion is precise: without an additional physical rule
selecting the compact second preparation, this first-packet specification
does not select the signs of these two second observables. This does not
disprove every blast model or change the validity of the original theorem
for its expressly chosen preparation. Comparisons with other models require
their own equations, initial conditions and observables.


## 9. Detailed whole-cell certificate and outward displays

A separately registered interval calculation fixes the same direction and target
before evaluation. It uses the authenticated N1024 background cover from the
preserved checkpoint, restricted to the causal initial segment: 241 continuous
cells (including the grouped flat tail) and one target. Fresh bump derivatives
through order three produce the new five initial components, their scaled
velocities and radial derivatives. Earlier prepared profile bounds are not
reused as bounds for this different profile. No new background or evolved-field
oracle query is made.

The detailed certificate has safe upper displays `N_h<4.5253e10` and
`E_h<1.6991e7`. Whole-cell positive weights and the complete feedback error give
kernel histories in `[2.1313e7,2.5493e7]` for volume and
`[3.0511e8,3.6022e8]` for expansion. Multiplication by the exact proper-time
factors and enclosed negative `F'_b/2` gives:

| Quantity | Safe outward enclosure |
|---|---|
| Scalar-mean slope `L_M[h]` | `[1.5393e9,1.7068e9]` |
| Relative-volume slope `L_v[h]` | `[-2.1094e-4,-1.5795e-4]` |
| Expansion slope `L_Theta[h]` | `[-22.9410,-17.4038]` |
| New `v2`, `lambda=-10^21` | `[1.0343e17,1.8442e17]` |
| New `Theta2`, `lambda=-10^21` | `[8.3831e21,1.8496e22]` |

The deterministic finite-lambda rule from section 8 was declared as algebraic
postprocessing before this evaluation. The sharper method also selects
`lambda=-10^21`; it changes no profile, partition, target or error gate.
Independent Fraction code reconstructs the whole-cell profiles, all five-vector
norms, histories, errors, slopes and witness from raw boxes without importing the
numerical producer. The elementary coarse route independently supplies weaker
sufficient bounds with the same signs, without quadrature. Both routes share the
retained analytic and background-enclosure premises; neither independently
solves their PDE or measures a physical state.

The positive scalar-mean slope is a sensitivity coefficient, not a claim that
the new absolute scalar mean reverses sign for this finite lambda. Volume and
expansion signs come from their complete time histories. Their strict bounds
apply at the one registered proper time and do not establish a signed numerical
time interval or a volume-acceleration sign.
