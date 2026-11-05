# A sharper odd scalar trace from velocity-only characteristic transport

4 November 2026, America/Los_Angeles. Conditional first-linear-field result
for the same registered timelike-shell packet. Standard characteristic,
positive-kernel, energy and interval-transfer arguments are used; no new
general method or cosmological discovery is asserted.

## 1. Retained problem, physical coefficient and conclusion

Retain the analytic branch, smooth linear Wentzell development, full dynamic
boundary law, reconstruction and uniqueness premises of the selected public
contact and past–future proofs. Their exact local copies are in `inputs/`.
For `ell=10^-6`, `T=3ell`, the registered master has displacement
`q(z)=Q(z/ell)` and velocity `h(z)=exp(-z)q(z)`, where

```
Q(x)=exp(4-1/((x-1)(2-x))) for 1<x<2, and zero otherwise.
```

It satisfies, in physical coordinate time,

```
w_tt-w_zz+W(t,z)w=0,
-w_z+mu*w_tt+beta(t)w=0 at z=0,
W(t,z)=V(z)-(35/4)tanh(t)^2,
beta(t)=beta0-(35/4)mu*tanh(t)^2,
beta0=gamma+(5/2)mu,   gamma=C_b/2.
```

The embedding-included first physical scalar coefficient is

```
Phi1(tau,Omega)=S(tau/rho_b)*Y(Omega),
Y=n1^2-n2^2,  tau=rho_b*t,
S(t)=G*cosh(t)^(-3/2)*D(t),
D(t)=w_z(t,0)-gamma*w(t,0),  G=3/(rho_b^3*p_b).
```

These statements use the same state, spatial normal, angular normalization
and anchored clock as the retained actual future waveform certificate. The
constant background shell scalar removes first-order clock/label amplitude
terms. The fields remain in the admitted collar. On the relevant causal
triangles the following deliberately coarse bounds follow from the selected
contact report and its underlying theorem:

```
mu>m=1/1000,  |W|<M=103000,
0<beta0<B=300,  |beta(t)|<B,
|gamma|<Gamma=300,
|G|<=Gstar=25000/129^3 < 3/250.
```

For example the raw positive mu/gamma intervals prove the first and boundary
bounds; the potential interval `[22153/4,102966]` lies inside `[-M,M]`.
The gain bound uses `rho_b>=129` and `|p_b|>=3/25000`. Evenness of the
complete problem supplies the positive elapsed-time continuation, rather
than silently enlarging the interval of the older past-only statement.

The new conclusions are

```
sup_{0<=t<=3ell}|S(-t)-S(t)| < 1/40,
sup_{(49/25)ell<=t<=3ell}|S(-t)-S(t)| < 1/40000.
```

Both are absolute differences of first coefficients in fixed model units.
They do not bound either individual amplitude by these numbers or supply a
relative error. Conditional on the retained actual future certificates,
these bounds recover two distinct past zero brackets and a signed past
post-pulse tail, as derived below.

## 2. Isolate a genuine velocity-only coupled solution

Let `w_plus` and `w_minus` have initial data `(q,+h)` and `(q,-h)` for the
same complete dynamic-boundary problem. Their smooth existence, compatibility
and uniqueness follow from the retained linear construction. Since W and
beta are even in time, the registered solution obeys

```
w_current(+t,z)=w_plus(t,z),
w_current(-t,z)=w_minus(t,z).
```

Set `psi=(w_plus-w_minus)/2`. This is an actual coupled solution with
initial data `(0,h)`, not a free approximation or a different background.
Thus `w_minus-w_plus=-2psi`. All initial shell corner jets vanish because
h is smooth and supported away from the shell. The physical spatial normal
does not reverse, and the cosh conjugation is even. Consequently

```
S(-t)-S(t)=-2G*cosh(t)^(-3/2)*D_psi(t),
D_psi=psi_z(t,0)-gamma*psi(t,0).
```

The factor two is essential. It is not a cancellation of second-order
sources or a comparison of two second-order preparations with a fixed first
tangent. The first coefficient and all its differences vanish at the shell
before the first front `t<=ell`. No nonlinear arrival time is inferred.

## 3. Velocity mass controls the full field and boundary speed

Write `y(t)=psi(t,0)`. The nonnegative initial velocity has

```
0<=h<=1,  integral_0^infinity h(r)dr <= L=ell.
```

The free Dirichlet solution `u_D` with data `(0,h)` is nonnegative. Its
retarded one-dimensional kernel has size at most one half, so
`|u_D|<=L/2<=L`. The actual Dirichlet representation is

```
psi(t,z)=u_D(t,z)+1_(t>=z)*y(t-z)+R[-Wpsi](t,z).
```

The half-line retarded kernel is positive and has total mass at most `t^2/2`.
All estimates below may be restricted to `t+z<=T`; those triangles stay
inside the validated collar and use the retained actual coefficients.
For `A=sup|psi|` there and `Vb=sup|y'|` on `[0,T]`,

```
A <= L+T*Vb+r*A,   r=M*T^2/2.
```

Incoming transport at the shell gives, exactly,

```
psi_t+psi_z=h(t)-I(t),
I(t)=integral_0^t W(s,t-s)*psi(s,t-s)ds,
mu*y''+y'=h-I-beta*y.
```

Integrate the last equation with its positive exponential kernel. Use
`integral_0^T h<=L`, `|I(t)|<=M*t*A`, and `|y(t)|<=t*Vb` to obtain

```
Vb <= L/m+(r/m)*A+s*Vb,   s=B*T^2/(2m).
```

The two inequalities have a positive inverse because

```
d=(1-r)*(1-s)-T*r/m > 0.
```

They give the explicit uniform bounds

```
A <= Astar=L*(1-s+T/m)/d,
Vb <= Vstar=L/(m*d).
```

Exact rational evaluation gives approximately `Astar=1.00300047e-6` and
`Vstar=.001000001815`. Decimals are explanatory only. The improvement
comes from the velocity mass L, not an assumed small displacement of the
registered original packet. The factor `1+T` appropriate to the older full
displacement-plus-velocity field must not be reused as the new source size.

## 4. Characteristic residual bound with all prior history retained

On an incoming ray, the velocity-only free Dirichlet field is

```
u_D(s,t-s)=(1/2)*integral_(t-2s)^t h_odd(r)dr >=0,
integral_0^t u_D(s,t-s)ds=(1/2)*integral_0^t r*h(r)dr.
```

Since `0<=h<=q`, its moment is bounded by
`integral_ell^(2ell) r dr=3ell^2/2`. Weighting by W therefore costs at most
`3Mell^2/4` on the complete interval, with no order-`Mt*q(t)` displacement
term. No spatial derivative bound or cancellation between two large
displacement fields is needed.

The delayed boundary lift and actual potential-scattering remainder obey

```
integral_0^t M*|y(2s-t)|*1_(2s>=t) ds <= M*Vstar*t^2/4,
integral_0^t M*|R[-Wpsi](s,t-s)| ds <= M^2*Astar*t^3/6.
```

Thus the global characteristic bound is

```
|I(t)| <= Istar,
Istar=3Mell^2/4+M*Vstar*T^2/4+M^2*Astar*T^3/6.
```

Let u be the free-bulk comparator with the same `(0,h)` data and the
constant complete boundary `-u_z+mu*u_tt+beta0*u=0`. For `b(t)=u(t,0)`,

```
mu*b''+b'+beta0*b=h,  b(0)=b'(0)=0.
```

Put `q_e=y-b`, `v_e=q_e'`. Subtracting the exact actual/comparator boundary
equations yields

```
mu*v_e'+v_e+beta0*q_e=-I+(35/4)mu*tanh(t)^2*y.
```

The actual boundary y, its coefficient discrepancy and the complete prior
forcing history all remain in this identity. The positive exponential
kernel, `|q_e(t)|<=t*sup|v_e|`, and `|y(t)|<=t*Vstar` imply

```
sup|v_e| <= Estar,
Estar=[T*Istar/m+(35/4)*Vstar*T^4/4]/(1-s),
sup|q_e|<=T*Estar.
```

Here `T*Istar` bounds the entire history because Istar is valid globally.
It is not a post-pulse bound substituted into an on-pulse convolution.
Incoming transport now gives the normal-trace identity

```
D_psi-D_u=-I-v_e-gamma*q_e,
|D_psi-D_u| <= Rstar=Istar+(1+Gamma*T)*Estar.
```

Exact arithmetic gives `Istar<7.749e-8`, `Estar<2.325e-10` and
`Rstar<7.772e-8`. These are characteristic trace estimates, not pointwise
normal traces inferred from an L2 energy norm.

## 5. Inertial comparator and the two explicit scalar ceilings

For the impulse response k,

```
mu*k''+k'+beta0*k=0,  k(0)=0,  k'(0)=1/mu.
```

Since `beta0>0`,
`[mu*(k')^2+beta0*k^2]'=-2*(k')^2` and its initial value is `1/mu`.
Therefore `|k'|<=1/mu<=1/m` and `|k(t)|<=t/m`. With `b=k*h`,

```
D_u(t)=h(t)-(k'+gamma*k)*h(t),
|(k'+gamma*k)*h(t)| <= Cmem=(1+Gamma*T)*ell/m.
```

This retains the local velocity input as well as the boundary-memory term.
Combining the comparator and actual-error estimates gives

```
|S(-t)-S(t)| <= 2Gstar*[h(t)+Cmem+Rstar].
```

Globally use `h<=1`. The derived upper bound is approximately
`.0233149976013`, and exact rational evaluation proves it is below `1/40`.

For `x>=49/25`, Q is decreasing until its flat support endpoint at 2, and

```
Q(49/25)=exp(-529/24)<exp(-22)<10^-9.
```

The last comparison needs no sampled exponential: `e>8/3` by its positive
Taylor series, and `(8/3)^22>10^9` is an exact integer comparison. Beyond
x=2 the bump vanishes. Hence `h(ell*x)<10^-9` for every `49/25<=x<=3`.
The localized derived upper bound is approximately `2.33144790525e-5`,
strictly below `1/40000` by exact arithmetic. The same bound holds pointwise
in angle for the first perturbation because `sup|Y|=1`.

## 6. Transfer of actual future witnesses: two distinct past zeros

The retained wider future certificate is for this precise packet, not the
different `ell=10^-10` tangent and not a free-comparator-only value. Its
`ell_S` interval is `ell*S`, so division by ell precedes widening.

At the two central witnesses, use the new uniform error `1/40`; at the
two late witnesses, use `1/40000`. The exact interval operations preserve

```
S(-(74999/50000)ell)<0,
S(-(149999/100000)ell)>0,
S(-(49/25)ell)>0,
S(-(197/100)ell)<0.
```

Continuity proves at least one zero in each disjoint open interval

```
(-149999ell/100000, -74999ell/50000),
(-197ell/100, -49ell/25).
```

Thus there are at least two distinct post-contact past amplitude zeros.
The central bracket has width `ell/100000`, which is `3125/64` times
narrower than the preceding `ell/2048` past bracket. This is a transfer of
the older fine future witness bracket, not a claim to improve its precision.
No uniqueness, total zero count or optimal localization follows. Signs of S
are coefficient signs; the full angular field contains the sign-changing Y.

## 7. Uniform signed past post-pulse tail

The retained old actual-H theorem is continuous on `2<=x<=3`. Its positive
interval, negative gain interval and `1-T^2<=cosh(t)^(-3/2)<=1` give a
future physical-S enclosure. Widen that enclosure by `1/40000`, which is
valid throughout this tail. Its new rational endpoints are approximately
`[-.000606648,-.000296086]`; the strict readable conclusion is

```
-1/1500 < S(-t) < -1/4000,    2ell<=t<=3ell.
```

This is a uniform nonzero sign for the actual conditional first scalar
coefficient over that past interval. It is not a globally signed angular
field, monopole, physical energy flux or finite-amplitude nonlinear tail.

## 8. A separate norm/cross-term transport

Define `H_plus(x)=ell*D_current(ell*x)` and
`H_minus(x)=ell*D_current(-ell*x)`. The difference is
`DeltaH=-2ell*D_psi(ell*x)` and vanishes for `x<=1`. The characteristic
estimate, not an L2-to-normal-trace inference, gives

```
||DeltaH||_L2(0,3)^2 <= 2*[2ell*(1+Cmem+Rstar)]^2 < (3e-6)^2.
```

The retained actual future integral has `I_plus=integral H_plus^2<9`, hence
`||H_plus||_2<3`. Keeping the full cross term yields

```
|I_minus-I_plus| < 6*(3e-6)+(3e-6)^2
                = 18000009/10^12 = 1.8000009e-5.
```

For the first-field quadratic scalar norm on the past proper-time window,
the background measure and conjugation are even. Since
`integral_(S3)Y^2 dOmega=pi^2/3`, the change from
`-3rho_b*ell<=tau<=0` to positive elapsed time has a positive integral,
and the same exact cancellation as in the future gives

```
C_minus=(K/ell)*I_minus,
K=3pi^2/(rho_b^2*p_b^2)>0.
```

Use the retained certified K and actual I_plus intervals, enlarged by the
explicit cross-term error. The resulting C_minus enclosure is approximately
`[6.0697590387e11,7.4275679288e11]`; exactly,

```
6.06e11 < C_minus < 7.44e11.
```

This is a quadratic norm of the admitted first field on a past interval.
It is not physical injected or Noether energy. Identifying it with the
epsilon-squared leading response of a finite-amplitude exact past family
requires the separate full past-family hypotheses and remainder control;
those are not supplied here.

## Scope, credit and reproduction

The analytic reflection, linear construction, characteristic representations,
kernel and norm identities are continuous mathematical arguments under the
stated retained premises. The old actual waveform and tail/integral interval
theorems remain conditional prior inputs; selected-file authentication does
not rerun or reprove their full campaigns. Exact Fraction guards supplement
the written proof and include every coefficient, trace conversion, interval
transfer and cross-term decision above. No old producer, numerical PDE,
background trajectory, bump sample or field/source oracle is invoked.

The previous project already contains conditional formal complete future
second Taylor constructions. They are not described as absent. The present
first-field results do not supply a quantitative particular past second-order
observable/baseline, a controlled nonlinear defect map, an effective
finite-amplitude interval, physical source selection, SI calibration,
radiation production or evidence for a higher-dimensional Big Bang.
Internal independent mathematical review is not external peer review or
formal proof-assistant certification; broader novelty/significance still
require specialist assessment.
