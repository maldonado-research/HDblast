# A direct characteristic estimate for the odd temporal scalar response

4 November 2026, America/Los_Angeles. Conditional first-linear-coefficient
result for the same registered packet, `ell=10^-6`.

This derivation works directly with the actual coupled difference, rather
than estimating an actual-minus-free comparator by the old displacement
packet's generic error. It establishes

```
sup_0<=t<=3ell |S_current(-t)-S_current(t)| < 1/40,
sup_(49/25)ell<=t<=3ell |S_current(-t)-S_current(t)| < 1/40000.
```

The second estimate includes the entire post-pulse interval `[2ell,3ell]`.
These are bounds on a difference of physical first scalar coefficients in
model units, not on either individual amplitude or on a nonlinear signal.

## 1. Retained problem, reflection and coarse constants

Retain the admitted branch, smooth linear Wentzell development, finite
propagation, compatible initial jets and physical reconstruction established
in the selected contact/comparison proofs. The registered data are
`w(0,z)=qell(z)=Q(z/ell)` and `w_t(0,z)=exp(-z)qell(z)`, where

```
Q(x)=exp(4-1/[(x-1)(2-x)]) on 1<x<2, zero elsewhere.
```

The complete equation and boundary law are even in time. Thus the actual
past minus future scalar response at matched times is reconstructed from a
single forward difference `u` on `0<=t<=T=3ell`, satisfying

```
u_tt-u_zz+W(t,z)u=0,
-u_z+mu*u_tt+beta(t)u=0 at z=0,
u(0,z)=0,  u_t(0,z)=h(z)=-2exp(-z)qell(z).
```

The factor two in `h` is essential. The past and future observation clocks
are `-rho_b*t` and `+rho_b*t`, anchored on the same section. The normal
derivative is spatial, so there is no additional reflection minus sign.

The selected contact enclosures and their even extension supply deliberately
coarse safe constants on the support:

```
mu > m=1/1000,  |W|<M=103000,
|beta|<B=300,  |gamma|<Gamma=300,
gamma=C_b/2,  |G|=|3/(rho_b^3*p_b)|<Gstar=3/250.
```

For example, the contact report directly encloses `mu`, `beta` and `gamma`;
its physical denominator bounds give `|G|<=25000/129^3<3/250`.
The original fields remain inside `z<=5ell`, within the verified collar
`z<=8ell`; the auxiliary extension has no effect there. The comparison
proof supplies uniqueness and the exact same-state identification. None of
these analytic premises is proved by authenticating the selected files.

## 2. The zero-displacement free Dirichlet term is small

Let `h_odd` be the odd extension of `h` to the full line. On the half-line,
the zero-boundary free wave is

```
f(t,z)=(1/2) integral_[z-t,z+t] h_odd(r) dr.
```

The function `h` is nonpositive on the positive half-line and supported in
`[ell,2ell]`; its odd extension is nonnegative on the negative half-line.
Each sign's total integral has magnitude at most `2ell`, because
`0<=Q<=1` and its support has length one in scaled coordinates. Therefore
an integral over any interval is the difference of two numbers in
`[0,2ell]`. Consequently `|f|<=ell`, without estimating `Q'` or sampling it.
Also

```
integral_0^T |h(t)| dt <= 2ell.
```

Let `b(t)=u(t,0)`, `v=b'`, and let `A` and `V` be the suprema of `|u|`
and `|v|` over the compact supported slab and its boundary. Smoothness and
support make both finite. The Dirichlet boundary lift is `b(t-z)` for
`t>=z`, zero otherwise. The Dirichlet Duhamel term is bounded by
`M*A*t^2/2`. Hence

```
(1-r)A - T V <= ell,    r=M*T^2/2.
```

## 3. Keep the complete boundary law

The incoming characteristic identity is

```
u_t(t,0)+u_z(t,0)=h(t)-I(t),
I(t)=integral_0^t W(s,t-s)u(s,t-s) ds,
|I(t)|<=M*t*A.
```

Combining this identity with the full Wentzell law gives exactly

```
mu*v'+v+beta(t)b=h-I,    b(0)=v(0)=0.
```

The integrating factor for `mu*v'+v` is positive. Its kernel is
`mu^-1 exp(-(t-s)/mu)`, bounded by `1/m`; no boundary spring or endpoint
inertia is omitted. Since `|b(s)|<=sV`, the support-integrated input bound
gives

```
-(r/m)A+(1-s)V <= 2ell/m,    s=B*T^2/(2m).
```

Put `d=(1-r)(1-s)-T*r/m`. Exact arithmetic below gives `d>0` and
`r,s<1`; the inverse of this two-by-two majorant matrix has nonnegative
entries. Thus the two inequalities imply

```
A <= Astar=ell*(1-s+2T/m)/d,
V <= Vstar=ell*(2-r)/(m*d).
```

This is a finite-supremum estimate for the retained smooth actual solution,
not a new existence theorem or a numerical PDE solve.

## 4. A pointwise physical trace, not an energy-to-normal-trace inference

The exact physical reconstruction uses

```
D_delta(t)=u_z(t,0)-gamma*b(t),
S_delta(t)=G*cosh(t)^(-3/2)*D_delta(t)
         =S_current(-t)-S_current(+t).
```

Incoming transport, rather than an L2 trace argument, gives

```
|D_delta(t)| <= 2exp(-t)Q(t/ell)+C,
C=M*T*Astar+(1+Gamma*T)*Vstar.
```

Using `cosh(t)^(-3/2)<=1`, `Q<=1`, and the exact constants yields

```
r=927/2000000000, s=27/20000000,
d=39999927404405029/40000000000000000,
Astar=40239946000/39999927404405029,
Vstar=79999981460000/39999927404405029,
C=80084415586628/39999927404405029,
Gstar*(2+C)=120119908836595029/4999990925550628625 < 1/40.
```

For `x>=49/25` on the remaining support `[49/25,2)`, `Q` decreases and
`Q(x)<=exp(-529/24)<exp(-22)`. The elementary exponential series gives
`exp(2)>109/15>7`; hence `Q(x)<7^-11<10^-9`. For `x>=2`, `Q=0`.
It follows uniformly on `49/25<=x<=3` that

```
|S_delta(ell*x)|
 < Gstar*(2/10^9+C)
 =120126743379724213215087/4999990925550628625000000000
 <1/40000.
```

All decimal summaries are non-decision aids. Both bounds are in the same
model units and normalization as the retained actual future scalar report.

## 5. Conditional fine past zeros and negative post-pulse coefficient

The older actual future waveform report concerns the same `ell=10^-6`
packet, background, nonzero registered velocity, full boundary law, physical
embedding and clock. Its `ell_S` entries enclose `ell*S`; divide endpoints
by `ell` before widening by the new difference bounds. The actual-H
enclosures already include comparator error. A comparator alone, or the
different `ell=10^-10` packet, would not supply this identification.

Widening the two central witnesses by `1/40` preserves their opposite signs;
widening the two late witnesses by `1/40000` also preserves their signs.
The resulting safe strict physical-S bounds are

```
S(-1.49998ell)<-3/2,    S(-1.49999ell)>1,
S(-1.96ell)>1/1000,     S(-1.97ell)<-29/100000.
```

Continuity therefore proves at least one zero in each of the disjoint past
intervals

```
(-1.49999ell,-1.49998ell),    (-1.97ell,-1.96ell).
```

The first has width `ell/100000`, narrower than the previous selected-cell
past bracket `ell/2048` by the exact factor `3125/64` (about 48.8). Neither
statement proves uniqueness or counts all zeros.

For the retained whole future tail `[2ell,3ell]`, multiply its actual-H
interval by the negative gain interval and the enclosure
`1-T^2<=cosh(t)^(-3/2)<=1`. Widening the resulting physical-S interval by
`1/40000` gives the safe strict **past** tail bound

```
-607/10^6 < S(-ell*x) < -296/10^6,    2<=x<=3.
```

This sign concerns the quadrupole coefficient. The angular factor changes
sign and has zero mean; no globally signed scalar field, radiation, physical
energy or entropy statement follows. The old actual waveform/tail theorem
remains a conditional retained premise; its original campaign is not rerun.

## 6. A past quadratic response coefficient with the cross term retained

Let `H_+(x)=ell*D_current(ell*x)` and
`H_-(x)=ell*D_current(-ell*x)`. Then
`delta H=H_--H_+=ell*D_delta(ell*x)`. The zero-displacement input obeys
`||h(ell*.)||_L2(0,3)<=2`, and the remaining trace history is uniformly
bounded by `C`. Consequently

```
||delta H||_L2(0,3) <= ell*(2+sqrt(3)*C)
                    < ell*(2+2C) < 3ell.
```

The retained actual future squared integral `I_+=integral_0^3 H_+^2 dx`
is less than 9, so `||H_+||<3`. Expanding the square, with the cross term
retained, gives

```
|I_--I_+| <= 2||H_+||*||delta H||+||delta H||^2
            < 18ell+9ell^2 = 18000009/10^12.
```

The background proper-volume integral of the squared first scalar has the
same even-cosh cancellation on the past interval as on the future interval.
Changing `t` to `-t` reverses the limits of the past integral, not its sign.
Thus the corresponding positive quadratic response functional is

```
C_-=K/ell * I_-,    K=3*pi^2/(rho_b^2*p_b^2)>0.
```

It is not physical injected or Noether energy. The old report encloses the
same background `K` and future `C_+`. Using its `K_upper` gives

```
|C_--C_+| < (K_upper/ell)*(18ell+9ell^2)
           =1021088587031956677959217/500000000000000000.
```

Widening the retained actual future `C_+` interval by this exact ceiling
still proves the convenient strict enclosure

```
6.06*10^11 < C_- < 7.44*10^11.
```

This defines and bounds the background-integrated square of the admitted
first scalar signal. Interpreting it as a leading coefficient for a common
exact nonlinear family on the entire past interval would require those
separate family/development hypotheses and is not established here.

## 7. Limits and verification

The new estimates are a direct zero-displacement application of established
Dirichlet Duhamel, characteristic transport and positive scalar integrating
factor methods. No new general principle or external novelty claim is made.
The selected files, exact rational majorants, witness transports and tail
product are checked by `verify.py`, with identical normal/optimized reports.
This arithmetic supplements the continuous argument and does not prove the
retained branch, reduction, smooth-development or old waveform theorems.

The result gives no nonlinear finite-amplitude interval, particular
second-order past baseline, physical incoming source, SI calibration,
observational fit or higher-dimensional Big Bang evidence. Earlier
conditional formal complete future second Taylor constructions are not
asserted absent. Internal independent review is not external peer review or
proof-assistant certification.
