# Conditional actual first-master derivatives through order four

For the same retained smooth actual master `U(x,r)`, all five scaled mixed
partials of total order four satisfy the strict bound **342395.408** on
`0<=x<=3, r>=0, x+r<=8`. The same shell trace obeys
`|y''''|<3.235927`. These are continuous all-point majorants conditional on
the retained smooth solution, background, boundary and initial-data premises.
They do not complete the gauge fourth-jet or complete bulk-source derivative
certificate. No second-response approximation is changed by this calculation.

The new step avoids a new `W_zzz` input. It applies the existing reflected
second-derivative identity to the actual function `w=U_xx`, then recovers
the remaining radial derivatives with the actual wave equation.

## Retained equation and temporal potential

Let `ell=1e-6`, `X=3`, and

```
U_xx-U_rr=-f,       f=ell^2 W(ell*x,ell*r) U,
U(0,r)=Q(r),       U_x(0,r)=ell V(r), V=exp(-ell*r) Q(r),
U(x,0)=y(x),       H=U_r(x,0)-ell*gamma*y(x),
y''=(ell/mu)H-A*y,
A=ell^2[5/2-(35/4)tanh^2(ell*x)].
```

The published input proof retains
`W(t,z)=W_static(z)-(35/4)tanh^2(t)`. Therefore
`W_tr=W_ttr=0`, and `|W_tt|<=35/2`. Indeed, with `z=tanh^2(t)` in `[0,1]`,
`(1-z)(1-3z)=1-4z+3z^2` lies in `[-1/3,1]`. Thus
`|A''|<=(35/2)ell^4`. These are derivatives of the registered analytic
potential, not derivatives of independently selected interval boxes.

Four pinned retained files supply actual master order-zero through three,
boundary `y,y',y'',y'''`, `H`, `mu`, `gamma`, bump derivatives through six,
the physical potential constants, and the proof of the separable temporal
term. The constants used below are
`M=102000`, `Mz=30000000`, `Mt=21/400000`, `Mtt=35/2`.
`Uj` denotes the authenticated upper bound for every scaled mixed partial
of order `j`. `Frr` and `Fx0` are retained exact actual forcing majorants.
Subscript `0` denotes shell restriction; it never denotes a different source.

## A fourth shell derivative without differentiating H cells

The exact incoming identity is

```
H=Q'+ell*V-y'-ell*gamma*y-I,
I(x)=int_0^x f(s,x-s) ds.
I''=f_x(x,0)+f_r(x,0)+int_0^x f_rr(s,x-s) ds.
```

Both moving-endpoint terms are retained. Consequently

```
Fr0=ell^3*Mz*|y|+ell^2*M*U1,
H2=Q3+ell*(Q2+2ell*Q1+ell^2)+|y'''|
    +ell*|gamma|*|y''|+Fx0+Fr0+X*Frr,
Y4=(ell/mu)_max*H2+(35/2)ell^4*|y|
   +2*|A'|*|y'|+|A|*|y''|.
```

The actual initial collar is flat near the shell. Finite propagation gives
`y=0` for `x<1`, so all required corner values vanish. No derivative of the
4,353 independent `H` interval cells is taken. The result is
`H2<6665.267481` and `Y4<3.235927`.

## Interior fourth jets through a differentiated wave

Ordinary product differentiation gives the following actual source bounds:

```
Fxx  =ell^4*Mtt*U0+2ell^3*Mt*U1+ell^2*M*U2,
Fxr  =ell^3*Mt*U1+ell^3*Mz*U1+ell^2*M*U2,
Fxxr =ell^4*Mtt*U1+2ell^3*Mt*U2+ell^3*Mz*U2+ell^2*M*U3,
Fxx0 =ell^4*Mtt*|y|+2ell^3*Mt*|y'|+ell^2*M*|y''|.
```

The factors of `ell` convert physical potential derivatives to scaled
coordinates. The vanishing mixed derivatives `W_tr,W_ttr` are used explicitly.

The smooth actual `w=U_xx` obeys

```
w_xx-w_rr=-f_xx,
w(0,r)=Q''-f(0,r),
w_x(0,r)=ell*V''-f_x(0,r),
w(x,0)=y''(x).
```

Its free-wave second derivatives, including the causal boundary lift, are
bounded by

```
Cfree=Q4+Frr+ell*(Q3+3ell*Q2+3ell^2*Q1+ell^3)+Fxr+Y4.
```

In particular the initial `f_rr` and `f_xr` terms must be included.
The initial profiles and velocities of `w` vanish near `r=0`; `y''` and its
required corner jets vanish there too. Thus the reflected free expressions
have compatible joins. The forcing `f_xx` and its first radial derivative
are actual continuous functions under the retained smooth premise, and its
corner trace vanishes. The physical odd-extension second-jet formula, with
its explicit crossing term, then yields

```
C=Cfree+X*Fxxr+Fxx0,
|w_rr|,|w_xr|<=C,     |w_xx|<=C+Fxx.
```

The remaining mixed derivatives follow from `U_rr=U_xx+f`:

| Derivative | Exact majorant formula | Approximate value for display |
| --- | --- | ---: |
| `U_xxxx` | `C+Fxx` | 342395.407974504 |
| `U_xxxr` | `C` | 342395.407948391 |
| `U_xxrr` | `C` | 342395.407948391 |
| `U_xrrr` | `C+Fxr` | 342395.407974504 |
| `U_rrrr` | `C+Frr` | 342395.407974505 |

The JSON retains exact rational values. For physical `t,z` derivatives of
total order four, divide the scaled cap by `ell^4`. These raw master
derivatives are not physical curvature tensors or observable forces.

## The smaller remaining gauge input

The inherited gauge system is
`V_xx-V_rr=P V_x+Q V_r+R V+s`, with zero Cauchy data and time-independent
mixed Dirichlet/Neumann rows. To extend its existing differentiated
contraction to order four, it is enough to bound **pure time** third
coefficient/source derivatives `P_xxx,Q_xxx,R_xxx,s_xxx`, together with
the retained mixed derivatives through order two. It is unnecessary to
introduce every mixed third coefficient derivative just for this route.

For `J=V_xxx`, the forcing after collecting the same principal operator is

```
3 P_x (V_xx)_x + 3 Q_x (V_xx)_r + 3 R_x V_xx
+3 P_xx (V_x)_x + 3 Q_xx (V_x)_r + 3 R_xx V_x
+P_xxx V_x + Q_xxx V_r + R_xxx V + s_xxx.
```

The stabilized initial data are

```
J(0)=P0*s0+s_x0,
J_x(0)=s_rr0 + P0*(P0*s0+s_x0) + 2P_x0*s0
       +Q0*s_r0+R0*s0+s_xx0,
J_r(0)=P_r0*s0+P0*s_r0+s_xr0.
```

These terms are not zero merely because the original gauge Cauchy data
are zero. With the same positive resolvent, the usual free-wave C1 bound
and commutator budget can control `N(J)`. The original PDE differentiated
twice then recovers the other fourth mixed derivatives using retained
coefficient/source jets through order two and lower field jets. The new
master C4 bound supplies the highest raw master derivative needed by
`s_xxx`, but the gauge coefficient/source third-time ledger and its complete
source-map composition still require implementation and independent review.
No numerical gauge C4 or complete actual bulk `F_x` cap is asserted here.

## Reproduction and limits

```
python -I -S -B check_master_fourth.py --root /path/to/HDblast --output /tmp/master4-normal.json
python -I -O -S -B check_master_fourth.py --root /path/to/HDblast --output /tmp/master4-optimized.json
```

The checker performs 87 exact polynomial identities for reflected second
derivatives, moving-endpoint integrals and the temporal derivative bound.
These finite tests support the displayed analytic all-point argument; they
do not replace it. Three deliberately wrong formulas are rejected in both
execution modes. The normal and optimized uncorrupted reports are byte equal.

This closes a specific quantitative first-master input, while preserving
all frozen candidates and sources. It does not prove the full differentiated
bulk source, a second-response C2 error, a nonlinear family, a cosmological
prediction, journal readiness or a higher-dimensional Big Bang origin.
