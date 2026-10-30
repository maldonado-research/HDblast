# Independent review of the actual first-master fourth-jet bridge

The bridge is mathematically consistent under the retained smooth actual
master premise. It supplies continuous fourth derivatives of the actual
first master on `0<=x<=3,r>=0,x+r<=8`, with `ell=10^-6`. It does not
differentiate interval selections or any P1 numerical profile. It does not
yet bound fourth derivatives of the reconstructed common gauge or the
complete differentiated quadratic forcing of the second response.

The review checker imports no producer code. It authenticates the previous
master and bump receipts, the source PDE and potential-separation proof
texts, and the producer's final `MASTER_FOURTH_READY` report. It independently
generates the product-rule terms, checks 49 incoming-integral polynomial
identities, and reconstructs the rational constants and all five mixed
fourth-derivative bounds. Normal and optimized output agree byte for byte.

## The differentiated equation closes without a spatial third potential jet

Write `f=ell^2 W(ell x,ell r) U`, so the actual master satisfies
`U_xx-U_rr=-f`. The potential is separable:

```
W(t,z)=W_static(z)-(35/4)tanh(t)^2.
```

Hence `W_tz=W_ttz=0`. The exact product derivatives are

```
f_xx = ell^4 W_tt U +2 ell^3 W_t U_x +ell^2 W U_xx,
f_xr = ell^3 W_t U_r +ell^3 W_z U_x +ell^2 W U_xr,
f_xxr= ell^4 W_tt U_r +2 ell^3 W_t U_xr
       +ell^3 W_z U_xx +ell^2 W U_xxr.
```

No W_zzz appears. The existing global actual U jets through order three,
W,W_z,W_t, and the new elementary W_tt bound therefore bound every term.
For `z=tanh(t)^2 in [0,1]`,

```
W_tt=-(35/2)(1-z)(1-3z),
-1/3 <=1-4z+3z^2<=1.
```

The lower inequality follows from
`1-4z+3z^2=3(z-2/3)^2-1/3`; the upper follows from
`1-(1-4z+3z^2)=z(4-3z)>=0`. Thus `|W_tt|<=35/2` globally.

## The dynamic shell trace needs both endpoint derivatives

With `y(x)=U(x,0)`, the retained incoming identity is

```
H=Q'+ell exp(-ell x)Q-y'-ell gamma y-I,
I=integral_0^x f(s,x-s) ds.
```

Classical Leibniz differentiation of that actual smooth integral gives

```
I''=f_x(x,0)+f_r(x,0)+integral_0^x f_rr(s,x-s) ds.
```

Both endpoint terms are necessary. The resulting absolute bound is

```
B_H2=B_Q3+ell(B_Q2+2ell B_Q1+ell^2 B_Q0)
     +B_y3+ell |gamma| B_y2+F_x,shell+F_r,shell+3 F_rr.
```

Here `F_r,shell=ell^3 M_z B_y+ell^2 M B_1`; using the global B_1
for U_r at the shell is conservative and valid. The equation
`y''=(ell/mu)H-A y` then yields

```
B_y4=(ell/mu)_max B_H2+ell^4(35/2)B_y
     +2 B_Aprime B_yprime+B_A B_ysecond.
```

This uses actual smooth differentiation of the incoming identity and
boundary ODE. It takes no derivative of the old interval H cells.

## Initial source terms and reflected jumps are retained

Set `V=U_xx`. It solves `V_xx-V_rr=-f_xx`, with

```
V(0,r)=Q''(r)-f(0,r),
V_x(0,r)=ell V_initial''(r)-f_x(0,r),
V(x,0)=y''(x),
V_initial(r)=exp(-ell r)Q(r).
```

The free reflected field for V has every second mixed derivative bounded by

```
B_free=B_Q4+F_rr+ell B_Vinitial3+F_xr+B_y4,
B_Vinitial3=B_Q3+3ell B_Q2+3ell^2 B_Q1+ell^3 B_Q0.
```

The `F_rr` and `F_xr` terms come from differentiating the initial source
values. Leaving them out would incorrectly treat U_xx's initial data as
the derivatives of the free packet alone.

For any smooth physical half-line source g, its odd extension has
`d_r(g_odd)=(g_r)_even+2g(x,0)delta_0`. The reflected principal-wave
identities are therefore

```
(Kg)_rr=1/2 integral[g_r(s,|a+|)-g_r(s,|a-|)]ds
         -1_(x>r)g(x-r,0),
(Kg)_xr=1/2 integral[g_r(s,|a+|)+g_r(s,|a-|)]ds
         +1_(x>r)g(x-r,0),
(Kg)_xx=g+(Kg)_rr.
```

In particular, at the shell `(Kg)_rr=-g` and `(Kg)_xx=0`, fixing the
jump signs and showing why a shell-source trace cannot be discarded.
Applying this to `g=f_xx` gives a common rr/xr correction bound
`3 F_xxr+F_xx,shell`; the xx correction additionally includes F_xx.
The boundary shell value is

```
F_xx,shell=ell^4(35/2)B_y+2ell^3 M_t B_yprime+ell^2 M B_ysecond.
```

The reflected rays stay in the same backward-closed domain. The open
flat initial/shell collar makes y and all corner source jets vanish near
the characteristic join. Thus the physical one-sided traces and reflected
formulas extend continuously to the whole triangle under the retained
smoothness premise. This is a bound on one actual smooth solution, not a
claim that arbitrary independent interval choices define such a solution.

## All five fourth derivatives

Put `B=B_free+3 F_xxr+F_xx,shell`. The differentiated V equation gives

```
|U_xxxx|<=B+F_xx, |U_xxxr|<=B, |U_xxrr|<=B.
```

Differentiating the original principal equation once in x and r, and twice
in r, gives exactly

```
U_xrrr=U_xxxr+f_xr,
U_rrrr=U_xxrr+f_rr.
```

Therefore the remaining two bounds are `B+F_xr` and `B+F_rr`. Each
independently reconstructed rational constant agrees exactly with the
producer report. In the registered scaled coordinates,

```
every total-order-four mixed partial of actual U has magnitude <342395.408,
|y''''|<3.235927.
```

Physical t,z derivatives of total order four have an additional factor
`ell^-4`. These are large absolute envelopes, not a measured error,
relative smallness claim, source sign, or cosmological observable.

The checker rejects eight deliberately omitted or sign-reversed terms.
The analytic derivations above justify the all-point estimates; finite
polynomial witnesses support the implementation without replacing the
underlying smooth actual-field and PDE premises. No common-gauge fourth-jet
certificate, full actual bulk F_x bound, second-response C2 accuracy,
natural Einstein residual, exact nonlinear family, or Big Bang-origin
evidence follows from this bridge alone.
