# Independent trace-recovery aggregation audit

This note checks the mathematical extraction used by a separate frozen
numeric campaign. It introduces no additional numerical result or premise
and does not alter the geometry checker's 21-identity count.

Let x be scaled background time and suppose the actual first shell trace
obeys, on 0<=x<=T,

```
u_xx + a(x)u_x + b(x)u = f(x),    u(0)=u_x(0)=0,
f=ell² rho_b² F'_b v/6,
a=3ell tanh(ell*x),
b=ell²[8sech²(ell*x)−4].
```

The exact trace identity and F'_b=2p Bcal<0 are proved in DERIVATION.md.
Here a_cap and b_cap bound the absolute coefficient values, while every
history source cell has a whole-cell enclosure of the actual f, including
the inherited first-field errors. Define

```
Fint=sum_cells cell_width * sup_cell |f|,
d=1−a_cap T−b_cap T²/2 > 0.
```

For V=sup_[0,T]|u_x|, the zero initial condition gives |u(s)|<=sV.
The once-integrated equation then gives

```
V <= Fint + [a_cap T+b_cap T²/2]V,
V <= Vcap=Fint/d,
|u| <= Ucap=T Vcap.
```

Thus the feedback acceleration −a u_x−b u has absolute value at most
G=a_cap Vcap+b_cap Ucap. Repeated integration gives

```
u_x(x) = ∫₀ˣ f(s)ds                  + error, |error|<=xG,
u(x) = ∫₀ˣ (x−s)f(s)ds              + error, |error|<=x²G/2,
∫₀ˣu(v)dv = ∫₀ˣ (x−s)²f(s)ds/2     + error, |error|<=x³G/6.
```

For source cell [l,h] and integer p=0,1,2, the exact nonnegative kernel
weight is zero when x<=l, and otherwise

```
w_p(x;l,h)=[(x−l)^(p+1)−(x−min(x,h))^(p+1)]/(p+1).
```

Each weight is nondecreasing in x. For a target cell [xlo,xhi], its interval
enclosure is exactly [w_p(xlo),w_p(xhi)]. Multiplication by the whole source
cell interval encloses the corresponding weighted integral: the kernel is
nonnegative and the kernel-weighted mean of f lies within the source
interval. This remains true for source cells truncated at the actual x.
Adding interval contributions can lose correlations but preserves
enclosure. This is continuous integration over whole cells, not point
quadrature.

The proper-time conversion and observable extraction are

```
P=rho_b ell ∫₀ˣu(v)dv,   u'=u_x/(rho_b ell),
A=rho_b cosh(ell*x),     H=tanh(ell*x)/rho_b, H'=1/A²,
E=u−HP, E'=u'−Hu−H'P,  L=P/A², theta1_amplitude=3E'+8L,
C2=10E²/(3A²), S2=80L²/9,
(Var theta)2=theta1_amplitude²/3,
Q2=2[theta1_amplitude²−40L²]/9.
```

All coefficients are second parameter derivatives. The strict signs of
the recovered intervals must come from the complete frozen numerical
report; the equations alone do not prove them. Charging source errors
through the entire history starting at x=0 requires no exact precontact
zero lemma. The independent old actual-master interval can provide an
overlap check; intersecting it is not needed to construct the recovered
signs.
