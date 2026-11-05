# Scalar-level area labels the retained preparation freedom

4 November 2026. Research checkpoint. Earlier published capsules are
preserved unchanged. This is a conditional result within the restricted class of the
4 November preparation-dependence theorem, not a new existence theorem or
physical preparation-selection rule.

## Setting and observable

The specified physical prepared Cauchy section has background data

```
gamma0=dy^2+rho(y)^2 kappa, phi0=phi0(y),
p=phi0_y!=0, a=rho_y/rho>0,
dy/dr=-ell rho, ell=10^-6, 1<=r<=2.
```

Here kappa is the unit three-sphere metric. The actual perturbed section is
not assumed rotationally symmetric. Retain the exact background, first seed,
parameter normalization, canonical momentum and scalar canonical second data,
weighted constraint inverse, conditional exact-family development, linear
uniqueness, physical boundary rows, anchored observables and proper-clock
target T=5/4 of the preserved
[preparation-dependence checkpoint](../HDBLAST_SHELL_PREPARATION_DEPENDENCE_20261104/DERIVATION.md).
In particular each exact family is C2 in its parameter into the retained lower
C4 topology on a common compact neighborhood, after family-dependent
parameter shrinkage. No effective parameter radius or radius uniform over the class is supplied.

For v=phi0(y), let S_e(v) be the connected scalar-level component in the
selected annular collar that continues the background sphere, and set

```
A_e(v)=Vol_3(S_e(v),gamma_e), A0(phi0(y))=|S^3|rho(y)^3,
D_h(y)=[d_e^2 A_e^h(phi0(y))|0-d_e^2 A_e^base(phi0(y))|0]/A0(phi0(y)).
```

The component qualification prevents summing other disconnected scalar
levels elsewhere on the full section. It can be dropped if global uniqueness
is independently known. The scalar label is fixed during parameter
differentiation. D_h uses actual second derivatives, with no extra factor 1/2.

## Finite regularity

Use a compact collar slightly larger than the closed interval of labels.
The uniform nonzero background p bound and C4 closeness make phi_e monotone
radially in this collar for sufficiently small parameter. Every continued
level is a unique graph F_e(y,omega)=(u_e(y,omega),omega), where

```
phi_e(u_e(y,omega),omega)=phi0(y), u_0=y.
```

This equation uses the fully nonspherical metric and scalar. Its pointwise
implicit solution and spatially differentiated equations first give continuity
of u_e in C4. Its first two parameter derivatives satisfy

```
u_dot=-phi_dot/phi_y,
u_ddot=-(phi_ddot+2phi_dot_y u_dot+phi_yy u_dot^2)/phi_y,
```

with right sides evaluated at the graph. Their spatial derivatives through
order two involve at most four spatial derivatives of the input fields and
uniformly nonzero denominators. The mean value formula applied to difference
quotients therefore gives u_e C2 into C2 with both parameter derivatives
continuous in that norm. This allows two spatial orders of loss and avoids
asserting differentiability of composition at unchanged regularity.

The induced metric is

```
sigma_e,AB=gamma_e,AB(F_e)
 +gamma_e,yA(F_e)partial_B u_e+gamma_e,yB(F_e)partial_A u_e
 +gamma_e,yy(F_e)partial_A u_e partial_B u_e.
```

Its first two parameter derivatives and their first spatial derivatives are
continuous: sigma_e is C2 into C1. Its determinant stays positive. Integrating
sqrt(det sigma_e) on the compact sphere gives A_e C2 in the parameter into
C1 in the background label. Thus the second area variation and the one label
derivative used below exist with exactly the retained regularity.

## Difference formula and shared-first cancellation

The two preparations have identical first metric and scalar fields. Write
their second-data differences as delta gamma2 and delta phi2. Their first
level-graph displacements coincide. Differencing the second implicit equation
gives delta u_ddot=-delta phi2/p. Because the background angular graph slope
and gamma0,yA vanish, the second induced-metric difference is exactly

```
delta sigma2,AB=delta gamma2,AB-(delta phi2/p)partial_y gamma0,AB.
```

Every term quadratic in the shared first metric, scalar and graph displacement
has cancelled. The determinant's first-quadratic terms coincide as well. Hence

```
D_h=average_kappa[(1/2)tr_sigma0(delta gamma2,tan)-3a delta phi2/p]. (1)
```

Only now use the prescribed second mean difference. If B_met is its angular
metric coefficient and w its scalar coefficient, (1) reads

```
D_h=(3/2)B_met-3a w/p.                                        (2)
```

Nonmean second differences vanish under this linear average. The absolute
area variations still contain the actual nonspherical first-field quadratics;
none has been set to zero or replaced by an orbit calculation.

In the raw representative the retained homogeneous preparation difference is

```
delta gamma2=h gamma0+Lie_(zeta partial_y)gamma0, delta phi2=0,
zeta=3(h_y-a h)/(2p^2)+9chi/(2rho^4 p^2), chi_y=rho^2 h.
```

Therefore B_met=h+2a zeta and

```
D_h=(3/2)(h+2a zeta).                                        (3)
```

After the retained pure second conditioning gauge, B_met=J=h and w=-p zeta.
Equation (2) gives the same answer. Equivalently,

```
D_h=-(3a/p)[w-p h/(2a)],
w=(p/(2a))(h-2D_h/3).                                        (4)
```

This checks the scalar pullback, sign and normalization of the conditioning.

## Gauge meaning

Scalar-level area is invariant under initial-section diffeomorphisms
transporting both metric and scalar. Formula (1) directly checks a pure
second generator X: the metric trace contributes 3a X^y+div_kappa X^tan,
the scalar term contributes -3a X^y, and the sphere integral of the remaining
divergence vanishes. X need not be radial. The conditioning cannot alter D_h.

The physical prepared Cauchy section is part of the observable. A spacetime
diffeomorphism transporting that section and the fields preserves it. Choosing
a different physical section is a different observable; no independence of
physical spacetime foliation has been proved.

## Reconstruction and injectivity

Retain the exact admissibility: h is smooth, supported in the packet, flat at
both endpoints, and integral_(yL)^(yU)rho^2 h dy=0, with yL=y(2), yU=y(1).
Define chi(y)=integral_(yL)^y rho^2 h ds. It vanishes on both exterior sides.
Let B_area=2D_h/3. Substituting zeta in (3) gives

```
h_y=(a-p^2/(3a))h-3chi/rho^4+(p^2/(3a))B_area.                (5)
```

Put b=ell rho a, g=ell rho p^2/(3a), z=chi/(ell rho^3).
Using rho_r=-ell a rho^2 and chi_r=-ell rho^3 h gives

```
z_r=3b z-h,
h_r=3ell^2 z+(-b+g)h-(2/3)g D_h,
z(1)=h(1)=0.                                                (6)
```

The initial conditions use both flatness of h and its exact zero moment.
The continuous-coefficient linear initial-value problem (6) has a unique
solution for any continuous D_h on a compact interval starting at r=1, by
the fundamental matrix or ordinary integral Gronwall inequality. Thus equal
profiles on [1,T] force equal h, chi and zeta there. Equal profiles on the
whole packet [1,2] force equal full admissible preparations. In particular
D_h=0 on the whole packet implies h=0.

A profile zero only on [1,T] imposes no condition on preparations outside
that causal segment. The map from admissible h to the full profile is
injective, not asserted onto arbitrary C1 profiles: endpoint and zero-moment
compatibility restrict its image. This is a result about the retained mean
freedom with fixed canonical second data, not unrestricted Cauchy data.

## Geometric class and old response estimates

The r derivative in a profile norm means the fixed background scalar-label
parametrization, not an arbitrary perturbed radial coordinate. Equivalently,

```
partial_r D_h=(-ell rho p)partial_v D_h, v=phi0(y(r)).
```

Here rho(v)=(A0(v)/|S^3|)^(1/3) is background level-area data. With this
convention define the geometric, background-calibrated class radius

```
R_A(h)=sup_(1<=r<=T)max(|D_h(r)|,|partial_r D_h(r)|).
```

This bounds a scalar-level area difference and its label derivative on the
specified physical section. Adoption of a particular radius and restriction to
the retained canonical data still need a physical preparation rationale.

The conditioned initial difference in (K,L,N,J,w) order obeys

```
W0=(h,0,0,h,w), W0_x=(0,V_L,0,0,0), W0_r=(h_r,0,0,h_r,w_r),
w=(p/(2a))(h-2D_h/3),
V_L=h_r+b h=3ell^2 z+g(h-2D_h/3),
w_r=(p/(2a))_r(h-2D_h/3)+(p/(2a))(h_r-2D_h,r/3).
```

For H=diag(1,2,1,3,1), set D0=sup|W0|_H, Dr=sup|W0_r|_H,
Dv=sup|W0_x|_H and N(h)=max((D0+T Dv)/3,Dr+Dv). These identities and (6)
reduce the old initial-coordinate norm bound to the geometric profile bound.
The old homogeneous coupled and observable bounds retain exactly their old
analytic hypotheses. The separate exact-rational report
[bounds/references/NORMAL.json](bounds/references/NORMAL.json) verifies
N(h)<=R_A(h)/2000 using the retained continuous background caps. Hence the
closed class R_A(h)<=2*10^33 satisfies N(h)<=10^30, strictly below the
preserved sufficient guard 1.1864*10^30 for both base coefficients to remain
negative. That report supplies the numerical coefficient check; equations
(1)--(6) supply its geometric identification. The normal and optimized
reports are byte-identical, with their finite checks and no new PDE,
background oracle or bump query.

For base+lambda*h, the actual difference profile is lambda*D_h and its radius
is |lambda|R_A(h). The unscaled direction does not bound arbitrary lambda.

This is conditional coefficient-level robustness at the retained target.
Retained C2 dependence yields the corresponding family-wise qualitative
small-parameter inequalities. It yields no uniform nonlinear parameter window
over the class, whose higher norms and remainders need not be bounded. It
does not select the base preparation, prove contraction, provide a finite
amplitude regime, or establish observational cosmology.


## Quantitative certificate for the declared class

The background caps below are verified against the preserved continuous annulus
and shell enclosures, including their initial-slice domain and raw first-jet
normalization. Reused boxes are shared premises, not fresh PDE observations.

```
129<=rho,rho_b<=131, 3/25000<=-p<=1/7000,
111/1000<=a<=14/125, |p_y|<=1/4000, |a_y|<=1/10000,
|Fprime_b|<=1/1000, q=1201/3200000.
```

## Rational ODE and conditioned-data bounds

In this section write `D=D_h`, `R=R_A(h)`, `rhoM=131`, `pM=1/7000`,
`aL=111/1000`, `pyM=1/4000` and `ayM=1/10000`.

Write `bM=ell*131*(14/125)` and
`gM=ell*131*(1/7000)^2/(3*(111/1000))`. Thus

```
bM=917/62500000, gM=131/16317000000000.
```

For the reconstruction matrix in equation (6), its infinity norm is bounded by

```
max(1+3bM,3ell^2+bM+gM)<=mu=20001/20000.
```

For `x=mu/4=20001/80000`, the exact rational Taylor calculation checks

```
exp(x)<=sum_(j=0)^6 x^j/j! + (x^7/7!)/(1-x/8)<4/3.            (Q1)
```

The tail bound follows because ratios after the first omitted term are at
most `x/8<1`. Variation of constants applied to (6), the exponential bound (Q1), and `Delta=1/4` give the
continuous, whole-segment bounds

```
|z|,|h|<=H R, H=Delta*(4/3)*(2/3)*gM=131/73426500000000,
|h_r|<=Hr R,
Hr=3ell^2 H+(bM+gM)H+(2/3)gM.
```

There is no mesh, numerical ODE solution or exponential oracle in this
argument. Write `s=w` and `k=p/(2a)`. The remaining exact identities are

```
k_r=-ell rho[p_y/(2a)-p a_y/(2a^2)],
s_r=k_r(h-2D/3)+k(h_r-2D_r/3), V_L=h_r+b h.
```

Use the rational absolute caps

```
kM=pM/(2aL), krM=ell rhoM[pyM/(2aL)+pM ayM/(2aL^2)],
S=kM(H+2/3), Sr=krM(H+2/3)+kM(Hr+2/3), V=Hr+bM H.
```

They give `|s|<=S R, |s_r|<=Sr R, |V_L|<=V R`. The retained conditioned five-mean
initial difference, in `(K,L,N,J,w)` order, is exactly

```
W0=(h,0,0,h,s), W0_r=(h_r,0,0,h_r,s_r),
W0_x=(0,V_L,0,0,0).
```

For the positive Gram matrix `diag(1,2,1,3,1)`, triangle inequalities give

```
D0<=(2H+S)R, Dr<=(2Hr+Sr)R, Dv<=2V R,
N=max((D0+T Dv)/3,Dr+Dv), T=5/4.
```

The script reconstructs the exact resulting constant

```
N<=c_N R<=R/2000,
c_N=798916710779559466175931203077/
    1861847711577000000000000000000000
   <0.000429099<1/2000.                                      (Q2)
```

The decimal in (Q2) is only an outward display; the reported rational and
the `1/2000` comparison use exact `Fraction` arithmetic. All five channels
and the mixed metric velocity are retained.

## Transfer to the two complete second observables

The retained homogeneous mean estimate is `X(W)<=N/(1-q)`. The initial
collar and finite propagation give the exact zero wedge `x+r<=1`, hence the
boundary scalar value obeys

```
|m(x)|<=(x-1)_+ N/(1-q).
```

The exact retained kernels satisfy `0<=k(T,x)<=T-x` and `0<=k_T(T,x)<=1`.
Integrating, with the exact response forcing factor `F'_b/2`, gives

```
|Delta v2|<=Cv N, Cv=(131ell)^2*(1/1000)*Delta^3/[12(1-q)]
                    =17161/767711760000000000,
|Delta Theta2|<=Ct N, Ct=131ell*(1/1000)*Delta^2/[4(1-q)]
                       =131/63975980000.                      (Q3)
```

For example `integral_1^T (T-x)(x-1)dx=Delta^3/6`; the scalar forcing's
`1/2` produces the denominator 12. These are complete coupled differences:
the identical actual first forcing, proper-clock and measure corrections and
first-quadratic Gauss terms cancel when subtracting families. Every such
term remains in the preserved absolute base intervals.

From the exact old negative upper endpoints and (Q3), the exact sufficient
threshold is their minimum ratio, reconstructed directly from the authenticated published base intervals. Its outward display is

```
Ncrit>1.18646887e30.
```

For the frozen **closed** class `R<=2e33`, (Q2) yields `N<=1e30<Ncrit`.
The actual rational bound is stronger, below `8.5820e29`, but the final
gates use the originally fixed `1e30` guard. This gives absolute difference
caps below `2.23535e16` for volume and `2.04765e21` for expansion.
Adding both signs of the exact difference bounds to the exact preserved
base intervals gives the predeclared outward enclosures

| Whole class at the fixed proper target | Outward enclosure |
| --- | --- |
| `v2` | `[-7.6900e16,-4.1600e15]` |
| `Theta2` | `[-1.1070e22,-2.3980e21]` |

Both upper endpoints are strictly negative. The complete exact interval
endpoints are in [bounds/references/NORMAL.json](bounds/references/NORMAL.json); the finite gates include both
outward containments and both strict signs.



## Compatibility with the earlier opposite-sign example

The comparison concerns the complete preparation increment. For the preserved
P1 direction and `lambda=-10^21`, independent rational reconstruction of the
old initial target boxes gives a safe outward area-profile enclosure

```
D_h(r=T) in [3.7621e12,4.5765e12],
D_h,r(r=T) in [-8.3575e13,-6.8702e13].
```

Therefore `|lambda D_h(T)|>3.7621e33>2e33`. The derivative gives an even larger
exclusion. Here `r=T` is an initial bulk scalar label, not the future shell
proper-time event. The previous sign reversal is outside the stated class;
it neither refutes this restricted bound nor becomes an admissible class member
by using the unscaled direction's norm.

The primary exact-rational calculation and separately written reconstruction
share the retained source enclosures. They do not independently solve the PDE,
prove the analytic assumptions, justify the geometric radius physically, or
establish a new general gauge-invariance principle. The initial section and
continued scalar-level component remain declared parts of the construction.
