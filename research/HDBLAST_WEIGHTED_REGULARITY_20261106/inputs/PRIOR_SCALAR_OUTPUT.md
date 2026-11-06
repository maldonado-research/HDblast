# A quantitative initial metric output bound for a specified weighted member

6 November 2026. Conditional initial-data result in the retained
Einstein–scalar model. All parameter coefficients are amplitude derivatives,
not ordinary Taylor coefficients.

## Statement

Retain the actual background, first packet, scalar/density-affine
prescription and `k=12` construction in
[REGISTERED_EXACT_FAMILY_BINDING.md](inputs/REGISTERED_EXACT_FAMILY_BINDING.md)
and [EXACT_INITIAL_THEOREM.md](inputs/EXACT_INITIAL_THEOREM.md). Specify the
following member of the admissible weighted constructions:

```
ell=10^-6, A={ell/2<z<5ell/2}xS3,
r(z)=(z−ell/2)(5ell/2−z)/ell²,
kappa=ell/100, d=kappa r, beta=kappa, a_w=2,
q=d², w=exp(−1/r).
```

Here d is a defining function, a_w is the integer weight parameter, and
the spatial dimension is four. The symbol r in this note is the displayed
dimensionless defining-function shape. The warp coefficient below is a.

Conditional on the retained exact-initial-data theorem and its background
premises, the selected metric second derivative `h_star=gamma''(0)` obeys

```
||h_star||_L2,prob <2.5e40.                                      (1)
```

The norm uses the fixed background tensor metric and the measure
`rho³|dy| dmu_kappa,prob=rho⁴dz dmu_kappa,prob`, with probability measure
on the unit sphere. The geometric spatial L2 norm is `pi sqrt2` times
this norm and consequently is less than `1.12e41`.

The scalar and scalar momentum-density second derivatives are zero in
this construction. The normal scalar momentum second derivative has the
known quadratic canonical conversion; it is not claimed zero.

This explicitly specifies a weight member. The numerical bound is not
an identification of an earlier construction whose defining function and
weight parameters were unspecified. For every such fixed admissible
radial weight, the same proof gives the separate inequality

```
||h_star||_L2,prob <5.58e39 C_H,
C_H=sqrt(max_A mu_H / min_K mu_H), K=[ell,2ell]xS3,
mu_H=w²q⁴.                                                     (2)
```

No H14 bound, evolved second-field estimate, nonlinear parameter radius,
or common physical observation cylinder follows from(1).

## Fixed operators and the minimum-energy property

On the initial slice,

```
gamma0=dy²+rho² kappa, rho_y=a rho,
a_y=−rho^−2−p²/3, p=phi0_y, p_y=V'−4ap,
R(gamma0)=p²+2V.
```

The canonical normalization in
[CANONICAL_CONSTRAINTS.md](inputs/CANONICAL_CONSTRAINTS.md) gives

```
P_g h=(div div h−Delta tr h−Ric:h+(dphi0 tensor dphi0):h)/2,
P_g*N=(Hess N−gamma0 Delta N−N Ric+N dphi0 tensor dphi0)/2.
```

All pairings use the fixed background volume. At the time-symmetric
background, the Hamiltonian metric correction and momentum correction
decouple. Put `f=−D²Chat_H[v_seed,v_seed]`. The selected metric solves

```
h_star=mu_H P_g*N, P_g h_star=f,
mu_H=w²q⁴=d^(4a_w)exp(−2beta/d).                              (3)
```

For tensors define `||h||_E²=integral_A mu_H^−1|h|² dmu0`.
Every finite-energy correction h with the same source and physical
zero-extension regularity satisfies

```
<h_star,h−h_star>_E=0,
||h_star||_E²+||h−h_star||_E²=||h||_E².                       (4)
```

To justify the pairing, choose compact smooth potentials N_j converging
to N in the dotted weighted domain. That convergence controls
`mu_H^(1/2)P_g*N_j` in L2. The distributional equation
`P_g(h−h_star)=0` gives `integral(P_g*N_j):(h−h_star)=0`.
Weighted Cauchy–Schwarz passes to the limit and proves(4). This uses no
unweighted boundary trace of N. Equivalently, the normalized equation
`w^−2 P_g(mu_H P_g*N)=w^−2 f` is tested with the weighted pairing w²;
its right side is exactly `integral Nf`.

Thus a compact particular correction gives a bound for the actual
selected output without estimating the inverse operator norm. It does
not identify that output with the particular correction.

## Explicit scalar-frozen compact corrections

For any smooth scalar N, divergence of the complete adjoint and the
background equations yield

```
div(P_g*N)=(p²/2)(N_y−aN)dy.                                 (5)
```

Indeed divergence of the Hessian and Ricci-gradient terms cancels,
leaving
`−N dR/4+[(grad N.grad phi0)+N Delta phi0]dphi0/2
 +N Hess(phi0)(grad phi0,.)/2`; substitute the background equations.
The scalar-gradient term in P_g* is essential.

Let `N=n(y)E_l`, `Delta_kappa E_l=−lambda E_l`. Its radial adjoint
component is

```
Q_yy=−(3a/2)(n_y−an)E_l+(lambda−3)nE_l/(2rho²), Q=P_g*N.
```

For lambda different from three, this and(5) imply

```
N=c Q_yy+b(div Q)_y,
c=2rho²/(lambda−3), b=6a rho²/[(lambda−3)p²].                  (6)
```

Taking the formal adjoint of(6) gives a compact correction for a
compact scalar source `F_l(y)E_l`:

```
h_part=cF_l E_l dy²−sym grad(bF_l E_l dy),
h_yy=[cF_l−(bF_l)_y]E_l,
h_ya=−(bF_l/2)D_aE_l,
h_ab=−abF_l E_l gamma0_ab,
P_g h_part=F_l E_l.                                          (7)
```

Here `sym grad alpha=(nabla_i alpha_j+nabla_j alpha_i)/2`.
No scalar second correction is present. The retained p is bounded away
from zero, so all coefficients are smooth. A flat compact source gives
a flat compact metric with the same support. For unit sphere norm of
E_l, the exact angular-integrated tensor norm squared is

```
[cF_l−(bF_l)_y]²+3a²b²F_l²+lambda b²F_l²/(2rho²).             (8)
```

For lambda8 and24 a simpler algebraic correction is available. Set
`T_l=Hess_kappa E_l+(lambda/3)E_l kappa`; then
`div T_l=−2(lambda−3)dE_l/3` and
`||T_l||²=2lambda(lambda−3)/3`. Choose

```
h_yy=h_ya=0,
h_ab=3rho⁴F_l T_l,ab/[lambda(lambda−3)].                       (9)
```

It has zero trace and zero Ricci and scalar-gradient contractions, and
its double divergence is `2F_l E_l`. Hence it also solves the same
scalar-frozen source equation. Its norm squared is

```
6rho⁴F_l²/[lambda(lambda−3)].                                (10)
```

This uses the standard scalar-derived trace-free Hessian inversion
already present in the localized constraint construction. It is not
an angular transverse-traceless tensor. The monopole uses(7); the
excluded eigenvalue three is the dipole and is absent by reflection
parity. Rotational invariance of the radial weighted normal operator
and uniqueness preserve these three source channels.

## Actual affine source and evaluated coefficients

The complete first datum and bump are fixed by
[REGISTERED_EXACT_FAMILY_BINDING.md](inputs/REGISTERED_EXACT_FAMILY_BINDING.md).
Writing `w_v=u_t` for the first velocity profile, the actual canonical
Hamiltonian Hessian is

```
D²Chat_H[v_seed,v_seed]=F Y²+J|DY|²,
F=−12u u_yy−12u_y²−144a u u_y
 −(192a²+96a_y+168/rho²)u²
 −(v_y+4up)²−V''v²−(6w_v²+v_t²)/rho²,
J=(18u²−v²)/rho².                                           (11)
```

These expressions are derived directly in the companion affine-source
calculation. The Hamiltonian source also equals half the physical
constraint Hessian in
[EXPLICIT_COMPACT_SECOND_INITIAL_JET.md](inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md),
because its first derivative has no K or Pi term at time symmetry.
No second correction prescribed in that earlier construction is used.

The actual background box in
[ACTUAL_POTENTIAL_JET_BOUNDS.json](inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json)
gives `129<=rho<=131`, `111/1000<=a<=14/125`,
`3/25000<=|p|<=7/50000`, and `|V'|<28/100000`. Substitution in(8) gives
for the monopole

```
|b|<2.67e11, |c|+|b_y|+sqrt3|ab|<1.6e12,
||h_part,0||<=1.6e12||F_0||+2.67e11||partial_y F_0||.          (12)
```

For example
`b=−2a rho²/p²` and
`b_y=−(2rho²/p²)[−rho^−2−p²/3+10a²−2aV'/p]`.
For unit-normalized degree2 and4 harmonics, (10) gives respectively
the unweighted gains6650 and1873. These are rational-safe ceilings.

The source calculation supplies the global bounds

```
|F|<1.22e22, |J|<4.1e6,
|F_y|<7.3e27, |J_y|<3.5e12.                                 (13)
```

They use the exact bump `Q=exp(4−1/[x(1−x)])`,
`x=(z−ell)/ell`, with flat extension, and rational background boxes.
The extra derivative is controlled by `|Q_zzz|<35000/ell³`,
`|u_yyy|<1.3e12`, `|v_yy|<3.3e16`, `|(v_t)_y|<1.1e11`.
The bump bound follows by differentiating in `b=1/[x(1−x)]>=4`
and maximizing each polynomial term times exp(4−b); no sampling is used.

Use harmonics

```
H0=1, H2=(n1²+n2²)/2−1/4, H4=Y²−(n1²+n2²)/2+1/12,
||H_l||²=(1,1/48,19/720), lambda_l=(0,8,24).
```

The source coefficients for `f=−D²Chat_H` are
`f0=−F/6−4J/3`, `f2=−F−4J`, `f4=−F+4J`.
Set the exact outward bounds

```
B0=1.22e22/6+4(4.1e6)/3,
B0y=7.3e27/6+4(3.5e12)/3,
B24=1.22e22+4(4.1e6).
```

The source and all compact comparison tensors are supported in
`K=[ell,2ell]xS3`. Since `dy=−rho dz`, their radial volume satisfies

```
integral_K rho³|dy|<=131⁴ell=294.499921,
Vhalf=sqrt(131⁴ell)=17.161.
```

Orthogonality of the three tensor harmonic sectors gives

```
||h_part||² <=Bpart²
 =[Vhalf(1.6e12 B0+2.67e11 B0y)]²
  +(Vhalf*6650*B24)²/48
  +(19/720)(Vhalf*1873*B24)²
 <(5.58e39)².                                               (14)
```

All terms, including the nonconstant modes, are retained in the exact
rational squared comparison. Equations(3),(4),(14) prove(2).

## Explicit weights, cover and final comparison

For the specified member, `max_A r=1` and `min_K r=3/4`.
The multiplier `mu_H=kappa^8 r^8 exp(−2/r)` is increasing in r, so

```
C_H=(4/3)^4 exp(1/3)<4.42.
```

The constant kappa cancels. An explicit Taylor sum with geometric tail
majorant bounds exp(1/3) in rational arithmetic. Combining this with(14)
gives `(4.42)²Bpart²<(2.5e40)²`, proving(1). Sphere conversion uses
`pi sqrt2<4.45`, which follows from `pi<22/7`, `sqrt2<99/70`.

The scaled cover is admissible. Let
`delta=min(z−ell/2,5ell/2−z)`. Then `0<delta<=ell`,
`r=delta(2ell−delta)/ell²<=2delta/ell`, and

```
q=kappa²r²<=(ell/2500)delta,
dist_gamma0(point,partial A)>=129delta.
```

Thus even radius4q balls stay inside A, in either the physical metric
or the product reference metric. Also `|grad_gamma0 r|<=2/(129ell)`.
On the physical-metric radius4q ball about x, relative r variation is at most
`epsilon=8kappa²/(129ell)<1e−11`, and

```
|log[w(x')/w(x)]|<=epsilon/(1−epsilon)<1e−11.
```

This verifies q and w comparability on the scaled balls. The smooth
quadratic defining function and exponential weight have the required
finite-order scaled multiplier bounds. If an equivalent coordinate norm
requires a smaller scale, simultaneous rescaling `d→t d,beta→t beta`
multiplies the scalar and momentum correction operators by constants
`t^(4a_w)` and `t^(4a_w−4)`. Fixed finite-order weighted norms are
equivalent, and inverse rescaling of the potential variables preserves
their physical ranges and the exact implicit-function curve. This
scaling does not identify different shape functions or beta/d ratios.

## Finite radial normal operator and limits

For completeness, write `D=partial_y`, `D^dagger=−D−3a`,
`R_t=3/rho²+p²/3−3a²`, and

```
A_l=−(3a/2)D+[3a²+(lambda−3)/rho²]/2,
B_l=−D²/2−aD+lambda/(3rho²)−R_t/2, J_l=D−a.
```

The exact scalar radial normal equation on the dotted weighted domain is

```
[A_l^dagger mu_H A_l+3B_l^dagger mu_H B_l
 +J_l^dagger(mu_H lambda/(2rho²))J_l
 +mu_H lambda(lambda−3)/(6rho⁴)]n=F_l.                        (15)
```

Its energy is the integral against rho³dy of
`mu_H[(A_l n)²+3(B_l n)²+lambda(J_l n)²/(2rho²)
 +lambda(lambda−3)n²/(6rho⁴)]`.
This records the complete finite scalar block; it does not impose
ordinary Dirichlet boundary conditions in place of the weighted domain.

The bound is deliberately coarse. A small nonzero scalar gradient can
make the monopole expensive: the lapse N=a rho has adjoint components
`Q_yy=a rho p²/2` and
`Q_ab=[rho p(V'−2ap)/3]gamma0_ab`, becoming an exact radial kernel in
the scalar-constant limit. Narrow annular width alone therefore does
not supply a small weighted output.

The identities, source joins, norm coefficients and final squared
arithmetic are finite verifications conditional on the stated analytic
theorems. Passing these checks does not establish a high-order inverse
estimate or propagate the bound to the physical second scalar, Weyl
stress, observer map or area. Those tasks remain separate.
