# Actual common short observer domain and geometric pullbacks

16 November 2026. New proof at immutable base
`5be99501c7f58e4400212f374abb6334c22cbca5`. All previous roots and
capsules remain read-only. The proper duration below is a qualitative
positive number. No historical target duration, evaluated observer
budget or new parameter restriction is asserted.

## 1. Actual inputs and conclusion

Use the current canonical family on its original compact parameter
interval `I`, the physical original-time slab `[-T,0]`,
`T=13ell/8`, `ell=10^-6`, and the actual fixed shell inclusion

```
iota(t,n)=(t,0,n),           h_e=iota^*g_e,
varphi_e=iota^*phi_e.
```

The retained physical join proves that this inclusion is an embedding
of the shell in the specified one-sided bulk chart, with all original
shell equations. The new parameter theorem gives the actual fixed-base
physical curve `C2` into `C X10`. The regularity and anchor companion
`COMMON_PROPER_CLOCK_AND_COMPOSITION.md` supplies
`h_e,varphi_e in C2(e; joint C7)` and their common Lorentzian
patch. In particular `K=partial_t` is future timelike; the function `t`
is temporal. The initial shell pair is the fixed background pair, so

```
h_e(0)=rho_b^2(-dt^2+kappa),
x_e(0,n)=(0,n),             u_e(0,n)=rho_b^(-1)partial_t.
```

The anchor and these initial velocities have zero first and second
parameter derivatives. The induced initial cut metric is exactly
`rho_b^2 kappa`, with `rho_b>0` fixed.

The actual companion `INTRINSIC_NORMAL_OBSERVER_CONSTRUCTION.md`
constructs, on a common interval
`0<=r<=r0`, all labels `n in S3`, and a compact parameter set `Kpar`
containing `I` in the interior of the analytical reduced extension,

```
x_e(r,n)=(t_e(r,n),omega_e(r,n)),
x_r=-u,             nabla^h_u u=0,
h(u,u)=-1,          u future,
t_e(r,n) in [-T/2,0],
(x,u) in C2(e; joint C4(r,n)),
(x,u,P,Q) in C2(e; joint C3(r,n)),
P=d_n x,            Q=d_n u.                              (1)
```

All derivatives here are actual compatible derivatives. The coordinate
and parameter dependence are jointly continuous on the compact common
domain. The initial angular map is the identity. We use this constructed
ODE result, not an assumed long observer tube or a formal variation
system. The physical equations are used only for parameters in `I`;
the analytical extension permits ordinary endpoint derivatives.

Then there is one common `0<r_*<r0`, independent of `e in I` and of
`n`, such that:

* `x_e:[0,r_*] x S3 -> [-T/2,0] x S3` is an embedding, including
  injectivity between different proper-clock cuts;
* every cut has metric `k_e=x_e(r,.)^*h_e` satisfying
  `(rho_b^2/2)kappa <= k_e <=(3rho_b^2/2)kappa`;
* the scalar, intrinsic curvature, full cut metric, area density and
  specified actual-frame curvature contractions have the actual `C2`
  parameter regularity and complete derivatives established below.

The physical shell observer map is `gamma_e=iota o x_e`. It is an
embedding into each chosen one-sided physical bulk chart because the
current `iota` is an embedding. These observers are geodesic within the
shell; they are generally accelerated in the ambient metric. The
proper past duration is `r=-tau`: `h(x_r,x_r)=-1`. The historical scaled
variable `q=r/rho_b` would instead give norm `-rho_b^2`. It is not
silently substituted for `r` in this proof.

## 2. One uniform positive cut gap from the actual initial gap

On the ODE domain (1), define the full pullback cut metric

```
k_e(r)(v,w)=h_e(x_e(r,n))(P_e v,P_e w),  v,w in T_n S3.
```

No time component of `P` is suppressed. The regularities in (1) and the
fixed-base metric regularity make `k` jointly continuous in
`(e,r,n)` as a symmetric bilinear form on the fixed label sphere.
At `r=0`, `P=(0,Id)` and `k_e(0)=rho_b^2 kappa` for every parameter.

Use the operator norm relative to `kappa` and set

```
m(a)=sup_(e in Kpar,0<=r<=a,n in S3)
                ||k_e(r,n)-rho_b^2 kappa||_kappa.
```

Uniform continuity on the compact product implies `m(a)->0` as
`a->0+`. Choose once `0<r_*<r0`, for example also `r_*<=r0/2`, so that
`m(r_*)<=rho_b^2/2`. This proves the common two-sided cut estimate in
Section 1. It shortens only the new observer duration, preserving the
entire original parameter interval and the original PDE time slab.
The fraction one-half is a convenient qualitative gap; no value of
`r_*` or evaluated Jacobi constant is claimed.

In particular every `P_e(r,n)` is injective: if `P v=0`, positive
definiteness of `k` gives `v=0`. Every full cut is a spacelike immersion.
The construction uses complete label spheres at every earlier `r`,
not labelwise stopped fragments.

## 3. Actual Gauss orthogonality and global angular inverse

Fix a parameter and write `U=x_r=-u`. For a label direction `v`, put
`P_v=d_n x[v]`. Geodesic normalization and the torsion-free pullback
connection give

```
nabla_U U=0,                 h(U,U)=-1,
nabla_U P_v=nabla_(P_v) U,
partial_r h(U,P_v)
 =h(nabla_U U,P_v)+h(U,nabla_U P_v)
 =(1/2)partial_v h(U,U)=0.
```

The last two equalities are identities along the two-parameter family;
they need no prior global injectivity. All derivatives exist by (1).
Initial normality gives the exact identity

```
h(U,P_v)=0                                    (2)
```

on the entire common cylinder. Thus Gauss orthogonality has been proved
for these actual normal geodesics, not assumed for a generic moving map.

For fixed `r`, let `A_(e,r)(n)=omega_e(r,n)`. If
`D_n omega[v]=0` for a nonzero label vector `v`, then
`P_v=(D_n t[v])K`. A nonzero multiple of timelike `K` has negative
square, and the zero vector has zero square. Both contradict positive
`k(v,v)`. Therefore `D_n omega` is invertible everywhere.

The smooth map `A_(e,r):S3->S3` is a local diffeomorphism. Its image is
open by local invertibility and closed by compactness, hence is the
whole connected sphere. A proper local diffeomorphism is a finite
covering. The actual maps `A_(e,s)`, `0<=s<=r`, give a homotopy from
`A_(e,0)=Id`. Thus `deg A_(e,r)=1`. On a connected oriented sphere
the local orientation sign of a local diffeomorphism is constant; its
degree is that sign times the number of sheets. Degree one forces one
sheet. Consequently every `A_(e,r)` is a global diffeomorphism.

Write `B_(e,r)=A_(e,r)^(-1)`. Finite-chart implicit differentiation and
uniqueness make `B` jointly regular; on the compact tube the inverse
angular differentials have one finite uniform bound. Parameter
differentiation uses the actual lower regularity of inverse composition:
from the `C2(e;C4)` angular maps one may safely take
`B in C2(e; joint C2(r,m))`. No unchanged top-order differentiability
of a same-order diffeomorphism-group operation is claimed. Only its
joint `C1` regularity is needed for the injectivity argument.

## 4. Injectivity between different proper-clock cuts

Use the global angular inverse to write each cut as a graph over the
fixed angular coordinate `m`:

```
Theta_e(r,m)=t_e(r,B_(e,r)(m)).
```

At fixed `m`, differentiating the identity for `omega` gives

```
V=partial_r[x_e(r,B_(e,r)(m))]
 =U-P (D_n omega)^(-1) omega_r
 =Theta_r K.
```

Contract with `U`. By (2),

```
Theta_r h(K,U)=h(V,U)=h(U,U)=-1.                    (3)
```

`K` is future timelike and `U` is past unit timelike, so their
Lorentzian inner product is strictly positive. Hence `Theta_r<0`
everywhere. Continuity and compactness even give a finite positive
uniform lower bound for `-Theta_r`; no numerical value is needed.

If `x_e(r1,n1)=x_e(r2,n2)=(t,m)`, then
`n_i=B_(e,ri)(m)` and `Theta_e(r1,m)=Theta_e(r2,m)`. Strict monotonicity
in (3), on the complete connecting interval, gives `r1=r2`; angular
injectivity gives `n1=n2`. This proves injectivity of the whole proper
cylinder, including points on different cuts. It uses neither a
coordinate-time inverse outside its common image nor an independent
small-angle/Schur guard.

The full differential is injective as well. If `a U+P_v=0`, contraction
with `U` and (2) gives `-a=0`, followed by `v=0` from cut positivity.
The domain is compact and the target Hausdorff, so the continuous
injection is a homeomorphism onto its image. The injective differential
and the established regularity make it an embedding, up to both
proper-clock boundary cuts. At `r>0`, (3) and `Theta(0,m)=0` also give
`t<0`; the ODE companion already supplied the common past-edge margin.

The fixed physical inclusion `iota` now transfers this result to
`gamma_e=iota o x_e`. Pullback metric equality alone would not prove
ambient injectivity for an arbitrary immersed moving shell; the actual
current embedded inclusion is a separate retained input used here.

## 5. Actual moving scalar and curvature compositions

All parameter derivatives below hold the actual proper duration `r`
and transported initial label `n` fixed. Put

```
Y=partial_e x,                Z=partial_e^2 x.
```

At the initial time edge, use the fixed bounded finite-regularity
extensions from the regularity/ODE companions when applying ordinary
open-domain composition rules. Equivalently extend the metric/scalar
component arrays across `t=0` by that parameter-independent linear
extension, retain the Lorentzian patch on a sufficiently small common
extension strip, and form the inverse, connection and curvature there.
Their one-sided jets through the stated orders agree at `t=0` with
the actual past fields. The reconstructed bulk components can be
extended in the same way for one-sided bulk curvature/normal arrays.
The actual observer retention keeps every evaluated point in `t<=0`,
and at `r>0` strictly in `t<0`; hence all final values and jets are
the original past ones, independent of the chosen extension. No future
field equation or physical continuation is asserted. This device is
needed for the Banach composition justification at the boundary, not
for an enlargement of the actual observer domain.

In a fixed shell chart, for any scalar or fixed-chart tensor component
field `f_e`, define `f0=f_e(x_e)` and its complete parameter derivatives

```
f1=(partial_e f_e)(x_e)+D f_e(x_e)[Y],
f2=(partial_e^2 f_e)(x_e)
       +2D(partial_e f_e)(x_e)[Y]
       +D^2 f_e(x_e)[Y,Y]+D f_e(x_e)[Z].             (4)
```

The explicit field derivatives, mixed derivatives, observer first
variation, Hessian and observer second variation all occur. At a
general parameter neither the last term nor the moving-evaluation
terms vanish. Smooth fixed chart transitions and their Hessians give
the same actual parameter curve; these local component formulas do not
declare coordinate components themselves invariant observables.

The current shell scalar is `C2(e;C7)` and `x` is `C2(e;C4)`. The
finite composition rule, with two state derivatives charged in (4),
therefore gives the actual pulled scalar

```
Phi_e(r,n)=varphi_e(x_e(r,n)) in C2(e; joint C4(r,n)). (5)
```

The intrinsic shell Riemann tensor and scalar curvature have fixed-base
regularity `C2(e;C5)`. Their moving values from (4) are safely
`C2(e;C3)` since `3+2=5`. The same statement holds for one-sided bulk
curvature restricted to this fixed shell and then evaluated at `x`.
All are on the actual common observer domain, now proved above.

For a genuine tensor pullback, include the differential of the map.
Let `T_0=U=-u`, `T_A=P_A`, and
`T1_a=partial_a Y`, `T2_a=partial_a Z`. For a covariant rank-four
tensor `R`, its full observer-coordinate pullback is

```
Robs_abcd=R0(T_a,T_b,T_c,T_d).
```

Its first derivative is `R1(T^4)` plus the four terms with one `T`
replaced by `T1`. Its second derivative is

```
R2(T^4)
 +2 sum_(one slot) R1(T,...,T1,...,T)
 +  sum_(one slot) R0(T,...,T2,...,T)
 +2 sum_(unordered distinct slot pairs)
                         R0(T,...,T1,...,T1,...,T).       (6)
```

Here `R1,R2` are the full moving field derivatives (4). This formula
keeps every velocity/label differential variation. As `Dx in C2 C3`,
the full tensor pullback is `C2(e;C3)`. It is a coordinate-independent
tensor on the fixed observer cylinder. The scalar curvature from (4)
is already a scalar and needs no such differential slots.

These statements do not claim unchanged top `C5` regularity for a
twice parameter-differentiated curvature composition. No clock inverse
is used in (4)--(6): `r` is already the actual proper parameter of the
constructed observer. Replacing it by a coordinate-time target would
introduce a separate common-image inverse and its derivatives.

## 6. Actual observer frames and curvature contractions

For an explicit physical frame within the shell, choose once a smooth
global orthonormal round frame `e_A(n)`, `A=1,2,3`, on `S3`; for
example the three unit-quaternion multiplication fields. This is a
smooth global frame, not an angular truncation or one global coordinate
chart. Define

```
E_A(0,n)=rho_b^(-1)(0,e_A(n)),
nabla^h_U E_A=0,
in coordinates: (E_A)_r=Gamma[h](x)[u,E_A].           (7)
```

The coefficients are those of the actual induced metric `h`. This is a
linear ODE along the already constructed common flow. Their finite
joint regularity and compactness give a unique solution on all
`[0,r_*]` with actual `C2(e;C4)` dependence; the initial parameter
derivatives are zero. Parameter differentiation of (7) gives the full
linear and quadratic coefficient/point/velocity/frame terms, rather
than holding a frame fixed at a moved point. Metric compatibility gives

```
h(E_A,E_B)=delta_AB,          h(u,E_A)=0.
```

For completeness, in a fixed COORDINATE bundle trivialization put
`L_e=Gamma[h_e](x_e)[u_e,.]`. In a global product/round-frame
representation this symbol must mean the full connection coefficient
in that frame, including its fixed frame-derivative terms. Equivalently,
with the fixed product connection `D0` and connection difference
`Aconn=nabla^h-D0`, the invariant equation is
`D0_r E=Aconn(u,E)`; writing it in components also differentiates the
fixed frame. No noncoordinate global frame is treated as one coordinate
chart. Smooth fixed transitions and uniqueness glue the local linear
equations. Their fixed additional terms have the same finite reserve.
The two-base-derivative reserve
`6=4+2` gives `L in C2(e; joint C4)`. The integral equation
`E(r)=E(0)+integral_0^r L(s)E(s) ds` has its ordinary uniformly
convergent linear Picard series on this finite interval. Differentiate
label variables through four in that equation. Each differentiated
row is triangular, with the same leading `L` and finite lower rows;
Gronwall bounds every row uniformly on the compact parameter/label
domain. Difference quotients in epsilon converge to the unique rows

```
(E1)_r=L E1+L1 E,
(E2)_r=L E2+2L1 E1+L2 E,
E1(0)=E2(0)=0.
```

Here, with `G=Gamma[h]` and its full moving derivatives `G1,G2`
defined by (4),
`L1=G1[u,.]+G[u_epsilon,.]` and
`L2=G2[u,.]+2G1[u_epsilon,.]+G[u_epsilon,epsilon,.]`.
The same integral estimates applied to remainders and continuous
coefficient differences prove actual `C2` dependence in the label
`C4` norm. Time derivatives in (7) recover all joint derivatives of
total order at most four and their parameter continuity. This proves
the claimed frame regularity without a further nonlinear PDE or an
extra state domain; the linear transport exists on the entire already
constructed interval.

Thus the intrinsic shell tidal components

```
E^h_AB=Riem[h](u,E_A,u,E_B)                            (8)
```

are actual `C2(e;C3)` scalar components in this chosen transported
orthonormal frame. Their first and second derivatives have (4) for
the moving curvature plus every first/second variation of the four
vectors, with the one-slot and two-slot product terms in (6).

One can also specify a one-sided ambient measurement frame without
equating the observers with ambient free fall. Let `n_out` be the actual
outward spacelike unit shell normal in the chosen copy. The completed
physical law, with the retained convention, is

```
K(V,W)=g(nabla^g_V n_out,W)=(sigma/6)h(V,W).
```

For the intrinsic unit geodesic `u`, the Gauss and Weingarten formulas
therefore give, in FUTURE proper time,

```
a^g=nabla^g_u u=(sigma/6)n_out,
nabla^g_u E_A=-K(u,E_A)n_out=0,
nabla^g_u n_out=(sigma/6)u.                            (9)
```

Consequently the four screen vectors `(E_1,E_2,E_3,n_out)` are
orthonormal and perpendicular to the actual future `u`, and satisfy
the Fermi--Walker equation

```
nabla^g_u V = u g(a^g,V)-a^g g(u,V).
```

For evolution in `r` its covariant left side is negated because
`U=-u`. These are per-copy one-sided identities. In particular the
normal need not be constant along the observer. The normal is a smooth
algebraic function of the current metric on the positive radial-normal
patch; it and both parameter variations composed with `x` are safely
`C2(e;C4)`.

The Israel law in this argument is used only for `e in I`. Frame
transport and its analytical parameter extension use the metric family
on the open reduced neighborhood, but impose no physical Israel law
outside `I`. When differentiating (9) or its Fermi--Walker consequence,
first work in the interior of `I` and pass to its endpoints by the
proved parameter continuity. No extension of a physical shell law
beyond its retained interval is assumed.

Thus, if desired, the explicitly defined one-sided ambient tidal
matrix `Riem[g](u,S_I,u,S_J)` in this frame is also `C2(e;C3)`, with
all curvature, velocity, frame and normal variations included. It is a
specified geometrical observable for these confined observers; it is
not a detector/noise calibration or a claim of smooth geodesic
transmission through the thin shell. The two copies keep their own
outward normal and one-sided curvature.

## 7. Full cut metric Hessian and area density

Let `H=h_e(x_e)` in fixed local components, and let `H1,H2` be its
complete moving derivatives from (4). Put

```
P1=d_n Y,                 P2=d_n Z.
```

These are supplied as actual mixed parameter/label derivatives by the
ODE and augmented-label theorem. The full pulled metric and its exact
parameter derivatives are

```
k=P^T H P,
k1=P^T H1 P+P1^T H P+P^T H P1,
k2=P^T H2 P+2P1^T H1 P+2P^T H1 P1
                 +P2^T H P+P^T H P2+2P1^T H P1.       (10)
```

In particular `P2`, the mixed field/observer derivative in `H2`, and
the two-label-first-variation term cannot be omitted. A spatial-only
metric restriction `h_AB(t,omega)` would not be the metric of these
moving cuts. Formula (10) and the stated regularities give

```
e -> k_e in C2(e; joint C3([0,r_*]xS3)).              (11)
```

Define the global area density relative to the FIXED round label
measure by

```
j_e=sqrt(det(kappa^(-1) k_e)),
dA_e=j_e dmu_kappa.
```

This determinant is an endomorphism determinant, invariant under label
coordinates; no global coordinate frame is required. The common cut
gap proved in Section 2 keeps it in a compact positive-definite patch.
The determinant and positive square root are smooth there, so

```
e -> j_e in C2(e; joint C3),
(rho_b^2/2)^(3/2) <= j_e <=(3rho_b^2/2)^(3/2).        (12)
```

For `A=k^(-1)k1`, the exact derivatives are

```
j1/j=(1/2)tr A,
j2/j=(1/2)tr(k^(-1)k2)
                   +(1/4)(tr A)^2-(1/2)tr(A^2).       (13)
```

The inverse and all mixed matrix products use the actual full cut
metric. In a component frame `A` need not be symmetric; the invariant
traces in (13) are still the correct identities. If operator norms are
later evaluated they must use the corresponding self-adjoint
conjugation, not a silently changed auxiliary matrix norm.

Integration over the fixed compact label sphere is a bounded linear
operation. Thus actual cut area and specified fixed-label weighted
scalar/curvature integrals formed from (5), (8) or the specified
one-sided contractions and `j` are `C2` curves, uniformly in this
short proper interval. Their derivatives differentiate both the moved
field and density. A normalized quotient additionally requires its
chosen denominator to be nonzero; no experimental template, response
rank or normalization is selected by this theorem.

## 8. Actual remainders and precise endpoint of the result

For any of these proved `C2` curves in its stated Banach target `B`,
compact parameter continuity of the second derivative gives

```
||F_(a+h)-F_a-h F'_a-(h^2/2)F''_a||_B
 <=(h^2/2) sup_(|c-a|<=|h|) ||F''_c-F''_a||_B
 =h^2 o(1).
```

Segments must stay in the admitted analytical parameter neighborhood;
reported physical conclusions concern the original interval `I`.
This is an actual qualitative local-observer remainder. It is not a
numerical response constant or a third-parameter-derivative/cubic bound.

The new conclusion is a common actual short intrinsic proper cylinder,
positive spacelike cuts, full global injectivity on that cylinder, and
actual twice parameter-differentiable scalar, curvature and area
pullbacks with the displayed frame/label terms. Its duration `r_*` is
positive but unspecified. It does not include the historical
`25ell/16` scaled-clock cylinder, the `6145ell/4096` target, any
prescribed long proper duration, ambient free fall, shell-crossing
transmission, evaluated source/detector/noise norms or a cosmological
preparation mechanism. None of those is inferred from compactness.

## 9. Source binding

`SOURCE_PINS.json` pins the
exact immutable inputs at the stated base and their limited roles.
The two new companion proofs provide actual regularity/anchor and
actual common intrinsic ODE construction; their original frozen hashes are
recorded in AUDIT.md, with exact original-to-public promotion. The present
proof supplies the positive-cut, global-degree/graph and moving-geometry
steps rather than assuming the older conditional observer theorem's
domain hypotheses. No ancestor producer or checker is imported or run.
