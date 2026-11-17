# A clock-and-degree route to a global observer domain

5 November 2026, UTC. Independent conditional argument. No PDE estimate,
amplitude, nonlinear constant, or physical family is evaluated here.

## 1. Scope and orientation

The input is an already existing regular observer map, including its parameter
dependence, on a closed proper-clock tube

```
X_e(q,n)=(t_e(q,n),omega_e(q,n)),
|e|<=a, 0<=q<=T_q, n in S^3,   q=-tau/rho_b, rho_b>0.
```

The target reference slab is the product cylinder `(-T_coord,0] x S^3`.
Coordinates and fields must be available on a neighborhood used for the
candidate flow, or the estimates must be applied on a stopped continuation
tube before its first reference-boundary exit. Such stopping means the first
common/global stopping clock: every earlier q-slice must still contain all
labels in S3. Labelwise stopped partial domains do not support the covering
argument below. This product description
uses the sphere as a manifold, not one angular chart.
The fixed initial anchor section and fixed transported labels are

```
X_e(0,n)=(0,n),       X_0(q,n)=(-q,n).
```

If the physical anchor is a different fixed embedded section, first choose
a fixed product identification putting it in this form, and charge that
identification in all bounds. A parameter-dependent reidentification is not
implicit here. A moving anchor requires its own extra terms and is outside
this stated theorem.

The physical future unit velocity is `u=partial_tau X`, so

```
X_q=-rho_b u,           h_e(X_q,X_q)=-rho_b^2.
```

Thus past integration requires `t_q` negative in the background and in the
clock guard below. Coordinate time is not declared equal to proper time.
All derivatives used as parameter inputs hold the actual `q` and transported
label fixed. The coordinate-clock inverse introduced later is a device for
this proof, not a replacement for the observable's proper-clock derivative.

Once Gauss normality is supplied, uniform C1 regularity suffices for the
global injectivity argument; use sufficient geodesic regularity to derive
that identity and C2 regularity for the inverse derivative formulas.
Existence of this map on an
open neighborhood of the closed tube, regular shell fields and their physical
identification remain upstream assumptions. The argument can certify that
the map stays inside the reference cylinder and is globally one-to-one.
On a stopped tube the strict containment guard rules out that boundary
exit; a compact regular state tube and ODE continuation must still be
provided. The argument does not create geodesics or replace that estimate.

The full shell metric is `h_e=iota_e^*g_e`, with the actual physical shell
map and all its parameter dependence. The injectivity proved here is that
of `X_e` INTO THE REFERENCE SHELL cylinder. For the physical composition
`iota_e o X_e` to be globally injective in ambient spacetime, the retained
physical shell identification `iota_e` must additionally be an embedding
on the containing reference-shell domain. A merely immersed shell map or
a fixed-coordinate restriction does not supply this premise; the pullback
metric identity alone does not imply global ambient injectivity.

## 2. Main physical theorem: normality and a timelike coordinate vector

There is a stronger route for the ACTUAL normal-geodesic observer map.
Assume throughout the already existing proper-clock tube that

```
h_e has signature (-,+,+,+),
K=partial_t is future timelike,
U=X_q=-rho_b u is past timelike, h_e(U,U)=-rho_b^2,
each induced cut metric k_q=X(q,.)^*h_e is positive definite.
```

The initial section is fixed and embedded, and the initial observer velocity
is its actual unit normal. Initial normality is not replaced by background
normal data. For a transported label variation `P=D_n X`, the geodesic
equation and the torsion-free connection give, in each label direction,

```
partial_q h_e(U,P)
 =h_e(nabla_U U,P)+h_e(U,nabla_U P)
 =0+h_e(U,nabla_P U)
 =(1/2)partial_n h_e(U,U)=0.
```

The anchor normality therefore implies `h_e(U,P)=0` on the entire tube.
This is the Gauss orthogonality identity. The factor `rho_b` is fixed in
the label, so the proper-norm derivative is zero. This proof concerns one
fixed physical e; no parameter differentiation is needed.

For each q, suppose `D_n omega(v)=0`. Then
`P(v)=(D_n t(v))K`. If this vector is nonzero, it is timelike; if zero,
its squared length is zero. Both possibilities contradict
`k_q(v,v)>0` for nonzero v. Thus `D_n omega` is invertible. Compactness,
the local diffeomorphism property and the homotopy from the anchored
identity imply that `A_q(n)=omega(q,n)` is a degree-one covering, hence
a diffeomorphism. This proof is global on S3 and uses no global chart.

Set `B_q=A_q^{-1}` and `Theta(q,m)=t(q,B_q(m))`. The derivative at fixed
angular coordinate is

```
V=partial_q X(q,B_q(m))
 =U-P(D_n omega)^(-1)omega_q
 =Theta_q K.
```

Gauss orthogonality now gives the exact identity

```
Theta_q h_e(K,U)=h_e(V,U)=h_e(U,U)=-rho_b^2.
```

The Lorentzian inner product of a future timelike and a past timelike
vector is positive. Hence `h_e(K,U)>0` and `Theta_q<0`. Each q-cut is
therefore a graph whose time height strictly decreases at fixed angular
coordinate. Equality of two reference-shell points forces equality of q by this
monotonicity and then equality of their labels by the angular
diffeomorphism. The ENTIRE finite proper-clock cylinder is injective,
including points from different proper-clock cuts.

This physical theorem needs no smallness of `D_n omega-I`, no Schur
inequality and no per-label sign condition on `t_q`. It needs the actual
normal, proper, geodesic conditions, positive cut metrics and a timelike
coordinate vector on the whole tube. Positive cut metric alone does not
state those other inputs. In particular a generic smooth-map counterexample
below does not contradict this theorem.

For a quantitative separation estimate, any certified
`0<h_e(K,U)<=B` gives `-Theta_q>=rho_b^2/B`. On the compact tube B is
finite when the fields are regular; evaluating B needs actual metric and
velocity budgets. For example an auxiliary positive norm yields the
sufficient bound `B<=||h_e||op ||K|| ||U||`. To certify K timelike from
a background margin, it suffices to have
`h_0(K,K)<=-gamma_0<0` and
`|(h_e-h_0)(K,K)|<=eta_h gamma_0` with `eta_h<1`, together with the
continuous future-orientation choice. No such constants are evaluated here.

Since `Theta(0,m)=0`, the physical theorem also gives `t(q,n)<0` for
all q>0. The displacement guard `delta_t<T_coord-T_q` below then places
the tube in the required reference slab. Flow existence, positive cut
metrics and timelike K must still be obtained by the prior ODE/metric
bootstrap; injectivity is not assumed in their place.

## 3. Quantitative alternative inputs, with manifold sphere comparisons

Use the unit round sphere embedded in R4. For `v in T_n S^3`, the notation
`I_n v=v` is its ambient inclusion. Suppose on the whole tube

```
|t+q| <= delta_t,
mu=1-eta_t >0,        mu <= -t_q <= 1+eta_t,
||D_n t|| <= alpha,   ||omega_q|| <= beta,
|D_n omega(v)| >= sigma |v|,       sigma>0.       (1)
```

Here `D_n omega:T_n S^3 -> T_omega S^3`; norms are the round norms or
their Euclidean ambient realizations. These definitions require no global
frame or angular chart. One sufficient angular input, often directly
available from the label-Jacobi variables, is

```
||D_n omega-I_n||_(T_n S^3 -> R4) <= epsilon <1,
sigma=1-epsilon.                                  (2)
```

The lower singular bound in (1) is the real input; (2) is only a convenient
sufficient budget. There is no requirement for a global ambient displacement
or a geodesic-to-chord factor in this route.

The two strict guards are

```
delta_t < T_coord-T_q,
alpha beta < mu sigma.                            (3)
```

Set `c=mu-alpha beta/sigma>0`. Units are consistent: `t,q,alpha,delta_t`
have the common coordinate units fixed by `t_0=-q`, `beta` has inverse
coordinate units, and `mu,sigma,c` are dimensionless. The conversion is
`[rho_b]=[tau]/[q]`; if tau and rho_b both have physical length units,
q (and the background t coordinate) is dimensionless. These coordinate
units must not be silently identified with physical proper duration.
The product norms include the chosen
reference scaling, atlas comparisons and embedding constants, if any.

## 4. Quantitative global injectivity across proper-clock cuts

Fix a parameter and a proper clock q. The angular map
`A_q:S^3 -> S^3`, `A_q(n)=omega(q,n)`, has invertible differential by
(1), so it is a local diffeomorphism. Compactness makes its image closed;
local invertibility makes it open. Its image is therefore all of the
connected sphere. A proper local diffeomorphism is a finite covering.

The continuous family `A_r`, `0<=r<=q`, is a homotopy from `A_0=Id`.
Consequently `deg A_q=1`. A local diffeomorphism of a connected oriented
manifold has one constant local orientation sign, and its degree is that
sign times its number of sheets. Degree one forces one sheet. Hence every
`A_q` is a diffeomorphism. This is a manifold argument; it does not identify
the entire sphere with one coordinate domain. The implicit function theorem
also makes `A_q^{-1}` jointly C1 in `(q,m)`.

Write each proper-clock cut as a graph over the fixed angular coordinate:

```
B_q=A_q^{-1},           Theta(q,m)=t(q,B_q(m)).
(B_q)_q=-(D_n omega)^(-1) omega_q,
Theta_q=t_q-D_n t (D_n omega)^(-1) omega_q.         (4)
```

The bound on the inverse is `||(D_n omega)^(-1)||<=1/sigma`. Thus

```
Theta_q <= -mu+alpha beta/sigma=-c<0.              (5)
```

If `X(q1,n1)=X(q2,n2)=(t,m)`, then
`n_i=B_(q_i)(m)` and `Theta(q1,m)=Theta(q2,m)`. Strict monotonicity in
(5) gives `q1=q2`; the angular diffeomorphism then gives `n1=n2`.
This proves injectivity of the FULL reference-shell map, not merely injectivity
on each proper-clock cut. Its derivative is nonsingular as well: the
Schur complement identity in any oriented local frames is

```
det DX=det(D_n omega)
        [t_q-D_n t (D_n omega)^(-1) omega_q].      (6)
```

For a compact tube the continuous injection is a homeomorphism onto its
image; the nonsingular derivative gives the corresponding smooth embedding
up to the boundary. When the given map is the physical normal-geodesic
map, this supplies a sufficient global noncrossing/label-domain condition.
It asserts neither an evaluated input nor a new field-equation theorem.

Finally `t_q<=-mu` and the fixed anchor give `t(q,n)<=-mu q<0` for q>0.
The residual bound and the first guard in (3) give
`t(q,n)>=-q-delta_t>=-T_q-delta_t>-T_coord`.
The entire proper-clock tube lies in the stated open-past reference slab.

## 5. Coordinate-clock inversion and the fixed-time angular degree argument

For each label, `q -> t(q,n)` is strictly decreasing. The inverse
`Q(t,n)` exists exactly when `t(T_q,n)<=t<=0`; its label domain can
depend on t. On every time slice where all labels are present, define

```
C_t(n)=omega(Q(t,n),n),
D_n Q=-D_n t/t_q,
D_n C_t=D_n omega-omega_q (D_n t)/t_q.            (7)
```

All quantities on the right are evaluated at `(Q(t,n),n)`. The rank-one
correction obeys

```
|D_n C_t(v)| >= [sigma-alpha beta/mu]|v|>0.        (8)
```

The same covering/degree proof now applies directly to `C_t`. The path
from t=0 to the chosen t stays in the common time domain, and `C_0=Id`;
therefore `deg C_t=1` and it is a diffeomorphism. This is the requested
fixed-coordinate-time angular-map proof on complete time slices.

There is an essential endpoint issue: the full proper-clock cylinder has
the generally nonconstant past edge `t(T_q,n)`. Below
`max_n t(T_q,n)`, `Q(t,n)` is only defined on a subset of labels. Degree
of a sphere map cannot be invoked there without an extension argument.
Proving (8) only on the common coordinate-time slab does not prove global
injectivity of the rest of the proper-clock tube. The graph argument
(4)-(5) closes exactly this gap, without extending the geodesic solution,
clipping an inverse, or silently charging no cost for such an extension.

A guaranteed common strip is

```
-(T_q-delta_t) <= t <=0,
```

since `t(T_q,n)<=-T_q+delta_t` for every n. This is a sufficient strip;
the actual common strip can be larger.

## 6. Exact slack for the requested tube

With `ell=10^-6`, use

```
T_coord=13 ell/8,   T_q=25 ell/16,
q_star=6145 ell/4096,
T_coord-T_q=ell/16=256 ell/4096,
T_q-q_star=255 ell/4096.                           (9)
```

Thus `delta_t=a Z_t<ell/16` is enough for reference-slab containment.
If one also needs the fixed coordinate time `t=-q_star` to be present for
every label, require the slightly stronger
`delta_t<255 ell/4096`. These two margins differ by `ell/4096`;
they must not be conflated. Full proper-clock injectivity itself needs no
target common-slice condition.

Anchoring and the derivative bound also imply `|t+q|<=eta_t q`, so one
may replace `delta_t` by the sharper of independent `a Z_t` and
`eta_t T_q`. If this derivative estimate is the only displacement input,
the exact sufficient clock guards are respectively

```
eta_t<1/25                       (slab containment),
eta_t<255/6400                   (target common slice).      (10)
```

These are conditional inequalities on unevaluated budgets, not a numerical
amplitude radius. To use the parameter variational estimates from the
observer-flow note, projections of the full parameter/label derivatives
can certify `a Z_t`, `eta_t`, `epsilon`, `alpha` and `beta`. The normal's
initial data and all chart/embedding conversions must be included. No
claim is made here that any of (1)-(3) has been checked for HDBLAST.

## 7. Inverse-clock and C2 costs

Let `f(e,q,n)=t_e(q,n)` and use a local label direction only to write the
following tensor identities. For a base variable `y=(e,t,n)`, put
`F(y,q)=f(e,q,n)-t`. At `q=Q(y)` one has

```
Q_i=-F_i/F_q,
Q_ij=-(F_ij+F_iq Q_j+F_jq Q_i+F_qq Q_i Q_j)/F_q. (11)
```

For label Hessians use the round covariant Hessian; the scalar implicit
identity is unchanged. In particular,

```
Q_t=1/f_q,  Q_e=-f_e/f_q,  D_n Q=-D_n f/f_q,
Q_ee=-(f_ee+2f_eq Q_e+f_qq Q_e^2)/f_q,
Q_tt=-f_qq/f_q^3.                                (12)
```

Hence, writing `|f_e|<=B_e`, `|f_ee|<=B_ee`,
`|f_eq|<=B_eq`, `|f_qq|<=B_qq`,

```
|Q_t|<=1/mu,  |Q_e|<=B_e/mu,  ||D_n Q||<=alpha/mu,
|Q_ee| <= B_ee/mu+2B_eq B_e/mu^2+B_qq B_e^2/mu^3. (13)
```

For `Omega(e,t,n)=omega_e(Q(e,t,n),n)`,

```
Omega_e=omega_e,par+omega_q Q_e,
Omega_ee=omega_ee,par+2omega_eq Q_e
          +omega_qq Q_e^2+omega_q Q_ee.           (14)
```

Thus bounds `C_e,C_ee,C_eq,C_qq` on the corresponding explicit omega
derivatives and `beta` on `omega_q` give

```
|Omega_e|<=C_e+beta B_e/mu,
|Omega_ee|<=C_ee+2C_eq B_e/mu+C_qq B_e^2/mu^2
             +beta[B_ee/mu+2B_eq B_e/mu^2+B_qq B_e^2/mu^3].
```

These denominators are real costs of coordinate-clock inversion. Equations
(11)-(14) must not be inserted as if they were derivatives at fixed proper
clock. Ambient second derivatives use the fixed sphere embedding; a chart
conversion requires the corresponding transition Hessians.

## 8. Explicit failures of insufficient hypotheses

On the flat real cylinder `R x S1` with `h=-dt^2+dtheta^2`, define the
smooth map

```
F(q,theta)=(-q(1-sin(theta)/2),theta+4q mod 2pi),
0<=q<=1.
```

Its anchor is `F(0,theta)=(0,theta)`, `-t_q` lies in `[1/2,3/2]`,
and every angular map at fixed q is a rotation diffeomorphism. The induced
cut metric is

```
k_q=(1-q^2 cos(theta)^2/4) dtheta^2 >=(3/4)dtheta^2.
```

Nevertheless there is the explicit collision

```
F(pi/4,pi/2)=F(pi/8,pi)=(-pi/8,3pi/2).
```

Both clocks lie in `[0,1]` because `0<pi<4`. This demonstrates that
positive cut metric, per-label monotone coordinate time, degree-one
angular maps on each proper-label cut, and a fixed embedded anchor still
do not by themselves imply global noncrossing. The missing Schur sign is
the substantive extra condition in (3). This is NOT a normal-geodesic or
proper-clock observer counterexample: its q-velocity is spacelike, since
`h(F_q,F_q)=16-(1-sin(theta)/2)^2>=55/4>0`. It makes no claim about
HDBLAST's field equations.

An even simpler fold-return map is `F(q,theta)=(-q(1-q),theta)`:
all cut metrics are exactly `dtheta^2`, while
`F(q,theta)=F(1-q,theta)`. This one also loses clock monotonicity.

There is a genuine unit-normal-geodesic return on a flat Lorentzian torus:
`F(q,theta)=(-q mod 1,theta)`. It returns after one unit with positive cut
metrics. Its periodic time coordinate violates the required real cylinder
and has closed timelike curves. It is only a topology caveat, not an
example in the reference domain.

For the actual timelike normal congruence with a future timelike coordinate
vector, section 2 proves the graph sign exactly. Positivity of the induced
metric in isolation does not state those other hypotheses. The quantitative
alternative instead uses explicit coordinate budgets.

## 9. Reproducibility and limitations

`check_clock_domain.py` uses only Python's standard library and exact
Fractions. It checks the universal Gauss-contraction algebra and clock
orientation scaling, the universal polynomial Schur determinant identity,
the inverse-clock second-jet and composed angular second-jet identities,
the rank-one perturbation guard algebra, the exact slab/target margins,
strict-boundary guard rejection, and exact rational multiples of pi in
the explicit collision. The supplied positive guard example is synthetic
algebra and is not an HDBLAST amplitude or evaluated flow budget.

It does not algorithmically prove the covering-space theorem, compactness,
analytic norm inequalities, inverse function theorem, geodesic existence,
or any PDE premise. Those are the mathematical arguments and conditional
inputs above. No nonlinear budget is estimated by this finite checker.

```
python3 check_clock_domain.py --output /ABS/FRESH/normal.json
python3 -O check_clock_domain.py --output /ABS/FRESH/optimized.json
cmp /ABS/FRESH/normal.json /ABS/FRESH/optimized.json
```

The mandatory output path must be absolute, fresh and non-symlink, with
an existing non-symlink parent chain. Exclusive creation forbids overwrite.
The complete canonical JSON is also written to stdout, identically in
normal and optimized Python. The manifest excludes itself and pins every
deliverable plus the two read-only source notes. No old file was modified,
legacy checker replayed, network request issued, or publication performed.
