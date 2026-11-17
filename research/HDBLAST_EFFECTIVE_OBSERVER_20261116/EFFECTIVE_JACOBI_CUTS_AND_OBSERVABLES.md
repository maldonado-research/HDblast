# Effective Jacobi cut gap and observer geometry through the full target cylinder

16 November 2026 (America/Los_Angeles; 17 November UTC). New proof at immutable base
`b4f2f59ab3918264c9832d87be1e39ed4dbfc178`. All earlier roots and
capsules remain sealed. The argument evaluates a sufficient cut/area
comparison for the actual family. It supplies no detector calibration,
arrival theorem or nonzero response lower bound.

## 1. Actual inputs, normalization and result

Keep the original admitted physical parameter interval `I`, the
original-time slab `[-T,0]`, and fixed background shell scale:

```
ell=10^-6,       T=13ell/8,       Q=25ell/16,
rho_b>0 fixed,
h_e=iota^*g_e,       iota(t,n)=(t,0,n),
barh_e=h_e/rho_b^2,
g0=-dt^2+cosh(t)^2 kappa,       p=dt^2+kappa.
```

Here `g0` denotes the normalized de Sitter SHELL metric in this proof;
it is not the full bulk metric. All tensor norms below are full
Frobenius norms for the positive product metric `p`, with every
covariant slot present. Matrix operator norms, when used, are labeled
`op`. The fixed reference connection is `D=D_p`.

The companion `EFFECTIVE_SHELL_METRIC_AND_CONNECTION.md` proves the
following for the actual
normalized metrics, uniformly on the whole existing shell slab and I:

```
k=barh-g0,
delta=sup_e sum_(m=0)^4 sup |D^m k|_F <10^-20,
|barh|op<=2,       |g0|op<=2,
|A0|_F<1,         |A|_F<2,
A0=nabla^(g0)-D,  A=nabla^(barh)-D,
|DeltaA|_F<=6delta,       |D DeltaA|_F<=48delta.        (1)
```

These are consequences of the actual common mixed solution bound and
the exact shell reconstruction; their pointwise hypotheses are not
postulated for a new solution. In particular delta is less than 1/4.
For the scalar the same companion proves the exact shell identity and
pointwise bound

```
varphi_e=phi_e|shell=phi_b+(v_e)_w,
sup |varphi_e-phi_b|<=10^7 Mfield,                     (2)
```

where `v_e` denotes the retained straightened sixteen-field solution,
and `Mfield` its common global `X14` bound. It is distinct from the
geodesic velocity V below.

The companion `EFFECTIVE_INTRINSIC_OBSERVER_CONTINUATION.md`
constructs the actual intrinsic normal
geodesics for every `0<=q<=Q`, with complete label spheres, strict
original-time retention, and actual parameter/label derivatives:

```
x_e(q,n)=(t_e,omega_e),       V=x_q,
barh(V,V)=-1,                 V past,
|V|_p<2,
-T+ell/32<t_e(q,n)<=-7q/8<0 for q>0,
(x,V) in C2(e; joint C4),      P=d_n x in C2(e; joint C3).
```

At q=0 the anchor is `(0,n)`, `V=-partial_t`, and the metric and its
initial time jet equal the background by the actual supported initial
collar. These facts, including parameter endpoint meanings, are the
phase companion's proved actual inputs. No cut gap is assumed there.

This q is the background-scaled proper past duration:

```
r=rho_b q,          x_r=V/rho_b=-u,
h(u,u)=-1,          u future.
```

Thus the effective cylinder below has physical proper duration
`0<=r<=rho_b Q`. It is not obtained by calling coordinate t proper
time. Normalized parallel screen vectors have barh-length one;
their physical h-unit versions are divided by rho_b.

We prove throughout `I x [0,Q] x S3` a normalized cut gap
`bar k_cut>=kappa/2`, global embedding of the whole cylinder, a relative
area error at most `2^18 delta q^2`, and the actual C2 geometric
continuation in the targets stated in Section 7. No new parameter
ceiling or PDE-slab restriction is imposed.

## 2. Curvature convention and the complete full-array error ledger

Fix the convention

```
R(X,Y)Z=nabla_X nabla_Y Z-nabla_Y nabla_X Z-nabla_[X,Y] Z.
```

For the background,

```
R0(X,Y)Z=g0(Y,Z)X-g0(X,Z)Y,       Ric(g0)=3g0.
```

This convention is checked below against the actual background Jacobi
field `cosh(q)I`, before assigning signs to any tidal matrix. Write
the mixed component array as `R^a_bcd`, meaning
`R(partial_c,partial_d)partial_b`, and lower its first index. Then

```
(R0)_abcd=(g0)_ac(g0)_bd-(g0)_ad(g0)_bc.
```

The reference product curvature cancels when subtracting the two
connection-difference formulas. Each antisymmetrization has two terms.
Contraction of full positive Hilbert arrays is bounded by the product
of their Frobenius norms. Thus

```
|R(barh)^(1,3)-R(g0)^(1,3)|_F
 <=2|D DeltaA|_F
       +2(2|A0|_F |DeltaA|_F+|DeltaA|_F^2)
 <=96delta+24delta+72delta^2
 <=138delta <=144delta.                             (3)
```

The factor two outside the quadratic expression charges the two
antisymmetric connection products; the inner factor two charges
`A0*DeltaA` and `DeltaA*A0`. The final inequality uses
`delta<=1/4`. No componentwise bound is relabeled as a full-array
bound, and no base-curvature term is discarded before subtraction.

In dimension four, `|Id|_F=2` and `|g0|_F<=4` follow from its operator
bound in (1). The constant-curvature formula therefore gives

```
|R0^(1,3)|_F<=2 |Id|_F |g0|_F<=16.
```

Lowering with the ACTUAL barh in the first term of a telescoping
difference gives

```
|R(barh)_(0,4)-R(g0)_(0,4)|_F
 <=|barh|op |DeltaR^(1,3)|_F+|k|op |R0^(1,3)|_F
 <=2*144delta+16delta=304delta.                     (4)
```

To compare with a constant-curvature tensor for the ACTUAL metric,
define `W(g)_abcd=g_ac g_bd-g_ad g_bc`. Since both full metric
Frobenius norms are at most four, direct telescoping of its two
products gives

```
|W(barh)-W(g0)|_F
 <=2(|barh|_F+|g0|_F)|k|_F
 <=16delta<=32delta.
```

As `R(g0)_(0,4)=W(g0)`, equations (3)--(4) yield the rounded interface

```
|R(barh)_(0,4)-W(barh)|_F
 <=336delta <=512delta.                            (5)
```

This sequence charges both the changed curvature array and the changed
metric used to lower the output. Comparing only to `W(g0)` would not
give the identity matrix in an actual barh-orthonormal screen.

## 3. Actual parallel screen and Jacobi coefficient

Choose a fixed smooth global round orthonormal frame `L_A` on S3,
for example quaternion multiplication fields. Along each constructed
geodesic solve actual normalized parallel transport

```
nabla^(barh)_q E_i=0,          E_i(0,n)=L_i(n),   i=1,2,3.
```

The equation in the product connection is
`D_q E_i=-A(V,E_i)`. It is a linear equation along the already retained
phase, not a new observer trajectory assumption. The p-compatible
connection and (1) give

```
|E_i(q)|_p<=exp(integral_0^q |A|_F |V|_p ds)
           <=exp(4q)<=exp(4Q)<2.                   (6)
```

For the last elementary bound, `4Q<1/2` and
`exp(x)<=1/(1-x)` for `0<=x<1` suffice. The initial p-norm is one.
Actual metric compatibility and normal initial data give

```
barh(E_i,E_j)=delta_ij,       barh(V,E_i)=0.
```

The full normalized congruence obeys Gauss orthogonality
`barh(V,P_A)=0`, either by its actual phase construction or by the
torsion-free/unit-normal proof in the retained observer theorem. Hence
the three-dimensional screen spans all label Jacobi fields, and

```
P_A=d_n x[L_A]=sum_i E_i S_iA,
S_iA=barh(E_i,P_A).
```

The labels need not be commuting angular coordinate generators. They
are fixed on the label sphere, commute with the q direction of the
parameter cylinder, and define the complete geometric differential.

Commuting the geodesic variation gives
`nabla_V^2 P_A=R(V,P_A)V`. In the parallel screen this is

```
S''=K S,
K_ij=barh(E_i,R(V,E_j)V)=R_(0,4)(E_i,V,V,E_j).
```

By (5), write `K=I+E`. Indeed
`W(barh)(E_i,V,V,E_j)=delta_ij` from unit timelike normalization and
screen orthogonality. Each error entry is bounded by
`512delta |E_i|_p |V|_p^2 |E_j|_p<=2^4*512delta`.
The full 3-by-3 Frobenius norm therefore has the bound

```
|E(q,n)|_F<=3*16*512delta<=2^15 delta=:eta.          (7)
```

Here three is the square root of the nine-entry count, sixteen is the
four vector-norm factors in (6), and 512 is the full tensor error (5).
This is a bound for all modes and all labels, not a diagonal tidal
approximation.

The actual initial pair gives

```
S(0)=I,             S'(0)=0.                        (8)
```

The second identity uses the initial TIME jet: the current shell metric
and first time jet equal the background, whose anchored cut is totally
geodesic at t=0. More explicitly, the supported normalized initial pair
has `v_e(0)=partial_t v_e(0)=0` on an open shell collar. Hence
`k(0)=D k(0)=0`; also `D g0(0)=0` since `cosh'(0)=0` and kappa is
parallel. Thus `A(0)=0`, and with the fixed initial phase velocity
`V=-partial_t`, `nabla_(L_A)V=0`. A prescribed background
initial position alone would not imply this identity.

As a sign check, for delta=0 the normal flow is
`x_0(q,n)=(-q,n)` and a parallel screen is
`E_i(q)=sech(q)L_i`. The actual differential is `P_A=L_A`, so
`S_0(q)=cosh(q)I`, with `S_0''=S_0`. This agrees with (7)--(8) and
rules out replacing the Jacobi equation by the opposite curvature sign.

## 4. Closed Volterra comparison and effective cut gap

The exact matrix equation is

```
S''-S=E S,
S(q)=cosh(q)I+integral_0^q sinh(q-s)E(s)S(s) ds.      (9)
```

First bootstrap `||S||op<=2`; do not assume `|S|_F<=2` in its place.
On any initial interval where this operator bound holds, use
`|ES|_F<=|E|_F||S||op` and
`sinh(q-s)<=2(q-s)` for `0<=s<=q<=Q<1/2`. Equation (9) gives

```
|S(q)-cosh(q)I|_F<=2^16 delta q^2=:d(q).            (10)
```

The elementary kernel bound follows, for example, by integrating
`cosh x<=2` on `[0,1/2]`. On the same interval
`cosh q<=1+q^2<=5/4`; this follows directly from its nonnegative
factorial series. With `delta<10^-20`,

```
d(q)<=2^16*10^-20/4<1/4,
||S(q)||op<=cosh q+d(q)<3/2<2.
```

The actual Jacobi matrix exists continuously throughout the retained
phase interval. A first point with operator norm two would contradict
this strict improvement. Thus the bootstrap closes on all `[0,Q]`,
and (10) holds there without a small unproved product constant.

The minimum singular value satisfies

```
smin(S)>=cosh q-||S-cosh(q)I||op
                >=1-d(q)>3/4.                     (11)
```

For the actual normalized cut metric, relative to the round label
frame,

```
bar k_cut=P^*barh=S^T S,
bar k_cut>=(9/16)kappa>=(1/2)kappa,
bar k_cut<=4kappa.                                 (12)
```

The physical cut metric is `k_cut=rho_b^2 bar k_cut`. This proves a
uniform effective positive cut gap for every q up to the SPECIFIED Q,
not merely for an unspecified shorter interval. The original I and
field-coordinate T remain unchanged.

## 5. Global embedding on the complete effective cylinder

All hypotheses of the retained Gauss/degree theorem are now actual
facts through Q: complete normal geodesics on full label spheres,
unit past normalized velocity V, Gauss orthogonality, the positive
cut gap (12), timelike `Kcoord=partial_t`, and identity angular anchor.
For completeness, the key implications are as follows.

If `D_n omega[v]=0` for nonzero label v, then
`P_v=(D_n t[v])Kcoord`; this has nonpositive barh-square, contradicting
(12). Hence each angular map `A_q:S3->S3` is a local diffeomorphism.
Compactness makes it a covering, and its homotopy from `A_0=Id` has
degree one, hence one sheet. It is therefore a global diffeomorphism.
Set `B_q=A_q^(-1)` and `Theta(q,m)=t(q,B_q(m))`.

At fixed m, `partial_q[x(q,B_q(m))]=Theta_q Kcoord` differs from V
by a label-tangent vector. Gauss orthogonality gives

```
Theta_q barh(Kcoord,V)=barh(V,V)=-1.                 (13)
```

The inner product is positive because Kcoord is future timelike and
V is past timelike. Thus `Theta_q<0` throughout the full connecting
interval. Two equal shell points must have the same q by (13), then
the same label by the angular diffeomorphism. This proves injectivity
between different clock cuts as well as on each cut.

The full differential is injective: contract `aV+P_v=0` with V to
obtain a=0, then use (12) for v=0. Compact-domain injection plus this
regular immersion gives an embedding up to both endpoint cuts. The
actual fixed shell inclusion iota is an embedding, so the physical
map `gamma=iota o x` is embedded in each one-sided bulk chart.

In physical h units the corresponding graph identity is
`Theta_q h(Kcoord,V)=-rho_b^2`. It is equivalent to (13), not a
different clock normalization. No coordinate-time inverse is clipped
or extended past a label-dependent clock image.

## 6. Effective relative area and scalar-mean bounds

Let `j` be the physical area density relative to the fixed round
label measure, and let the background comparison area be

```
j=rho_b^3 sqrt(det(S^T S))=rho_b^3 |det S|,
j_bg(q)=rho_b^3 cosh(q)^3,
A(q)=integral_S3 j dmu_kappa,
A_bg(q)=2pi^2 rho_b^3 cosh(q)^3.
```

The determinant of S is positive by (11), continuity and det S(0)=1,
so equivalently `j=rho_b^3 det S` on this entire cylinder.
For `H=S/cosh(q)`, (10) gives `||H-I||op<=d(q)`. Every singular
value of H lies in `[1-d,1+d]`, so

```
(1-d(q))^3<=j/j_bg<= (1+d(q))^3.
```

For `0<=d<=1/4`, both cubic deviations from one are at most `4d`:
`3d+3d^2+d^3<=4d`, and `1-(1-d)^3<=3d<=4d`.
Integrating the positive densities yields the explicit uniform bound

```
|A(q)/A_bg(q)-1|<=4d(q)<=2^18 delta q^2.             (14)
```

For the actual observed scalar, use the typed maps

```
varphi=iota^*phi,      Phi=varphi o x=phi o gamma,
B(q)=integral Phi j dmu_kappa,
Mobs(q)=B(q)/A(q).
```

The positive cut gap gives A>0. The exact scalar identity (2) holds
at every moved point of the retained physical shell, and hence

```
|Phi-phi_b|<=10^7 Mfield,
|B-phi_b A|<=10^7 Mfield A,
|Mobs-phi_b|<=10^7 Mfield.                          (15)
```

No derivative of the moving observer is needed for this VALUE bound;
the complete moving derivatives are still required for parameter
regularity. Bounds (14)--(15) are upper bounds only. They prove no
nonzero scalar, area or tidal response, no arrival of the initial
perturbation and no experimental sensitivity.

## 7. Actual C2 continuation and frames through the full Q

The phase companion gives actual `C2(e; joint C4)` phase and
`C2(e; joint C3)` label derivatives on all `[0,Q]`, using a retained
compact phase tube, finite continuation and Volterra/jet arguments.
Uniqueness of the same intrinsic geodesic initial-value problem
identifies this solution with the retained short observer map on their
overlap after `r=rho_b q`; the anchors and transported labels agree.
Parallel-screen uniqueness makes the same identification for frames.
This is continuation of those actual observers and their pullbacks.
This construction does not assume cut positivity. Conversely our
Jacobi proof uses that phase and now supplies the independent global
geometric gap (12). Its strict margin persists on a qualitative open
parameter neighborhood of I by compact continuity, so the inherited
ordinary endpoint parameter derivatives remain meaningful. No
numerical enlargement of I is reported.

The same finite linear transport argument as in the retained observer
proof extends the screens throughout Q. In fixed coordinates its
coefficient is the full connection evaluated on the actual phase;
in noncoordinate round/product frames the fixed frame-derivative
terms must also be included. Equivalently the invariant equation is
`D_q E=-A(V,E)`. The connection is `C2(e;C6)` and the phase is
`C2(e;C4)`, so the two-base-derivative reserve `6=4+2` gives the
coefficient `C2(e;C4)`. Linear Volterra differentiation supplies
actual first and second parameter derivatives in label C4, and the
transport equation recovers the joint q rows. Bound (6) is a finite
uniform screen bound over the whole effective cylinder.

All moving-evaluation and tensor/label product formulas proved in the
retained observer capsule now apply on this FULL domain, with the same
safe strong targets:

```
phase, future velocity and screens:      C2(e; joint C4),
Phi=varphi o x:                         C2(e; joint C4),
curvature evaluated/pulled back:         C2(e; joint C3),
full cut metric and j:                  C2(e; joint C3),
A, B and Mobs:                         C2(e; C3([0,Q])).
```

For example if `Y=x_epsilon,Z=x_epsilon,epsilon`, the second moved
field derivative contains
`f_ee+2D f_e[Y]+D^2 f[Y,Y]+D f[Z]`. The cut-metric Hessian also
contains `d_nY,d_nZ` and the two-label-first-variation term. Every
velocity/frame slot of a curvature contraction varies. Curvature
composition consumes two of the fixed-base C5 spatial derivatives,
leaving C3; an unchanged top curvature target is not asserted.

In physical units set `u=-V/rho_b` and `Ephys_i=E_i/rho_b`.
For each copy's actual outward unit normal, the physical Israel law
gives

```
a^g=nabla^g_u u=(sigma/6)n_out,
nabla^g_u Ephys_i=0,       nabla^g_u n_out=(sigma/6)u.
```

Thus the actual one-sided orthonormal screen
`(Ephys_1,Ephys_2,Ephys_3,n_out)` remains Fermi--Walker along the
confined observer, with equation
`nabla_u W=u g(a^g,W)-a^g g(u,W)`. Its q or past-r evolution carries
the corresponding negative/scale factor. The metric-normal expression
is a smooth actual state function on the retained positive radial
patch. Consequently the full one-sided physical tidal components,
including normal, velocity and frame variations, remain `C2(e;C3)`
through Q. The identities use the physical shell law only on I;
differentiate in its interior and extend to endpoints by continuity.
No physical law outside I, ambient free fall, or smooth crossing spray
is being assumed.

The original-time edge uses the accepted fixed finite-regularity
coefficient extensions only to justify open-domain ODE/composition
rules. The actual strict inward retention makes final values and jets
the original past ones. All conclusions are independent of an auxiliary
future extension. The strengthened domain does not require an X15
estimate or strong top-norm general-data dependence.

For every displayed C2 curve the Banach Taylor remainder is
`h_parameter^2 o(1)` uniformly on compact admitted parameter segments.
No third parameter derivative, cubic bound or evaluated numerical
second-derivative modulus follows from this statement.

## 8. Source boundary

`SOURCE_PINS.json` records the exact
immutable source inputs and their limited roles. The new effective
metric and phase proofs supply the actual quantitative interfaces (1)
and full-Q phase existence; their frozen hashes are recorded with the
final joint source ledger. The previous observer capsule supplies
already proved finite-regularity composition, proper normalization,
physical-frame and global Gauss/degree arguments. Its formerly
qualitative proper radius is not substituted for Q here: equations
(3)--(15) and the effective phase proof establish the new sufficient
whole-Q domain.
