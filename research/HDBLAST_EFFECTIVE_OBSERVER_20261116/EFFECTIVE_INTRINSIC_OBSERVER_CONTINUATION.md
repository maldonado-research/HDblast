# Effective continuation of the actual intrinsic observers through Q

16 November 2026 (America/Los_Angeles; 17 November UTC). New proof at immutable base
`b4f2f59ab3918264c9832d87be1e39ed4dbfc178`. All earlier private roots
and public capsules are read-only. Source identities and companion
bindings are recorded in `SOURCE_PINS.json`.
This proof continues the already constructed actual normal observers;
it does not replace their short branch by a formal variation solution.

## 1. Actual inputs, normalization and conclusion

Retain the selected original canonical interval and coordinate slab

```
I=[-a_iter,a_iter], a_iter=10^(-10^500209),
ell=10^-6, T=13ell/8,
Q=25ell/16, q_star=6145ell/4096,
T-Q=ell/16, Q-q_star=255ell/4096.                    (1)
```

On the actual fixed physical shell `[-T,0] x S3` write

```
h_e=iota^*g_e, iota(t,n)=(t,0,n),
gbar_e=h_e/rho_b^2, gbar_0=-dt^2+cosh(t)^2 kappa,
p=dt^2+kappa, D=D_p,
A_e=nabla^(gbar_e)-D, A0=nabla^(gbar_0)-D.           (2)
```

The constant rho_b is the retained fixed shell scale, independent of e.
The connection is unchanged by this constant metric normalization.
All tensor norms in this proof are full p-array norms, with every
input, output and derivative slot included. No angular truncation is
made. The original one-copy or identical-copy convention is retained.

The metric companion proves, from the actual admitted field bound and
exact nonzero-state shell reconstruction, a uniform number delta with

```
0<=delta<10^-20,
sup_(e in I) sum_(m=0)^4 sup_(t,n) |D^m(gbar_e-gbar_0)| <=delta,
|A_e-A0|<=6delta, |D(A_e-A0)|<=48delta,
|A0|<12T<1, |D A0|<14.                             (3)
```

In particular `|A_e|<2`, `|D A_e|<15`. The sufficient inequality
delta<10^-20 is a consequence of the unchanged registered amplitude,
not a new restriction imposed on e. The companion derives its stronger
bound from `delta<=10^11 M` and the retained common M. No general
reduced-PDE coefficient envelope is assigned to this geometric
connection.

PR147 supplies an actual common short intrinsic normal flow in physical
past proper duration r, with future unit velocity u,
`x_r=-u`, `nabla_r u=0`, actual anchor `(0,n0)` and
`u(0)=rho_b^-1 partial_t`. Its unaugmented phase is C2 in e into
joint C4, and its actual label augmentation is C2 into joint C3.
The current field theorem supplies an open reduced-family parameter
extension and `h_e` C2 in e into joint C7. Consequently A_e is C2
into C6. These qualitative regularity premises are retained separately
from the evaluated small field estimates (3).

Set

```
q=r/rho_b=-tau/rho_b,
V=x_q=-rho_b u=(a,w),
gbar_e(V,V)=-1.                                    (4)
```

We prove that the actual branch continues for every
`(e,n0) in I x S3` through the full common interval `0<=q<=Q`.
It remains in the original past slab and satisfies

```
|a+1|+|w| <=48delta q,
|t+q| <=24delta q^2,
|n(q)-n0|_(R4) <=24delta q^2,
|a+1|+|w|<1/8,
-T+ell/32<t(q)<=-7q/8<0 for q>0.                    (5)
```

The lower t bound also holds at q=0. The phase is an actual C2 parameter
curve in joint C4(q,n0), and its actual label derivative table is C2 in
joint C3 on this full interval. Ordinary two-sided derivatives at the
endpoints of I are supplied by a qualitative open parameter neighborhood
with the same Q. No numerically enlarged parameter ceiling is evaluated.

The physical common proper duration is therefore `rho_b Q`; q_star
is included with its strict interior gap in (1). This proves phase
continuation without assuming positive later cuts, absence of focal
points, global injectivity or a nonzero response. Those geometric and
observable conclusions belong to the separate Jacobi/cut proof.

## 2. Exact global sphere equation and constraints

Use once the unit sphere embedding `S3 subset R4`. Along a true shell
phase let `n in S3`, `w in T_n S3`. The actual normalized geodesic
equation, with V past directed, is

```
t_q=a,                  n_q=w,
a_q=-A_e^0(V,V),
w_q=-|w|^2 n-A_e^S(V,V),
(t,n,a,w)|q=0=(0,n0,-1,0).                         (6)
```

The angular output A^S is embedded in R4. The sphere term has the
displayed minus sign because `n_q=w`. This convention differs from
the r system written using future u. It follows either from (4) and
the accepted phase equation or directly from
`D_q V=-A_e(V,V)`.

For an ambient state `(t,omega,a,W)` near the constraint, put

```
n=omega/|omega|, Pi_n=Id-n tensor n, v=(a,Pi_n W),
F_e(t,omega,a,W)=
 (a,W,-A_e^0(t,n)[v,v],-|W|^2 omega-A_e^S(t,n)[v,v]). (7)
```

The field is smooth enough on the fixed annulus `1/2<|omega|<3/2`
and the actual open past time domain. All angular descriptions are
global; there is no single angular chart or uncharged chart transition.
The same formula applies near the anchored time using PR147's fixed
finite-regularity coefficient extension when a local open neighborhood
is needed. Effective estimates below are used only on t<=0, where
these are the actual field coefficients.

The A angular output in (7) is orthogonal to omega. Thus, for
`c1=|omega|^2-1`, `c2=omega.W`, the exact constraint equations are

```
(c1)_q=2c2,       (c2)_q=-|W|^2 c1.                 (8)
```

Their initial values vanish. Elementary uniqueness for this homogeneous
linear pair gives c1=c2=0. Therefore all continued solutions of (7)
from the actual branch remain in the true tangent-bundle phase space,
where (6) is precisely the actual geodesic equation.

The normalized background terms are explicit, with no curvature sign
convention needed:

```
A0^0(V,V)=cosh(t)sinh(t)|w|^2,
A0^S(V,V)=2tanh(t) a w.                            (9)
```

These are the time/angular difference coefficients for
`gbar_0=-dt^2+cosh(t)^2 kappa` relative to p.

## 3. Stopped estimates for the actual past phase

Shorten the already proved common short branch if necessary so the
following stops are initially strict. Continue it as an actual local
solution until Q or the first exit from

```
-T<t<0 for q>0,
|a+1|<1/4,       |w|<1/4.                          (10)
```

No label-Jacobi or cut condition is included in this phase stop. Every
solution is considered with all of its original labels, and Section 5
will show that the same common horizon exists for their whole compact
set. On (10),

```
|V|_p^2=a^2+|w|^2 <=25/16+1/16=13/8<2.
```

Since a<=-3/4, time points strictly inward. Equation (3) bounds each
time/angular perturbing acceleration by `12delta`. For the upper
right Dini derivative of `e(q)=|a+1|+|w|`, use (9) and
`|cosh(t)sinh(t)|<=2T`, `|tanh(t)|<=T` on `|t|<=T`.
For w!=0 the sphere term in (6) is orthogonal to w, and consequently
makes no contribution to the derivative of |w|. At w=0 its value is
zero and the upper Dini bound remains valid. Thus

```
D^+ |a+1| <=2T|w|^2+12delta,
D^+ |w| <=2T|a||w|+12delta,
D^+ e <=3T|w|+24delta <=e+24delta.                 (11)
```

The coefficient `3T` charges `2T(1/4)+2T(5/4)`; `3T<1`.
All contractions in the perturbing terms use the full norm of A;
projecting a time or angular output block costs at most one. The sum
in (11) charges both output blocks, rather than dropping either.

The continuous e is locally absolutely continuous. Its integral
Gronwall inequality, or the elementary Dini comparison obtained after
multiplying by exp(-q), gives

```
e(q)<=24delta(exp(q)-1)<=48delta q,      0<=q<=Q.    (12)
```

For the last step, Q<2ell<1/2 and exp(Q)<2, so
`exp(q)-1=int_0^q exp(s) ds<=2q`. No smallness of a generic
uncomputed ODE coefficient is used. Integrating t_q+1=a+1 and
n_q=w yields

```
|t(q)+q|<=24delta q^2,
|n(q)-n0|_(R4)<=24delta q^2.                        (13)
```

The entire estimates are valid before a proposed first stop. They
therefore improve that stop strictly. From delta<10^-20 and Q<2ell,

```
48delta Q <96*10^-26 <1/8,            since 768<10^26,
24delta Q^2 <96*10^-32 <ell/32,       since 3072<10^26. (14)
```

In particular `(a,w)` lies in the improved bounds in (5). They imply
`a<=-7/8` and

```
t(q)<=-7q/8,
t(q)>=-Q-ell/32=-T+ell/32.                         (15)
```

The strict gap T-Q=ell/16 is charged in the last equality. Thus no
velocity stop or past coordinate-time exit can occur at or before Q.
The future time face cannot be an exit for q>0 because a remains
negative. This calculation is independent of the angular label and e.

## 4. Explicit point-state and Lipschitz bounds in one fixed global norm

For clarity we give a fully charged ambient local continuation bound.
Use the dimensionless q phase and the fixed block-sum norm

```
|(dt,domega,da,dW)|_oplus=|dt|+|domega|+|da|+|dW|,
```

where angular norms are the ambient Euclidean norms. This is one fixed
positive auxiliary norm; no physical measurement norm is inferred.
Consider the larger compact range

```
1/2<=|omega|<=3/2, |a|<=3/2, |W|<=1/2, -T<=t<=0.
```

Set H=A_e(t,n)[v,v], with n,v as in (7). Then |v|<2,
|v|^2<=5/2<3, `|A|<2`, and `|D A|<15` by (3). We claim

```
|dH| <=45|dt|+110|domega|+8|da|+8|dW|.             (16)
```

Here H is the full time/angular vector in its ambient representation.
This derivative is not identified with D A alone. Its terms are:

* Time coefficient variation costs at most `15*3=45`.
* For angular variation, `|dn|<=2|domega|`. The covariant
  coefficient derivative costs `15*3|dn|<=90|domega|`.
* Covariantly differentiating a projected fixed ambient input W on the
  sphere gives `D_(dn)(Pi_n W)=-(n.W)dn`. Its norm is at most
  `|W||dn|<=|domega|`. The two input slots therefore cost at most
  `2*2*2|domega|=8|domega|`.
* The angular output's ambient normal derivative costs at most
  `|H^S||dn|<=2*3*2|domega|=12|domega|`.
* A velocity-input variation da or dW costs at most
  `2|A||v|` times that variation, bounded by eight. Pi_n is an
  orthogonal projection, so no additional dW factor is used.

The angular sum is `90+8+12=110`. Every state normalization,
input slot, output slot and sphere representation derivative is
therefore counted. Summing the two acceleration output blocks costs
at most a further factor two, since `|H^0|+|H^S|<=2|H|`.
The derivatives of the sphere term are bounded by
`(1/4)|domega|+(3/2)|dW|`. Adding the two kinematic blocks in (7)
gives

```
|D F_e|_(oplus->oplus)
 <=max(90,220+1/4,17,16+3/2+1)<256,
|F_e|_oplus<3/2+1/2+12+3/8<16.                    (17)
```

These are bounds on local ambient balls contained in the displayed
range. Convex balls of sufficiently small fixed radius around true
phases stay in the annulus; a straight segment between arbitrary
antipodal sphere points is not used. Equations (16)-(17) supply an
explicit finite local Lipschitz and speed bound for continuation.
There is no unexplained generic constant such as 10^100.

## 5. Actual common continuation through Q

First consider the actual short branch from PR147, in q units. By its
uniform continuity on the compact parameter/label set, choose a common
`0<q0<=Q/2` within that branch with (10) strict. Let qbar be the supremum of
common clocks up to which every e in I and every initial label has
the actual continued phase satisfying (10), capped at Q. We have
qbar>=q0>0. Estimates (12)-(15) hold on all earlier clocks.

Suppose qbar<Q. On clocks `q>=q0/2` near qbar, the actual state values
lie in the single compact set

```
-T+ell/32<=t<=-7q0/16,
|omega|=1, |a+1|<=1/8, |W|<=1/8, omega.W=0.         (18)
```

Its time values are strictly inside (-T,0), velocities are strictly
inside the larger bounds in Section 4, and angular values have a fixed
annular margin. The choice q0<=Q/2 makes this a nonempty relevant
time range whenever continuation is still needed.

The speed bound sixteen makes all phase values uniformly Cauchy as
q approaches qbar from below. Their limits lie in (18). Choose a common
positive radius d so every block-norm ball of radius d around that compact
set stays in the original open time/annular/velocity domain. The finite
local contraction argument for (7) then works with, for example, any
step theta>0 satisfying

```
16theta<d/2,       256theta<1/2.                    (19)
```

It applies to every endpoint limit and every parameter/label. Thus the
actual phases extend for a common positive step beyond qbar, matching
the earlier solutions by local uniqueness. The improved stop (14)-(15)
persists for a sufficiently short part of this step, contradicting the
definition of qbar. This proves qbar=Q and actual continuation through
Q. At Q the same compact-set argument also gives local continuation
past Q as a phase, if needed; no larger target is asserted.

For completeness, the extended point flow is jointly continuous in
parameter, initial label and q. On each of the uniform local balls,
subtracting its integral equations gives the usual integral inequality
with Lipschitz constant 256 and the uniform coefficient-parameter
modulus. The initial map is smooth in the label. Gronwall gives
continuous dependence, and uniqueness agrees on overlaps. The common
uniform step and compactness yield continuity of the full glued flow.
In particular the endpoint limits used above are a common whole-label
map, not merely separate labelwise continuations. The global sphere
equation and (8) ensure no angular chart domain is lost.

This proves the phase domain without any inverse label map, cut metric
positivity, Jacobi inverse or noncrossing assumption. A focal point does
not itself prevent the phase ODE from continuing; cut positivity is
proved separately.

## 6. Actual C2 dependence and label regularity on the full interval

The effective constants above control phase continuation. The parameter
argument first applies at interior points of I, with continuous one-sided
limits at its endpoints; Section 7 supplies the open neighborhood for
ordinary two-sided endpoint derivatives. Higher
coefficient/parameter derivatives needed for qualitative C2 dependence
are finite on the compact phase range, by the retained metric C2 C7
theorem. They need not be numerically evaluated to prove regularity,
and do not change (12)-(17).

At the anchored time face one may reuse PR147's fixed bounded C7
coefficient extension for the open-neighborhood calculus. Its future
derivatives are only qualitative finite bounds, not inputs to the
effective past phase estimates. All actual continued curves have t<=0;
the resulting derivatives and flow are independent of this auxiliary
extension by actual local uniqueness. At all later clocks the phase
range lies in the actual open past domain.

Here are explicit finite-regularity arguments for the continued targets.
On each compact local range, F_e is C2 in e into C6 in state. Its
Nemytskii map on `E4=C4(S3;R10)` is C2: four label derivatives and
two Frechet state derivatives require at most six state derivatives.
The ordered finite product/Taylor proof is exactly the one proved for
the actual short branch. No higher PDE solution-flow theorem is an input.

Actual label derivatives of orders at most four initially exist from
that branch. Whenever local continuation has been made, differentiating
the actual equation in ordered fixed label directions gives a triangular
system for each label jet J_m:

```
(J_m)_q=(D_Y F_e)J_m+R_m,
```

where R_m is a finite sum involving lower label jets and state derivatives
of F of order at most m. The reference label connection is fixed;
q commutes with it, and angular derivative orders are not commuted.
On the compact phase range these coefficients are uniformly bounded.
For m=1 the remainder has no higher jet; induction and the integral
Gronwall inequality give finite uniform bounds for all m<=4 on [0,Q].
Their q derivatives are also bounded by the same equations. Hence these
actual jets have limits at any continuation endpoint, so continuation
does not lose C4 label regularity. This proves a continuous whole-label
solution in `C([0,Q];E4)`.

The same finite product bounds give finite Nemytskii derivative bounds
on a common bounded E4 neighborhood of the actual graphs. For nearby
parameters subtract the actual integral equations on [0,Q]. Gronwall
in E4 first proves parameter continuity and a local Lipschitz estimate.
The exact first secant equation has the segment coefficients

```
p_h=int_0^q [B_h p_h+b_h] ds,
B_h=int_0^1 D_Y mathcalF(e+theta h,Y_e+theta DeltaY)dtheta,
b_h=int_0^1 D_e mathcalF(e+theta h,Y_e+theta DeltaY)dtheta. (20)
```

At the anchored face the common coefficient extension makes these open
segments legitimate; between the actual member graphs they remain in
a common compact annular neighborhood. The E4 coefficients converge
uniformly to `B=D_Y mathcalF(e,Y_e)`, `b=D_e mathcalF(e,Y_e)`.
The linear Volterra operator has a unique solution on any finite
interval: iterated integrals are bounded by `(L4 Q)^j/j!` for a finite
operator bound L4, so its Neumann series converges. Subtracting the
two integral equations and applying Gronwall identifies p_h with a
limit U, the actual first derivative. No assumption `L4 Q<1` is made.

The coefficients B,b are then C1 curves, and actual secants of U converge
by the same Volterra estimate to the unique second derivative Z:

```
U_q=F_Y U+F_e,
Z_q=F_Y Z+F_YY[U,U]+2F_Ye U+F_ee,
U(0)=Z(0)=0.                                      (21)
```

Continuity of U,Z follows by comparing these integral equations and
using the continuous second derivatives of the Nemytskii map. This
constructs actual C2 parameter dependence into C_q E4 on the continued
domain. The finite compact-state and lower-jet bounds justify every
local patch used in this argument; no qualitative compactness argument
has been used to identify the former short duration with Q.

Recover joint q/label regularity from the exact ODE. Define
`G1=F`, `G_(p+1)=D_Y G_p F`, p<=3. Actual time jets obey
`partial_q^p Y=G_p(e,Y)` for p<=4. Requesting 4-p label derivatives
and two parameter derivatives needs at most

```
(p-1)+(4-p)+2=5<=6
```

state derivatives of F. Each row is therefore C2 into
`C_q C_n^(4-p)`, and the actual integral identities identify the compatible
joint C4 table. Bounded label differentiation yields the actual C2 joint
C3 augmentation `(x,V,d_nx,d_nV)` on [0,Q] x S3. Changing back to
physical r=rho_b q and u=-V/rho_b is a fixed linear rescaling, preserving
these regularity conclusions.

## 7. Endpoint parameter neighborhood and preserved geometric identities

The original-I phase estimates have strict uniform margins: the actual
phase lies inside the larger compact annular/velocity range, and its past
coordinate edge has gap ell/32. The underlying metric family is already
C2 on an open parameter interval J. Ordinary local continuous dependence
and the uniform whole-label compact range therefore give, around every
e in I, an open parameter neighborhood whose flows exist through the
same Q with slightly weaker strict margins. Compactness of I selects
finitely many such neighborhoods. Their union contains an open interval
J_eff around I; shrink it qualitatively within J if necessary.

The C2 argument of Section 6 applies on that open interval. Initial
normal data stay fixed there by the retained initial collar. Negative
a preserves inward t<0 for q>0. This supplies ordinary two-sided
derivatives at both endpoints of I, without a new numerical epsilon
ceiling. No new full Einstein–scalar physical claim outside I follows
from the use of this reduced-family extension.

The actual connection is metric compatible. Since (6) is its geodesic
equation, `gbar(V,V)=-1` remains true from its exact initial value.
Torsion-free commutation of q and the transported labels gives

```
partial_q gbar(V,d_nx)
 =gbar(V,nabla_(d_nx)V)
 =(1/2)d_n[gbar(V,V)]=0.
```

The initial contraction is zero, so Gauss orthogonality holds throughout
[0,Q]. These identities are propagated on the actual continued flow;
they do not depend on an assumed embedded cut or invertible label map.
The physical proper-clock identity is `h(x_r,x_r)=-1`, and the future
unit velocity remains `u=-V/rho_b`. The observer stays intrinsic to
the shell and is generally accelerated relative to ambient bulk free
fall, as in PR147.

## 8. Exact scope

The continuation covers every original parameter and full S3 label
through the prescribed Q and hence q_star, on the same original
coordinate field slab and amplitude interval. Its effective point-state
estimate is supplied by the literal normalized shell metric and its
fully charged connection, not by an invented large generic ODE bound.
The separate intrinsic Jacobi comparison may now prove positive cuts,
global embedding, full moving curvature/area/scalar pullbacks and
transported-frame identities on this domain.

This theorem does not prove that a disturbance reaches either target,
that a scalar or curvature response is nonzero, or that an observable
moment has a sign or a zero. It supplies no detector/noise/nuisance
calibration, physical energy normalization, third parameter derivative
or cosmological origin mechanism. The former inner causal shadow and
all separate scientific limitations remain in force.
