# Explicit prepared family: all five initial harmonic means

1 November 2026. Portable conditional derivation. The new initial curve is the
reviewed explicit B-mean preparation [D], not the old weighted-only curve or
the older explicit physical scalar correction. Public sources are read only;
their exact bytes and the reviewed preparation inputs are pinned in
`SOURCE_PINS.json`. The operator and forcing are retained, while the initial
mean displacement and velocity below are new.

**Result.** Equations (8)–(9) give all ten initial amplitudes, and (10) gives
their scaled radial displacement derivative. The exact-family identification
in section 7 applies simultaneously to the five coupled means. It does not
evaluate their subsequent response or import an old initial norm bound.

## 1. Coordinates, normalization, and the mean projector

The retained background and first seed [S1, sections 1–2] are

```
g0=dy²+rho(y)²[-dt²+cosh(t)² kappa_S3],
a=rho_y/rho, p=phi0_y, dy=-rho dz,
ell=10^-6, x=t/ell, r=z/ell,
s=ut/rho, j=vt/rho, Y=n1²-n2², <Y²>=1/6.
```

Thus `partial_x=ell partial_t`, `partial_r=-ell rho partial_y`, and `t`
is the unscaled de Sitter coordinate time. It is not shell proper time;
background shell proper time is `rho_b t`. In particular the initial velocity
below is `partial_x U`, with no extra `rho_b` division. The background-frame
radial derivative conversion is `partial_y=-rho^-1 partial_z`. The opposite
display `dy=-rho^-1 dz` in [S5] is a prose typo, contradicted by the common
operator and its derivative conversion; it is not used here.

All parameter expansions are `X_e=X0+e X1+e² X2/2+o(e²)`. The canonical
data, physical data, metric `k=g2`, and scalar `w=phi2` therefore denote
second parameter derivatives. Every extraction of an `e²` coefficient below
includes multiplication by two. A mean is integration against the probability
measure on the unit S3. In tensor language it is Haar averaging under the
natural SO(4) action, not averaging arbitrary coordinate components.

In the radially parallel background frame

```
e_y=partial_y, e_0=rho^-1 partial_t,
e_a=(rho cosh t)^-1 E_a,
P0 k=K0 dy²+2rho L0 dy dt+rho²N0 dt²
       +rho²cosh(t)²J0 kappa,       P0 w=w0.                (1)
```

The registered mean vector is `U=(K0,L0,N0,J0,w0)` with Gram diagonal
`(1,2,1,3,1)` [S3, section 1; S4, retained problem]. `N0=k(e_0,e_0)` is
the covariant time-time component; there is no Lorentzian sign reversal in
this definition. All invariant angular one-form and angular trace-free tensor
means vanish. The single Dirichlet amplitude is `L0`.

## 2. New canonical mean and constraints

The reviewed preparation [D,R] fixes `b=Q`, the same compact flat bump as
the first packet. With increasing outward y and `yL=y(2ell)<yU=y(ell)`,

```
I=integral_yL^yU rho⁴ QH dy,  J=integral_yL^yU rho² b dy>0,
f=2I b/(9J),
C[f]=-3(f_yy+3a f_y)/2-V f,
xi=-(rho⁴p²)^-1 integral_yL^y rho⁴(QH+C[f])dy,
A=f+2xi_y,                    B=f+2a xi,
A_y=f_y+2xi_yy,               B_y=f_y+2a_y xi+2a xi_y.      (2)
```

Here `A,B` are metric profiles, not the correction map or bulk matrices.
The canonical mean is
`<gamma2>=diag(A,B,B,B)`, `<pi2/rho0>=diag(0,D,D,D)` with
`D=-QM/(6a)`, and `phi2=varpi2=0` in all angular sectors. The volume density
`rho0` is distinct from the radius `rho`.

For reference, the source formulas are fixed functionals of the unchanged
first data, with every prime denoting `partial_y` [W, section 2]:

```
QH=-2u_y²-24a u u_y-2u u_yy
   -(16a_y+32a²+8p²/3)u²-4u²/rho²-s²
   -(j²+v_y²+V''v²)/6-(4/3)p u v_y-(4/3)v²/rho²,
QM=-4u s_y-12a u s-2s u_y+j v_y/3
   =-2s u_y+j(v_y+4pu)/3.
```

Their constraint normalization is `QH=<D²H[v_seed,v_seed]>` and
`QM=<D²(M/rho0)_y[v_seed,v_seed]>`. The linearized metric and momentum rows
are

```
Pg(A,B)=-3B_yy/2+3a A_y/2-6a B_y
          +(3a_y+6a²+p²/2)A-3B/rho²,
Pm(C,D)=-2[C_y+3a(C-D)].
```

The retained radial identities give
`Pg(f gamma0)=C[f]` and
`Pg(Lie_(xi partial_y)gamma0)=rho^-4(rho⁴p²xi)_y`.
Since `integral rho⁴C[f]dy=-(9/2)integral rho²fdy=-I`, the primitive
in (2) has zero total integral. It extends flatly by zero at both endpoints.
Thus `Pg(A,B)=-QH` and `Pm(0,-QM/(6a))=-QM`, exactly as required by
`P x2+D²Chat[v_seed,v_seed]=0`. An extra half on the right would be wrong.
The retained bounds `rho>=129`, `|p|>=3/25000`, `a>=111/1000` justify all
displayed divisions on the support. These are initial constraints, not
Einstein accelerations prescribed independently of the reduced evolution.

## 3. Canonical-to-physical conversion

The action and momentum normalization are those of [S1, section 2]:

```
K_ij=(2/rho_gamma)[gamma_ia gamma_jb pi^ab
                       -gamma_ij gamma_ab pi^ab/3],
Pi=varpi/rho_gamma,             K=Lie_n gamma/2.
```

At `pi0=Pi0=0`, no `gamma2` enters the second physical velocity. The
quadratic inverse-map terms depend only on the unchanged first tangent:

```
K2=2[tau2-I tr(tau2)/3]+diag(36us,-6us,-6us,-6us)Y²,
Pi2=-2uj Y².                                                  (3)
```

Consequently, writing `R2=<K2_yy>` and `S2=<K2_angular>` in the fixed frame,

```
R2=QM/(3a)+6us,      S2=-us,       <Pi2>=-uj/3,
R2+3S2=QM/(3a)+3us.                                           (4)
```

This uses `tau2=pi2/rho0` with contravariant components expressed in the
fixed background orthonormal frame; it is not `pi2/rho_gamma` at the
perturbed metric. Confusing the two removes the quadratic terms in (3).

The selected lapse and shift **values** are `N_e=rho(1+e uY)` and
`beta_e=0`. Therefore

```
<k_tt>=-rho²u²/3, <k_tI>=0, <k_yy>=A,
<k_ij>=rho²B kappa_ij, phi2=0,
<dot gamma2_yy>=2rho[QM/(3a)+16us/3],
<dot gamma2_angular>=-(4/3)rho us,
dot phi2=0.                                                   (5)
```

The angular velocity in (5) is the fixed-frame coefficient, to be multiplied
by `rho² kappa` for coordinate angular components. The scalar cancellation
is pointwise: `dot phi2=rho Pi2+2rho uY Pi1=0`, since
`Pi1=jY`. Thus the scalar zero initial pair is exact in every angular sector.

## 4. Exact relative-harmonic completion and its second mean

Let `C^A=g^{BC}(Gamma[g]^A_BC-Gamma[g0]^A_BC)` and
`R^A=g^{BC}Gamma[g0]^A_BC`. The reference is fixed but its contraction uses
the perturbed inverse metric. At zero shift the exact prescription [S2]
is

```
dot N=N² tr_gamma K-N³R^t,
dot beta^I=N²(Gamma_gamma^I-R^I)-N gamma^IJ partial_JN,
dot gamma=2NK,                 dot phi=N Pi.                  (6)
```

Here `R^t=0` on the initial slice. For transparency, take general second
diagonal spatial profiles `A2,B2` and physical curvature `R2pt,S2pt`, with
first fields `U1=uY`, `T1=utY`. Set
`Aeps=1-4eU1+e²A2/2`, `Beps=1+2eU1+e²B2/2`, `P=1+eU1`.
The exact radial shift expression is

```
dot beta^y=rho²P²[Aeps_y/(2Aeps²)-3Beps_y/(2Aeps Beps)
              -4a/Aeps-P_y/(Aeps P)+a/P²+3a/Beps].
```

Expanding with the factorial convention gives

```
dot N1=rho T1,             dot beta1^y=-6rho²(U1_y+4aU1),
dot beta1^i=0,
dot N2=rho²(R2pt+3S2pt)-24rho U1 T1,
dot beta2^y=rho²[A2_y/2-3B2_y/2+4aA2-3aB2
                              -74U1 U1_y-194aU1²],
dot beta2^i=22U1 D^iU1-(D^iA2+D^iB2)/2.
```

For general non-diagonal second tensors the omitted linear terms are the
corresponding linear harmonic operator. Their invariant projection is zero
unless they lie in the radial/isotropic block. Equivariance of (6) at the
background therefore lets us use the displayed diagonal formulas to compute
the full mean, without asserting that the full new second datum is diagonal.
With `H=A_y/2-3B_y/2+4aA-3aB`, (4) implies

```
<dot N2>=rho²[QM/(3a)-us],
<dot beta2^y>=rho²[H-(37/3)u u_y-(97/3)a u²],
P0(dot beta2^angular)=0.                                     (7)
```

Metric conversion must retain the first-order products:
`dot g_tt=-2N dot N` and `dot g_ty=gamma_yy dot beta^y`.
Their second derivatives contain `-4N1 dot N1` and
`2gamma1_yy dot beta1^y`. In particular

```
<dot k_tt>=-2rho³QM/(3a)+(4/3)rho²u ut,
<dot k_ty>=rho²[H-(13/3)u u_y-(1/3)a u²].
```

Omitting the mixed metric product would incorrectly leave coefficients
`-37/3,-97/3`; importing the old completion would omit all `A` terms.

## 5. Registered five-mean initial API

At `t=0`, the time derivative of every angular `cosh(t)` frame factor is
zero. Radial frame factors depend only on y. Dividing the mixed coordinate
component by `rho`, the time-time component by `rho²`, and the angular
component by `rho²` therefore gives

```
U0=(A, 0, -u²/3, B, 0),                                     (8)

V0=partial_x U|x=0
  =ell*(
       2rho QM/(3a)+(32/3)u ut,
       rho*[A_y/2-3B_y/2+4aA-3aB-(13/3)u u_y-a u²/3],
       -2rho QM/(3a)+(4/3)u ut,
       -(4/3)u ut,
       0).                                                   (9)

partial_r U0=-ell rho*(A_y,0,-2u u_y/3,B_y,0).                (10)
```

`mean_initial.py` implements these formulas without importing any old
initial-profile evaluator. It accepts exact or interval arithmetic inputs.
The two shape values and their y derivatives are enough for the map. To
construct them from (2) one needs `xi` through two y derivatives and `f`
through one; differentiating the primitive directly uses `QH_y` and `f_yyy`.
Stable equivalent primitive formulas can avoid numerical cancellation, but
must preserve (2) and its zero endpoint condition.

As an algebraic cross-check independent of `QM`, the physical trace identity
from (9) is
`V0_K+V0_N+6V0_J=4ell u ut`; it includes the lapse/time-sign conventions.
Every component of (8)–(10) vanishes on the initial open shell collar and
outside the closed packet support `[1,2]` in r. The full nonlinear family
need not have this smaller support: its unprojected corrections and higher
parameter derivatives can extend across the larger correction annulus.

## 6. Initial gauge, constraints, and corner implications

The new polynomial trial `x0+e v_seed+e²Bmean/2+Acor U(e)` is treated by
the retained weighted IFT exactly as in [D,R]. Its derivative in U is the
same invertible operator. Haar averaging commutes with the radial operators
and weights, preserves the dotted closures, and maps the reflection-invariant
domain into itself. Since `P Bmean=-P0 D²Chat[v_seed,v_seed]`, injectivity
gives `P0 U''(0)=0` and `P0 x''(0)=Bmean`. This reasoning needs no claim
that every individual rotation preserves the reflection-invariant subspace.

The exact canonical curve satisfies all constraints. Equations (6) impose
`C=0` for every parameter on the initial slice, with no modification of
physical data. The reduced Einstein accelerations and the physical
Hamiltonian/momentum constraints imply `n(C)=0` there: the constraints are
the normal projections of `E-g tr(E)/2`, not the single equation `E_nn=0`
[S1, section 4; S2, section 5]. Differentiating twice gives `C2=dot C2=0`;
lower-order normal/frame variations multiply already vanishing subsidiary
data. These accelerations are supplied by the reduced PDE, not additional
initial inputs. The exact open zero shell collar supplies every corner
compatibility order needed by the retained finite-regularity theorem.

This establishes the initial constraints and harmonic gauge analytically.
The finite exact checks below corroborate the conversion and Christoffel
algebra; they do not replace the retained weighted-IFT or shell-existence
premises.

## 7. Conditional exact-family identification of all five means

Retain the same background branch, first master/gauge fields, common
fixed-shell chart, identical copies, and prescribed-slab theorem as [S1].
Adding Bmean changes only smooth compact initial second data and the local
IFT constants. After a possibly family-dependent parameter shrinkage, the
same strict common domain arithmetic holds: `Lsupport=5ell/2`,
`Tphys=3ell/2`, `Z=6ell`, `c*=2`, and `Lsupport+c*Tphys<Z`.
The initial gauge completion is (6). The entire first metric/scalar jet is
unchanged, including its mixed gauge velocity, so the retained first linear
uniqueness identifies the exact family's first derivative with the actual
common-harmonic first pair throughout the common causal domain.

The finite-regularity parameter theorem gives a C2 exact family into X10
with sufficient lower C4 spacetime control [S1, section 6]. Parameter
differentiation, the necessary spacetime derivatives, sphere integration,
and physical-side boundary traces therefore commute. Differentiating the
exact reduced bulk and shell equations twice gives

```
Lbulk(k,w)=-D²Ered[(q,eta),(q,eta)],
Lshell(k,w)=-D²B[(q,eta),(q,eta)],
```

with exactly the retained actual-first-field quadratic sources. The
linearized background operators and the linear shell/gauge rows are
SO(4)-equivariant. Their natural action commutes with Haar averaging on
arbitrary regular tensor sections. Thus the five invariant slots form a
closed subsystem even if the full nonlinear family or its second derivative
has additional angular components. All invariant second metric tensors and
scalars have precisely the form (1); no vector or tensor slot is missing.

Projecting the complete retained operator [S3, equations (1)–(8)] gives

```
U_xx-U_rr=A0 U_x+B0 U_r+D0mean U+F0,
L0|r=0=0,   P_N U_r=M0 P_N U+G0,                            (11)
```

where `F0=P0(-2ell²rho² Qred,ell²rho² Qphi)` and `G0` is the same
projected physical/gauge boundary forcing as in [S3]. The factors two and
signs are those of its reduced residual convention. The trace gauge source
is still `GN=-2tr_gamma QR-2(Q_C)_y`. Neither these sources nor these
operators contain the changed second initial data. In mean order the
homogeneous outward Robin derivative rows are

```
K_y=-4aK+2a(-N+3J)+(4sigma'/3)w,
N_y=-aK-sigma'w/3,     J_y=aK+sigma'w/3,
w_y=pK/2-sigma''w/2,
```

and `U_r=-ell rho U_y`; `L` has its Dirichlet row. The forcing parts are
`K_y:GN`, `N_y:-2QR00`, `J_y:-(2/3)sum_a QRaa`, `w_y:-QP` before
conversion to r. These are frame components. In particular no coordinate
`2a` self-coupling is added to the tangential frame rows.

Construct the unique mean response with (8)–(9), the actual retained
operators/sources, and the reflected free wave `Ufree_new`:

```
Unew=Ufree_new+K(A0 Unew_x+B0 Unew_r+D0mean Unew+F0)
                      +L(M0 P_N Unew|r=0+G0).               (12)
```

The reflection is odd for L and even for K,N,J,w. The retained contraction
and smoothness proof in [S3] applies to any such smooth compatible finite
initial pair; the contraction constant depends on the homogeneous operator,
not the size of the chosen data. On the smaller backward cone of [S4], the
mean zero wedge `x+r<1` is preserved, and its sharpened bound
`3409/8000000<1` applies to the mean operator. This is a second-derivative
statement; it does not claim that all higher orders of the exact family
vanish on that wedge.

Both the projected exact family and (12) solve the same five-component
linear initial/boundary problem. Subtracting cancels the quadratic sources
and all initial data. Retained mixed linear uniqueness, localized under the
same strict domain-of-dependence margins, gives

```
P0(partial_e²(g_e,phi_e)|e=0)=Unew                            (13)
```

in all five slots throughout the common domain. This is a new application
of old uniqueness to new data. It asserts no equality with the old total
five-mean solution. The auxiliary zero-initial forcing response remains the
same, while its homogeneous-initial contribution changes. Old initial
caps, the failed old accuracy gate, and any total bound based on the old
scalar correction cannot be transferred without a fresh estimate.

## 8. Checks and scope

`checks/check_initial_mean.py` forms inverse metric and all Christoffels
directly at three nonspecial rational one-jets with independently varied
radial/tangential second values and derivatives. It checks all five
relative-harmonic components through order two, including angular gauge
velocities, then checks the canonical conversion, frame map and scalar
cancellation. Mutation controls deliberately remove the radial second
profile, lapse-square value, metric mixed product, or parameter factorial.
The checker also verifies the listed source hashes. Repository-independent authentication uses the explicit relative paths and
hashes, without requiring Git or a particular checkout HEAD.
Normal and optimized reports record the same exact mathematical checks.

This work gives new explicit initial data, a harmonic completion, and a
conditional exact-family/linear-response binding. Formula (12) defines the
evolved mean response but does not numerically evaluate it. No shell-volume
sign, heating history, matter-production mechanism, radiation, force, or
finite-parameter remainder is inferred. Internal AI-assisted checking is
not external specialist review.

Source keys: D = reviewed explicit mean design; R = its independent review;
W = prior canonical mean/source derivation, used for source algebra only;
S1 = exact-family TT transfer (20261031); S2 = explicit initial harmonic
completion (20261015); S3 = common-gauge full response construction (20261016);
S4 = total-response preflight manuscript (20261023); S5 = old initial-gate
proof, used solely to locate the registered mean ordering, never its data.
