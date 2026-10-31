# Conditional exact-family identification of the registered postcontact TT response

31 October 2026. Conditional analytic checkpoint. This note identifies a specified
second derivative, conditional on the retained finite-regularity analytic
theorems and actual-solution premises. It neither computer-certifies those
theorems nor changes the already published numerical calculation.

The new statement is that the actual scalar TT solution used in the
postcontact calculation is the selected angular projection of the second
parameter derivative of an exact family with a precisely selected initial
curve, gauge and common domain. The initial curve need not have the same
scalar second data as the registered explicit second jet. The initial TT
pair does agree, and the linear TT projection decouples from those other
second data. This suffices for the registered slicing observer.

## 1. Fixed model, branch, coordinates and derivative convention

Use the five-dimensional Einstein–scalar residuals and the one-copy outward
shell laws

```
E_AB = Ric_AB - phi_A phi_B - (2/3)V(phi) g_AB,
S = Box_g phi - V'(phi),
K_out,ab = sigma(phi) h_ab/6,  n_out(phi) = -sigma'(phi)/2.
```

There is one common tension and two identical interior copies, each with its
own outward normal. Put

```
Wsuper(phi)=1-phi+phi^3/3,
V=(Wsuper')^2/2-(2/3)Wsuper^2,
sigma=2Wsuper+delta(1+c phi),
delta=1/1000, c=96714321425692/161839258930731.
```

Writing eta=phi-1 gives exactly the retained polynomial
`V=-2/27+(14/9)eta^2+(50/27)eta^3-eta^4/6-(4/9)eta^5-(2/27)eta^6`.
Thus the explicit current potential and the shell theorem's model coincide;
this is not a change of scalar branch or action normalization.

Select the same near-phi=+1 regular member used by the current actual-jet
certificate, with its fixed shell data, not an arbitrary member of the
model. Its background is

```
g0=dy^2+rho(y)^2[-dt^2+cosh(t)^2 kappa_S3], phi0=phi0(y),
a=rho_y/rho, p0=phi0_y,
rho_y=a rho, p0_y=V'-4a p0,
a_y=-rho^-2-p0^2/3, V=p0^2/2-6a^2+6rho^-2.
```

The branch theorem and wider-collar first-exit enclosure are retained
analytic/background inputs. On `0<=t<=3ell, 0<=z<=8ell` they give smooth
coefficients, `129<=rho<=131` and
`-14/10^5<=p0<=-12/10^5`; in particular the divisions below are legitimate.
Here `ell=10^-6`, `dy=-rho dz`, `x=t/ell`, and `r=z/ell` increases inward.
The unscaled t is the de Sitter coordinate time; the notation Tphys below
distinguishes it from scaled x and does not denote elapsed proper time.
The symbol `r` below always denotes this scaled inward coordinate. In the
ambient initial chart write `R=|X|=b exp(-z)` instead, so `s=R^2`,
`rho=Omega R`, and `gamma0=Omega(s)^2 delta` on the initial section.

All parameter coefficients are derivatives:
`g_e=g0+e q+e^2 k/2+o(e^2)`, `phi_e=phi0+e eta1+e^2 zeta/2+o(e^2)`.
No additional one-half is attached to `k` or to the screen second derivative.
Use the ambient quadrupole matrix `diag(1,-1,0,0)`. The old notation
`Y=n1^2-n2^2` and the current zero-based notation `Y=n0^2-n1^2` describe
these same ordered first two axes; no rotation or orientation reversal is
performed. Keep the registered outward sphere orientation.

## 2. The complete actual first canonical seed is admissible

Let `Q(z)` be the smooth flat bump

```
Q=exp(4-ell^2/((z-ell)(2ell-z))) for ell<z<2ell,
Q=0 otherwise,
alpha=(rho/rho_b)^(3/2) p0/p_b,
kappa_r=alpha_z/alpha=rho(5a/2-V'/p0),
u=alpha Q/rho^2, ut=exp(-z)u,
v=3alpha(Q_z+kappa_r Q)/(rho^3 p0),
vt=3alpha exp(-z)(Q_z+(kappa_r-1)Q)/(rho^3 p0).
```

These are the actual master solution's initial data, not its comparison
function. In a background orthonormal spatial frame, with radial slot first,
the physical first Cauchy datum is

```
gamma1=diag(-4uY,2uY,2uY,2uY),
K1=diag(-2S1,S1,S1,S1), S1=ut Y/rho,
phi1=vY, Pi1=vt Y/rho.
```

The time extrinsic curvature convention is `K=Lie_u gamma/2` for the future
normal `u`, distinct from `K_out` at the timelike shell. Use the checkpoint
action `R/2-|dphi|^2/2-V` and canonical normalization

```
pi^{ij}=rho_gamma (K^{ij}-gamma^{ij} tr K)/2,
varpi=rho_gamma Pi,
K_ij=(2/rho_gamma)[gamma_ia gamma_jb pi^{ab}
                  -gamma_ij gamma_ab pi^{ab}/3].
```

Thus `pi1/rho0=diag(-3S1/2,0,0,0)` in that frame. Here a density is first
divided by the fixed background volume form `rho0`; its coordinate
component contains the sphere volume determinant and is not itself a radial
scalar. The shell gauge companion instead uses the action multiplied by two,
and hence `pi_comp=2pi`, `varpi_comp=2varpi`. Its allowed consistent
checkpoint normalization produces exactly the same physical data. No
companion conversion formula is used without this rescaling.

There is an explicit match to the earlier initial theorem's ADM seed, not
merely an agreement of angular matrices. That theorem uses

```
q_master=G(s)(X1^2-X2^2), p_master=Kprofile(s)(X1^2-X2^2),
h=2q_master(s delta-3XX^T),
delta phi=-3(D q_master+4q_master)/(2Omega^2 Phi_s),
delta pi=-3p_master XX^T/(2Omega),
delta varpi=-3Omega(D p_master+5p_master)/(2Phi_s), D=R partial_R.
```

Choose

```
G=Omega^2 u/s^2=alpha Q/s^3,
Kprofile=Omega^2 ut/s^(5/2)=G/b.
```

Since `R=b exp(-z)`, `D=-partial_z` and
`Phi_s=p0 rho/(2R^2)`, direct substitution gives all four displayed physical
first sectors, including
`delta varpi=Omega^4 Pi1` in the ambient chart. In particular the current
exp(-z) velocity factor maps to the constant ambient profile ratio `1/b`.
The first linearized constraints hold by the current all-point physical
binding, using `v=-3(u_y+2au)/p0` and its time derivative. The seed is smooth,
flat, compactly supported in `[ell,2ell]`, and invariant under each ambient
coordinate reflection in every canonical sector. Scalar and density
transformations preserve that invariance.

The exact localized theorem quantifies over **every** such smooth invariant
linear-constraint datum. Its concrete consequence explicitly permits
arbitrary smooth radial `G,Kprofile`, and states that no additional frozen
radial-profile condition is required. We apply that theorem anew to the
displayed seed. We do not identify an old unspecified nonlinear curve with
this new one.

## 3. Select the weighted curve and a common shell-interaction domain

Choose the correction annulus and domain

```
A'={ell/2<z<5ell/2} x S3,
Lsupport=5ell/2, Tphys=3ell/2, Z=6ell, c*=2.
```

The annulus is smooth and concentric also in the ambient initial chart,
contains the entire seed support strictly, and has closure away from shell
and vertex. Select radial defining function and radial weights as in the
exact initial theorem, rescaling the defining function so its scaled cover
balls stay inside the annulus. At the single finite choice `k=12`, let

```
x_e=x0+e v_seed+A u(e),
A_H N=w^2 qweight^4 P_g^*N,
A_M Zeta=rho0 w^2 qweight^2 P_m^*Zeta.
```

Here the scalar and scalar-momentum sectors are prescribed affinely and
`u(0)=u'(0)=0`. Retain the theorem's local adjoint classification,
invariant-kernel removal, scalar and symmetric-divergence API estimates,
scaled-cover estimates, dotted-space inverse, smooth nonlinear mapping and
zero-extension conclusions. They yield a smooth canonical exact constraint
curve in `H14 x H13`, with full correction support contained in `closure A'`.
The initial scalar and scalar-density second derivatives are zero for this
curve. Its second data need not agree with all the registered local physical
second sectors; their equality is unnecessary for the TT result below.

The required causal arithmetic is strict:

```
Lsupport+c*Tphys=11ell/2<6ell=Z,
(t_dagger,z_dagger)=(5ell/4,ell/8),
z_dagger+c*t_dagger=21ell/8<Z,
x<=3/2<3 and x+r<=3/2+6=15/2<8.
```

Thus the whole proposed rectangle fits both the retained background collar
and the actual first-field triangle. At background radial speed is one;
on this compact rectangle, positivity of rho and cosh(t) gives strict
hyperbolicity, positive auxiliary wave forms and a timelike shell. These
properties, and the strict speed bound `<2`, persist in a sufficiently small
state neighborhood. The exact boundary state vector field is smooth near
zero; its local flow can be cut off in a fixed collar of width
`0<h<ell/4`, further decreased to fit its local flow interval. This collar
lies inside the exact initial zero collar `z<ell/2`. Keep this fixed h and
state neighborhood for all parameters, and shrink the parameter interval
using the prescribed-slab theorem to retain them on all `[0,Tphys]`.

This is a qualitative small-parameter selection; no numerical lower bound
on the allowed parameter radius or on the state-flow interval is claimed.
The shell theorem's analytic hypotheses and construction are retained.
Smooth coefficient extension and any artificial outer construction have no
influence on the observation's domain of dependence under the strict causal
margin. The parameter-dependent cones need not coincide with background
cones. A common background-labelled `(t,z,Omega)` chart and fixed shell
`z=0` are used on this entire rectangle. In that chart `dy=-rho0 dz` uses
the fixed background function, not a parameter-dependent proper-distance
coordinate.

## 4. Complete the gauge, constraints and first solution

For this exact canonical curve choose initial lapse and shift

```
N_e=rho(1+e uY), beta_e=0.
```

Use the exact local relative-harmonic completion with reference connection
`Gamma[g0]`. Write `R^A=g_e^{BC}Gamma[g0]^A_BC`. The companion's formulas are

```
(partial_t-beta.d)N=N^2 tr_gamma K-N^3 R^0,
(partial_t-beta.d)beta^i
 =N^2(Gamma_gamma^i-R^i-beta^i R^0)-N gamma^{ij} partial_j N,
partial_t gamma=2NK+Lie_beta gamma,
partial_t phi=N Pi+beta.d phi.
```

They give `C=0` initially exactly and preserve the open zero shell collar.
The exact physical constraints are the projections of
`H=E-g tr(E)/2`. Given the reduced equations, they imply
`partial_t C=0` initially. Neither `E_nn=0` alone nor an arbitrary prescribed
second acceleration substitutes for this implication. The exact Bianchi
identity with `S=0` gives `(Box_g+Ric_g)C=0`. Codazzi and both physical shell
laws give

```
C(n_out)=0, (nabla_n_out C)_a-K_out,a^b C_b=0,
or C_y=0, partial_y C_a-N sigma C_a/3=0.
```

The factor `1/3` includes the connection contribution. The retained linear
mixed energy uniqueness therefore propagates `C=0`. The initial zero collar
supplies the finite corner identities by the local reduced recursion. On
identical copies the same construction and uniqueness give matching induced
fields and the copy-even branch. These are exact-family premises and
consequences of the retained analytic shell theorem, not inferences from
intervals containing zero or from a finite list of source identities.

At first order this initial completion gives the actual full harmonic
one-jet, including

```
dot q_ty=-6rho^2(u_y+4au)Y.
```

Indeed the registered gauge wave has `f=b=d=f_t=b_t=d_t=0` as initial
functions, but `f_tt=-6rho^2(u_y+4au)`, `b_tt=d_tt=0`. Its nonzero mixed
velocity must be retained. The canonical fields and lapse/shift values
agree, and the five linearized harmonic equations determine uniquely the
remaining five gauge velocities. This is the completion permitted by the
shell theorem; the older precontact evolution with background lapse for all
parameters is a different representative.

For first-field identification on the slab, retain the actual smooth master
and gauge solutions, including their exact equations as neighborhood
identities. The current natural master binding proves the longitudinal pair
solves `DE=DS=0`. The dynamic shell trace as well as the master equation is
needed for the physical scalar junction. The registered gauge vector obeys
`f|shell=0, b_r|shell=d_r|shell=0`; it is shell preserving, gives the four
mixed coordinate rows, preserves the physical linearized junctions, and
its PDE gives `DC[q]=0`. Hence the registered pair and the exact family's
first derivative solve the same linearized reduced bulk and physical/gauge
boundary system with the same full initial one-jet. Linear uniqueness in
the retained mixed energy class identifies them on their common causal
domain, in particular on the entire backward region needed at Pdagger.
This applies locally if the auxiliary construction uses a different outer
wall: neither wall reaches that backward region. Equality of initial metric
values alone would not establish this step.

## 5. Initial TT displacement and velocity of the exact weighted curve

Let `P_TT` denote angular projection on the closed round S3. In the actual
radial time-symmetric background,

```
P_g^*N=(Hess N-gamma0 Delta N-N Ric0)/2
       +N dphi0 tensor dphi0/2,
P_m^*Zeta=Lie_(Zeta sharp) gamma0, with raised indices as appropriate.
```

The lowered angular block of the first is a sphere Hessian plus angular
trace. The radial scalar-gradient term has no angular block, and the
background angular Ricci block is isotropic. The normalized momentum
adjoint angular block is a sphere symmetric gradient plus trace. For every
smooth TT tensor T,

```
integral T^{ab} D_aD_b N=0,
integral T^{ab} D_(a V_b)=0,
integral T^{ab} f kappa_ab=0.
```

These identities follow by integration by parts, using `div T=tr T=0`.
Radial weights commute with the projection; the identities extend from
smooth potentials to the finite Sobolev weighted closure. Since
`x''(0)=A u''(0)` lies in the fixed range, they prove

```
P_TT gamma2=0, P_TT(pi2/rho0)=0.
```

The physical K2 requires the quadratic canonical conversion, not only its
linear part. Let `U1=uY`, `S1=utY/rho`. Differentiating the full inverse
canonical map at `pi0=0`, in an orthonormal frame, gives

```
K2=2[pi2_lower/rho0-I tr(pi2_lower/rho0)/3]
   +diag(36 U1 S1,-6 U1 S1,-6 U1 S1,-6 U1 S1).
```

No gamma2 enters because it multiplies pi0. Therefore the angular quadratic
correction is pure trace and `P_TT K2=0`. The specified lapse and shift give

```
dot gamma2=2rho K2+4rho U1 K1,
P_TT dot gamma2=0.
```

The second term is again an angular trace. The fixed background frame
rescaling and multiplication by ell preserve the zero initial pair. At
t=0 the derivative of cosh(t) vanishes; more generally the derivative of
any radial background frame scale multiplies a zero TT displacement.

This construction directly selects the exact family in the registered
chart; no unconstructed second coordinate transformation is required. The
related pullback lemma remains valid: if X and X_t vanish as initial spatial
functions, the angular value and time derivative of
`2Lie_X q+Lie_X^2 g0` vanish even if X_tt is nonzero; `Lie_Y g0` has zero
angular TT for any smooth second generator. That lemma alone would not
identify moved hypersurfaces or the shell embedding, which is why the direct
initial completion is used here.

## 6. General angular projection and the exact second variation

The shell theorem and its finite-regularity addendum give
`e -> (g_e,phi_e)` of class C2 into X10 on the common fixed chart, with the
nonlinear fields in X14. In four spatial dimensions this controls a lower
piecewise spacetime C4 topology. Consequently parameter differentiation,
the finite number of spacetime derivatives below, angular integration and
physical-side shell traces commute at the needed orders. The second
variation is an actual derivative in this topology, not a formal solution
defined only by differentiated equations.

For clarity, finite closure of the registered 23-component calculation alone
would not prove decoupling from arbitrary angular second data. The extension
used here is the rotational invariance of the **linearized background**
bulk and boundary operators. Scalar-derived, vector-derived and TT tensor
harmonics form their usual orthogonal S3 decomposition. Scalar/tensor
cross terms contain traces, sphere Hessians or symmetric gradients and
vanish against T by the preceding integration-by-parts identities. The
rough sphere Laplacian is self-adjoint and preserves the TT eigenspaces.
The angular curvature endomorphism of this warped round background acts
algebraically on a TT tensor, and time/radial coefficient functions commute
with angular projection. Shell rows have the same isotropic background
structure. Thus each TT eigenspace is an autonomous linear channel for
arbitrary second data, with no requirement that the nonlinear solution
retain a finite angular expansion. Equivalently, the lowest tensor
representations `(2,0)` and `(0,2)` of SO(4) occur in the angular TT block
and not in scalar or vector slots.

More explicitly, write a pure projected angular perturbation as
`k_ab=Rphys^2 H(t,y) T_ab` in coordinate angular slots, where
`Rphys=rho cosh(t)` and `Delta_kappa T=-lambda_T T`. The TT projection of
the linearized Einstein–scalar equation on a spherical warped product is,
up to the common nonzero principal factor used in the forcing adapter,

```
Box_base H+3(dRphys/Rphys).dH-(lambda_T+2)H/Rphys^2,
base=-rho^2 dt^2+dy^2.
```

Here is the direct tensor calculation, retaining the mixed connection
terms. A pure angular TT k has zero full spacetime trace and divergence,
so the linear harmonic covector vanishes. In the curvature convention of
the registered source, the angular natural residual is

```
DE_ab=-Box_0 k_ab/2-Riem_aC bD k^CD.
```

Indeed the two Ricci endomorphism terms in the trace-free Ricci variation
combine to `(2V/3)k_ab` by the angular background Einstein equation and
cancel the linear potential term. The radial scalar-gradient variation has
no TT block. If `J=|dRphys|_base^2`, the warped connection gives, after
dividing angular coordinate components by Rphys^2,

```
(Box_0 k)_ab/Rphys^2
  =[Box_base H+3(dRphys/Rphys).dH-(lambda_T+2J)H/Rphys^2]T_ab,
Riem_aC bD k^CD/Rphys^2=-(1-J)H T_ab/Rphys^2.
```

The `-2J` in the rough wave comes from the angular connection acting on
each tensor slot through mixed base/angular components; setting those
covariant derivatives to zero because the mixed components themselves
vanish would lose it. Substitution cancels J and gives precisely
`DE_T=-[Box_base H+3(dRphys/Rphys).dH-(lambda_T+2)H/Rphys^2]/2`.
Thus the base wave contributes `H_yy+aH_y-rho^-2 H_tt`, the three angular
directions supply `3aH_y-3tanh(t)rho^-2 H_t`, and the total angular term
is `-(lambda_T+2)H/Rphys^2`. The `2` is the unit-sphere curvature
contribution; with `lambda_T=6`, the coefficient is eight. A reduced
gauge-gradient term has angular symmetric-gradient/trace form and zero TT
projection. In the fixed inward coordinate,
`rho^2 H_yy=H_zz+rho a H_z` and `rho^2 aH_y=-rho aH_z`. Thus the universal
homogeneous projected equation is exactly

```
H_tt-H_zz=-3tanh(t)H_t-3rho aH_z-8H/cosh(t)^2.
```

In the background orthonormal angular frame the exact principal
normalization is

```
DE_T=[H_tt-H_zz+3tanh(t)H_t+3rho aH_z+8H/cosh(t)^2]/(2rho^2).
```

Consequently `DE_T=-QE_T` gives exactly the registered scaled forcing
`-2ell^2 rho^2 QE_T`. The TT projection of the reduced gauge gradient is
zero, so this also equals the reduced source convention `QEred_T`.

For the shell, linearizing `K_out,ab=sigma h_ab/6` gives in TT slots
`(partial_y k_ab-2a k_ab)/2=0`. Since
`partial_y Rphys^2=2a Rphys^2`, this is `H_y=0`, hence `H_r=0`.
Scalar variations of the tension, lapse and normal multiply background
traces, and mixed-vector terms project to zero. The fixed first field adds
the retained quadratic inhomogeneous source. After division by Rphys^2,
the linear shell residual is `H_y/2`; hence `B'_T[k]+QKab_T=0` gives
`H_r=2ell rho QKab_T` with the positive inward r convention. These formulas identify the
projected full PDE and boundary law directly; the finite 23-component TT
rows then fix precisely the same amplitude and source normalization.

The registered `Tplus,Tminus` have `div=tr=0`, rough eigenvalue `-6`, curl
eigenvalues `+3,-3`, equal squared probability-sphere norms `1/6`, and zero
cross inner product. Define, in the fixed background angular frame,

```
U_exact,plus(x,r)=<k_ab,Tplus>/<Tplus,Tplus>,
U_exact,minus(x,r)=<k_ab,Tminus>/<Tminus,Tminus>.
```

No parameter-dependent TT projector or perturbed-volume inner product is
used in this definition. The exact second-variation equations contain the
linear operator applied to k and the complete quadratic forcing of the
already identified first pair. In the registered normalization each selected
coefficient therefore satisfies

```
U_xx-U_rr=A_T U_x+B_T U_r+D_T U+F_T,
U_r(x,0)=G_T(x), U(0,r)=U_x(0,r)=0,
A_T=-3ell tanh(t), B_T=-3ell rho a,
D_T=ell^2[-6/cosh(t)^2
  +rho^2(-2a^2+2(tanh(t)/rho)^2-2bcurv-8dcurv+4V/3)],
bcurv=-(a_y+a^2), dcurv=rho^-2-a^2.
```

Both background identities in section 1 reduce D_T to
`-8ell^2/cosh(t)^2`. The unsimplified coefficient is precisely the one used
in the frozen numerical enclosure; that enclosure need not be rerun with
a different interval dependency pattern. The TT shell Robin term is zero.
The retained exact source maps give the same two 78-term forcing rows and

```
F_T=-2ell^2 rho^2 (QEred)_T,
G_T=2ell rho QKab_T
   =-(4/3)ell rho [a A-(2p0/3) scalar1] M.
```

The shell source is generally nonzero; it cannot be dropped. Here A,M and
scalar1 are the registered actual first frame amplitudes. `QEred` uses the
full-metric-lowered harmonic convention with `C0=C1=0`; its factor and sign
are those of the second derivative, not its half-coefficient. These physical
bulk, shell and gauge-source identifications retain the audited all-point
source formulas and actual first-equation premises. They do not assume that
arbitrary interval selections are coherent on-shell jets.

Let U be the complete actual scalar formal response already enclosed by
the postcontact checkpoint. Its reflected integral problem has the same
F_T,G_T, operator and zero initial TT pair. On

```
Cdagger={0<=s<=5/4, max(0,s-9/8)<=v<=11/8-s}
```

the retained contraction bound is strictly below one. Subtracting the two
integral equations gives `X(U_exact-U)<=q X(U_exact-U)` with q<1, so

```
U_exact,plus=U_exact,minus=U on Cdagger,
and their x derivatives agree there, in particular at Pdagger=(5/4,1/8).
```

This uniqueness proof closes the TT coefficient match even though other
second initial sectors differ. It uses the reflected physical shell
history through `x-r=9/8`, not a precontact cone.

## 7. Conditional theorem and registered screen consequence

**Theorem.** Accept the retained smooth actual background/master/gauge
premises and their physical boundary/corner conditions; the exact localized
weighted inverse theorem on the explicit annulus of section 3; the finite
mixed-boundary shell evolution theorem, subsidiary uniqueness and its
finite-regularity parameter addendum in the selected small state
neighborhood; and the retained physical source and registered screen
identities. Select the exact weighted initial curve and initial harmonic
completion specified above. Then for sufficiently small |e| this curve has
an exact copy-even Einstein–scalar evolution through the observation on a
common controlled chart. Its full first derivative is the registered actual
harmonic first pair, and both selected angular TT components of its second
derivative equal the published scalar U on Cdagger. In particular the same
published interval encloses the registered screen second derivative for
this exact family.

**Proof.** Sections 1–2 identify the model, branch and all four admissible
first seed sectors. Section 3 applies the retained initial and shell
theorems with a strict support/causal margin. Section 4 supplies the full
gauge match and propagates the exact constraints, then uses linear
bulk/boundary uniqueness to identify the first fields. Section 5 proves
the exact family's zero angular TT initial pair. Section 6 justifies
actual differentiation, unrestricted angular projection, the same forced
TT equations and their unique solution. The registered kinematic screen
identity then applies to this actual differentiable family. No equality of
the remaining second components is used. QED.

For that identity, at every e use the future normal to the fixed t slices,
the positive unit spatial normal to the fixed y spheres, and the screen
tangent to their intersection, with the positive symmetric inverse-square-
root frame and the registered sphere orientation. This rule varies the
observer and frame algebraically with the metric exactly as the original
screen computation prescribes; it is not a rule for geometrically
transporting clocks or worldlines.

With `Rphys=rho cosh(t)` and the original complete first-field quadratic
map Fquad,

```
Pi_pair H''(0)=h(Tplus-Tminus),
h=3U_x/(2ell rho^2 cosh(t))+4Fquad/Rphys.
```

The map includes the initial first field's lapse, shifts, connection,
screen-volume and observer variations. It applies because the first pair
is identical and the linear screen contribution from any other second
scalar/vector sector has zero selected TT projection. Therefore the prior
primary outward rational interval is inherited without a new numerical
evaluation:

```
[-1049651791613786574257886764672101/2596148429267413814265248164610048,
 -18954240643899648600527883239121/1298074214633706907132624082305024].
```

It lies strictly below zero. This is a conditional nonzero second
coefficient for the stated exact family and registered prescription.

**Qualitative small-parameter corollary.** Apply the registered screen rule
at the same coordinate event at every parameter, and set

```
J(e)=3 integral_S3 <H_e(Pdagger),Tplus-Tminus> dmu_kappa,
```

where the probability-sphere measure and tensor components use the fixed
registered frame identification. The factor three uses
`||Tplus-Tminus||^2=1/3`. Thus `J''(0)=h`. The background screen tensor
vanishes, and its complete first variation vanishes for the actual scalar
first pair, so `J(0)=J'(0)=0`. The lower C2 parameter topology and smooth
screen construction give J of class C2. Let h_hi denote the strict negative
upper rational endpoint above. Continuity gives an unspecified delta>0
within the common parameter interval such that `J''(e)<=h_hi/2` for
`|e|<delta`. The exact integral Taylor formula, valid for either parameter
sign,

```
J(e)=e^2 integral_0^1 (1-s) J''(s e) ds <= h_hi e^2/4 < 0
       for 0<|e|<delta,
```

proves a qualitative small-amplitude sign for this selected family and
registered angular coefficient. No effective delta or bound at a specified
nonzero parameter follows. This elementary C2 consequence adds no claim
about a transported observer, another nonlinear completion or a detector.

### The selected nonlinear completion is essential

The same first tangent alone does not determine this TT history. To make
the limitation concrete, let `C=Tplus+Tminus`. The retained all-point angular
identity `Pi_TT(DY tensor DY)=-(2/3)C` shows that C is nonzero, with squared
probability-sphere norm `1/3`, and is invariant under the whole selected
reflection group. Indeed Y and its polar tensor square are invariant, and
the round-sphere TT projector commutes with every isometry. Individual
chiral templates need not themselves be reflection invariant.

For a nonzero smooth radial f compactly supported inside A', set
`B_ab=rho^2 f C_ab`, `B_yy=B_ya=0`, with every other canonical sector zero.
This spatial metric direction has zero full spatial trace and divergence:
the radial divergence is minus a times its angular trace, and the angular
divergence is the sphere divergence of C. Its Ricci contraction and radial
scalar-gradient contraction are zero. Hence the linear Hamiltonian
variation vanishes, while the momentum variation has no metric term at
the time-symmetric background. Thus `D Chat[B]=0`.

Replace the weighted IFT trial path by
`x0+e v_seed+e^2 B/2+A u_tilde(e)`. The same unknown derivative is the
same weighted isomorphism L, and the mapping/support arguments are unchanged.
The IFT therefore gives another exact constraint curve with the same first
tangent and `u_tilde''(0)=u''(0)`, since `D Chat[B]=0`. Its angular TT
second metric displacement differs by `rho^2 f C`. The same local completion
and prescribed-interval estimates apply after parameter shrinkage. Its TT
evolution differs by a homogeneous solution with this nonzero initial datum.
This does not assert a different numerical endpoint for an unspecified f;
it establishes that common first tangent and reflection symmetry do not
force a common TT initial pair or history. The theorem above accordingly
selects the fixed weighted-range curve, or permits another curve only after
its completed TT initial pair is separately matched.

### Primary context for the universal tensor operator

The generic warped-product tensor equation is established literature.
Kodama, Ishibashi and Seto, *Brane World Cosmology — Gauge-Invariant
Formalism for Perturbation*, Physical Review D **62**, 064022 (2000),
[arXiv:hep-th/0004160v3](https://arxiv.org/pdf/hep-th/0004160), section II.C.1,
equations (19)–(23), gives the tensor equation for an arbitrary Lorentzian
base and constant-curvature fiber with angular-isotropic background matter.
Its metric convention is `h_ij=2R^2 H_T T_ij`; the present binding is
`H_T=H/2`, fiber dimension 3, curvature 1 and rough eigenvalue 6. Its
anisotropic-stress source uses a mixed angular index. The canonical scalar
background has no linear tensor anisotropic stress; in covariant indices
the pressure times metric-variation term cancels when raising an index.

Cai and Cao, *Generalized Formalism in Gauge-Invariant Gravitational
Perturbations*, Physical Review D **88**, 084047 (2013),
[arXiv:1306.4927](https://arxiv.org/pdf/1306.4927), section III.B.1 equation
(3.10), gives the tensor equation before spectral expansion, and section
IV.A equations (4.1)–(4.3) explicitly display that covariant pressure
correction. These primary results corroborate the generic identity and its
normalization; neither is being used as a nonlinear shell-evolution theorem.
The present advance is the source-bound application and selected exact-family
coefficient match, not a new universal tensor perturbation equation.

## 8. Premise and verification ledger

The new proof closes the finite matching joins: actual model and branch,
all four canonical seed sectors and profile normalization, new weighted
curve selection, quadratic canonical-to-physical TT conversion, permitted
initial lapse/shift completion, causal arithmetic, first-field uniqueness
premises, angular TT projection outside a finite truncation, and scalar TT
uniqueness. The analytic existence, local adjoint/API/weighted estimates,
mixed-boundary energy, actual master/gauge regularity and parameter
theorems remain explicitly retained mathematical inputs. Source
authentication and executable identities do not establish those theorems.

This note does not claim an equality with every nonlinear completion of the
same first tangent, a transported-observer invariant, an effective epsilon
radius, a quantitative remainder or sign at a specified finite amplitude, detector or
radiation energy response, stability, or a cosmological-origin result.
The full recovered 51-member quasilinear theorem capsule is authenticated by
its externally pinned manifest, alongside the selected formula sources.
These checks retain the earlier published actual-field/numerical premise
chain; they do not constitute external specialist review.

`SOURCE_PINS.json` records the exact inspected source bytes. The local
verification script checks specified finite algebra and domain identities
and rejects actual formula mutations in fresh normal and optimized
processes. It does not import a numerical producer, evolve a PDE or certify
the analytic argument. Detailed execution receipts are separate files.

The source paths in that manifest are relative to a complete checkout's
`research` directory. The principal claim-to-source mapping is:

| Claim or input | Exact retained source |
|---|---|
| Generic weighted curve, fixed range, support and second derivative | `HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md`, sections “Conditional inputs”, “Fixed weighted spaces”, “Physical zero extension”, “Exact fixed-regularity conclusion”, “Concrete commuting quadrupole consequence” |
| Adjoint formulas, actual canonical normalization, inherited ambient seed | `HDBLAST_SECOND_ORDER_LOCALIZED_20261014/manuscript/SECOND_ORDER_LOCALIZED.md`, sections 1–3 |
| Prescribed finite interval and common-parameter theorem | `HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/EXACT_FINITE_REGULARITY_SHELL_EVOLUTION.md`, sections 1–5 |
| Exact gauge completion, subsidiary system and conservative parameter topology | Same capsule, `gauge-parameter/INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md` and `FINITE_REGULARITY_VARIATION_COMPLETION.md` |
| Zero-fixing boundary flow and prescribed-interval analytic estimates | Same capsule, `inputs/NONLINEAR_BOUNDARY_STRAIGHTENING.md`, sections 1–2, and `astra/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md`, sections 1–7 |
| Actual packet and all physical first data | `HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_COMPACT_SECOND_INITIAL_JET.md`, section 1 |
| Exact initial lapse/shift and first mixed velocity | Same capsule, `inputs/EXPLICIT_INITIAL_HARMONIC_COMPLETION.md`, sections 1–4, and `first-fields/FIRST_HARMONIC_FIELDS_UNPUBLISHED.md`, sections 1–3 |
| Current initial constraint and complete ADM implementation joins | `HDBLAST_REGISTERED_INITIAL_CONSTRAINTS_20261030/theory/PHYSICAL_INITIAL_BINDING_NOTE.md` and `REGISTERED_ADM_BINDING_NOTE.md` |
| Current first natural equation, dynamic shell trace and physical boundary | `HDBLAST_NATURAL_FIELD_BINDING_20261030/theory/NATURAL_FIRST_MASTER_BINDING.md` and `FIRST_SHELL_MASTER_BINDING_NOTE.md` |
| Current harmonic PDE and complete first one-jet | `HDBLAST_BULK_SHELL_CONSTRAINT_BINDING_20261030/theory/FIRST_WAVE_BINDING_NOTE.md` |
| Actual/reduced source roles and analytic subsidiary premises | `HDBLAST_REDUCED_SOURCE_BINDING_20261031/theory/CONSTRAINT_PROPAGATION_CHAIN.md` |
| Current common actual first fields and full domain | `HDBLAST_GEOMETRIC_BULK_SOURCES_20261016/numeric/ACTUAL_FIRST_JET_AND_SOURCE_PROOF.md` |
| Actual wider branch collar | `HDBLAST_INTERIOR_MASTER_JETS_20261015/potential/inputs/BACKGROUND_DERIVATIVE_AUDIT.md`, section 1, and `HDBLAST_SECOND_INITIAL_NORMS_GAUGE_20261015/background/inputs/ACTUAL_POTENTIAL_JET_BOUNDS.json` |
| Frozen scalar equation, actual forcing, reflected uniqueness and numerical result | `HDBLAST_POSTCONTACT_SCREEN_CURVATURE_20261031/DERIVATION.md`, `main/REFLECTED_TT_PROOF.md`, `main/POSTCONTACT_NORMAL.json` |
| Original screen, first-field quadratic map and orientation | `HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/theory/DERIVATION.md` |
| Reflection-invariant nonzero TT template for completion dependence | Same capsule, `review/symmetry/SYMMETRY_AUDIT.md` |

With SymPy 1.14.0, the local commands are

```
python -I -B check_binding.py --repo-root CHECKOUT --run-controls --output NORMAL.json
python -I -B -O check_binding.py --repo-root CHECKOUT --run-controls --output OPTIMIZED.json
```

The report contains 40 finite identities/inequalities or authenticated
report joins, ten rejected formula mutations, and two rejected copied-source
byte mutations. One source control changes a selected input; the other
changes a member of the full analytic capsule outside the selected-source
list. These authentication rejections are distinct from algebraic formula
mutation rejections. Both modes must produce byte-identical reports. The
inherited numerical endpoint and contraction arithmetic checks read the
authenticated old report; they are not a new evaluation of its source
cells. The local universal-TT coordinate conversion check does not replace
the tensor derivation in section 6.
