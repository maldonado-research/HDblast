# An anchored wave-map realization of the registered screen coefficient

31 October 2026. This corollary applies to the conditional exact family selected
in [the exact-family derivation](DERIVATION.md), retaining all of its analytic,
first-field and numerical premises. It specifies an additional reference
structure.

## 1. Result and its precise domain

An independently anchored physical-to-background wave map can give the
registered coefficient exactly. On a conditional exact family already known
to satisfy the registered relative harmonic gauge and mixed coordinate shell
rows, existence is witnessed by the identity map. A scalar-principal
difference-energy argument gives uniqueness on its controlled domain of
dependence. Consequently the construction extends naturally to the family's
boundary-preserving diffeomorphism orbit. It is not necessary to assume a
standalone wave-map existence theorem for an arbitrary nearby metric to make
this restricted statement.

The result is reference dependent. It depends on the target metric, initial
anchor map and velocity, angular templates, orientation and measure. Changing
those data can change the observer and coefficient. Diffeomorphism covariance
of the complete anchored object does not imply independence of these choices.

The [exact-family binding](DERIVATION.md) establishes the common first tangent
and second TT initial data conditional on its stated hypotheses. This lemma
retains those hypotheses and the numerical error charges. Its conclusion is
that the second derivative of the functional in Section 5 is exactly the
previously registered `h`.

## 2. Independently specified reference data and the correct wave-map sign

Fix the target manifold with boundary

```
M0 = [0,T] x [0,Z] x S3,
g0 = rho(Y)^2 [dZ^2-dT^2+cosh(T)^2 kappa], dY=-rho dZ,
B0={Z=0}.
```

Here `Y` is the old outward radial coordinate, while `Z` is inward. This
notation uses `Y` only for the target radius; below `U_init` denotes the full
registered initial lapse perturbation, including its angular harmonic. The
target also carries its time orientation, the registered outward ambient
sphere orientation, the round metric `kappa`, its probability measure
`dmu_kappa/Vol(S3)`, and the fixed TT templates `Tplus,Tminus` with squared
normalized L2 sphere norms `1/6` and normalized cross inner product zero.

Before evolving, prescribe an initial physical hypersurface `Sigma`, its
future unit normal `u`, an anchor map `f:Sigma -> {T=0}`, and a target tangent
vector `V` along `f`. Seek a map `Phi:(M,g)->(M0,g0)` with

```
Phi|Sigma=f,                 dPhi(u)|Sigma=V,
tau(Phi)^A = g^{mu nu}(partial_mu partial_nu Phi^A
                -Gamma[g]^lambda_mu nu partial_lambda Phi^A
                +Gamma[g0]^A_BC(Phi) partial_mu Phi^B partial_nu Phi^C)=0.
```

The plus sign on the target connection is essential. In the registered
common chart,

```
tau(id)^A = -g^{mu nu}(Gamma[g]^A_mu nu-Gamma[g0]^A_mu nu) = -C^A.
```

Thus the retained condition `C=0` makes the identity a wave map in the
proposed physical-to-target direction. Reversing the map direction is a
different equation.

The initial reference rule must match the completed lapse and shift. Write
the physical initial ADM normal as

```
u_e=N_e^{-1}(partial_t-beta_e^i partial_i).
```

On the fixed initial chart take `f_e(t=0,z,Omega)=(0,z,Omega)` and prescribe
the target velocity by the *initial* rule

```
V_e=N_e^{-1}(partial_T-beta_e^i partial_{X^i}) along f_e.
```

For the retained completion `N_e=rho(1+e U_init)`, `beta_e=0`, this gives

```
V_e=[rho(1+e U_init)]^{-1} partial_T,
V_0=rho^{-1} partial_T,
V'_0=-U_init/rho partial_T,
V''_0=2 U_init^2/rho partial_T.
```

The rule is prescribed from the initial packet and chosen completion, before
the solution is known. It is a fixed rule for the family, although its values
vary with `e`. Fixing `V_e=rho^{-1}partial_T` for every parameter would fail
already at first order. Neither the lapse nor the reference velocity is
determined by the canonical physical data alone; the chosen reference data
must remain part of the definition. Under a physical diffeomorphism the
anchor map and velocity along it transform with the initial hypersurface.

## 3. Boundary law, uniqueness, and the limited existence claim

Use the geometric boundary law

```
Phi(B) subset B0,
dPhi(n_out) is g0-orthogonal to T B0.
```

In the fixed target product chart this is exactly

```
Phi^Z=0,                 n_out(Phi^a)=0, a=T,Omega.
```

It is invariant under physical diffeomorphisms. A target coordinate change
must transform the target metric and boundary projector as well; arbitrary
untransformed component conditions are not the geometric definition.

For `Phi=id`, the second boundary condition is

```
n_out^a = +/- g^{aZ}/sqrt(g^{ZZ})=0.
```

This is equivalent to the retained mixed coordinate rows `g_Za=0` on the
timelike shell. The Israel and scalar shell laws do not themselves imply
these coordinate rows. The condition `C(n)=0` is a separate primal reduced
Einstein boundary row; the subsidiary argument proves `C=0` throughout and
therefore supplies the wave-map bulk identity. No additional incoming
Einstein boundary conditions are imposed here.

The initial data equal the background on an open shell collar, including
the lapse and shift. Hence the anchor data there agree with the background
identity data, supplying the compatible corner wedge on the existing
family. This is stronger than checking only the zeroth boundary trace.

For completeness the uniqueness argument needed here is short. Fix one
physical metric satisfying the declared hyperbolic and timelike-shell
margins. If `Phi` and `Psi` are two sufficiently regular maps with the same
anchors and boundary data, write their difference in target product charts.
On a compact region where their first derivatives are bounded and the target
coefficients are smooth, the difference satisfies

```
Box_g w^A = A^{A mu}_B partial_mu w^B + B^A_B w^B,
w^Z|B=0,                    n_out w^a|B=0,
w|Sigma=0,                  u(w)|Sigma=0.
```

The coefficient matrices are bounded by the mean-value formula for the
semilinear target-connection term. Choose a smooth future timelike energy
multiplier tangent to the physical shell. In the registered chart the slice
normal has this property because `g^{Za}=0`. For each component the boundary
flux is a constant multiple of `(multiplier w)(n_out w)`: it vanishes on
Dirichlet components since the multiplier is tangential, and on Neumann
components by the normal condition. Add `|w|^2` to the positive auxiliary
energy and apply Cauchy-Schwarz to the lower-order terms, obtaining

```
E(t) <= C integral_0^t E(s) ds.
```

Gronwall gives `w=0`. A finite atlas and local uniqueness patch this argument
on the sphere; alternatively use the standard embedding `S3 subset R4`,
whose angular components all have homogeneous normal derivative and a smooth
semilinear wave-map equation in an ambient neighborhood. Equality initially
and continuation keep the two maps in common target neighborhoods. At an
artificial outer boundary use the same outer data or restrict to the
backward-closed region protected by the exact family's strict *perturbed*
speed margin. No unspecified outer boundary is silently omitted.

This uses only bounded coefficients and sufficient classical regularity
(piecewise spacetime `C4` is more than enough for this uniqueness argument).
It is not a finite-`C4` existence theorem. For the conditional exact harmonic
family selected in [the derivation](DERIVATION.md), the identity is an explicit global candidate satisfying the
independently prescribed data, so existence and invertibility on that family
are immediate. On its diffeomorphism orbit the candidate is the corresponding
coordinate map. This proves uniqueness of the independently specified
reference on the restricted class; it does not define the reference by
declaring it equal to the answer after solving.

For an arbitrary metric outside that class, existence on the whole slab,
corner regularity, parameter dependence, and invertibility remain additional
gates. A nonzero Jacobian gives local invertibility only; global event
uniqueness requires an injectivity/covering argument or an explicit common
inverse domain.

## 4. The full observer, radial normal, screen, orientation and event

Let `tau=T o Phi`, `zeta=Z o Phi`. Require `d tau` timelike and its level
sets future oriented. Define

```
L=sqrt(-g^{-1}(d tau,d tau)),
U=-g^{-1}(d tau)/L,
P^{mu nu}=g^{mu nu}+U^mu U^nu,
v=P(d zeta),
n=-v/sqrt(g(v,v)),
S^{mu nu}=P^{mu nu}-n^mu n^nu.
```

The minus sign makes `n` outward: `Z` increases inward, whereas the old
positive `y` normal increases outward. The radial normal is projected within
the time slice; replacing it by the full spacetime gradient of `zeta` would
change the observer construction when a time shift is present. The screen
is `ker(d tau) intersect ker(d zeta)`, with induced positive metric `s`.
Transport its orientation by the target sphere map, retaining the registered
outward ambient sphere orientation.

For each target angular label `omega`, the event is

```
p_e(omega)=Phi_e^{-1}(Tdagger,Zdagger,omega),
Tdagger=5 ell/4, Zdagger=ell/8, ell=10^-6.
```

Thus the old `(x,r)=(5/4,1/8)` denotes a labeled observation sphere for the
angular coefficient, rather than a single unspecified angular event. For
`Phi=id`, these events and all the displayed observer/screen formulas agree
exactly with the registered prescription. On the shell `n_out(tau)=0`, so
`U` is tangent to the shell. At interior events the nondegeneracy conditions
are open and hold for a sufficiently small common parameter interval by
the exact family's `C4` continuity and the background margins.

This time is a wave-map label. `U(tau)=L`, so proper time along this observer
satisfies `d proper_time/d tau=1/L`. This lapse need not be one, and the
observer need not be geodesic. At the background shell it is `rho_b`, as in
the retained source. This construction supplies no new proper-time clock or
geodesic transport assertion.

## 5. Tensor, frame prescription and complete angular functional

Use the retained curvature sign and define the screen tensor, including its
observer and its induced epsilon tensor, by

```
H_AB=(1/2) epsilon(s)_(A^{CD} Weyl(g)_{B) U C D},
```

where parentheses are normalized symmetrization. The source derives the
equivalent ADM Codazzi expression and includes its complete lapse, shift,
radial, connection, epsilon and frame variations.

Here is a frame-free definition of the source's positive symmetric inverse
square root prescription. Write `Psi=Phi^{-1}` and restrict it at each fixed
target `T,Z` to the target sphere. Put

```
R0(T,Z)=rho(Y(Z)) cosh(T),
sbar=Psi^*g restricted to T(S3),
A=kappa^{-1}(sbar/R0^2),
B=A^{-1/2},                 J=(1/R0) dPsi o B.
```

The endomorphism `A` is positive and `kappa`-self-adjoint. Its unique positive
square root makes `J` an orientation-preserving isometry from the round
target tangent space to the physical screen. Define the target tensor

```
Hhat(v,w)=H(Jv,Jw).
```

This specifies every axis without choosing an additional rotation. In the
registered gauge it is exactly the source's inverse-square-root screen
frame relative to the background frame. Merely pulling back a covariant
two-tensor without this normalization would not reproduce that finite-
parameter convention.

Define the fixed target projection and scalar functional by

```
a_plus(e)=6 integral_S3 <Hhat_e,Tplus>_kappa dmu_kappa/Vol(S3),
a_minus(e)=6 integral_S3 <Hhat_e,Tminus>_kappa dmu_kappa/Vol(S3),
mathfrak_h(e)=(a_plus(e)-a_minus(e))/2,
Pi_pair Hhat=a_plus Tplus+a_minus Tminus.
```

Equivalently the physical templates and integration measure are transported
from these target data by the same sphere map. The measure is the chosen
reference probability measure, not an unspecified physical-area measure.
All angular fields, normalizations and the orientation are therefore part
of the functional. Finite-parameter `Hhat` need not lie in this two-mode
space.

On the existing common-gauge family, `Phi_e=id` for every `e`. Therefore
`Hhat_e` is exactly the previously defined registered tensor for every
parameter; there are no additional wave-map observer corrections to append
to its existing variation formula. Using the conditional exact-TT identification
in [the derivation](DERIVATION.md), the source gives

```
Pi_pair Hhat''(0)=h(Tplus-Tminus),
mathfrak_h''(0)=h,
h=3 U_Tx/(2 ell rho^2 cosh(t)) + 4 F_quad/(rho cosh(t)).
```

The derivative convention is `g_e=g0+e q+e^2 k/2+o(e^2)`; thus
`mathfrak_h(e)=e^2 h/2+o(e^2)` since the full defined `H0,H1` vanish.
There is no extra one-half in `mathfrak_h''(0)`. The previously enclosed
negative interval applies to this second derivative conditional on all the
same exact-family, first-field and numerical premises. It is not a new
numerical calculation.

There is also a qualitative finite-parameter consequence. Let
`u=-0.014601815851`, the outward upper endpoint of the retained decimal
enclosure, so `mathfrak_h''(0)=h<=u<0`. The full defined background and
first tensor variations vanish, hence `mathfrak_h(0)=mathfrak_h'(0)=0`.
The `C2` regularity established below gives a possibly smaller parameter
interval `|e|<delta` with `mathfrak_h''(e)<=u/2`. For either sign of `e`,

```
mathfrak_h(e)=e^2 integral_0^1 (1-s) mathfrak_h''(s e) ds
             <= (u/4)e^2 = -0.00365045396275 e^2 < 0,
             0<|e|<delta.
```

Here `delta>0` is existential and ineffective; it also lies inside the
common exact-family interval and all domain and nondegeneracy margins.
No specified nonzero parameter window or computable Taylor remainder
constant is supplied. This elementary consequence is a sign statement for
the selected angular projection of the reference-dependent tensor. It is
not a pointwise assertion at every angular event or a detector prediction.

## 6. Exact naturality, including events and templates

Let `F: Mtilde -> M` be an orientation- and time-orientation-preserving
diffeomorphism carrying the initial hypersurface and shell to their original
ones. Pull back the metric and scalar, and set

```
gtilde=F^*g,  Sigmatilde=F^{-1}(Sigma), Btilde=F^{-1}(B),
ftilde=f o F|Sigmatilde,     Vtilde=V o F|Sigmatilde.
```

The vector `V` takes values in the *target* tangent bundle, so this last
formula is composition along the transformed anchor, not a source-vector
pushforward. Wave-map and boundary covariance show that `Phi o F` solves
the transformed anchored problem. Uniqueness gives

```
Phitilde=Phi o F,           Psitilde=F^{-1} o Psi.
```

It follows successively that the clock and radius pull back, the observer
and radial normal transform by `dF^{-1}`, the screen and its oriented metric
transform naturally, `Htilde=F^*H` on that screen, the event becomes
`ptilde=F^{-1}(p)`, and `Jtilde=dF^{-1} o J`. Consequently

```
Hhat[gtilde,ftilde,Vtilde](T,Z,omega)=Hhat[g,f,V](T,Z,omega),
mathfrak_h[gtilde,ftilde,Vtilde]=mathfrak_h[g,f,V].
```

This is an exact identity of the whole target functional. It includes the
event, frame, angular fields and measure. It is stronger than checking an
isolated selected coordinate component under a first-order gauge change.
It applies equally to parameter-dependent diffeomorphisms when the relevant
parameter derivatives exist. Changing the reference data while keeping the
physical metric fixed is not the operation proved invariant here.

## 7. Finite regularity and parameter derivative accounting

The retained exact theorem provides a family `C2` in the lower `X10`
topology, hence `C2` into a lower piecewise spacetime `C4` norm. In the
registered representative `Phi_e=id`, inverse maps and target events are
parameter constant. Curvature loses two spacetime derivatives, so the
fixed-target curvature functional and sphere integral are `C2` in `e`.
The positive square root and all normalizations are smooth algebraic maps
on the uniformly nondegenerate neighborhood. Thus no unchanged-top-norm
flow differentiability or spatially smooth common family is required.

For a proposed general moving-reference theorem, one must account for the
loss explicitly. A sufficient *conditional* input for the scalar functional
is `g_e C2 into C4` and `Phi_e C2 into C3`, on a compact common inverse domain
with uniform nondegeneracy. At a physical event the functional depends
pointwise on the two-jet of `g` and one-jet of `Phi`. Two parameter derivatives
of evaluation at `p_e=Phi_e^{-1}(p0)` use spatial derivatives through order
four of `g` and order three of `Phi`, including the mixed derivatives.
The inverse event's derivatives follow from

```
p' = -(D Phi)^{-1} Phi',
p''=-(D Phi)^{-1}[Phi''+2 D Phi' p'+D^2 Phi[p',p']].
```

Compactness makes these formulas uniform for the angular integration. This
direct jet-level argument avoids assuming that inverse-map composition is
`C2` in an unchanged high spatial norm. Writing a pulled-back metric and
asserting `C2 into C4` for it without a loss calculation would be unjustified.
The present sources have not proved that an arbitrary metric's wave-map
solution enjoys these stated parameter spaces. For the restricted harmonic
family they are automatic because the map is constant.

Physical shell traces remain one-sided. The selected event has `Z>0`; no
distributional curvature on the shell is being evaluated. The construction
stops at the common-domain endpoint or loss of its time/radial margins or
inverse domain. On the known identity/diffeomorphism orbit there is no
additional caustic theorem to prove; on a larger class these are genuine gates.

## 8. Primary comparison and remaining claims

The primary text of Bruni, Matarrese, Mollerach and Sonego,
*Perturbations of spacetime: gauge transformations and gauge invariance at
second order and beyond*, gr-qc/9609040v2 (22 May 1997), was reread at Section
4.1 Proposition 1 and Section 4.2 Proposition 2, equations (4.8)-(4.9):

```
delta T_new-delta T_old = Lie_xi1 T0,
delta2 T_new-delta2 T_old = (Lie_xi2+Lie_xi1^2)T0+2 Lie_xi1 delta T_old.
```

These formulas concern a well-defined tensor with the required pullback
law. They do not supply a wave-map boundary existence theorem and do not
equate different observer prescriptions. Here exact naturality was proved
for the full reference-dependent functional before using its vanishing
background and first variation. No shortcut from a selected coordinate
component is needed. Bibliographic source:
<https://arxiv.org/abs/gr-qc/9609040v2>, DOI
<https://doi.org/10.1088/0264-9381/14/9/014>. No fresh external theorem or
unbounded literature search is invoked. Raw third-party text is not a runtime
dependency or part of this capsule.

The result is a diffeomorphism-covariant realization of the registered
coefficient relative to a specified anchored background wave map on the
conditional exact family's orbit. The time is a wave-map label, and the
observer is its unit-normal congruence. A proper-time/geodesic prescription
requires its own comparison and second-order terms. No observer independence,
detector signal, radiation theorem, or cosmological-origin conclusion follows.

## 9. Verification and sources

`check_reference_identities.py` performs bounded exact algebra checks of the
reference formulas, with named wrong-formula rejection controls. Its source
receipt authenticates the retained repository files before symbolic imports.
These identities and finite fixtures test formula implementation; they do not
prove nonlinear evolution, wave-map existence outside the stated orbit,
regularity estimates, or the retained numerical enclosure. Normal and optimized
runs are required to produce identical scientific JSON reports.

The exact-family selection and TT match are in [DERIVATION.md](DERIVATION.md).
Retained public sources are:

- [Postcontact coefficient and full observation convention](../HDBLAST_POSTCONTACT_SCREEN_CURVATURE_20261031/DERIVATION.md).
- [Prior next-step statement distinguishing the proper-time/tubular prescription](../HDBLAST_POSTCONTACT_SCREEN_CURVATURE_20261031/NEXT_STEP.md).
- [Complete screen tensor, frame, quadratic variation and TT norms](../HDBLAST_SCREEN_MAGNETIC_CURVATURE_20261026/theory/DERIVATION.md).
- [Finite-regularity exact shell theorem](../HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/EXACT_FINITE_REGULARITY_SHELL_EVOLUTION.md).
- [Harmonic completion, subsidiary constraints and parameter budget](../HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/gauge-parameter/INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md).
- [Finite-regularity variation completion](../HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/gauge-parameter/FINITE_REGULARITY_VARIATION_COMPLETION.md).
- [First harmonic gauge operator and inward/scaled coordinates](../HDBLAST_FIRST_GAUGE_WAVE_20261016/operator/GAUGE_OPERATOR_AND_CONSTRUCTION.md).
- [Explicit initial lapse and shift completion](../HDBLAST_COMMON_GAUGE_FORCING_20261016/inputs/EXPLICIT_INITIAL_HARMONIC_COMPLETION.md), Section 1 equation (2) and Section 2.
- [Scalar-principal mixed-boundary energy audit](../HDBLAST_QUASILINEAR_SHELL_EVOLUTION_20261015/sol/linear_audit/LINEAR_VARIABLE_WAVE_AUDIT.md), hypotheses and basic energy. Only its low-order energy mechanism is used for map uniqueness here.
- [Bounded prior primary-source comparison](../HDBLAST_POSTCONTACT_SCREEN_CURVATURE_20261031/literature/METHODS_SCOPE.md), including Bruni et al. Section 4.1 and Section 4.2.

Source hashes and bounded primary-body inspection metadata are recorded in
`PUBLIC_SOURCE_RECEIPT.json`. No raw primary text is needed by the checker.
