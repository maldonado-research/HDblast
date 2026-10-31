# Domain and common-seed gate for the exact TT comparison

31 October 2026. Supplemental source audit. The result below applies the retained conditional analytic theorems;
the arithmetic and symbolic checks do not computer-certify those theorems.

**Verdict:** the proposed annulus and prescribed slab are admissible for a
newly selected weighted exact initial curve with the actual registered first
seed. Its complete harmonic initial data agree with the actual common-gauge
first solution. Accepting the retained weighted-inverse and prescribed-time
mixed-evolution theorems, qualitative parameter shrinkage suffices. There is
no additional effective epsilon or shell-arrival-transversality requirement
for this coordinate coefficient comparison. This does not identify the old
unspecified exact curve, or a transported observer.

## Coordinates, background, and the target

Distinguish three radial variables and three time variables. In the current
certificate, `x=t/ell` and `r=z/ell`, with `ell=1/1000000`. Here `t` is the
unscaled de Sitter coordinate time, not unit-normal proper time. The positive
conformal radius is denoted `R` in this note:

```
R=b exp(-z),  s=R^2,  rho=Omega(s) R,
dy=Omega(s)dR=-rho dz,
ambient_T=R sinh(t),  ambient_R=R cosh(t).
```

The common initial section `t=0` is exactly `ambient_T=0`. On it
`partial_t=R partial_ambient_T`, while `partial_y=Omega^-1 partial_R`.
In scaled variables `partial_t=ell^-1 partial_x` and
`partial_y=-(ell rho)^-1 partial_r`. The label `T_phys` in the earlier
domain note therefore means the unscaled `t` interval; it is not a proper
clock duration. The background future unit normal is `rho^-1 partial_t`.
The metric conversion and ambient coordinate equations are source-bound by
`COVARIANT_FLUX_AND_PHYSICAL_DOMAINS.md:258-269` and
`FIRST_HARMONIC_FIELDS_UNPUBLISHED.md:89-97,150-163`.

The target is the current registered coordinate event
`(x,r)=(5/4,1/8)`, hence `(t,z)=(5ell/4,ell/8)`.
`DERIVATION.md:189-213` gives its complete reflected linear backward set,
initial footprint `[0,11ell/8]`, and shell history through `9ell/8`.
The current actual first fields, including the dynamic shell lift and
reflected branch, are supplied on

```
D = {0<=x<=3, r>=0, x+r<=8}.
```

This is the global actual-field certificate, not the precontact exterior
restriction `r>=x`; see `ACTUAL_FIRST_JET_AND_SOURCE_PROOF.md:21-43,110-128`.

The actual-background audit proves the strict first-exit continuation on
`0<=z<=8ell`; its potential-jet report explicitly enlarges to the rectangle
`0<=t<=3ell, 0<=z<=8ell`. The retained ranges are
`129<=rho<=131`, `-7/50000<=p0=phi_y<=-3/25000`, and
`111/1000<=a<=14/125`. Thus `rho>0` and `p0` never vanishes.
The proper-distance change is at most `131/125000<9<=T_b`, so this is the
actual regular branch. The historical `8e-10` collar is not its final extent.
These are retained analytic branch and rational enclosure results, not a
new trajectory solve. See `BACKGROUND_DERIVATIVE_AUDIT.md:34-62` and
`ACTUAL_POTENTIAL_JET_BOUNDS.json:domain,background_intervals`.

## Full common first seed, including velocity

The actual packet is the smooth flat bump `Q=q_ell(z)` on `[ell,2ell]`,
with `Y=n1^2-n2^2`, `alpha=(rho/rho_b)^(3/2)p0/p0_b` and
`k=alpha_z/alpha`. The source fixes, without a free velocity choice,

```
u=alpha Q/rho^2,             u_t=alpha exp(-z) Q/rho^2,
v=3alpha/(rho^3 p0)(Q_z+kQ),
v_t=3alpha exp(-z)/(rho^3 p0)(Q_z+(k-1)Q).
```

It also fixes `gamma1_yy=-4uY`, `gamma1_ab=2uY gamma0_ab`,
`K1_yy=-2u_tY/rho`, `K1_ab=u_tY gamma0_ab/rho`,
`phi1=vY`, and `Pi1=v_tY/rho`; mixed spatial components vanish.
These are all four physical Cauchy sectors, and the retained reconstruction
solves their linearized constraints. See
`EXPLICIT_COMPACT_SECOND_INITIAL_JET.md:66-109`.

Here is the explicit join to the old conformal ADM seed. Denote the old
Cartesian master displacement and velocity by `q_c,p_c`, to distinguish
`p_c` from the background scalar derivative `p0`. The inherited equations
are in `SECOND_ORDER_LOCALIZED.md:67-72`:

```
h=2q_c(s I-3XX^T),
delta_phi=-3(Dq_c+4q_c)/(2Omega^2 Phi_s),
delta_pi=-3p_c XX^T/(2Omega),
delta_varpi=-3Omega(Dp_c+5p_c)/(2Phi_s),  D=R partial_R.
```

With `H=X1^2-X2^2=sY`, choose

```
q_c=G(s) H,  p_c=K(s) H,
G=Omega^2 u/s^2=alpha Q/s^3,
K=Omega^2 u_t/s^(5/2)=G/b.
```

The last equality uses the actual velocity `u_t/u=exp(-z)=R/b` and the
initial clock conversion `partial_t=R partial_ambient_T`. In particular,
`p_c=q_c/b`; replacing `p_c` with `u_t` would give the wrong seed.
Both profiles are smooth and supported strictly inside the selected
annulus. They even obey the old frozen-profile relation with `omega=1/b`;
the general exact initial theorem independently permits arbitrary smooth
profiles (`EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md:278-301`).

All sectors match. For clarity, `D=-partial_z` and
`Phi_s=p0 rho/(2R^2)` give

```
q_c=alpha Q/R^4 Y,
Dq_c+4q_c=-alpha/R^4 (Q_z+kQ)Y,
Dp_c+5p_c=-alpha/(bR^4) (Q_z+(k-1)Q)Y.
```

The inherited scalar reconstruction is exactly `vY`; the inherited
scalar momentum density is exactly `Omega^4 Pi1` in the Cartesian chart.
The metric and its unit-normal extrinsic curvature also match the displayed
physical values. The checkpoint convention is
`pi=rho_gamma(K^ij-tr(K)gamma^ij)/2`, `varpi=rho_gamma Pi`.
The shell gauge companion's action convention instead doubles both
momenta; its lines 30-44 explicitly permit the checkpoint normalization
when constraints and inverse conversion are changed consistently. This
constant factor does not change the seed or the support after conversion.

The old one-based `Y=n1^2-n2^2` is exactly the same polynomial as a
zero-based implementation's `n0^2-n1^2` after index relabeling. The matrix
is `diag(1,-1,0,0)` in that fixed order in both displacement and velocity.
No angular rotation, amplitude rescaling, or change of family parameter is
being made.

## Selecting the whole correction support and gauge

Apply the exact initial theorem afresh on

```
A'={ell/2<z<5ell/2} x S3.
```

In conformal Cartesian coordinates this is the smooth concentric annulus
`b exp(-5ell/2)<R<b exp(-ell/2)`. It stays away from the shell `R=b`
and the vertex `R=0`; `p0!=0` on a neighborhood of its closure. The
reflection group preserves the radial background and the complete seed.
The retained local adjoint classification therefore has the same kernel
and zero invariant part. The scaled-cover/weighted-isomorphism and
nonlinear mapping estimates remain analytic premises of the retained
theorem, applied at this fixed annulus and `k=12`; no uniform bound as
annulus width tends to zero is needed.

The exact curve `x(e)=x0+e v+A u(e)` has its **entire** nonlinear metric
and momentum corrections supported in `closure(A')`. The correction need
not stay within the first bump's `[ell,2ell]` annulus. The theorem's
`H0^14 x H0^13` zero-extension result prevents an artificial distributional
boundary layer, and gives a fixed open background initial collar `z<ell/2`
and a background inner region `z>5ell/2`.
See the exact initial theorem at lines 209-269.

Choose the permitted initial gauge directly:
`N_e=rho(1+e uY)`, `beta_e=0`, and determine the lapse/shift time
velocities by the exact local relative-harmonic formulas. These formulas
use only the data and first spatial derivatives. Thus all completed
reduced initial perturbations have support in `closure(A')` too; affine
lapse support is smaller, `[ell,2ell]`.
`INITIAL_GAUGE_SUBSIDIARY_AND_PARAMETER.md:62-129` proves this local support
preservation and uniqueness of the remaining five velocities.

This matches the actual common-gauge first one-jet: its first generator
and its first time derivative vanish initially, but
`f_tt=-6rho^2(u_y+4au)` is retained, producing the nonzero mixed `yt`
metric velocity. The first shift **value** is zero; its velocity is not.
The initial scalar and angular metric velocities remain the stated ones.
Simply taking the old precontact evolution with background lapse would
not give this completed tangent. See
`FIRST_HARMONIC_FIELDS_UNPUBLISHED.md:180-196`.

For the exact boundary state map, the full recovered source now supplies
the local flow proof, not just the main theorem's reference. Its vector
field is smooth where `N>0` and `h` has the required Lorentzian signature,
and it vanishes at zero perturbation. Compactness of the chosen shell-time
interval gives a common positive local flow range in a smaller state ball.
Choose once and for all

```
0<h<min(that uniform local flow range,ell/4),
theta(z)=z chi(z/h),  |theta|<=h,
```

with `chi=1` near zero and supported before `z=h`. The resulting map is
invertible, equals the identity outside that collar, and fixes the zero
state. Its initial collar lies strictly inside `z<ell/2`, where the exact
data are background, so it creates no initial support. No numerical value
for the uniform flow range or for `h` is asserted. See
`NONLINEAR_BOUNDARY_STRAIGHTENING.md:27-80`.

## Exact domain inequalities and prescribed-time smallness

Set `L=5ell/2`, `T_t=3ell/2`, `Z=6ell`, and `c_star=2`. The full
correction and completed gauge support, together with the fixed flow
collar, obey

```
L+c_star T_t=11ell/2 < 6ell=Z,    margin=ell/2.
T_t<3ell, Z<8ell,
max_slab(x+r)=3/2+6=15/2<8,      triangle margin=ell/2.
target speed-2 backward reach=ell/8+2(5ell/4)=21ell/8<Z.
```

Thus the whole selected slab lies in the registered global first-field
domain and background rectangle; the target and its histories are
strictly internal. At time `T_t` the perturbation is still zero near the
artificial outer wall. The full variable-wave proof permits homogeneous
Dirichlet data there and proves that this construction has no physical
effect (`QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md:15-20,414-430`). The
physical shell has the original mixed physical and harmonic boundary
rows, including all postcontact history.

On this compact slab the background metric in `(t,z)` is
`rho^2(-dt^2+dz^2+cosh(t)^2 kappa)`. Its radial characteristic speed is
one, `-g0^tt>0`, `g0^zz>0`, and the full spatial quadratic form is
positive relative to a fixed product metric. The induced shell metric is
Lorentzian. The required energy positivity is stronger than timelikeness
of `dt` alone; it holds here and persists in a small state neighborhood.
Compactness allows uniform positive margins and speed below two. The
full boundary flow remains defined in a smaller common neighborhood.

Let `R_state>0` be the chosen fixed state/energy ball and `C_seed` bound
the completed initial perturbation norm by `C_seed |e|`. The retained
whole-interval proof closes once

```
C0 C_seed |e| exp(C1 T_t) < R_state.
```

All constants are finite and depend on this fixed branch, annulus, collar,
slab and regularity. Intersecting this positive parameter interval with
the exact initial IFT interval and the lapse-positivity interval supplies
a common nonzero interval. The theorem expressly permits this shrinkage.
Neither an effective numerical epsilon radius nor a lower bound uniform
in `ell` is required for the existence of the second derivative at zero.
See main theorem lines 42-103 and 198-231; variable-wave proof lines
316-340 and 414-430. No positive claim about a specified finite amplitude
follows from this argument.

## Common first solution and remaining analytic scope

The exact family's first variation and the current retained actual first
pair now have the same background, common initial section, four physical
Cauchy sectors, lapse and shift values, and all harmonic velocities. They
satisfy the same relative-harmonic linear bulk equations, shell coordinate
rows, linearized Israel/scalar laws and harmonic normal row, on both
identical copies. The current master/gauge construction supplies those
physical rows; the actual first gauge has `f=0,b_y=d_y=0` at the shell.
Linear uniqueness for this completed mixed problem identifies the fields
on the needed causal region, independently of an auxiliary outer
construction. To assert agreement on the whole selected slab, apply the
linear finite-propagation theorem to each completed first field: their
initial support is `[ell,2ell]`, their background radial speed is one,
and the physical shell rows are homogeneous at first order. Each is zero
for `z>2ell+t`, hence zero near `Z=6ell` throughout this slab. They
therefore have the same zero perturbation data at the artificial wall,
so whole-slab uniqueness also applies. This uses the linear
energy/uniqueness and finite-propagation theorems, not equality of a formal
tangent in different time variables.

Analytic inputs retained explicitly: the actual regular branch; local
adjoint classification, API and scaled-cover weighted inverse; the
whole-interval quasilinear mixed solver, constraint propagation, parameter
regularity, finite propagation and linear uniqueness; and the current
master/gauge existence and boundary premise. The recovered 51-member
shell capsule now contains the relevant written proofs, and its external
manifest pin authenticates those bytes. Their finite symbolic fixtures do
not independently certify the PDE theorem.

There is no shell-arrival transversality gate for this coordinate
comparison: the current first packet's reference edge time is `ell`, but
the exact nonlinear correction can extend toward `ell/2`, so actual
nonlinear first arrival is not claimed to equal the reference time.
The fixed event has `z>0`. Geometric observer transport and its own
transversality, any finite-amplitude sign, and the general TT projection
argument are separate gates.

`SOURCE_PINS.json` authenticates the 15 scientific source inputs with
checkout-relative paths resolved under an explicit `--repo-root`. Private
planning notes and recovery metadata are not runtime dependencies. The new
checker reauthenticates its externally pinned input manifest and all 51
members of the retained analytic capsule, checks the displayed exact
inequalities, and verifies nine seed identities including both scalar
sectors. It rejects three nonzero wrong seed formulas. Normal and `-O`
outputs agree byte for byte. Separate subprocess controls mutate actual
source copies and reject them before algebra; numerical edge-case controls
reject loss of each strict margin. These checks establish the stated
arithmetic/algebra facts and source binding, not an analytic theorem.
