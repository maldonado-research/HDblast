# Next physical gate: an observer-defined screen magnetic Weyl coefficient

The most concrete next target is a symmetric magnetic Weyl tensor on the
three-dimensional angular screen. Its linear TT contribution uses the
already certified, strictly negative `u_x`, rather than the unresolved
second-derivative signs. The complete quantity also contains quadratic
first-field terms that have not been evaluated. This is a feasibility
derivation, not a curvature sign, gauge-invariant observable, or origin claim.

## Define the quantity before computing it

For a future unit timelike observer U and an orthogonal unit spacelike
radial vector n, let `s=g+U⊗U−n⊗n` be their screen projector. Orient the
screen with the spacetime orientation and the ordered pair `(U,n)`, and
denote its three-form by epsilon. Define the full screen tensor

```
H_ab = (1/2) epsilon_(a^{cd} C_{b)Ucd},
```

with all displayed spatial indices screen-projected and parentheses
denoting normalized symmetrization. This is a particular observer/radial
projection of the five-dimensional Weyl tensor. The full observer magnetic
Weyl has three indices on a four-dimensional rest space; it is not the
usual four-dimensional spacetime magnetic two-tensor. A four-dimensional
FLRW vanishing theorem cannot be imported here.

At background level choose `U0=rho^-1 partial_t`, `n0=partial_y`, and
`e_a=(rho cosh t)^-1 E_a` with the registered sphere orientation. A concrete
route to extend this prescription geometrically is to take n from the
normalized, spacelike scalar gradient, with its sign fixed to n0; specify
a clock on the physical shell and extend it by `n·d tau=0`; then set U to
the future unit normal of tau. The background shell clock must be scaled
to give `tau0=t`. This requires a nonvanishing scalar gradient, timelike
clock gradient, a valid shell-to-point transport collar, and specified
shell clock/embedding perturbations. These are proposed additional
anchoring data, not premises already established by the TT certificate.
A covariantly prescribed initial clock and its shell/radial transport must
be shown to preserve U0; a geodesic proper-time clock generally changes
that background congruence. The transport needs its own first-variation
ODE/integral error ledger, which has not been constructed.

Anchor the target sphere by `tau=t_apex` and the background apex scalar
value, and transport angular labels from an explicitly marked shell
sphere. Fix an oriented triad by the chosen screen metric and labels.
This specifies what observer, radius, time, axes, and angular projection
mean. The simpler registered-coordinate U,n prescription is immediately
usable for a gauge-fixed diagnostic, but cannot silently stand in for
this geometric construction.

## Two separate vanishing checks

The background metric `dy^2+rho(y)^2 ds_dS4^2` is conformally flat for
positive rho. With `dz=dy/rho`, remove the conformal factor and consider
`dz^2+ds_dS4^2`. Its Ricci tensor is `Ric_zz=0`, `Ric_ij=3h_ij`, and its
scalar curvature is 12. In five dimensions its Schouten tensor is
`P_zz=-1/2`, `P_ij=h_ij/2`; direct substitution gives `Riem=g wedge P`,
including the mixed z components. Hence the full background Weyl is zero.

Full first-order magnetic Weyl does **not** vanish merely because a
perturbation is scalar on S3. For example, on a time slice a local scalar
extrinsic-curvature perturbation
`delta K_ab=mY gamma_ab`, `delta K_yy=-2mY`, `delta K_ya=0`
has zero angular momentum Ricci component but Codazzi gives a nonzero
`delta C_0abc=m(gamma_ab D_cY−gamma_ac D_bY)` for nonconstant Y, up to the
stated curvature convention. This kinematic example is not an assertion
about the registered on-shell first solution; it rules out the generic
dimensional shortcut.

The **entire defined symmetric screen tensor** H has a stronger scalar
selection rule. For scalar harmonics, angular extrinsic-curvature terms
are combinations of `Y kappa_ab` and `D_aD_bY`, with mixed radial terms
proportional to `D_aY`. Their angular Codazzi curls, including radial
connection terms, reduce to metric times a scalar gradient. Commuting
derivatives of the scalar Hessian uses the unit-sphere constant curvature.
The epsilon dual is antisymmetric in the two remaining angular indices,
so symmetrization kills it. Weyl Ricci-subtraction terms have the same
property. Thus `H^(1)=0` as a complete screen tensor for the registered
scalar first sector, not only after TT projection. Background Weyl zero
also removes first-order changes of U,n from this statement.
The scalar sector here includes all registered gauge-generated scalar
metric amplitudes and the scalar/radial matter perturbation through its
metric response. A field-equation substitution must retain those matter
terms; no fluid or scalar equation is used to erase them in this kinematic
Codazzi argument.
Explicitly, `D_c Hess_bd Y−D_d Hess_bc Y=kappa_bc D_dY−kappa_bd D_cY`.
Contracting the right side with `epsilon_a^{cd}` gives
`2 epsilon_ab{}^d D_dY`, whose symmetric part is zero. The same reduction
applies to the scalar trace term and mixed radial-connection terms.

This gives a promising, but conditional, Stewart–Walker route: if U,n and
the anchoring are genuinely transforming geometric data and the complete
first-order identity is independently certified, then the usual second
gauge transformation terms `L_xi2 H0 + L_xi1^2 H0 + 2 L_xi1 H1` vanish.
No such invariant certificate is claimed here. Different physical observer
prescriptions can still give different H2 through their first-order boosts
acting on nonzero full C1. Vanishing of a TT projection alone would not
suffice for this argument.

## A usable TT formula with no second x/r derivative

Use `R^a_bcd=partial_c Gamma^a_db−partial_d Gamma^a_cb+...` and positive
time-slice extrinsic curvature `K=gamma_dot/(2N)`. For a pure linear TT
perturbation with registered orthonormal coefficient `h_ab=u(y,t)T_ab`,
the coordinate angular metric is `(rho cosh t)^2(kappa_ab+u T_ab)`.
Writing `c=cosh t`, its extrinsic curvature is

```
K_ab = (tanh t/rho) gamma_ab + (rho c^2/2) u_t T_ab.
delta C_0abc = [u_t/(2 rho^2 c)](D_c T_ab−D_b T_ac),
delta H_ab = [u_t/(2 rho^2 c)] curl(T)_ab,
curl(T)_ab = epsilon_(a^{cd} D_c T_{b)d}.
```

The background expansion times u cancels in the all-angular Codazzi
derivative; TT transversality removes the mixed Ricci term. These
identities were derived independently twice. Reversing the curvature or
orientation convention changes the corresponding sign. With scaled time,
`u_t=u_x/ell`. Thus the known first derivative supplies the linear
second-response contribution without `u_xx`, `u_xr`, or `u_rr`.
The factors are explicit: `delta^2 g_AB=rho^2 c^2 u T_AB`, the time normal
is `rho^-1 partial_t`, and each angular orthonormal derivative contributes
`(rho c)^-1 D_kappa`. Codazzi therefore converts its coordinate coefficient
`rho c^2 u_t/2` by three angular frame factors to `u_t/(2rho^2 c)`.
The equation uses the second parameter derivative k itself, so no further
factor of one-half is inserted from `epsilon^2 k/2`.

SO(4) equivariance and inequivalent scalar/vector/TT representation types
mean that TT projection of this *linear* operator sees only the TT part
of the second metric. The required next audit is the exact curl matrix
on the frozen `T+,T−` basis, including normalization and orientation.
Its eigenvalues must be derived in that basis, not assumed from labels.
Equal TT amplitudes need not make the spatial curl tensor vanish.

## Precise next deliverable and bottlenecks

Export and independently check a finite angular coefficient map for

```
Pi_TT H''(0)
 = Pi_TT L_H[k] + Pi_TT Q_H[q,q; U1,n1,screen1],
g(epsilon)=g0+epsilon q+epsilon^2 k/2+... .
```

The quadratic map must retain the complete Weyl Hessian, its Ricci and
scalar subtractions, and first observer/screen/volume corrections acting
on C1. Second observer variations multiply C0 and vanish. Existing actual
first metric two-jets, hence inherited master/gauge C3, suffice for the
quadratic curvature terms. The TT linear term needs only u_x and
background values. No full 23-channel second-derivative certificate is
indicated by this representation argument. The chosen observer transport
may require new first-order radial integrals and shell data; that is a
distinct input bottleneck, not a demand for fourth metric derivatives.

The first milestone should therefore be an exact map and independent
all-point checks of `C0=0`, full `H1=0`, TT block isolation, curl
normalization, and the observer correction terms. Then propagate actual
interval inputs at the already declared apex into the complete projected
H2. Accept a bounded interval even if it contains zero. The signed u_x
alone cannot establish the sign or nonvanishing of the complete quadratic
curvature coefficient. This would move the work from a metric amplitude
toward a specified observer's local curvature diagnostic; detector
response, energy, and cosmological origin remain further questions.

The local literature note describes Compere–Hoque, arXiv:2610.26623, as
linear localized-source perturbation theory on four-dimensional spatially
flat FLRW. Its distinction between gauge, origin/axes, and measurement is
useful here. None of its field equations or magnetic-sector claims were
imported. This feasibility note used local files only and reran no source
or response campaign.

The adjacent standard-library `check_screen_angular_identities.py` checks
all 625 product Schouten/Riemann components, the reduced scalar-curl
identity, and the Codazzi-to-symmetric-curl sign/factor on all 18 independent
angular tensor derivative components. Isolated normal and optimized
outputs agree. It does not substitute for the full registered coefficient
map or the proposed observer transport proof.
