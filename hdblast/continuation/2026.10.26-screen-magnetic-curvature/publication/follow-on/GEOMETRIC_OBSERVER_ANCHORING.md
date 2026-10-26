# Follow-on: what must be added to anchor the observer geometrically

This note is separate from the frozen screen-magnetic coefficient capsule. It
changes neither the accepted metric coefficient nor its registered-coordinate
interpretation. The argument below is the standard fixed-vector obstruction
for an equivariant construction at symmetric data. It is not a new physical
no-go theorem, and it does not rule out clocks, observers, or measurements that
include additional physical preparation data.

The present background data do **not** select the registered observer
`U0=rho^-1 partial_t` by a local covariant finite-jet construction. A marked
initial section of the shell, with a clock/congruence construction and its
perturbations, is one concrete way to supply the missing datum. An intrinsic
shell clock extended along the scalar normal can preserve the desired `U0`;
full-spacetime geodesic transport generally does not.

## Precise hypothesis and obstruction

Consider the background in a collar of an isotropic shell,

```
g0=dy^2+rho(y)^2 h_dS4,
phi0=phi0(y),
Shell0={y=y_s},
h_dS4=-dt^2+cosh(t)^2 kappa_S3.
```

Assume `rho>0`. The supplied structure may include spacetime orientation,
time orientation, the shell side, and the tensors naturally induced on the
shell. It does **not** include a marked shell event, marked initial spacelike
section, preferred congruence, non-invariant clock field, or asymmetric boundary
conditions. The statement concerns points where a smooth covariant construction
is supposed to be defined on these exact background data.

A local finite-jet observer rule is a map from a finite jet of those data at a
point to a tangent vector, equivariant under the allowed local diffeomorphisms.
Its covariance implies invariance under every automorphism of the supplied
background that fixes that point. Finite differential order is not essential
to the following pointwise argument, but it specifies the class at issue.

All de Sitter isometries lift to isometries of `g0` by acting trivially on `y`.
They preserve `phi0`, the shell, its induced metric and second fundamental form,
and their covariant derivatives. At a point, their proper orthochronous
stabilizer is `SO^+(1,3)`. Its action on the five-dimensional tangent space is

```
T_p M = span(n0) direct-sum V_(1,3),
n0=partial_y,
SO^+(1,3) acts trivially on span(n0) and by the Lorentz vector representation
on V_(1,3)=T_p dS4.
```

The common fixed space of that representation is exactly `span(n0)`. To see
this directly, spatial rotations force the three spatial components of a
fixed vector in `V_(1,3)` to vanish. A nontrivial proper orthochronous boost
then forces its time component to vanish. Thus covariance requires

```
U[g0,phi0,Shell0](p) = a(p) n0.
```

Such a vector is spacelike or zero, never a unit timelike vector. Therefore
**no covariant vector rule defined at these symmetric data can produce a
unit timelike observer there**. In particular, no smooth local finite-jet
rule can reproduce `rho^-1 partial_t` from the listed inputs alone.

The obstruction already follows from equivariance at the background point;
smoothness is not needed for that step. Smoothness matters when asking for a
rule that extends regularly to this point from perturbed metrics.

Time orientation does not remove the obstruction. Proper orthochronous boosts
preserve the future cone while moving each future unit timelike vector. They
select the future sheet of the unit hyperboloid, not a point on that sheet.
Similarly, spacetime orientation and a chosen shell side still leave the same
proper orthochronous stabilizer.

An isotropic constant-`y` shell provides its spacelike normal and tensors
proportional to its de Sitter metric, such as
`chi_shell=(rho_y/rho) h_shell`. It supplies no preferred timelike tangent
direction. This remains true of all its covariant jets. If the shell is given
only on a region whose boundary or endpoint is itself marked, that marking
may reduce the symmetry; the theorem then applies to the actual augmented
inputs, not to an artificially enlarged symmetry group. If `dphi0` vanishes,
it can supply less directional information, not a timelike direction that
escapes the same stabilizer argument.

The coordinate `t` in the displayed global de Sitter slicing already encodes
additional structure: the chosen spacelike `S3` section at `t=0` and the
normal congruence it determines. General de Sitter isometries do not preserve
this slicing. Calling `t` a scalar under a coordinate change does not make it
a scalar constructible from the invariant metric and matter data alone.

## Why perturbative eigenvectors do not evade the hypothesis

A perturbation can break the background symmetry and split a degenerate
curvature or matter eigenvalue. An eigenvector may then be chosen smoothly
along one specified one-parameter family, with a chosen limiting direction.
That observation is compatible with the obstruction.

A universal covariant rule that is continuous at the symmetric background
must have the same value on every approaching family. Applying a background
stabilizer transformation to a family rotates or boosts its preferred
limiting eigenvector while leaving its limiting background unchanged. A
unique continuous limit would therefore have to be fixed by the stabilizer,
which cannot be timelike. Normalizing a vector that vanishes at the background,
or dividing perturbative eigenvalue differences by the perturbation parameter,
can fail regularity or produce a direction-dependent limit.

Supplying the first perturbation `q` itself as an additional independently
specified tensor datum is a different problem: it can reduce the stabilizer.
It is not a rule of the instantaneous `g,phi,shell` jets alone at the degenerate
background. A pathwise choice can be useful, but its path dependence and
regularity assumptions must be stated. None of this excludes a physical
clock field, prepared observer, selected boundary congruence, or marked initial
surface as valid extra data.

## A concrete additional datum that preserves the desired background

One minimal direction-selecting datum is a marked spacelike three-section
`Sigma_s` of the physical shell, together with its smooth perturbation and the
existing time orientation. Its future unit normal **within the shell** selects
an initial observer. In a collar before intrinsic Gaussian-coordinate caustics,
extend those normals along intrinsic timelike shell geodesics and let `sigma`
be signed shell proper time measured from `Sigma_s`.

For the unperturbed section `Sigma_s0={y=y_s,t=0}`, this gives
`sigma0=rho_s t`, where `rho_s=rho(y_s)`. A fixed constant rescaling
`tau_shell=sigma/rho_s` therefore has `tau_shell0=t`. The rescaling only fixes
clock units; even without it the normalized background observer below is the
same. The marked section and its perturbed embedding are extra preparation
data, not consequences of shell isotropy.

Assume in a shell-to-target collar that the physical scalar gradient is
nonvanishing and spacelike. Orient its unit normal `n_phi` to agree with
`n0=partial_y`, and solve the transport problem

```
n_phi(tau)=0,                    tau|_shell=tau_shell.
```

Require the scalar-normal flow to meet the shell transversely and uniquely in
the collar, and require `d tau` to remain timelike. Define

```
U_geo=-grad_g(tau)/sqrt(-g^-1(d tau,d tau)),
n_geo=n_phi.
```

These vectors are future/unit and orthogonal by construction. At the background,
`n_phi=partial_y` and `tau0=t`, so `U_geo0=rho^-1 partial_t`, exactly the
registered observer. This construction is covariant **for the augmented data**
that include the marked section and the specified transport rule. This is a conditional construction on a smooth field family; it does not
construct that family or its second metric. It is
nonlocal across the collar and is not a local finite-jet rule of the original
symmetric fields. At the perturbed shell, the resulting observer need not
equal the intrinsic shell geodesic congruence unless `n_phi` is the shell
normal or an equivalent compatibility condition is supplied. The assertion
here is exact background matching, not that additional perturbed matching.

Do not replace the intrinsic shell construction by full-spacetime freely
falling observers without checking the background congruence. For the desired
fixed-`y` observer,

```
nabla_(U0) U0 = (rho_y/rho) n0.
```

It is generally accelerated in the five-dimensional spacetime, even though its
shell projection is an intrinsic shell geodesic. A full-spacetime geodesic
clock therefore need not reproduce `U0`.

The geometric target sphere and its axes must also be specified, for example
by a clock value, a scalar value, and marked angular labels transported from
`Sigma_s`. This note supplies a candidate prescription and its first-variation
inputs; it does not claim that this augmented construction has already been
implemented or certified.

## Finite first-variation data for this candidate

Use the accepted frame conventions `e0=rho^-1 partial_t`, `ey=partial_y`,
`ea=R^-1 E_a`, `R=rho cosh(t)`, and write the actual scalar perturbation
`phi1=s(y,t)Y(Omega)`. Put `p=partial_y phi0`, and assume `p` is bounded away
from zero on the chosen collar. The sign of the normalized gradient is chosen
so that its background vector is `ey`: `n_phi=sgn(p) grad(phi)/|grad(phi)|`.
The formulas below use the signed `p`, whose sign is constant on a connected
nonvanishing-gradient collar.
The registered metric amplitudes include
`q_yy=A Y`, `q_y0=B Y`, `q_ya=C D_aY`, `q_0a=E D_aY`.

The scalar-normal first variation is

```
(n_phi)1 = -(A Y/2) ey
           + (B-s_t/(rho p))Y e0
           + (s/(pR)-C) D^aY ea.                              (4)
```

This follows by differentiating both `g^-1 dphi` and its normalization; it
must not be replaced by the spatial normal to `y` in the registered time slice.
Writing `tau=t+lambda tau1+O(lambda^2)`, the radial transport equation gives

```
partial_y tau1 = [-B/rho+s_t/(rho^2 p)]Y.                     (5)
```

Let `beta_pull(t,Omega)` denote the first variation of the boundary clock
pulled back by the specified perturbed shell embedding, and let `zeta^t` be
that embedding's first coordinate-time displacement. Define the effective
Eulerian boundary datum `beta=beta_pull-zeta^t`, because `tau0=t`. This
subtraction includes the time-identification term; it is unnecessary only
for an explicitly time-preserving shell identification. The boundary datum
also depends on the marked section and intrinsic-clock variations. Then,
at this order,

```
tau1(y,t,Omega) = beta(t,Omega)
  + integral_(y_s)^y [-B(y',t)/rho(y')
                     +s_t(y',t)/(rho(y')^2 p(y'))]Y(Omega) dy'. (6)
```

The boundary evaluation and its identification have to be specified as part of
the augmented data. A displacement of the radial starting point alone produces
no extra first-order term here because `partial_y tau0=0`.

Compared at the same registered event with the accepted coordinate-normal
pair `(U_reg,n_reg)`, whose zeroth-order pair agrees, define

```
Delta U1=b n0+v^a ea,
Delta n1=b U0+w^a ea.
```

Equations (4)–(6) give

```
b=(B-s_t/(rho p))Y,
v_a=-(rho/R)D_a tau1,
w_a=s D_aY/(pR).                                             (7)
```

The common radial-boost coefficient `b` is required by differentiating
`g(U,n)=0`; it is also equal to `-rho partial_y tau1` by transport. These are
candidate observer-correction inputs, not evaluated enclosures.

The inherited actual first-field two-jets suffice to construct the complete
linear Weyl tensor `C1` kinematically. Equations (5)–(7) additionally need only
`B,s,s_t`, background `rho,p`, the boundary-clock variation and its angular
derivative, and a certified radial integral on the selected collar. In the
standing quadrupole sector the angular differentiations are exact finite
harmonic operations. The existing differentiated master/gauge bounds through
order three are sufficient for those local metric jets; no fourth metric
derivative is indicated. A lower bound on `|p|`, collar existence, boundary
clock data, and a transport/integral error ledger are additional premises or
work, not automatically supplied by the frozen TT response certificate.

## Exact observer-change identity to apply before making a sign claim

Here compare the same metric family and the same event, with two smooth
oriented observer/screen frames that agree at `lambda=0`. Let
`C(lambda)=lambda C1+lambda^2 C2/2+...`; the full background Weyl is zero.
The difference between the first frame variations is a Lorentz-algebra
variation. After matching the first angular rotation, write

```
Delta U1=b n0+v^a ea,
Delta n1=b U0+w^a ea,
Delta(ea)1=v_a U0-w_a n0.                                    (8)
```

In an adapted oriented orthonormal frame, epsilon has constant numerical
components. It therefore already includes its metric/volume variation through
the varying frame. For the same definition
`H_ab=(1/2)epsilon_(a^{cd}C_{b)Ucd}`, differentiation gives

```
Delta H''_ab = epsilon_(a^{cd} [
    C1(Delta(e_{b)})1,U0,ec,ed)
   +C1(e_{b)},Delta U1,ec,ed)
   +C1(e_{b)},U0,Delta(ec)1,ed)
   +C1(e_{b)},U0,ec,Delta(ed)1) ].                            (9)
```

The displayed `a,b` symmetrization is normalized. The overall coefficient is
one because the second parameter derivative contributes a factor two and the
definition of `H` has a factor one-half. Equivalently,
`Delta H''=2 Delta(projection)_1(C1)`.

There is no `C2` difference: the metric is held fixed. Second observer/frame
variations multiply `C0=0` and vanish. An additional relative first angular
rotation acts on `H1`; for the accepted scalar first sector the complete
`H1=0`, so that rotation contributes zero. Nonetheless the other terms in (9)
need not vanish merely because this particular projection of `C1` vanishes.
One must retain the full `C1`, whose other components need not be zero.

Equation (9) is the required correction to evaluate and project before carrying
a registered-coordinate coefficient sign over to a new physical observer.
It does not compute that correction or assert its sign. For a moving target
sphere, its embedding and angular-label pullback must also be specified and
audited. A covariance or gauge-invariance claim requires the full augmented
construction to transform with its data; the symmetry obstruction and the
formula alone provide no such certificate.

## Bounded proof checks and scope

The adjacent proof checks use a rational proper orthochronous boost with
`cosh(chi)=5/3`, `sinh(chi)=4/3`, along three independent screen axes,
and signature `(-++++)`. They verify metric preservation, determinant one, time orientation,
and that the common fixed space is radial. A nontrivial boost sends the
coordinate time unit vector to a different future unit timelike vector.
These finite matrices illustrate the stabilizer step; the complete argument
above identifies why the continuous background isometry group contains them.

The independent audit also checks the Lorentz-algebra form (8) and the factor
in (9). An explicitly synthetic algebraic Weyl example satisfies the curvature
symmetries, Bianchi identity and tracefree condition, has this screen `H1=0`,
and has a nonzero observer correction. It disproves the unrestricted inference
`H1=0` implies observer independence of `H2`; it is not an actual-field example. These are mathematical consistency checks
of the follow-on proposal,
not a new evolution computation, new metric, complete second-order solution,
or external specialist review.
