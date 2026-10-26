# Independent audit of the observer anchoring follow-on

The proposed obstruction and first observer-change formula are correct under
the stated hypotheses. This is a standard isotropy/naturality argument, not a
new impossibility theorem about physical clocks. The checked scripts use no
actual first-field values and establish no new sign for the registered result.

## Exact scope of the obstruction

Let the input consist of a metric `dy²+rho(y)² g_dS4`, scalar fields depending
only on `y`, a constant-`y` shell with de Sitter invariant intrinsic/extrinsic
geometry, and compatible orientation and time orientation. No marked section,
worldline, clock, congruence, symmetry-breaking boundary data, or perturbation
tangent is part of this input. At any bulk point, the connected de Sitter
isotropy group fixing that point acts on its tangent space as

```
T_p M = span(radial) direct-sum T_p dS4,
G_p = SO^+(1,3),
Fix(G_p)=span(radial).
```

The same group fixes the constant-radius shell as an embedded set and fixes
all its invariant geometric jets. In particular, a shell normal supplies a
radial spacelike direction; it does not supply a shell timelike rest frame.
The shell induced metric and isotropic second fundamental form are invariant
under the full tangent Lorentz group. A specified material rest frame would
be an additional datum and would change these premises.

Every finite jet of the background inputs is fixed by these actual
isometries. Equivariance of a natural vector assignment therefore forces its
value at the point to be fixed by `G_p`. The only possible output is radial
(or zero), hence cannot be unit timelike. This pointwise argument needs no
smoothness assumption; smoothness becomes relevant when asking for a regular
assignment on a neighborhood of the background. The proper orthochronous
subgroup preserves both supplied orientations, so time orientation does not
remove the obstruction. It picks a component of the timelike cone, not a
particular vector inside it.

The argument applies to the indicated input category. It does not exclude
clocks defined from additional symmetry-breaking data or chosen transport
paths. Even a deterministic equivariant nonlocal rule using only the *same*
exactly symmetric input would still face its stabilizer, but no general
nonlocal-clock impossibility is required here.

A curvature/scalar-gradient eigenvector selected after perturbation can have
a smooth limit along one chosen family. At the symmetric background its
relevant tangential eigenspace is Lorentz degenerate. Limits can depend on the
family's first tangent `q`; normalization of an `O(lambda)` quantity can also
be singular as a functional of the original fields. Giving `q` or choosing an
eigenvector branch at the limiting background adds data. It does not yield a
single regular equivariant unit timelike functional of the unaugmented
background jets.

## Relative observer and screen variation

Compare the same metric family, same Weyl tensor, and same event, using two
oriented adapted orthonormal frames that agree at parameter zero. Write the
*relative* first variations as

```
delta U1 = b n + v^a e_a,
delta n1 = b U + w^a e_a,
delta e_a1 = v_a U - w_a n,                   (zero relative screen rotation)
```

where all vectors on the right are background vectors. These are precisely
the seven independent variations allowed by first orthonormality. A relative
SO(3) rotation can be included if the axes are to change; its contribution is
just the corresponding rotation of the first magnetic tensor and disappears
when that complete tensor vanishes.

In these adapted components the numerical epsilon symbol is constant. If
`C(lambda)=lambda C1+lambda² C2/2+...`, then direct multilinearity gives

```
Delta H''_ab = Sym_ab epsilon_a^{cd}
  { C1(delta e_b1,U,e_c,e_d) + C1(e_b,delta U1,e_c,e_d)
   +C1(e_b,U,delta e_c1,e_d) + C1(e_b,U,e_c,delta e_d1) }.
```

Here `Sym_ab` denotes normalized symmetrization. There is one epsilon prefactor because the second-derivative cross factor `2`
cancels the definition's factor `1/2`. No extra free-`a` slot is missing:
normalized symmetrization includes its counterpart. Radial-screen changes
enter through `delta e_a1`; the epsilon/volume changes are absorbed by using
adapted orthonormal components. In a fixed-coordinate expression, their
explicit variations must be retained instead. Second frame variations multiply
`C0=0` and vanish. This is a same-event comparison; altered point/clock/radius
anchoring requires its own displacement bookkeeping.

## Concrete additional datum and finite inputs

A marked spacelike `S³` section of the shell is sufficient locally, together
with future time orientation, to choose the shell normal geodesic congruence.
Construct the shell clock from its timelike eikonal equation with zero value
on that section. Before caustics, transport it off the shell along a specified
radial collar. A nonzero spacelike scalar gradient supplies one possible
collar direction; a Gaussian shell-normal collar is another. Neither is
silently identified with the registered-coordinate normal away from the
background.

When the marked background section is the registered `t=0` sphere, the shell
clock is proportional to `t`; radial transport preserves this property and
its normalized normal is `rho^-1 partial_t`. A positive clock calibration is
needed to identify the numerical apex time, although it does not affect the
unit observer. Existence of a valid timelike clock and a non-caustic collar,
the perturbed section/embedding data, and the chosen transport prescription
must be stated.

For the first anchoring correction one needs the linearized clock and radial
transport equations, their initial/boundary data, and their certified first
variations. The local Weyl input `C1` needs the metric first perturbation and
its covariant derivatives through order two. Scalar-gradient transport needs
the scalar first perturbation through first derivatives; differentiating the
transport equation as needed remains within the available smooth first-field
jet control. Existing C3 first-field data can supply local coefficients, but
cannot invent missing clock/section data or certify new transport integrals
without additional work. No second metric solution or curvature sign is
computed here.

## Exact algebra controls

`check_anchor_witness.py` verifies with rational arithmetic:

- three proper orthochronous boosts with `cosh=5/3`, `sinh=4/3` preserve the
  five-dimensional metric and fix the radial vector;
- the stacked fixed-vector equations have rank four, leaving exactly the
  radial line;
- all seven independent relative adapted-frame generators obey first
  orthonormality;
- a synthetic algebraic Weyl tensor is tracefree and satisfies all curvature
  symmetries and Bianchi, has `H1=0`, and acquires a nonzero second magnetic
  change under a first observer boost.

The last example only prevents the invalid inference `H1=0 => observer
independence of H2`. Its displayed rational entry is not an evaluation or
sign claim about the registered actual solution.
