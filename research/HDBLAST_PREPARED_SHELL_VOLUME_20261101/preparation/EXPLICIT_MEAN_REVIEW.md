# Independent review of the explicit second monopole preparation

1 November 2026. **ACCEPT, conditional on the retained weighted-IFT,
background and shell-evolution premises.** No material error was found in
the explicit design. This review is separate from the accepted publication
and changes no public or frozen source file.

Reviewed design SHA-256:
`d2db3ac342a9225fc6eef54461f2829bd52cf6cd2e8024aa346fafe6146ab01f`.
The cited public sources were independently authenticated against their
recorded hashes and byte counts at the retained source version. The prior
weighted-mean derivation also matches its earlier reviewed hash. Current portable source records and finite checks are in the adjacent
`../gauge/SOURCE_PINS.json` and `../gauge/references/` files.

## Explicit Hamiltonian and momentum data

The normalization is correct: the second constraint equation is
`P x2=-D^2 Chat[v,v]`, since the trial uses `epsilon^2 B/2`. For increasing
outward y, `yL=y(2ell)<yU=y(ell)`, and

```
J=integral_(yL)^(yU) rho^2 b dy
 =integral_ell^(2ell) rho^3 Q dz >0.
```

The metric identity and its moment give

```
Pg(f gamma0)=C[f],
Pg(Lie_(xi d_y)gamma0)=rho^-4(rho^4 p^2 xi)',
integral rho^4 C[f]dy=-(9/2)integral rho^2 fdy.
```

Thus `f=(2I/(9J))b` makes the total integral of
`rho^4(QH+C[f])` exactly zero. Differentiating the chosen negative primitive
then yields `Pg h=C[f]-(QH+C[f])=-QH`, with the correct sign. The earlier
strict moment result gives `I<0`, hence `f<0` in the bump interior, but this
sign is not required for the right-inverse construction.

The source and bump are smooth and flat at the packet endpoints. Their
zero-total-integral primitive vanishes on both exterior sides, and all its
endpoint derivatives vanish. The lower bound on `|p|` makes division by
`rho^4 p^2` smooth near the support; extend xi by zero outside rather than
requiring a nonzero p on an unrelated larger region. The resulting h is
smooth and supported strictly within the correction annulus. All required
finite Sobolev norms are finite even if their numerical sizes are large.

For the normalized contravariant momentum tensor in the fixed frame,

```
Pm(C,D)=-2[C'+3a(C-D)].
```

Therefore `C=0`, `D=-QM/(6a)` is correct. The proposed alternative
`D=QM/(3a)` would produce `+2QM`, so the stated correction is essential.
The retained `a>=111/1000` makes the chosen division legitimate. QM is
smooth and flat, so this momentum datum has the same compact support.
An invariant tangential one-form on S3 is zero; no angular momentum
monopole remains missing. Together these identities give `P B=-P0 Q` in
the full scalar/one-form constraint range.

## Averaging and the conditional IFT

The added fixed compactly supported polynomial term changes only the
prescribed seed constants in the retained nonlinear estimates. It does not
alter the derivative in U at the background. The same finite-k IFT therefore
applies after shrinking the parameter interval, and gives

```
L U''=-w^-2(Q+P B), U(0)=U'(0)=0.
```

The averaging argument is sound. Define P0 on the full weighted spaces by
Haar averaging the natural SO(4) action. Radial background coefficients and
weights make P, its adjoints, A and L commute with this averaging. Background
covariant Sobolev norms make rotations isometries; equivalent weighted norms
give a bounded projector. Averaging compactly supported smooth fields keeps
their radial support, so it extends to the dotted closures.

The proof correctly avoids claiming that each rotation preserves the
coordinate-reflection invariant subspace. What is needed is only that P0
maps the relevant invariant domain back into itself. Its scalar and one-form
outputs are radial, hence reflection invariant. A general alternate
justification is that SO(4) is normal in O(4): any selected reflection
commutes with the Haar projector, so reflection-invariant input retains that
invariance after averaging. This also handles the tensor/density outputs
with their natural pullback and the fixed background density.

Consequently `L P0 U''=0` in the reflection-invariant dotted domain. The
retained injectivity on that domain implies `P0 U''=0`; no unrestricted
full-space inverse or individual-rotation invariance of the domain is needed.
Commutation with A gives `P0 x''=B`. Scalar and scalar-density second data
remain zero in every angular sector because the correction acts only on
metric and metric momentum.

The independence from remaining radial weight choices is exactly limited
to this newly prescribed second monopole. It does not assert weight
independence of the complete nonlinear family, higher parameter derivatives,
or the other second angular components.

## Physical conversion and TT/evolution applicability

Substitution of `C=0`, `D=-QM/(6a)` into the already independently checked
canonical conversion gives exactly the displayed mean K2 and coordinate
metric velocities. In particular the radial linear contribution is
`QM/(3a)`, while the angular linear contribution is zero. The scalar-density
conversion still gives `<Pi2>=-uj/3`; its first-lapse cross term cancels in
the coordinate scalar velocity, so `phi2=partial_t phi2=0` initially.

The added radial metric and momentum have isotropic angular blocks and zero
angular TT projections. The inherited weighted-adjoint correction also has
zero angular TT projection. At `pi0=0`, the quadratic canonical conversion
depends only on the unchanged first field, and its angular term is pure
trace. The first-lapse cross term is also pure angular trace. Thus the
second TT displacement and velocity remain zero.

The support stays inside the same fixed annulus, leaving the physical shell
collar exactly background. The retained relative-harmonic completion and
corner argument therefore apply to this new initial curve. Conditional on
the same prescribed-slab theorem, the previous causal margins persist after
family-dependent parameter shrinkage. The full first field is unchanged by
the retained linear uniqueness argument. The first-only shell coefficients
and the TT second-response identification can therefore be reapplied.

This does not identify all second data with the prior weighted-only family,
produce an effective parameter radius, or compute the evolved mean. It is
a new mathematical preparation of initial data, not a physical mechanism
generating the packet or a radiation/volume-acceleration result.

## Checks and limitations

The original independent preparation checker verified 14 finite algebraic/source
claims, including the integrated coefficient, primitive cancellation,
momentum factor and its deliberate wrong alternative, physical conversion,
parameter factorial, and nonvanishing source-box denominators. It makes
zero oracle calls. Haar averaging, weighted closures, the IFT and evolution
applicability are the analytic arguments reviewed above; the finite checker
does not certify those retained analytic theorems. This is internal
AI-assisted review, not external specialist review.

Portability note: the reviewed-design hash above identifies the original
reviewed text. The portable design makes documentation/path changes only;
its current bytes are authenticated in gauge/ESSENTIAL_PINS.json. The
original preparation-only checker is not a runtime dependency of this
component; the gauge checkers reproduce the stated conversion tests and
add direct five-dimensional harmonic and conditioned-gauge checks.
