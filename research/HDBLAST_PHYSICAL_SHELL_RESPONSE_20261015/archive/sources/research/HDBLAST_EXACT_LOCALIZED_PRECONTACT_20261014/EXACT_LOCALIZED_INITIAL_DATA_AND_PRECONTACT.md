# Exact localized initial data and a short precontact family at fixed regularity

14 October 2026. Separate analytic continuation with completed internal
cross-review of the weighted initial-data argument. It does not modify
the frozen capsules, certify novelty, claim spatial `C∞` regularity on
a common parameter interval, or solve nonlinear shell contact. All
regularity claims concern one fixed, sufficiently high finite class.

## Theorem and short-time corollary

**Theorem (exact localized initial-data family).** Under the smooth
background, support and reflection hypotheses stated below, fix an integer
`k≥12`. Every smooth compactly supported `Γ`-invariant datum `v`
solving the linearized Einstein–scalar constraints is tangent to a
parameter-smooth family solving the exact bulk initial constraints. Its
metric and contravariant metric momentum density have regularities
`H^{k+2}` and `H^{k+1}`, respectively; its scalar and scalar momentum
density can be prescribed as their background values plus their first
parameter derivatives `εψ_v` and `εϖ_v`. The corrections equal
zero outside one fixed closed annulus, and the initial fields equal the
background on open neighbourhoods of the shell and vertex. The theorem
holds independently in both interior copies.

**Corollary (exact sufficiently short precontact family).** If the fixed
background also exists smoothly on a regular spacetime neighbourhood of
the initial two-copy section and its timelike shell, the initial family
has an exact Einstein–scalar evolution for a common sufficiently short
time interval. It is `C²` in the parameter in a piecewise spacetime
`C⁴` norm and equals the background on an open shell collar. Keeping
the shell embeddings equal to their background values therefore satisfies
all matching, metric/scalar junction and normal-force equations exactly
on this interval. The time is chosen using the actual perturbed causal
cones and a strict support margin. No particular long background-cone
contact time is asserted.

## Conditional inputs already supplied by the checkpoint

Let `𝒜′` be a smooth concentric annulus away from shell and vertex, with
a smooth radial time-symmetric Einstein–scalar background. Its spatial
metric is positive definite and its scalar has `Φ′≠0` on a neighbourhood
of the closed annulus. Let `Γ≅(Z₂)^4` act by coordinate reflections.
Use the checkpoint's local adjoint classifications and the corrected
(API) arguments for the scalar metric operator and symmetric divergence.

Let `v=(h_v,ψ_v,π_v,ϖ_v)` be smooth, `Γ`-invariant, compactly
supported strictly inside `𝒜′`, and solve the linearized constraints.
Then the following argument produces exact constraint data tangent to
`v`, subject to the standard scaled-cover estimates from Delay Appendix
10 used explicitly below.

## Fixed weighted spaces and linear inverse

Write `d` for a radial boundary defining function and choose

```
q=d²,    w=d^{2a−4} exp(−β/d),    a∈N, β>0.
X_k=(Ḣ^{k+4}_{q,w} ⊕ Ḣ^{k+2}_{q,w})^Γ,
Y_k=(Ḣ^k_{q,w} ⊕ Ḣ^k_{q,w})^Γ,     integer k≥12.
```

The dot denotes closure of compactly supported smooth sections in the
indicated weighted norm. Normalize the momentum density using the fixed
background density, so the range pairing is a fixed-background one-form
pairing. Denote this constraint map by `Ĉ=(H,M/√detγ₀)`.
Rescale the defining function by a fixed small constant if needed so
`B(p,q(p))` stays inside the annulus, as in Delay's equation (10.4).

Write `ρ₀=√detγ₀` for the fixed background density. A momentum
correction is represented as `δπ^{ij}=ρ₀τ^{ij}`, with `τ`
an ordinary symmetric contravariant tensor. The fixed linear operator is
`(P_mτ)_i=−2γ₀,ik D⁰_jτ^{jk}`. Its formal adjoint, for a
one-form potential `Z` and vector `Z♯=γ₀⁻¹Z`, is
`P_m*Z=(L_{Z♯}γ₀)^{ij}`, with both indices raised using `γ₀`.
This states the normalization and pairing explicitly; its kernel is the
rotation Killing space and its (API) is the symmetric-divergence one.

For scalar potential `N` and one-form potential `Z`, define fixed linear
correction operators

```
A_HN=w²q⁴ P_g* N,
A_MZ=w²q² P_m* Z.
```

Convert the second output into the contravariant momentum density by
multiplication by `ρ₀`. These smooth fixed factors do not change any
estimate. The resulting `A(N,Z)` corrects only metric and metric momentum,
leaving scalar and scalar momentum prescribed.

The derivative of `F(ε,u)=w⁻²Ĉ(x₀+εv+Au)` at the origin is
the direct sum

```
L_H=w⁻² P_g w²q⁴ P_g*,
L_M=w⁻² P_m w²q² P_m*.
```

Delay Theorem 3.8 identifies each block as a bounded isomorphism on
the complement of its adjoint kernel. Both kernels have zero invariant
part, and the weights are invariant. Every invariant source is weighted
orthogonal to the full kernel; uniqueness of the inverse implies that
its solution is invariant too. Thus `L:X_k→Y_k` is a bounded
isomorphism. This step does not require classifying KIDs of the perturbed
metric.

## Weighted nonlinear mapping lemma at fixed k

Use a scaled covering by coordinate balls `B_i` of radius comparable
to `q_i=q(p_i)`, with bounded overlap of enlarged balls and comparable
`q,w` on each enlarged ball. Delay's Appendix 10 supplies these scaling
properties for `q=d²` and exponential `w`; weighted derivatives of
`q,w` also obey its uniform multiplier estimates. In dimension four,
the weighted norm is equivalent to the sum of squares of local scaled
norms multiplied by `w_i q_i²`.

Set

```
V_H,i=w_i q_i² ||N||_{H^{k+4}(scaled enlarged ball)},
V_M,i=w_i q_i² ||Z||_{H^{k+2}(scaled enlarged ball)}.
```

Their squared sums are controlled by `||u||²_X` and conversely, with
the usual enlarged-cover constants. Uniform differentiation of the fixed
correction multipliers gives

```
||h||_{H^{k+2}(scaled ball)} ≤ Cw_i V_H,i,
||Π||_{H^{k+1}(scaled ball)} ≤ Cw_i q_i⁻¹ V_M,i.
```

Here `Π` is the physical momentum density correction. In particular a
small `X_k` ball gives a uniformly small physical metric in `C⁰`,
because `H^{k+2}` embeds into `C⁰` on the unit balls and `w` is
bounded. The inverse metric and reciprocal determinant consequently
exist uniformly and are smooth local Sobolev multipliers. The rescaled
background and its derivatives have uniform bounds on the closed
annulus.

After subtracting the background value and its full linearization, the
curvature terms have orders

```
h D²h, (Dh)², h Dh Dγ₀, h² D²γ₀,
```

with smooth inverse-metric coefficient factors. The Hamiltonian momentum
terms have order `Π²`; the scalar momentum term has order `ϖ²`
with prescribed compact support. The scalar gradient remainder has
orders `h²|Dφ₀|²`, `h Dφ₀ Dψ`, and `(Dψ)²`.
The potential `U(φ₀+εψ_v)` is a smooth, compactly supported parameter
contribution after subtracting its background and linear terms. No
derivatives of an unknown scalar correction occur.

There is an especially useful exact coordinate formula for the momentum
constraint with contravariant symmetric weight-one density `Π^{jk}`:

```
M_i(γ,Π,φ,ϖ)
 = −2γ_{ik}∂_jΠ^{jk} −2Π^{jk}∂_jγ_{ki}
   +Π^{jk}∂_iγ_{jk} +ϖ∂_iφ.
```

It follows by expanding the density divergence; the trace-connection
terms cancel. Its nonlinear metric/momentum remainder is exactly of
orders `h DΠ` and `Π Dh`, plus the prescribed scalar term
`ϖ Dψ`. Thus there is no hidden loss of momentum derivatives and
no unknown inverse-metric coefficient in this density formulation.

Since `k>2`, the ordinary unit-ball `H^k` algebra and composition
estimates apply. Each nonlinear remainder has local norm bounded by

```
||Q_H||_{H^k(scaled)}+||Q_M||_{H^k(scaled)}
 ≤ C w_i² q_i⁻² (V_H,i+V_M,i+|ε|V_seed,i)².
```

The seed constants include all four ADM sectors of `v` and are
supported on a compact subannulus; only finitely many balls of a locally
finite scaled cover meet that support, so their squared sum is finite.
Lower-order background terms satisfy this bound after
enlarging its constant, because `q` is bounded. Normalize by `w⁻²`
and apply the range norm factor `w_i q_i²`. The local range contribution
is at most `Cw_i(V_H,i+V_M,i+|ε|V_seed,i)²`. Hence

```
||w⁻²Q||_Y ≤ C(||u||_X+|ε|)²,
||w⁻²(Q(ε,u)−Q(ε,z))||_Y
 ≤ C(|ε|+||u||_X+||z||_X)||u−z||_X.
```

For the sum, use
`ΣV_i⁴ ≤ (ΣV_i²)²` and `sup_i w_i<∞`.
All higher Fréchet derivatives have the same bounded local products;
additional metric factors introduce additional bounded factors of `w_i`
and local potential norms. Matrix inverse and determinant composition are
smooth on the uniformly positive metric neighbourhood. The prescribed
scalar potential term is smooth in `ε` on its compact support. These
observations establish `F:R×X_k→Y_k` as a smooth map on a small
neighbourhood of the origin, in the fixed finite Sobolev topology.

For the dotted range explicitly, first define the map on compactly
supported smooth potentials. There the physical corrections and the
constraint differences have compact support, so their normalized values
belong to the dotted range. The estimates and derivative bounds extend
this map smoothly to the dotted domain by density, retaining its values
in the closed dotted range `Y_k`.

## Physical zero extension at fixed finite regularity

The same cover controls the unweighted physical norms. Taking the highest
derivative in the scaled balls gives

```
||h||_{H^{k+2}(B_i)} ≤ Cw_i q_i⁻k V_H,i,
||Π||_{H^{k+1}(B_i)} ≤ Cw_i q_i⁻k V_M,i.
```

The lower derivatives obey the same bound after adjusting a constant.
For every fixed finite `k`, `sup_i w_i q_i⁻k<∞` because the
exponential dominates every power at the boundary. Summing squares gives
continuous maps

```
A_H:Ḣ^{k+4}_{q,w}→H₀^{k+2}(𝒜′),
A_M:Ḣ^{k+2}_{q,w}→H₀^{k+1}(𝒜′).
```

To obtain the subscript zero, approximate each potential in its defining
weighted closure by compactly supported smooth potentials, averaging over
`Γ` if necessary. Their physical corrections are compactly supported
smooth tensors. The unweighted estimates show convergence in the physical
Sobolev spaces; the limit belongs to their compact-support closures.
Consequently extension by zero gives global `H^{k+2}` metric and
`H^{k+1}` momentum corrections without a distributional boundary layer.
The corrections vanish on fixed open neighbourhoods of the physical
shell and vertex because the entire annulus closure stays away from both.

## Exact fixed-regularity conclusion

The Banach implicit function theorem now yields, for this fixed `k`,
a small interval and a smooth map `u(ε)∈X_k` with `F(ε,u(ε))=0`.
Since `Pv=0`, differentiating at zero gives `u(0)=u′(0)=0`.
Thus

```
x(ε)=x₀+εv+A u(ε)
```

solves the **exact** initial constraints, is tangent to `v`, and has
parameter-smooth dependence in physical `H^{k+2}×H^{k+1}`.
The prescribed scalar and scalar momentum are smooth compact-support
parameter polynomials. Smallness keeps the metric positive definite.
The data equal the background near the shell, so its exact initial
junction conditions can be retained unchanged.

The actual family selects a particular second derivative given by

```
x″(0)=−A L⁻¹[w⁻² D²Ĉ|_{x₀}[v,v]].
```

This is a bounded weighted-operator formula, with the same fixed-density
constraint normalization used throughout. Its scalar and scalar-momentum
second derivatives are zero because those sectors were prescribed
linearly in `ε`. The formula solves the full second-order initial
constraint equation but may choose a different particular correction and
gauge representation from the explicit correction in the frozen localized
capsule. No equality with that specific frozen `w` is asserted.

In four dimensions `H^{k+1}` embeds into `C^{r,α}` whenever
`r+α<k−1`; choose large finite `k` to meet any specified finite
classical regularity and parameter-differentiability requirement. For
example the fixed choice `k≥12` supplies ample finite regularity beyond `C⁴`.
This proves a fixed-regularity implication; it does not manufacture a
single `C∞` spatial family by intersecting infinitely many IFT intervals.

## Concrete commuting quadrupole consequence

Take `A=B=diag(1,−1,0,0)` and smooth radial profiles `G,K`
compactly supported strictly inside the selected annulus. The inherited
master initial data are

```
q=G(s)(x₁²−x₂²),    p=K(s)(x₁²−x₂²),    s=|x|².
```

Their inherited ADM reconstruction solves the linearized constraints and
is invariant under each independent coordinate reflection, in all four
ADM sectors. The theorem therefore gives an exact initial-data family
tangent to this concrete quadrupole, and the corollary below gives its
exact evolution on a uniformly sufficiently short protected slab.
Identical sequences on the two copies realize the even family. To
realize opposite first-order copy signs, use `x_+(ε)` and
`x_-(ε)=x_+(-ε)`; the shell collar is background in either case.

More generally, commuting real symmetric quadrupole matrices can be
simultaneously orthogonally diagonalized, and the same argument applies
using the coordinate-reflection group in that fixed basis. No additional
radial frozen-profile condition is required for this sufficient symmetry
construction.

This consequence directly advances the admissible-quadrupole target, but
it is an analytic existence construction. It supplies neither a closed
formula for the evolved metric nor a numerical spacetime evolution, and
it does not construct evolution through shell contact.

## Exact evolution on a uniformly short protected slab

Let the initial closed annulus be contained in
`R_-≤R≤R_+`, with `0<R_-<R_+<b`. These radii enclose the
union of all initial corrections, not just the first-order packet.
Choose `B` strictly between `R_+` and `b`, and work inside a smooth
regular neighbourhood of the initial ball of radius `B`. The background
is the fixed actual solution on this neighbourhood and near the shell;
no regular branch extension is inferred from the initial constraint proof.

Use a background-adapted generalized harmonic reduction with fixed smooth
reference connection equal to the connection of the actual background
on the controlled patch, so the background satisfies the same reduced
equations. Extend that reference connection only for the auxiliary
construction below. Choose initial lapse and shift equal to the
background. Set their time derivatives by the local harmonic Cauchy
gauge construction, so the harmonic vector vanishes initially. This
construction is local in the canonical fields and their first spatial
derivatives, so its parameter-dependent gauge data also agree with the
background outside the initial closed annulus. The exact initial
Hamiltonian and momentum constraints set the normal derivative of the
gauge vector to zero. The scalar equation and contracted Bianchi identity
give the homogeneous gauge propagation equation.

The reduced Einstein–scalar system is a smooth quasilinear vector-bundle
wave system on the regular Lorentzian neighbourhood. Its initial metric
and scalar fields have `H^{s}` regularity and velocities have
`H^{s−1}`, with `s=k+2≥14`. The high-regularity local Cauchy
theorem gives a common short time for a bounded sufficiently small
parameter interval. Positivity of the initial spatial metric and the
Lorentzian lapse remain uniform by shrinking this interval and time.

The smooth initial parameter curve has two differentiable solution
variations in lower Sobolev orders. Differentiating the reduced wave
equations first gives a linear hyperbolic system with the actual solution
metric as principal metric. A source of order `g₁∂²g` is in
`H^{s−2}`, giving a first variation of order `H^{s−1}` by the
energy estimate. The quadratic second-variation source
`g₁∂²g₁` is in `H^{s−3}`, giving a second variation of
order `H^{s−2}`. These estimates describe the formal variation
values; they alone do not assert differentiability of the solution map
at those highest variation orders.

For actual parameter differentiability, compare finite difference
quotients with the variation equations. A coefficient difference
multiplying `∂²g₁` lies only in `H^{s−3}` in this comparison.
Further comparing first variations and their quotients introduces
`∂²` of a first-variation difference, controlled in a still lower
source norm. The same linear hyperbolic energy estimates give quotient
convergence and continuity in lower orders. We deliberately allow four
total Sobolev losses and record only `C²` parameter dependence in
fields `H^{s−4}` and velocities `H^{s−5}`, on the common short
interval. This includes additional comparison margin beyond the formal
variation-value estimates and uses no sharp top-norm smooth-flow claim.

With `s≥14`, taking four more time/spatial derivatives leaves
`H^{s−8}` with `s−8≥6`, above the four-dimensional continuity threshold
`4/2=2`. Thus the parameter derivatives have piecewise spacetime
`C⁴` regularity on a compact slab. This is finite classical regularity,
not an assertion of spatial analyticity or `C∞` regularity.

For an explicit local patching construction, embed the radius-`B` ball
in a coordinate chart of a sufficiently large compact four-torus.
Extend the reduced Cauchy data and reference connection smoothly outside
that ball, keeping the initial metric uniformly Lorentzian and spatial
metric positive. The parameter-dependent corrections have support inside
the annulus, so this extension can be common to all parameters outside
it. Compactness gives finite global Sobolev norms without imposing decay
on the background. Its constraints outside the controlled ball need not
be used. The background auxiliary reduced solution agrees with the
actual background on the domain of dependence of the controlled ball.

Choose a fixed coordinate speed `c>1` bounding the spatial speeds of
all null directions of the actual solution metrics in the controlled
neighbourhood, for both time directions and all sufficiently small
parameters. Such a finite common bound follows from the regular
background and the continuous dependence just established. In the
conformal ambient coordinates the background speed is one, but the
perturbed cones are not assumed to be identical to the background ones.
Choose a positive time `τ` below the Cauchy existence time and with

```
2cτ < B−R_+,
cτ < R_-/2,
```

and, if only a particular controlled background neighbourhood is used,
shrink `τ` further so the entire controlled causal region remains there.
The constraints and gauge propagate in `R<B−c|T|`, because their
initial values vanish in the radius-`B` ball. The difference of each
auxiliary reduced solution from the auxiliary background solves a linear
homogeneous hyperbolic system with coefficients built from the two
solutions and compactly supported initial differences. Finite propagation
therefore confines that difference to

```
max(0,R_-−c|T|) ≤ R ≤ R_++c|T|.
```

There is consequently a nonempty common annular overlap between
`R_++cτ` and `B−cτ`. On it the auxiliary solution, auxiliary
background and actual background agree. Choose one fixed patching radius
in that overlap. Use the evolved solution inward and the exact actual
background outward. Their agreement on an open overlap gives an exact
Einstein–scalar solution throughout each physical interior, without a
new interface or a dependence on the arbitrary auxiliary extension.

The shell lies in the protected background region for this sufficiently
short time; its background trajectory is `R_b(T)=√(b²+T²)≥b`.
Both copies therefore retain the exact common induced fields, outward
extrinsic curvatures and scalar normal derivatives there. Keeping the
shell embeddings fixed to the background trajectory solves every exact
matching and Israel/scalar junction equation, and the exact normal-force
identity then applies. The protected inner neighbourhood likewise remains
background, by the second strict inequality above. This proves the stated
short precontact corollary.

This argument is a local high-regularity Cauchy construction. It claims
neither a long slab up to the background-cone contact bound, nor any
nonlinear evolution once a perturbation reaches the physical shell.

## Completed review and remaining scope

The scaled-cover nonlinear mapping and physical Sobolev zero-extension
arguments have been checked independently against the fixed-density
constraint conventions and Delay's weighted closure/isomorphism
hypotheses. The substantive step beyond algebraic surjectivity is the
fixed-`k` Banach estimate and its unweighted extension, not the formal
series induction. Their accepted finite-regularity conclusion does not
depend on completing an infinite smoothness bootstrap.

Spatial `C∞` regularity on one common parameter interval would require
the additional weighted nonlinear regularity argument. No convergence
claim about an arbitrary formal coefficient sequence is made. The
initial implicit-function solution and the short exact Cauchy evolution
are actual solutions in the stated finite classes; their extension
through shell contact remains a separate unproved problem.

Primary linear source: E. Delay, *Smooth compactly supported solutions of
some underdetermined elliptic PDE, with gluing applications*,
[arXiv:1003.0535v4](https://arxiv.org/html/1003.0535v4), Theorem 3.8,
operator-specific (API) arguments, and Appendix 10.

Primary local hyperbolic sources:

- T. J. R. Hughes, T. Kato and J. E. Marsden, *Well-posed quasilinear
  second-order hyperbolic systems with applications to nonlinear
  elastodynamics and general relativity* (1977), Theorem III in Section
  3.1, [author PDF](https://www.cds.caltech.edu/~marsden/bib/1977/02-HuKaMa1977/HuKaMa1977.pdf).
  Its common-neighbourhood short-time Sobolev existence and continuous
  dependence apply to the smooth uniformly hyperbolic reduced system.
  The parameter argument above is derived by differentiating, comparing
  quotients and applying energy estimates with four reserved Sobolev
  losses; top-norm smoothness of the solution map is not assumed.
- H. Ringström, *Future stability of the Einstein–non-linear scalar field
  system*, Invent. Math. 173 (2008), 123–208,
  [DOI](https://doi.org/10.1007/s00222-008-0117-y),
  [author PDF](https://people.kth.se/~hansr/invgnls.pdf).
  The local Einstein–scalar reduction in Section 1.3, causality arguments
  in Sections 3.2–3.3 and homogeneous gauge propagation in Section 4
  supply the machinery used here. Section 4.1 explicitly permits the
  constraints to be imposed on only a subset of the initial slice.
  The paper's global stability conclusion and its potential assumptions
  are not invoked.

**Independent cross-review.** The literature-review agent checked the
fixed-density momentum formula, Delay's dotted-space isomorphism,
Appendix 10 scaling/multiplier hypotheses, all four seed sectors, the
dimension-four correction factors, and the unweighted zero-extension
exponents. It found no substantive fixed-finite-regularity gap, requesting
the defining-function rescaling explicitly added above. This is internal
agent review, not external peer review.

The same independent reviewer checked the short precontact corollary:
actual perturbed cone bounds, the union of correction supports, the
background reference connection, local homogeneous gauge propagation,
uniform existence time, the open-overlap background patch, and the four
reserved Sobolev parameter losses. It accepted the argument after the
compact auxiliary extension was specified above. The primary Cauchy and
Einstein–scalar source passages were inspected for those local uses.
