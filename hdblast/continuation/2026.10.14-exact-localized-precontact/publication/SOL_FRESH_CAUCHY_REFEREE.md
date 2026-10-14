# Fresh internal referee report: short precontact Cauchy construction

14 October 2026. Independent review of the exact short precontact corollary in
`publish-exact-precontact-20261014/research/HDBLAST_EXACT_LOCALIZED_PRECONTACT_20261014/EXACT_LOCALIZED_INITIAL_DATA_AND_PRECONTACT.md`,
especially lines 308–427, at the publication state associated with PR 49,
head `29a3d374f613393232660cef4ef1aebc3b5982d2`.

## Judgment

**No substantive gap found in the stated fixed-finite-regularity, sufficiently
short precontact corollary, conditional on its exact initial-data theorem.**
The four reserved Sobolev losses suffice for actual twice continuously
differentiable parameter dependence. This conclusion does not rely on treating
formal variations as derivatives of an already smooth nonlinear flow. The
compact auxiliary construction, local recovery of the full Einstein equations,
actual-cone support bound, and open-overlap patch are also consistent.

This is fresh internal agent review, not external peer review. I checked the
primary source passages and independently reconstructed the estimates; earlier
agent acceptance was not used as a premise. A separately delegated check of the
quotient identities agreed with the reconstruction below. I changed no tracked
publication files.

## Regularity and the local existence input

Let `s=k+2≥14`, and denote the complete reduced metric/scalar unknown by `U`.
The supplied initial curve is smooth into `H^s×H^(s−1)`. Conversion from the ADM
momentum densities to extrinsic curvature and scalar velocity is a smooth
Sobolev composition on the positive-metric neighborhood. The local harmonic
Cauchy gauge prescription uses first spatial derivatives of the metric, so its
metric velocities have exactly the available `H^(s−1)` regularity. It agrees
with the background wherever the canonical data do. No spatial `C∞` claim is
needed for these operations.

Hughes–Kato–Marsden, §3.1, assumptions (a1)–(a3), Theorem III, supplies a common
time for a neighborhood of data in `H^(r+1)×H^r`, together with uniqueness and
continuity uniformly in time. Here `r=s−1≥13` and spatial dimension is four.
Both its general threshold `r>n/2+1` and its weaker threshold for principal
coefficients independent of field derivatives are comfortably satisfied. The
reduced Einstein–scalar principal coefficients depend on the metric, while
the scalar potential and reference-connection terms are lower order. A compact
regular coefficient neighborhood gives bounded derivatives and uniform
hyperbolicity; the smooth potential need only be defined on a neighborhood of
the controlled scalar range. In the adapted background coordinates the initial
time direction and slices have the required timelike/spacelike character, and
this persists on a common smaller slab.

Theorem III is stated on `R^n`, not literally on `T^4`. Passing its local energy
and iteration argument to a finite compact cover is legitimate: partitions of
unity introduce lower-order terms, and the same compact Sobolev product and
energy estimates apply. Equivalently, one can use local extensions and finite
propagation. The torus removes the nondecay and zero-source integrability
problem. The packet is appropriately using this local adaptation rather than
pretending the original theorem was stated on the torus. HKM explicitly warns
that its top-space continuous solution map need not have stronger smoothness;
the parameter claim therefore requires the separate quotient proof.

“Extend the reduced Cauchy data smoothly outside that ball” is read as a smooth
exterior extension of the background fields, retaining the finite-Sobolev
corrections inside. Since those corrections vanish on a neighborhood of the
extension region, this does not assert smoothness of the full perturbed data.
A positive spatial metric can be extended using its convex positive cone;
lapse/shift variables give a Lorentzian extension with a uniform time slicing.
Connections admit extension because their differences are tensors.

## Actual parameter derivatives: quotient estimates

For clarity, write `E^m=C_T H^m ∩ C_T^1 H^(m−1)` on a common compact time
interval. On the compact auxiliary space, suppress the fixed smooth coordinate
and reference-connection dependence and write the reduced system as

```
P(U)=a(U)∂²U−F(U,∂U)=0,
L_U=DP_U,                 B_U=D²P_U.
```

The notation includes all spacetime second derivatives. Uniform `E^s` bounds
and strict hyperbolicity give the linear energy estimate, for the lower orders
used below,

```
||w||_(E^m) ≤ C( ||Cauchy data(w)||_(H^m×H^(m−1))
                   + ||L_U w||_(L_T^1 H^(m−1)) ).
```

Lower-order terms such as `Da(U)[w]∂²U` belong to the source norm and can be
absorbed by Gronwall. Their coefficients have `H^(s−2)` regularity. The wave
equation, with uniformly invertible time principal coefficient, also gives
`∂_T²w∈C_T H^(m−2)` at these orders. Thus the estimates do not silently replace
time derivatives by spatial derivatives without controlling them.

First define candidate variations by linear equations, without assuming
differentiability of the nonlinear flow:

```
L_U v=0,
L_U z=−B_U(v,v),
```

with the first and second differentiated initial data. They satisfy uniform
`v∈E^(s−1)` and `z∈E^(s−2)` bounds. The term `v∂²v` has order `H^(s−3)`;
the remaining quadratic terms have at least that regularity.

Fix a parameter `θ`, put `d=U_(θ+h)−U_θ`, and `q=d/h`. The exact difference
equation, with the principal operator taken from one actual solution, gives

```
||d||_(E^(s−1)) ≤ C|h|,       ||q||_(E^(s−1)) ≤ C.
```

This uses the smooth initial curve and a source of the form
`[a(U_(θ+h))−a(U_θ)]∂²U_θ`; it does not divide a merely continuous top-norm
bound by `h`. Exact Taylor expansion of the differential expression then gives

```
L_U(q−v)=−h ∫₀¹(1−t) B_(U+td)(q,q) dt.
```

Its source is `O(|h|)` in `H^(s−3)`, and the initial quotient error is also
`O(|h|)`. Therefore

```
||q−v||_(E^(s−2)) ≤ C|h|.
```

Let `ṽ=v_(θ+h)`. Comparing the candidate first-variation equations gives the
exact identity

```
L_U(ṽ−v)=−h ∫₀¹ B_(U+td)(q,ṽ) dt,
```

and hence `||ṽ−v||_(E^(s−2))≤C|h|`. For `p=(ṽ−v)/h`, subtraction of the
candidate second-variation equation yields a residual equal to

```
−B_U(q−v,ṽ) − B_U(v,ṽ−v)
− ∫₀¹ [B_(U+td)−B_U](q,ṽ) dt.
```

The first two terms are `O(|h|)` in `H^(s−4)`: two derivatives of an
`E^(s−2)` comparison difference have exactly that order. The last term is
`O(|h|)` in `H^(s−3)`, by the smooth coefficients and the `E^(s−1)` difference
bound. Consequently

```
||p−z||_(E^(s−3)) ≤ C|h|.
```

Finally, subtracting the second-variation equations gives

```
L_U(z_(θ+h)−z_θ)
 = −(L_(U_(θ+h))−L_U)z_(θ+h)
   −[B_(U_(θ+h))(ṽ,ṽ)−B_U(v,v)].
```

The source is `O(|h|)` in `H^(s−4)`, so
`||z_(θ+h)−z_θ||_(E^(s−3))≤C|h|`. The corresponding initial errors follow
from the smooth fixed-space initial curve. All estimates are uniform on a
smaller closed parameter interval. At the minimal `s=14`, the lowest source
space in this chain is `H^10`, well above the dimension-four multiplication
threshold. In fact this reconstruction reaches actual `C²` into `E^(s−3)`;
the manuscript claims only `E^(s−4)`, which has additional margin. No sharp
topology of the nonlinear flow is inferred.

For the stated spacetime conclusion, the wave equation recursively expresses
the next time derivatives as smooth Sobolev products. Actual `C²` into
`E^(s−4)` therefore implies actual `C²` into
`C_T^j H^(s−4−j)` for `0≤j≤4`. Any total of four time/spatial derivatives
leaves at least `H^(s−8)`, with `s−8≥6>4/2`. Sobolev embedding on compact
patches gives parameter `C²` dependence in the piecewise spacetime `C⁴` norm,
including convergence of the differentiated quotients in that norm.

## Gauge propagation and causal patching

The fixed reference connection equals the connection of the actual background
on the controlled patch. Thus the actual background has zero harmonic gauge
vector and solves the same reduced system there. The harmonic-compatible
Cauchy construction sets that vector to zero initially; the exact Hamiltonian
and momentum constraints set its normal derivative to zero in the radius-`B`
ball. The scalar equation makes its stress tensor divergence free, and the
contracted Bianchi identity gives a homogeneous wave equation for the gauge
vector. The high finite regularity is more than sufficient for this identity
and for its local linear energy uniqueness argument.

Ringström §4, equations (55)–(57), exhibits this mechanism. Section 4.1
explicitly allows constraints on only an initial subset `S`; its conclusion is
the full Einstein–scalar equations on `D(S)` defined by the solution metric.
The auxiliary extension may consequently violate constraints outside the ball.
No gauge or constraint claim outside the physical domain is required.

A common bound `c` for actual causal coordinate speeds follows from uniform
Lorentzian regularity and the common-time Sobolev continuity. The bound is used
inside the controlled ball; arbitrary auxiliary behavior outside cannot enter
`R<B−c|T|` during the chosen interval. The exact constraints/gauge therefore
propagate throughout that region. The difference between a perturbed reduced
solution and the auxiliary background satisfies a homogeneous linear equation:
coefficient differences in the principal part become terms linear in the
difference, multiplying derivatives of the other solution. Finite propagation
places its support inside

```
R_-−c|T| ≤ R ≤ R_++c|T|.
```

The inequalities `2cτ<B−R_+` and `cτ<R_-/2` give respectively a nonempty fixed
annular overlap and an inner protected neighborhood. Choose a radius strictly
between `R_++cτ` and `B−cτ`. A fixed open radial neighborhood of that radius
belongs to both the exact physical domain and the region where the auxiliary
solution agrees with the actual background, uniformly for all parameters and
all times in the slab. Patching across that overlap creates no interface
distribution, and its result is independent of the arbitrary extension. Using
a fixed radius also preserves the parameter norm statement.

Outside that radius the physical field is the actual background. The shell
trajectory is therefore surrounded by an open background collar in both
copies. Every induced field and every derivative entering the metric/scalar
junction and force equations is exactly its background value. Retaining the
background shell embeddings satisfies those equations exactly. This step
requires only the assumed exact background shell solution; it is not a
timelike-boundary existence theorem. The opposite-sign parameter choice on the
second copy has the same protected collar.

## Source precision and scope

One minor citation label should be corrected separately: Ringström Theorem 6
is a zero-data domain-of-dependence **uniqueness** theorem, not an existence
theorem. The ledger's phrase “linear-wave existence/uniqueness statement” is
overbroad. Its actual use here is uniqueness, so this does not weaken the
corollary. Its smooth statement is being implemented through the same local
energy uniqueness argument at the available high finite regularity.

Ringström's special assumptions on `V(0)`, `V′(0)` and `V″(0)` are hypotheses
for that paper's stability setting. None is needed for the local reduction,
Bianchi identity, or the compact finite-Sobolev Cauchy construction with an
arbitrary smooth potential on the controlled range. No global stability
conclusion is transferred.

The reviewed conclusion remains conditional on the supplied exact initial
family and on the actual smooth regular background neighborhood. It supplies
a sufficiently short exact two-copy family before perturbations can reach the
shell. It supplies neither one common spatial `C∞` parameter family, a long
background-cone contact time, through-contact evolution, global stability,
Big Bang conclusions, nor a novelty claim.

Primary local snapshots inspected:

- `cauchy-sources/HKM1977_PRIMARY.md`, §3.1, particularly lines 378–423.
- `cauchy-sources/RINGSTROM_PRIMARY.md`, §1.3, §3.1, Theorem 6 in §3.3,
  §4 equations (55)–(57), and §4.1/its domain-of-dependence conclusion.
- The publication packet's `CAUCHY_SOURCE_LEDGER.md`.
