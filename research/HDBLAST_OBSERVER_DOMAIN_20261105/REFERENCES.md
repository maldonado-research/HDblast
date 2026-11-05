# Primary method references for a common observer domain

5 November 2026 UTC. These three author-written mathematical sources are a
focused comparison for embedding stability, normal-tube obstructions, and ODE
continuation. The indicated passages were read; this is not a complete proof
audit or an exhaustive originality search. The first and third sources are
monographs presenting classical methods; the second is a research preprint.

1. **Morris W. Hirsch**, *Differential Topology*, Graduate Texts in
   Mathematics **33**, Springer, 1976,
   [doi:10.1007/978-1-4684-9449-5](https://doi.org/10.1007/978-1-4684-9449-5).
   Read Chapter 2, §1, Lemma 1.3 and Theorem 1.4, printed pp. 36–38,
   in the [book-text mirror](https://www.academia.edu/143469285/Hirsch_DifferentialTopology);
   the [publisher page](https://link.springer.com/book/10.1007/978-1-4684-9449-5)
   verifies the edition. Theorem 1.4 gives openness of embeddings in the
   strong `C^r` topology, `r>=1`; on compact source manifolds the preceding
   topology discussion identifies strong and weak topologies. Its proof
   combines local immersion/injectivity with separation of distant points.
   Also read Chapter 5, §1, printed pp. 122–125, for degree as a signed
   count and homotopy invariance. These support the qualitative topology;
   they do not provide the observer-flow constants. The explicit sphere
   chord estimate below is an elementary specialization proved here.

2. **O. C. Durumeric**, *Thickness Formula and C¹-Compactness for C¹,¹
   Riemannian Submanifolds*,
   [arXiv:math/0204050v1](https://arxiv.org/abs/math/0204050v1),
   3 April 2002. Read the [primary PDF](https://arxiv.org/pdf/math/0204050v1),
   Theorem 1 on printed p. 2 and §2, Definitions 2, 5, 6 on printed
   pp. 5–6. The normal injectivity radius is defined by global one-to-one
   behavior of the normal exponential map. For compact boundaryless
   `C^{1,1}` submanifolds of a complete connected smooth **Riemannian**
   manifold, Theorem 1 expresses it as the minimum of geometric focal
   distance and half the minimal double-critical distance. This precisely
   exhibits a local/focal obstruction and a global normal-intersection
   obstruction. It is a comparison source for the distinction; its metric
   balls, minimizing distances, completeness assumptions, and thickness
   formula are not imported into the Lorentzian shell problem. A positive
   cut-metric determinant or a nondegenerate label differential alone
   does not establish global noncrossing.

3. **Gerald Teschl**, *Ordinary Differential Equations and Dynamical
   Systems*, Graduate Studies in Mathematics **140**, American
   Mathematical Society, 2012. Read the authorized
   [author's preliminary PDF](https://www.mat.univie.ac.at/~gerald/ftp/book-ode/ode.pdf),
   Chapter 2, Theorems 2.10–2.11, printed pp. 46–48, and Lemma 2.14,
   Corollaries 2.15–2.16, printed pp. 52–53. The
   [author's book page](https://www.mat.univie.ac.at/~gerald/ftp/book-ode/)
   verifies publication and PDF authorization. The first passages justify
   higher parameter/initial-data derivatives and variational equations
   under the stated `C^k` hypotheses. The latter give extension when the
   solution remains in a compact subset of the open ODE domain. The
   relevant use is local existence followed by an a priori confinement
   estimate and continuation, followed by parameter differentiation on
   the resulting common interval. This source supplies no shell PDE
   family, physical embedding, or numerical coefficient bound.

The following are explicit specializations of standard methods, rather than
additional external theorems.

## A quantitative cylinder embedding comparison

Let `T_q` be a coordinate-time duration and set the **proper duration**
`S=rho_b*T_q`. On `C=[0,S]×S³`, use the fixed reference Riemannian metric
`ds²/rho_b²+g_round` and the auxiliary ambient realization

```text
F_0(s,Omega)=(-s/rho_b,Omega) in R×R⁴.
```

This auxiliary Euclidean geometry only measures maps; it is not the shell
Lorentzian metric or the physical bulk-shell embedding. If `theta` is the
round spherical distance between two unit vectors, then `0<=theta<=pi` and

```text
|Omega-Omega'|=2 sin(theta/2) >= (2/pi) theta,
|F_0(p)-F_0(q)| >= (2/pi) d_C(p,q).
```

The first inequality follows from concavity of sine on `[0,pi/2]`; the
second follows by combining it with the unchanged time-coordinate
difference in the Euclidean product norm. Suppose the complete observer
map exists on `C` and its ambient representation is `F_e=F_0+E_e`, with
`sup_C ||D E_e||<=eta<2/pi` in these fixed norms. Integrating `D E_e`
along a minimizing product geodesic gives

```text
|F_e(p)-F_e(q)| >= (2/pi-eta) d_C(p,q).
```

Thus `F_e` is globally injective. Since `D F_0` is an isometry on each
tangent space, `|D F_e v|>=(1-eta)|v|` gives immersion as well. Compactness
then gives an embedding, including the time-boundary cuts. If
`sup ||partial_e D F_e||<=B` on the admitted parameter/cylinder box, the
guard `a B<2/pi` is sufficient. This coefficient depends on the stated
normalization; changing reference scales changes the comparison constants.
The argument proves global noncrossing after existence and the derivative
estimate are established. Spacelike-cut positivity needs its own metric
comparison.

An alternative is to prove a uniform strict coordinate-clock derivative
bound and uniform endpoint coverage, invert the clock label by label,
and examine the resulting angular map `G_t:S³→S³` on a common coordinate
interval. If each `G_t` is a local diffeomorphism homotopic to the identity,
compactness makes its fibers finite, its image open and closed, and its
orientation sign constant. Its degree is one, so the signed fiber count
has one member; `G_t` is therefore a diffeomorphism. Local invertibility
without the global argument is insufficient. The inverse-clock derivative
bounds and endpoint coverage must be proved before using this route.

## Continuation from a fixed coefficient tube

Fix an open regular ODE domain and a closed radius-`R` coefficient tube
around the known background states `X_0(s,Omega)`, `0<=s<=S`, contained
compactly in that domain. Use a fixed finite-dimensional representation
and norm. Assume local existence/uniqueness and, throughout the entire
fixed tube and `|e|<=a`, bounds

```text
||D_X F_e||<=L,   ||partial_e F_e||<=P,
||X_e(0,Omega)-X_0(0,Omega)||<=B_0 |e|.
```

The field bounds are hypotheses about coefficients at candidate states,
not hypotheses about an already existing full perturbed flow. Before a
first exit, Gronwall gives

```text
||X_e(s,Omega)-X_0(s,Omega)||
 <= |e|[exp(Ls) B_0 + phi_L(s) P],
phi_L(s)=(exp(Ls)-1)/L if L>0, and phi_0(s)=s.
```

The strict guard `a[exp(LS)B_0+phi_L(S)P]<R` prevents a first exit.
Compact-set continuation then excludes a maximal time at or before `S`.
The tube includes all labels and is inside an open regular domain, so the
local existence and continuation argument can be uniform in `Omega`.
Only after this step does parameter-flow differentiation on the common
closed cylinder become available. Higher variational estimates require
the corresponding uniform coefficient derivatives on that fixed tube.
Teschl's quoted smooth-dependence theorem assumes joint `C^k` regularity;
if only continuous clock dependence with `C^2` parameter/state dependence
is used, the corresponding integral-equation argument must be supplied
rather than claiming it is the literal stated theorem.

## Retained source interface and remaining assumptions

The retained [nonlinear bridge][bridge], §§1–3, already assumes a
conditional exact shell family and an unevaluated admissible parameter
interval. Its past-domain application requires the reflected exact
constraint curve and complete corresponding gauge data. Its arithmetic
slab margins do not certify the evolution, inverse, regularity, or
parameter-radius premises. Time reversal is tensor pullback, including
the mixed metric components and gauge jets; negating every coordinate
time derivative is invalid.

The bridge §3 anchors observers to a specified common **physical initial
shell section**, with its round labels transported. A different physical
section can change the observable. Its smaller-common-tube statement is
qualitative and family dependent. The retained `OBSERVER_FLOW_BUDGET.md`,
§1, assumes the whole admitted compact proper-clock tube, spacelike cuts,
noncrossing, and valid label maps; its derivative bounds do not establish
those domain assumptions. The connecting interval from the initial
section to the target cuts must be included, not only the short zero
bracket.

Both constructions need the complete pulled-back shell fields
`h_e=iota_e^*g_e` and `varphi_e=iota_e^*phi_e`, and the full initial section,
unit normal, label derivatives, and parameter derivatives. Restricting bulk
coordinates to the background shell omits physical-embedding terms.
Observer-budget §3 also requires fixed scaled norms and the second
derivative costs of chart transitions; its auxiliary embedding is distinct
from `iota_e`. Joint diffeomorphism transport of fields and physical
anchors preserves the chosen observable; changing its anchors does not.

These methods can turn effective fixed-tube coefficient and initial-data
inputs into a conditional common-domain statement. They do not evaluate
those inputs, prove a new Einstein–scalar PDE theorem, establish a
numerical amplitude radius, or supply a completion-uniform domain.

[bridge]: https://github.com/maldonado-research/HDblast/blob/fe1d84d420eefcad60156bdfeff1c3cbde9f66a5/research/HDBLAST_CHARACTERISTIC_ODD_TRACE_20261104/nonlinear/NONLINEAR_BRIDGE.md
