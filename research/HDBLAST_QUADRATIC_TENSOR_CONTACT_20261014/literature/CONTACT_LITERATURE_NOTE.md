# Focused literature update for sourced tensor evolution at shell contact

Prepared 14 October 2026 (America/Los_Angeles). This is an internal research aid, not external peer review. The current model is five-dimensional Einstein gravity with a canonical bulk scalar, a finite timelike shell with scalar-dependent tension, and two bulk copies. The exact localized initial-data/short-precontact checkpoint and material junction supplement are the controlling model definitions.

The useful new literature concerns **causal Robin wave propagation** and **boundary constraint propagation**. It supplies building blocks and design ideas; none of the sources checked here directly proves nonlinear contact evolution for the registered Einstein–scalar shell.

## What is newly useful

| Primary source and verified scope | What to use in the next calculation | What does not follow |
|---|---|---|
| Costeri, Dappiaggi and Juárez-Aubry, [arXiv:2609.24465v1](https://arxiv.org/html/2609.24465v1), 21 September 2026. Proposition 3.5, Theorem 3.6 and Appendix B. Fixed-background Klein–Gordon equation on globally hyperbolic spacetimes with timelike boundary, in every dimension \(d\ge2\). | A causal Green-operator construction for a scalar wave with homogeneous Robin data, including Neumann as \(\kappa=0\). Appendix B uses a weak form, the epsilon trace inequality, a **bulk** coercivity shift, Galerkin approximation and time/space regularity. These are useful checks for a forced tensor amplitude after its exact scalar reduction has been established. | The operator is fixed-background scalar Klein–Gordon, not the coupled Einstein–scalar tensor/constraint/gauge system. The stated presentation takes real constant \(\kappa\) and an **inward** normal. Smooth variable-\(\kappa\) extensions are discussed, rather than furnished as a separate model-specific theorem. Positivity of a quantum ground state, boundary convexity and finiteness of reflected rays are additional questions, not prerequisites for the classical causal-existence statement. |
| Stevens and Valiente Kroon, [arXiv:2604.24931v1](https://arxiv.org/html/2604.24931v1), 27 April 2026. Sections 2.1, 5.1–5.2, 6 and 7. Four-dimensional **vacuum** conformal Einstein equations at a finite timelike boundary. | Analyze the boundary principal symbol and choose a boundary-adapted reduction so the subsidiary constraint system has no incoming transverse modes. The evolution boundary law relates incoming/outgoing Weyl components; fully reflective data form a carefully selected compatible subclass. | The boundary must be ruled by conformal geodesics and satisfy the restrictions of that gauge. Proposition 2 assumes all-order corner compatibility; the paper explicitly does not solve that corner problem. Scalar matter, dimension five and an Israel/scalar interface are absent. Its geometric boundary restrictions must not be imposed on HDBLAST without deriving equivalence to the physical shell law. |
| Souêtre, [arXiv:2507.20661v1](https://arxiv.org/html/2507.20661v1), 28 July 2025. Four-dimensional vacuum, negative cosmological constant, timelike **conformal** boundary. | Its tensorial conformal formulation and auxiliary boundary system provide a further example of deriving geometric reflective/Robin data from gauge-dependent data, instead of assuming that a suggestive boundary name guarantees well-posedness. | An asymptotically AdS conformal boundary is different from the registered finite scalar shell. Its Neumann condition refers to boundary stress, not automatically to the rescaled physical tensor amplitude. |

The September 2026 paper is particularly close to the immediate **linear sourced-wave** problem, but its quantum/Hadamard sections need not be imported. The current task can use its classical PDE argument without changing the theory or adding quantum assumptions.

## Established methods that remain more direct

1. Ginoux and Murro, [arXiv:2007.02544v6](https://arxiv.org/html/2007.02544), 8 July 2024, Theorem 6.8 and Section 6.3, prove well-posedness for a normally hyperbolic vector-bundle wave operator with an admissible boundary space and compatible smooth initial/source data. Equation (6.7) and the following calculation show that covariant Neumann data are admissible. This gives a more general bundle framework than scalar Klein–Gordon alone. Application still requires the actual tensor operator, a suitable positive fiber metric and a verified boundary-domain equivalence. Its Section 7.1 Robin discussion uses a different symmetric-positive reduction and restricted coefficient signs; it is not a blanket arbitrary-Robin result.

2. Araneda and Dotti, [arXiv:1611.03534v3](https://arxiv.org/html/1611.03534), Section II.5, make the Robin boundary energy and its sign explicit. On \(x\in(-\infty,0)\), outward normal \(+\partial_x\), with
   \[
   u_{tt}-u_{xx}+Vu=F,\qquad u_x(0)=\gamma u(0),
   \]
   the conserved homogeneous stiffness energy contains
   \[
   E_\gamma=\frac12\int_{-\infty}^0
   (u_t^2+u_x^2+Vu^2)\,dx-\frac{\gamma}{2}u(0)^2.
   \]
   Thus an outward-form condition \(\partial_nu+\kappa u=0\) contributes \(+\kappa u(0)^2/2\). This identity follows directly by integration by parts and fixes the convention independently of a paper's terminology. An unfavorable boundary coefficient can produce instability even for nonnegative \(V\). Finite-time well-posedness is a different conclusion from global positive energy or stability.

3. Kudoh and Tanaka, [arXiv:hep-th/0112013v3](https://arxiv.org/html/hep-th/0112013), Section II.1–II.2, compute second-order metric/scalar bulk sources and Gaussian-normal gauge transformations in a scalar-stabilized five-dimensional Randall–Sundrum model. Equations (27)–(30) exhibit quadratic displacement/evaluation corrections; equations (35)–(38) keep boundary sources alongside the bulk sources. These provide a valuable bookkeeping template. Their background has Minkowski slices and two distinct branes, with static axial/spherical restrictions and copy symmetry. Their coefficient expansion also differs from a derivative/factorial convention. In particular, the authors expressly warn after equation (7) that the variables called “TT” at second order are not automatically the fully trace-free nonlinear metric perturbation. No source formula, cancellation or mode spectrum should be transplanted into HDBLAST.

4. Bridgman, Malik and Wands, [arXiv:astro-ph/0107245](https://arxiv.org/abs/astro-ph/0107245), supply bulk/brane gauge-invariant cosmological perturbation and Gaussian-normal reconstruction procedures. Ruiz, Rinne and Sarbach, [arXiv:0707.2797](https://arxiv.org/abs/0707.2797), supply a generalized-harmonic constraint/gauge/radiation boundary hierarchy. The former's master-variable simplifications depend on its background; the latter studies vacuum outer/absorbing boundaries. Neither directly replaces the physical second-order Israel/scalar data.

Mars's second-order hypersurface perturbation formulas and the previously checked vacuum umbilic/finite-boundary results remain relevant from the earlier literature checkpoint. The new papers do not remove their dimensional, matter or boundary-condition limitations.

## Consequences for the calculation now underway

The material equation gives a necessary Neumann condition for the **rescaled dimensionless tensor**, not for its warped metric component. Once a sourced tensor mode and that boundary law are derived independently, a causal solution can be organized as a compatible homogeneous initial solution plus a retarded/Duhamel response. The total source must include quadratic geometry and scalar terms before projection. A nonzero matter-only projection neither proves nor rules out cancellation of the total source.

For Robin or conjugated tensor variables, keep the actual boundary potential in the energy. Use the trace inequality to absorb an unfavorable boundary term into a sufficiently large bulk \(L^2\) shift. Adding an arbitrary positive boundary \(u^2\) control term is not justified by the basic bulk \(H^1\times L^2\) energy: its derivative requires the uncontrolled trace of \(u_t\). An inhomogeneous boundary source needs a proper lifting or a separately proved trace estimate, with its time derivatives and corner jets included.

The compact-support shell collar is useful: initial perturbations and their recursively derived jets vanish near the initial corner. At later contact the physical source and boundary forcing must come from the same first-order solution and material pullback; compatibility cannot be inferred merely from initially vanishing boundary values.

The remaining full-model question is whether the geometric Israel/scalar conditions can be completed with compatible generalized-harmonic gauge and constraint boundary conditions, and whether the resulting boundary symbol and energy domain satisfy the requirements of a coupled Einstein–scalar theorem. A solved tensor amplitude alone does not settle the scalar/vector/gauge sectors, exact nonlinear contact continuation, global stability or cosmological observational predictions.

## Prior local material and search coverage

The companion [PRIOR_METHODS.md](PRIOR_METHODS.md) and [ARCHIVE_GREEN_METHODS.md](ARCHIVE_GREEN_METHODS.md) record useful local material actually read: warped tensor rescaling and Neumann laws, Sturm–Liouville/Wronskian and Prüfer checks, fixed-domain Robin energy, recursive corners and matched graph geometry. Their source manifests and verification receipt pin 21 files; six duplicate active/archive copies matched byte-for-byte. No old scientific script or opaque cache was executed and no old background numbers were adopted. The January RS2 notes contain sign/rescaling cautions explicitly recorded in the companion inventory.

Coverage was bounded: twelve targeted paper-index searches, one citation-graph search, relevant indexed passages, and six primary full-text HTML retrievals, concentrating on timelike-boundary waves, second-order scalar-supported branes, and 2025–2026 Einstein boundary work. The local search covered the accessible HDblast active/archive/work directories, not an inaccessible personal computer. The September 2026 Robin abstract endpoint failed once; its primary full HTML and separate metadata index succeeded and were checked. This is not an exhaustive web search or a novelty certificate. No external correspondence or publishing was performed by this literature subtask.

Machine-readable search/source records are in SOURCE_LEDGER.json. Downloaded full-paper text under sources/ is internal reference material; a public checkpoint can publish this synthesized note and source ledger without republishing those texts.
