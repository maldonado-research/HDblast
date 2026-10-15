# Literature for angular-uniform contact estimates

Prepared 14 October 2026 in America/Los_Angeles; retrieval workspace dated
15 October 2026 UTC. Bounded primary-source method comparison. This note
does not assert external novelty, external peer review, or an exact
nonlinear shell-contact theorem supplied by another paper.

## Finding

The normal change of field variables proposed in this continuation has
a closely related published antecedent: Li–Lu–Rao convert a scalar
Robin channel into Neumann form by an exponential normal multiplier.
Their published theorem first diagonalizes a Robin matrix assumed
similar to a real symmetric matrix. The full nonsymmetric twelve-field
HDBLAST block therefore needs a separate proof; their theorem cannot be
quoted as already establishing that extension.

For the homogeneous transformed boundary problem, the decisive structure
is the common scalar wave principal part and vanishing normal–tangential
inverse-metric coefficients at the shell. Shibata's primary linear
Neumann theorem permits the zero-dissipation case, and the waveguide
literature explains the boundary-flux condition directly. These support
smooth frozen linear solves and energy bookkeeping. A componentwise mixed
Dirichlet/Neumann quasilinear lemma, finite Sobolev regularity budget,
initial harmonic completion and continuation to a prescribed contact
interval must still be proved for the actual formulation.

The earlier exact localized initial-data theorem is a legitimate
inherited input: at one fixed sufficiently high finite class it provides
an exact constraint family tangent to the reflection-invariant aligned
seed and equal to the background on an open initial shell collar. Its
short precontact corollary does not supply a prescribed contact lifespan
or a common spatial C-infinity parameter interval.

## Primary comparisons

| Source actually read | Exact location | Useful conclusion | Material restriction |
| --- | --- | --- | --- |
| Fournodavlos–Smulevici, arXiv:2104.08851v1 | Theorem 1; Section 2.2 (16)–(18); Section 4.2 Proposition 2; Appendix A Proposition 3/(29) | Complete the umbilic metric boundary problem using a normal harmonic row and derive tangential subsidiary laws; scalar homogeneous Robin energy uses a bulk L2 shift and trace inequality | Four-dimensional vacuum, constant umbilic factor; Appendix A is scalar and does not furnish a standalone mixed nonsymmetric matrix finite-Hm theorem |
| Li–Lu–Rao, arXiv:2202.00760v1 | Theorems 3.1–3.2; Theorem 4.1, (4.3)–(4.9); (2.6) | Scalar normal exponential conjugation converts Robin into Neumann plus first-order bulk terms | Their coupled Robin theorem assumes B similar to a real symmetric matrix; autonomous Euclidean wave/control problem |
| Shibata, Tsukuba J. Math. 12 (1988), 149–209 | Section 2, (A.1)–(A.6), Theorems 2.2–2.3; Sections 5–7 | Smooth linear vector Neumann existence and finite-order regularity, including non-strict boundary dissipation | All-Neumann theorem; explicit symmetry/coercivity and coefficient/source regularity; not by itself the nonlinear mixed-component theorem |
| Metcalfe–Sogge–Stewart, arXiv:math/0411513v2 | Section 5, (5.1)–(5.2), citing the energy estimate (4.3) | Homogeneous Neumann energy retains its form when normal–tangential principal flux vanishes | Scalar perturbative wave/Klein–Gordon problem in homogeneous waveguides; global dimensional and smallness assumptions differ |
| Kreiss–Reula–Sarbach–Winicour, arXiv:0807.3207v1 | Section I.1, Theorem 1, (1)–(4); Sections II–III | Common metric principal part on a vector bundle and a positive auxiliary fiber metric are compatible with a quasilinear wave IBVP | Boundary class contains a strictly timelike oblique/Sommerfeld vector and small derivative couplings; pure Neumann is not obtained by silently taking a singular limit |
| Kreiss–Reula–Sarbach–Winicour, arXiv:0707.4188 | Section 2, (1)–(3); Appendix 1 | Strong boundary stability includes control of boundary derivatives for the stated dissipative condition | Its condition has alpha>0 and beta_1^2+beta_2^2<alpha^2; it is not the homogeneous reflecting energy claim used here |

Links: [Fournodavlos–Smulevici](https://arxiv.org/abs/2104.08851v1),
[Li–Lu–Rao](https://arxiv.org/abs/2202.00760v1),
[Shibata primary PDF](https://core.ac.uk/download/pdf/56641220.pdf),
[Metcalfe–Sogge–Stewart](https://arxiv.org/abs/math/0411513v2),
[Kreiss et al. 2008](https://arxiv.org/abs/0807.3207v1),
[Kreiss et al. 2007](https://arxiv.org/abs/0707.4188).

## The boundary-input distinction

No estimate of an arbitrary inhomogeneous Neumann source by its bare
L2 spacetime boundary norm should be inferred from the homogeneous
energy. The boundary term contains the boundary trace of the velocity,
which the basic H1 bulk norm does not directly control.

Li–Lu–Rao Section 2 and Theorem 4.1 make this concrete. For their stated
L2 boundary controls the conclusion is in H^alpha times H^(alpha-1),
with alpha=3/5-epsilon on a general smooth domain in their stated result
and alpha=3/4-epsilon on a parallelepiped. These are the regularities
in that paper, not a claim that no stronger theorem exists: later
microlocal work improves some of the general-domain estimates.

Paunonen–Seifert–Vanspranghe, arXiv:2412.14786, Section 6.1 Theorem 6.4,
was inspected through relevant primary passages. Their transfer-function
analysis distinguishes eta-admissibility of Neumann control from ordinary
energy-space admissibility; L2 input produces H^(1-eta) times H^(-eta)
in the cited sharp regularity framework. The Schrödinger Theorem 6.11
is a different result and has not been substituted for the wave theorem.

An explicit smooth lifting is a valid way to formulate the present
higher-regularity forced problem. For example, in a product collar with
boundary z=0, a cutoff times z times a smooth Neumann datum has the
prescribed normal derivative. Subtracting it produces homogeneous
boundary conditions and an explicit bulk source. The estimate must
charge the derivatives of this lifting actually used by the wave
operator. It must not be described as a lossless L2 boundary-control
theorem.

## Current 2025–2026 comparisons

An–Anderson's series supplies genuinely new geometric vacuum boundary
results, but each boundary prescription matters:

- [Part I, arXiv:2505.07128](https://arxiv.org/abs/2505.07128) was screened
  by canonical metadata. Its Dirichlet theorem has a Brown–York
  Lorentz-signature/convexity hypothesis.
- [Part II, arXiv:2507.15567v3](https://arxiv.org/abs/2507.15567v3)
  was retrieved in full primary HTML and checked in the relevant sections.
  Theorem 1.1 concerns twisted Dirichlet
  data; Sections 3–4 allow derivative loss, and the nonlinear conclusion
  uses Nash–Moser. It does not justify discarding derivative-loss
  bookkeeping in another boundary problem.
- [Part III, arXiv:2503.12599v3](https://arxiv.org/abs/2503.12599v3)
  was retrieved in full primary HTML and checked in the relevant sections;
  its version header is dated
  18 June 2026. Theorems 1.1–1.2 establish injectivity and dense range
  under their smooth bounded-geometry assumptions, not general smooth
  surjectivity or a full nonlinear theorem. Section 3.4 Remark 3.10,
  (3.29), explicitly discusses the half-derivative issue for inhomogeneous
  scalar Neumann data. Its conformal-mean-curvature failure mechanisms
  cannot be relabeled as failures of umbilic Israel data.

[Hou–Yin–Yuan, arXiv:2602.14019v1](https://arxiv.org/abs/2602.14019v1)
was retrieved in full primary HTML and checked in the relevant sections.
Equation (1.3) explicitly imposes the
nonlinear Neumann admissibility condition; (1.4) states the corner
hierarchy. Their Theorems 1.1–1.2 concern small-data two-dimensional
exterior scalar waves, with additional dimensional, structural and
regularity hypotheses. This is current support for tracking normal
principal flux, not a direct five-dimensional matter-shell theorem.

[Ivanovici–Sánchez Sanz, arXiv:2609.40114](https://arxiv.org/abs/2609.40114)
was screened by canonical metadata and primary passages including the
introduction and Theorem 1.5. It constructs the glancing Neumann
parametrix for strictly convex fixed domains and obtains dispersive/
Strichartz conclusions. Its caustic-induced dispersive loss is a different
question from a uniform H1 homogeneous bulk energy bound. Glancing
does not alone falsify that energy estimate.

[Stevens–Valiente Kroon, arXiv:2604.24931](https://arxiv.org/abs/2604.24931)
was canonically rechecked; the previous primary-source ledger is reused
for the restriction to a conformal-geodesic finite timelike boundary and
a compatible maximally dissipative subset of vacuum conformal equations.
It is not an Einstein–canonical-scalar Israel theorem.

[Aigner–Skrepek, arXiv:2505.08602](https://arxiv.org/abs/2505.08602)
was checked by canonical metadata and primary passages. Its boundary
triple framework allows split displacement/velocity/stress boundary
conditions under positivity/self-adjointness/dissipativity assumptions.
It is a useful fixed linear operator comparison, not automatic coverage
of arbitrary nonsymmetric time-dependent matter-shell boundary matrices.

## Exact applicability of the proposed transformation

For a fixed linear Robin row -U_z+M U=g, choose an invertible smooth
collar matrix T with T_z=M T at the boundary and T(0)=I. Writing U=T V
makes that row -V_z=g at z=0. No diagonalization or positive symmetric
weight for M is needed for this algebraic identity. The common scalar
principal operator is unchanged after multiplication by T inverse;
derivatives of T enter lower-order bulk terms. A positive fiber norm on
V pulls back to a positive norm on U through T inverse, independently
of whether the Lorentzian metric-component pairing is positive.

For a nonlinear boundary row U_z=B(U), the flow identity
U=Psi(z,V), Psi(0,V)=V and Psi_z=B(Psi) similarly gives
U_z-B(U)=D_V Psi V_z. This is an algebraic reduction when the flow
exists and its Jacobian is invertible. Field-value dependence of Psi
preserves a common scalar principal symbol, but generates new lower-order
nonlinearities. The whole Dirichlet subbundle must be preserved, and the
physical metric used in the wave principal part must still satisfy the
boundary normal–tangential block condition.

These statements describe the method and its proof obligations. They
are not a theorem quoted from Li–Lu–Rao, Shibata, or a vacuum Einstein
paper. A complete proof must include the mixed high-order estimate,
existence/uniqueness and convergence, corner compatibility, subsidiary
propagation and the finite regularity lost under parameter
differentiation. The registered potential, scalar-dependent tension and
two-copy physical shell equations remain the model.

## Coverage and exclusions

This subtask performed ten targeted paper-index searches, thirteen
ordinary web searches, one bounded citing-paper expansion from
arXiv:2104.08851 (nine candidates returned), twelve canonical metadata
inspections, eleven indexed primary-passage retrievals, eight full
primary arXiv HTML retrievals with targeted section reading, two full
primary PDF parses with targeted theorem reading, and one
publisher metadata/reference-page reading. Search results and some
metadata are only discovery evidence, not full-text theorem evidence.

The Shibata–Nakamura 1989 quasilinear paper,
[DOI 10.1007/BF01180683](https://doi.org/10.1007/BF01180683),
was bibliographically verified on the publisher page. Its theorem body
was not available in that retrieval, so it has not been invoked with
unchecked hypotheses. The earlier Ikawa 1969 full primary theorem
was read and excluded as a direct pure-Neumann citation because its
stated main theorem uses a strict dissipative inequality.

No directly applicable published five-dimensional canonical-scalar,
variable-tension finite-shell contact theorem was established by this
bounded search. That is a coverage statement, not a claim that none
exists anywhere. No unrelated quantum or string sector is needed for
this boundary analysis. No external messages, Zenodo record or APS
submission were created. Raw full texts and tool responses remain
internal reference material; only the synthesis, source ledger and
verification notes are candidates for the research packet.
