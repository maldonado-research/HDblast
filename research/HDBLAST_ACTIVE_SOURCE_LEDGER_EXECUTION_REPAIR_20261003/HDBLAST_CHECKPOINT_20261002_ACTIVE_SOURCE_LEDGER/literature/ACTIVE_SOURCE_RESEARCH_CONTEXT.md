# Primary-literature context for the next active-source conservation test

Reviewed October 2, 2026 in America/Los_Angeles. This is a bounded internal
AI-assisted methods review. It evaluates no physical modes, response integrals,
acceptance thresholds or trajectories. The completed source-free result remains
limited to its registered interval; all earlier metric failures remain failures.

Eight first-page INSPIRE, arXiv and Crossref searches returned 109 metadata
records including duplicates. One narrow older-source follow-up returned 15
records; a selected-ID request checked four version entries. Four primary PDFs
were retrieved and selected passages read. A fifth primary article was read in
the full-text transcription embedded in its INSPIRE record. Three previously
retrieved numerical-method papers remain inherited context; selected passages
in two of those were reread. These counts describe access, not five complete
paper verifications or an exhaustive survey. Initial sandbox sockets were
blocked; authorized read-only network execution then completed the requests.

## Immediate priorities supported by the literature

**Construct the stress and source observable in one convention.** For a real
scalar coupled to an external scalar potential, Stephen A. Fulling, Thomas E.
Settlemyre and Kimball A. Milton give

\[
\partial_\mu T^\mu{}_{\nu}+\tfrac12(\partial_\nu V)\phi^2=0.
\]

Their [2018 article](https://doi.org/10.3390/sym10030054), §3.4, Eq. (34),
checks cancellation between the regularized stress and field-square terms.
In particular, a local term proportional to `V g_mu_nu` need not be conserved
by itself: it participates in the source-exchange identity together with the
corresponding field-square term. Their worked potential is static and flat
spacetime is used. This supplies a clear provenance for testing paired
stress/source subtraction and matching; it supplies no coefficient for the
present curved, time-dependent metric response. Translate signature, index
placement and source normalization from the actual HDBLAST action before
using an energy-density version of this identity. Do not change finite
couplings or contacts after seeing a conservation defect.

**Linearized conservation includes the variation of the connection.** Tomislav
Prokopec's [2512.22958v2](https://arxiv.org/abs/2512.22958v2), first submitted
December 28, 2025 and revised March 7, 2026, was read in selected passages on
printed pp.23–28 and46–48. Its full-metric classical, quantum and geometric
counterterm stresses obey their background identities. After expanding the
metric, the separated gravitational response kernels generally have
non-transverse terms tied to background one-point functions. The complete
linearized equation cancels these terms when the background equation holds.
Equations (95)–(99) explicitly retain first and second variations before the
paper sets the scalar perturbation to zero. The discussion explains the
connection to the in-in formalism; the paper primarily derives the in-out
case. The practical application is to derive the complete HDBLAST response
identity, including connection/background terms and any independently present
source-current term. A zero divergence demanded of an isolated separated
kernel is not automatically the correct check. Neither this source nor its
de Sitter examples remove the need to derive HDBLAST's own causal kernel,
state, source contacts and finite-cutoff convention.

**Keep quadrature accuracy separate from trajectory accuracy.** The inherited
[Iserles–Maierhofer 2404.11448v2](https://arxiv.org/abs/2404.11448v2) uses a
weight obeying `w'=G w` and a nonoscillatory auxiliary equation whose solution
reduces an oscillatory integral to endpoints. For a chosen constant scalar
carrier `exp(i Omega t)`, let an independently constructed envelope approximation
be `p` and define `r=f-(p'+i Omega p)`. Then the exact integration error of the
endpoint expression is `integral r exp(i Omega t) dt`, bounded by
`integral |r| dt` if the latter bound is independently justified. This is the
standard residual identity, not a new method or a supplied continuum
certificate. Values of `r` at saved samples do not bound its unsampled values.
It is useful to record such a residual separately from the mode-flow defect
and from the conservation residual. A ledger defined from measured final
density would destroy the independence of the check.

[Lee–Appelö 2607.13580v1](https://arxiv.org/abs/2607.13580v1), Theorem3.1 and
Remarks3.6–3.8, bounds Filon envelope-interpolation error using derivatives of
the envelope and warns that a poor carrier choice makes those derivatives
large. The source-active interval needs its own envelope and derivative
contract; increasing an arbitrary frequency cannot manufacture accuracy.
The bound is for its stated controlled-system quadrature, not a transferred
error constant for renormalized HDBLAST stress. The inherited
[Domínguez 2503.08169v3](https://arxiv.org/abs/2503.08169v3) remains useful for
stable polynomial-exponential moments; no new passage reading or new method
implementation from that paper is claimed here. Interpolation error,
oscillatory-moment arithmetic, mode-flow error and momentum truncation/refinement
should remain separately observable under the predeclared old tolerances.

## Additional primary results and their limits

| Primary source, access and date | Useful application | Limit |
| --- | --- | --- |
| Matías I. Caruso, Javier Fernández, Cora Tori and Marcela Zuccalli, [Variational integrators using forced discrete Hamiltonian systems, 2607.02694v1](https://arxiv.org/abs/2607.02694v1), July2, 2026; selected pp.2–3 and13–15, Proposition4.10/Remark4.11 | Derive the discrete balance/evolution law for the actual Hamiltonian and forcing before assuming an invariant. The paper explicitly gives the forcing-dependent evolution of the canonical symplectic form and identifies the closed-force exception. | Geometric preservation alone does not certify the stress-energy ledger, state accuracy or renormalization. HDBLAST's externally time-dependent conservative potential may be represented by a nonautonomous Hamiltonian; it is not automatically the dissipative-force example. No integrator replacement is selected on paper performance claims. |
| Beatrice Costeri, Claudio Dappiaggi and Michele Goi, [Conservation Law and Trace Anomaly for the Stress Energy Tensor of a Self-Interacting Scalar Field, 2411.07109v2](https://arxiv.org/abs/2411.07109v2), first submitted November11, 2024, revised April8, 2025, [journal DOI](https://doi.org/10.1007/s00023-025-01580-0); selected pp.1–4 and19–20, Remark3.3/Theorem3.4 | Its discussion explains why a classically vanishing EOM term can contribute after the composite operator is defined and why local covariance by itself does not guarantee quantum conservation. This reinforces an independent contact/operator audit. | The authors work in perturbative algebraic QFT and establish the interacting result under stated assumptions and orders. Their `eta=1/3` free Wick-ordering coefficient is prescription-specific; it is not a coefficient to copy into the current finite-band adiabatic benchmark. |
| Anamitra Paul and Sonia Paban, [Stress-Energy Tensor of a Scalar Field on a Product Spacetime with a Time-Dependent Compact Dimension, 2603.12444v3](https://arxiv.org/abs/2603.12444v3), first submitted March12, 2026, revised September18, 2026; [JHEP DOI](https://doi.org/10.1007/JHEP09(2026)188); selected printed pp.2–3,7,11,15–17 | Provides a directly relevant future dimensional-stress comparison: a varying compact scale carries an independent internal pressure and dimensionally reduced scalar dynamics. Their Eq.(4.22) checks the higher-dimensional divergence; Section6 reconstructs a reduced action and checks gravitational components. | The authors propose an approximate WKB prescription, test stated finite/conservation/conformal limits, and defer cosmological applications. It is a homogeneous product spacetime with a compact circle, not a moving shell or blast. It supplies no HDBLAST energy-transfer rate, shell junction, thermal history or validated Big Bang mechanism. |

The brane/energy search returned 11 metadata records, mainly late-time dark
energy, flux-vacuum or membrane results. None was selected as a repair for the
four-dimensional active-source calculation. This is a relevance decision for
this bounded sample, not a claim that no recent brane-energy-transfer paper
exists. Modified-dispersion renormalization, spin1/vector stress and unrelated
engineering oscillatory methods were screened at metadata/abstract level;
changing dispersion, field content or regulator would change the benchmark.

## Other owner projects and provenance

The public UTOE default head was freshly checked at
`c284afaa99c9a1e1fbde92a5547892f0ff9a2ce4`; its pinned
[`docs/RESEARCH_STATUS.md`](https://github.com/maldonado-research/Unified-Theory-of-Everything/blob/c284afaa99c9a1e1fbde92a5547892f0ff9a2ce4/docs/RESEARCH_STATUS.md)
has SHA256 `63cbed15d066ed61797d172cedb63d63d1faf1a159582f235fd3f867b1995536`.
Unlike the earlier ledger check, this note includes an October operator audit.
It preserves source and curvature contacts and warns that EOM/basis
cancellations cannot determine physical loop coefficients. That is a useful
bookkeeping principle. Its pulse remains prescribed in flat spacetime with
expansion, source dynamics, backreaction and complete energy accounting open;
no coefficient or physical validation is imported. No other owner repository,
private archive, local computer folder or D-Blast3 directory was searched in
this targeted context check.

The next scientific task should fix the source-active interval, stress/current
operators, subtraction/matching conventions, time and momentum refinements,
quadrature construction and acceptance gates before evaluating physical data.
The literature informs that design; it is not evidence that the resulting
experiment passes, an extra dimension was observed, or a new fundamental law
was discovered. Raw third-party PDFs and their extracted text stay outside the
public payload; public files contain metadata, hashes and original review.
