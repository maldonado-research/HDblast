# Sources and scope of the nonlinear limit

Sources checked on 14 November 2026 UTC. Lower-energy convergence,
individual top-order time continuity, strong top-order convergence of
approximations and top-order initial-data continuity are distinct claims.
This capsule proves the first two by separate arguments. Its numerical
difference and factorial bounds come from the evaluated actual system;
no external qualitative constant is promoted to a numerical bound.

## Primary method passages

| Versioned source | Supporting passage and limit of use |
|---|---|
| Ifrim–Tataru, *Local well-posedness for quasilinear problems: a primer*, [arXiv:2008.05684v2](https://arxiv.org/abs/2008.05684v2), revised 23 April 2022 | [Section 4.1, PDF pp.12–14, especially Remark 4.2](https://arxiv.org/pdf/2008.05684v2#page=12), separates a uniformly high-bounded, lower-convergent iteration from strong continuity at the top Sobolev order. [Section 5.4, p.19, (5.11)](https://arxiv.org/pdf/2008.05684v2#page=19), adds frequency-tail control for strong top approximation. [Section 6.2, pp.20–21, (6.4)–(6.5)](https://arxiv.org/pdf/2008.05684v2#page=20), explains Kato/Bona–Smith comparison and places the smoother solution in the factor carrying the extra derivatives. Their full-space symmetric-hyperbolic model supplies a method guide, not this mixed-boundary theorem or its constants. |
| Ifrim–Pineau–Tataru–Taylor, *Sharp Hadamard local well-posedness, enhanced uniqueness and pointwise continuation criterion for the incompressible free boundary Euler equations*, [arXiv:2309.05625v3](https://arxiv.org/abs/2309.05625v3), revised 26 March 2025 | [Theorem 1.3, PDF p.8](https://arxiv.org/pdf/2309.05625v3#page=8), has the Euler Taylor-sign and problem-specific control-norm hypotheses. [Section 9, pp.104–109](https://arxiv.org/pdf/2309.05625v3#page=104), separately regularizes data, proves uniform bounds and lower convergence, then uses frequency-envelope tails for strong Sobolev convergence and continuity. This recent revision is a related method comparison; its moving fluid domain, pressure equation, norms and constants do not transfer to the scalar/Israel mixed system. |
| Velázquez-Mendoza–Sandoval-Romero–Diaz Carlos, *An Intrinsic Sobolev Trace Theory for Riemannian Vector Bundles with Applications to Quasilinear Elliptic Equations*, [arXiv:2609.26744v1](https://arxiv.org/abs/2609.26744v1), 22 September 2026 | [Theorems 4.6–4.8, PDF pp.11–13](https://arxiv.org/pdf/2609.26744v1#page=11), give continuous integer-order covariant spatial traces and a zero-trace kernel characterization on compact smooth metric bundles with compatible connection. They support ordinary boundary passage at available trace orders. They give neither a temporal compactness theorem nor an ordinary Neumann trace for a top commuted section that is only H1. Previously inspected passages were reused. |
| An–Anderson, *Well-posed geometric boundary data in General Relativity, II: twisted Dirichlet boundary data*, [arXiv:2507.15567v3](https://arxiv.org/abs/2507.15567v3), revised 30 May 2026 | Theorem 1.1 is smooth tame Fréchet vacuum theory for different boundary data. [Remark 5.2, PDF p.22](https://arxiv.org/pdf/2507.15567v3#page=22), proposes a finite-derivative extension with loss three but explicitly leaves its theorem hypotheses unchecked. It supplies no finite-Sobolev convergence or compatibility theorem for the present sixteen-field mixed system. Previously inspected passages were reused. |

Dates identify the exact arXiv versions inspected. The 2025 Euler paper
is a revision of a 2023 preprint. Journal publication status was not
independently established from the retrieved records. Selected passages,
rather than complete papers and proofs, were reviewed. Authors' novelty
claims are not adopted as a novelty determination for this capsule.

## Three separate top-order statements

| Statement | What is established here |
|---|---|
| One solution has a continuous mixed X14 jet table in time. | Established by the additional weak linear energy, mixed-boundary reflection and interior normal-recovery proof. |
| The original Picard or supported-data approximations converge strongly in C X14. | Not established. The proved lower convergence and uniform top bound do not imply it. |
| The initial-data-to-solution map is continuous into C X14. | Not established. Evaluated low-energy stability and qualitative continuity into C X13 are distinct from this claim. |

The lower limit first uses the evaluated consecutive-iterate source and
an exponentially weighted Volterra inequality. Ordered integration gives
the factorial denominator on the full prescribed interval. The seed
needs its known starting bound, not a PDE residual. Spatial compactness
upgrades the field/velocity pair below top order; the exact shifted
equation then supplies the other temporal jets. Temporal derivatives are
not interpolated as unrelated spatial data.

The first statement in the table requires an **additional PDE argument**.
[WEAK_LINEAR_ENERGY_REGULARITY.md](WEAK_LINEAR_ENERGY_REGULARITY.md)
proves a closed-manifold linear energy lemma directly, using spatial
Friedrichs commutators, Lipschitz symmetrizable principal coefficients,
an L1 L2 source and weak time traces. Its energy identity yields
continuity of the energy; positivity plus weak continuity then gives
strong time continuity. The qualitative regularization constants change
none of the numerical admission or factorial bounds.

[TOP_ORDER_TIME_CONTINUITY.md](TOP_ORDER_TIME_CONTINUITY.md) verifies the
missing boundary interface before applying that lemma. Highest
time/rotation-commuted sections satisfy a weak mixed conormal identity;
an unavailable classical Neumann trace is not assigned to them. Parallel
fiber reflections and the actual zero cross-coefficient boundary traces
extend the linear commuted equation to a radial circle of circumference
4Z. This preserves Lipschitz coefficients, positivity and the weak
equation. The resulting continuous tangential energies, followed by the
accepted interior normal recurrence, give individual strong C X14.

Thus X15 smoothing is **not logically required for this individual
time-continuity conclusion**. The Ifrim–Tataru passages explain why a bound
alone is insufficient and provide routes to stronger approximation and
data-continuity statements; they do not exclude the separate weak-energy
argument used here. Strong top approximation and top initial-data
continuity still need their own estimates, such as an appropriate
smoothing comparison or substitute.

The integrated
[reduced evolution and causal-support theorem](REDUCED_EVOLUTION_AND_CAUSAL_SUPPORT.md)
also proves finite propagation from the exact initial support using the
full radial flux. Its speed ceiling is c=20/19 in the registered elapsed
and radial coordinates, giving the outer support cap 80ell/19 and shell
nonarrival before 19ell/40. These are coordinate statements and do not
assert that a signal arrives at that endpoint, provide proper-observer
calibration, or establish physical evidence. Subsidiary, harmonic and
Israel-condition propagation, evolved parameter derivatives and observable
remainders remain separate tasks.

## Bounded discovery and archive coverage

The connector was unavailable for collection, so direct primary metadata
and versioned PDFs were used. A focused title/author query located the
primer and related recent records; the 2025 Euler revision received
targeted convergence-passage inspection. The 2026 trace comparison and
An–Anderson limitation reuse exact earlier passages. This is a bounded
method screen, not an exhaustive recent survey or a claim that no other
relevant method exists.

A targeted search of visible D-Blast 3 Markdown/text files for factorial
Picard convergence, Bona–Smith, weak H14/H13, Rellich and top-continuity
terms returned five filenames, four being related handoff copies. The
inspected original supplement bounds an older two-by-two forcing
majorant for a different source/response problem. Handoff contexts repeat
that result and warn against transferring parabolic smoothing to a
hyperbolic Goursat problem. No archived convergence constant or theorem
was transferred. No archive producer was executed, and no hidden files,
ZIP contents, unmounted storage or whole-computer coverage is claimed.

The resulting theorem is a conditional mathematical construction for the
retained reduced system. This literature review supplies no empirical
confirmation, observationally calibrated amplitude or higher-dimensional
origin claim.
