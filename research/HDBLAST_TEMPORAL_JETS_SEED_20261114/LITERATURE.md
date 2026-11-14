# Method sources and scope of the compatible seed

Sources checked on 14 November 2026 UTC. Semigroup smoothing is an
established way to construct trace right inverses. This capsule evaluates
a finite compatible lift in its specific mixed Sobolev hierarchy; no
external qualitative constant is used as a numerical seed bound. This
focused review establishes neither mathematical novelty nor physical
evidence for a higher-dimensional blast.

## Primary texts inspected

| Versioned primary source | Supporting passage and applicability |
|---|---|
| Meyries–Schnaubelt, *Interpolation, embeddings and traces of anisotropic fractional Sobolev spaces with temporal weights*, [arXiv:1202.3870v1](https://arxiv.org/abs/1202.3870v1), 17 February 2012 | [Lemma 4.3, PDF p.20](https://arxiv.org/pdf/1202.3870v1#page=20), constructs an initial-trace right inverse using `R_Omega exp[-t(1-Delta)^m] E_Omega`. It assumes a Banach space of class HT, `1<p<infinity`, `1/p<mu<=1`, integer `m>=1`, and the stated smooth-domain geometry. Its Besov input and weighted maximal-Lp output give a single parabolic trace, rather than this multiple-jet continuous mixed norm. General surjectivity after Theorem 4.2 is only suggested; the proved special Lemma 4.3 is the relevant precedent. ArXiv comments report publication in *Journal of Functional Analysis* 262 (2012), 1200–1229. |
| Grubb, *Regularity of spectral fractional Dirichlet and Neumann problems*, [arXiv:1412.3744v4](https://arxiv.org/abs/1412.3744v4), revised 28 June 2015; [related DOI](https://doi.org/10.1002/mana.201500041) | [Definition 2.1 and Theorem 2.2, PDF p.3](https://arxiv.org/pdf/1412.3744v4#page=3), identify positive-power domains of smooth elliptic Dirichlet/Neumann realizations with Sobolev spaces satisfying all meaningful iterated boundary conditions `B A^l u=0`, including a special condition at critical indices. Arbitrary H^r data cannot be substituted for these domains. The actual cylinder, tensor and graph-norm constants are evaluated separately here. |
| Velázquez-Mendoza–Sandoval-Romero–Diaz Carlos, *An Intrinsic Sobolev Trace Theory for Riemannian Vector Bundles with Applications to Quasilinear Elliptic Equations*, [arXiv:2609.26744v1](https://arxiv.org/abs/2609.26744v1), 22 September 2026 | [Theorems 4.6–4.8, PDF pp.11–13](https://arxiv.org/pdf/2609.26744v1#page=11), give integer-order covariant spatial traces for finite-rank metric bundles with compatible connection on compact smooth manifolds: `Tr_j:W^(m,p) -> W^(m-j-1,p)`, `j<=m-1`, for `m>=1`, `1<=p<infinity`. They identify the closure of smooth interior-supported sections as the kernel of these traces. This recent preprint supplies related spatial trace language, not a temporal jet lift, mixed D/N spectral realization or numerical seed constant. Journal status was not independently established. |
| An–Anderson, *Well-posed geometric boundary data in General Relativity, II: twisted Dirichlet boundary data*, [arXiv:2507.15567v3](https://arxiv.org/abs/2507.15567v3), revised 30 May 2026 | Theorem 1.1 is a smooth tame Fréchet vacuum IBVP result for different boundary data. [Remark 5.2, PDF p.22](https://arxiv.org/pdf/2507.15567v3#page=22), proposes a finite-derivative Nash–Moser extension with loss three but explicitly leaves verification of that theorem's hypotheses unchecked. It supplies no completed finite-Sobolev theorem for the present Israel/scalar rows. Previously inspected passages were reused. |

Dates identify the exact arXiv versions read, rather than internal PDF
dates. Selected passages were inspected; complete papers and proofs are
not claimed audited. The versioned preprint text is distinguished from
any related journal publication.

## What is evaluated in this application

The retained [mixed evolution proof](inputs/QUASILINEAR_VARIABLE_MIXED_EVOLUTION.md),
Section 3.1, uses a smooth Borel lift, while its Section 5 describes a
qualitative short precontact solution as a compatible seed. Neither
argument supplies the finite numerical seed constant by itself. The
[current continuation interface](inputs/CURRENT_CONTINUATION.md) instead
requires a bound for jets `J_j` in H^(14-j), followed by a compatible seed
on the prescribed elapsed-time slab.

For the auxiliary positive self-adjoint A, the corrected finite heat
polynomial uses terms `A^r J_(j-2r)`. Their derivative count matches the
spectral scale `D(A^((14-j)/2))`. The positive heat counterterms come from
a finite formal Taylor identity; no backward heat operator is applied to
rough data. Gaussian multiplier bounds and spectral dominated convergence
must establish strong initial traces, including the top L2 time trace.
A bound for positive time alone would not establish those traces.

The common initial zero collars provide the iterated boundary
compatibility. The explicit spectral/frame comparison, all derivative
slots, and cutoff multiplier costs remain part of the evaluated proof.
A fixed radial cutoff equal to one on the initial support removes heat
tails. Applying the natural finite reflection average afterward preserves
the required symmetry without assuming that a component Laplacian
commutes with that tensor action.

The degree-fourteen construction preserves interior jets through order
fourteen. Its zero collars also give the smooth boundary-corner zero
conditions. These are distinct assertions; an order-fifteen smoothing
argument needs its own additional interface. The auxiliary heat lift
implies no physical finite-propagation or nonlinear evolution theorem.

## Bounded collection and archive relevance

The connector was unavailable for collection, so this review used direct
primary arXiv metadata and versioned PDFs. Nine bounded API queries
covered semigroup trace lifting and recent Sobolev jet/extension work.
The recent mathematical trace query returned 45 records; only its five
newest were screened, and the bundle-trace paper above received passage
inspection. This is not an exhaustive survey or an originality finding.
No screened recent result replaces the application-specific spectral
domain, strong-trace or numerical norm calculations.

Prior mounted D-Blast 3 relevance reports were reused. Two new searches
of visible Markdown/text files returned fourteen broad semigroup/lift
matches and zero matches for a narrower heat-lift pattern. Inspection of
the promising scopes found a Beta–Jacobi flow for PTA residuals, a
characteristic adjoint trace-discrepancy lift, and a proposed empirical
frequency-profile extension. These concern different operators and
questions; no archived numerical constant was transferred. No archive
solver was executed, and no whole-computer coverage is claimed.

The contribution here is a quantitative application of established
analysis to the retained initial-data construction. Nonlinear spacetime
energy, continuation, parameter estimates and proper-observer tube bounds
remain separate mathematical tasks. This literature check supplies no
observational calibration or empirical confirmation.
