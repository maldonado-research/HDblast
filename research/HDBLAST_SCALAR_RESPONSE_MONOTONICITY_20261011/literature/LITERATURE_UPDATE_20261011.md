# Targeted literature update for the regular scalar-boundary branch

Research session: 11 October 2026. Primary-source readbacks were made at approximately 07:53–07:56 UTC; the accompanying JSON receipts retain exact timestamps, URLs, byte counts and SHA-256 hashes. This is a bounded literature update, not an exhaustive search or a priority certificate.

**Integration update, later on 11 October:** the priorities below record the
literature review at the start of this calculation. The accompanying
[completed result](../RESULTS.md) now supplies the sharper local theorem and
exact bounds through `delta = 0.002`. The public 27 September checkpoint already
contains gauge-invariant coupled perturbation equations and numerical spectral
diagnostics for the corresponding near-`phi=1` branch. The next stability task
is to validate those existing operators and results over the certified interval,
including the small-detuning limit and precise parameter convention; it is not
necessary to restart that derivation. References below to needing a HDBLAST
perturbation problem mean a fully validated calculation in the present scope,
not an assertion that no prior project perturbation calculation exists.

**The immediate research priority remains a sharper, explicitly controlled nonlinear boundary-response theorem.** The most directly useful additional precedent is the established Hamilton–Jacobi treatment of curved Einstein–scalar walls. Two 2026 physics preprints supply useful future perturbation and observable calculations, and a 2026 revision of a mathematical paper sharpens the distinction between full boundary data and a few response coefficients. None supplies evidence that a higher-dimensional blast caused the Big Bang.

## What was searched and read

Six keyword queries were made through the arXiv API, INSPIRE and Crossref. They returned 82 metadata entries before removing duplicates, irrelevant results and overlaps. This is not 82 papers read. The first 20 results of each of two larger arXiv searches were inspected; those searches reported 138 and 62 total matches respectively. Other searches returned 7 arXiv entries, 5 and 15 INSPIRE entries, and 15 Crossref entries. Three additional identifier lookups checked source metadata. Exact queries are retained in `SEARCH_LOG.json`.

The queried concepts were de Sitter branes, curved domain walls, holographic self-tuning, Einstein–scalar curvature/boundaries, Hamilton–Jacobi methods, and Dirichlet-to-Neumann maps. Six author PDFs were fetched and the specific sections identified below were read. The seven-paper 8 October comparison and its cached primary texts were also consulted. No claim is made that all later citing literature or every relevant mathematical-relativity paper was searched. The 2026 papers below were newly incorporated into this project review; they were not released after the 8 October checkpoint.

Initial sandbox requests returned `Operation not permitted`. The supported approved network retry succeeded with inherited proxy settings and verified TLS. The runtime reported current observations but policy readiness `unknown`; successful HTTP responses, rather than an assumed readiness state, establish these particular accesses. No proxy or TLS bypass was used. Third-party PDF/text bodies remain in the private source cache and should not be included in the public research packet.

## 1. Direct mathematical precedent: curved-wall Hamilton–Jacobi theory

K. Skenderis and P. K. Townsend, *Hamilton-Jacobi method for curved domain walls and cosmologies*, **Physical Review D 74, 125008 (2006)**, [doi:10.1103/PhysRevD.74.125008](https://doi.org/10.1103/PhysRevD.74.125008), [hep-th/0609056v2](https://arxiv.org/abs/hep-th/0609056v2). Read §§1–2 and §4, especially Eqs. (1.1)–(1.4), (2.1)–(2.4), (4.7)–(4.15), and (4.23).

The paper already formulates curved single-scalar domain walls using a constrained radial Hamiltonian. It explicitly includes de Sitter slices. In their notation the lapse imposes the Hamiltonian constraint, and the momenta are derivatives of a Hamilton–Jacobi characteristic function. Equation (4.12) contains an integrated term proportional to the slice-curvature parameter. Their flat superpotential formula cannot simply be reused without this contribution on a curved branch. The authors also state that the first-order construction for de Sitter-sliced walls is not a BPS construction.

The normalization can be translated at the bulk-action level: their scalar is \(\sigma_{ST}=\sqrt{2}\phi_H\), their potential is \(V_{ST}=2U_H\), and \(e^{\beta_{ST}\varphi_{ST}}=\rho_H\). Their curvature symbol \(k_{ST}=+1\) specifies unit de Sitter slices and is **not** HDBLAST's inverse AdS radius \(k_H\). This comparison does not identify their complex generating function with HDBLAST's chosen real \(W\) on a curved solution.

The already reviewed Ghosh et al. paper [1711.08462v3](https://arxiv.org/abs/1711.08462v3), §2.2, Eqs. (2.15)–(2.20), makes the curvature mismatch especially explicit:

\[
S^2-SW'+\frac{2}{d}T=0.
\]

Its text immediately after Eq. (2.20) explains that the difference between \(S\) and \(W'\) measures curvature. Appendix B relates that formulation to Skenderis–Townsend. These functions are evaluated along one radial solution using the scalar as a coordinate; that construction requires care at a scalar turning point. It does not directly equal differentiation with respect to the endpoint Dirichlet value across a family of regular interiors.

**Impact on the next calculation.** Exploiting the Hamiltonian constraint before taking absolute-value estimates is well motivated. A regular-center Jacobi or Wronskian identity that preserves a small curvature factor may substantially improve a sufficient monotonicity estimate. Its possible contribution is the exact identity with stated endpoint hypotheses, a nonvanishing response denominator, and an evaluated consequence for this particular boundary problem. Hamilton–Jacobi reduction, scalar/metric coupling and curvature-dependent mismatch are established methods, not new physical principles.

The reviewed sources do not establish that the proposed particular response identity is absent from all literature. Any final paper should say what the new estimate proves and by how much it improves the certified domain; it need not claim priority for the underlying formalism. In particular, radial branch uniqueness is not Lorentzian dynamical stability.

The earlier effective-action comparison also remains relevant: Charmousis, Kiritsis and Nitti, [1704.05075v3](https://arxiv.org/abs/1704.05075v3), Appendix C, Eqs. (C.9)–(C.13), integrates out the two bulk regions and obtains both interface matching conditions by varying a boundary action. Their fixed UV data and integration parameter are not additional variables to minimize. Equality of mixed derivatives or a response reciprocity relation is natural in such a framework, but neither that fact nor a negative static susceptibility establishes a positive dynamical Hessian in the constrained gravitational theory.

## 2. A 2026 observable calculation worth developing after the static theorem

S. Aoki, *Primordial Correlators from a Kaluza–Klein Graviton Continuum*, [arXiv:2608.01762v1](https://arxiv.org/abs/2608.01762v1), submitted 3 August 2026. No journal DOI was present in the arXiv metadata read here. Read §§2.3–2.4, §§3.1–3.5, especially Eqs. (3.41)–(3.55), (3.71)–(3.79), (3.82)–(3.90), and §4.

This preprint studies an RS2-like de Sitter brane with two reflection-symmetric AdS sides and a brane-localized inflaton. It derives the transverse-traceless graviton Schrödinger potential and the brane Robin condition. In that model the potential approaches \(9H^2/4\), with a graviton zero mode and a continuum beginning at \(m=3H/2\). The brane wavefunction sets the continuum spectral weight. Near \(\mu=\sqrt{m^2/H^2-9/4}=0\), the paper finds a weight proportional to \(\mu\). Integrating over the continuum dephases the usual isolated-mass logarithmic clock signal in a primordial four-point seed function.

The paper itself limits the displayed exchange calculation to helicities \(\pm2\) and identifies the other massive-spin-2 helicities, realistic slow roll, conversion to curvature perturbations, and observational templates as further work. These qualifications prevent treating its plotted seed as a ready-made observational test of HDBLAST.

**Useful next project:** after fixing a dynamically consistent HDBLAST background and perturbation problem, calculate its normalized brane-to-brane tensor spectral measure and compare it to this known constant-scalar RS2/de Sitter baseline. Boundary scalar backreaction could alter the interior potential or spectral weight, but that alteration has not been calculated here. The mass gap and threshold law above belong to Aoki's specified model and must not be imported without deriving HDBLAST's operator, shell condition and asymptotic domain.

This also reinforces the existing identifiability issue: a four-dimensional continuum with the same spectral measure and couplings can reproduce the same exchange correlator. One correlator or one spectral weight is not, by itself, evidence uniquely selecting a higher-dimensional origin, much less a blast origin. A useful distinguishability calculation must specify the comparison class and additional jointly constrained observables.

## 3. A 2026 perturbation result: spectra require their own analysis

H.-L. Jia, W.-D. Guo, Y.-T. Gu and Y.-X. Liu, *Real Part Emergence in Purely Imaginary Quasinormal Modes in Perturbed de Sitter Braneworlds*, [arXiv:2606.04424v1](https://arxiv.org/abs/2606.04424v1), submitted 3 June 2026. No journal DOI was present in the arXiv metadata read here. Read §II.A–B, Eqs. (2)–(16), and the conclusion.

The example is a smooth **thick** brane, with an unperturbed warp factor \(A_0(z)=-n\log\cosh(\beta z)\), \(0<n<1\). The analysis concerns transverse-traceless gravitational perturbations. It uses a reduced potential \(V_{re}=V-9\alpha^2/4\), and its standard frequency-domain equation applies at zero spatial momentum; the nonzero-momentum equation retains the explicitly time-dependent \(e^{-2\alpha t}p^2\) term. A double-kink deformation of the background scalar changes the warp factor and the potential. Their discussion notes that the deformation changes the asymptotic bulk geometry as well.

The authors report that modes which were purely imaginary acquire real parts after the deformation, with transient oscillatory time-domain signatures. Their conclusion also says the effects are generally limited in the time domain and largely negligible with present observational capabilities. This is a statement about the characteristic spectrum of their model. It is not a theorem that every scalar-supported de Sitter brane is unstable or that a static curvature susceptibility determines stability.

**Impact:** HDBLAST needs its own coupled scalar, metric and brane-displacement perturbation problem, not just a copied tensor potential. A tensor-sector check would be an informative first step, but would not settle scalar ghosts, tachyons, boundary kinetic terms or full nonlinear stability. The current radial existence/error theorem should retain that scope boundary.

## 4. A recent mathematical inverse-boundary result clarifies identifiability

A. Enciso, G. Uhlmann and M. Wrochna, *The Dirichlet-to-Neumann map on asymptotically anti-de Sitter spaces and holography*, [arXiv:2512.03587v2](https://arxiv.org/abs/2512.03587v2), first submitted 3 December 2025; version read dated 14 March 2026. No journal DOI was present in the arXiv metadata read here. Read §1.1–1.3 and Theorems 1.1–1.3, together with their stated hypotheses.

This is a mathematical result for a **linear Klein–Gordon field on a fixed asymptotically AdS spacetime**. It identifies the forward Dirichlet-to-Neumann map microlocally with a fractional boundary wave operator modulo lower-order terms. Theorem 1.2 says that equality of the full maps determines the boundary Taylor series of the metric outside a possible countable set of mass parameters. Real analyticity, or the specified Einstein and null-convexity assumptions, yields further isometry conclusions. The maps use asymptotic data at conformal infinity and forward-supported spacetime fields.

HDBLAST's scalar response \(p(T,\beta)\), or the single linear ratio \(F'(T)/F(T)\), is a symmetry-reduced static endpoint response at a finite shell. It is not the full spacetime Dirichlet-to-Neumann operator in this theorem. The HDBLAST bulk is backreacted by the scalar and cut off by an interface, so the theorem does not automatically apply.

**Impact:** the paper offers a useful conceptual route from sufficiently rich boundary-response data to geometry *within a specified geometric model class*. It does not overcome HDBLAST's continuum-bath equivalence with a single spectral or curvature coefficient. Nor does uniqueness among metrics exclude nongravitational realizations outside that class. A serious next inverse problem would first define the available measurements, boundary, admissible bulk class and equivalence relation.

## 5. Nonlinear geometric Dirichlet-to-Neumann maps are also established

S. Blitz, A. R. Gover, J. Kopiński and A. Waldron, *The Dirichlet-to-Neumann Map for Poincaré–Einstein Fillings*, [arXiv:2307.08470v2](https://arxiv.org/abs/2307.08470v2), first submitted 17 July 2023; version read dated 23 October 2025. No journal DOI was present in the arXiv metadata read here. Read the introduction and the setup/statements of Theorems 1.1–1.4.

The nonlinear map here concerns Riemannian Poincaré–Einstein fillings, conformal boundary metrics, and a boundary tensor related to variation of renormalized volume. The paper constructs conformal hypersurface invariants; odd bulk dimensions require additional boundary qualifications in the stated theorems. The four-dimensional example involves a Cotton-tensor boundary quantity, not a scalar junction susceptibility.

**Impact:** calling a boundary response “nonlinear Dirichlet-to-Neumann” is useful standard terminology, not itself an originality claim. This paper does not supply HDBLAST's finite-shell scalar branch or its explicit detuning radius. Its boundary-at-infinity, Riemannian, pure-Einstein assumptions differ from a doubled Lorentzian Einstein–scalar interior.

## 6. A recent no-go theorem that should not be misapplied

F. Niedermann and A. Padilla, *Higher co-dimension de Sitter branes*, [arXiv:2506.19515v1](https://arxiv.org/abs/2506.19515v1), submitted 24 June 2025. INSPIRE also returned publication DOI [10.1103/9423-j3hl](https://doi.org/10.1103/9423-j3hl); the text read was the arXiv v1 author PDF, not a separate journal readback. Read §§II–IV, including Eqs. (1)–(15) and the discussion of thin-source limits and finite thickness.

The authors extend integrated-curvature no-go arguments to **codimension two and higher**, with compact internal space, specified bulk stress assumptions, and controlled treatment of localized sources. Their conclusions constrain de Sitter solutions for codimension two and require negative effective tension contributions in higher codimension under the stated assumptions. The discussion identifies finite thickness and localized curvature terms as possible changes to the setup.

HDBLAST's present four-dimensional shell in five dimensions has codimension **one**. This theorem is therefore not a direct exclusion of the current static branch. Conversely, the existence of that codimension-one branch does not evade or refute the cited higher-codimension theorem. Any future extension to more extra dimensions must reconsider the hypotheses rather than carry over the five-dimensional proof unchanged.

## Research priorities identified by this review

The [accompanying research checkpoint](../README.md) records the completed response calculation and evaluated bounds selected by items 1–2. This report supplies its literature context; it does not independently certify that calculation.

1. Finish the exact regular-center response identity and independent sign/denominator review; use it to replace a loose sufficient estimate only after every altered hypothesis and constant has been checked.
2. Evaluate the improved certified interval, curvature remainder and comparison margin on the already registered models. Distinguish proof coverage from successful floating-point trajectories and preserve the previous accepted certificate.
3. Keep dynamical stability as a separate next calculation with the coupled scalar/metric/interface action and boundary conditions stated explicitly.
4. Build an eventual observable comparison from the brane-to-brane spectral measure and a declared alternative-model class. Aoki's 2026 calculation is a concrete baseline, not evidence for this hypothesis.
5. Cite the old Hamilton–Jacobi and curved-flow formalisms alongside any improved estimate. The contribution should be the particular theorem and usable quantitative gain. This limited review establishes no “first-ever” or APS-significance claim.

This review does not assess journal acceptance or establish a cosmological origin. The accompanying checkpoint records mathematical results; publication and deployment history are recorded separately.
