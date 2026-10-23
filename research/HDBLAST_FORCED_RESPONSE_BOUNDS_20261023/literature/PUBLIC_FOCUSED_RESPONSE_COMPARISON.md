# Focused literature update for the forced-response certificate

Checked 22 October 2026 in America/Los_Angeles (23 October UTC). This is a bounded update of the 16 October comparison, not an exhaustive search or a novelty assessment. It supplies no replacement theorem for the coupled five-dimensional Einstein–scalar material shell.

The four closest previously pinned numerical/boundary references have **unchanged primary arXiv histories**. Four additional papers were inspected in full primary HTML. They support careful boundary, constraint and residual accounting, but none directly certifies the registered 23-amplitude response in its pointwise first-derivative norm. The retained conditional integral-operator certificate remains the applicable route.

## Previously pinned references: no revision change found

| Primary record | Latest version shown on the retrieved primary record | Consequence |
| --- | --- | --- |
| [Dong, Mascotto and Wang, arXiv:2411.03264](https://arxiv.org/abs/2411.03264) | v3, 3 May 2026 | The earlier continuous-reconstruction comparison remains current. |
| [Dong, Georgoulis, Mascotto and Wang, arXiv:2509.08537](https://arxiv.org/abs/2509.08537) | v2, 23 April 2026 | The earlier time/space/transfer residual accounting remains current. |
| [An and Anderson, Part II, arXiv:2507.15567](https://arxiv.org/abs/2507.15567) | v3, 30 May 2026 | No later revision was shown. Its vacuum geometric boundary conditions still differ from the material shell. |
| [An and Anderson, Part III, arXiv:2503.12599](https://arxiv.org/abs/2503.12599) | v3, 18 June 2026 | The weakened injective/dense-range comparison recorded previously must not be upgraded into a general shell inverse. |

Version statements above use the primary submission-history blocks, rather than search-index update dates. Their full texts were not reread for this update; the 16 October proof-comparison remains the substantive source for them.

## Additional primary papers inspected

These papers are newly considered in this update; they were already posted before the 16 October checkpoint. No claim of a new paper appearing after that checkpoint is made.

### Robin fundamental solutions

Beatrice Costeri, Claudio Dappiaggi and Benito Alberto Juárez-Aubry, [arXiv:2609.24465v1](https://arxiv.org/html/2609.24465v1), posted 21 September 2026, *The Hadamard parametrix on globally hyperbolic spacetimes with Robin boundary conditions: Fundamental solutions and Hadamard states*.

Theorem 1.1/3.6 constructs advanced and retarded Klein–Gordon fundamental solutions on globally hyperbolic spacetimes with timelike boundary. Equation (3.3) uses an inward unit normal and a real constant Robin parameter. Remark 3.1 discusses smooth parameter extensions; that discussion is not an already checked matrix-valued coupled-shell theorem. Convexity and finite-reflection assumptions enter the later microlocal claims, not this basic Green-operator existence result. Sections 3.1 and Appendices A–B separate energy propagation, Galerkin existence and causal support.

**Use here:** independently check normal orientation, boundary-parameter signs and the distinction between existence, causal propagation and numerical error. **Transfer gap:** this is a scalar linear operator, not the coupled metric/scalar Robin system, dynamic material-shell rows, subsidiary constraints or a computable HDBLAST inverse norm. No Green kernel or numerical certificate is imported.

### Coupled heat–wave residual estimates

Pelin G. Geredeli, Mohammad Mahabubur Rahman and Deepanshu Verma, [arXiv:2610.03575v1](https://arxiv.org/html/2610.03575v1), posted 2 October 2026, *Error Estimates in Physics-Informed Neural Networks to Heat-Wave Fluid-Structure Interactive System*.

Theorem 3.3 and equations (21)–(22) estimate mixed energy/integrated errors using separate initial, PDE, outer-boundary and interface residuals. Theorem 3.7 assumes fourth-order regularity and adds quadrature-error contributions to training error. Its proof, equation (82), bounds an integral by quadrature plus a term depending on the integrand's second derivatives; sample residuals alone do not supply that bound. The theorem is stated with big-O constants and a different two-/three-dimensional heat–wave interface system.

**Use here:** maintain distinct initial, bulk, boundary and interface error budgets. **Transfer gap:** these integrated/energy estimates do not bound the registered pointwise first-derivative norm, and the assumptions/constants have not been verified for the shell problem. A PINN fit would add regularity and quadrature obligations; it would not certify a waveform merely by producing a small training loss.

### Conformal Einstein boundary constraints

Chris Stevens and Juan A. Valiente Kroon, [arXiv:2604.24931v1](https://arxiv.org/html/2604.24931v1), posted 27 April 2026, *Wellposedness of the initial boundary value problem for the conformal field equations*.

Section 5.2, Proposition 2 gives local existence for the reduced conformal evolution equations with smooth maximally dissipative data and compatible corners. Section 6, Proposition 3 then treats propagation of the constraints. Section 7 separately checks consistency of the boundary-adapted conformal-geodesic gauge and geometry. The paper explicitly distinguishes reduced evolution solvability from satisfaction of the full equations.

**Use here:** keep reduced response existence, subsidiary/constraint propagation, corner compatibility and geometric shell conditions as separate checks. **Transfer gap:** this is a different conformal Einstein system and boundary-adapted gauge, without the scalar-dependent material-shell junction law or the required numerical error estimate. Its propositions do not complete the nonlinear HDBLAST problem.

### A Robin-specific estimator in a different PDE class

Patrick Ciarlet, Minh-Hieu Do, Mario Gervais and François Madiot, [arXiv:2604.02890v1](https://arxiv.org/html/2604.02890v1), posted 3 April 2026, *A posteriori error estimates for mixed finite element discretization of the multigroup Neutron Simplified Transport equations with Robin boundary condition*.

Equation (4.9) defines a distinct Robin boundary estimator. Theorem 4.1 combines cell residual, flux and boundary contributions in its local reliability estimate. This makes the omitted-boundary-error failure mode particularly explicit.

**Use here:** retain the actual Robin defect as an error contribution rather than treating a bulk fit as sufficient. **Transfer gap:** this is a multigroup transport/diffusion mixed formulation with its own coefficient weights and norms, not the hyperbolic shell response. Neither its estimator constants nor its reliability theorem transfer without a fresh derivation.

## Consequence for the current implementation

The applicable existing model-specific identity is still

    R = V - S - T V,
    ||W - V||_X <= ||R||_X/(1-q).

For approximate inputs/operators, an exact bookkeeping decomposition is

    R = (V - S_tilde - T_tilde V)
        + (S_tilde - S) + (T_tilde - T)V.

The first term needs a continuous numerical residual bound, including integration and reconstruction. The second needs the actual source-approximation error; the retained bounds give at most `3||F_tilde-F||_C0(H) + ||G_tilde-G||_C0(H)` for this source term. The last term needs a bound for operator/coefficient approximation in the retained norm. These statements use the already accepted conditional operators, not a theorem imported from the new references.

Construct the candidate in the specified compatible function space, preserve every coupled Robin and gauge row, and check the actual initial/corner prescription. A finite Neumann partial sum need not satisfy the exact Robin law. Its source/iteration truncation bound does not by itself certify an evaluated numerical candidate. Continuous source approximation and the retained C0-to-X integral bounds avoid importing the PINN paper's additional source-derivative assumptions; the full-cell input error and coherent reconstruction still have to be supplied.

## Search and verification limits

Two Firecrawl paper-index queries covered recent wave-error/Robin and Einstein/timelike-boundary formulations. Two additional ordinary web queries restricted to primary arXiv pages in the 16–23 October window returned no results. This limited search found no directly applicable later theorem, which is not evidence that none exists.

The primary arXiv records for nine candidates were retrieved. Four candidates above were read in full primary HTML; the inverse-spectral paper [arXiv:2610.13629v1](https://arxiv.org/abs/2610.13629), posted 9 October, was screened only by abstract/metadata and excluded because it addresses wave-speed recovery from spectral data, not the present forward residual certificate. Other paper-index candidates were screened by abstracts only and not used substantively.

Primary-source metadata, exact versions and retrieval hashes accompany this bounded comparison. They document source inspection, not external peer review, originality or physical confirmation.
