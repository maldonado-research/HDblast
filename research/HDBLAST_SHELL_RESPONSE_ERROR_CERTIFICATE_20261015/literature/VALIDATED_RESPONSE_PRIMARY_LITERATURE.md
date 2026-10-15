# Literature comparison for a validated shell-response calculation

Research date: 15 October 2026. This is a bounded search: three paper-index searches, two web searches, three primary content reads and five canonical metadata lookups. It is not an exhaustive review. Search excerpts, abstracts and primary passages are distinguished in the accompanying ledger. No downloaded code was executed, and no raw paper bodies are included here.

## Finding

The closest sources support two useful choices: retain the boundary mass in the energy space, and bound the residual of a continuous reconstruction. None of the sources read supplies a directly transferable, fully explicit interval certificate for all of

\[
w_{tt}-w_{zz}+W(t,z)w=0,\qquad
-w_z(t,0)+\mu w_{tt}(t,0)+\beta(t)w(t,0)=0,
\qquad \mu>0,
\]

with a certified error in the normal-derivative shell output. This is a statement about this search, not a novelty claim or a claim that such literature does not exist.

## Primary comparisons

**Boundary mass and the force sign.** [Hipp and Kovács (2019), §2.2.1, equations (8)–(10)](https://arxiv.org/html/1901.01792#S2.E8) treat a bulk wave coupled to a boundary wave with positive mass. Their variational mass form includes both the bulk integral and a positive boundary contribution. At the left endpoint of a one-dimensional collar, the outward normal derivative is \(\partial_n=-\partial_z\); their equation \(\mu u_{tt}+\kappa u+\partial_nu=f_\Gamma\), without a boundary Laplacian, has the same normal-force sign as the displayed HDB equation. This motivates retaining \(\mu |w_t(t,0)|^2\) in a combined energy.

The correspondence has limits. Their boundary-Laplacian coefficient called \(\beta\) is **not** the HDB boundary stiffness \(\beta(t)\); their \(\kappa\) is the analogous stiffness. Their prototype assumes positive constant stiffness and a coercive, time-independent bilinear form; the abstract hypotheses appear in §2.1, equations (1)–(3). A time-dependent, possibly indefinite \(W\) or \(\beta\) requires a separate argument. Their error analysis does not establish the present derivative-output certificate.

**Residuals must belong to a continuous reconstruction.** [Chaumont-Frelet and Ern (2024)](https://arxiv.org/abs/2403.12954) derive reliable damped energy-norm error estimates for finite elements in space and leapfrog in time, using a globally \(C^2\), piecewise quartic time reconstruction. The continuous residual and reconstruction correction both enter the estimate; see §3, §4 and equations (1.7)–(1.8). The leading displayed estimate has higher-order terms. Remark 4.5 discusses computable bounds for those terms: they cannot simply be omitted from a numerical certificate.

The paper's norm is a damped time integral of bulk velocity and spatial-gradient errors, with Dirichlet boundary conditions, stated source regularity and an initial-data setup different from this shell problem; see §2.1. Its fixed mesh and constant time-step analysis has a CFL condition for the leapfrog results. We use it as a methodological comparison, not as a theorem covering HDB. A residual that vanishes only at grid points, or a sequence of apparently converging grids, does not provide the continuous residual bound needed here.

**A reliable estimator can still leave unevaluated constants.** [Chaumont-Frelet, Gimperlein, Labarca-Figueroa and Nick (2025)](https://arxiv.org/abs/2509.04217) prove residual-based estimates for exterior acoustic Dirichlet scattering, using boundary integral equations and convolution quadrature. Their §1.1 specifies that different problem; Theorem 1 in §4.3 retains a temporal error term and constants depending on geometry, final time and data regularity. These are valid analytical estimates in their setting, but they are not numerical enclosures for our shell trace. An unspecified constant or a neglected temporal remainder cannot be promoted to a rigorous error bar.

## Recent candidates, read only at abstract/metadata level

Two 2026 results were found and canonically inspected:

- [Hu, Geng and Wang, arXiv:2603.00594](https://arxiv.org/abs/2603.00594), *An adaptive integrating factor midpoint method for second order evolution equations*. The abstract reports energy-based residual estimates and quadratic continuous time reconstruction. Its body was not read, so dynamic-boundary applicability is unresolved here.
- [Cao, Dörich and Li, arXiv:2607.18320](https://arxiv.org/abs/2607.18320), *Finite element exponential integration for rough solutions of nonlinear wave equations. Part II: Dynamic boundary conditions on curved domains*. The abstract reports low-regularity convergence of a displacement-velocity pair in weak norms. This does not by itself bound our normal derivative or certify roundoff.

The 2024 and 2025 primary reads above provide the more directly inspected recent methodological support. None is evidence for a higher-dimensional origin of the Big Bang.

## Consequences for the proposed HDB certificate

The following are requirements of a new HDB analysis and implementation, rather than conclusions imported from these papers:

1. Specify a smooth approximate field on continuous cells and evaluate both its bulk residual and its dynamic boundary residual. Include initial errors and any auxiliary boundary contribution affecting the observation region.
2. Enclose the actual retained background, including every required derivative, and propagate the coefficient errors in \(W,\mu,\beta\) into the residual. A floating-point background solve without a proven enclosure is insufficient.
3. Prove stability in a norm containing the positive boundary mass. If the potential or stiffness can be negative, use a justified shifted/coercive energy or a separate explicit growth bound; do not assume positivity from \(\mu>0\) alone.
4. Prove an estimate for the **output actually used**. A basic energy controlling \(w_t,w_z\) in bulk \(L^2\) and the boundary velocity need not control the pointwise normal derivative \(w_z(t,0)\). The boundary identity expresses it through \(\mu q_{tt}+\beta q\), where \(q=w(t,0)\); an additional time-regularity estimate or a characteristic boundary estimate is needed to bound its error.
5. If differentiating the equation in time to obtain that regularity, include \(W_t w\), \(\beta_t q\), differentiated residuals, compatible initial derivatives, and coefficient-enclosure errors. No derivative may be assigned a numerical bound solely because the undifferentiated quantity is bounded.
6. Enclose every residual supremum or integral between sampling points and include arithmetic, quadrature and reconstruction errors. Directed rounding or exact rational inequalities must cover the final numerical certificate, not merely the plotted waveform.
7. Only a lower bound for the response integral that remains positive after the proven output-error bound establishes a numerical positive-response certificate. Otherwise report a conditional analytic bound or an uncertified numerical experiment accurately.

The sibling mathematical derivations must stand on their own assumptions and proofs. This literature file neither proves their estimates nor verifies an implementation.

## Reproducibility and limitations

[VALIDATED_RESPONSE_SOURCE_LEDGER.json](VALIDATED_RESPONSE_SOURCE_LEDGER.json) records the five compared sources, exact locators, evidence level and transfer limits. [BOUNDED_SEARCH_RECEIPT.json](BOUNDED_SEARCH_RECEIPT.json) records all 13 research operations and discovery results. Full paper passages were inspected through Firecrawl for the 2019 and 2024 papers; the 2025 arXiv HTML was read directly through Firecrawl. Raw passages and full HTML are intentionally absent from this public-ready folder. No DOI, peer-review status or theorem applicability is inferred when it was not verified.
