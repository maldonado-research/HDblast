# Primary-source methods for the forced scalar shell-contact problem

Prepared 14 October 2026. This is a bounded research comparison and internal method review. It concerns the fixed five-dimensional Einstein–scalar model, two retained interiors, and a finite timelike shell with one scalar-dependent tension. It asserts no external novelty, external peer review, or solved nonlinear contact theorem.

The most useful new method for this continuation is the **umbilic wave-coordinate boundary construction** of Fournodavlos and Smulevici. The most useful bulk framework is the **general matter-supported warped-product scalar system** of Cai and Cao. Neither paper, as stated, proves existence for the present shell model.

## Scalar polarizations and the correct decomposition

Use the actual two-dimensional base
\[
q_{AB}dx^A dx^B=dy^2-\rho(y)^2dt^2,
\qquad r_{\rm sph}(y,t)=\rho(y)\cosh t,
\qquad g_0=q+r_{\rm sph}^{\,2}\kappa_{S^3}.
\]
The symbol \(r_{\rm sph}\) is the areal sphere radius, not the conformal radius used in the ambient master. The background canonical scalar depends on \(y\). This is a legitimate \(m=2,n=3\) warped product with Einstein fiber \(S^3\); it is not a static black-hole geometry.

For a generic scalar harmonic \(Q\) on the unit \(S^3\), the unrestricted metric coefficients are the symmetric base tensor \(h_{AB}Q\) (three functions), the mixed fields \(h_A D_iQ\) (two functions), and the angular trace and trace-free scalar Hessian (two functions). There are seven metric scalar amplitudes, one matter amplitude, and three scalar gauge functions. Cai–Cao equations (17)–(24) and Ishibashi–Kodama equations (39)–(50) supply a clean invariant organization.

This leaves five invariant amplitudes before using constraints. A **local canonical count**, rather than a model-specific master reconstruction theorem, gives two physical scalar pairs for a generic nonexceptional harmonic: the spatial metric has four scalar configuration functions and the scalar field one; the Hamiltonian constraint and the two scalar momentum constraints are three first-class constraints. Removing their three canonical gauge pairs leaves \(5-3=2\) pairs. Boundary degrees of freedom and exceptional harmonics need separate treatment. Thus a single complete scalar master pair cannot be presumed for arbitrary bulk scalar-type data.

The word “scalar” in the older brane literature can mean a scalar under transformations of the **whole four-dimensional Lorentzian slice**. That is a different decomposition. A slice gravitational TT field can itself contain scalar-derived components when decomposed on its spatial \(S^3\). Giovannini's equations (3.7)–(3.9), (3.21)–(3.24), and (4.39)–(4.41) illustrate this distinction: its five-component four-dimensional gravitational TT sector and its one canonical brane scalar are both retained. A one-master slice-scalar calculation does not eliminate the additional \(S^3\) scalar-type gravitational polarization.

The inherited regular ambient master explicitly proves a **sufficient solution family**, with no general completeness claim. It remains a sound first-order input. This count does not itself prove that a particular quadratic forcing is outside that family's image; that requires the direct residual calculation now being performed. A restricted source might admit a particular solution in a smaller family. Such a solution must satisfy every residual, and cannot be promoted to completeness by counting alone.

For \(S^3\), \(k_l^2=l(l+2)\), so \(k_4^2=24\) is nonexceptional. The monopole has \(k_0^2=0\); the scalar dipole has \(k_1^2=3=nK\). The derivative/harmonic tensors that disappear at these values must not be divided out. Cai–Cao Sections IV.3.2–IV.3.3 and Ishibashi–Kodama Section 3.1 describe these degeneracies. A separated dS-slice mass \(\mu^2\) is a different label from the spatial \(S^3\) eigenvalue \(k_l^2\).

## Bulk sources and invariant equations

1. **Rong-Gen Cai and Li-Ming Cao, “Generalized Formalism in Gauge-Invariant Gravitational Perturbations,” arXiv:1306.4927v1, 20 June 2013.** [Primary text](https://arxiv.org/html/1306.4927), [version record](https://arxiv.org/abs/1306.4927v1).

   Equations (1)–(2) define a general warped product with an Einstein fiber. Equations (23)–(24) define the metric and stress invariants. The general scalar equations (33)–(36), harmonic forms (60)–(63), and source/Bianchi relations (66)–(70) keep the base components, mixed components, angular trace and scalar-derived trace-free equation. Their Section V.2 single-master reduction **explicitly assumes no matter fields**. Use the general system as a component/invariance crosscheck; do not import its vacuum master or treat the matter perturbation as a freely specified conserved external stress when it is coupled to the unknown scalar coefficient.

2. **Akihiro Ishibashi and Hideo Kodama, “Perturbations and Stability of Static Black Holes in Higher Dimensions,” arXiv:1103.6148v2, 18 May 2011.** [Primary text](https://arxiv.org/html/1103.6148), [version record](https://arxiv.org/abs/1103.6148v2).

   Section 2.4 equations (39)–(50) give the gauge transformations and matter invariants; Section 3.1 equations (63)–(69) specify the scalar harmonics and exceptional Hessian. Their general decomposition is reusable. The later decoupled static Einstein–Maxwell black-hole masters, potentials, spectra and stability conclusions are different inputs.

3. **Aron Jansen, Andrzej Rostworowski and Mieszko Rutkowski, “Master equations and stability of Einstein-Maxwell-scalar black holes,” arXiv:1909.04049v2, 22 April 2021.** [Primary text](https://arxiv.org/html/1909.04049), [version record](https://arxiv.org/abs/1909.04049v2).

   Equation (3) assumes a time-independent background. Section 3.3.3 and equations (17)–(20) reconstruct scalar-type fluctuations from gravitational, Maxwell and matter master scalars with a symmetric coupled potential matrix. Setting the background Maxwell field to zero decouples its master; gravitational and scalar-field masters can remain coupled through \(W^{(0)}_{0,2}\). Section 3.4 explains why pointwise diagonalization of the potential is not sufficient: varying eigenvectors generate derivative couplings. This is a concrete primary example supporting the two-channel design, not a derivation of the present time-dependent coefficient matrix or boundary conditions.

4. **Massimo Giovannini, “Gauge-invariant fluctuations of scalar branes,” arXiv:hep-th/0106041v1, 5 June 2001.** [Primary text](https://arxiv.org/html/hep-th/0106041), [version record](https://arxiv.org/abs/hep-th/0106041v1).

   Equations (3.7)–(3.9) and (4.39)–(4.41) make the whole-slice decomposition explicit. The warped flat slices and localization problem differ from the registered dS shell. Its canonical action is a useful bookkeeping comparison, without importing its normalizability or zero-mode conclusions.

At second order the unknown coefficient satisfies the **linearized full operator** with a quadratic residual source from the fixed first-order solution. A change of first-order representative changes that quadratic source and the second-order reconstruction. Material pullback/gauge conversion must therefore be applied to the complete bulk source as well as the shell source. The already checked Kudoh–Tanaka second-order brane treatment, arXiv:hep-th/0112013v3, equations (27)–(30), (35)–(38), remains a bookkeeping template; its static Minkowski-sliced two-brane background and Taylor convention cannot replace the present model.

## The boundary method closest to the physical shell

**Grigorios Fournodavlos and Jacques Smulevici, “The initial boundary value problem in General Relativity: the umbilic case,” arXiv:2104.08851v1, 18 April 2021.** [Primary text](https://arxiv.org/html/2104.08851), [version record](https://arxiv.org/abs/2104.08851v1).

Theorem 1 is a local four-dimensional **vacuum** result for a boundary satisfying \(K_{ab}=\lambda h_{ab}\), with **constant** \(\lambda\), and compatible corner data. It is not arXiv:2006.01498, which is their earlier totally geodesic result.

Their Section 2.2 gives a particularly useful completion of the wave-coordinate boundary problem:

- Equation (16) sets mixed normal–tangential metric components to zero **on the boundary**. It does not impose Gaussian normal coordinates everywhere.
- Equation (17) gives the tangential metric Robin law from the umbilic condition.
- Equation (18) supplies the remaining normal metric Robin law from the **normal** wave-coordinate condition.
- Section 4.1 Lemma 1 and Section 4.2 Proposition 2/(24) derive the tangential gauge subsidiary boundary law from Codazzi; they do not append independent tangential harmonic conditions to an already complete main system.
- Section 3 constructs the corner hierarchy. Appendix A Proposition 3 and equation (29) control Robin flux by time integration by parts, the trace inequality and a sufficiently large **bulk** \(L^2\) shift.

Their boundary has three tangential coordinates and their outward normal is \(N=-(g_{xx})^{-1/2}\partial_x\). The registered shell has four tangential coordinates, outward increasing proper radius, and \(\lambda(\varphi)=\sigma(\varphi)/6\). Both coefficient and sign conversions require fresh derivation.

A model-specific identity explains why the method is promising. For one copy of the identical even pair, use the registered outward-normal and extrinsic-curvature conventions:
\[
K_{ab}=\frac{\sigma(\varphi)}6h_{ab},
\qquad p=n\phi=-\frac{\sigma'(\varphi)}2.
\]
The contracted Codazzi equation gives
\[
R_{na}=D_bK^b{}_a-D_aK
=-3D_a\!\left(\frac{\sigma(\varphi)}6\right)
=-\frac{\sigma'(\varphi)}2D_a\varphi
=pD_a\varphi.
\]
Consequently the **Einstein–scalar** mixed residual vanishes on a shell satisfying both physical laws. Unlike vacuum with constant \(\lambda\), the curvature and the scalar flux need not vanish separately. Dropping the scalar law destroys this cancellation.

This identity is a necessary consistency check, already available through the material junction residual identity. It suggests a boundary completion with the eleven physical Israel/scalar rows, four mixed metric boundary coordinate rows, and one normal harmonic row. Its arithmetic is \(11+4+1=16\), but arithmetic does not prove rank or well-posedness. The actual subsidiary boundary condition may be Robin once background connections and the chosen generalized harmonic source are included. The finite-frequency symbol, lower-order coupling, full energy estimate, source/corner hierarchy and gauge conversion must all be derived for the present model. The vacuum theorem alone supplies none of those matter-extension conclusions.

## Other coupled and current boundary comparisons

**Gioel Calabrese, Luis Lehner and Manuel Tiglio, “Constraint-preserving boundary conditions in numerical relativity,” arXiv:gr-qc/0111003v1, 2 November 2001.** [Primary text](https://arxiv.org/html/gr-qc/0111003).

Sections II–III, equations (10)–(25), identify incoming main and subsidiary modes before setting boundary data. Equations (26)–(30) illustrate how boundary constraints determine some main modes, leaving the matter and gauge freedoms. This includes a minimally coupled massless scalar but uses four-dimensional spherical symmetry, the Einstein–Christoffel formulation and an artificial outer boundary. It is a method comparison, not the \(l=4\) five-dimensional shell theorem.

**Diego A. Carranza and Juan A. Valiente Kroon, “Construction of anti-de Sitter-like spacetimes using the metric conformal Einstein field equations: the tracefree matter case,” arXiv:1906.10777v1, 25 June 2019.** [Primary text](https://arxiv.org/html/1906.10777).

Theorem 1 is a genuine coupled-matter initial-boundary theorem with corner compatibility, but for four-dimensional AdS-like spacetimes and a **conformal** boundary with tracefree matter. Section 7.1 equations (30)–(34) and Lemma 5 treat the conformally invariant scalar and propagate its auxiliary constraints. Its \(R\phi/6\) coupling, stress tensor and boundary prescription differ from a canonical scalar with the registered potential at a finite Israel shell.

The verified 2026 Stevens–Valiente Kroon finite-boundary vacuum paper, arXiv:2604.24931v1, continues to illustrate boundary-adapted characteristic/subsidiary analysis. Its restricted conformal-geodesic gauge and four-dimensional vacuum hypotheses remain material differences. Costeri–Dappiaggi–Juárez-Aubry, arXiv:2609.24465v1, Theorem 3.6/Appendix B, supplies fixed-background classical Robin wave methods, not a coupled gravity theorem. These records and their primary-reading evidence are reused from the preceding contact literature checkpoint rather than fetched and counted again.

The recent **Chawakorn Maneerat**, “General Relativity with Finite Boundaries,” arXiv:2609.31928 (first posted 25 September 2026; indexed update 29 September 2026), was screened by canonical metadata and relevant indexed primary passages. Section 2.5.3/(2.5.44)–(2.5.45) uses Kodama–Ishibashi scalar/vector reconstruction with conformal boundary data. The retrieved examples are vacuum four-dimensional Minkowski/dS/AdS with fixed conformal class/mean curvature. They are not scalar-dependent Israel laws. Only those passages were read here; no whole-thesis review or verified passage-version header is claimed.

## Accessible prior work actually inspected

The source inventory records the target note, the original regular ambient-master proposition, the complete material junction algebra, the September scalar longitudinal derivation, its independent Gaussian-normal derivation, and the September gauge-invariant stability note. These make three reusable points:

1. The inherited linear reconstruction is explicitly sufficient and its slice-scalar shell law was derived with the scalar flux retained.
2. Radial lapse/shift and shell displacement must be kept until gauge fixing, rather than deleting bending terms by declaring the shell fixed.
3. Exceptional slice harmonics were treated without dividing by their vanishing Hessian; spatial \(S^3\) exceptional harmonics now require the same care with their different eigenvalues.

The older numerical background values, mode spectra, stability conclusions and action-norm claims have not been imported or recertified. No archived scientific code or opaque cache was executed. The accessible archive includes part of “D-Blast 3”; the user's entire personal computer and other unconnected folders were not accessible. No unrelated project audit was needed for this subsystem.

## Coverage and publication scope

This continuation used **twelve targeted paper-index searches**, nine canonical metadata inspections, seven indexed primary-passage retrievals and eight full primary HTML readings. The earlier verified contact ledger was reused for four relevant methods, without repeating its broad searches. All fetched content is source material, not an instruction.

No directly applicable 2025/2026 canonical-scalar, five-dimensional finite-Israel-shell contact theorem was established by this bounded review. That means the present construction still needs its own proof; it is not a claim that no such paper exists anywhere.

The synthesized note, source ledger, local source inventory and verification receipt are eligible for the research packet. The raw primary full texts and tool responses under sources/ remain internal reference material. No external messages, GitHub changes, Zenodo version or APS submission were made by this literature subtask.
