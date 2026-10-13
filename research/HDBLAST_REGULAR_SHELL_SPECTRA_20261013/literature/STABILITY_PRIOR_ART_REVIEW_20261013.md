# Bounded primary-literature review of the proposed HDBLAST stability bridge

Review performed 13 October 2026, with primary retrievals approximately 06:50–07:07 UTC (late 12 October through early 13 October in America/Los_Angeles). This is a targeted primary-source comparison, not an exhaustive literature survey or a priority certificate. It evaluates the proposed contribution's relationship to existing work; it does not independently certify the new perturbation calculation.

## Scope and reproducible evidence

The review read the October 11 public literature report and the private `PROSPECTIVE_ZERO_MODE_STRATEGY.md`. It then made bounded arXiv, INSPIRE, and Crossref queries covering Einstein–scalar de Sitter branes, radion stability/ghosts, eigenparameter boundary conditions, and relevant authors. Exact query strings, result caps, timestamps, URLs and response hashes are in `SEARCH_LOG_20261013.json`, `FOLLOWUP_SEARCH_METADATA.json`, and the source receipt files. A large query's total match count is not a paper-reading count. Only returned metadata pages were inspected; no citation graph was exhausted.

Two freshness checks returned zero entries: an INSPIRE query with `date > 2026-10-10` and brane/domain-wall plus scalar/stability/radion terms, and an arXiv `lastUpdatedDate:[202610110000 TO 202610132359]` scalar/brane query. This says only that those particular API searches returned no entries at the recorded time. Indexing delays, different terminology, mathematical papers, and omitted sources remain possible. Three 2026 physics papers discussed below predate October 11; they are newly inspected here, not newly released since that date.

Twenty-two author PDFs were fetched to the private scratch cache and converted with `pdftotext -layout`. Selected passages from twenty were read; two additional mathematical PDFs were retrieved but not substantively read. The reading register distinguishes these cases and records the exact sections used. It does not claim twenty complete-paper readings. Version dates come from the arXiv API; some old PDFs display later TeX-generation dates on their first pages. All retained third-party PDF/text bodies are excluded from the proposed public deliverables.

The cloud-runtime networking instructions and `/etc/codex/network-policy.json` were consulted. The approved shell accesses used inherited proxy/CA configuration, permitted destinations, and verified TLS. Runtime observations were current but HTTP policy state was `unknown`; successful HTTP responses establish these specific readbacks, not a claim of globally enforced policy or credential readiness. No proxy/TLS bypass, publication action, or repository mutation was performed.

## Principal assessment

The coupled scalar/metric equations, the de Sitter curvature shift in the cited Einstein–scalar radion/master formulations, soft scalar junctions with the eigenvalue in the boundary condition, enlarged Hilbert spaces, zero-mode comparison, ground-state transformations, and tensor supersymmetric factorization are established. The October 11 Jacobi construction may provide a useful **specific** bridge from certified background signs to an explicitly nodeless zero-mass solution and a strict endpoint inequality on the registered regular-cone branch. This search did not identify that exact identity in the selected literature, but does not establish its absence elsewhere.

The contribution should therefore be stated as the particular identity and its proved spectral consequence under the derived physical junction/domain hypotheses. A positive mathematical norm is not by itself a physical ghost exclusion; parity and exceptional harmonics must be derived from the unreduced equations and action. A strict inequality for each fixed positive detuning is not automatically a uniform quantitative gap as detuning tends to zero.

The exact branch being compared is

\[
W(\eta)=\tfrac13+\eta^2+\tfrac13\eta^3,\quad
U=\tfrac12(W')^2-\tfrac23W^2,\quad
\sigma=2W+\delta(1+c+c\eta),\quad
c=2/(1.0357712571566784)-4/3,
\]

with the decimal interpreted rationally. The harmonic convention is \(\Box_\gamma Y=\mu^2Y\), \(m_4^2=\mu^2H_b^2\), and \(\Lambda=\mu^2+4\). Literature coordinates/actions must be translated before transferring a formula. In particular, a bulk mass parameter, de Sitter curvature, radial logarithmic derivative, and shifted spectral \(\Lambda\) are different quantities.

## Closest scalar precedents and apparent challenges

**Frolov and Kofman, “Can Inflating Braneworlds be Stabilized?”,** [hep-th/0309002v1](https://arxiv.org/abs/hep-th/0309002v1), [Phys. Rev. D 69, 044021](https://doi.org/10.1103/PhysRevD.69.044021). Read §§II–V, especially Eqs. (3)–(9), (13)–(20), and the discussion of brane displacement before the scalar junction. Their generalized longitudinal-gauge coupled Einstein–scalar system has \(\Psi=-\Phi/2\), and the radial equations contain \(m^2+4H^2\). The scalar junction is obtained before taking the stiff-potential limit; their quantitative Rayleigh estimate uses the stiff limit. Their static interbrane separation can be fixed while the radion remains tachyonic. Their statement that possible global brane displacements are set aside is a scope choice, not an exclusion of every exceptional displacement in another model. The HDBLAST master equation is not new by virtue of repackaging these equations.

**Lesgourgues and Sorbo, “Goldberger–Wise variations: stabilizing brane models with a bulk scalar”,** [hep-th/0310007v2](https://arxiv.org/abs/hep-th/0310007v2), [Phys. Rev. D 69, 084010](https://doi.org/10.1103/PhysRevD.69.084010). Read §§I–II and the zero-mode/stability discussion in §III. Eqs. (16)–(17) already give scalar junctions of the form \(g_\pm[\Psi'+2H_{rad}\Psi]+m^2\Psi=0\); the soft, eigenparameter-dependent boundary structure thus predates the 2006 work below. Their single-brane instability result is an important apparent challenge. Its setup is **flat Minkowski slices**, a compact orbifold, and a smooth second fixed point with \(a(0)>0\), \(\phi'_0(0)=H_{rad}(0)=0\), and \(V(\phi_0(0))=0\), Eqs. (18)–(19). HDBLAST instead has a vanishing cone warp factor, \(\rho\sim y\), \(a_{rad}\sim1/y\), de Sitter slices and an infinite conformal throat. Their theorem does not directly exclude this branch; the difference must be documented, rather than treating “one brane” as a sufficient dictionary.

**Konikowska, Olechowski and Schmidt, “Radion stabilization with(out) Gauss–Bonnet interactions and inflation”,** [hep-th/0603014v2](https://arxiv.org/abs/hep-th/0603014v2), [Phys. Rev. D 73, 105018](https://doi.org/10.1103/PhysRevD.73.105018). Read §§3–4, Eqs. (21)–(30), (31)–(44). Setting the Gauss–Bonnet coefficient to zero leaves the relevant Einstein–scalar precedent. Their Eqs. (23)–(25) give \(-(pQ')'+qQ=\lambda pQ\) with \(\lambda=m^2+4H^2\) and endpoint conditions \(\partial_nQ-(\lambda/b_i)Q=0\). Eqs. (27)–(29) explicitly put endpoint terms \(pQ^2/b_i\) in the **denominator** of the Rayleigh quotient. They state that inflating-brane stability requires \(\lambda_0>4H^2\). Positivity of the shifted eigenvalue alone is insufficient. Their analytic necessary/sufficient sign conditions are subsequently developed for \(H=0\), regular coefficients on a finite two-brane interval. They explicitly warn that zeros of \(\phi'\) make the divided formulation singular. The HDBLAST cone and the \(\delta=0\) constant-scalar limit require their own domain treatment.

**Olechowski, “k-stabilization in brane models”,** [0801.1605v2](https://arxiv.org/abs/0801.1605v2), [Phys. Rev. D 78, 084036](https://doi.org/10.1103/PhysRevD.78.084036). Read §§3.1–3.2 and §4, especially Eqs. (20)–(37) and (38)–(43). The paper treats generalized scalar kinetic terms and more complicated eigenvalue-dependent junctions. §3.2 constructs modified scalar products with boundary values/derivatives, reducing for ordinary kinetic terms to the endpoint-mass structure relevant here. §4 gives an explicit zero-mass reference solution and uses its derivative and boundary signs to locate massless or tachyonic modes. Neither enlarged inner products nor a zero-mode boundary comparison is a new method.

**Olechowski, “Stability of multibrane models”,** [2408.15343v2](https://arxiv.org/abs/2408.15343v2), [Nucl. Phys. B 1011 (2025), 116807](https://doi.org/10.1016/j.nuclphysb.2025.116807). Read the setup and §§4.1–4.2, Eqs. (45)–(50), with the introductory discussion of intermediate-brane junctions. It generalizes exact zero-mode and variational criteria to multibrane Einstein–scalar models. Equation (46) contains brane/interface terms in both the energy numerator and norm denominator. The additional scalar variable and junction structure of an intermediate brane differ from an orbifold fixed boundary. Its finite, flat-sliced setup supplies no ready-made regular-cone theorem, but is a strong warning against silently treating a doubled interior as only the orbifold-even problem.

## Physical kinetic norm, bending, and exceptional harmonics

**Kofman, Martin and Peloso, “Exact identification of the radion and its coupling to the observable sector”,** [hep-ph/0401189v1](https://arxiv.org/abs/hep-ph/0401189v1), [Phys. Rev. D 70, 085015](https://doi.org/10.1103/PhysRevD.70.085015). Read §III.B, Eqs. (18)–(27), together with its flat two-brane setup. The quadratic action is diagonalized using a canonical bulk variable, while boundary total derivatives are retained. A short exact passage is: “Total derivatives in y cannot be dropped”. The canonical norm in Eq. (26) contains a boundary subtraction and then becomes a positive bulk expression in the metric/scalar variables. This is direct primary precedent for distinguishing a formal Sturm–Liouville normalization from the physical kinetic coefficient. It motivates, but does not replace, the de Sitter Einstein–Hilbert/Gibbons–Hawking/junction-action calculation on the HDBLAST branch. A derived relation such as \(\mathcal K\propto\Lambda\mathcal N_{SL}\), if independently established in the present conventions, is additional physical information.

**Gen and Sasaki, “Radion on the de Sitter brane”,** [gr-qc/0011078v2](https://arxiv.org/abs/gr-qc/0011078v2), [Prog. Theor. Phys. 105, 591](https://doi.org/10.1143/PTP.105.591). Read §III, Eqs. (3.5)–(3.11), §IV, and Appendix A's two-sided thin-wall setup. Their scalar-free two-brane radion corresponds to relative displacement and has the \(-4H^2\) shift; the scalar-free single positive-tension brane has no radion. Appendix A starts from two manifolds with one common boundary and explicitly does not impose \(\mathbb Z_2\) symmetry for most of the discussion. It is a relevant source for displaced-wall and two-sided bookkeeping. It does not identify which odd modes survive the dressed HDBLAST shell. One common induced geometry, orientation signs, common versus relative bending, continuity and scalar junctions have to be derived, even if the background has an exchange symmetry.

**Kodama and Ishibashi, “A master equation for gravitational perturbations of maximally symmetric black holes in higher dimensions”,** [hep-th/0305147v3](https://arxiv.org/abs/hep-th/0305147v3), [Prog. Theor. Phys. 110, 701](https://doi.org/10.1143/PTP.110.701). Read §2.1, Eqs. (2.2)–(2.10), and Appendix B. Constant harmonics and the first nonconstant scalar harmonics make the derivative-generated harmonic structures absent or degenerate. Their generic gauge-invariant variables and equations therefore require separate treatment in these sectors. Their vacuum black-hole conclusions about gauge/parameter modes are **not** transferable conclusions about a scalar-dressed de Sitter shell. The transferable lesson is the degeneracy of the reduction, not the physical classification. HDBLAST's constant sector and translation-type \(\mu^2=-4\) sector must be checked unreduced.

**Pilo, Rattazzi and Zaffaroni, “The Fate of the Radion in Models with Metastable Graviton”,** [hep-th/0004028v2](https://arxiv.org/abs/hep-th/0004028v2), [JHEP 07 (2000), 056](https://doi.org/10.1088/1126-6708/2000/07/056). Read §§2–3, especially Eqs. (2.14)–(2.17). Their radion kinetic coefficient is computed from the action; it is negative in the GRS limit. Giving a ghost a mass does not repair its negative kinetic term. This is a concrete counterexample to identifying a stabilized modulus or real mass spectrum with ghost freedom. Its negative-tension/multiple-region geometry differs from HDBLAST.

**Gorbunov, Koyama and Sibiryakov, “More on ghosts in DGP model”,** [hep-th/0512097v1](https://arxiv.org/abs/hep-th/0512097v1), [Phys. Rev. D 73, 044016](https://doi.org/10.1103/PhysRevD.73.044016). Read the introduction and §III, Eqs. (19)–(26), including the brane-bending/spin-2 mixing at the self-accelerating degeneracy. It demonstrates why physical kinetic/Hamiltonian analysis is needed when scalar bending mixes with a tensor helicity. DGP has an induced brane Einstein term and different bulk/junction dynamics; its ghost conclusion is not a HDBLAST exclusion.

## Mathematical tools versus the particular Jacobi identity

**Shahriari, Jodayree Akbarfam and Teschl, “Uniqueness for Inverse Sturm–Liouville Problems with a Finite Number of Transmission Conditions”,** [1203.6608v3](https://arxiv.org/abs/1203.6608v3), [J. Math. Anal. Appl. 395 (2012), 19–29](https://doi.org/10.1016/j.jmaa.2012.04.048). Read §5, Eqs. (5.1)–(5.6), Lemmas 5.1–5.2. It places eigenparameter-dependent boundary data into \(L^2(w)\oplus\mathbb C^2\) with positive endpoint weights and defines the trace-linked operator domain. This directly demonstrates that augmented Hilbert-space realization is standard. Its regular finite interval is not the singular cone domain; checking HDBLAST's closed form, trace identification and limit-point endpoint remains necessary.

**Pinchover and Tintarev, “A ground state alternative for singular Schrödinger operators”,** [math/0411658v1](https://arxiv.org/abs/math/0411658v1). Read §1's operator/form assumptions and §2, Lemma 2.4, Eq. (2.6). Their familiar ground-state identity turns the quadratic form of \(\psi v\) into a positive derivative integral weighted by \(\psi^2\), under the stated boundary hypothesis. Calling the corresponding one-dimensional calculation a Picone identity or ground-state transform does not confer novelty. Retaining the nonzero shell term and proving its sign for this model is the substantive application.

**Kostenko, Sakhnovich and Teschl, “Commutation Methods for Schrödinger Operators with Strongly Singular Potentials”,** [1010.4902v2](https://arxiv.org/abs/1010.4902v2), [Math. Nachr. 285 (2012), 392–410](https://doi.org/10.1002/mana.201000108). Read §§3.1–3.2, Theorem 3.1, Hypothesis 3.2 and Lemma 3.3. Factorization, equivalence of nonzero partner spectra, and the need to match operator domains/endpoint conditions are established. Their assumptions cannot simply be imported to a boundary-mass Sturm–Liouville problem; they provide the mathematical comparison, not the missing physical domain proof.

The companion **“Weyl–Titchmarsh Theory for Schrödinger Operators with Strongly Singular Potentials”**, [1007.0136v3](https://arxiv.org/abs/1007.0136v3), [IMRN 2012, 1699–1747](https://doi.org/10.1093/imrn/rnr065), §2 setup, also explicitly separates the differential expression from its self-adjoint realization and specified endpoint conditions. Only this setup was inspected, not a claimed theorem proving HDBLAST limit-point behavior. The project's local reduction-of-order calculation at \(\Lambda=4\), with \(P_0\sim y^4\) and the second solution \(\sim y\), is the appropriate direct endpoint evidence; the Weyl alternative and form-representation machinery are standard.

The specific proposed HDBLAST identities are instead

\[
X_0=h+p\rho^2j,\quad\Psi_0=v+a\rho^2j=-\rho^2\mathcal D/12,
\quad M_0=hJ+B_\partial p\rho^2j,
\]

and \(M_0/X_0>J\ge w/2\), with independently proved pointwise Jacobi signs. These give a nodeless \(\Lambda=4\) reference with a strict shell remainder. The background junction determinant alone does not prove nodelessness, physical kinetic positivity or parity completeness. The audited scalar ground-state identity can then supply a spectral exclusion, but the method itself is standard.

## Tensor factorization and de Sitter threshold

**DeWolfe, Freedman, Gubser and Karch, “Modeling the fifth dimension with scalars and gravity”,** [hep-th/9909134v4](https://arxiv.org/abs/hep-th/9909134v4), [Phys. Rev. D 62, 046008](https://doi.org/10.1103/PhysRevD.62.046008). Read §6, Eqs. (62)–(72). The transverse-traceless graviton equation has potential \(9A_z^2/4+3A_{zz}/2\), factors into first-order partners, and has zero mode \(e^{3A/2}\). The authors explicitly verify the norm required for the positivity argument and state that they have not thereby excluded coupled scalar tachyons.

**Garriga and Sasaki, “Brane-world creation and black holes”,** [hep-th/9912118v2](https://arxiv.org/abs/hep-th/9912118v2), [Phys. Rev. D 62, 043523](https://doi.org/10.1103/PhysRevD.62.043523). Read §II, Eqs. (2.12)–(2.14) and the immediately following spectrum discussion. For doubled AdS interiors and a de Sitter brane they find the normalizable massless graviton and a volcano-potential tail \(9/4\), with KK continuum beginning at \(m=3H/2\). Neither the massless graviton nor this threshold is new.

The possible HDBLAST-specific tensor step is the exact background substitution

\[
V_-=S^2-S'=\tfrac94+\rho^2(9p^2/16-U/8)\ge\tfrac94
\quad (U\le0),\qquad S=\tfrac32(\log\rho)_z,
\]

and its application throughout the certified tube. A claim that only the even massless graviton lies below threshold also needs the actual even/odd shell domains: under the proposed map, even Neumann tensor modes map to a Dirichlet partner, whereas odd Dirichlet modes map to a Robin partner. The sign of the Robin form term must be verified in the chosen outward-normal/conformal-coordinate convention. This search found the underlying factorization and gap, not an explicit prior proof of this exact \(U\le0\) branch comparison. This is no priority certificate for the inequality or parity argument.

## Recent work and an interpretation discrepancy

**Dong, Guo, Wei and Lu, “Complex scalar field thick branes: stability of linear perturbation and evolution of scalar Kaluza–Klein modes coupled with gravity”,** [2609.21421v1](https://arxiv.org/abs/2609.21421v1), submitted 18 September 2026; no DOI in retrieved arXiv metadata. Read §§3.1, 3.1.2 and 3.3, Eqs. (3.18)–(3.35), (3.55)–(3.60). It studies smooth complex-scalar thick branes, obtains a coupled scalar master problem, and numerically reports a positive lowest scalar eigenvalue for sampled de Sitter parameters. It repeats the standard tensor factorization and localized zero-mode result. Its absence of a thin shell and its two-component scalar background prevent transfer of either spectrum or shell norm to HDBLAST. It is evidence that scalar stability on de Sitter branes is an active, established model-dependent problem, not evidence for a universal no-go or universal stability theorem.

**Karmakar and SenGupta, “Tachyonic (In)stability in Randall–Sundrum Braneworld Scenarios”,** [2605.18403v3](https://arxiv.org/abs/2605.18403v3), submitted 18 May 2026, version updated 17 September; [10.1103/xv43-tr6m](https://doi.org/10.1103/xv43-tr6m). Crossref identifies the publication as *APS Open Science* 1, 000133, 15 September 2026, with the journal title “Tachyonic instability in Randall-Sundrum braneworld scenarios”. The text read was arXiv v3, not a separate journal full-text readback. Read §4.5.1, Eqs. (105)–(112), in the two-brane low-energy scalar–tensor context. Equation (108) is

\[
\delta\ddot\phi+3H\delta\dot\phi+\mu_{eff}^2\delta\phi=0.
\]

The paper rescales \(\delta\phi=a^{-3/2}u\); its Eq. (110) gives \(\ddot u+(\mu_{eff}^2-9H^2/4)u=0\) for constant \(H\). A short exact context quotation is: “If \(\Omega^2<0\), then \(u\sim e^{\Omega t}\) and the perturbation grows exponentially”. It then calls \(\mu_{eff}^2<9H^2/4\) an instability condition, Eq. (112). The equation context presents a convention/interpretation discrepancy: growth of the rescaled variable does not by itself establish growth of the original perturbation.

Direct algebra from their Eq. (108) gives

\[
\delta\phi\propto\exp\!\left[\left(-\tfrac32\pm
\sqrt{\tfrac94-\mu_{eff}^2/H^2}\right)Ht\right].
\]

For \(0<\mu_{eff}^2<9H^2/4\), both physical-field exponents are negative. The exact example \(\mu_{eff}^2=2H^2\) gives \(\delta\phi\propto e^{-Ht},e^{-2Ht}\) while \(u\propto e^{Ht/2},e^{-Ht/2}\). Growth of this original homogeneous field requires \(\mu_{eff}^2<0\); at zero mass the constant solution is marginal. A different physical observable or norm would require its own justification. This review does not declare the whole paper invalid. Its stated threshold must not be imported as a criterion excluding every positive-mass complementary-series scalar. The HDBLAST \(\mu^2>0\) result need not exclude all masses below \(9/4\).

**Rubin, “Continuum Landscape of stable extra-dimensional metrics”,** [2610.08006v1](https://arxiv.org/abs/2610.08006v1), submitted 6 October 2026; no DOI in retrieved arXiv metadata. Read §§3–5, especially Eqs. (27)–(37). It treats quadratic-curvature modified gravity and a breathing-mode low-energy action. It separates potential stationarity, positive kinetic coefficient and positive second derivative, and uses the implicit function theorem to propagate strict inequalities in a family of metrics. Its own §4 limits the conclusion to radion stability; it is not a proof of the entire higher-derivative spectrum. It supports careful scope language, not a direct comparison theorem for a five-dimensional two-derivative Einstein–scalar cone.

The mathematical query also returned 2026 papers [2603.14817v4](https://arxiv.org/abs/2603.14817v4) and [2606.31623v3](https://arxiv.org/abs/2606.31623v3), the latter updated 4 October, concerning root functions/minimality for linear eigenparameter boundary conditions, and [2502.18289v3](https://arxiv.org/abs/2502.18289v3) on inverse-problem stability. These were inspected as **metadata/abstracts only**, not used as theorem sources. They do not justify saying that no recent mathematical work exists. They also illustrate that inverse-spectral “stability” and Lorentzian dynamical stability are different uses of the word.

## Defensible claim and remaining boundary

A defensible formulation is: “For the specified rational-coupling Einstein–scalar regular-cone branch at each positive detuning in the certified interval, background Jacobi identities construct a nodeless zero-mass scalar reference and a strictly positive shell remainder. Together with the independently derived perturbation junctions and self-adjoint domain, this yields a positive scalar spectral bound in the covered sector.” If the independent action and odd/exceptional analyses succeed, add their precise results explicitly. Cite the older equations, soft-boundary norms and zero-mode methods, and state which part of the doubled physical system is covered.

This bounded review supplies no global novelty, first-ever, nonlinear stability, cosmological-origin, experimental-evidence, tunnelling, uniform-detuning-gap, or journal-significance claim. The October 11 review's Hamilton–Jacobi and inverse-boundary scope distinctions continue to apply. The useful advance is the model-specific proved identity and audited consequences, with its actual hypotheses and interval, rather than the names of the mathematical tools.
