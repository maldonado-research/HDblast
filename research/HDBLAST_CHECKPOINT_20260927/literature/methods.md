# HDBLAST literature sweep 3: mathematical and numerical methods

Research round of 27 September 2026, written 28 September 2026. **Work in progress.** This is an annotated bibliography with a methods map and one small validated-numerics pilot. Nothing here is a claim of the HDBLAST programme until it appears in a dated checkpoint or a Zenodo record.

## Read this first: the web search did not run

**This sweep could not be carried out as specified (procedural negative result).**

- The task required at least 15 distinct web searches. **None was executed.** All four attempted queries were refused by the tool: the session's web-search budget (200 of 200 calls) had already been used by earlier work in the same session. The refusal notice asked to continue with information already gathered, so no other search route was tried.
- A page fetch was also blocked: `journals.aps.org` returned `EGRESS_BLOCKED`.
- Attempts, verbatim:

- `Z4 formulation constraint damping Gundlach Calabrese Hinder Martin-Garcia 2005`: refused: session web-search budget used (200 of 200)
- `conformal and covariant Z4 CCZ4 formulation Alic Bona-Casas Bona Rezzolla Palenzuela`: refused: session web-search budget used (200 of 200)
- `Z4c formulation Bernuzzi Hilditch constraint-preserving spherical symmetry`: refused: session web-search budget used (200 of 200)
- `generalized harmonic formulation constraint damping Pretorius Lindblom Scheel Kidder 2006`: refused: session web-search budget used (200 of 200)
- fetch `https://journals.aps.org/prd/abstract/10.1103/PhysRevD.44.1007`: EGRESS_BLOCKED (journals.aps.org blocked by the network egress proxy)

Because of this, the bibliography below is assembled only from sources that are already on record. It does **not** reflect a fresh search of the 2020–2026 literature. Two kinds of source were used, and every entry says which one it is:

1. **Sweep-1 search records.** The URL was returned by a web search earlier in this session, during sweep 1 (theory), and is stored in `theory_checks/theory_sources.json`. The finding here re-reads that record's snippet paraphrase from a methods angle. The evidence is therefore **second-hand search-snippet evidence**.
2. **Programme records.** The URL is written verbatim in a file of the programme's own archive (`new-files/D-Blast 3/…`). An earlier project agent reports having opened it on the stated date. The finding paraphrases that record only. **This sweep did not open, search or re-verify these sources.**

The builder (`methods_checks/build_methods_md.py`) machine-checks each URL against its claimed origin. It refuses to build if a URL is not byte-identical to its sweep-1 record, or not present verbatim in the named programme file.

A third list, **leads**, gives standard methods references from general knowledge. The task names several of them: CCZ4, Z4c, generalized harmonic, hyperboloidal compactification, the CDL negative-mode papers, CosmoLattice, and validated ODE libraries. They carry **no URL**, were **not verified**, may contain bibliographic errors, and **are not cited as sources**. They are paired with 20 prepared queries for a later run.

**Totals:** 42 sources with URLs (26 from sweep-1 search records, 16 from programme records); 23 unverified leads without URLs; 20 planned queries (none run); implications: constrains 5, context 6, method 31; published 2020-2026: 17, before 2020: 25.

No paper was read in full in this sweep. Every access label is either "search snippet only" (second-hand) or "programme record only".

## Question

Which mathematical and numerical methods would most directly advance the programme's open problems? The problem list below is taken from the 22 September checkpoint (Zenodo 22922928), Chat 14, the registered-shell controls, and this round's workstreams.

| Code | Open problem | Sources below |
|---|---|---|
| P1 | Rigorous existence and local uniqueness of the static +1 branch (certification) | C2, C3, C7, C10, C11, C12, E1, E2, E3, E4 |
| P2 | Stability beyond the computed sectors: negative modes, gauge modes, coupled matter | B7, C1, C2, C3, C4, C5, C6, C8, C9, C10, C13, C14, E2, E4 |
| P3 | Long 5D evolutions: constraint growth, far-boundary contamination, attraction to the +1 branch | A1, A2, A3, A4, A5, B1, B6 |
| P4 | The collapse fate: unconverged terminal collapse; horizon or crunch | B3, B4, B5 |
| P5 | Shell/junction numerics: boundary closure, moving shell, thin versus thick wall | A3, A4, B1, B2, B3, B6, B7 |
| P6 | Matter and radiation: non-perturbative particle production with consistent backreaction | C8, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12 |

## Main findings

Each finding carries a label. **exact** means verified symbolically in this folder. **numerical** means floating-point or multiprecision with stated evidence. **conditional** means it depends on a stated assumption. **judgement** means a methods recommendation, not a result.

1. **Certifying the +1 branch is feasible with methods the programme already uses (conditional, judgement).** The programme already holds a computer-assisted existence-and-local-uniqueness proof for a static shell. That proof is PHYS-M462, 17 Aug 2026, in `D-Blast 3/untitled folder 139/.../M462_…CERTIFIED_FIXED_T_CONTINUUM_BVP_ROOT.md`. It covers the original shell at t = 0.001, not the +1 branch. Its ingredients were:
   - a certified regular endpoint germ;
   - a pinned `kv` v0.4.62 validated Taylor integrator (order 12, outward-rounded binary64);
   - a defect chart that subtracts the flat-wall direction;
   - a strict Krawczyk inclusion in a two-parameter shooting box.

   Porting this pipeline to the +1 branch is the most direct route to P1. The pilot below identifies what the port needs:
   - **A local germ.** It is explicit: the regular linear cone solution is exactly η = a·C₁₄⁽²⁾(cosh u)/680, with u = y/9 (**exact**).
   - **An outward direction that is well conditioned in the linear regime (conditional).** The unwanted, cone-singular solution decays like e^{−18u} while the regular one grows like e^{14u} (**exact**). Over [0.25, u_b] the unwanted direction is suppressed by about 10^43.3, with u_b = 3.364076 at δ = 0.001.
   - **The η-formulation, which is mandatory.** The germ amplification is G(u_b) = 6.288479e+18 (10^18.799). In binary64, 1+η_h rounds to exactly 1. A direct φ-formulation would need about 39 decimal digits to carry η_h to 16 significant digits. The checkpoint already uses η = φ − 1, which removes this problem.
   - **A shooting scale.** By analogy with M462's scaled shooting variable, a natural variable of order δ⁰ is â = η_h·G(u_b)/δ. At δ = 0.001 it is -0.08406, against the leading-order value −9c/64 = -0.08404.
   - **A good centre.** The closed-form leading-order (LO) model misses the nonlinear solution only by O(δ) relative amounts (**numerical**, table below). Either the LO model, or better the O(δ⁸) series and 31–41-digit BVP solutions of this round's `analytic_structure` workstream, can serve as the approximate solution around which a Krawczyk or Newton–Kantorovich proof is built.

   Precedent: the Euclidean +1 branch is a symmetry-reduced (cohomogeneity-one) Riemannian Einstein–scalar metric on a ball with a regular centre. Wang (2026, E1) built such metrics by exactly this approximation → residual bound → fixed-point route. What is still missing is the enclosure itself: a validated integrator (none is installed here), the nonlinear remainder bound, and the Krawczyk step. **No existence proof is claimed.**

2. **Spectral statements can in principle be certified too (judgement).** In the Frolov–Kofman variables (C3) the scalar problem is a self-adjoint Sturm–Liouville problem. The Chat 9 file reproduced the registered eigenvalue μ² = −7.7178716 with them. Rigorous eigenvalue enclosures of the kind covered by verification-method texts (E2) could turn this round's floating-point result "no scalar mode with μ² < 9/4 on the +1 branch" into a certificate. Verified determinant bounds (E4) could certify the nonsingular junction Jacobian that gives local uniqueness. Separately:
   - Lorentzian linear stability does **not** settle the **Euclidean negative-mode count**. That count is what matters if the static branch is read as a Garriga–Sasaki-type creation or tunnelling saddle (C2, C10). The CDL negative-mode literature could not be retrieved here; it is listed under leads.
   - Brane matter coupled to φ can make tachyonic ladders normalizable (C8). Any matter extension therefore needs its own coupled negative-mode check.
   - The vector sector is untreated in this round. The 2026 thick-brane paper (C9) is a current example that includes it.

3. **The far boundary of the evolutions can be removed rather than pushed out (exact observation, judgement on the remedy).** In the evolution chart, dz = dy/ρ. On the φ = +1 AdS background, z = ln tanh(y/18) + const → −∞ at the regular cone, and the (t, z) part of the metric is −ρ²dt² + dy² with ρ/y → 1 (**exact**, pilot part 1). So the far end of the bulk is a Rindler-type horizon of the de Sitter slicing, which the conformal chart reaches only as z → −∞. Any finite far boundary z = −L (as in `registered_solver.py`, with L = 6 by default) is therefore an artificial truncation. Signals from such a boundary and from the seed taper could reach the shell before the latest Chat 14 values, which is why the 22 September checkpoint does not treat those values as settled.
   - The standard remedies are a horizon-penetrating (hyperboloidal-type) slice or a characteristic chart across that horizon. The task-named hyperboloidal and characteristic literature is listed under leads.
   - He–Tian–Zhang (A2) warn that convergence can be lost at a compactified edge. Any such map needs its own convergence test there.

4. **Constraint growth has two standard complements to this round's fix (judgement).**
   - The `evolution_constraints` workstream removed the stall by projecting the initial data onto the discrete constraint.
   - The remaining growth is geometric: the outgoing weighted density e^{3A}(C_H+2C_M) is conserved while the warp factor falls. Constraint-damped formulations (Z4, CCZ4 and Z4c, and generalized harmonic; leads) add explicit decay.
   - A posteriori and adaptive local-time-stepping wave methods (A3, A4) target the thin near-shell region where the error is generated.
   - These are recommendations. Nothing was implemented here.

5. **External calibrations exist and have not yet been used (judgement).** Three published or exact benchmarks would test the programme's evolution and perturbation codes against answers they did not produce:
   - one published BraneCode test (B1), for the nonlinear 5D brane-plus-bulk-scalar evolution;
   - Kraus's exact moving AdS shell (B2), for the moving-boundary closure with the scalar frozen;
   - the exactly solvable Koyama–Takahashi model (B7), for perturbations with both junctions.

   The thin-wall versus thick-wall cross-validation of B3 is the check needed before replacing the Israel shell by a resolved wall.

6. **The collapse fate needs horizon diagnostics, which need new output (judgement).**
   - Colliding-wall and AdS-bubble simulations (B4, B5) end in spacelike singularities behind apparent horizons.
   - In the 1+1 reduction a trapped region is where (∂_t+∂_z)A < 0 and (∂_t−∂_z)A < 0 both hold. Evaluating this needs velocity-field snapshots, which the Chat 14 archive lacks. New runs should store them.
   - B3 gives an explicit expansion-versus-collapse criterion to evaluate on the trajectories.

7. **Particle production with consistent backreaction: benchmarks exist; a verified 2020–2026 lattice reference does not, here (judgement, negative for the search goal).**
   - The records supply exact calibration problems: conformal non-production (D4); the solvable sech pulse (D2); Pöschl–Teller non-production points (D3); and the instant-preheating limit (D1), which this round's `preheating` workstream uses.
   - They supply one fully self-consistent reheating computation, into a holographic 5D sector (D10). This is the closest template for the bulk-coupled feedback that the preheating workstream left unsolved.
   - They supply warnings about bulk-leakage bookkeeping (D6, D7).
   - **No 2020–2026 lattice paper with renormalized backreaction could be verified in this sweep.** CosmoLattice and the adiabatic-regularization literature appear only as leads.

8. **"Gregory–Ruth" was not identified.** The task names it among the negative-mode references. No record in this session or in the scanned programme files matches that pair of names. Without search it could not be resolved, so it is not cited.

## Validated-numerics pilot (theme E, problem P1)

Script: `methods_checks/certification_pilot.py` → `methods_checks/certification_pilot.json`. It uses SymPy 1.14 and mpmath 1.3 (`mpmath.iv` at 60 digits, point checks at 100 digits), runs in 0.54 s on one core, and has script SHA-256 prefix `35362a0d7032`. The input is the checkpoint's `static_branch/PLUS_BRANCH_RESULTS.json`, read-only with its SHA-256 checked.

**What is enclosed.** The *leading-order truncated model* only:
- the thin-brane metric junction k·coth u_b = σ(1)/6, with k = 1/9;
- the linear scalar germ on pure AdS;
- the linearised scalar junction η_y = −2η_b − δc/2.

This is **not** an enclosure of the nonlinear +1 branch.

| Check | Result | Label |
|---|---|---|
| U''(1)/k² = 252, so the linearized equation is (x²−1)η'' + 5xη' = 252η with x = cosh u | holds | exact |
| C₁₄⁽²⁾(cosh u) solves it; C₁₄⁽²⁾(1) = 680; exponential form Σ(j+1)(15−j)e^{(14−2j)u} | holds | exact |
| Growth exponents at large x: 14 (regular) and −18 (unwanted) | holds | exact |
| LO scalar junction gives η_b = −(9c/64)δ at leading order | holds | exact |
| Evolution chart: dz/dy = 1/ρ, z → −∞ and ρ/y → 1 at the cone (Rindler-type horizon) | holds | exact |
| 100-digit point values inside all 30 interval enclosures | all contained | interval |
| LO metric-only H² vs checkpoint `metric_only_H2` | agree to ≤ 6×10⁻¹⁴ relative | numerical |
| LO η_b vs checkpoint `eta_finite_curvature_linear` | agree to ≤ 4×10⁻¹⁶ relative (same linear model, implemented independently) | numerical |
| Controls: Gegenbauer degree 13 or 15, index 3/2, 4D-like damping coefficient | all fail, as required | control |
| Control: degree-13 amplification vs the checkpoint's η_b/η_h | misses by > 90 % | control |
| Control: c perturbed by +1 % | detected at the three smallest δ (deviation > 10× baseline) | control |
| Control: wrong-sign scalar junction | detected: η_b = -6.715e-04, about 8× too large (see note) | control |
| Control: containment test on a deliberately shifted point | rejected, as required | control |

Comparison of the LO enclosures with the checkpoint's nonlinear floating-point solutions (numerical):

| delta | G(u_b) at checkpoint u_b | (eta_b/eta_h)/G - 1 | eta_b: LO/ckpt - 1 | eta_h: LO/ckpt - 1 | H^2: LO/ckpt - 1 | LO vs ckpt `metric_only_H2` |
|---:|---:|---:|---:|---:|---:|---:|
| 0.0003 | 2.853884e+22 | -1.762e-05 | -1.480e-05 | 5.738e-07 | 4.717e-06 | -5.6e-14 |
| 0.001 | 6.288479e+18 | -5.876e-05 | -4.934e-05 | 1.769e-06 | 1.574e-05 | 2.3e-14 |
| 0.003 | 2.937632e+15 | -1.765e-04 | -1.481e-04 | 4.076e-06 | 4.734e-05 | 4.3e-15 |
| 0.01 | 6.921797e+11 | -5.912e-04 | -4.949e-04 | -6.768e-07 | 1.592e-04 | 2.8e-15 |
| 0.03 | 3.902683e+08 | -1.798e-03 | -1.495e-03 | -1.220e-04 | 4.897e-04 | 1.1e-15 |
| 0.1 | 1.711731e+05 | -6.292e-03 | -5.103e-03 | -1.732e-03 | 1.762e-03 | 4.4e-17 |

- The transport ratio (η_b/η_h)/G(u_b) − 1 measures only the nonlinear correction to the germ between cone and shell, since it uses the checkpoint's own u_b. It scales as -0.0587·δ; the log–log slope over the three smallest δ is 1.0008.
- The LO error in η_b is -0.0493·δ (slope 1.0004).
- The LO error in H² is 0.01572·δ (slope 1.0015). This is the known −c²δ²/384 scalar-profile term of the checkpoint's expansion: the predicted coefficient is 27c²/(384(1+c)) = 0.01572, and the observed value at δ = 3×10⁻⁴ is 0.01572.
- The LO error in η_h has a small O(δ) coefficient that changes sign near δ ≈ 0.01. The coefficients by increasing δ are +0.0019, +0.0018, +0.0014, -0.0001, -0.0041, -0.0173, and the three-point slope is 0.8528. The O(δ) and O(δ²) terms partly cancel there, so this slope is less clean than the others.
- The differences between the checkpoint's default and refined solver settings are ≤ 2×10⁻¹² relative for η_h. That is far below every deviation in the table, so the table measures the model truncation, not solver noise.

**Note on a control.** The wrong-sign-junction control was first written to expect a sign flip of η_b. That expectation was wrong: with k·G_u/G ≈ 14/9 < 2 the wrong-sign denominator is also negative, so η_b keeps its sign and becomes about 8 times too large. The first run recorded "not met". The criterion was changed to a magnitude test, and both outcomes are stored in `control_design_note` in the JSON.

**Scope.** Items 1–4 of the pre-declared expectations in the script docstring were met. The pilot shows that the linear germ and conditioning are favourable. It proves nothing about the nonlinear branch: no existence, no uniqueness, no stability.

## Annotated bibliography

Every entry is based on a search snippet (second-hand, through the sweep-1 record) or on a programme record. **None was read in full in this sweep.** Authors are given only as they appear in the originating record.

### A. Constraint-damped formulations, outer boundaries and error control

**A1. Robustness of slow contraction to cosmic initial conditions** (2020)
- Authors: A. Ijjas, W. G. Cook, F. Pretorius, P. J. Steinhardt, G. N. Davies. arXiv:2006.04999.
- URL: <https://arxiv.org/abs/2006.04999>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `C7` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: numerical-relativity simulations vary all freely specifiable initial data (shear, curvature, scalar-field and velocity profiles), including data far outside the perturbative regime, and find slow contraction robust. A follow-up (arXiv:2104.12293) reports robustness enhanced with several modes and reduced symmetry. The evolution formulation and its constraint-damping choices are not stated in the record.
- Relevance to HDBLAST (P3): A template: test whether the two Chat 14 fates survive non-perturbative disturbances rather than the narrow balanced two-bump family. Arbitrary data first need a constraint solve. The programme's discrete-constraint projection (evolution_constraints workstream) is the analogue of that step.
- Implication: `method`

**A2. Characteristic evolution of conformal scattering: I. Scalar Waves in Minkowski Spacetime** (2026)
- Authors: He, Tian and Zhang (surnames as recorded). arXiv:2608.15729v1 (16 Aug 2026).
- URL: <https://arxiv.org/html/2608.15729v1>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 141/FORWARD_ERROR_CONTROL_COMPONENTS_20260905/METHODS_AND_LITERATURE_20260905.md` (primary preprint page, checked 5 Sept 2026 by an earlier agent).
- Finding: Per the programme record: studies compactified double-null scattering with alternative characteristic stencils, and reports loss of convergence near spatial infinity. Its stability, convergence and Richardson-extrapolation evidence are numerical, not interval certificates.
- Relevance to HDBLAST (P3): A compactified or characteristic outer region is one way to stop the far-boundary signals that contaminate the late Chat 14 values (H/H0 near 0.628690 arrives after taper and boundary signals could reach the shell). The reported loss of convergence at the compactified edge is a warning to test any such map with a convergence study at the new boundary.
- Implication: `method`

**A3. A posteriori error analysis and adaptivity of a space-time finite element method for the wave equation in second order formulation** (2026)
- Authors: Dong, Georgoulis, Mascotto and Wang (surnames as recorded). arXiv:2509.08537v2; Numerische Mathematik, DOI 10.1007/s00211-026-01561-3 (per programme record).
- URL: <https://arxiv.org/html/2509.08537v2>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 141/FORWARD_ERROR_CONTROL_COMPONENTS_20260905/METHODS_AND_LITERATURE_20260905.md` (preprint and institutional record, checked 5 Sept 2026 by an earlier agent).
- Finding: Per the programme record: explicit a posteriori error bounds from temporal and spatial reconstructions, including mesh-change discontinuities, for a scalar wave equation on a bounded domain.
- Relevance to HDBLAST (P3, P5): A route from 'grid convergence' to computable error estimates for wave evolutions. The HDBLAST system (coupled, nonlinear, with a shell boundary) would need its own constants. At most it can guide an error indicator near the steep wall.
- Implication: `method`

**A4. Adaptive FEM with explicit time integration for the wave equation** (2025)
- Authors: M. J. Grote, O. Lakkis, C. S. Santos. arXiv:2507.11193v3; J. Comput. Appl. Math. 481, 117272 (2026) (per programme record).
- URL: <https://arxiv.org/abs/2507.11193v3>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 136/M489G_LITERATURE_AND_PUBLICATION_REVIEW_20260903.md` (arXiv full text (html v3), read 3 Sept 2026 by an earlier agent).
- Finding: Per the programme record: combines a posteriori error indicators, evolving spatial meshes and local explicit time steps for the wave equation.
- Relevance to HDBLAST (P3, P5): The constraint workstream of this round found the outgoing constraint front is generated within about 0.02 of the shell in the first ~0.05 time units. Local refinement with local time stepping in exactly that region is the natural use. Not an error certificate.
- Implication: `method`

**A5. Damped energy-norm a posteriori error estimates for fully discrete approximations of the wave equation using C2-reconstructions with the leapfrog scheme** (2024)
- Authors: T. Chaumont-Frelet, A. Ern. arXiv:2403.12954v2 (20 Dec 2024).
- URL: <https://arxiv.org/abs/2403.12954v2>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 145 APS JOURNAL/HDBLAST-v24-APS-review-package/manuscript/HDBLAST-v24-manuscript.md` (preprint metadata and abstract, checked 9 Sept 2026 by an earlier agent).
- Finding: Per the programme record: scalar finite-element/leapfrog a posteriori estimates in a damped energy norm; the programme noted that they concern a different formulation and do not supply its constants.
- Relevance to HDBLAST (P3): A second example of energy-norm error control for explicit wave schemes. Relevance is indirect for the method-of-lines finite-difference solvers used in HDBLAST.
- Implication: `context`

### B. Thin shells, junction conditions, moving boundaries and collapse diagnostics

**B1. BRANECODE: A Program for Simulations of Braneworld Dynamics** (2004)
- Authors: J. Martin, G. Felder, A. Frolov, M. Peloso, L. Kofman (per sweep-1 record). arXiv:hep-ph/0404141; companion 'Brane world dynamics with the BraneCode', arXiv:hep-th/0309001.
- URL: <https://arxiv.org/abs/hep-ph/0404141v1>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `E12` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: a C++ code for fully nonlinear evolution of 5D braneworlds with bulk scalar fields, described as the first of its type, used to confirm the tachyonic instability of 'stabilized' de Sitter branes. The programme's Chat 9 prior-art file (which read the papers through ar5iv on 16 Sept 2026) gives hep-th/0309001 as the BraneCode reference; the two identifiers may be a program paper and a dynamics paper, which this sweep could not check.
- Relevance to HDBLAST (P3, P5): The closest published analogue of rolloff5d (nonlinear 5D evolution of a brane with bulk scalar and junction conditions). Reproducing one published BraneCode test would give the programme an external calibration that its internal controls cannot provide.
- Implication: `method`

**B2. Dynamics of Anti-de Sitter Domain Walls** (1999)
- Authors: P. Kraus. JHEP 12 (1999) 011; arXiv:hep-th/9910149.
- URL: <https://arxiv.org/abs/hep-th/9910149v1>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `D5` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: moving domain walls in the Randall-Sundrum universe with the bulk patched from AdS5 black-hole solutions; the junction equations fix the wall motion, which observers on the wall see as cosmological expansion or contraction.
- Relevance to HDBLAST (P5): An exact moving-shell benchmark. A thin shell in pure AdS5 (scalar frozen at a vacuum) evolved by the HDBLAST moving-boundary code must reproduce the exact Friedmann-type trajectory, including the a^-4 bulk-mass term when a black-hole patch is used. This tests the shell closure separately from the scalar physics.
- Implication: `method`

**B3. When do colliding bubbles produce an expanding universe?** (2003)
- Authors: J. J. Blanco-Pillado and colleagues (per sweep-1 record). Phys. Rev. D 69, 103515 (2004); arXiv:hep-th/0306151.
- URL: <https://arxiv.org/abs/hep-th/0306151>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `C10` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: a symmetric collision leaves a collapsing universe unless the bulk expansion rate before the collision is large compared with the momentum transfer in the fifth dimension; thick-wall numerical collisions confirm the thin-wall result.
- Relevance to HDBLAST (P5, P4): A published thin-wall versus thick-wall numerical cross-validation, the check HDBLAST needs if the Israel shell is ever replaced by a resolved wall. P4: an explicit expansion-versus-collapse criterion to evaluate along the Chat 14 trajectories.
- Implication: `method`

**B4. Dynamics of colliding branes and black brane production** (2007)
- Authors: Y. Takamizu, H. Kudoh, K. Maeda. Phys. Rev. D 75, 061304; arXiv:gr-qc/0702138.
- URL: <https://arxiv.org/abs/gr-qc/0702138>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `C13` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: self-gravitating colliding domain walls (initial data from a BPS wall in 5D supergravity) generically form a spacelike curvature singularity covered by a horizon, a black brane with trapped walls; non-relativistic weak-field cases can pass through or bounce.
- Relevance to HDBLAST (P4): A concrete hypothesis and a method for the unconverged Chat 14 collapse: track trapped surfaces. In the 1+1 reduction with metric e^{2B}(-dt^2+dz^2)+e^{2A}dx_3^2, a trapped region is where (d_t+d_z)A and (d_t-d_z)A are both negative. Evaluating this needs the velocity-field snapshots the Chat 14 archive lacks, so it must be added to new runs.
- Implication: `method`

**B5. Crunch from AdS bubble collapse in unbounded potentials** (2024)
- Authors: not captured in snippet. arXiv:2411.07692.
- URL: <https://arxiv.org/abs/2411.07692>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `B17` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: an AdS bubble nucleated in a Minkowski false vacuum, in a potential with an infinitely deep true vacuum, has an interior that collapses to a spacelike singularity behind an apparent horizon; trapped surfaces form once kinetic energy dominates.
- Relevance to HDBLAST (P4): The diagnostic set for the collapse fate: apparent-horizon location, trapped-surface formation, and the sign of the local energy density. The registered U is also unbounded below (sweep 1 checked U ~ -(2/27) phi^6 exactly).
- Implication: `method`

**B6. Simulating the universe(s) II: phenomenology of cosmic bubble collisions in full General Relativity** (2014)
- Authors: Johnson, Peiris, Lehner and collaborators (per query context in sweep 1). arXiv:1407.2950 (also arXiv:1112.4487, arXiv:1508.03641).
- URL: <https://arxiv.org/pdf/1407.2950>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `C17` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: full-GR simulations of bubble collisions that compute cosmological observables directly from the scalar Lagrangian, giving relativistic predictions that differ significantly from non-relativistic approximations.
- Relevance to HDBLAST (P3, P5): Established practice for evolving resolved scalar walls in full GR in a symmetry-reduced setting. It is an alternative to a thin Israel shell, and a later route from a violent 5D event to observables.
- Implication: `method`

**B7. Exactly solvable model for cosmological perturbations in dilatonic brane worlds** (2003)
- Authors: K. Koyama, K. Takahashi. arXiv:hep-th/0307073.
- URL: <https://arxiv.org/pdf/hep-th/0307073>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `E9` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: a bulk scalar with exponential potential and exponential coupling to the brane tension gives power-law inflation on a single brane; the backreacted background and the perturbations obeying the junction conditions are solved exactly.
- Relevance to HDBLAST (P2, P5): An exact benchmark with a bulk scalar, a moving de Sitter-like brane and both junction conditions. It is the cheapest external test of the full perturbation and junction machinery.
- Implication: `method`

### C. Gauge-invariant perturbations and negative modes of de Sitter-sliced walls

**C1. Perturbations on domain walls and strings: A covariant theory** (1991)
- Authors: J. Garriga, A. Vilenkin. Phys. Rev. D 44, 1007.
- URL: <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.44.1007>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `E1` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: a covariant perturbation formalism for walls and strings, including nucleated bubbles and walls in de Sitter space; the perturbation field often has a tachyonic mass and a nonminimal coupling to worldsheet curvature, and, without gravity, one scalar describes wall translations.
- Relevance to HDBLAST (P2): The standard way to separate wall-translation (gauge) modes from physical instabilities. This round's gauge-invariant stability workstream already classifies the mu^2 = -4 harmonics on the +1 branch as pure gauge in the bulk.
- Implication: `method`

**C2. Brane-world creation and black holes** (2000)
- Authors: J. Garriga, M. Sasaki. Phys. Rev. D 62, 043523; arXiv:hep-th/9912118.
- URL: <https://arxiv.org/abs/hep-th/9912118>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `B1` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: an inflating brane-world can be created from nothing with its AdS bulk; the de Sitter brane instanton has AdS inside; 5D and 4D actions and entropies agree when the instanton is much larger than the AdS radius; thermal instantons with a bulk AdS black hole may describe creation of a hot universe.
- Relevance to HDBLAST (P1, P2): The Euclidean static branch (an S^4 shell bounding a hyperbolic 5-ball, Z2-doubled) is this instanton with a bulk scalar. If the branch is to be read as a tunnelling or creation saddle, its negative-mode count must be computed. The Lorentzian stability result of this round does not settle that count.
- Implication: `method`

**C3. Can Inflating Braneworlds be Stabilized?** (2003)
- Authors: A. Frolov, L. Kofman. arXiv:hep-th/0309002.
- URL: <https://arxiv.org/abs/hep-th/0309002v1>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `E3` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: for de Sitter branes the radion mass squared is typically negative, confirmed with BraneCode numerics. The programme's Chat 9 file (full text read via ar5iv, 16 Sept 2026) adds: scalar perturbations in generalized longitudinal gauge reduce to a self-adjoint Sturm-Liouville problem with a Rayleigh bound m^2 <= -4H^2 + (2/3) int(dw/a) / int(dw/(a phi'^2)), and the programme reproduced its registered eigenvalue mu^2 = -7.7178716 with these equations.
- Relevance to HDBLAST (P2, P1): Because the problem is self-adjoint Sturm-Liouville, rigorous eigenvalue enclosures of the kind found in verification-method texts (E2) apply in principle. The +1-branch statement 'no normalizable scalar mode below 9/4' could be certified rather than only computed.
- Implication: `method`

**C4. Radion on the de Sitter brane** (2000)
- Authors: U. Gen, M. Sasaki. arXiv:gr-qc/0011078.
- URL: <https://arxiv.org/pdf/gr-qc/0011078>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `E2` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: with two de Sitter branes of opposite tension the radion has negative mass squared proportional to the brane curvature; with a single positive-tension de Sitter brane there is no radion.
- Relevance to HDBLAST (P2): Fixes which scalar modes can exist with one shell and a regular cone.
- Implication: `context`

**C5. Thick Brane Worlds and Their Stability** (2001)
- Authors: S. Kobayashi, K. Koyama, J. Soda. Phys. Rev. D 65, 064014 (2002); arXiv:hep-th/0107025.
- URL: <https://arxiv.org/abs/hep-th/0107025>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `E5` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: non-singular thick Poincare, de Sitter and AdS branes with a dilaton have positive-definite effective potentials in their scalar master equations, hence stability.
- Relevance to HDBLAST (P2): A master-variable (Schrodinger-form) method. A positive-definite effective potential would be an analytic stability proof, stronger than a numerical spectrum.
- Implication: `method`

**C6. Fake supergravity and domain wall stability** (2004)
- Authors: D. Z. Freedman, C. Nunez, M. Schnabl, K. Skenderis. Phys. Rev. D 69, 104027.
- URL: <https://ui.adsabs.harvard.edu/abs/2004PhRvD..69j4027F/abstract>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `E6` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: a generalized Witten-Nester spinor argument proves stability of flat domain walls without supersymmetry in any dimension, extended to AdS-sliced walls.
- Relevance to HDBLAST (P2): Covers the balanced flat wall (delta = 0) only. The record shows no analogous positive-energy theorem for de Sitter-sliced walls, so +1-branch stability has to be computed, or proved by a new argument.
- Implication: `constrains`

**C7. Hidden Supersymmetry of Domain Walls and Cosmologies** (2006)
- Authors: K. Skenderis, P. K. Townsend. Phys. Rev. Lett. 96, 191301.
- URL: <https://www.osti.gov/etdeweb/biblio/20777232>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `E7` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: every domain-wall solution with Minkowski or AdS worldvolume admits first-order equations with a solution-determined superpotential; FLRW cosmologies obey similar equations through pseudo-Killing spinors.
- Relevance to HDBLAST (P1): A first-order reformulation reduces the static BVP to a lower-order system, which is often the key simplification in a computer-assisted proof. The de Sitter-sliced case needs the curved version.
- Implication: `method`

**C8. Coupled bulk and brane fields about a de Sitter brane** (2006)
- Authors: A. Cardoso, K. Koyama, A. Mennim, S. S. Seahra, D. Wands. Phys. Rev. D 75, 084002 (2007); arXiv:hep-th/0612202.
- URL: <https://arxiv.org/abs/hep-th/0612202>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `E8` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: a bulk scalar coupled linearly to a scalar on a de Sitter brane has zero, one or two bound states plus resonances; with brane test-field matter an infinite ladder of discrete tachyonic modes can become normalizable.
- Relevance to HDBLAST (P2, P6): Any matter extension that couples a shell field to phi needs a coupled-sector negative-mode check. Adding brane fields can create normalizable tachyons that the vacuum analysis cannot see.
- Implication: `constrains`

**C9. Complex scalar field thick branes: stability of linear perturbation and evolution of scalar Kaluza-Klein modes coupled with gravity** (2026)
- Authors: not captured in snippet. arXiv:2609.21421 (18 Sept 2026).
- URL: <https://arxiv.org/abs/2609.21421>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `E11` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: Minkowski, de Sitter and AdS thick branes from a complex scalar show no instability in the scalar and vector sectors; a factorized tensor equation excludes tachyonic tensor modes.
- Relevance to HDBLAST (P2): A current sector-by-sector template, including the vector sector, which this round's stability workstreams did not treat.
- Implication: `method`

**C10. Revisiting Coleman-de Luccia transitions in the AdS regime using holography** (2021)
- Authors: J. K. Ghosh, E. Kiritsis, F. Nitti, L. T. Witkowski. JHEP (2021); arXiv:2102.11881.
- URL: <https://arxiv.org/abs/2102.11881>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `B10` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: CDL AdS-to-AdS decays in Einstein-scalar theories are interpreted as vev-driven holographic RG flows on de Sitter space; such flows do not exist for generic potentials, and existence is tied to exotic RG flows.
- Relevance to HDBLAST (P1, P2): Existence of de Sitter-sliced regular solutions is non-generic, so an existence proof for the +1 branch is not a formality. P2: gives the tunnelling interpretation in which a negative-mode count matters.
- Implication: `method`

**C11. Searching for Coleman-de Luccia bubbles in AdS compactifications** (2022)
- Authors: G. Dibitetto, N. Petri. Phys. Rev. D 107, 046020 (2023); arXiv:2207.02172.
- URL: <https://arxiv.org/abs/2207.02172>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `B11` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: fully backreacted instantons are obtained by numerically integrating first-order Hamilton-Jacobi equations for smooth de Sitter-foliated walls between AdS vacua.
- Relevance to HDBLAST (P1): An independent solver architecture (first-order Hamilton-Jacobi) for the same class of solutions; a cross-check on the +1 branch that shares no code with the checkpoint's shooting solvers.
- Implication: `method`

**C12. Holographic RG flows on curved manifolds and quantum phase transitions** (2017)
- Authors: J. K. Ghosh, E. Kiritsis, F. Nitti, L. T. Witkowski. JHEP 05 (2018) 034; arXiv:1711.08462.
- URL: <https://arxiv.org/abs/1711.08462>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `F1` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: Einstein-dilaton flows with dS, AdS or sphere slicing and general potential; UV and IR asymptotics; new flows at finite curvature with no flat counterpart.
- Relevance to HDBLAST (P1): The regular-endpoint (IR) expansions classified there are the local cone solutions a certification needs. In this problem the linear germ is the explicit polynomial C_14^(2)(cosh u) (pilot, part 1).
- Implication: `method`

**C13. Generalized surface tension bounds in vacuum decay** (2017)
- Authors: A. Masoumi, S. Paban, E. J. Weinberg. arXiv:1711.06776.
- URL: <https://arxiv.org/abs/1711.06776>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `B16` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: defines a tension for thin walls whose curvature radius changes across the wall and derives CDL-like bounds for bounces with Minkowski or AdS false vacua.
- Relevance to HDBLAST (P2): Context for tunnelling interpretations of a shell with a nonconstant scalar profile.
- Implication: `context`

**C14. Gravitational Effects on and of Vacuum Decay** (1980 (per sweep-1 caveat, year from general knowledge))
- Authors: S. Coleman, F. De Luccia. Phys. Rev. D (1980).
- URL: <https://www.semanticscholar.org/paper/Gravitational-Effects-on-and-of-Vacuum-Decay-Coleman-Luccia/a60b863946ca2d85da9d8b48b72257e8ec05e3ed>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `B15` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: gravity can prevent decay of a Minkowski or AdS false vacuum; in the thin-wall limit this happens above a tension bound.
- Relevance to HDBLAST (P2): Origin of the CDL setting in which the 'negative-mode problem' arises. The negative-mode literature itself (Lavrelashvili, Lee-Weinberg, Koehn-Lavrelashvili-Lehners and others) could not be retrieved in this sweep; see LEADS.
- Implication: `context`

### D. Non-perturbative particle production and backreaction

**D1. Instant preheating (title and authors not written in the programme record; see caveat)** (1998 (arXiv number))
- Authors: not stated in the programme record. arXiv:hep-ph/9812289v2.
- URL: <https://arxiv.org/abs/hep-ph/9812289v2>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 142/00_READ_FIRST.md` (cited as the established benchmark in the 6 Sept 2026 read-first).
- Finding: Per the programme record: the established instant-preheating benchmark, adapted there to a brane field with m_eff^2 = m_chi^2 + g_b^2 (phi - phi_*)^2 crossing linearly: n_k = exp[-pi (k^2 + m_chi^2)/mu^2], n_chi = mu^3/(8 pi^3) exp[-pi m_chi^2/mu^2], mu^2 = g_b |phi_dot_b|, in the locally flat adiabatic in/out limit. Backreaction, the quantum state and thermalization were left open.
- Relevance to HDBLAST (P6): The mechanism used by this round's preheating workstream, which found it insufficient on the fixed archived trajectory. Its open item is the self-consistent bulk-coupled backreaction.
- Implication: `method`
- *Caveat:* To this sweep's general knowledge the identifier is Felder, Kofman and Linde, 'Instant preheating'. This attribution was not verified here.

**D2. Exact Bogoliubov coefficients for a sech-profile mass pulse (paper title not recorded)** (2019)
- Authors: Alves and Camilo (as recorded). Eur. Phys. J. C (2019), article s10052-019-6581-2.
- URL: <https://link.springer.com/article/10.1140/epjc/s10052-019-6581-2>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 144/HDBLAST_CHAT8_INFORMATION_AND_QUANTUM_20260906/quantum/EXACT_SCALAR_PULSE_PRODUCTION.md` (cited with equation numbers in the 6 Sept 2026 exact-pulse memo).
- Finding: Per the programme record: equations (28)-(32) give the massless profile and Bogoliubov coefficients of the solvable sech mass quench, which the programme reproduced as a benchmark; integer indices give zero asymptotic production.
- Relevance to HDBLAST (P6): An exact calibration for any mode-function particle-production code, alongside the Landau-Zener limit used for instant preheating.
- Implication: `method`

**D3. Poschl-Teller non-production points (paper title not recorded)** (2022 (arXiv number))
- Authors: Ahmadiniaz et al. (as recorded). arXiv:2205.15946.
- URL: <https://arxiv.org/pdf/2205.15946>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 144/HDBLAST_CHAT8_INFORMATION_AND_QUANTUM_20260906/quantum/EXACT_SCALAR_PULSE_PRODUCTION.md` (cited (section V) in the 6 Sept 2026 exact-pulse memo).
- Finding: Per the programme record: section V discusses Poschl-Teller non-production points in an electric-field construction with additional momentum dependence.
- Relevance to HDBLAST (P6): A second exact benchmark family for particle-production codes (reflectionless profiles).
- Implication: `method`

**D4. Parker (1969), Physical Review 183, 1057 (title not written in the programme record)** (1969)
- Authors: L. Parker (surname as recorded). Phys. Rev. 183, 1057.
- URL: <https://journals.aps.org/pr/abstract/10.1103/PhysRev.183.1057>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 144/HDBLAST_CHAT8_INFORMATION_AND_QUANTUM_20260906/quantum/EXACT_SCALAR_PULSE_PRODUCTION.md` (cited for the conformal non-production check).
- Finding: Per the programme record: the established conformal non-production result (beta_k = 0 for a conformally coupled massless field in the conformal vacuum), used as a control.
- Relevance to HDBLAST (P6): A mandatory zero-production control for any gravitational particle-production calculation.
- Implication: `method`

**D5. Probing inflationary particle production with the CMB power spectrum** (2026)
- Authors: Jense, Abu El-Haj, Hill and Philcox. arXiv:2606.26823 (25 June 2026).
- URL: <https://arxiv.org/html/2606.26823v1>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 144/HDBLAST_CHAT8_INFORMATION_AND_QUANTUM_20260906/review/FRESH_SCOPE_AND_REFERENCES.md` (abstract, model setup and methods sections read 6 Sept 2026 by an earlier agent).
- Finding: Per the programme record: a particle-production burst from an inflaton-spectator mass crossing is turned into CMB templates normalized by the inverse-covariance norm (their eqs. 26-27, 30); the weak joint-data preference is not an established detection.
- Relevance to HDBLAST (P6): If a mass-crossing burst ever happens during an inflating stage of an HDBLAST variant, this gives the observable-template method. It does not address reheating.
- Implication: `method`

**D6. Massive Graviton Dark Matter from a Gapped Continuum** (2026)
- Authors: Megias, Prieto and Quiros. arXiv:2607.07295 (v2 3 Aug 2026).
- URL: <https://arxiv.org/html/2607.07295v2>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 144/HDBLAST_CHAT8_INFORMATION_AND_QUANTUM_20260906/review/FRESH_SCOPE_AND_REFERENCES.md` (abstract and sections 2.1, 6, 7 inspected 6 Sept 2026 by an earlier agent).
- Finding: Per the programme record: in a linear-dilaton braneworld, brane-to-bulk graviton leakage scales as the eighth power of temperature, and the hidden fluid later scales as matter.
- Relevance to HDBLAST (P6): Consistent backreaction must include energy leaving the brane into the bulk. Not every bulk channel acts as dark radiation.
- Implication: `constrains`

**D7. RS2 gravitational-wave calculation with bulk Kaluza-Klein leakage (title not recorded)** (2006)
- Authors: Hiramatsu (as recorded). arXiv:hep-th/0601105v2.
- URL: <https://arxiv.org/abs/hep-th/0601105v2>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 142/00_READ_FIRST.md` (cited in the 6 Sept 2026 read-first).
- Finding: Per the programme record: in established RS2 calculations, bulk Kaluza-Klein leakage can cancel an expansion enhancement in their radiation examples.
- Relevance to HDBLAST (P6): A numerical precedent for brane-bulk energy bookkeeping in a radiation era; it limits any expansion-only inference.
- Implication: `constrains`

**D8. Particle Creation from Entanglement Entropy** (2026)
- Authors: not stated in the programme record. PTEP 2026, 013A01.
- URL: <https://academic.oup.com/ptep/article/2026/1/013A01/8400342>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 144/HDBLAST_CHAT8_INFORMATION_AND_QUANTUM_20260906/review/FRESH_SCOPE_AND_REFERENCES.md` (abstract and methods read by an earlier agent).
- Finding: Per the programme record: the main calculation is restricted to low-velocity moving mirrors; the programme treated it as an analogy, not a brane energy source.
- Relevance to HDBLAST (P6): Moving-mirror (dynamical Casimir) production is the closest analogue of production by a moving shell. The low-velocity restriction does not cover the Chat 14 roll-off.
- Implication: `context`

**D9. Collision of Domain Walls and Reheating of the Brane Universe** (2004)
- Authors: Y. Takamizu, K. Maeda. Phys. Rev. D 70, 123514; arXiv:hep-th/0406235.
- URL: <https://arxiv.org/abs/hep-th/0406235>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `C12` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: particle production when two domain walls collide in 5D, with a reheating estimate (formula garbled in the snippet); a follow-up treats asymptotically AdS5.
- Relevance to HDBLAST (P6): A wall-collision production estimate against which a moving-shell production code can be compared.
- Implication: `method`

**D10. A dynamical inflaton coupled to strongly interacting matter** (2023)
- Authors: C. Ecker, E. Kiritsis, W. van der Schee. Phys. Rev. Lett. 130, 251001; arXiv:2302.06618.
- URL: <https://arxiv.org/abs/2302.06618>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `F6` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: the Einstein-inflaton equations are coupled self-consistently to a holographic strongly coupled QFT, giving inflation, reheating and a thermal QFT-dominated universe.
- Relevance to HDBLAST (P6): The only item in this session's records that computes reheating with full self-consistent backreaction; the heat bath is a 5D (holographic) sector. Methodologically the closest template for bulk-coupled feedback, which this round's preheating workstream did not solve.
- Implication: `method`

**D11. Gravitational reheating at strong coupling** (2023)
- Authors: A. Buchel. JHEP 07 (2023) 159; arXiv:2304.11195.
- URL: <https://arxiv.org/abs/2304.11195>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `F7` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: estimates the maximal reheating temperature of strongly coupled theories after a rapid exit from de Sitter; most efficient when H greatly exceeds the conformal-breaking scale.
- Relevance to HDBLAST (P6): An upper-bound method for gravitational production during a change of de Sitter rate.
- Implication: `constrains`

**D12. Reheating the Universe in Braneworld Cosmological Models with bulk-brane energy transfer** (2008)
- Authors: not captured in snippet. arXiv:0805.1792.
- URL: <https://arxiv.org/abs/0805.1792>
- Access: search snippet only (second-hand, through the sweep-1 record). Origin of the URL: sweep-1 record `D11` (`theory_checks/theory_sources.json`).
- Finding: Per the sweep-1 snippet record: analyses the reheating era in 5D braneworlds with brane-bulk energy exchange.
- Relevance to HDBLAST (P6): A precedent for the energy-exchange law of the programme's matter extension.
- Implication: `context`

### E. Validated numerics and computer-assisted proofs

**E1. Computer-assisted construction of SU(2)-invariant negative Einstein metrics** (2026)
- Authors: Qiu Shi Wang. Ann. Global Anal. Geom. 69, 18 (2026); arXiv:2504.21644v2.
- URL: <https://arxiv.org/abs/2504.21644v2>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 136/M489G_LITERATURE_AND_PUBLICATION_REVIEW_20260903.md` (publisher full text and versioned preprint, read 3 Sept 2026 by an earlier agent).
- Finding: Per the programme record: high-accuracy approximate solutions, rigorous residual bounds and a fixed-point argument yield actual Riemannian Einstein metrics, with explicit constants and reproducible computations.
- Relevance to HDBLAST (P1): The closest structural precedent found. The Euclidean continuation of the +1 branch is a symmetry-reduced (cohomogeneity-one) Riemannian Einstein-scalar metric on a ball with a regular centre and an S^4 shell boundary. The Riemannian approximation -> residual -> fixed-point architecture therefore applies more directly to the static branch than to any Lorentzian evolution. The static Lorentzian problem is the same ODE BVP.
- Implication: `method`

**E2. Numerical Verification Methods and Computer-Assisted Proofs for Partial Differential Equations** (2019)
- Authors: M. T. Nakao, M. Plum, Y. Watanabe. Springer book, DOI 10.1007/978-981-13-7669-6.
- URL: <https://doi.org/10.1007/978-981-13-7669-6>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 135/HDBLAST_PHYS_M488G_FIXED_P0_CAUSAL_EXCISION_AND_GOURSAT_ADJOINT_RESEARCH_CHECKPOINT_20260902/lanes/lane_methods/PHYS_M488G_CHARACTERISTIC_PDE_AND_VALIDATED_NUMERICS_METHODS_ROUTE.md` (listed as a reference for the computer-assisted-proof gate, 2 Sept 2026).
- Finding: Per the programme record: title, authors, year and DOI only; the record lists it as a primary reference for its next computer-assisted-proof step and does not summarise the content.
- Relevance to HDBLAST (P1, P2): The standard reference class for verified solutions of elliptic/ODE BVPs and verified eigenvalue enclosures. The latter would certify spectral statements such as 'no scalar mode below 9/4 on the +1 branch'. Content not checked in this sweep.
- Implication: `method`

**E3. Verification Methods: Rigorous Results Using Floating-Point Arithmetic** (2010)
- Authors: S. M. Rump. Acta Numerica 19 (2010), 287-449.
- URL: <https://doi.org/10.1017/S096249291000005X>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 135/HDBLAST_PHYS_M488G_FIXED_P0_CAUSAL_EXCISION_AND_GOURSAT_ADJOINT_RESEARCH_CHECKPOINT_20260902/lanes/lane_methods/PHYS_M488G_CHARACTERISTIC_PDE_AND_VALIDATED_NUMERICS_METHODS_ROUTE.md` (listed as a reference for the computer-assisted-proof gate, 2 Sept 2026).
- Finding: Per the programme record: bibliographic entry only (title, venue, pages, DOI).
- Relevance to HDBLAST (P1): The survey behind interval-Newton/Krawczyk-type existence tests with outward-rounded floating point, the step used by the programme's earlier certified BVP root (M462). Content not checked in this sweep.
- Implication: `method`

**E4. Verified Bounds for the Determinant of Real or Complex Point or Interval Matrices** (2020)
- Authors: S. M. Rump. J. Comput. Appl. Math. 372 (2020), 112610.
- URL: <https://doi.org/10.1016/j.cam.2019.112610>
- Access: programme record only (not opened, searched or re-verified in this sweep). Origin of the URL: programme file `D-Blast 3/untitled folder 135/HDBLAST_PHYS_M488G_FIXED_P0_CAUSAL_EXCISION_AND_GOURSAT_ADJOINT_RESEARCH_CHECKPOINT_20260902/lanes/lane_methods/PHYS_M488G_CHARACTERISTIC_PDE_AND_VALIDATED_NUMERICS_METHODS_ROUTE.md` (listed as a reference for the computer-assisted-proof gate, 2 Sept 2026).
- Finding: Per the programme record: bibliographic entry only (title, venue, DOI).
- Relevance to HDBLAST (P1, P2): Certifying that a junction Jacobian is nonsingular is the local-uniqueness half of a Krawczyk argument. This round's stability workstream reports floating-point condition numbers of 15-17 for the static junction Jacobian. Content not checked in this sweep.
- Implication: `method`


## Leads to verify (general knowledge; no URL; not verified; not cited)

These are the standard references for the methods the task names. They are listed so that a later run with search enabled can verify them. The bibliographic details may be wrong.

*A. Constraint-damped formulations, outer boundaries and error control*

- Gundlach, Calabrese, Hinder, Martin-Garcia (2005): constraint damping in the Z4 formulation and harmonic gauge. Why: Damping terms that drive constraint violations to zero at a tunable rate. The direct remedy for the outgoing constraint front (P3). (planned query 1)
- Alic, Bona-Casas, Bona, Rezzolla, Palenzuela (2012): conformal and covariant Z4 (CCZ4). Why: Standard damped formulation used in spherical and higher-dimensional codes (P3, P4). (planned query 2)
- Bernuzzi, Hilditch (2010): Z4c; and later work on constraint-preserving boundary conditions for Z4c. Why: Constraint-preserving outer boundaries, the other half of the far-boundary problem (P3). (planned query 3)
- Pretorius (2005) and Lindblom, Scheel, Kidder, Owen, Rinne (2006): generalized harmonic evolution with constraint damping. Why: Alternative damped formulation; widely used for AdS and cosmological problems (P3, P4). (planned query 4)
- Zenginoglu (2008): hyperboloidal foliations and scri-fixing; later hyperboloidal spherical-symmetry work. Why: P3: in the evolution chart the far end of the bulk (z -> -infinity) is the Lorentzian continuation of the regular cone. That is a Rindler-type horizon of the de Sitter slicing, reached only asymptotically in coordinate time. Horizon-penetrating or hyperboloidal-type slices would let outgoing signals leave instead of reflecting from a finite boundary. (planned query 5)
- Sarbach, Tiglio (2012), Living Reviews in Relativity: continuum and discrete initial-boundary value problems and Einstein's equations. Why: Well-posed boundary conditions and summation-by-parts (SBP) discretizations with energy estimates (P3, P5). (planned query 8)
- Summation-by-parts operators with simultaneous-approximation-term (SAT) boundary and interface closures (review literature on SBP-SAT). Why: P5: Chat 14 found that the order of the one-sided shell closure dominated the error at delta=0.001. SBP-SAT closures give provably stable boundary/interface treatments of a chosen order. (planned query 8)
- Higher-dimensional numerical relativity with symmetry reduction (modified cartoon method, Shibata-Yoshino; 5D CCZ4 codes). Why: Established 5D formulations if the programme leaves the 1+1 reduction (P3). (planned query 4)
- Chesler, Yaffe (2014): characteristic formulation for asymptotically AdS gravitational dynamics. Why: Characteristic (null) evolution in AdS5 avoids an artificial timelike outer boundary (P3, P4). (planned query 4)

*B. Thin shells, junction conditions, moving boundaries and collapse diagnostics*

- Israel (1966): singular hypersurfaces and thin shells; Ida (2000): brane-world cosmology from junction conditions. Why: Foundations for the moving-shell benchmarks (P5). (planned query 7)

*C. Gauge-invariant perturbations and negative modes of de Sitter-sliced walls*

- Lavrelashvili (2006) and Lee, Weinberg (2014): negative modes of Coleman-De Luccia bounces. Why: The CDL negative-mode problem: counting and gauge-invariant definition of the negative mode (P2). (planned query 10)
- Koehn, Lavrelashvili, Lehners (2015): towards a solution of the negative-mode problem in quantum tunnelling with gravity. Why: Gauge-invariant Euclidean perturbation variables for O(4)/O(5)-symmetric bounces (P2). (planned query 10)
- Tanaka, Sasaki (1992); Khvedelidze, Lavrelashvili, Tanaka (2000); Dunne, Wang (2006); Gratton, Turok (2000): negative modes and fluctuations of gravitational instantons. Why: Methods for counting negative modes of the Euclidean +1 and -1 branches (P2). (planned query 11)
- 'Gregory-Ruth' (named in the task): NOT identified. No record in this session or in the scanned programme files matches this pair of names. Why: Left open; not cited. (planned query 12)

*D. Non-perturbative particle production and backreaction*

- Kofman, Linde, Starobinsky (1997): towards the theory of reheating after inflation. Why: Broad-resonance preheating; background for P6. (planned query 14)
- Figueroa, Florio, Torrenti, Valkenburg (2021): CosmoLattice; and lattice reviews of early-universe field dynamics (2020-2023). Why: Lattice tools with full classical backreaction for P6; no 2020-2026 lattice paper with renormalized backreaction could be verified in this sweep. (planned query 16)
- Kofman, Linde, Liu, Maloney, McAllister, Silverstein (2004): 'Beauty is attractive', moduli trapping at enhanced symmetry points. Why: Backreaction of production at a mass crossing on the rolling field, the self-consistency issue flagged by the preheating workstream (P6). (planned query 14)
- Parker, Fulling (1974): adiabatic regularization; later semiclassical-backreaction work on renormalized stress tensors. Why: Renormalized energy density needed for any backreacted production calculation (P6). (planned query 15)
- Durrer, Ruser (2007): dynamical Casimir effect in braneworlds (graviton production by a moving brane). Why: Particle production by a moving brane in AdS5, the closest physical analogue for P6. (planned query 15)

*E. Validated numerics and computer-assisted proofs*

- Krawczyk (1969); Moore (1966); Tucker (2011, 'Validated Numerics'): interval Newton/Krawczyk existence tests. Why: The existence/uniqueness step used by the programme's M462 certified BVP root (P1). (planned query 18)
- Kashiwagi's kv library (validated Taylor ODE integration; used by the programme's M462 release as kv v0.4.62); CAPD::DynSys (Kapela et al. 2021); Arb/FLINT (Johansson). Why: Validated ODE integrators for a shooting-type proof of the +1 branch (P1). None is installed in this environment. (planned query 17)
- Radii-polynomial approach: Day, Lessard, Mischaikow (2007); van den Berg, Lessard (2015). Why: Newton-Kantorovich proofs for BVPs in Chebyshev/Taylor bases; an alternative to shooting (P1). (planned query 17)
- Lohner (1987): the wrapping effect in interval ODE integration. Why: Controls overestimation in validated shooting; the pilot suggests the transverse direction contracts strongly outward (P1). (planned query 18)


## Planned queries (prepared, not run)

1. `constraint damping Z4 formulation spherical symmetry Einstein scalar field evolution`
2. `CCZ4 spherical symmetry scalar field collapse constraint damping parameters`
3. `Z4c spherical symmetry constraint-preserving boundary conditions scalar field`
4. `generalized harmonic spherically symmetric scalar field anti-de Sitter numerical evolution`
5. `hyperboloidal compactification spherical symmetry scalar field evolution 2023 2024 2025`
6. `horizon-penetrating hyperboloidal slice Rindler horizon de Sitter slicing AdS numerical`
7. `numerical relativity thin shell Israel junction conditions moving boundary evolution`
8. `summation-by-parts simultaneous approximation term interface conditions numerical relativity`
9. `braneworld numerical evolution bulk scalar brane junction conditions 2020 2026`
10. `negative mode Coleman-De Luccia bounce gravity Lavrelashvili Lee Weinberg`
11. `negative modes de Sitter brane instanton Garriga Sasaki bulk scalar Euclidean`
12. `Gregory Ruth negative mode bubble brane (identify the intended reference)`
13. `gauge-invariant perturbations de Sitter sliced domain wall bulk scalar master variable 2020 2026`
14. `instant preheating lattice simulation backreaction 2021 2022 2023 2024 2025`
15. `non-perturbative particle production renormalized stress tensor backreaction adiabatic regularization 2022 2026`
16. `CosmoLattice preheating tachyonic lattice 2024 2025`
17. `computer-assisted proof boundary value problem ODE radii polynomial Einstein metric`
18. `validated numerics interval arithmetic existence proof domain wall soliton gravity`
19. `computer-assisted eigenvalue enclosure Sturm-Liouville stability rigorous`
20. `rigorous numerics self-similar solution Einstein scalar computer-assisted proof`

## Limitations

- **No new search.** The 2020–2026 methods literature was not surveyed in this sweep. Absence of an item here says nothing about whether it exists.
- **Second-hand evidence.** Sweep-1 entries rest on snippet paraphrases made by another agent. Programme entries rest on earlier project records that were not re-verified.
- **Recommendations are judgements.** The methods map is untested. No damped formulation, hyperboloidal chart, SBP-SAT closure, apparent-horizon finder or eigenvalue enclosure was implemented.
- **The pilot is not a proof.** It encloses a truncated closed-form model and compares it with floating-point nonlinear solutions.

## Reproduction

From `literature/`, with python3, SymPy 1.14 and mpmath 1.3:

```
python3 methods_checks/certification_pilot.py        # ~1 s, 1 core; writes methods_checks/certification_pilot.json
python3 methods_checks/build_methods_md.py           # validates records, writes methods.md and methods_checks/methods_sources.json
python3 methods_checks/builder_negative_control.py   # corrupted records must be refused; writes builder_negative_control.json
```
