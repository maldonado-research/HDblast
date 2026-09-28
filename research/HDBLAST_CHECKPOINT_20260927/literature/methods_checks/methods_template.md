# HDBLAST literature sweep 3: mathematical and numerical methods

Research round of 27 September 2026, written 28 September 2026. **Work in progress.** This is an annotated bibliography with a methods map and one small validated-numerics pilot. Nothing here is a claim of the HDBLAST programme until it appears in a dated checkpoint or a Zenodo record.

## Read this first: the web search did not run

**This sweep could not be carried out as specified (procedural negative result).**

- The task required at least 15 distinct web searches. **None was executed.** All four attempted queries were refused by the tool: the session's web-search budget (200 of 200 calls) had already been used by earlier work in the same session. The refusal notice asked to continue with information already gathered, so no other search route was tried.
- A page fetch was also blocked: `journals.aps.org` returned `EGRESS_BLOCKED`.
- Attempts, verbatim:

{{ATTEMPTED}}

Because of this, the bibliography below is assembled only from sources that are already on record. It does **not** reflect a fresh search of the 2020–2026 literature. Two kinds of source were used, and every entry says which one it is:

1. **Sweep-1 search records.** The URL was returned by a web search earlier in this session, during sweep 1 (theory), and is stored in `theory_checks/theory_sources.json`. The finding here re-reads that record's snippet paraphrase from a methods angle. The evidence is therefore **second-hand search-snippet evidence**.
2. **Programme records.** The URL is written verbatim in a file of the programme's own archive (`new-files/D-Blast 3/…`). An earlier project agent reports having opened it on the stated date. The finding paraphrases that record only. **This sweep did not open, search or re-verify these sources.**

The builder (`methods_checks/build_methods_md.py`) machine-checks each URL against its claimed origin. It refuses to build if a URL is not byte-identical to its sweep-1 record, or not present verbatim in the named programme file.

A third list, **leads**, gives standard methods references from general knowledge. The task names several of them: CCZ4, Z4c, generalized harmonic, hyperboloidal compactification, the CDL negative-mode papers, CosmoLattice, and validated ODE libraries. They carry **no URL**, were **not verified**, may contain bibliographic errors, and **are not cited as sources**. They are paired with 20 prepared queries for a later run.

**Totals:** {{COUNTS}}.

No paper was read in full in this sweep. Every access label is either "search snippet only" (second-hand) or "programme record only".

## Question

Which mathematical and numerical methods would most directly advance the programme's open problems? The problem list below is taken from the 22 September checkpoint (Zenodo 22922928), Chat 14, the registered-shell controls, and this round's workstreams.

| Code | Open problem | Sources below |
|---|---|---|
{{PROBLEM_MAP}}

## Main findings

Each finding carries a label. **exact** means verified symbolically in this folder. **numerical** means floating-point or multiprecision with stated evidence. **conditional** means it depends on a stated assumption. **judgement** means a methods recommendation, not a result.

1. **Certifying the +1 branch is feasible with methods the programme already uses (conditional, judgement).** The programme already holds a computer-assisted existence-and-local-uniqueness proof for a static shell. That proof is PHYS-M462, 17 Aug 2026, in `D-Blast 3/untitled folder 139/.../M462_…CERTIFIED_FIXED_T_CONTINUUM_BVP_ROOT.md`. It covers the original shell at t = 0.001, not the +1 branch. Its ingredients were:
   - a certified regular endpoint germ;
   - a pinned `kv` v0.4.62 validated Taylor integrator (order 12, outward-rounded binary64);
   - a defect chart that subtracts the flat-wall direction;
   - a strict Krawczyk inclusion in a two-parameter shooting box.

   Porting this pipeline to the +1 branch is the most direct route to P1. The pilot below identifies what the port needs:
   - **A local germ.** It is explicit: the regular linear cone solution is exactly η = a·C₁₄⁽²⁾(cosh u)/680, with u = y/9 (**exact**).
   - **An outward direction that is well conditioned in the linear regime (conditional).** The unwanted, cone-singular solution decays like e^{−18u} while the regular one grows like e^{14u} (**exact**). Over [0.25, u_b] the unwanted direction is suppressed by about 10^{{SUPP}}, with u_b = {{UB}} at δ = 0.001.
   - **The η-formulation, which is mandatory.** The germ amplification is G(u_b) = {{AMP}} (10^{{LOGAMP}}). In binary64, 1+η_h rounds to exactly 1. A direct φ-formulation would need about {{DIGITS}} decimal digits to carry η_h to 16 significant digits. The checkpoint already uses η = φ − 1, which removes this problem.
   - **A shooting scale.** By analogy with M462's scaled shooting variable, a natural variable of order δ⁰ is â = η_h·G(u_b)/δ. At δ = 0.001 it is {{AHAT}}, against the leading-order value −9c/64 = {{M9C64}}.
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

Script: `methods_checks/certification_pilot.py` → `methods_checks/certification_pilot.json`. It uses SymPy 1.14 and mpmath 1.3 (`mpmath.iv` at 60 digits, point checks at 100 digits), runs in {{RUNTIME}} s on one core, and has script SHA-256 prefix `{{PILOT_SHA}}`. The input is the checkpoint's `static_branch/PLUS_BRANCH_RESULTS.json`, read-only with its SHA-256 checked.

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
| 100-digit point values inside all {{NCONT}} interval enclosures | all contained | interval |
| LO metric-only H² vs checkpoint `metric_only_H2` | agree to ≤ 6×10⁻¹⁴ relative | numerical |
| LO η_b vs checkpoint `eta_finite_curvature_linear` | agree to ≤ 4×10⁻¹⁶ relative (same linear model, implemented independently) | numerical |
| Controls: Gegenbauer degree 13 or 15, index 3/2, 4D-like damping coefficient | all fail, as required | control |
| Control: degree-13 amplification vs the checkpoint's η_b/η_h | misses by > 90 % | control |
| Control: c perturbed by +1 % | detected at the three smallest δ (deviation > 10× baseline) | control |
| Control: wrong-sign scalar junction | detected: η_b = {{WRONGSIGN}}, about 8× too large (see note) | control |
| Control: containment test on a deliberately shifted point | rejected, as required | control |

Comparison of the LO enclosures with the checkpoint's nonlinear floating-point solutions (numerical):

{{PILOT_TABLE}}

- The transport ratio (η_b/η_h)/G(u_b) − 1 measures only the nonlinear correction to the germ between cone and shell, since it uses the checkpoint's own u_b. It scales as {{COEF_RATIO}}·δ; the log–log slope over the three smallest δ is {{SLOPE_RATIO}}.
- The LO error in η_b is {{COEF_ETAB}}·δ (slope {{SLOPE_ETAB}}).
- The LO error in H² is {{COEF_H2}}·δ (slope {{SLOPE_H2}}). This is the known −c²δ²/384 scalar-profile term of the checkpoint's expansion: the predicted coefficient is 27c²/(384(1+c)) = {{H2EXP}}, and the observed value at δ = 3×10⁻⁴ is {{H2OBS}}.
- The LO error in η_h has a small O(δ) coefficient that changes sign near δ ≈ 0.01. The coefficients by increasing δ are {{COEF_ETAH_LIST}}, and the three-point slope is {{SLOPE_ETAH}}. The O(δ) and O(δ²) terms partly cancel there, so this slope is less clean than the others.
- The differences between the checkpoint's default and refined solver settings are ≤ 2×10⁻¹² relative for η_h. That is far below every deviation in the table, so the table measures the model truncation, not solver noise.

**Note on a control.** The wrong-sign-junction control was first written to expect a sign flip of η_b. That expectation was wrong: with k·G_u/G ≈ 14/9 < 2 the wrong-sign denominator is also negative, so η_b keeps its sign and becomes about 8 times too large. The first run recorded "not met". The criterion was changed to a magnitude test, and both outcomes are stored in `control_design_note` in the JSON.

**Scope.** Items 1–4 of the pre-declared expectations in the script docstring were met. The pilot shows that the linear germ and conditioning are favourable. It proves nothing about the nonlinear branch: no existence, no uniqueness, no stability.

## Annotated bibliography

Every entry is based on a search snippet (second-hand, through the sweep-1 record) or on a programme record. **None was read in full in this sweep.** Authors are given only as they appear in the originating record.

{{BIBLIOGRAPHY}}

## Leads to verify (general knowledge; no URL; not verified; not cited)

These are the standard references for the methods the task names. They are listed so that a later run with search enabled can verify them. The bibliographic details may be wrong.

{{LEADS}}

## Planned queries (prepared, not run)

{{PLANNED}}

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
