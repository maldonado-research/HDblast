# HDBLAST: the late-time brane is stable, and the missing hot Big Bang is now a sharper problem

Research checkpoint, 27 September 2026 (calculations and independent audits completed 28 September 2026). Ricardo Maldonado's HDBLAST program; prepared with AI assistance. New calculations, an independent audit of each of the six computational workstreams (the literature sweeps and the program map are self-checked only), and snippet-level literature sweeps. No publication or external peer review is claimed.

**Hypothesis under test:** a higher-dimensional "blast", a five-dimensional gravitational event, sparked our Big Bang. The registered model is unchanged from the 22 September checkpoint (Zenodo record 22922928):

$$W=1-\phi+\tfrac13\phi^3,\qquad U=\tfrac12W_\phi^2-\tfrac23W^2,\qquad \sigma=2W+\delta(1+c\phi),$$

with δ = 0.001, c = 0.5975949350280132, a doubled (Z2) bulk and one shell. The junction conditions are ρ′/ρ = σ/6 and φ′ = −σ′/2.

## In plain words

1. **The strongest new result: the proposed late-time state is stable.** On 22 September we found a corrected static solution, the "+1 branch": an empty, inflating (de Sitter) brane with H/H₀ = 0.6067. This round tested whether small disturbances grow on it. Two independent methods say they do not: Method A covers every de Sitter harmonic of the scalar sector plus the tensor sector; Method B covers only the homogeneous (FRW-symmetric) scalar sector. Both methods first reproduced the known instability of the original shell, and both detected an instability we planted on purpose. Status: **numerical, independently verified**, but not a mathematical proof, and only for small (linear) disturbances.
2. **The most important negative result: nothing we tried turns the blast into a hot Big Bang.** We computed particle production ("instant preheating") along the actual simulated roll-off. We also screened six other routes to heat (dark radiation from the bulk, dark-bubble cosmology, a dissipative shell coupling, sudden conversion of tension, gravitational particle production, and an ekpyrotic-type collapse), plus two ways of retuning the model. None produces a radiation-dominated era in the registered model. The obstacle is a leftover vacuum energy on the final brane: radiation would have to beat it for about 22 e-folds, and the best budget allows about 0.6. Status: **conditional negative** (it holds under stated assumptions).
3. **Two leads remain open, and both change the model.** (a) Adding one quadratic term to the tension can cancel the leftover vacuum. Radiation then reaches 69% of the expansion budget in a pilot run at δ = 0.1 (100 times the registered detuning, within a budget with large cancellations), but "dark radiation" from the fifth dimension stays larger than observations allow, so far. (b) The independent auditor found that the tension parameter c can be continued to the vacuum-cancelling value after all (shown at δ = 0.1 only). That gives a new, different starting shell whose fate is not yet known; one exploratory single-grid run exists. Status: **inconclusive**.
4. **Three earlier claims in this round were wrong and are corrected here** (see §4). Nothing in this checkpoint is a discovery, a proof, a hot Big Bang mechanism, or observational support for the hypothesis.

What this means for the hypothesis: in its registered form, the 5D event ends as a stable, empty, inflating brane. That is a consistent universe, but not ours: ours had a hot, radiation-filled phase. To stay viable, HDBLAST now needs an extra ingredient or a model change that (i) removes or overwhelms the leftover vacuum and (ii) keeps dark radiation below the N_eff limits. §6 lists the tests that decide this.

## 1. What changed since 22 September

| Question left open on 22 Sept | Answer now | Status |
|---|---|---|
| Is the +1 branch stable? | Linearly stable (scalar with brane bending and both junctions, and tensor) at six δ values from 0.0003 to 0.1 (Method B: homogeneous sector at δ = 0.001, 0.003, 0.01); slowest decay is e^{−3H₊τ/2} | numerical, verified |
| Exact small-δ structure? | φ_b, H², ρ_b exact in c through δ⁸, and the leading (δ⁸) cone displacement; δ⁹ of φ_b and H² exact (auditor) | exact-verified |
| Why a degree-14 polynomial? | n = Δ − 4 with Δ = 3W″/W = 18 at φ = +1 | exact-verified |
| Can instant preheating give a radiation era? | No, in the scanned range at the M₅ cutoff (≤1.7×10⁻⁴ of the threshold for screen-passing couplings at the data end; ≤5.9×10⁻⁴ over the stored profile) on the fixed archived trajectory | conditional negative, verified |
| Constraint stall in the 5D evolutions? | Cause found (initial-data defect near the shell) and removed (order ≈3 restored) | numerical, verified |
| Any other route to heat? | Six screened: none works in the registered model; two model-changing leads open | mixed, see §3.4 |
| What does 2022–2026 literature say? | Strong prior art for "universe from a 5D event"; every hot-origin model adds an ingredient | snippet-level |

## 2. The headline result: stability of the +1 branch

### 2.1 Method A — gauge-invariant perturbation theory (`stability_gauge_invariant/`)

Perturbations are expanded in de Sitter harmonics, □Y = μ²Y, so the 4D mass is m² = μ²H². A mode is unstable if it is normalizable with μ² < 0; the continuum starts at μ² = 9/4. From the full linearized 5D Einstein–scalar equations, symbolic algebra (55/55 checks, 4/4 wrong-formula controls) gives a gauge-invariant first-order system in X and Z = Ψ/φ′. The shell condition is

$$M(\mu^2)=B\,X+\frac{3(\mu^2+4)}{\rho^2}Z=0\ \text{ at } y_b,\qquad B=\frac{\phi''}{\phi'}+\frac{\sigma''}{2}.$$

| Check | Result | Status | Audit |
|---|---|---|---|
| Calibration: original shell, δ = 0.001 | one bound state, μ² = −7.717871625260, growth 1.657193631259 (Chat 9: −7.717871625176) | numerical | confirmed |
| +1 branch, six δ from 0.0003 to 0.1 | **no scalar eigenvalue with μ² < 9/4**; normalized mismatch ≥ 0.99274 (a root needs 0) on a 360-point real scan of [−400, 2.2499]; below −400 only the conditional B > 0 bound applies | numerical | confirmed |
| Complex eigenvalues | none inside [−60,2]×[−30,30] and [−400,2]×[−200,200] | numerical | partially confirmed (subset re-run) |
| Analytic bound | if B > 0 the spectrum is real and μ² > −4; on the +1 branch B = 3.479–3.555 (original shell: B = −2.68×10⁻⁴) | conditional | confirmed |
| ℓ = 1 and ℓ = 0 | no physical mode | numerical | confirmed |
| Tensor | only the massless graviton; gap m ≥ 3H/2 | μ² ≥ 0 exact; no bound states numerical | partially confirmed (label split) |
| Planted instability (sign of σ″ flipped) | μ² = −28555 at δ = 0.001 | numerical control (tests only a strong tachyon, far outside the main scan) | confirmed |

The better measure of detection power is the auditor's margin test: a real mode with μ² < 0 would need the shell coefficient B inside [−0.659, 0.0081] (the range over all tested δ), while the actual B is 3.479–3.555 (`stability_gauge_invariant/verify/indep_spectrum.json`, key `Bstar_control`). Stability also rests on an assumed criterion: modes must be regular (normalizable) at the cone, which is standard in the dS-brane literature.

### 2.2 Method B — time domain and discrete operator (`stability_time_domain/`)

This method uses the unchanged registered PDE solver, linearized about the +1 branch, in the homogeneous (FRW-symmetric) sector. It computes eigenvalues directly and also evolves constraint-satisfying disturbances in time.

- Calibration: the time-domain rate on the original shell is 1.657193368 (reference 1.6571936312). **Status: numerical, confirmed.**
- On the +1 branch, every mode with support on the shell lies on Re λ = −3/2, the de Sitter continuum. None lies above it. **Status: numerical, confirmed.**
- The measured decay rate is −1.500 ± 0.001 in units of H₊, i.e. −0.910 H₀. **Status: numerical, partially confirmed:** the precision was overstated, and the "shell-localised" test profile was really a bulk bump.
- Positive rates exist only for (i) a pure-gauge mode (λ → 1), (ii) far-boundary artefacts (≈ 0.6/L), and (iii) grid-scale oscillations. **Correction from the audit:** the semi-discrete operator also has grid-scale modes with Re λ ≈ +0.245 near the shell. They do not shrink under refinement and are damped only by the RK4 time stepper. They are numerical artefacts, not physics, but the original headline wording was too strong.
- **Weak-mode sensitivity (audit caveat).** The evolved signal barely reaches the shell (peak |f_b| ≈ 8.5×10⁻⁹ of the bulk amplitude), so the time-domain runs cannot exclude a weakly bound, weakly excited shell mode. The no-bound-state conclusion of Method B rests on its spectral (eigenvalue) evidence, not on the time-domain decay.
- Planted instability (σ″ flipped): λ = 167.48925. **Status: numerical control, confirmed**; like Method A's control, it shows only that a strong tachyon is detected.
- Whole-domain constraints do not converge in the far stretched region. **Status: negative, confirmed.**

### 2.3 Do the two methods agree? (`synthesis/RECONCILIATION.json`, 20/20 checks, 4 deliberate wrong-formula controls fail as required)

| Overlap | A | B | Agreement |
|---|---|---|---|
| Original-shell growth rate | 1.657193631259 | 1.6571936573 (eigenvalue), 1.657193368 (time domain) | 1.6×10⁻⁸ and 1.6×10⁻⁷ relative |
| Planted tachyon (δ = 0.001) | rate 167.4892463 | 167.4892463 | 3.7×10⁻¹⁰ relative |
| Bound states below 9/4 on the +1 branch | none | none | agree |
| Slowest decay | continuum edge ⇒ −3/2 | −1.500 ± 0.001 (Richardson; −1.4994 to −1.5007 depending on the assumed order) | agree within the supported 1×10⁻³ |
| The λ ≈ 1 mode | ℓ = 1 harmonic (μ² = −4) is pure gauge | residual gauge mode, λ → 1 | λ = 1 ⇔ μ² = −4; consistent |

**Why they agree, and where they differ.** The scalar eigenvalue condition depends only on μ². Every μ² has a homogeneous representative with time dependence e^{λH₊τ}, where μ² = −λ(λ+3). Method B's homogeneous sector therefore tests the same scalar condition as Method A. Method A adds the tensor sector, complex eigenvalues and the analytic bounds. Method B adds an independent code path, the actual nonlinear right-hand side at small amplitude, and a real-time decay rate. Their positive discrete rates have no counterpart in A because they belong to the truncated grid, not the physics.

One further new cross-check: A's shell coefficient B extrapolates to 3.5555554 as δ → 0. The exact limit is 32/9 = 14/9 (the growth rate of the regular Δ = 18 mode) + 2 (σ″/2 at φ = 1). They agree to 3×10⁻⁸.

**Verdict (numerical, independently verified).** Within the scalar sector (all harmonics in Method A, homogeneous sector in Method B) and the tensor sector (Method A only), at the tested δ, the +1 branch is a linearly stable, empty, inflating brane. The only marginal mode is the massless graviton.

**Limits.** This is floating-point numerics, not interval-certified. The vector sector was not computed. Real μ² < −400 and complex μ² outside the two winding contours are excluded only by the conditional B > 0 theorems, and the large winding contour was not re-run by the auditor. Nonlinear stability, the basin of attraction, quantum tunnelling of the +1 AdS vacuum (U(+1) = −2/27 lies above U(−1) = −50/27), and whether the Chat 14 roll-off actually reaches this branch are all open. Stability of an empty brane says nothing in favour of a radiation era.

## 3. Results by workstream

### 3.1 Exact analytic structure (`analytic_structure/`)

- **Series through δ⁸, exact in c** (exact-verified; confirmed). The 22 Sept coefficients are reproduced. New examples:
  - h₃ = −c²(11767c + 11362)/2826240;
  - a₃ = −27c(60469c² + 49568c + 4048)/331612160 for φ_b − 1.

  All coefficients are in `SERIES_COEFFICIENTS.json`. Two independent multiprecision solvers (40–65-digit working precision) agree with each other to 10⁻³⁹–10⁻⁴¹ relative at δ = 10⁻³ to 10⁻² and with the series. For example, at δ = 0.001:
  - H² = 5.9240147943288795996237552780×10⁻⁵;
  - φ_b = 0.999915947316913353875784077511.
- **Correction (refuted claim).** The workstream said the δ⁹ coefficient needs global matching to the cone. The auditor showed it is local and computed it exactly: at registered c, the δ⁹ coefficient of φ_b − 1 is 3.249147254269×10⁻⁶ and that of H² is −2.759694789×10⁻⁷. Only the cone value η_h gets a δ¹⁶ ln δ term. The drop-out of the resonance constant was shown exactly at degree 9 only; higher resonances were not examined.
- **Why degree 14** (exact-verified). For superpotential potentials, the bulk mass at a vacuum gives the conformal dimension Δ = 3W″/W. At φ = +1, Δ = 18 and the regular static solution is the Gegenbauer polynomial C₁₄⁽²⁾(cosh y/9) with n = Δ − 4 = 14. At φ = −1 no polynomial exists. Closed forms are C₁₄⁽²⁾(cosh u) = Σ(j+1)(15−j)e^{(14−2j)u}; and for any bulk mass the solution is the elementary f_ν = (1/sinh u)·d/du[sinh(νu)/sinh u].
- **Leading-order tensor spectrum** (exact at leading order in δ). The spectrum is fixed by (μ − 3/2)P^{−μ}_{1/2}(cosh u_b) = 0, which gives a massless graviton plus a continuum from m = 3H/2. The gap itself is known literature (Garriga–Sasaki); the compact quantization form is written out here for this model.
- The decoupled bulk scalar has no bound state for H/k < 16.65 (conditional: it ignores the scalar–metric mixing that Method A includes).

### 3.2 Instant preheating on the archived roll-off (`preheating/`)

This is an explicit mode-function calculation for a shell field χ with m_χ² = m₀² + ḡ²(φ_b − φ*)², run on the archived Chat 14 +1 trajectory. It scans the coupling G = 10 to 10⁶, the crossing point φ* = 0.25, 0.5, 0.75, 0.9, and a decay Yukawa y ≤ 1.

- **Calibration** (numerical; confirmed, but the two wrong-formula controls are fixed by Gaussian algebra and so test the comparison metric, not the solver). The flat-space instant-preheating spectrum is reproduced to ≤ 8.4×10⁻⁷. An independent solver agrees on the real trajectory to ≤ 2.2×10⁻⁸.
- **Produced energy versus what a radiation era needs** (conditional; confirmed). This assumes the fixed archived trajectory (no backreaction), χ masses below the 5D gravity scale M₅ and instant full conversion. Then:
  - the produced energy is ≤ 1.7×10⁻⁴ of the leftover-vacuum threshold for couplings that pass the local screen (evaluated at the data end);
  - the maximum over the stored profile is 5.9×10⁻⁴;
  - over the whole scan it is ≤ 2.8×10⁻³;
  - with the generous "post-crossing" mass cutoff it is ≤ 4.0×10⁻³ for screen-passing couplings.

  Reaching the threshold needs χ masses ≥ 18 M₅ (≥ 6.3 M₅ with a generous relaxation), and stronger coupling makes it worse (R ∝ 1/√G).
- **Self-consistency** (conditional; confirmed). Wherever the channel could matter, its backreaction on the scalar junction is 10–79%. The fixed-trajectory calculation therefore stops being valid exactly there. The loop correction to the tension is 1.1–3.7 times σ′, so keeping the registered σ would need tuned counterterms.
- **Correction:** the decay-channel maximum is 2.34×10⁻⁴, not "≤ 1.3×10⁻⁴" (a stale README number). A derivation note was also corrected (the turning-point term is 15/(2π)). Neither changes the conclusion.
- **Verdict:** a conditional negative. The fully coupled bulk problem (all three modified junctions) remains unsolved.

### 3.3 Constraint control in 5D evolutions (`evolution_constraints/`)

This is a numerical-method result; it does not bear directly on the hypothesis.

- **Cause of the stall** (numerical; confirmed). The Hamiltonian-constraint stall of 22 Sept had an observed order of 0.24 at the finest grids. It comes from a sampled initial-data defect near the shell that does not shrink with the grid spacing. An exact transport law then carries it outward and amplifies it geometrically.
- **Remedy A: discrete-constraint projection of the initial data** (numerical; confirmed). It restores convergence, with orders ≈ 2.9–3.1 at the finest pair, and lowers H_max(t = 0.5) from 6.38×10⁻³ to 9.35×10⁻⁴. Caveat from the audit: the projection correction itself (~10⁻⁹) does not converge below h = 10⁻⁴.
- **Remedy B: outgoing-characteristic damping** (exact transport law, verified). Combined with Remedy A it gives 2.28×10⁻⁵. **Correction:** the reduction factor is about 25–200 at t = 0.5–1 (≈ 8 at t = 0.25), not "13–100".
- **Fourth-order convergence** is shown on smooth calibration data. For the actual seed it is not reached (inconclusive): φ_t and B_t converge at about 2, on the shell-seeded front.

### 3.4 Mechanism screen for a radiation era (`mechanisms/`)

| Candidate | Result | Status after audit |
|---|---|---|
| Bulk Weyl ("dark radiation") made by the blast | The blast does make it: peak 0.0806–0.0809 H₀² (≈ 15% of H²) on 3 grids; late value not converged. It cannot be the hot Big Bang, and N_eff caps it at a few percent of the radiation (N_eff inputs are snippet-level) | numerical (confirmed); negative (partially confirmed) |
| 4D energy budget with graviton zero-mode Planck masses | Leftover vacuum H_vac²/H₀² → 0.368 as δ → 0 (0.372 at δ = 0.01, 0.406 at δ = 0.1); radiation can dominate for ≤ 0.59 e-folds (0.47 at δ = 0.1) vs ≈ 21.8 needed | conditional (confirmed) |
| Dissipative shell coupling, 5D pilot (δ = 0.1) | captures 0.37–2.2% of the tension drop; R/R_crit ≤ 0.038 | negative (confirmed) |
| Sudden tension→radiation conversion | gravitationally screened: the Weyl term must cancel 80–99.9% | exact-verified |
| Gravitational particle production | short by 85–122 orders of magnitude | negative, order of magnitude |
| Collapse branch as an ekpyrotic phase | w_eff = 0.76–0.85 over the ~1.1 resolved e-folds | inconclusive (resolved range only) |
| Dark-bubble cosmology | not available in the registered model | conditional (rests on a cited theorem and the Z2 topology) |
| Tuning the vacuum with c | **Refuted as first stated.** The family does continue past c = 0 to c* = −0.99307 (δ = 0.1). The new starting shell has φ_b = −1.7287 and, in one exploratory run (single grid, one seed sign, no radiation), rolls across φ = −1 to φ_b = 0.80 within the reliable window | corrected; endpoint open |
| Tuning with a quadratic tension term d* ≈ −3.1 (model change; pilots at δ = 0.1) | radiation reaches 69% of H², inside a budget with large cancellations (vacuum −117%, Weyl +121%); Weyl/radiation is 1.76 at the stop and 1.23 at the end of the extended reliable window, about 12× the ΔN_eff-derived limit of 0.1. After the chart freeze the ratio is not usable (it falls through 0.6 and below zero while proper time barely advances, then swings back as H turns negative; `CRITIC_CHECKS.json`). Convergence rests on one grid pair | inconclusive (partially confirmed; stop cause corrected) |

A further conditional result: if the leftover vacuum is identified with today's dark energy, the blast e-folds on a 6.4 Gyr timescale and δ ≈ 1.1×10⁻⁶² to 7.6×10⁻⁶⁸ (for AdS radii from 38.6 µm down to 0.1 µm). That is the brane-world form of the cosmological-constant problem, not an explanation of it.

### 3.5 Literature, 2022–2026 (`literature/`, merged in `LITERATURE_2022_2026.md`)

Scope: 162 distinct sources, 82 dated 2022–2026. **Every entry is snippet-only.** Paper sites were blocked, and the methods sweep could not run new searches.

- **Prior art** (context). "Our universe as a wall born in a 5D event" is an established published line: brane-world creation, ekpyrotic collisions, holographic white holes, and the 2018–2026 dark-bubble programme, including black-hole-catalysed creation. HDBLAST must not claim the concept; its content is the registered model and its computed results.
- **The pattern of hot-origin models** (context). Every published hot-origin model adds an ingredient: a bulk black hole, colliding walls, a hot holographic sector, or matter on the wall. This round's negative results for a pure-tension wall fit that pattern.
- **Registered model versus standard formulas** (exact; self-checked with controls, no separate auditor). The balanced tension 2W equals the critical Randall–Sundrum/Karch–Randall tension at both vacua. At O(δ²) the thin-brane H² differs from the +1 branch by exactly −c²δ²/384; this term predicts −7.849 ppm of the measured −7.869 ppm shift in H, and higher orders supply the rest.
- **Observations** (constrain, do not test). The model has no radiation era and no perturbation spectrum, so CMB, BBN, DESI, JWST and H₀ data are future requirements, not tests.
  - The registered PTA knee sits at 0.9972 × (1/yr), on the pulsar position/proper-motion fitting blind spot. That blind spot is standard in PTA work but new to this project.
  - No 2025–2026 release tests the knee at likelihood level. Noise reanalyses trend towards γ = 13/3, which is unfavourable to the frozen shape but not a rejection.
  - Tabletop gravity bounds would put the static branch at a TeV-scale vacuum, with its horizon-scale gravitational-wave signal in the sub-mHz band, not the nHz band (conditional).
  - The pulsar-knee work has still never been derived from the 5D model (program map).

### 3.6 Program map (`program_map/PROGRAM_MAP.md`)

The map is a dated chronology from August 2025 to 23 September 2026 with a claims ledger (established / negative / withdrawn / open) and 63 key paths. Its self-verifier passes 34/34 checks. It was not separately audited. The version and file list of record 22922928 are unverified, because zenodo.org is blocked here.

## 4. Corrections made in this round

1. **δ⁹ series coefficient:** the claim that it is non-local is refuted. It is local, with exact values given in §3.1.
2. **Vacuum tuning with c:** the claim that the family ends at c → 0⁺ is refuted. The family continues to c*, where a new starting shell exists; its fate is open.
3. **Preheating decay maximum:** the stale value ≤ 1.3×10⁻⁴ is replaced by 2.34×10⁻⁴.
4. **Smaller corrections:**
   - the Method B headline wording and the precision of its decay rate;
   - the tensor-claim label;
   - the ρ_b agreement (8.5×10⁻¹⁴);
   - the damping factors;
   - the normalisation of the transport gain;
   - the dark-bubble exclusion relabelled as conditional;
   - the collapse-branch test limited to the resolved range;
   - the stop cause of the d* pilot.

   All are itemised in `CLAIM_REVIEW.md`.
5. **Critic review of this report (28 September 2026).** A final completeness and overclaim pass over the top-level documents tightened wording that was stronger than the audits support, added missing caveats (the δ = 0.1 pilots, the scan range and cone-regularity assumption of Method A, the weak shell excitation of Method B's time-domain runs, the fixed-trajectory condition of the preheating bound), withdrew an unsupported post-chart-freeze number for the d* pilot, completed the list of unverified citations, and removed AI product names. The changes and the checks behind them are listed in `CLAIM_REVIEW.md` (section "Critic review") and `CRITIC_CHECKS.json`.
6. **Cross-workstream reconciliation.** The theory sweep says a pure-tension brane has no dark radiation. That is true for static configurations. The dynamical roll-off does make a positive Weyl term (§3.4).

## 5. What this means for the hypothesis

- **What stands.**
  - The registered 5D model has an unstable starting shell (certified earlier).
  - It has a linearly stable, empty, inflating end state (numerical, scalar and tensor sectors, this round).
  - A roll-off connects them in the Chat 14 runs, as a tendency, not a proven attractor.

  This is a consistent five-dimensional story of a "blast" that ends in de Sitter expansion.
- **What fails.** Within the registered model and every mechanism screened so far, that story does not contain a hot, radiation-dominated era. The central obstruction is quantitative: the leftover vacuum H_vac²/H₀² ≈ 0.37 is too large for any screened heat source to overcome.
- **What is open.** Surviving requires a model change that cancels the leftover vacuum without creating too much dark radiation, or an added ingredient such as a bulk black hole or a collision. Both are standard in the literature and would need to be justified, not just assumed. The tuned-vacuum pilot is the most promising lead, and it is currently failing the dark-radiation test.
- **Observational status.** HDBLAST currently makes no prediction that 2023–2026 data can test. The registered PTA knee is a separate phenomenological template, not derived from the 5D model.

## 6. The discriminating next tests (details in `NEXT_TESTS.md`)

1. **Tuned-vacuum roll-off past the chart freeze.** Evolve the d* model in a proper-clock gauge with constraint projection and damping. Decision: HDBLAST with d* survives only if the Weyl-to-radiation ratio settles at ≤ 0.1 while radiation dominates (the normalisation of this N_eff-derived limit and its inputs must be fixed and checked against the papers first; see `NEXT_TESTS.md` A1).
2. **The c* shell.** Continue the signed shell family to c* at δ = 10⁻² and 10⁻³ (so far it exists only at δ = 0.1), then run the c* starting shell on several grids to find where it ends.
3. **Nonlinear attraction.** Does the Chat 14 +1 run relax onto the +1 branch? If it relaxes linearly, deviations should decay at least as fast as e^{−0.91 H₀τ} (e-folding time 1.10 in H₀τ).
4. **Bulk black hole.** Add a bulk black-hole (Weyl) parameter. This is the cheapest literature-backed hot ingredient, and it must pass the N_eff bound.
5. **Collapse fate.** Search the Chat 14 collapse runs for trapped surfaces and an apparent horizon.
6. **Completing the stability result:**
   - vector sector;
   - interval certification (the margin is large: the normalized mismatch is ≥ 0.99);
   - Euclidean-action comparison of the two static branches;
   - tunnelling of the +1 vacuum.
7. **Observations.** Pre-register how the PTA knee handles the 1/yr degeneracy before NANOGrav 20-yr and IPTA DR3. Derive a perturbation spectrum before any CMB claim.

## Package map and reproduction

| File | Content |
|---|---|
| `RESULTS_SUMMARY.json` | 66 claims with producer label, verifier verdict, status and evidence paths |
| `CLAIM_REVIEW.md` | Claim-by-claim table, corrections and the stability reconciliation |
| `LITERATURE_2022_2026.md` | Merged, deduplicated annotated bibliography |
| `NEXT_TESTS.md` | Prioritised next calculations and observational tests |
| `REPRODUCE.md`, `requirements.txt` | Exact commands per workstream |
| `PUBLIC_SUMMARY.md` | One-page summary for GitHub visitors |
| `synthesis/` | Reconciliation, quoted-number checks (56/56 plus 3 controls), builders |
| `critic_corrections.py`, `CRITIC_CHECKS.json` | Post-build critic corrections to `RESULTS_SUMMARY.json` and `CLAIM_REVIEW.md`, plus checks of the numbers, citations and wording added in the critic review |
| `MANIFEST.sha256.json` | SHA-256 of every file in this folder |
| workstream folders | Each of the six computational workstreams has `README.md` (producer) and `VERIFICATION.md` plus `verify/` (independent audit); `literature/` and `program_map/` have self-checks only. Workstream READMEs keep their original wording; where an audit corrected them, the corrected statement is the one in this file and in `CLAIM_REVIEW.md` |

Every number above is traceable to a JSON output. `synthesis/check_quoted_numbers.py` checks the key ones automatically, and `critic_corrections.py` checks the numbers added in the critic review (`CRITIC_CHECKS.json`). Original research folders under `new-files/` were only read, never modified.
