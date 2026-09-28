# Next tests for HDBLAST, in priority order (after the 27 September 2026 checkpoint)

Prepared with AI assistance for Ricardo Maldonado's HDBLAST program.

The tests are ordered by how much each outcome would change the status of the hypothesis. Each one has a pass/fail rule fixed *before* it is run, where possible. The cost estimates assume the current machine (4 cores; runs of ≤ 20 min each on ≤ 2 cores).

The starting point is the checkpoint's status:
- a certified unstable starting shell;
- a linearly stable (numerical; scalar and tensor sectors), empty, inflating end state (the +1 branch);
- no radiation era in the registered model;
- two open, model-changing leads.

## A. Decisive calculations (could revive or retire the hot-Big-Bang reading)

### A1. Tuned-vacuum roll-off past the chart freeze (highest priority)

**Model.** σ = 2W + δ(1 + cφ + dφ²/2), with registered c and d = d*(δ) from `mechanisms/M8_QUADRATIC_TENSION_TUNING.json`. This is an explicit change to the registered model.

**Why it matters.** It is the only screened configuration where shell radiation becomes a large share of H² (69% in δ = 0.1 pilots, inside a budget with large cancellations: vacuum −117%, Weyl +121%; `mechanisms/verify/V2_PILOT_CONSISTENCY.json`). But the bulk Weyl ("dark radiation") term stays larger than the radiation:
- 1.76× at the original stop;
- 1.23× at the end of the extended reliable window.

Beyond the chart freeze the ratio is not usable: proper time barely advances while the ratio falls through 0.6 and below zero and then swings back as H turns negative (`CRITIC_CHECKS.json`, key `post_freeze_dstar`). That is why a proper-clock gauge is needed.

The N_eff bound needs ≤ 0.03–0.1 in the normalisation ρ_DR/ρ_SM at production (mapping in `mechanisms/M5_OBSERVATIONAL_SCREEN.json`; the input N_eff values are snippet-level, one of them, arXiv:2603.13226, has no record in the literature files, and the ACT DR6 value is quoted differently in two sweeps). The observations sweep quotes ρ_dr/ρ_γ ≲ 0.016 (ACT DR6) at recombination, a different normalisation. Fix the normalisation and verify the inputs against the papers before running.

**Code.** Start from `mechanisms/pilot_5d/rolloff5d_matter.py`, which already has all three modified junctions, the radiation ledger and `--d`. Then:
1. Remove the φ_b > 0.995 stop, since the tuned endpoint is at φ_b ≈ 1.036.
2. Replace the near-shell conformal slicing with the proper-clock gauge from folder 152.
3. Add the missing dR/dt and source terms to `neumann_t`.
4. Use Remedy A (discrete-constraint projection) and Remedy B (switched-off characteristic damping) from `evolution_constraints/lab.py` to control the constraint error.

**Runs.**
- δ = 0.1 with Y ∈ {0.3, 1, 3}, on two shell spacings and two seeds.
- δ = 10⁻³ on the Chat 14 refined grid (all pilots so far are at δ = 0.1, and d* convergence rests on one grid pair).

**Decision rule, fixed now:**
- **Pass:** 𝒲a⁴ plateaus with 𝒲/(radiation terms) ≤ 0.1 (≤ 0.03 for the combined bound) while Ω_r → 1. The tuned model then has a candidate radiation era, and the next step is to replace the friction closure with the derived χ source.
- **Fail:** 𝒲 plateaus above that ratio. The d* model is then excluded by N_eff, and no screened mechanism survives.

Either result is publishable as a checkpoint.

**Caveat.** A pass would still leave the tuning 1 + c + d/2 ≈ 0 unexplained: it is the cosmological-constant problem restated.

### A2. Fate of the c* starting shell

The auditor found that the registered linear tension *can* be tuned (`mechanisms/verify/V3B_CONTINUE_TO_CSTAR.json`). At c* = −0.99307 (δ = 0.1) a static shell exists with φ_b = −1.7287. In a single-grid Y = 0 pilot it rolls across φ = −1 to φ_b = 0.80 (`V5_CSTAR_PILOT_SIGNED.json`).

**What to do:**
1. Repeat the continuation at δ = 10⁻² and 10⁻³, with the signed cone parametrisation.
2. Certify the shell root, reusing the M462-type pipeline.
3. Compute its tachyon.
4. Run the roll-off on three grids and with both seed signs, with Y = 0 and Y > 0.

**Decision:** does it end on a vacuum-free brane? If so, compute the Weyl-to-radiation ratio as in A1.

### A3. Bulk black-hole (Weyl) ingredient

This is the cheapest literature-backed "hot" ingredient (Garriga–Sasaki thermal instanton; Kraus; Savonije–Verlinde; dark-bubble black-hole catalysis; all snippet-level). Add an AdS–Schwarzschild mass behind the shell to the static and roll-off solvers. Then test:
- whether the +1 branch persists;
- which a⁻⁴ term appears on the brane;
- whether the BBN/CMB dark-radiation bounds (ρ_dr/ρ_γ ≲ 0.016–0.1, snippet-level) can be met while a Standard-Model radiation sector dominates.

This is a new model branch, and it should be registered before any runs.

## B. Consolidating what is already established

### B1. Nonlinear attraction to the +1 branch

Does the Chat 14 positive-side run approach the +1 branch? The linear spectrum predicts that deviations decay at least as fast as e^{−(3/2)H₊τ} = e^{−0.910 H₀τ} (e-folding time 1.10 in H₀τ units; `synthesis/RECONCILIATION.json`, conditional on linear relaxation). The late H/H₀ ≈ 0.6386 at H₀τ = 6.9 (22 Sept report) is still 5% above 0.6067.

**What to do:** re-run with constraint projection, a compactified or characteristic outer region (the far end is a Rindler-type horizon: `literature/methods.md`), and full-state snapshots.

**Test:** fit the decay of H/H₀ − 0.6067. Also check that the cone region actually moves from near φ = −1 to near φ = +1, which requires following the departing wall.

### B2. Complete the linear stability result

- **Vector sector.** Graviphoton-type modes were not computed.
- **Interval certification.** Certify that the scalar mismatch has no zero on −4 < μ² < 9/4. The margin is large (normalized mismatch ≥ 0.99; B ≈ 3.5 against an instability threshold B* ∈ [−0.66, 0.008]), so a coarse validated ODE enclosure should suffice. Use the η = φ − 1 formulation: the cone germ is amplified by 6.3×10¹⁸ at δ = 10⁻³ (`literature/methods_checks/certification_pilot.json`).
- **Existence and local uniqueness.** Certify the +1 branch, with the exact series through δ⁸ (δ⁹ from the audit) as the proof centre.

### B3. Quantum and global stability

- Compare the Euclidean on-shell actions of the original shell and the +1 branch. This decides which configuration a nucleation event prefers.
- Estimate Coleman–de Luccia / brane-nucleation decay of the +1 AdS vacuum in the presence of the detuned shell. U(+1) = −2/27 lies above U(−1) = −50/27, and positive-energy protection has only been shown without branes.

### B4. Collapse fate

Search the Chat 14 collapse runs for trapped surfaces and an apparent horizon, which the literature expects to form (black brane or crunch). Evaluate the Blanco-Pillado et al. criterion comparing the expansion rate with the fifth-dimensional momentum transfer. If the full state is missing from the archive, re-run with full snapshots.

## C. Matter-sector consistency (needed before any radiation claim)

- **C1. Coupled preheating.** Where instant preheating could matter (R ≈ 1), backreaction is 10–79% (`preheating/SUMMARY.json`). Solve the coupled problem, feeding all three modified junctions back into the bulk, with an adiabatically renormalized stress tensor. Address the Coleman–Weinberg tension shift (1.1–3.7 × σ′).
- **C2. Equivalence principle and fifth force.** The χ mass depends on φ, so brane matter couples directly to the bulk scalar. A 2025 dark-bubble paper (snippet) reports severe equivalence-principle violation from such couplings. This needs an estimate.
- **C3. External solver calibration.** Reproduce the exactly solvable Koyama–Takahashi dilatonic brane model, or a BraneCode test. This would be the first external calibration of the HDBLAST 5D solver.

## D. Observational tests

The registered 5D model currently has no radiation era and no perturbation spectrum. The items below are either requirements or tests of the separate PTA-knee template.

- **D1. PTA knee: pre-register the 1/yr degeneracy.** The knee is at 0.9972 × (1/yr), on the sky-position/proper-motion blind spot. Before NANOGrav 20-yr or IPTA DR3 data appear, fix in writing how the analysis propagates timing-model absorption near 1/yr. Also fix the likelihood-level test to be run: free-spectrum posteriors, not power-law summaries. Otherwise the knee's position is not identifiable.
- **D2. Derive, don't assume, the link to 5D.** The knee has never been derived from the 5D model. Under tabletop scale-setting, the static branch's horizon-scale gravitational-wave signal would be in the sub-mHz band, 1,600–2,600 × f_knee (conditional; `literature/observations_checks/knee_and_scale_checks.json`). Either derive a mechanism that puts a 5D signal at nHz, or treat the knee as an independent phenomenological test.
- **D3. N_eff and dark radiation.** Any surviving HDBLAST variant with a bulk Weyl remnant must satisfy CMB and BBN dark-radiation bounds. Current and next-generation CMB experiments (for example ACT, SPT-3G and Simons Observatory) will tighten these. Verify the exact bounds from the papers before use.
- **D4. Short-range gravity.** The AdS radius ℓ is bounded by torsion-balance tests (Yukawa range < 38.6 µm, snippet-level). Any HDBLAST prediction of micron-scale gravity changes should be stated as a falsifiable number once physical scales are fixed.
- **D5. Primordial spectrum.** Derive scalar and tensor perturbations for any surviving variant before comparing with n_s, r, JWST or H₀. Until then those data are requirements, not tests.

## E. Housekeeping

- Open the key papers (Garriga–Sasaki; Frolov–Kofman; Kiritsis–Nitti et al.; the dark-bubble series; ACT DR6; NANOGrav 15-yr papers) and replace snippet-level citations with checked ones.
- Run the 20 prepared methods queries (`literature/methods_checks/methods_sources.json → queries_planned_not_run`) once a search budget is available.
- Restore `analytic_structure/MANIFEST.sha256.json` by deleting its 25 `verify/…` entries. The 29 original entries still match the current files (`synthesis/QUOTED_NUMBER_CHECKS.json`).
