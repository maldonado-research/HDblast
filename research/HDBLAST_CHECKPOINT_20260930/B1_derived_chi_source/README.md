# B1: the derived χ source in place of the friction closure (tuned model, δ = 0.1)

Ricardo Maldonado's HDBLAST program, workstream B1, 30 September – 1 October 2026.
Status: internal, not peer reviewed. Every number below comes from a saved script with JSON output in this folder.
Labels used: **exact-verified**, **numerical**, **conditional**, **negative**, **inconclusive**.

## Post-audit corrections (1 October 2026)

This README was written before the independent audit. The sentences below were corrected in place to agree with `AUDIT.md` and with the final labels in `../CLAIM_REVIEW.md`. The registered verdict (INCONCLUSIVE at λ_c = 1) is unchanged. No run was repeated for these edits.

- **"Essentially final" withdrawn.** The b_cut values r = 94 to 8.5×10³ (y = 1) are run-end values, not registered plateau values. The cell φ\* = 0.9, G = 100, y = 1 is not fully decayed at the end (Ω_χ/Ω_r ≈ 0.095, Ra⁴ still rising 16% over the last Δln a ≈ 0.12), and its r is about 1.9× high because the field-space ramp under-produces (§6.2, §6.4).
- **y = 0.1 b_cut values are snapshots, not lower bounds.** Undecayed χ carries 0.6 to 44 times the radiation energy at the end. A linear "all χ → R" estimate would bring r down to about 52–4.4×10³ (§6.2). They are still ≫ 0.1.
- **Constraint bound.** "≤ 1.2×10⁻⁴" holds on the fine grid only. The coarse-grid maximum is 4.2×10⁻⁴.
- **Production-model controls were mis-paired (§6.5).** The values −338, −1149 and 22.5 are r at each control's *reliable end* (H₀τ = 1.85, 1.57, 2.62), while the baseline 0.17 is r at the *run end* (H₀τ = 5.19). Compared at the run end on both sides, the controls give r = 0.316 (ramp 1.5), 0.082 (instant) and 0.169 (no D) against 0.170: about 2× sensitivity, and no effect from D. The corrected values come from `audit/AUDIT_TIMESERIES.json`. **`B1_RESULTS.json` (`production_model_controls`, `last_r`) still holds the mis-paired values and has deliberately not been edited.** All run-end values lie beyond the runs' reliable ends.
- **Label of the main conclusion** changed from "negative (conditional)" to **"negative (conditional/descriptive)"**; the registered outcome is INCONCLUSIVE.
- **Registration note timing.** The timing of the second dated note relative to the start of the grid batch cannot be verified from timestamps (git was not used).

See `AUDIT.md` for the checks and `audit/*.json` for the recomputed values.

## Short answer

- **Registered verdict: INCONCLUSIVE at λ_c = 1** (registration §4).
  - All 16 grid cells ran at both spacings. Not one reached the A1 plateau, which needs Δln a = 0.5 with 𝒲a⁴ and Ra⁴ steady.
  - Two different things stop the runs:
    - **b_cut cells:** the solver chart ends the evolution at H₀τ ≈ 7.4, only Δln a ≈ 0.2–0.4 after the roll-off and the decay.
    - **b_need cells:** near-shell constraint residuals rise above 0.05 shortly after production, and the runs then stop non-finite at H₀τ ≈ 4.4–5.6.
  - No rule was relaxed.
- **What the evidence says, short of a registered classification.**
  - At the M5 cutoff (λ_c = 1) the derived channel makes far too little radiation. The roll-off's own Weyl (dark-radiation) term dominates by 2 to 5 orders of magnitude.
  - The 5D runs at b_cut give r = 𝒲/rad = 94 to 8.5×10³ for y = 1, and 3×10² to 1.1×10⁵ for y = 0.1. Ω_r ≤ 0.012 throughout.
  - The two spacings agree to better than 0.1%. Constraints stay ≤ 1.2×10⁻⁴ on the fine grid (≤ 4.2×10⁻⁴ on the coarse grid).
  - These are run-end values, not registered plateau values: the longest window meeting the A1 tolerance is Δln a = 0.394 < 0.5. For three of the four y = 1 cells χ has decayed and 𝒲a⁴ and Ra⁴ are close to flat at the end. The cell φ\* = 0.9, G = 100 is not fully decayed (Ω_χ/Ω_r ≈ 0.095, Ra⁴ still rising 16%), and its r is about 1.9× high from ramp under-production. The y = 0.1 values are snapshots, not lower bounds (undecayed χ; r could fall to about 50 in the most favourable cell). All values remain ≫ 0.1.
  - The fixed-background reduced estimate (conditional) gives FAIL-Weyl in all four y = 1 cases (r = 97–8.8×10³). It gives recollapse (FAIL-no-radiation-era) in all four y = 0.1 cases.
  - In this test-field regime the reduced estimate matches the 5D runs to 3–4% in r (φ\* = 0.5 and φ\* = 0.9 with G = 1000). The one exception is a factor of 1.9 from a production-ramp artefact (see Limitations).
- **What it would take.**
  - The reduced estimate needs λ_c ≈ 11–50 for r ≤ 0.1 and ≈ 17–74 for r ≤ 0.03, with the cutoff taken over the whole trajectory. That means χ masses 10 to 75 times above M₅.
  - At those b the backreaction is O(1) (R_j = 0.35–1.9). The 5D runs there are unreliable, and they show the reduced estimate is not valid there: R is 5–25 times lower than the fixed-background value.
  - One exploratory, unregistered 5D cell at λ_c ≈ 5.4 (b = 10⁻³; φ\* = 0.5, G = 100, y = 1) is reliable throughout. It settles at r = 0.282 with Ω_r = 0.82: 𝒲a⁴ is flat to 10⁻⁴ and Ra⁴ to 0.5% over its last Δln a ≈ 0.18. That is still above the conservative bound.
- **Conclusion (negative, conditional/descriptive, within this model at δ = 0.1; the registered outcome is INCONCLUSIVE).** With χ masses kept below M₅, the derived χ channel cannot replace the friction closure: the roll-off's own dark radiation dominates by a large margin. Whether some b above the cutoff gives r ≤ 0.1 is not settled. Reaching r ≈ 0.28 already needs λ_c ≈ 5, and the strong-backreaction regime beyond that could not be evolved reliably.

## 1. Question

A1 (28 Sept) found a radiation era with r = 0.073 in the tuned model d = d\*(0.1), but only with the phenomenological friction closure κ₅²j = Yv at Y = 1. B1 replaces that closure with the χ-particle source derived in the 22 Sept matter extension, in which the bulk-scalar-dependent χ mass is

  m_χ² = m₀² + G²(φ_b − φ\*)².

χ is produced at each crossing of φ\* by instant preheating and decays into radiation at rate Γ = y²m/(8π).

The question is whether this channel yields a radiation era with r = 𝒲/rad ≤ 0.1 (and ≤ 0.03) by the A1 rule, with masses kept below the M₅ cutoff, b = κ₅²H₀³ ≤ (λ_c/Λ̂)³.

In the tuned model the residual vacuum energy is ≈ 0, so the competitor is the roll-off's own Weyl term, not the residual-vacuum threshold that killed the 27 Sept preheating channel.

## 2. Restart and reuse

This round ran in three sessions; two were interrupted by usage limits.

**Reused from the earlier sessions, after review:**
- `a1_copy/`: byte-identical copies of the A1 code.
- `shell_matter.py`, `evolve_b1.py`.
- `b1_symbolic.py` → `B1_SYMBOLIC_CHECKS.json`.
- `b1_reduced.py` → `B1_REDUCED.json`.
- `b1_analyze.py`.
- `REGISTRATION.md` with its first dated note. It was written in the second session, before any grid run.
- Job lists `jobs_cal.txt` and `jobs_grid.txt`, and `queue.sh`.

**Kept but not used for any result:**
- `dev/`: development runs. Their outcomes are disclosed in the registration notes.
- `superseded_time_ramp/`: the first grid cell made with the time-ramp production, plus the K1/K2 runs made with an older solver build.

**Done in this session:**
- Second and third dated notes in `REGISTRATION.md`.
- Recalibration with the final build (K1, K2).
- The full registered grid, the production-model controls, the analysis, and the cross-checks `compare_reduced_5d.py` and `diag_blowup.py`.
- The exploratory cell (`jobs_explore.txt`, `explore_analyze.py`) and this README.

## 3. Registration summary (`REGISTRATION.md`; three dated notes appended; the audit found the body consistent with being unedited, but this cannot be verified without version control, nor can the timing of the second note)

**Model change B1.**
- Tuned A1 model: σ = 2W + δ(1 + cφ + dφ²/2), with δ = 0.1, c = 0.5975949350280132 and d = d\* = −3.106933495673783.
- The friction closure is removed (Y = 0). The shell instead carries a χ particle gas, with m₀ = 0, decaying into radiation R.
- Seed: dc = 10⁻² only.

**Rule.** The A1 rule, applied verbatim, with R taken as the decay radiation only:
- r = 𝒲/rad, with rad = σR/18 + R²/36.
- Plateau: over Δln a = 0.5, H > 0 and both 𝒲a⁴ and Ra⁴ steady to 5%.
- PASS-conservative: Ω_r ≥ 0.9 and |r| ≤ 0.1. PASS-combined additionally requires |r| ≤ 0.03.
- FAIL-Weyl, FAIL-no-radiation-era and INCONCLUSIVE as in A1.
- Reliability: both spacings, dz_f = 10⁻³ and 5×10⁻⁴, give the same class and agree on r, and near-shell constraints stay < 0.05.

**Grid.** φ\* ∈ {0.5, 0.9} × G ∈ {100, 1000} × y ∈ {0.1, 1} × b ∈ {b_cut(λ_c = 1), b_need}.
- b_need comes from the reduced estimate. For y = 0.1 the y = 1 b_need is used (first dated note).

**Aggregate.**
- B1-PASS needs a reliable PASS in some b_cut cell.
- B1-FAIL needs a reliable FAIL in every completed b_cut cell, and the reduced estimate failing for all eight.
- Anything else is INCONCLUSIVE.

**Dated notes.**
1. b values fixed; y = 0.1 cells run at the y = 1 b_need.
2. Restart disclosure, and production inserted by a C² ramp in *field space* (with a stall switch) instead of in time. Reason: the time ramp's force j_prod = P/v diverges when the scalar stalls.
3. Informative additions, announced before they were run: the reduced-versus-5D comparison, the blow-up diagnostic, and the exploratory cell at b = 10⁻³.

## 4. The closed shell system (exact-verified: `b1_symbolic.py` → `B1_SYMBOLIC_CHECKS.json`, 22 identities and 8 controls, all as required)

All quantities are in H₀ units, with s = H₀τ and κ₅²ρ/H₀ = b ρ̂.

**χ gas.** Per momentum node i (24 generalised Gauss–Laguerre nodes per cohort):
- n_i = N_i a⁻³, p_i = k_i/a, ω_i = √(p_i² + m_χ²).
- ρ_χ = Σ n_iω_i and p_χ = Σ n_ip_i²/(3ω_i).
- **j = −∂L_m/∂φ = Σ n_i ∂ω_i/∂φ = G²(φ − φ\*) Σ n_i/ω_i** (C1, C1b).
- Ward identity: ρ̇_χ + 3H(ρ_χ + p_χ) = jv (C2).
- Wrong-sign control: j → −j leaves the residual 2jv (C2x).

**Decay.** Ṅ_i = −Γ(m/ω_i)N_i, time-dilated per mode, with Γ = y²m/(8π).
- Q = Γ m Σn_i.
- Ledgers: ρ̇_χ + 3H(ρ_χ + p_χ) = jv − Q and Ṙ + 4HR = Q. Q cancels in the total (C3, C3b).
- Controls: without the time dilation the ledger is unbalanced (C3x).

**Production.** At each crossing of φ\*:
- n_k = (1 + D/q) exp(−πk²/q), with q = G|dφ_b/ds| at the crossing.
- D is the O(1/q) curved/expanding correction of `HDBLAST_CHECKPOINT_20260927/preheating`, from quartic fits of the shell history over the preceding 0.2 H₀⁻¹. It is capped at |D/q| ≤ 0.3; the realised |D/q| was 0.002–0.036.
- The cohort is filled by a C² ramp in |φ − φ\*| of length 3|v\*|/√q.
- The work is paid by the scalar through the bounded force j_prod = P/v. The Ward identity then holds exactly, including the stall switch (C4c, C4d).
- Control: instant insertion leaves unpaid energy P (C4x).

**Junctions** (22 Sept matter extension), with ρ = ρ_χ + R and p = p_χ + R/3:
- n·A = (σ + κ₅²ρ)/6
- n·B = (σ − κ₅²(2ρ + 3p))/6
- **n·φ = −(σ′ + κ₅²(j + j_prod))/2**

Checks on the junctions:
- They give the signed shell budget d(λ + ρ)/dτ + 3H(ρ + p) = −2(n·φ)v/κ₅² (C5b).
- The + sign, or dropping the factor ½, breaks that budget (C5x, C5y).
- The Weyl balance holds with the χ source (C6). Controls: 4H → 3H, and dropping J (C6x, C6y).

**Cutoff.** Λ̂ = max(m̂_χ over the evolution, √q) and b_cut = (λ_c/Λ̂)³. The "post" variant uses masses after the crossing only.

## 5. Method

**Solver.** `evolve_b1.py` is the A1 solver (`evolve_a1.py`) with these settings:
- bounded chart, 4th-order differences, RK4;
- Remedy A projection on Z ≥ −15.9 and Remedy B damping κ = 10;
- L = 16, T_f = 14;
- x_c from a coarse pre-run with the same matter.

The only change is that the shell matter is the pluggable sector in `shell_matter.py`:
- `friction`: the A1 code path;
- `toy`: a friction-equivalent fluid decaying into R;
- `chi`: the derived source.

**Calibration (numerical, `B1_RESULTS.json` → `calibration`), with the final build:**

| Check | Requirement | Result |
|---|---|---|
| K1: χ path on with G = 0, against A1 `main_dstar_Y0_dc1e-2_dzf1e-3` | φ_b, H, 𝒲 within 10⁻⁸; same stop | **PASS**: identical to 0.0 on 305 records; same stop (non-finite at H₀τ = 7.3526) |
| K2: friction-equivalent toy (Y = 1, decay rate 10 H₀), against A1 `main_dstar_Y1_dc1e-2_dzf1e-3` | same class; r within 20% or 0.02 | **PASS**: PASS-conservative in both; r = 0.0738236 in both; φ_b, H, 𝒲 within 1.5×10⁻¹¹ |

**Reduced estimate (conditional, `B1_REDUCED.json`).**
- The shell ODE runs on the archived A1 tuned Y = 0 trajectory, with no backreaction.
- Beyond the archive, φ, 𝒲a⁴ and the residual vacuum are frozen.
- Archive robustness: the dzf = 10⁻³ and dc = 10⁻⁴ archives change r(b_cut) by ≤ 0.9% and b(r = 0.1) by ≤ 0.1% (`reduced_grid_spread`).

## 6. Results

### 6.1 Reduced estimate (CONDITIONAL; primary archive, "lowE" variant)

| φ\* | G | y | q | b_cut (λ_c = 1) | class at b_cut, r | b for r ≤ 0.1 → λ_c (full / post) | b for r ≤ 0.03 → λ_c (full) | R_j max at b_need |
|---|---|---|---|---|---|---|---|---|
| 0.5 | 100 | 1 | 128 | 6.5×10⁻⁶ | FAIL-Weyl, 97 | 0.010 → 11.5 / 11.5 | 0.032 → 17.0 | 0.69 |
| 0.5 | 1000 | 1 | 1283 | 6.5×10⁻⁹ | FAIL-Weyl, 697 | 5.6×10⁻⁵ → 20.5 / 20.5 | 1.8×10⁻⁴ → 30.1 | 0.51 |
| 0.9 | 100 | 1 | 75 | 1.4×10⁻⁶ | FAIL-Weyl, 1.7×10³ | 0.032 → 28.4 / 4.3 | 0.10 → 41.6 | 4.2 |
| 0.9 | 1000 | 1 | 746 | 1.4×10⁻⁹ | FAIL-Weyl, 8.8×10³ | 1.8×10⁻⁴ → 50.4 / 7.7 | 5.6×10⁻⁴ → 74.0 | 3.3 |
| all | all | 0.1 | | | FAIL-no-radiation-era (recollapse) at every b scanned | none | none | |

In the "post" column the cutoff ignores the χ mass before the crossing; for φ\* = 0.9 that is m = Gφ\* at φ = 0.

### 6.2 5D grid (NUMERICAL; `B1_RESULTS.json` → `runs`, `convergence`, `aggregate`)

Registered class: **every one of the 32 runs is INCONCLUSIVE (no plateau within the reliable evolution), at both spacings.** Per the registration, no cell is reliably classified, and the aggregate is **INCONCLUSIVE at λ_c = 1**.

The descriptive values below are **not classifications**.
- "End" values are at the last record, which is also the last reliable record for every b_cut run.
- dz = 10⁻³ and dz = 5×10⁻⁴ agree to ≤ 0.1% in r at b_cut and to ≤ 6% at the run end for b_need.

| Cell | b | max R_j | Reliable to H₀τ | r at end of reliable range | Ω_r there | Run end H₀τ | r at run end | Ω_r at run end |
|---|---|---|---|---|---|---|---|---|
| φ\*=0.5 G=100 y=1, b_cut | 6.5×10⁻⁶ | 5.7×10⁻⁴ | 7.44 (whole run) | **94.0** | 0.012 | 7.44 | 94.0 | 0.012 |
| φ\*=0.5 G=1000 y=1, b_cut | 6.5×10⁻⁹ | 1.6×10⁻⁴ | 7.35 | **674** | 0.0016 | 7.35 | 674 | 0.0016 |
| φ\*=0.9 G=100 y=1, b_cut | 1.4×10⁻⁶ | 3.6×10⁻⁴ | 7.35 | **3.56×10³** | 3×10⁻⁴ | 7.35 | 3.56×10³ | 3×10⁻⁴ |
| φ\*=0.9 G=1000 y=1, b_cut | 1.4×10⁻⁹ | 5.8×10⁻⁵ | 7.35 | **8.47×10³** | 1.3×10⁻⁴ | 7.35 | 8.47×10³ | 1.3×10⁻⁴ |
| y = 0.1, b_cut (four cells) | as above | ≤ 1.7×10⁻³ | 7.35–7.44 | 306 to 1.1×10⁵ (snapshots, not lower bounds: χ not yet decayed, Ω_χ ≤ 0.02 but 0.6–44× the radiation energy) | ≤ 0.002 | same | same | same |
| φ\*=0.5 G=100 y=1, b_need | 0.010 | 0.35 | 2.8–3.1 | 4.4–11 | 0.17–0.36 | 5.19 | 0.17 | 0.88 |
| φ\*=0.5 G=1000 y=1, b_need | 5.7×10⁻⁵ | 0.47 | 1.67 | −128 | 0.003 | 4.58 | 0.91 | |
| φ\*=0.9 G=100 y=1, b_need | 0.032 | 1.85 | 2.46 | −9 | 0.03 | 5.54 | −0.09 | |
| φ\*=0.9 G=1000 y=1, b_need | 1.8×10⁻⁴ | 0.60–0.67 | 2.04 | −20 to −39 | 0.003 | 4.46 | 0.058 | |
| y = 0.1, b_need (four cells) | as above | 0.35–1.9 | 1.7–2.6 | not meaningful | < 10⁻³ | 4.6–5.6 | −5 to 1.0 | |

- **b_cut cells (test-field regime).**
  - The constraint residuals stay ≤ 1.2×10⁻⁴ throughout on the fine grid; on the coarse grid the maximum is 4.2×10⁻⁴ (audit).
  - For three of the four y = 1 cells χ has decayed (Ω_χ ≤ 1.3×10⁻⁷) and Ra⁴ is flat to about 1% over the post-decay range, but that range is only Δln a = 0.24–0.37 long, shorter than the registered 0.5. The longest window meeting the A1 tolerance in any run is Δln a = 0.394.
  - **Exception (audit):** φ\* = 0.9, G = 100, y = 1 is not fully decayed at the end. Ω_χ = 3×10⁻⁵ is about 0.095 of Ω_r, and Ra⁴ rises 16% over the last Δln a ≈ 0.12. Its r = 3.56×10³ is also about 1.9× high from the ramp under-production (§6.4).
  - These are run-end values read descriptively, not final or plateau values. The conclusion r ≫ 0.1 does not depend on that.
  - **y = 0.1 (audit):** the values are snapshots, not lower bounds. Under a linear "all χ → R" estimate r would fall to about 52 (φ\* = 0.5, G = 100), 194, 2.5×10³ and 4.4×10³; massive χ redshifts more slowly than radiation, so the true values could be lower still. All remain ≫ 0.1.
- **b_need cells (strong backreaction).**
  - The near-shell constraint residuals exceed 0.05 within H₀τ ≈ 0.2–1.6 after the crossing. Every run then stops non-finite (H₀τ = 4.4–5.6) at the same coordinate time (to 0.01) at both spacings.
  - The run-end values are therefore outside the registered reliability window. They are listed only to show that they straddle 0.1 and do not form a plateau.

### 6.3 Internal consistency (numerical; per run in `B1_RESULTS.json`)

| Check (registration K3) | Requirement | Result |
|---|---|---|
| Shell-matter ledger, total, median relative residual | ≤ 10⁻³ | b_cut: 3.6–5.1×10⁻⁵. b_need: 1.8×10⁻⁴–1.3×10⁻³ for y = 1. **Exceeded** for φ\* = 0.9, G = 100 (2.4–2.9×10⁻³ for y = 1; 1.7–2.1×10⁻² for y = 0.1) and φ\* = 0.5, G = 1000 (1.1–1.3×10⁻³) |
| Wrong-sign control (j → −j), at b_need | ≥ 10× worse | **PASS**: 0.10–0.71, i.e. 20–2700× worse |
| Weyl identity with matter, median | ≤ 10⁻³ | **PASS**: 1.0–4.5×10⁻⁴ |
| Control 4H → 3H | ≥ 10× worse | **PASS**: 0.95–1.8×10⁻², i.e. 24–160× worse |
| Control "j dropped from w", at b_need | ≥ 10× worse | Passes for φ\* = 0.5, G = 100 (0.029–0.044, 180–290×). **Fails** for the G = 1000 and φ\* = 0.9 cells (about 1×). There j_prod acts only briefly, so dropping it barely moves the median. The control does not discriminate there; this is not evidence against the identity |

### 6.4 Reduced estimate against 5D (`compare_reduced_5d.py` → `B1_REDUCED_VS_5D.json`; numerical cross-check)

The reduced ODE is rerun with each 5D run's own production event, on the archived Y = 0 background.

- **b_cut, y = 1.**
  - φ\* = 0.5 (G = 100 and 1000) and φ\* = 0.9, G = 1000: R_5D/R_red = 1.03–1.04 and r_5D/r_red = 0.966 at the end.
  - φ\* = 0.9, G = 100: R_5D/R_red = 0.54. The field-space ramp needs |φ − φ\*| = 3|v\*|/√q ≈ 0.28, but φ_b only reaches 1.036, so production stops at about half its target (see Limitations).
- **b_cut, y = 0.1.** Agreement to ≤ 0.6% in r, except the same φ\* = 0.9, G = 100 cell (factor 1.9).
- **b_need.** R_5D/R_red = 0.04–0.21 at the end (≤ 0.32 at any time). The scalar backreaction (R_j = 0.35–1.9) cuts the produced radiation by a factor of 5 to 25, and r is not reproduced: the 5D end value is 0.9–12 times the reduced value, or has the opposite sign. **The reduced b_need and λ_c-needed values in 6.1 are therefore not reliable as quantitative thresholds.** They mark where backreaction becomes O(1).

### 6.5 Production-model controls (coarse spacing, φ\* = 0.5, G = 100, y = 1, b_need; registration step 4)

**Corrected after the audit.** The original table paired each control's r at its *reliable end* (−338, −1149, 22.5) with the baseline's r at the *run end* (0.17). The corrected values below are from `audit/AUDIT_TIMESERIES.json`. `B1_RESULTS.json` → `production_model_controls` → `last_r` still holds the mis-paired values; that JSON was not edited.

| Control | Reliable end H₀τ | r at reliable end | Run end H₀τ | r at run end | Change vs baseline at run end |
|---|---|---|---|---|---|
| ramp 1.5/√q | 1.85 | −338 | 4.76 | 0.316 | +86% |
| instant insertion (unpaid energy) | 1.57 | −1149 | 4.55 | 0.082 | −52% |
| D switched off | 2.62 | 22.5 | 5.22 | 0.169 | +0.3% (no effect) |
| baseline | | | 5.19 | 0.170 | |

No run has a plateau, so the registered "r changes by more than 20%" test cannot be evaluated. All run-end values lie beyond the runs' reliable ends. Read descriptively, the end state at b_need is sensitive to the ramp and to instant insertion at the level of about 2×, not orders of magnitude; the O(1/q) D correction has no effect. This sensitivity is still a sign that this regime is at the edge of the instant-preheating production model.

### 6.6 Why the runs stop (`diag_blowup.py` → `diag/*_blowup.json`; informative)

- **b_need** (φ\* = 0.5, G = 100, y = 1, dz = 10⁻³):
  - After production, a bulk region next to the shell (Z ≈ −0.4) crosses φ = −1 at T ≈ 1.85. φ = −1 is the critical point of W, beyond which the bulk potential is unbounded below.
  - φ_min falls to −1.33 by T = 2.85 while that region moves inward. The run stops at T = 2.89.
- **Y = 0 A1-equivalent run:** φ_min stays ≥ −1.007 until T = 3.0. The stop at T = 3.04 comes with growth of the shell lapse (B_b ≈ 4.4), consistent with the bounded chart's known limitation for a slowly expanding end state (A1 registration, second note).
- The b_need stop occurs at the same coordinate time at both spacings (6.2), so it is not a resolution artefact. Whether the b_need bulk excursion is a physical singularity or a gauge/chart effect is not determined here.

### 6.7 Exploratory cell (NOT registered; third dated note; `explore_analyze.py` → `B1_EXPLORE.json`)

Cell: φ\* = 0.5, G = 100, y = 1, b = 10⁻³ (λ_c = 5.36).

- Both spacings are reliable over the whole run: max R_j = 0.083, near-shell constraint residuals ≤ 1.9×10⁻³.
- The runs end non-finite at H₀τ = 5.78.
- From H₀τ ≈ 4 on, 𝒲a⁴ = 7.653 (constant to 10⁻⁴) and Ra⁴ = 11.9 (constant to 0.5%), so **r = 0.282** (the spacings agree to 4×10⁻⁴). Ω_r = 0.82 and Ω_vac = −0.05.
- This is not an A1 plateau: the steady interval is Δln a ≈ 0.18 < 0.5, and Ω_r < 0.9. Formally it is INCONCLUSIVE; read descriptively, it lies above the conservative bound.
- Compared with the b_cut run of the same (φ\*, G, y), r falls from 94 to 0.28 as b rises by 154×, roughly r ∝ b^−1.15.
- Production also lowers the Weyl term itself: 𝒲a⁴ = 7.65 here, against about 16 in the Y = 0 run.

## 7. Verdict and interpretation

| Statement | Label |
|---|---|
| Closed shell system (Ward identity, decay ledger, energy-paid production, junctions, Weyl balance) | **exact-verified** (sympy, with controls) |
| Calibration: the χ code path reproduces A1 Y = 0 exactly, and the friction-equivalent toy reproduces A1 Y = 1 (r = 0.0738) | **numerical** |
| Registered aggregate at λ_c = 1 | **inconclusive** (no cell reached the A1 plateau; nothing reliably classified) |
| At λ_c = 1 the roll-off's dark radiation dominates the χ-produced radiation: run-end r = 94 to 8.5×10³ (y = 1), and the radiation share is ≤ 1.2% | **numerical (descriptive)** (5D b_cut runs, two spacings, constraints ≤ 1.2×10⁻⁴ fine grid / 4.2×10⁻⁴ coarse); run-end values, not registered plateau values; φ\* = 0.9, G = 100 not fully decayed; y = 0.1 values are snapshots, not bounds |
| Reduced estimate: FAIL-Weyl at b_cut for all y = 1 cells, recollapse for all y = 0.1 cells | **conditional** (fixed background, validated against 5D to 3–4% in r where the ramp completes) |
| Can the derived χ channel replace the friction closure with χ below M₅? | **negative (conditional/descriptive)**: not within the scanned (φ\*, G, y); the registered outcome is INCONCLUSIVE |
| λ_c needed for r ≤ 0.1 (or 0.03) | **inconclusive**: λ_c ≳ 5 is required (exploratory r = 0.28 at λ_c = 5.4). The reduced values 11–50 (17–74) lie in a regime with O(1) backreaction, where the 5D runs are unreliable and the reduced estimate fails |

What this changes: the A1 PASS depends on the friction closure. The derived particle source does not reproduce it with χ masses below the 5D gravity scale. Near the cutoff the needed energy transfer is 10²–10⁴ times larger than the χ channel supplies.

What it does not change (as stated in the registration):
- the fine tuning d = d\* (the cosmological-constant problem restated);
- the restriction to δ = 0.1;
- the transient radiation era, because the residual vacuum is negative;
- the untreated thermalisation, rescattering, Bose enhancement and Coleman–Weinberg matching.

## 8. Limitations

- **No plateau anywhere.** The registered plateau needs Δln a = 0.5. The solver's bounded chart ends the evolution at H₀τ ≈ 7.4 (Δln a ≈ 0.25–0.4 after decay) even in the clean b_cut runs. A chart or solver that reaches later times is the main missing piece; the A1 Y = 1 run reached H₀τ ≈ 29 only because radiation dominated.
- **Strong-backreaction regime unreliable.** At b_need, production drains O(1) of the scalar's energy, constraint residuals exceed 0.05, and the bulk scalar crosses φ = −1 next to the shell. The b_need results, and the reduced λ_c values derived there, are not trustworthy as thresholds.
- **Production model.**
  - Instant preheating with the O(1/q) correction is computed for a prescribed trajectory. At b_need the trajectory is changed by the production itself.
  - The C² field-space ramp is a method choice. When the end point lies closer to φ\* than the ramp length (φ\* = 0.9, G = 100), production stops at about 50%.
  - The production-model controls show about 2× sensitivity of the run-end r at b_need to the ramp and to instant insertion, and none to D (corrected after the audit; see §6.5).
  - The cutoff test uses Λ̂ = max(m_χ, √q) over the whole run. The "post" variant (masses after the crossing) relaxes the φ\* = 0.9 requirements by 4–7×.
- **Seed.** Only one seed, dc = 10⁻², is used. A1 found seed differences of 5×10⁻⁴ in r.
- **Physics not included:** thermalisation of the decay products, rescattering, Bose enhancement at repeated crossings (only one crossing occurred in every run), and Coleman–Weinberg terms.
- **Control that does not discriminate.** The "j dropped" Weyl-identity control does not discriminate in 6 of 8 b_need cells (6.3), and the total-ledger tolerance is exceeded in the φ\* = 0.9, G = 100 and φ\* = 0.5, G = 1000 b_need cells.
- **N_eff thresholds** are the snippet-level values carried over from A1. They were not re-verified, since paper sites are blocked.

## 9. Files

**Scripts and solver:**
- `shell_matter.py`: matter sectors (friction / toy / χ).
- `evolve_b1.py`: 5D solver.
- `static_w.py`: copy of the A1 module.
- `b1_symbolic.py`: symbolic checks.
- `b1_reduced.py`: reduced estimate.
- `b1_analyze.py`: registered rule, calibration and ledgers.
- `compare_reduced_5d.py`, `diag_blowup.py`, `explore_analyze.py`: cross-checks, diagnostic and exploratory analysis.

**Outputs:**
- `B1_SYMBOLIC_CHECKS.json`, `B1_REDUCED.json`, `B1_RESULTS.json`, `B1_REDUCED_VS_5D.json`, `B1_EXPLORE.json`.
- `diag/*_blowup.json`.
- `b1_analyze.out`: printed table.

**Runs:**
- `runs/cal`: K1, K2.
- `runs/pre`: x_c pre-runs.
- `runs/main`: 32 grid runs.
- `runs/ctl`: production-model controls.
- `runs/explore`: exploratory cell.
- Each run has `*_summary.json` (parameters, code hashes, production events) and `*_timeseries.npz`.

**Kept, not used for results:** `dev/`, `superseded_time_ramp/`, `a1_copy/`.

**Run tolerances:**
- RK4 with Δt = 0.5 min ΔZ;
- reduced ODE with LSODA, rtol 10⁻⁹;
- 24 momentum nodes per cohort;
- convergence evidence from two spacings for every 5D run.

## 10. Reproduction (from this folder; 1 core, about 75 min for the 5D part)

```bash
export OMP_NUM_THREADS=1
python3 b1_symbolic.py                      # -> B1_SYMBOLIC_CHECKS.json
python3 b1_reduced.py                       # -> B1_REDUCED.json (fixed archived background)
./queue.sh jobs_cal.txt                     # K1, K2 -> runs/cal
./queue.sh jobs_grid.txt                    # pre-runs, 32 grid runs, 3 controls -> runs/{pre,main,ctl}
./queue.sh jobs_explore.txt                 # exploratory cell -> runs/explore
python3 -W ignore b1_analyze.py > b1_analyze.out   # -> B1_RESULTS.json
python3 -W ignore compare_reduced_5d.py     # -> B1_REDUCED_VS_5D.json
python3 -W ignore explore_analyze.py        # -> B1_EXPLORE.json
./chain_after.sh                            # (re)runs the explore queue and the two blow-up diagnostics -> diag/
```

`queue.sh` skips jobs whose summary already exists; delete a summary to rerun that job.
