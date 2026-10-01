# B1 pre-registration: the derived χ source in place of the friction closure (tuned model, δ = 0.1)

Registered 30 September 2026, before the reduced estimate, the calibration runs and any grid run. This file is not edited after registration; changes are appended below as dated notes with reasons.
Ricardo Maldonado's HDBLAST program. Status of everything that follows: numerical or conditional, internal, not peer reviewed.

## 0. Disclosure: what was done or seen before this registration

- **Restart.** A first attempt at this workstream was interrupted by a usage limit. It left copied code (`a1_copy/`, byte-identical to the A1 files, checked by `diff`), a draft solver `evolve_b1.py` and matter module `shell_matter.py`, and one development run `dev/smoke_chi` (G = 100, φ\* = 0.05, y = 1, b = 10⁻³, coarse grid dz = 2×10⁻³, L = 12, T = 3; its end state was seen: 𝒲/H₀² = 3.8×10⁻³, rad/H₀² = 9.9×10⁻⁵ at H₀τ = 7.08, i.e. r ≈ 39). No REGISTRATION.md existed.
- **This attempt, before registering:**
  - Reviewed and kept the first-attempt code; changed the production step from instantaneous insertion to the energy-consistent ramp of §1.3 (the instantaneous insertion creates χ energy that no junction pays for).
  - `b1_symbolic.py` → `B1_SYMBOLIC_CHECKS.json`: 20 exact identities and 8 wrong-sign or wrong-factor controls, all as required.
  - One development run `dev/dev2_ramp_G100_ps0.5_b1e-3` (φ\* = 0.5, G = 100, y = 1, b = 10⁻³, coarse grid, L = 12, x_c = 1.9, T = 3.5). Seen: the crossing at H₀τ = 1.487 (q = 128.3, D/q = −0.020), an unpaid energy of 1.6×10⁻⁹ H₀, and at H₀τ = 5.32 (end, no plateau yet) 𝒲/H₀² = 1.87×10⁻³ and rad/H₀² = 6.4×10⁻³ (r ≈ 0.29). This run is not used for any verdict.
  - The archived A1 tuned time series were inspected (read-only): the φ\* crossings of the Y = 0 run (φ\* = 0.5 at H₀τ = 1.486, v = 1.28, h = 0.645; φ\* = 0.9 at H₀τ = 1.819, v = 0.80, h = 0.285) and 𝒲a⁴ ≈ 16 at its end (ln a ≈ 2.1, residual vacuum −3.4×10⁻⁴ H₀²).
  - An order-of-magnitude estimate from these numbers suggested that at the λ_c = 1 cutoff the χ energy is roughly 10³ times too small for r ≤ 0.1. This estimate motivated including b_need (§3) in the grid.

## 1. Model (MODEL CHANGE B1, labelled)

### 1.1 Background model (unchanged from A1)

- 5D Einstein–scalar, κ₅² = 1 in the bulk: W = 1 − φ + φ³/3, U = ½W_φ² − (2/3)W², Z2-doubled bulk, one shell.
- **A1 model change, kept:** σ = 2W + δ(1 + cφ + dφ²/2), δ = 0.1, c = 0.5975949350280132, d = d\*(0.1) = −3.106933495673783 (`M8_QUADRATIC_TENSION_TUNING.json`).
- Seed dc = 10⁻² only (the pilot/A1 seed). A1 found the two seeds agree to 5×10⁻⁴ in r, so one seed is used to save cost. This is a registered reduction.

### 1.2 Shell matter: the B1 model change

The phenomenological friction closure κ₅²j = Yv is **removed** (Y = 0). The shell carries the 22 Sept matter-extension χ field instead. Its mass is m_χ² = μ₀² + ḡ²(φ_b − φ\*)², with μ₀ = 0. All quantities below are in H₀ units: G = ḡ/H₀, s = H₀τ, b = κ₅²H₀³ = (H₀/M₅)³, and κ₅²ρ/H₀ = b ρ̂.

- **χ gas.** The gas is represented by 24 momentum nodes per produced cohort (generalised Gauss–Laguerre in k²).
  - Number and momentum per node: n_i = N_i a⁻³ and p_i = k_i/a, with ω_i = (p_i² + m_χ²)^{1/2}.
  - Density, pressure and scalar source: ρ_χ = Σn_iω_i, p_χ = Σn_ip_i²/(3ω_i), and j = −∂L_m/∂φ = ḡ²(φ_b − φ\*)Σn_i/ω_i.
  - **Ward identity:** ρ̇_χ + 3H(ρ_χ + p_χ) = (j + j_prod)v − Q.
- **Decay.** χ → ψψ̄ into massless fermions, with rest-frame rate Γ = y²m_χ/(8π), time-dilated per mode: Ṅ_i = −Γ(m/ω_i)N_i.
  - Energy transfer: Q = Γm_χΣn_i.
  - Radiation ledger: Ṙ + 4HR = Q.
  - R is the decay radiation. Its thermalisation is not addressed.
- **Junctions (all three, exact, from the matter extension §2):**
  - n·A = (σ + κ₅²ρ)/6
  - n·B = (σ − κ₅²(2ρ + 3p))/6
  - n·φ = −(σ′ + κ₅²(j + j_prod))/2
  - Here ρ = ρ_χ + R and p = p_χ + R/3.

### 1.3 Production at each crossing of φ\*

- **Number.** Instant preheating, n_k = (1 + D/q)exp(−πk²/q), with q = G|dφ_b/ds| at the crossing.
- **Correction D.** D is the O(1/q) curved/expanding correction in the closed form of `HDBLAST_CHECKPOINT_20260927/preheating` (README §5, `CORRECTION_FIT.json`):
  - It is evaluated from a quartic fit to the preceding 0.2 H₀⁻¹ of the shell history.
  - It is applied when |D/q| ≤ 0.3. Otherwise it is capped at ±0.3 and the run is flagged "production formula outside validity".
- **Energy-consistent insertion (method choice, B1).**
  - The cohort is filled with a C² ramp lasting 3/√q in H₀τ, starting at the crossing.
  - The production power P = Σ Ṅ_i|_prod ω_i a⁻³ is paid by the scalar through j_prod = P/v in the scalar junction.
  - As a result, the Ward identity and the shell budget d(λ + ρ)/dτ + 3H(ρ + p) = −2(n·φ)v/κ₅² hold at all times (`B1_SYMBOLIC_CHECKS.json` C4, C5b).
- **Assumptions (stated, not tested):**
  - adiabatic vacuum before the crossing;
  - no Bose enhancement at repeated crossings (flagged if more than one crossing occurs);
  - no rescattering;
  - Coleman–Weinberg terms absorbed into σ (27 Sept R8 notes that this requires tuning);
  - perturbative decay.

### 1.4 Effective-theory cutoff

- The cutoff is Λ̂ = max(m̂_χ over the evolution from τ = 0, √q) ("full").
- b_cut(λ_c) = (λ_c/Λ̂)³, with **λ_c = 1 primary**.
- For the pre-set b values, Λ̂ is taken from the archived Y = 0 trajectory (computed by `b1_reduced.py`). It is re-checked after each run from the realised trajectory. If the realised Λ̂ exceeds the planned one by more than 2%, the run's actual λ_c = Λ̂_real b^{1/3} is reported.
- The "post" variant (masses after the crossing only) is reported alongside.

## 2. Quantities and the rule (A1 rule, verbatim, with R the decay radiation)

- **Shell identity:** H² = (σ + κ₅²ρ)²/36 + v²/12 − (σ′ + κ₅²(j + j_prod))²/48 + U/6 + 𝒲.
- **Radiation terms:** rad = σR/18 + R²/36, where R = κ₅²ρ_R is the decay radiation only. χ is **not** counted as radiation.
- **Weyl ratio and radiation share:** r = 𝒲/rad and Ω_r = rad/H².
- **Plateau (A1):** over the last Δln a = 0.5, H > 0 and both 𝒲a⁴ and Ra⁴ vary by at most 5% of their value. While χ is still decaying, Ra⁴ grows, so a plateau requires the decay to be essentially complete.
- **Per-run classes (A1, unchanged):**
  - PASS-conservative: plateau with Ω_r ≥ 0.9 and |r| ≤ 0.1. PASS-combined additionally requires |r| ≤ 0.03.
  - FAIL-Weyl: plateau with |r| > 0.1.
  - FAIL-no-radiation-era: either recollapse (H/H₀ < −0.05) before any time with Ω_r ≥ 0.9 and |r| ≤ 0.1, or a plateau with Ω_r < 0.5.
  - INCONCLUSIVE: anything else within the reliable evolution.
- **Reliability (A1, unchanged):**
  - The two spacings (dz_fine = 10⁻³ and 5×10⁻⁴, coarse spacing refined by the same factor) give the same class, and r agrees to 20% or 0.02 (for a recollapse, the H = 0 time agrees to 0.05 H₀⁻¹).
  - The near-shell relative constraint residuals are below 0.05 up to the classification time on the finer grid.
- **Run settings (A1 main batch):** L = 16, T_final = 14, κ = 10, projection window 15.9, bounded chart, x_c from a coarse old-chart pre-run with the same matter (A1 rule), wall limit 25 min.

## 3. Steps and grid

1. **Reduced estimate (CONDITIONAL, `b1_reduced.py` → `B1_REDUCED.json`).**
   - The shell ODE of §1.2–1.3 is solved on the archived A1 tuned Y = 0 trajectory, primary `main_dstar_Y0_dc1e-2_dzf5e-4`, with the dzf = 10⁻³ and dc = 10⁻⁴ archives as controls.
   - The background is fixed: no scalar backreaction and the archived 𝒲.
   - Beyond the archive, φ is frozen and 𝒲a⁴ and vac are frozen. H² = vac + 𝒲 + σ_f(R + X)/18 + (R + X)²/36 (variant "lowE"). A variant "fixed" omits the matter terms from H² and is used for r only.
   - The A1 rule is applied to the reduced series on a b grid 10⁻¹⁰ … 10^{−0.5} (39 log-spaced values).
   - **b_need** is the smallest grid b giving PASS-conservative in the lowE variant on the primary archive, rounded up to 2 significant digits. If there is none, b_need = 10^{−0.5} and it is labelled "maximum scanned".
   - Also reported: λ_c needed = Λ̂ b_need^{1/3} (full and post), the same for r ≤ 0.03, and R_j = b|ĵ|/|r_bσ′| along the archive.
2. **Calibration (must pass before grid runs are trusted).**
   - **K1:** chi code path switched off (`--source chi --G 0`) with the parameters and x_c of A1 `main_dstar_Y0_dc1e-2_dzf1e-3` (x_c = 1.9). Requirement: max |Δφ_b|, |ΔH/H₀| and |Δ𝒲/H₀²| against the A1 time series ≤ 10⁻⁸ over the common range, and the same stop.
   - **K2:** friction-equivalent toy source (`--source toy --Y 1 --Gt 10`: a radiation-like fluid fed by j = Yv that decays into R, so the total obeys the A1 ledger) with the parameters and x_c of A1 `main_dstar_Y1_dc1e-2_dzf1e-3` (x_c = 3.1).
     - Formal requirement (A1 reliability tolerances): the same class as A1 (PASS-conservative), with r within 20% or 0.02.
     - Stronger expectation, reported: max |Δφ_b|, |ΔH/H₀| ≤ 10⁻⁶.
   - **K3 (ledgers on every χ run, from the records):**
     - The shell-matter ledger d(X + R)/dτ + 3H(X + R + P) − κ₅²(j + j_prod)v = 0 has median relative residual ≤ 10⁻³. The wrong-sign control (j → −j) must be ≥ 10× worse where the source matters (b_need runs).
     - The Weyl identity with matter, Ẇ + 4H𝒲 = (1/6)[4k_s w v − 4Hv² − v v̇ + w ẇ − U′v], has median relative residual ≤ 10⁻³. The controls (4H → 3H, and j dropped from w) must be ≥ 10× worse; the j-dropped control applies to the b_need runs only.
3. **5D grid (MODEL CHANGE B1).**
   - Axes: φ\* ∈ {0.5, 0.9} × G ∈ {100, 1000} × y ∈ {0.1, 1} × b ∈ {b_cut(λ_c = 1), b_need}, where b_need is set per (φ\*, G, y) by step 1.
   - Every cell runs at both spacings.
   - Queue order (time-boxed at about 3.5 h in total): (i) y = 1 at b_need; (ii) y = 1 at b_cut; (iii) y = 0.1 at b_need; (iv) y = 0.1 at b_cut.
   - A cell not completed at both spacings within the time box is reported as NOT RUN or UNRELIABLE (single resolution). It is never classified.
4. **Production-model controls (coarse spacing only, informative).** For φ\* = 0.5, G = 100, y = 1 at b_need:
   - ramp length 1.5/√q instead of 3/√q;
   - instant insertion (unpaid energy);
   - D switched off.
   - A change of r by more than 20% flags the verdict as sensitive to the production model.

## 4. Verdicts and what they mean

- **B1-PASS at λ_c = 1:** at least one (φ\*, G, y) at b_cut(λ_c = 1) gives a reliable PASS-conservative.
  - Meaning: within this model at δ = 0.1, the derived χ channel can replace the friction closure with χ masses below M₅. PASS-combined is reported separately.
- **B1-FAIL at λ_c = 1:** every completed b_cut cell gives a reliable FAIL (either kind), and the reduced estimate gives r > 0.1 at b_cut for all eight (φ\*, G, y).
  - Meaning: with χ masses below M₅ the derived channel cannot replace the friction closure for the scanned parameters. The dark radiation made by the roll-off itself dominates.
- **INCONCLUSIVE at λ_c = 1:** otherwise.
- **Above the cutoff (conditional statement):** the b_need cells report whether a reliable PASS exists at all, and the λ_c it needs.
  - A PASS there with λ_c_needed > 1 means the channel works only with χ masses above the 5D gravity scale, i.e. outside the stated effective theory.
  - A FAIL there as well means the channel fails even when the cutoff is relaxed to b_need.
- **Backreaction.** R_j (max over the run) is reported for every run. The shell-matter ledger and the Weyl identity are the internal consistency checks.
- **What any outcome does not change:**
  - the tuning d = d\* (the cosmological-constant problem restated, A1);
  - the δ = 0.1 restriction;
  - the transient character of the radiation era (negative residual vacuum);
  - the untreated thermalisation, rescattering, Bose enhancement and Coleman–Weinberg matching.

---

### Dated note, 30 September 2026 (after the reduced estimate `B1_REDUCED.json`, before any grid run)

1. **b values fixed by the registered rule (primary archive, lowE variant, rounded up to 2 significant digits).**
   - b_need(y = 1): 0.010 (φ\* = 0.5, G = 100), 5.7×10⁻⁵ (0.5, 1000), 0.032 (0.9, 100), 1.8×10⁻⁴ (0.9, 1000).
   - b_cut(λ_c = 1, full): 6.49×10⁻⁶ (0.5, 100), 6.49×10⁻⁹ (0.5, 1000), 1.39×10⁻⁶ (0.9, 100), 1.39×10⁻⁹ (0.9, 1000).
2. **Change for the y = 0.1 cells (reason given).** For y = 0.1 the reduced estimate finds no PASS at any scanned b. The reason is not too little energy: the χ decay (Γ = 0.02–0.2 H₀) completes only when H has dropped so far that the negative residual vacuum turns the expansion around before an Ra⁴ plateau can form.
   - The registered fallback b = 10^{−0.5} would put κ₅²ρ_χ far above the tension (b ρ̂ ≫ σ̂, R_j ≫ 1). That is a different, high-energy regime which the fixed-background estimate cannot represent and which says nothing about replacing the friction closure.
   - The y = 0.1 cells are therefore run at the y = 1 b_need of the same (φ\*, G), i.e. the b at which the produced energy suffices when the decay is prompt. The b = 10^{−0.5} runs are not made.
   - The pass/fail rule, thresholds, reliability criteria, queue order and all other settings are unchanged.

### Second dated note, 1 October 2026 (restart; before the grid batch with the final code)

1. **Restart disclosure.** The previous attempt (30 September) was interrupted by a usage limit after the first dated note. It left, after the registration: the K1/K2 calibration runs (made with an older solver build, sha 888cfc92…; now moved to `superseded_time_ramp/cal_oldcode/` and rerun with the final build), the first grid cell with the time ramp of §1.3 (`superseded_time_ramp/`), and development runs `dev/dev3_*`, `dev/dev4_*`, `dev/dev5_*`. All of these outcomes were seen before this note:
   - Time-ramp cell φ\* = 0.5, G = 100, y = 1, b = 0.010: both spacings stopped non-finite shortly after production (H₀τ = 2.86 at dz_f = 10⁻³, 2.24 at 5×10⁻⁴), no classification.
   - Field-space ramp, same cell, dz_f = 10⁻³, x_c = 2.2 (`dev3`, `dev4`): ran to H₀τ ≈ 5.2–5.3 and then stopped non-finite; at the end 𝒲/H₀² ≈ 1.7×10⁻³, rad/H₀² ≈ 9.9×10⁻³ (r ≈ 0.17, still decreasing), Ω_r ≈ 0.88, no plateau. The bulk scalar minimum near the shell (Z ≈ −0.5) drifts below −1 from H₀τ ≈ 2 onward before the stop.
   - Chi sector off (Y = 0) with x_c = 3.1 and 4.1 (`dev5`): stopped non-finite at H₀τ = 6.49 and 6.29 (the A1 Y = 0 run with x_c = 1.9 stops at 7.35).
2. **Method change (production insertion; model, thresholds, rule, reliability criteria and grid unchanged).** The C² ramp in *time* of §1.3 is replaced by the C² ramp in *field space* (cohort progress u advances with |φ_b − φ\*|/Δφ_r, Δφ_r = 3|v\*|/√q, i.e. the same 3/√q at the crossing velocity), with a C¹ stall switch that stops production when the scalar stops moving away from φ\* (`shell_matter.ChiGas` docstring). Reason: in the time ramp the force j_prod = P/v diverges when v → 0 while production is still running, which happens at b_need where the backreaction slows the scalar; the field-space ramp keeps j_prod bounded. Energy is still paid exactly (Ward identity, `B1_SYMBOLIC_CHECKS.json` C4c, C4d). The `--ramp` control (1.5 instead of 3) now refers to the field-space ramp.
3. **Expectation stated now (so it is not read as a later choice).** From item 1, the 5D runs may stop before the A1 plateau window (Δln a = 0.5 after decay) is reached. Such runs are classified INCONCLUSIVE by the registered rule; no rule is relaxed. Values of r and Ω_r at the last reliable time are reported as descriptive numbers only (not classifications).

### Third dated note, 1 October 2026 (during the grid batch; additions only, nothing registered is changed)

Seen so far (y = 1 cells): every b_need run stops non-finite at H₀τ ≈ 4.5–5.5 at both spacings, with near-shell constraint residuals above 0.05 after production (so UNRELIABLE by §2); the b_cut runs reach H₀τ ≈ 7.3 with Δln a ≈ 0.24 after the decay (no A1 plateau), r ≈ 94 at the end for φ\* = 0.5, G = 100. Added, **informative only, not used for any registered verdict**:
1. `compare_reduced_5d.py` → `B1_REDUCED_VS_5D.json`: the fixed-background reduced ODE recomputed with each 5D run's own production event and compared with the 5D R and r over the common range (validates the reduced estimate in the test-field regime; measures backreaction at b_need).
2. `diag_blowup.py` → `diag/`: dense-snapshot rerun of one coarse b_need cell and the A1-equivalent Y = 0 case to locate where the evolution becomes non-finite.
3. Exploratory 5D cell at an intermediate b (φ\* = 0.5, G = 100, y = 1, b = 10⁻³, i.e. λ_c ≈ 5.4), both spacings, to see the trend of r with b between b_cut and b_need. Reported as exploratory.
