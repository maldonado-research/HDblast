# B2: the A1 tuned-vacuum test at the registered detuning δ = 10⁻³

Ricardo Maldonado's HDBLAST program; prepared with AI assistance. Internal working note, not peer reviewed. Every number below comes
from a script in this folder that writes JSON (file named in each row). Status labels: exact-verified / certified / numerical /
conditional / negative / inconclusive.

## Post-audit corrections (1 October 2026)

This README was written before the independent audit. The sentences below were corrected in place to agree with `AUDIT.md` and
with the final labels in `../CLAIM_REVIEW.md`. The registered verdict (INCONCLUSIVE at δ = 10⁻³) is unchanged. No run was
repeated for these edits.

- **Breakdown mechanism restricted.** "Every Y = 1 run ends 1–2.5 H₀⁻¹ after the crossing" is false as a general statement. It
  holds for the dc = 10⁻² runs, the x_c = 18.7 sensitivity set and the S3′/C3′/E runs. The registered dc = 10⁻⁴ pair ends
  *before* F_inf, about 5–6 H₀⁻¹ before the crossing located by its κ = 10 pre-run (H₀τ ≈ 18.7–18.8), because the registered
  x_c = 11.8 is a κ = 0 artifact. C7 (κ = 0) ends at T = 6.79.
- **x_c is not ruled out.** Moving x_c from 11.3 to 12.3 moves the run end from H₀τ 12.47 to 13.22 (reliable end 11.44 → 12.70),
  a larger gain than any resolution doubling. Only the qualitative failure near F_inf is x_c-independent.
- **"Plain refinement cannot fix this / ~10⁶× cost"** is relabelled **conditional** (a linear extrapolation from three points
  of one seed, one x_c and one chart family). The fixed-T residual converges by about 10× per doubling.
- **Amendment 2 not formally accepted.** As written, criterion (ii) fails (`monotone: false`) and criterion (iii) fails (C2d
  misses the late estimator by −2.2×10⁻⁴). The earlier text reinterpreted both after the data were seen. The seed-offset
  explanation is supported only by post-hoc estimators, including an independent offset-free derivative estimator from the audit
  (1.657184–1.657187 on C2a8/C2b; 1.657048 on C2d).
- **C2 post-hoc values.** 1.65718–1.65723 holds on the ρ-grid runs C2a–C2c only. On C2d (tanh grid) the mode fit gives
  1.65705–1.65709, about 1.4×10⁻⁴ low. The late-estimator passes are marginal (spreads 1.8 and 1.9×10⁻⁴ against 2×10⁻⁴).
- **Five dated notes, not four.** S3′ was launched about a minute before note 4, and D6 before note 3 (both diagnostic only).
- **A1 stop at 11.0 (dz = 2×10⁻⁴):** only its loss of reliability (T ≈ 7) is clearly numerical; its end (H₀τ 10.97) is near the
  physical crossing (11.36).
- **Y = 0 baseline:** the S1 Y = 0 runs lose reliability at or before their crossing; "the solver works for a fast roll" holds
  at S2 only.
- **C5 is a weak check here:** the 4H → 3H control discriminates only at the median, and "R dropped" is indistinguishable.
- **Minor:** δ = 0.1 crossing H/H₀ is 0.20 (D7), not 0.26; max Ω_r ≤ 0.003 holds within the reliable part only (0.013 over all
  records, C3′); x_c = 18.7 for dc = 10⁻⁴ is single-resolution; quoted end H₀τ values are the run-summary values (e.g. 12.47),
  while `B2_RESULTS.json` `H0tau_end_run` is the last record (e.g. 12.37); C7's residual at the stop is O(1), not 10⁻².

See `AUDIT.md` and `audit/*.json` for the checks.

## Bottom line

- **Registered verdict at δ = 10⁻³: INCONCLUSIVE.** No Y = 1 run (two seeds × two spacings, plus every control and
  exploratory variant) reaches a plateau, a recollapse, or any other classifiable state. The dc = 10⁻² runs, the x_c = 18.7
  sensitivity set and the S3′/C3′/E runs end 1–2.5 H₀⁻¹ after the shell crosses the light cone of the static vertex, with
  φ_b ≤ 0.18 and Ω_r ≤ 0.003 within the reliable part (0.013 over all records, C3′). The registered dc = 10⁻⁴ pair ends *before*
  F_inf, about 5–6 H₀⁻¹ before its own crossing, because of the registered x_c = 11.8 artifact; C7 (κ = 0) ends at T = 6.79. The
  A1 δ = 0.1 PASS therefore still says nothing about the registered detuning.
- **Why (numerical diagnosis, restricted to the dc = 10⁻² runs and the sensitivity set; see Post-audit corrections):**
  - At δ = 10⁻³ with the Y = 1 friction closure, the shell reaches the vertex light cone **before it rolls**: the old-chart proper
    time freezes at H₀τ ≈ 11.36, with φ_b ≈ 0.094 and H/H₀ ≈ 1.002 (D6/D7). At δ = 0.1 the crossing came at the *end* of
    the roll (φ_b ≈ 1, H/H₀ ≈ 0.20 at F_inf, D7).
  - Past the crossing, the bounded chart's shell lapse grows at d ln(lapse)/dT ≈ 3.5–4 (δ = 0.1: 1.24), and a near-shell
    momentum-constraint violation (z ≈ −0.13) blows up within 0.3–0.9 chart units.
  - The far-region chart saturation is not the cause: s_max = 200 and the asinh chart end at the identical point. (These
    charts coincide with the default near the shell, so they do not test the near-shell gauge.) Neither is the time step
    (halving it moves the end by 0.001 in T), nor the finite-difference error of the static background in full-field mode
    (D5: ≤ 2×10⁻⁷).
  - **x_c is not ruled out:** moving x_c from 11.3 to 12.3 (C3′) moves the run end from H₀τ 12.47 to 13.22, more than a
    resolution doubling gains. What is robust is that every such run dies within about 0.3–0.9 chart units of F_inf.
  - Resolution helps only marginally. Each doubling of near-shell resolution adds about 0.65 H₀⁻¹ (run end H₀τ 12.47 → 13.10
    → 13.77 for S1 → S2 → S3′) and less to the reliable end (11.44 → 12.12 → 12.55), at 4× the cost per doubling.
- **Growth calibration C2: fails as registered, explained.** The registered windows give 1.65630 with spread 2.4×10⁻³ (limit
  2×10⁻⁴). The cause is a neutral-mode offset in the seed, whose effect on the window rates decays as e^{−λτ} (measured gap ratio
  0.42–0.45 per 0.5 H₀τ; predicted e^{−λ/2} = 0.437).
  - A supplementary late-window estimator passes on both registered spacings: 1.657086 and 1.657071. It was defined in dated
    note 1, after the registered C2a result had been seen but before C2a8 and C2b were run, and is labelled post-registration.
  - A three-term mode fit (post hoc diagnostic) gives λ = 1.65718–1.65723 on the ρ-grids (C2a–C2c), against the Chat 14
    target 1.65719 (D4). On the tanh-grid run C2d it gives 1.65705–1.65709, about 1.4×10⁻⁴ low.
  - By Amendment 2's own criteria as written, the seed-transient interpretation is **not formally accepted** ((ii) and (iii)
    fail). The explanation rests on post-hoc estimators; the audit's offset-free derivative estimator supports it independently.
  - The A1 C2 failure (1.698) was mostly this offset, made worse by the coarse A1 tables.
- **x_c rule fixed and pre-registered:** (a′) a κ = 10 pre-run. It gives x_c = 11.3 for dc = 10⁻², the same at two
  resolutions (11.2 and 11.3). The dc = 10⁻⁴ value (18.7) rests on one resolution. The undamped pre-runs are numerically unstable at δ = 10⁻³. That instability is why the A1
  pre-run returned x_c = ∞.

## Question

Does the A1 candidate radiation era (δ = 0.1, d = d\*, Y = 1: PASS-conservative, r = 0.073) also exist at the registered
detuning δ = 10⁻³ with d = d\*(10⁻³) = −3.1942416958680835? The task also asked:
- why A1 could not finish δ = 10⁻³;
- whether the x_c rule can be fixed;
- whether the growth calibration C2 can be passed.

## Registration (REGISTRATION.md; five dated notes appended; the original body cannot be verified as unedited without a snapshot)

- **Status of the file.** `REGISTRATION.md` was written on 30 Sept by the first, interrupted attempt of this round. The restart
  instructions said none existed; it does, so it was treated as binding.
  - **Note 1** (written before any further run) discloses everything seen since registration and adds:
    - amendment 1 (x_c fallback (a′): a κ = 10 pre-run);
    - amendment 2 (C2 runs to T_final = 8; a supplementary late-window estimator, labelled post-registration).
  - **Note 2** (after the pre-runs, before any main run) keeps the registered x_c for dc = 10⁻⁴ (11.8, which comes from an
    unstable κ = 0 pre-run). It adds a sensitivity set at x_c = 18.7 that can only downgrade the result.
  - **Notes 3, 4 and 5** add exploratory chart variants and one finer-resolution diagnostic. None enters the aggregate.
    D6 was launched before note 3, and S3′ about a minute before note 4 (both diagnostic only; corrected after the audit).
- **Model (MODEL CHANGE as in A1):**
  - σ = 2W + δ(1 + cφ + dφ²/2), with δ = 10⁻³, c = 0.5975949350280132 and d = d\*.
  - Shell radiation R with p = R/3, and friction closure κ₅²j = Yv (Y = 1; Y = 0 as baseline).
  - Seeds: static shell with c + dc, dc ∈ {10⁻², 10⁻⁴}.
- **Rule:** A1 §3, unchanged.
  - Plateau: Δln a = 0.5 window with 5% constancy of 𝒲a⁴ and Ra⁴.
  - PASS-conservative: Ω_r ≥ 0.9 and |r| ≤ 0.1.
  - FAIL-Weyl: |r| > 0.1. FAIL-no-radiation-era: recollapse, or Ω_r < 0.5 at the plateau.
  - INCONCLUSIVE: anything else, including runs that end first.
  - Reliability: the two spacings agree, and the near-shell residuals stay < 0.05.
  - Aggregate: PASS if both seeds give a reliable PASS; FAIL if both give a reliable FAIL; INCONCLUSIVE otherwise.

## Method

- **Solver.** `evolve_a1.py` is a copy of the A1 solver, with the first attempt's numerical-only changes:
  - static-reference tables at dx = 2.5×10⁻⁵ (A1: 2.5×10⁻⁴; the table junction error falls from 1.8×10⁻¹⁰ to 3×10⁻¹⁴, see D1);
  - a ρ-adapted grid (a fixed number of points per local wall width), checkpoint/restart, and caches;
  - this round added `--smax` (exploratory only).
- **Spacings.** S1: dz_f = 3×10⁻⁴ (42.7 points per wall 1/ρ_b). S2: dz_f = 1.5×10⁻⁴ (85.4 points per wall).
- **Run settings.** L = 24, κ = 10 (Remedy B), no projection, bounded chart, deviation form until T_switch = F_inf − 1, RK4.
- **x_c selection:** `select_xc.py` → `XC_SELECTION.json`.
- **Analysis:** `b2_analyze.py` → `B2_RESULTS.json`. It is the A1 classification, plateau, Weyl-identity and ledger code
  unchanged, plus the supplementary C2 estimators.
- **Diagnostics:** `diag/*.py` → `diag/D*.json`.

## Results

### Diagnosis of the A1 δ = 10⁻³ failure (D2, D7; numerical)

| A1 symptom | Cause found |
|---|---|
| Pre-run never met d b_b/dT < −0.9, so x_c = ∞ | The undamped (κ = 0) old-chart run is numerically unstable at δ = 10⁻³: the lapse *rises* from T ≈ 1 and the near-shell residual grows. This happens on the A1 tanh grid, on the B2 ρ-grid (non-finite at T = 6.41) and at S1 (\|B_b\| > 40 at T = 6.79). With κ = 10 the lapse decays smoothly and the criterion is met (x_c = 11.2 at 32 points per wall, 11.3 at 42.7; residual ≤ 4×10⁻⁶). For dc = 10⁻⁴, κ = 0 meets the criterion at resolution-dependent times (11.8 and 13.7) driven by a growing residual at z ≈ −1; κ = 10 gives 18.7. |
| Spacings stopped at H₀τ 3.1 (dz 1.5×10⁻⁴) and 11.0 (2×10⁻⁴) | Two different numerical losses of reliability (the 2×10⁻⁴ run's *end*, H₀τ 10.97, is close to the physical crossing at 11.36 and is not clearly numerical). At 1.5×10⁻⁴ the Newton projection (Remedy A) diverged and the run started from constraint-violating data. At 2×10⁻⁴ a residual at the abrupt tanh fine-to-coarse transition (z ≈ −0.05, ratio 20) grew until breakdown. |
| C2 drift (1.698 vs 1.65719) | (i) The A1 tables violated the φ junction by 1.8×10⁻¹⁰, 35× the dc = 10⁻⁸ seed mismatch (D1); the control rerun with A1 tables reproduces early window rates of 1.88 and 1.75. (ii) Both A1 and B2 C2 runs carry a seed offset (neutral mode) whose effect decays as e^{−λτ}; see C2 below. Not a window, seed-size or run-length problem alone, and not grid coarseness: 75 and 134 points per wall agree to 2×10⁻⁵. |

### Growth calibration C2 (registered model d = 0, Y = 0, dc = 10⁻⁸; target 1.65719 within 2×10⁻⁴)

| Run | Registered estimator (windows 2–5): rate, spread | Late windows 4.5–7 (supplementary): rate, spread | Mode fit λ (D4) |
|---|---|---|---|
| C2a (75 pts/wall, T_f 6.5) | 1.65630, 2.4×10⁻³ (fail) | 1.65711, 0.9×10⁻⁴ | 1.65718–1.65722 |
| C2a8 (same, T_f 8) | 1.65630, 2.4×10⁻³ | **1.657086, 1.8×10⁻⁴ (pass)** | 1.65718–1.65722 |
| C2b (134 pts/wall) | 1.65632, 2.4×10⁻³ (fail) | **1.657071, 1.9×10⁻⁴ (pass)** | 1.65719–1.65723 |
| C2c (dc 3×10⁻⁸) | 1.65691, 7.9×10⁻⁴ | 1.65695, 4.2×10⁻⁴ (nonlinear onset) | 1.65719–1.65725 |
| C2d (Chat 14 tanh grid, fine tables) | 1.65615, 2.5×10⁻³ | 1.65697, 1.8×10⁻⁴ (diff −2.2×10⁻⁴) | 1.65705–1.65709 |
| A1 configuration, A1 tables (control) | 1.67100, 3.8×10⁻² | n/a | 1.65743–1.65773 |

- **Registered C2: FAIL** (numerical). The registered window set includes the offset transient: window rates rise
  1.6437 → 1.6571 with successive gaps shrinking by 0.42–0.45 per 0.5 H₀τ. That is the e^{−λ/2} = 0.437 expected from a
  constant offset B in dev = Ae^{λτ} + B.
- **Supplementary late estimator: PASS** on C2a8 and C2b (post-registration estimator, defined after C2a was seen and before
  these runs).
- **Amendment-2 test (corrected after the audit):**
  - Criterion (i) holds, marginally (spreads 1.8 and 1.9×10⁻⁴ against the limit 2×10⁻⁴).
  - Criterion (ii) **fails as written**: `B2_RESULTS.json` reports `monotone: false`, because the late windows turn down by
    ≤ 1.8×10⁻⁴ (plausibly the nonlinear onset captured by the C term of the mode fit).
  - Criterion (iii) **fails as written**: C2d fails the late estimator (diff −2.2×10⁻⁴), although it shows the same gap ratios.
  - An earlier version of this README read (ii) as holding "up to the 5–6 window" and (iii) as holding on the gap ratios. That
    was a reinterpretation of the pre-stated criteria after the data were seen. By the amendment's own rule the seed-transient
    interpretation is **not formally accepted**. Practical effect: none, since no tuned run was classified.
- The seed-independence check C2a vs C2c on the registered estimator differs by 6×10⁻⁴ (expected 10⁻⁴), because the larger seed
  turns nonlinear inside the window. The mode fit agrees to 4×10⁻⁵.
- **Conclusion (numerical, post hoc):** the solver reproduces the Chat 14 growth rate to ≤ 4×10⁻⁵ on the ρ-grids once the seed
  offset is modelled; on the C2d tanh grid both the mode fit and the audit's derivative estimator are about 1.4×10⁻⁴ low. The registered C2 rule, as written, cannot pass with this seed. The downgrade clause of §5 is moot, because no tuned
  run reached a classification.

### Pre-runs and x_c (`XC_SELECTION.json`)

| (Y, dc) | Chain | x_c used |
|---|---|---|
| (1, 10⁻²) | registered ρ4e-4 κ0: non-finite T 6.41; (a) S1 κ0: \|B_b\| > 40 at T 6.79; **(a′) S1 κ10: 11.3** (H₀τ 11.26) | 11.3 |
| (1, 10⁻⁴) | **registered ρ4e-4 κ0: 11.8** (artifact, see above); for the record (a) 13.7, (a′) 18.7 | 11.8 (sensitivity set: 18.7) |
| (0, 10⁻²) | registered: 1.6 | 1.6 |
| (0, 10⁻⁴) | registered: 2.7 | 2.7 |

### Main runs, controls and exploratory runs (`B2_RESULTS.json`, `diag/D7_CHART_END.json`)

| Run | x_c | Ends at T (F_inf) | H₀τ end | Reliable to H₀τ | max φ_b | max Ω_r | Class |
|---|---|---|---|---|---|---|---|
| Y1 dc1e-2 S1 | 11.3 | 11.58 | 12.47 | 11.44 | 0.118 | 0.002 | INCONCLUSIVE |
| Y1 dc1e-2 S2 | 11.3 | 11.63 | 13.10 | 12.12 | 0.146 | 0.002 | INCONCLUSIVE |
| Y1 dc1e-4 S1 | 11.8 | 11.60 (before F_inf) | 12.89 | 12.09 | 0.001 | 0.000 | INCONCLUSIVE |
| Y1 dc1e-4 S2 | 11.8 | 11.70 (before F_inf) | 13.56 | 12.82 | 0.001 | 0.000 | INCONCLUSIVE |
| Y0 dc1e-2 S1 / S2 (baseline) | 1.6 | 3.20 / 2.93 | 7.81 / 7.85 | 2.46 / 7.49 | 1.000 | 0 | INCONCLUSIVE (no radiation by construction) |
| Y0 dc1e-4 S1 (baseline) | 2.7 | 4.25 | 8.87 | 3.41 | 1.000 | 0 | INCONCLUSIVE |
| Y0 dc1e-4 S2 (baseline) | 2.7 | 3.97 | 8.92 | 8.19 | 1.000 | 0 | INCONCLUSIVE |
| sens Y1 dc1e-4 S1 / S2 | 18.7 | 18.97 / 18.70 | 19.83 / 20.47 | 18.76 / 19.32 | 0.012 / 0.015 | 0.000 | INCONCLUSIVE |
| C3′ Y1 dc1e-2 S1 | 12.3 | 13.20 | 13.22 | 12.70 | 0.151 | 0.003 | INCONCLUSIVE (class unchanged) |
| C7 Y1 dc1e-2 S1 κ = 0 | 11.3 | 6.79 (deviation mode) | 8.29 | 6.93 | 0.029 | 0.000 | INCONCLUSIVE (κ = 0 unstable; near-shell residual O(1) at the stop) |
| E1 s_max = 200 (exploratory) | 11.3 | 11.580 | 12.4736 | 11.44 | 0.118 | 0.002 | identical to default |
| E3 asinh chart (exploratory) | 11.3 | 11.580 | 12.4736 | 11.44 | 0.118 | 0.002 | identical to default |
| E4 time step halved, cfl 0.25 (exploratory) | 11.3 | 11.579 | 12.45 | 11.40 | 0.118 | 0.002 | same end: time-step instability excluded |
| S3′ dz_f 7.5×10⁻⁵, 170 pts/wall (diagnostic) | 11.3 | 11.68 | 13.77 | 12.55 | 0.182 | 0.003 | INCONCLUSIVE |

- **Aggregate (registered): INCONCLUSIVE.** The class is the same at both spacings for every pair, but none is decisive.
- **Checks during the evolved part (numerical):**
  - C5 Weyl identity: median relative residual 7×10⁻⁸ to 1.1×10⁻⁶ (≤ 10⁻³, pass).
  - 4H → 3H control: 1.7–1.9×10⁻⁵, 16–250× worse (pass at the median only; at p95 the identity residual is ≥ the 3H control, so
    C5 is a weak check in this regime).
  - "R dropped" control: indistinguishable, because R is negligible in this phase. As in A1, it cannot discriminate.
  - C6 radiation ledger: median 5–8×10⁻⁷ (≤ 10⁻⁴, pass).
  - The Y = 0 baseline rolls to φ_b ≈ 1 (H/H₀ ≈ 0.14 at F_inf and 0.057 at the end) and runs 1.3–1.6 chart units past F_inf.
    However, both S1 Y = 0 runs lose reliability (residual > 0.05) at or before their crossing (H₀τ 2.46 and 3.41). Only at S2
    does the machinery work reliably for a fast roll at δ = 10⁻³.
- **Last reliable state of the Y = 1, dc = 10⁻² runs (H₀τ ≈ 11.4):** 𝒲/H₀² ≈ −7.9×10⁻⁴ and rad/H₀² ≈ 1.6×10⁻³. This is a
  transient, not a plateau, and is not classified.

### Why the post-crossing evolution fails, and what it would need (D5, D6, D7; numerical)

- **Background error is not the cause.** The full-field finite-difference error of the exact static background at T_switch is
  2×10⁻⁷ (S1) and 4×10⁻⁸ (S2). That is smaller than in the A1 δ = 0.1 runs (2×10⁻⁵ to 3×10⁻⁴), which ran fine.
- **The crossing comes before the roll.**
  - Old chart (D6): d B_b/dT → −0.98 and b_b + T → 11.36 (so x_c = 11.3 is consistent). The shell crosses the light cone at
    H₀τ ≈ 11.36 while still near de Sitter (φ_b 0.094, H/H₀ 1.002).
  - After the crossing, the chart lapse grows at 3.5–4 per unit T. The momentum-constraint violation at z ≈ −0.13 reaches 0.29
    within ≈ 0.2 chart units, and the run stops. At fixed T this residual converges by about 10× per doubling: a steepening
    feature that each resolution resolves a little longer, not an O(1) continuum blow-up at a fixed time.
- **What it would take to finish:**
  - The roll at Y = 1 is friction-limited, at 0.33 per H₀τ in log deviation (D2). Reaching φ_b ~ 1 from φ_b ≈ 0.1 needs
    ≳ 7 H₀τ more past the crossing, followed by a plateau window (Δln a = 0.5).
  - Each doubling of resolution buys about 0.65 H₀τ of run (0.4–0.7 of reliable run) at 4× the cost. S3′ took 44 CPU-min to
    reach H₀τ 13.8. A linear extrapolation says the ≳ 7 H₀τ still needed would take ~10 more doublings, about 10⁶ times the
    S3′ cost. **Conditional:** this extrapolates three points of one seed at one x_c on one chart family, so "plain refinement
    cannot fix this" is not established.
  - The likely requirement is a gauge adapted to the post-crossing region: a second chart stage, re-mapping the null coordinates
    with the measured post-crossing lapse, or a proper-time-type gauge for the shell. That is a solver development task, not a
    longer run.
  - The x_c rule itself is fixed (a′).

## Limitations

- **Registration timing.**
  - The registration and the x_c amendment were written by this program during the work: note 1 came after the first attempt's
    κ = 10 development pre-run (x_c = 11.2) had been seen.
  - The supplementary C2 estimator was defined after the registered C2a result was seen.
  - Both are labelled as such.
- **The dc = 10⁻⁴ main pair uses x_c = 11.8** from a κ = 0 pre-run that is a numerical artifact. Using it was required by the
  registered chain. The x_c = 18.7 sensitivity set fails the same way.
- **The friction closure is phenomenological.** At δ = 10⁻³, "Y = 1" is a 10× stronger friction in Hubble units than at δ = 0.1
  (Y ρ_b = 79), which is why the roll is slow. This was not rescaled (as registered).
- **Out of scope:** the derived χ source and other Y values were not run.
- **Exploratory runs** E1, E3 and E4 cover only S1. E2 was dropped after E1 matched the default to 5 digits.
- **Data files.** Time series (`*_timeseries.npz`) and checkpoints/shell caches (`*.pkl`, `shell_cache/`, about 1.1 GB,
  regenerable) are excluded by the checkpoint folder's `.gitignore`. The JSON summaries and result files carry every quoted
  number.

## Reuse of the first attempt (30 Sept)

Inspected and reused after the checks noted below:
- `evolve_a1.py`, `static_w.py`: copies with numerical-only changes. `a1_analyze.py` is identical to the A1 file (sha256
  checked).
- `b2_analyze.py`: extended here.
- `diag/D1`–`D3`, `TIMING_GRID`, the dev runs (`runs/dev`) and the 30 Sept pre-runs (`runs/pre/pre_Y0_*`, `pre_Y1_dc1e-2`).
- C2a (`runs/cal/C2a_rho1.7e-4`) and the diagnostic C2 control (`runs/diag`).
- `shell_cache/`.

The interrupted `pre_Y1_dc1e-4` log was renamed `*_attempt1_interrupted.log` and the run was redone. Nothing was deleted.

## Reproduction

```bash
cd research/HDBLAST_CHECKPOINT_20260930/B2_delta1e-3
export OMP_NUM_THREADS=1
P=2 ./queue.sh jobs_r2_stageA.txt      # pre-runs (x_c chain) and C2a8/b/c/d
python3 select_xc.py                    # -> runs/pre/xc_*.json, XC_SELECTION.json
P=2 ./queue.sh jobs_r2_main.txt         # main runs, sensitivity set, C3', C7, Y0 (the first 7 job lines were run through this queue; the rest via slotq.sh)
./slotq.sh jobs_r2_rest.txt             # D6 old-chart continuation + remaining controls (jobs_r2_q3/q4 = the order actually used)
./run_chain.sh runs/explore E1_smax200_Y1_dc1e-2_S1 --delta 0.001 --dstar --Y 1 --dc 1e-2 --grid rho --dzf 3e-4 --dzc 3e-2 --L 24 --Tf 45 --kappa 10 --xc 11.3 --smax 200 --wall 1450
./run_chain.sh runs/explore E4_cfl0.25_Y1_dc1e-2_S1 --delta 0.001 --dstar --Y 1 --dc 1e-2 --grid rho --dzf 3e-4 --dzc 3e-2 --L 24 --Tf 45 --kappa 10 --xc 11.3 --cfl 0.25 --wall 1450
./run_chain.sh runs/explore E3_asinh_Y1_dc1e-2_S1 --delta 0.001 --dstar --Y 1 --dc 1e-2 --grid rho --dzf 3e-4 --dzc 3e-2 --L 24 --Tf 45 --kappa 10 --xc 11.3 --chart asinh --wall 1450
MAXCPU=5400 ./run_chain.sh runs/res S3p_Y1_dc1e-2_dzf7.5e-5 --delta 0.001 --dstar --Y 1 --dc 1e-2 --grid rho --dzf 7.5e-5 --dzc 1.5e-2 --L 24 --Tf 45 --kappa 10 --xc 11.3 --wall 1450
python3 b2_analyze.py                   # -> B2_RESULTS.json
python3 diag/c2_mode_fit.py             # -> diag/D4_C2_MODE_FIT.json
python3 diag/d5_fullmode_truncation.py  # -> diag/D5_FULLMODE_TRUNCATION.json
python3 diag/d7_chart_end.py            # -> diag/D7_CHART_END.json
```

Wall times (one core): S1 Y = 1 runs take about 6 min, S2 about 19 min, C2b about 5 min, and S3′ about 50 min.
