# B3: r(Y) in the tuned model at δ = 0.1 (friction closure)

Ricardo Maldonado's HDBLAST program; prepared with AI assistance, 1 October 2026. Status of every number: **numerical** (floating-point PDE evolutions, internal, not peer reviewed), unless labelled otherwise. Pre-registration: `REGISTRATION.md` (written before any scan run; five dated notes and a closing note appended; that the original text was never edited cannot be proven without version history, see `AUDIT.md`).

## Post-audit corrections (1 October 2026)

This README was written before the independent audit. The sentences below were corrected in place to agree with `AUDIT.md` and with the final labels in `../CLAIM_REVIEW.md`. No registered class changed. No run was repeated for these edits.

- **Y = 3 is conditional, not "reliable".** Its value r ≈ 0.010 is robust, but its PASS-combined class rests on the tertiary resolution level added after the coarser outcomes were seen; the Y = 3 outcome was known before registration; and the dc = 10⁻⁴ repeat was not run. "Minimum r with a reliable class" now reads "smallest r with a (conditional) class".
- **Tertiary convergence order.** The quoted orders 3.8–4.0 (and 3.8–4.1 for the trigger pairs) are measured at the last common reliable time, H₀τ = 43.8 and 52.4, near or after the turnaround (H = 0 at H₀τ ≈ 47.7), not at the plateau epoch. At H₀τ 20–26 the 1.25×10⁻⁴ drift is −0.09% to −0.15% (opposite sign to the coarser drifts), so the apparent order there is undefined. The drift is far below the 5% tolerance.
- **Y = 2 headline qualified.** The pass rests on the pre-registered secondary pair; the 10⁻³ level has no plateau. The A1 r tolerance max(20%, 0.02) is about 100% of r near 0.02, so the same-class requirement does the work. Ω_r = 1.36 at the plateau is helped by Ω_vac = −0.39. It is a numerical, model-internal pass of a snippet-level threshold.
- **Step 1 (D1) values.** At settling D1 gives Ω_r = 1.09–1.12 and Ω_vac = −0.10 to −0.13 (REGISTRATION §0 quotes 1.05 / −0.06, which is wrong). The `conclusions` strings in `D1_Y3_DIAGNOSIS.json` ("~0.75 e-folds", "lapse 10 -> 50") disagree with its own computed fields (0.995–1.016 e-folds friction-off → turnaround; 0.57–0.64 settling → turnaround). The lapse amplification of the drift is a correlation only.
- **E-folds after the plateau, stated per Y.** Measured: 0.383 (Y = 1.5), 0.319 (Y = 2, dc = 10⁻²), 0.227 (Y = 3), 0.093 (Y = 5, unreliable). For Y = 0.5, Y = 0.7 and Y = 2 at dc = 10⁻⁴ the runs stopped before the turnaround, so 0.50 / 0.48 / 0.32 are estimates only (within-run lower bounds 0.17, 0.26, 0.23). "Increasing Y cannot make the radiation phase longer" was an extrapolation; it now reads "in the scanned range a larger Y shortens the expansion left after the plateau".
- **Power-law fit is descriptive.** It includes the conditional Y = 3 point (p = 1.63 without it), and the "mechanism" restates the ledger identity Ra⁴ = Y ∫v²a⁴dτ through the fitted r(Y); it is not an independent explanation.
- **Y = 5 stiffness:** "|λ|dz grows ~ Y" holds only for Y ≳ 3.
- **Minor undisclosed items:** `secondary_trigger.py` and `tertiary_order.py` were modified after all runs (immaterial); the pre-registration "dev check of `b3_analyze.py`" cannot be traced, since that version was not kept.

See `AUDIT.md` and `audit/*.json` for the checks.

## Question

In A1, the tuned-vacuum model (σ = 2W + δ(1 + cφ + d\*φ²/2), δ = 0.1) with the phenomenological friction closure κ₅²j = Yv gave r = 𝒲/rad = 0.463 at Y = 0.3 (FAIL-Weyl), 0.073 at Y = 1 (PASS-conservative) and no plateau at Y = 3. Two questions:

1. Why did the Y = 3 runs never reach a plateau?
2. Does any Y give r ≤ 0.03 (PASS-combined, the combined CMB+BBN+BAO level of the A1 registration)? What is the trend r(Y)?

## Short answer

- **Y = 3 diagnosis (numerical).**
  - The friction source is off long before the runs end.
  - The Weyl term settles at r ≈ 0.010.
  - The A1 runs then lose the plateau to a **4th-order discretisation drift** of 𝒲a⁴. The measured order is 3.8–4.5, and the drift grows while the shell lapse of the bounded chart grows (a correlation; no other chart was run to show causation).
  - The physics also leaves only ≈ 0.6 e-folds between the Weyl term settling and the turnaround that the negative residual vacuum causes. The registered 0.5 e-fold window therefore has to end where the lapse is ≈ 20–25.
  - Finer resolution near the shell late in the run fixes this. I obtained it with validated restarts.
  - The Y = 3 PASS needed a third resolution level, added in a dated note after the coarser results had been seen (see Limitations).
- **Main result (numerical, within the friction closure, δ = 0.1, d = d\*).**
  - **Y = 2 gives a reliable PASS-combined for both seeds (numerical, model-internal):**
    - r = 0.0213 at dc = 10⁻² and r = 0.0213 at dc = 10⁻⁴;
    - the pass rests on the pre-registered secondary pair (5×10⁻⁴, hybrid 2.5×10⁻⁴); the 10⁻³ level has no plateau;
    - the A1 r tolerance max(20%, 0.02) is about 100% of r near 0.02, so the same-class requirement, not the r tolerance, does the work;
    - Ω_r = 1.36 at the plateau is helped by Ω_vac = −0.39 (the negative vacuum lowers H²);
    - the Richardson-extrapolated r is 0.02128–0.02129;
    - this is a pass of a snippet-level threshold inside a toy model, not observational evidence.
  - **Y = 3 gives PASS-combined at dc = 10⁻² on the tertiary resolution pair: conditional.** r = 0.0101, Ω_r = 1.66, Ω_vac = −0.68. The value is robust; the class rests on a resolution level added after the coarser outcomes were seen, the outcome was known before registration, and dc = 10⁻⁴ was not run.
  - Y = 5 is UNRELIABLE (only the finest run has a plateau, r = 0.0044).
- **Answer to question 2 (numerical, model-internal):** yes, r ≤ 0.03 is reached for Y ≳ 1.6–2, established at Y = 2 (both seeds).
- **Trend:** r falls monotonically. A power law r ≈ 0.069 Y^(−1.67) fits Y = 0.3–3 to within 11% (descriptive; it includes the conditional Y = 3 point, and without it p = 1.63); the local slope steepens from 1.42 to 1.84.
- **Strong caveats:**
  - Every plateau sits only a fraction of an e-fold before the shell stops expanding (measured: 0.38 at Y = 1.5, 0.32 at Y = 2, 0.23 at Y = 3, 0.09 at Y = 5 (unreliable); 0.50 and 0.48 at Y = 0.5 and 0.7 are estimates), because the residual vacuum that d\* leaves (≈ −3.5×10⁻⁴ H₀²) turns the expansion around. Ω_r > 1 at the plateau because Ω_vac < 0.
  - This is a "radiation era" in the sense of the A1 rule, not a long radiation-dominated history.
  - The fine tuning of d, the phenomenological Y and the open δ = 10⁻³ case are untouched.

## What was reused from the first (interrupted) attempt

The output folder did not exist when this round started, so nothing from an earlier attempt was found or reused. The solver, static-shell module and analysis are copies of the A1 files:

- `a1_analyze.py` is copied verbatim; its `classify` is the registered rule.
- `evolve_a1.py` has only optional, default-off additions:
  - `--save_T` / `--restart`: state save, and finer-grid or same-grid continuation;
  - `--stop_recollapse`;
  - `--order`, which is exposed but not used.

  With default options it is the A1 solver. C11 confirms this: the new runs reproduce A1 to ≤ 5×10⁻⁹.
- A1 run outputs (Y = 0.3, 1, 3) are read in place, read-only, as context.

## Method

- **Model:** exactly A1 (`REGISTRATION.md` §1). The scan runs use the A1 solver settings, with three numerical changes, each registered with its reason:
  1. **L = 10 instead of 16.** The shell is exact for T < 20. C11 shows L = 10 and L = 16 agree to 5×10⁻⁹.
  2. **Stop at H/H₀ < −0.05.** This never triggered.
  3. **cfl 0.25 for Y = 5 only.** The friction boundary term is stiff, with |λ|dz ≈ 9.8 at Y = 5 (`D2_FRICTION_STIFFNESS.json`); a large real friction eigenvalue exists only for Y ≳ 3. RK4 with Δt = 0.5 dz is unstable for Y ≳ 4: the Y = 5 attempts at cfl 0.5 blew up within 14 steps and are kept in `runs/unstable_cfl0.5/`.
- **Classification:** the A1 rule, unchanged (plateau over Δ ln a = 0.5 with 5% constancy of 𝒲a⁴ and Ra⁴; PASS-combined needs Ω_r ≥ 0.9 and |r| ≤ 0.03). The A1 reliability criterion is also unchanged: same class at two spacings, r agreeing to max(20%, 0.02), near-shell constraint residuals < 0.05.
- **Resolution pairs:**
  - **Primary** pair: dz_fine (10⁻³, 5×10⁻⁴).
  - **Secondary** pair, registered before use: if the primary pair was unreliable and the late drift converged at order 3–5 (`SECONDARY_TRIGGER.json`), the 5×10⁻⁴ run's full state at the last saved T with shell lapse ≤ 3 was moved by quintic-spline interpolation onto a 2.5×10⁻⁴ grid and continued (a "hybrid 2.5×10⁻⁴" run). The pair (5×10⁻⁴, hybrid 2.5×10⁻⁴) was then assessed with the same criterion.
  - **Tertiary** pair (third and fifth dated notes): the same construction one level finer (1.25×10⁻⁴, coarse 10⁻³), for Y = 3 and 5.
- **Restart validation (C10):** a 10⁻³ → 5×10⁻⁴ restart reproduces the full 5×10⁻⁴ run's late drift to 1–3% at H₀τ = 18–24: **PASS** (`C10_RESTART_CHECK.json`).
- **Supplementary (not decisive):** Richardson-extrapolated series on the finest pair, with p = 3.5, 4 and 4.5.

## Results

### Step 1: why A1's Y = 3 had no plateau (`D1_Y3_DIAGNOSIS.json`)

| Finding | Evidence |
|---|---|
| Friction is switched off | Yv²/(4HR) < 10⁻³ after H₀τ ≈ 8.6 (dc 10⁻²) and 13.5 (dc 10⁻⁴); Ra⁴ is then constant to 10⁻⁵ |
| 𝒲a⁴ settles, then drifts upward | Settled r ≈ 0.0101–0.0102 (D1: Ω_r ≈ 1.09–1.12, Ω_vac ≈ −0.10 to −0.13; REGISTRATION §0's 1.05 / −0.06 is wrong) |
| The drift is numerical, at 4th order | Drift(10⁻³)/drift(5×10⁻⁴) = 14–23, order 3.8–4.5; the Weyl identity fails in the late window (p95 0.3) while its v-terms are ~10⁻⁷ |
| The drift grows while the bounded chart's shell lapse grows (correlation only) | Lapse ≈ 10 at settling, ≈ 23 at the earliest possible plateau end, ≈ 35 at H₀τ = 30. No run with a different chart was made, so causation is not shown. (The `conclusions` strings in `D1_Y3_DIAGNOSIS.json`, "~0.75 e-folds" and "lapse 10 -> 50", disagree with its computed fields; the fields are used here.) |
| The plateau window is physically tight | Only 0.57–0.64 e-folds from settling to turnaround |
| 'Non-finite' stops are irrelevant | They occur after H < 0; the A1 "H₀τ 31–45" is the near-shell constraint reliable end |

### Step 2–3: r(Y) scan (`B3_RESULTS.json`; dc = 10⁻²; plateau values from the finest reliable run)

| Y | Class (registered rule) | Pair | r | Ω_r | Ω_vac | e-folds plateau → turnaround | Richardson r (p = 4) |
|---|---|---|---|---|---|---|---|
| 0.3 (A1) | FAIL-Weyl | primary | 0.463 | 0.785 | −0.148 | 0.51 | — |
| 0.5 | FAIL-Weyl | primary | 0.2237 | 0.943 | −0.154 | 0.50 (estimate; run stopped before turnaround, measured ≥ 0.17) | 0.2236 |
| 0.7 | FAIL-Weyl | primary | 0.1326 | 1.036 | −0.173 | 0.48 (estimate; measured ≥ 0.26) | 0.1326 |
| 1 (A1) | PASS-conservative | primary | 0.0734 | 1.124 | −0.207 | 0.44 | — |
| 1.5 | PASS-conservative | primary | 0.0360 | 1.236 | −0.280 | 0.38 | 0.0360 |
| **2** | **PASS-combined** | secondary | **0.0213** | 1.357 | −0.386 | 0.32 | 0.0213 |
| **2, dc = 10⁻⁴** | **PASS-combined** | secondary | **0.0213** | 1.356 | −0.385 | 0.32 (estimate; measured ≥ 0.23) | 0.0213 |
| 3 | PASS-combined (**conditional**: post-hoc tertiary level) | tertiary | 0.0101 | 1.659 | −0.675 | 0.23 | 0.0101 |
| 5 | UNRELIABLE (1.25×10⁻⁴: PASS-combined; 2.5×10⁻⁴, 5×10⁻⁴, 10⁻³: no plateau) | none | (0.0044) | (3.22) | (−2.23) | (0.09) | none |

- **Smallest r with a (conditional) class:** r = 0.0101 at Y = 3 (dc = 10⁻² only; tertiary level added after seeing the coarser outcomes; the dc = 10⁻⁴ repeat was not run). The smallest r with a class that does not depend on a post-hoc level is r = 0.0213 at Y = 2.
- **Is it ≤ 0.03?** Yes, at Y = 2, for both seeds.
- **Fit (numerical):**
  - Over the classified points Y = 0.3–3 (including the conditional Y = 3 and the two A1 points), r ≈ 0.069 Y^(−1.67), with a maximum deviation of 11% (`fit.power_law_reliable`). Without Y = 3 the exponent is 1.63. The fit is an empirical description only.
  - The local slopes are 1.42, 1.55, 1.66, 1.75, 1.83 and 1.84, from low Y to high Y.
  - The fit gives r = 0.03 at Y ≈ 1.64.
  - The unreliable Y = 5 value (0.0044) continues the trend, with a slope of 1.63 from Y = 3.
- **Re-expression of the trend (not an independent explanation; it restates the ledger identity through the fitted r(Y)):**
  - The ledger makes Ra⁴ = Y ∫v²a⁴dτ exactly; the check ratio is 1.0000.
  - So, at low energy (R/σ ≈ 10⁻³), r = (18/σ)(𝒲a⁴/∫v²a⁴dτ)/Y.
  - The Weyl yield per unit a⁴-weighted kinetic integral itself falls as Y^(−0.67) (`fit.trend`): 0.040 at Y = 0.3, 0.0088 at Y = 3.
  - Hence r ∝ Y^(−1.67).
  - In words (an interpretation): friction both diverts the roll's kinetic energy into R (∝ Y) and slows the roll, which reduces the bulk (Weyl) emission per unit kinetic integral.
  - The share of the released energy that ends up in R, 1/(1+r), rises from 0.68 (Y = 0.3) to 0.99 (Y = 3).
- **The price:** the plateau comes later in Y, so Ω_vac at the plateau grows more negative (−0.15 at Y = 0.5, −0.39 at Y = 2, −0.68 at Y = 3, −2.2 at Y = 5), and the expansion left after the plateau shrinks (0.50 at Y = 0.5 is an estimate only; measured 0.319 at Y = 2, 0.227 at Y = 3 and 0.093 at Y = 5 (unreliable); where a turnaround was reached, the estimate ¼ ln((Ω_r+Ω_W)/|Ω_vac|) matches it to 1%). In the scanned range a larger Y shortens the expansion left after the plateau; nothing beyond Y = 5 was tested. The recollapse comes from the residual vacuum and is set by d, not Y. At Y ≳ 5 the registered 0.5 e-fold window barely fits before the turnaround, which is why Y = 5 could not be classified reliably.

### Controls

| Control | Result |
|---|---|
| C5 Weyl identity, correct form | Median residual 5×10⁻⁷ – 4×10⁻⁵ on all finest runs (≤ 10⁻³: pass) |
| C5 4H → 3H | 124–2000× worse (pass) |
| C5 Yv dropped from the scalar junction | 330–156000× worse (pass) |
| C5 R dropped, over **all** records | 2.7–138× worse. **Fails 10× for Y ≤ 1.5 and for Y = 2 at dc = 10⁻⁴** (2.7–8.9×), as A1 found at Y = 1 |
| **C5 R controls (dropped, 3R, −R) on the R-sensitive records** (weight of the R term ≥ 10⁻³; 60–165 records) | **R dropped: 40–675× worse; 3R and −R: 75–1300× worse. Pass on every run, including A1's Y = 0.3, 1, 3** |
| C6 radiation ledger | Median 3×10⁻⁶ – 1.1×10⁻⁴. Pass except Y = 0.5 (1.07×10⁻⁴, just above 10⁻⁴; A1's Y = 0.3 was similar) |
| C10 restart validation | PASS (drift reproduced to 1–3%) |
| C11 L = 10 vs 16 | PASS (5×10⁻⁹) |
| Near-shell constraints up to classification | ≤ 1.6×10⁻⁵ (H), ≤ 1.4×10⁻⁴ (M) on the finest reliable runs; Y = 5 finest run 3.3×10⁻³ / 2.1×10⁻² (still < 0.05) |
| Restart triggers (convergence order of the late drift) | 3.78–4.15 for the primary pairs (`SECONDARY_TRIGGER.json`). For 5×10⁻⁴ vs 2.5×10⁻⁴ and 2.5×10⁻⁴ vs 1.25×10⁻⁴ the quoted 3.8–4.0 (`TERTIARY_ORDER*.json`, "order_at_last_common_reliable") is measured at H₀τ = 43.8 / 52.4 (Y = 3), near or after the turnaround, not at the plateau epoch. At H₀τ 20–26 the 1.25×10⁻⁴ drift is −0.09% to −0.15%, opposite in sign to the coarser drifts, and the order there is undefined. Its size is far below the 5% tolerance |

**Why the A1 "R dropped" control was weak, and why no all-record R control can be ≥ 10× worse at small Y:**

- The R term enters the identity only through k = (σ+R)/6 in 4kwv/6.
- Its weight in the identity's scale is median ~2×10⁻⁵ and at most 0.5%, because R/σ ≤ 3% during the roll and v is small whenever R is large.
- So on most records an error in the R coefficient changes the residual by less than the numerical floor (10⁻⁵ – 10⁻⁴), and the all-record median cannot separate the wrong factor from the correct one.
- On the records where the R term carries ≥ 0.1% of the identity's scale, every R-factor error is 40–1000× worse than the correct identity. That is the discriminating control registered here, and it passes on all runs. It was registered after its values on the A1 runs had been seen; it was blind on the new runs.

## Limitations

- **Phenomenological closure.** All results are for the friction closure κ₅²j = Yv, not the derived χ source. Y is a free parameter here: the result says which Y would be needed, not that the model has it.
- **Fine tuning.** d = d\* (1 + c + d/2 ≈ 0 tuning) is unchanged. A1 showed a 1% change in d destroys the radiation era.
- **Very short "radiation era".** All plateaus end with Ω_vac < 0 and a turnaround 0.09–0.38 e-folds later where measured (0.5 for Y = 0.5 and 0.7 is an estimate). Ω_r > 1 reflects the negative vacuum, not extra radiation. Nothing here is a realistic thermal history.
- **δ = 0.1 only.** The registered-scale δ = 10⁻³ is not addressed.
- **Thresholds.** The N_eff thresholds and the r → ΔN_eff mapping are those of the A1 registration: snippet-level literature values, background-only.
- **Hybrid runs.** The decisive Y = 2 classes come from the secondary pair, a registered method change. The hybrid 2.5×10⁻⁴ runs share their pre-restart history (lapse ≤ 3) with the 5×10⁻⁴ parent. The primary pair at Y = 2 was unreliable only because the 10⁻³ run's 4th-order drift destroys its plateau; the Richardson values agree with the hybrid to 10⁻⁴.
- **Outcomes seen before registering.** The Y = 3 outcome (r ≈ 0.010) was known before this registration (disclosed in §0). The scan values Y = 0.5, 0.7, 1.5, 2 and 5 were blind. The pre-registration "dev check of `b3_analyze.py`" cannot be traced (that version was not kept), and `secondary_trigger.py` and `tertiary_order.py` were modified after all runs (not disclosed at the time; immaterial, since the trigger orders are far from the 3/5 bounds).
- **Y = 3 tertiary pair.** Both tertiary runs restart from the same 5×10⁻⁴ state and keep dz_coarse = 10⁻³, so the pair does not test the coarse-region or pre-restart error (the audit estimates the latter at ~10⁻⁵).
- **Not done (budget):**
  - Y = 3 at dc = 10⁻⁴. It would need the same three-level restart chain. The A1 L = 16 Y = 3 dc = 10⁻⁴ pair is unreliable for the same numerical reason.
  - The dc = 10⁻⁴ repeat for Y = 1.5 (PASS-conservative). In every completed case (A1 Y = 0.3 and 1; B3 Y = 2), dc = 10⁻⁴ only shifts the time axis and changes r by < 10⁻³ relative.
- **Y = 5.** It is unreliable even with three resolution levels. Its plateau falls 0.09 e-folds before the turnaround, at shell lapse ≈ 40, where only the 1.25×10⁻⁴ run holds 𝒲a⁴ to 5%.
- **Resolution levels were added after the registration.** The Y = 3 and Y = 5 classes needed resolution levels that were added in dated notes after the registration, once the coarser outcomes had been seen. The Y = 3 PASS-combined therefore rests on a method extension decided after seeing data. The Y = 2 PASS-combined (secondary pair, registered up front) does not.

## Files

- `REGISTRATION.md`: pre-registration and dated notes.
- `README.md`: this file.
- `B3_RESULTS.json`: all classifications, pairs, Richardson, C5, the fit and trend.
- `D1_Y3_DIAGNOSIS.json`, `d1_y3_diagnosis.py`: step 1.
- `D2_FRICTION_STIFFNESS.json`, `d2_friction_stiffness.py`: Y = 5 stability.
- `SECONDARY_TRIGGER.json`, `secondary_trigger.py`; `TERTIARY_ORDER.json`, `TERTIARY_ORDER_Y5.json`, `tertiary_order.py`: restart triggers.
- `C10_RESTART_CHECK.json`, `c10_restart_check.py`; `C11_L_CHECK.json`, `c11_L_check.py`.
- `evolve_a1.py`, `static_w.py`, `a1_analyze.py`: solver and rule.
- `b3_analyze.py`: B3 analysis.
- `make_restart_job.py`: the registered restart-time rule.
- `queue.sh`, `chain_b.sh`, `chain_c.sh`, `chain_d.sh`, `jobs_*.txt`: job lists. `jobs_chain_e.txt` holds the Y = 5 restarts. `jobs_dc4_main_Y3_NOT_RUN.txt` is the dropped job. `b3_analyze.out` is the printed table.
- Runs:
  - `runs/pre`, `runs/main`, `runs/fine`, `runs/ctl`: summaries, time series, logs.
  - `runs/main/states`: saved full-field states.
  - `runs/unstable_cfl0.5`: failed Y = 5 cfl 0.5 attempts.
  - `dev/`: the one disclosed development run.
  - `diag/`: helper scripts and `DEV_COARSE_REGION_L_CHECK.json`.

## Reproduction (from this folder; one core, `OMP_NUM_THREADS=1`)

```
python3 d1_y3_diagnosis.py                         # step 1 (reads A1 runs)
python3 d2_friction_stiffness.py
./queue.sh jobs_pre.txt; ./queue.sh jobs_main.txt  # scan pairs (Y = 5 lines fail at cfl 0.5; see runs/unstable_cfl0.5)
./queue.sh jobs_chain2.txt; ./queue.sh jobs_dc4_pre.txt; ./queue.sh jobs_dc4_main.txt; ./queue.sh jobs_y5_main.txt
python3 secondary_trigger.py; python3 make_restart_job.py <parent_tag> runs/fine <new_tag> 2.5e-4 1e-3   # hybrids (see jobs_chain_c.txt, jobs_chain_d.txt)
./queue.sh jobs_chain_c.txt; ./queue.sh jobs_chain_d.txt; ./queue.sh jobs_chain_e.txt
python3 c10_restart_check.py; python3 c11_L_check.py; python3 tertiary_order.py 3; python3 tertiary_order.py 5
python3 b3_analyze.py                              # -> B3_RESULTS.json
```

Single-run example: `python3 evolve_a1.py --out runs/main --tag main_Y2_dc1e-2_dzf5e-4 --dstar --Y 2 --dc 1e-2 --dzf 5e-4 --dzc 2e-3 --L 10 --Tf 20 --kappa 10 --project 9.9 --stop_recollapse --save_T 0.25 --xc_from runs/pre/pre_dstar_Y2_dc1e-2_summary.json` (≈ 7 min on one shared core).
