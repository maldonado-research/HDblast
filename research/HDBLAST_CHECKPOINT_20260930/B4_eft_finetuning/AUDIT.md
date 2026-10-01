# Independent audit of B4 (EFT cross-check and fine-tuning measure)

Audit date: 1 October 2026. The audit had one core. I wrote only this file and the `audit/` subfolder. I modified no producer file.

## What was done

1. **Re-runs.** All four producer scripts were copied into `audit/` and re-run from scratch there (`OMP_NUM_THREADS=1`):
   - `static_lambda.py`
   - `b4_validate.py`
   - `b4_predictions.py`
   - `b4_finetuning.py`

   Every output matches the producer's JSON. `B4_VALIDATION.json` primary, `B4_PREDICTIONS.json` and `B4_FINETUNING.json` are identical apart from the script hash. Logs: `audit/static.log`, `audit/validate.log`, `audit/pred.log`, `audit/fine.log`.
2. **Independent extrapolation check (`audit/audit_extrap.py` → `audit/audit_extrap.json`).**
   - Λ_res(d) was refit on four subsets (q ≥ 0.80, 0.85, 0.90 and 0.94), with polynomial degrees 2, 3 and 4.
   - d*_exact was also obtained by inverse interpolation d(H²) evaluated at H² = 0.
3. **Independent read of the A1 5D time series (`audit/audit_a1_series.py` → `audit/audit_a1_series.json`).** This used A1's own `plateau_index` and read the series only. It checked:
   - the late vacuum value and φ_s;
   - the plateau values of r and Ω_vac;
   - W a⁴ at φ_b = 0.95 compared with the plateau;
   - the e-folds with |Ω_vac| ≤ 0.1 after the scalar stops.
4. **Algebra re-derived by hand and with sympy.**
   - The Λ < 0 radiation solution: u = a², (u′/2)² = K − L u², u = √(K/L) sin(2√L t). This gives H = 0 at π/(4√L) and the crunch at π/(2√L).
   - The M8 series root d* = −3.106933 and the series slope 0.122558 H0².
   - The BBN tolerance arithmetic for the 10 MeV and 10¹⁶ GeV rows.
5. **File times (birth and modification).**
   - `REGISTRATION.md` was born at 10:40:18. That is after `dev/hj_effective_theory_chat9_copy.py` (10:37:53, a disclosed pre-registration check) and before `static_lambda.py` (10:40:50), `eft_reduced.py` (10:51:49) and the first `B4_VALIDATION.json` (10:54:32).
   - Its last modification, at 10:56:10, follows the final validation write (10:55:51). That fits the appended dated notes.
   - The original body text cannot be verified as unedited, because there is no snapshot and git use was not allowed. I saw nothing in the text that suggests an edit; the targets match A1_RESULTS and the code.

## Verdicts on claims

| # | Claim | Verdict | Basis |
|---|---|---|---|
| 1 | C1: EFT growth rate 1.65750 vs 1.65719; μ² = −7.71980 | confirmed | Re-run reproduces it. Note that the EFT rate does not depend on δ, so this is a check of the EFT functions against a δ = 1e-3 5D number, not of the δ = 0.1 model |
| 2 | C2: static solver at d = 0 reproduces M2 | confirmed | Re-run gives 0.4061685 and φ = 0.9914347 |
| 3 | Λ_res(d*) = −3.399e-4 H0², d*_exact = −3.10415, slope 0.12209, Λ(1.01 d*) and Λ(1.05 d*) | confirmed | Re-run is identical. Alternative fits give Λ(d*) between −3.3992e-4 and −3.3993e-4 (degree ≥ 3) and d*_exact between −3.104148 and −3.104149. A1 Y = 0.3 late vac is −3.4000e-4 and φ_end = 1.036102. Caveat: A1 Y = 1 late vac drifts to −3.50e-4 (dz 5e-4) and −4.63e-4 (dz 1e-3) by H0τ ≈ 29, so the agreement holds for the stable plateau segment and the Y = 0.3 run |
| 4 | Registered verdict NOT VALIDATED: T1–T4 fail, T5–T8 pass | confirmed | Re-run is identical: T1 1.6705, T2 0.9147, T3 −0.0407, T4 −0.0315, T5 14.00, T6 19.14, T7 7.37, T8 8.61. Under the registered rule, T3/T4 failing alone forces NOT VALIDATED, so the seed note (dated note 1) cannot have changed the verdict. The A1 targets match A1 series (plateau r 0.0734 / 0.4629, Ω_vac −0.207 / −0.148) |
| 5 | Controls: φ_m systematics, rtol, C3 wrong-formula sensitivity | partially confirmed | Numbers reproduce. C3 discriminates: T2 goes 0.915 → 1.27 and 6.99. However, the registered rtol pair was 1e-8 vs 1e-11 and the code used 1e-10 vs 1e-12, with no dated note. This is harmless, since the differences are 8e-13, but it is an undisclosed deviation |
| 6 | Post-registration: the failure is in the matching (recipe on 5D gives r = 0.94 / 1.60; EFT tracks 5D to 1–3% up to φ ≈ 0.7; R 15% low at 0.95) | partially confirmed | r = 0.940 and 1.595 reproduce. Two corrections: (a) at φ = 0.7 the R and v errors are 4.5%, not 1–3%; at φ = 0.95, v is 37% low and H 9.5% low. (b) **"W a⁴ falls about 15× during the final stop" is not in any saved JSON** (the field is `null`). From the A1 series, W a⁴(φ_b = 0.95)/W a⁴(plateau) = **33.5 at Y = 1** and **11.2 at Y = 0.3**. The qualitative conclusion stands |
| 7 | Post-registration variants (Einstein dump r = 2.57 / 6.15; radiation-only: T3, T5–T8 pass, T4 fails) | confirmed | Re-run is identical; labelled post-registration |
| 8 | r(Y) scan: no usable prediction | confirmed | Re-run is identical; correctly labelled |
| 9 | Roll durations and radiation production vs Y | confirmed (as conditional) | Re-run is identical |
| 10 | Tuning window d/d* ∈ [1.000, 1.0075], [0.9975, 1.015], [0.9975, 1.0175] | partially confirmed | Numbers reproduce, but they come from the **post-registration radiation-only variant**, chosen after it passed the timing targets. Section 3.2 registered the window for model E. **The registered model's windows are [1.015, 1.0225] (Y = 0.3), [1.0125, 1.0275] (Y = 1) and [1.015, 1.03] (Y = 3). They exclude d* itself, and Y = 3 is cut off at the grid edge.** These are in `B4_PREDICTIONS.json` but are left out of the README and the summary. Also, the second window criterion was registered as "≥ 1 e-fold with Ω_r ≥ 0.9 and H > 0". The code (`rd_efolds`) adds \|Ω_vac\| ≤ 0.1 without a dated note. The half-width ≈ 2e-4 reproduces |
| 11 | Late Λ < 0 radiation solution (exact-verified) | confirmed | Derived independently. The sympy residual is 0 and the ODE agrees to 4e-11 |
| 12 | Why Ω_vac < 0 at d*; N_RD = 0.52; H = 0 at ≈ 44.2, crunch at ≈ 86.8; 13.7 / 6.9 at 1.01 / 1.05 d* | partially confirmed | Measured directly from A1, the e-folds after the stop with \|Ω_vac\| ≤ 0.1 are **0.48** (Y = 1) and 0.46 (Y = 0.3), so "≈ 0.5 e-fold" stands. Re-estimating from the plateau state (τ = 12.92) gives H = 0 at ≈ 43.9, consistent. These times are closed-form extrapolations beyond the end of the A1 run (28.8), so they are estimates. **Incorrect statement:** "Ω_vac ≈ −0.011 when the scalar stops (A1, H0τ ≈ 4.3)" is Λ_res/H², not A1's Ω_vac. A1's actual Ω_vac at that record is −0.160, because the local vac (−4.8e-3) has not yet relaxed to Λ_res, which takes until H0τ ≈ 7. A1's Ω_vac passes through about −0.03 near τ ≈ 5.3 and then grows to −0.21 |
| 13 | One-e-fold tolerance \|Δd/d*\| ≤ 1.3e-4 | confirmed (numerical, with inputs at the A1 "stop" record) | Arithmetic re-done. rad_s is taken at a stop record where vac has not relaxed. Its effect is O(1) on Λ_tol and only logarithmic in N |
| 14 | BBN tolerance 7e-7 … 7e-80; observed Λ 5e-41 … 5e-114 | confirmed as conditional | Arithmetic is correct (ρ_c = 8.0992e-47 h² GeV⁴). Correctly labelled conditional and framed as a restatement of the CC problem. Note: the numerical d*_exact has a fit spread of about 1e-6 relative (degree 2 vs degree ≥ 3 inverse fits), which is larger than the 10 MeV tolerance. The statement is only about distance from the exact zero |
| 15 | Restart note: no files from a first attempt | unverifiable | The earliest file is 10:37:53 and the parent folder is 10:33. Nothing contradicts the note |

## Registration

- **Respected for the verdict.** It was written before any reduced-model trajectory (birth times above) and discloses prior exposure to A1. The verdict follows the registered rule, and post-registration material is labelled.
- **Undisclosed deviations, all minor:**
  - the rtol pair;
  - the extra |Ω_vac| ≤ 0.1 condition in the e-fold window criterion;
  - headlining the post-registration radiation-only window while omitting the registered model's window, which excludes d*.
- The first validation output, from before dated note 1, was overwritten and not kept.

## Corrections requested

1. Replace "W a⁴ falls about 15×" with the measured factors (33.5× at Y = 1, 11.2× at Y = 0.3, from φ_b = 0.95 to the plateau), and save them in JSON.
2. Fix "Ω_vac ≈ −0.011 when the scalar stops": A1's Ω_vac at H0τ = 4.28 is −0.16. The −0.011 is Λ_res/H².
3. Report the registered model E's tuning windows ([1.015, 1.0225], [1.0125, 1.0275], [1.015, ≥1.03]) next to the radiation-only ones, and add a dated note on the variant choice and the changed e-fold criterion.
4. Add a dated note on the rtol pair (1e-10 / 1e-12 used vs 1e-8 / 1e-11 registered).
5. "EFT tracks 5D to 1–3% up to φ ≈ 0.7" should read about 1–5% (R and v are 4.5% off at 0.7). Also state that v is 37% low at the match.

## Overall

The core results stand and reproduce exactly:
- the NOT VALIDATED verdict;
- the exact residual vacuum Λ_res(d*) ≈ −3.40e-4 H0² and d*_exact = −3.10415, robust to fit choice;
- the exact Λ < 0 late solution;
- the ~0.5 e-fold transient radiation era;
- the fine-tuning restatement.

Several secondary descriptive statements are wrong or unsupported (items 6, 10 and 12), and the tuning-window presentation favours a post-registration variant. None of these changes the headline.
