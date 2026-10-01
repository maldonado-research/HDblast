# Critic notes: HDBLAST checkpoint 30 September 2026

Date: 1 October 2026. This was the final consistency review of the top-level files (00_READ_FIRST.md, CLAIM_REVIEW.md, NEXT_TESTS.md, RESULTS_SUMMARY.json, REPRODUCE.md). I compared them with each workstream's README.md, REGISTRATION.md and AUDIT.md, and with the JSON outputs (B1_*.json, B2_RESULTS.json, XC_SELECTION.json, diag/D4_C2_MODE_FIT.json, B3_RESULTS.json, B4_*.json, and the audit/*.json files).

I ran no new computation. I edited only the top-level files and added this file. Workstream folders were not modified.

## Checks made

| Check | Outcome |
|---|---|
| AI model names in any .md/.py/.json/.txt/.sh/.out/.log file | none found (the phrase "prepared with AI assistance" appears; no model or vendor is named) |
| Personal information (e-mail, address, phone) | none found; the program owner's name appears as attribution only |
| "proves", "proof", "discovery", "confirms the hypothesis" | none used as a claim (only in negations, such as "no claim of discovery or proof", and in audit wording such as "does not prove") |
| Registrations edited after the fact | Birth times precede the first decisive runs in all four workstreams (per the audits, re-checked with `stat`). Body text cannot be proven unedited (git was not used). Known post-hoc items are already listed in CLAIM_REVIEW.md: B2 Amendment 2 reinterpretation; B3 tertiary level; B4 rtol pair and e-fold criterion; B1 second-note timing unverifiable. The B3 REGISTRATION.md mtime (12:33) is later than its README birth (12:06), which is consistent with the appended closing note |
| Top-level numbers vs JSON | Spot-checked all headline numbers: B3 r per Y, Richardson, fit p and Y(r=0.03), Omega values; B1 K2 r, b_cut r range, explore r; B2 C2 rates and x_c values; B4 T1-T8, Lambda_res, d*_exact, tolerances, BBN rows. They match, except for the items corrected below |

## Changes made to the top-level files

1. **B2 C2 post-hoc values (00_READ_FIRST, RESULTS_SUMMARY).**
   - The old wording was "post-hoc fits recover the target 1.65719 to about 3e-5" and "1.65718-1.65722".
   - Those values hold only on the rho-grid runs C2a, C2b and C2c.
   - On C2d (tanh grid), D4 `lam_range` is 1.657047-1.657089, and the audit derivative estimator gives 1.657048, about 1.4e-4 low.
   - The supplementary late-window estimator gives 1.65707-1.65711.
   - The wording now gives all of these and says that the explanation is not formally accepted under Amendment 2.
2. **B3 Y = 2 headline.**
   - Added that the pass rests on the pre-registered secondary pair and that the 1e-3 level has no plateau.
   - Added that the A1 r tolerance max(20%, 0.02) is about 100% of r near 0.02, so the same-class requirement does the work.
   - Added that Omega_r = 1.36 at the plateau is helped by Omega_vac = -0.39.
   - Added that this is a model-internal pass of a snippet-level threshold.
   - "Yes, a Y passes" is now labelled numerical and model-internal.
3. **B3 Y = 3.** Added that its outcome (r ~ 0.010) was known before the B3 registration, and that dc = 1e-4 was not run.
4. **B3 e-folds.**
   - The old wording was "Each radiation era lasts only 0.1 to 0.5 e-folds" and "Plateau 0.09-0.38 e-folds (measured)".
   - These are now stated per Y, as the expansion left after the plateau time.
   - Y = 5 is marked unreliable and Y = 3 conditional. The Y = 0.5 and 0.7 values are marked as estimates.
   - The garbled CLAIM_REVIEW row (0.383 attributed to "Y = 0.5 estimate") is fixed: 0.383 is Y = 1.5.
5. **B4.**
   - "r wrong by about 12x" applies to Y = 1 only; at Y = 0.3 it is 3.6x (B4_VALIDATION T1, T2).
   - "The error comes from the matching step, not from the 4D roll" is softened. These are post-registration checks, the EFT tracks 5D to about 1-5% up to phi_b ~ 0.7, and v is 37% low at the match point (B4 AUDIT).
   - The one-e-fold tolerance now shows its input dependence: 1.3e-4 with Y = 1 inputs, 8.7e-5 with Y = 0.3 inputs.
   - The BBN tolerance now notes that the d*_exact fit spread is about 1e-6.
6. **B1.**
   - "r is at least 94" is now "run-end r = 94 to 8.5e3 (y = 1), not a registered plateau value".
   - The corrected b_need control values (0.316, 0.082, 0.169 vs 0.170) are now marked as coming from the B1 audit JSON (`audit/AUDIT_TIMESERIES.json`). They lie beyond each run's reliable end. B1_RESULTS.json still holds the mis-paired values.
7. **B2.1.** Added that the dc = 1e-4 value x_c = 18.7 is single-resolution.
8. **RESULTS_SUMMARY.json.**
   - Added the missing B2.2 entry.
   - Added a note that the "registered_model" block (delta = 1e-3) differs from the delta = 0.1 / MODEL CHANGE A1 setting of B1, B3 and B4.
   - Added caveat fields (pair used, Omega_vac, fit scope, tolerance range).
   - Pointed to this file.
9. **NEXT_TESTS.md.** "A larger Y cannot lengthen the era" was an extrapolation beyond Y = 5. It now reads "in the scanned range a larger Y shortens the expansion left after the plateau".
10. **REPRODUCE.md.** Added the B2 diagnostic scripts d1/d2/d3, whose JSON is cited in the B2 README, and listed CRITIC_NOTES.md.

## Unresolved issues (not fixable from the top level)

- **Workstream READMEs were not updated after their audits.** The top-level files follow the audits, but the producer READMEs still contain the audited errors:
  - **B1 README:** "essentially final" for all y = 1 cells; the mis-paired control values (-338, -1149, 22.5).
  - **B2 README:** "four dated notes" (there are five); "x_c ruled out"; "every Y = 1 run breaks down after the crossing"; the Amendment 2 reinterpretation; "cannot fix" labelled numerical.
  - **B3 README:** Y = 3 called a "reliable PASS-combined" and "minimum r with a reliable class" (the top level labels it conditional); "order 3.8-4.0" measured near turnaround; the Omega values in REGISTRATION section 0 and the D1 conclusion strings.
  - **B4 README:** "W a^4 falls about 15x" (the audit measures 33.5x / 11.2x; this is in B4 audit/audit_a1_series.json, not in a producer JSON); "EFT tracks 1-3%"; "Omega_vac ~ -0.011 at the stop" (A1 gives -0.160); the registered model E tuning windows that exclude d* are omitted.

  A reader who opens a workstream README without its AUDIT.md will see stronger statements than the checkpoint supports.
- **Time series are not stored.** *.npz and *.pkl are git-ignored. Every PDE claim can be re-checked only by re-running the solvers (REPRODUCE.md). The audits' independent re-analyses depended on local files that will not travel with the checkpoint.
- **Single-seed or single-resolution items.**
  - B2 x_c = 18.7 (dc = 1e-4).
  - B3 Y = 3 (dc = 1e-2 only; tertiary pair from a shared restart state and a shared dz_coarse).
  - B3 Y = 1.5 (dc = 1e-2 only).
  - B1 exploratory cell (one b, two spacings, unregistered).
- **The B3 power-law fit** includes the conditional Y = 3 point and two A1 points. Its exponent is descriptive, and the local slopes vary from 1.42 to 1.84.
- **Inherited and unverified.** The N_eff thresholds and the r -> Delta N_eff mapping come from A1 at snippet level and were not re-verified in this round (NEXT_TESTS item 7).
- **Unverifiable registration details**, as listed by the audits: the B1 second-note timing; that the original registration bodies were never edited; the B3 "dev check of b3_analyze.py" before registration, whose earlier version was not kept.
- **Repeated attribution line.** "Prepared with AI assistance" appears in the workstream READMEs and registrations. It names no model. The author may wish to keep or reword it.
