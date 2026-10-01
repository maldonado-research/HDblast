# HDBLAST checkpoint, 30 September 2026 (written 1 October 2026)

This checkpoint follows the A1 tuned-vacuum verdict
(`research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/A1_VERDICT.md`). It has four workstreams, B1 to B4. Each one was produced and then independently audited. Each subfolder holds a README.md, REGISTRATION.md (pre-registration, with dated notes appended) and AUDIT.md. Where a producer and its audit disagree, this checkpoint follows the audit. The reconciliation is in CLAIM_REVIEW.md.

**Status.** Nothing here is observational evidence. Nothing here has been peer reviewed. All results are numerical or algebraic statements about a toy 5D model, and in most cases about the tuned variant (MODEL CHANGE A1: d = d*). Most of them also use a phenomenological friction closure.

## Plain-language summary

The hypothesis is that a five-dimensional "blast" left behind a hot, radiation-filled universe. The model has a 4D world (the "shell") moving through a 5D space. A scalar field on the shell rolls downhill. Two things then compete:
- **Real radiation.** This is made when the rolling field loses energy to matter.
- **A "dark radiation" term (the Weyl term).** It is left over from the 5D geometry. Observations allow only a little of it.

The figure of merit is **r = Weyl / radiation**. The A1 rules call r <= 0.1 "conservative pass" and r <= 0.03 "combined pass". A1 also found that the leftover vacuum energy must be tuned through a coefficient d, or there is no radiation era at all.

Results of the four workstreams:

1. **B1: can a real matter field supply the radiation?** A1 put the energy loss in by hand as a "friction" with strength Y. B1 replaced that with a derived particle-production channel (chi field).
   - The registered outcome is **INCONCLUSIVE**, because no run lasted long enough to reach the plateau window.
   - Descriptively, the derived channel makes far too little radiation when the chi particles are lighter than the 5D cutoff: at the run end (not a registered plateau) r = 94 to 8.5e3 for y = 1, against a target of 0.1 or less.
   - Heavier chi helps: one exploratory cell at 5.4x the cutoff gives r = 0.28. The regime that might reach r <= 0.1 has strong backreaction, and there the simulations are not reliable.
   - **The derived chi source cannot yet replace the friction closure.**
2. **B2: the realistic detuning delta = 1e-3.**
   - B2 found and fixed the reason A1's pre-run failed: an undamped gauge was unstable.
   - B2 also proposed an explanation of the failed growth calibration: a transient from the seed offset. Post-hoc estimators support it on the rho-grid runs C2a/C2b/C2c (three-term mode fit 1.65718-1.65725; audit derivative estimator 1.657184-1.657195), against the target 1.65719. On the C2d (tanh) run both give 1.65705-1.65709, about 1.4e-4 low, and the supplementary late-window estimator gives 1.65707-1.65711 on C2a/C2b. By Amendment 2's own acceptance criteria, as written, the explanation is not formally accepted, and the registered C2 still FAILS.
   - Every main run still stops before the field finishes rolling. **delta = 1e-3 remains INCONCLUSIVE.** A new gauge is needed for late times.
3. **B3: how strong must the friction be?** B3 scanned Y at delta = 0.1.
   - **Y = 2 passes the combined bound (r = 0.0213) for both seeds** (numerical). This is the main registered positive result of the round. It rests on the pre-registered secondary resolution pair (5e-4 and a restarted 2.5e-4 run); the coarsest level (1e-3) shows no plateau, although its r at the plateau time (0.024) is also below 0.03. It is a model-internal pass of a snippet-level threshold under a phenomenological friction closure, not an observational test.
   - Omega_r >= 0.9 is met at the plateau partly because the negative residual vacuum lowers H^2 (Y = 2: Omega_r = 1.36, Omega_vac = -0.39); it is not a long radiation-dominated expansion.
   - Y = 3 gives r = 0.0101, but only conditionally, because it relies on a resolution level added after the coarser outcomes were seen, and the Y = 3 value (r ~ 0.010) was already known from A1 data before B3 was registered.
   - Overall r falls roughly as Y^-1.67 and crosses 0.03 near Y ~ 1.6.
   - After the plateau time only 0.09 to 0.38 e-folds of expansion remain (measured, Y = 1.5 to 5; about 0.5 estimated at Y = 0.5) before the negative leftover vacuum turns the expansion around.
4. **B4: a cross-check with a simpler 4D theory, and how much tuning is needed.**
   - The registered simple model was **NOT VALIDATED**: it gets r wrong by about 12x at Y = 1 and 3.6x at Y = 0.3. Post-registration checks place the main error in the sudden matching step: the 4D roll tracks 5D to about 1-5% up to phi_b ~ 0.7, but v is 37% low at the match point, so the roll is not error-free either.
   - The exact static calculation explains A1's leftover vacuum: A1's d* came from a series, and the exact zero is 0.09% away.
   - A radiation era lasting down to nucleosynthesis needs d tuned to a relative precision of about 7e-7 (reheating at 10 MeV) down to 7e-80 (1e16 GeV) (conditional; measured from the exact zero, whose own numerical fit spread is ~1e-6).
   - **This is the cosmological-constant problem restated, not solved.**

**In short:** the model can meet the Weyl bound only if (a) the vacuum coefficient is extremely fine-tuned and (b) an energy-loss channel with an effective strength Y of about 1.6 or more exists. No derived channel tested so far provides (b). Even when both hold, the radiation era is short-lived at delta = 0.1.

## Results table

| ID | Result | Final label |
|---|---|---|
| B1.1 | Closed shell system with chi source (junctions, Ward identity, ledgers, Weyl balance): 22 identities, 8 controls | exact-verified |
| B1.2 | Calibrations K1 (chi off = A1 Y=0) and K2 (toy = A1 Y=1, r = 0.0738) | numerical |
| B1.3 | Registered verdict at lambda_c = 1: no plateau in any of the 32 runs | inconclusive |
| B1.4 | At the M5 cutoff with y=1, run-end r = 94 to 8.5e3 and Omega_r <= 0.012 (no plateau reached). The phi*=0.9, G=100 cell is not fully decayed. The y=0.1 values are snapshots, not bounds (r could fall to ~50) | numerical (descriptive) |
| B1.5 | The reduced estimate gives FAIL-Weyl or recollapse at b_cut; it agrees with 5D to 3-4% in end-value r | conditional |
| B1.6 | lambda_c needed for r <= 0.1 is ~11-50; unreliable (O(1) backreaction) | conditional |
| B1.7 | Exploratory cell at lambda_c = 5.36: r = 0.282, not a plateau | numerical (exploratory, unregistered) |
| B1.8 | The derived chi channel cannot replace the friction closure with chi masses below M5 | negative (conditional/descriptive) |
| B1.9 | b_need regime unreliable; production-model sensitivity is about 2x at the run end (not orders of magnitude). These run-end values lie beyond each run's reliable end and come from the audit (B1 audit/AUDIT_TIMESERIES.json); B1_RESULTS.json still holds the mis-paired values | inconclusive |
| B2.1 | x_c = inf in A1 came from an unstable kappa=0 pre-run; kappa=10 gives x_c = 11.3 (dc=1e-2); the dc=1e-4 value 18.7 is single-resolution | numerical |
| B2.2 | A1 spacing stops were numerical failures in reliability. The end of the dz=2e-4 run lies near the physical crossing | numerical (partially confirmed) |
| B2.3 | Registered C2 calibration fails (1.65630, spread 2.4e-3) | negative |
| B2.4 | The seed-offset explanation of the C2 failure; post-hoc mode fit 1.65718-1.65725 and audit derivative estimator 1.657184-1.657195 on C2a/C2b/C2c, but 1.65705-1.65709 on C2d; late-window estimator 1.65707-1.65711 (C2a/C2b); target 1.65719. Amendment 2's own criteria are not met as written, so it is not formally accepted | numerical (post hoc) |
| B2.5 | Registered tuned test at delta = 1e-3 | inconclusive |
| B2.6 | Breakdown about 1-2.5 H0^-1 after the light-cone crossing (dc=1e-2 runs and sensitivity set only); x_c is NOT ruled out | numerical (diagnosis) |
| B2.7 | "Plain refinement cannot fix it" (+0.65 H0tau per doubling) | conditional (3-point extrapolation) |
| B3.1 | Y=0.5 FAIL-Weyl r=0.224; Y=0.7 FAIL-Weyl r=0.133; Y=1.5 PASS-conservative r=0.036 | numerical |
| B3.2 | **Y=2 PASS-combined, both seeds, r = 0.0213** (registered secondary pair; r at the plateau time converges at 4th order across 4 levels, but the coarsest level has no plateau; the A1 r-agreement tolerance max(20%, 0.02) is ~100% of r here, so the same-class requirement does the work) | numerical |
| B3.3 | Y=3 PASS-combined at dc=1e-2 only, r = 0.0101 (post-hoc tertiary level; outcome known before registration; dc=1e-4 not run) | conditional |
| B3.4 | Y=5: finest-level r = 0.0044, not converged | inconclusive |
| B3.5 | r ~ 0.069 Y^-1.67 (p = 1.63 without Y=3); r = 0.03 at Y ~ 1.64-1.66 | numerical (fit) |
| B3.6 | Expansion left after the plateau time before turnaround: 0.38 (Y=1.5), 0.32 (Y=2), 0.23 (Y=3, conditional), 0.09 (Y=5, unreliable) e-folds, measured; Y=0.5/0.7 values are estimates only; it shrinks with Y in the scanned range | numerical |
| B4.1 | Registered reduced EFT model reproduces A1: T1-T4 fail (r 3.6x high at Y=0.3, 12x at Y=1; Omega_vac off by >0.1), T5-T8 pass | negative (NOT VALIDATED) |
| B4.2 | Failure is in the sudden matching step (recipe on exact 5D gives r = 0.94) | numerical (post-registration) |
| B4.3 | Lambda_res(d*_M8) = -3.40e-4 H0^2; d*_exact = -3.10415 (0.09% from d*) | numerical |
| B4.4 | Lambda<0 radiation solution: H=0 at pi/(4 sqrt|L|), crunch at pi/(2 sqrt|L|) | exact-verified |
| B4.5 | A1 d* PASS has ~0.5 e-fold with |Omega_vac| <= 0.1; H=0 ~44 and crunch ~87 H0^-1 are extrapolations | numerical / conditional |
| B4.6 | Tolerance for 1 radiation e-fold: |dd/d*| <= 1.3e-4 (A1 Y=1 inputs; 0.87e-4 with Y=0.3 inputs) | numerical |
| B4.7 | Tolerance down to BBN: 7e-7 (10 MeV) to 7e-80 (1e16 GeV) | conditional |
| B4.8 | Tuning windows in d: the radiation-only variant includes d*; the REGISTERED model's windows exclude d* | conditional |

## What changed relative to the A1 verdict

- **Minimum r over Y.** A1 had Y=1 PASS-conservative at r = 0.073 and Y=3 with no plateau. B3 now gives a registered PASS-combined at Y=2 (r = 0.0213, both seeds). The minimum reliable registered r is 0.0213 at Y=2; the conditional minimum is 0.0101 at Y=3. Y=5 (0.0044) is not converged. **Yes, a Y passes the combined bound (r <= 0.03), with the friction closure only, at delta = 0.1 and d = d* (numerical, model-internal; Omega_r >= 0.9 there is helped by Omega_vac < 0).** A1's missing Y=3 plateau is diagnosed as discretisation drift. Its link to the bounded-chart lapse is a correlation only.
- **Friction closure.** In A1, Y was a free phenomenological parameter. B1 tested the derived chi source. At chi masses at or below M5 it gives r >= 94, i.e. an effective Y far below the ~1.6 needed. **It cannot replace the friction closure as it stands.** The registered B1 outcome is INCONCLUSIVE, so this negative is descriptive and conditional.
- **delta = 1e-3.** It is still INCONCLUSIVE. The A1 failure modes (x_c = inf, C2 drift, the lapse stops) are diagnosed. The new limit is a breakdown after the light-cone crossing in the friction-limited Y=1 roll. At delta = 1e-3, Y=1 is a 10x stronger friction in Hubble units than at delta = 0.1.
- **Fine tuning.** A1 found that a 1% change in d destroys the radiation era. B4 quantifies this:
  - A1's d* leaves Lambda_res = -3.4e-4 H0^2, which explains the A1 recollapse and the negative Omega_vac.
  - One radiation e-fold needs |dd/d*| <~ 1e-4.
  - A radiation era lasting to BBN needs 7e-7 to 7e-80.
  - This is the cosmological-constant problem.
- **Reduced model.** The 4D EFT cannot stand in for 5D runs when predicting r (B4 NOT VALIDATED).

The files are CLAIM_REVIEW.md, NEXT_TESTS.md, RESULTS_SUMMARY.json, REPRODUCE.md and CRITIC_NOTES.md (final consistency review of the top-level files against the workstream JSON, registrations and audits; it lists the corrections made and the open issues).
