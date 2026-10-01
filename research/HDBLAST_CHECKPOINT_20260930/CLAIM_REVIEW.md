# Claim review: HDBLAST checkpoint 30 September 2026

Rule: the audit verdict wins over the producer label. A claim that was "partially-confirmed" keeps its label only where the confirmed part supports it, and the stated correction is applied. A claim that was "unverifiable" is not upgraded. The full reasoning is in each subfolder's AUDIT.md.

## Cross-workstream points resolved here

- **B1 versus B3: effective Y.** B3 needs Y ~ 1.6 or more for r <= 0.03. B1's derived chi channel at lambda_c = 1 gives r >= 94, far worse than A1's Y = 0.3 (r = 0.46). So the derived source is far below the required effective Y in the scanned (phi*, G, y). This is inferred from labelled results of both workstreams, and is not a new computation. Label: negative (conditional/descriptive).
- **B4 versus A1/B3: Omega_vac < 0.** The negative residual vacuum that ends every B3 plateau is Lambda_res(d*_M8) = -3.4e-4 H0^2 from B4 (numerical, reproduced). This is consistent with B3's "about -3.5e-4 H0^2".
- **The B4 r(Y) predictions** are not usable against B3 (NOT VALIDATED). No comparison is made.

## B1: derived chi source (MODEL CHANGE B1)

| Claim | Producer | Audit | Final |
|---|---|---|---|
| Closed shell system with chi source, 22 identities, 8 controls | exact-verified | confirmed | exact-verified |
| K1 / K2 calibrations | numerical | confirmed | numerical |
| Registered aggregate INCONCLUSIVE at lambda_c = 1 | inconclusive | confirmed (longest A1-tolerance window Delta ln a = 0.394 < 0.5) | inconclusive |
| b_cut r = 94-8.5e3 (y=1), 3e2-1.1e5 (y=0.1), Omega_r <= 0.012 | numerical | partially confirmed | numerical (descriptive). Corrections: phi*=0.9, G=100 not fully decayed (Omega_chi/Omega_r 0.095, R a^4 still rising 16%; its r is ~1.9x high from underproduction); the y=0.1 values are snapshots, not lower bounds (r may fall to ~50); constraints <= 1.2e-4 hold on the fine grid only (coarse 4.2e-4) |
| Reduced estimate FAIL-Weyl / recollapse at b_cut | conditional | confirmed (end values; transient R ratio 0.74-0.92) | conditional |
| lambda_c needed 11-50 (r<=0.1) | conditional | confirmed as unreliable | conditional (not a threshold) |
| Exploratory b = 1e-3: r = 0.282 | numerical | confirmed (not a plateau, window 0.251) | numerical, exploratory and unregistered |
| Chi cannot replace friction below M5 | negative | partially confirmed | **negative (conditional/descriptive)**; the registered outcome is INCONCLUSIVE |
| b_need unreliable; controls change r to -338, -1149, 22.5 | inconclusive | partially confirmed: values mis-paired | inconclusive. Corrected: at the run end the controls give 0.316, 0.082, 0.169 vs 0.170 (~2x sensitivity; noD no effect) |
| K3 consistency checks | numerical | confirmed, exceptions disclosed | numerical |

Registration: respected as far as the evidence allows. The timing of the second note relative to the grid start is unverifiable.

## B2: delta = 1e-3

| Claim | Producer | Audit | Final |
|---|---|---|---|
| kappa=0 pre-run unstable; kappa=10 gives x_c = 11.3 | numerical | confirmed (dc=1e-4 value 18.7 single resolution) | numerical |
| A1 stops at 3.1 / 11.0 were numerical | numerical | partially confirmed | numerical for loss of reliability. The dz=2e-4 end (10.97) is near the physical crossing (11.36), so its end is not clearly numerical |
| Registered C2 fails | negative | confirmed | negative |
| Seed-offset explanation; post-hoc estimators give 1.65720 +- 3e-5 | numerical | partially confirmed; independent derivative estimator gives 1.657184-1.657187 | numerical (post hoc). Amendment 2 criteria (ii) monotone=false and (iii) C2d fail as written, so the transient interpretation is not formally accepted |
| Registered tuned test INCONCLUSIVE | inconclusive | confirmed (bit-identical rerun) | inconclusive. C5 is weak here: the 3H control discriminates only at the median, and R-dropped is indistinguishable |
| Breakdown after light-cone crossing in every Y=1 run; x_c ruled out | numerical | partially confirmed; two parts refuted | numerical (diagnosis), restricted to the dc=1e-2 runs and the sensitivity set. The registered dc=1e-4 pair ends BEFORE F_inf (x_c=11.8 artifact), and C7 ends at T=6.79. **x_c is NOT ruled out** (11.3->12.3 moves the end 12.47->13.22) |
| +0.65 H0tau per doubling; refinement cannot fix | numerical | partially confirmed | end times numerical; "cannot fix" is **conditional** (3-point extrapolation; fixed-T residual converges ~10x per doubling) |
| Y=0 baseline: solver works for fast roll | numerical | partially confirmed | numerical at S2 only; S1 Y=0 runs lose reliability at or before their crossing |
| C7 kappa=0 main run unusable | numerical | confirmed (residual at stop is O(1), not 1e-2) | numerical |

Registration: largely respected. There are five dated notes, not four. S3' was launched about a minute before its note, and D6 before note 3 (both diagnostic only). The Amendment 2 criteria were reinterpreted after the fact (no practical effect).

## B3: r(Y) scan (friction closure, delta = 0.1, d = d*)

| Claim | Producer | Audit | Final |
|---|---|---|---|
| A1 Y=3 no-plateau diagnosis | numerical | partially confirmed | numerical for the 4th-order drift and the 0.57-0.64 e-fold window. The amplification by the lapse is a correlation only. D1 conclusion strings and the Omega values in section 0 are wrong (D1 gives 1.09-1.12 / -0.10 to -0.13) |
| Y=5 cfl 0.5 stiffness | numerical | confirmed (grows ~Y only for Y>=3) | numerical |
| Y=0.5 FAIL-Weyl 0.2237 | numerical | confirmed | numerical |
| Y=0.7 FAIL-Weyl 0.1326 | numerical | confirmed | numerical |
| Y=1.5 PASS-conservative 0.0360 | numerical | confirmed | numerical |
| **Y=2 PASS-combined both seeds, r = 0.0213** | numerical | confirmed (4 levels, 4th order, bit-identical rerun) | **numerical** |
| Y=3 PASS-combined r = 0.0101 | conditional | partially confirmed (value robust; class rests on post-hoc level; the tertiary order is measured near turnaround) | conditional |
| Y=5 UNRELIABLE, r = 0.0044 | inconclusive | confirmed | inconclusive |
| Power law Y^-1.67, Y(r=0.03) ~ 1.64 | numerical | partially confirmed (p = 1.63 without Y=3; "mechanism" restates the ledger identity) | numerical (empirical fit only) |
| e-folds to turnaround 0.50->0.09 | numerical | partially confirmed | numerical: measured 0.383 (Y=1.5), 0.319 (Y=2, dc=1e-2), 0.227 (Y=3), 0.093 (Y=5, unreliable). The values for Y=0.5, Y=0.7 and Y=2 dc=1e-4 are estimates only (runs stopped before turnaround) |
| Strengthened C5 controls | numerical | partially confirmed (all-record R-dropped also fails at Y=2 dc=1e-4, 8.6x) | numerical; window chosen after seeing A1 |
| C10 / C11 | numerical | confirmed | numerical |
| C6 ledger | numerical | confirmed (Y=0.5 marginal fail) | numerical |

Registration: respected. Minor undisclosed items: two trigger scripts were modified after the runs (immaterial), and the dev check of b3_analyze.py cannot be traced.

## B4: EFT cross-check and fine tuning

| Claim | Producer | Audit | Final |
|---|---|---|---|
| C1 growth rate 1.65750 vs 1.65719 | numerical | confirmed (delta-independent check) | numerical |
| C2 static solver reproduces M2 | numerical | confirmed | numerical |
| Lambda_res(d*) = -3.399e-4; d*_exact = -3.10415 | numerical | confirmed (12 alternative fits) | numerical |
| Registered reduced model NOT VALIDATED | negative | confirmed | negative |
| Controls/systematics | numerical | partially confirmed (rtol pair differs from registered, no note) | numerical |
| Failure is in matching | numerical (post-reg.) | partially confirmed | numerical (post-registration). Corrected: W a^4 falls 33.5x (Y=1) / 11.2x (Y=0.3), not ~15x; EFT tracks 5D to ~1-5% up to phi 0.7; v is 37% low at the match |
| Post-registration variants | conditional | confirmed | conditional |
| r(Y) prediction unusable | inconclusive | confirmed | inconclusive |
| Tuning window in d (1-2.3%) | conditional | partially confirmed | conditional, and these are post-registration radiation-only windows. The registered model's windows [1.015,1.0225], [1.0125,1.0275], [1.015,>=1.03] exclude d* |
| Lambda<0 radiation solution | exact-verified | confirmed (independent derivation) | exact-verified |
| N_RD ~ 0.5; H=0 ~44, crunch ~87 | numerical | partially confirmed | numerical for N_RD (0.48 direct count). The times are **conditional** (closed-form extrapolation past the A1 run end 28.8). "Omega_vac -0.011 at stop" is wrong: A1's value is -0.160 |
| One-e-fold tolerance 1.3e-4 | numerical | confirmed | numerical |
| BBN tolerance 7e-7 ... 7e-80 | conditional | confirmed | conditional (relative to the exact zero; the d*_exact fit spread ~1e-6 exceeds the 10 MeV tolerance) |
| No first-attempt files to reuse | - | unverifiable | unverifiable |

Registration: largely respected. Undisclosed: the rtol pair, the added |Omega_vac| <= 0.1 e-fold condition, and headlining the post-registration window.

## Final critic pass (1 October 2026)

A last consistency pass compared the top-level files with the workstream JSON outputs, registrations and audits (details and open issues in CRITIC_NOTES.md). Corrections applied to the top-level files, none of which changes a registered verdict:
- B2 C2 post-hoc values: the quoted 1.65718-1.65722 holds on the rho-grid runs only; on C2d both the mode fit and the derivative estimator give 1.65705-1.65709 (about 1.4e-4 low). "Recover the target to about 3e-5" was removed.
- B3 Y=2: the pass rests on the pre-registered secondary pair; the 1e-3 level has no plateau; Omega_r >= 0.9 at the plateau is helped by Omega_vac = -0.39. B3 Y=3: noted that the outcome was known before registration.
- B3 e-folds after the plateau: stated per Y, with the Y=0.5/0.7 values marked as estimates.
- B4: "r wrong by about 12x" is Y=1 only (3.6x at Y=0.3); "the error is not from the 4D roll" softened (v is 37% low at the match).
- B1: the r >= 94 values are run-end values, not plateau values; the corrected b_need control values come from the B1 audit JSON and lie beyond the runs' reliable ends.
- RESULTS_SUMMARY.json: the missing B2.2 entry was added, and scope/caveat fields were added.

## Global caveats

All results are for a toy model with phenomenological pieces: the friction closure, a fine-tuned d, and N_eff thresholds carried over from A1 at snippet level. Nothing is observational evidence. Nothing is peer reviewed. No claim of discovery or proof is made.
