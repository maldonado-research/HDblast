# B4 pre-registration: reduced (shell / 4D EFT) cross-check of A1 and a fine-tuning measure

Registered 1 October 2026, before any reduced-model trajectory was computed and before any comparison of reduced-model numbers with the A1 5D numbers.
Ricardo Maldonado's HDBLAST program. This file is not edited after registration; changes are appended as dated notes with reasons.

## 0. Disclosure (what was seen before registering)

- The full A1 verdict, registration and `A1_RESULTS.json` (all target numbers below), and printed time histories (phi_b, H, v, R, Weyl, vac, rad, kin, fric versus H0 tau) of the A1 runs d*/Y=1, d*/Y=0.3, d*/Y=0, d*/Y=3 (dc = 1e-2, dz = 5e-4) and 1.05 d*/Y=1. From these I noted: the scalar stops at phi_b = 1.0361 with strong damping; the late vacuum term is -3.40e-4 H0^2 at d*; the Weyl term during the roll is largely a bookkeeping offset of the local vacuum term and only becomes a ~a^-4 fluid after the scalar stops.
- Development check (no comparison with A1): the Chat 9 EFT functions I, f, Z, Z_E evaluated at phi = 0 ... 0.9999 (`dev/hj_effective_theory_chat9_copy.py`). They are singular at phi -> 1 (f' ~ (1-phi)^(-8/9), Z ~ (1-phi)^(-17/9)); the canonical field reaches phi = 1 at finite distance.
- Back-of-envelope (not decisive, recorded for honesty): with H^2 = Lambda + C/a^4 and Lambda < 0, H = 0 at tau = pi/(4 sqrt|Lambda|) after the effective origin; inserting the A1 late vacuum values gives numbers of the size of the A1 H_zero times. This informed the choice of the late-time sector below.
- The B3 (5D Y scan) outputs were not read and will not be read.

## 1. Reduced model "E" (labelled: a reduced model of the A1 model change, not a new 5D model)

5D model being reduced: A1's sigma = 2W + delta(1 + c phi + d phi^2/2), delta = 0.1, registered c, kappa5^2 = 1, friction closure kappa5^2 j = Y v, Z2 bulk, one shell.

**Phase I (roll; Chat 9 two-derivative EFT, brane/Jordan frame = shell proper time tau and shell scale factor a):**
S4 = int sqrt(-g)[ f R/2 - Z (d phi)^2/2 - V ] + radiation, with f = 2I(phi), Z = 2W(1 - 2WI/3)/W_phi^2, V = sigma - 2W = delta(1 + c phi + d phi^2/2) (t -> 0 coefficients, no O(delta) recalibration; same d as the 5D run).
Integrated in the Einstein frame (g_E = f g, t_E, a_E = sqrt(f) a): 3 H_E^2 = Z_E phi'^2/2 + V/f^2 + rho_E, Z_E = Z/f + (3/2)(f'/f)^2, rho_E = R/f^2 (R = kappa5^2 rho on the brane);
Z_E(phi'' + 3 H_E phi') + Z_E' phi'^2/2 + (V/f^2)' = -Gamma_E phi', Gamma_E = Y f^(-3/2); d(rho_E a_E^4)/dt_E = Gamma_E phi'^2 a_E^4.
(This reproduces the 5D ledger dR/dtau + 4HR = Y v^2 exactly in the brane frame.) Brane-frame H = sqrt(f)(H_E - f' phi'/(2f)).
Initial data: the EFT stationary point of the seed tension (c + dc, same d) with zero velocity, evolved with c; dc in {1e-2, 1e-4}.

**Matching (Weyl fluid fed by the shell identity):** Phase I runs until phi = phi_m (primary phi_m = 0.95). There the brane-frame H and R are kept continuous; the scalar is replaced by its exact late static value and the identity

H^2 = Lambda_res(d) + rad + Weyl,   rad = sigma_s R/18 + R^2/36 (scaled as in A1),

defines Weyl_m = H_m^2 - Lambda_res(d) - rad_m: all scalar energy not yet in radiation is assigned to the bulk Weyl (dark-radiation) fluid, as the 5D runs do when the shell scalar stops.

**Phase II (exact once v = 0):** R a^4 and Weyl a^4 constant (exact: transport law and ledger with v = 0), H^2 = Lambda_res(d) + rad(R) + Weyl, integrated through H = 0 into recollapse (H' from the same relation).

**Lambda_res(d) and phi_s(d)** (not taken from the 5D runs): the exact static '+1' branch (copy of the M8 solver) for d/d* in {0.90, 0.93, 0.96, 0.98, 0.99, 0.995} (where H^2 > 0), fitted by a cubic in d and extrapolated to d/d* = 1, 1.01, 1.05 (and the tuning window). Extrapolation error = difference between cubic and quadratic fits; the O(delta^2) series slope is a cross-check.

**Units.** H0 = initial Hubble rate of each model's own initial shell. Lambda_res and phase II are in 5D units with the 5D H0^2 = 1/rho_b^2; phase I times are converted with the EFT's own H0 (the EFT and 5D initial H0^2 differ by ~1% at delta = 0.1); this is part of the tested approximation.

**Diagnostic (not part of the verdict):** the exact Weyl transport law with matter (V2 form, as in `a1_analyze.weyl_identity`) integrated along the phase-I EFT trajectory, compared with the identity-defined Weyl at phi_m; it measures how far the EFT trajectory is from a solution of the 5D shell relations.

## 2. Validation rule (fixed now)

Targets (A1, `A1_RESULTS.json`; seed dc = 1e-2 unless stated, finest pair averaged):

| # | Quantity | A1 value | Tolerance |
|---|---|---|---|
| T1 | r at plateau, Y = 0.3 | 0.463 | relative 30% (0.324-0.602) and r > 0.1 |
| T2 | r at plateau, Y = 1 | 0.0734 | relative 30% (0.0514-0.0954) and r <= 0.1 |
| T3 | Omega_vac at plateau, Y = 0.3 | -0.148 | absolute 0.05 |
| T4 | Omega_vac at plateau, Y = 1 | -0.207 | absolute 0.05 |
| T5 | H = 0 time, d = 1.01 d*, Y = 1 | 13.54 | relative 15% |
| T6 | recollapse (H/H0 < -0.05) time, 1.01 d* | 18.85 | relative 15% |
| T7 | H = 0 time, d = 1.05 d*, Y = 1 | 6.93 | relative 15% |
| T8 | recollapse time, 1.05 d* | 8.20 | relative 15% |

The plateau and recollapse are found by applying A1's rule (copied `plateau_index`, H/H0 < -0.05) to the reduced time series.

- **VALIDATED cross-check:** all eight pass with the primary phi_m = 0.95 and both seeds.
- **PARTIAL (late sector validated):** T3-T8 pass, T1 or T2 fails. Then the timing / tuning-window results are reported as conditional-supported and the r(Y) predictions as unvalidated (conditional, not usable as a prediction for B3).
- **NOT VALIDATED:** any of T3-T8 fails. All reduced-model outputs are then conditional and only the exact-algebra fine-tuning statements (section 3) stand.
- Systematics reported (do not change the verdict): phi_m in {0.90, 0.98}; ODE rtol 1e-8 vs 1e-11 (results must agree to 1e-4 relative, otherwise the numerics are fixed first and a dated note is added).

## 3. Predictions and fine-tuning measure (computed whatever the verdict; labelled by it)

1. r(Y) and Omega_r at plateau for Y in {0.3, 0.5, 0.7, 1, 1.5, 2, 3, 4, 5} at d* and both seeds.
2. Tuning window: the set of d/d* for which the A1 per-run rule gives a radiation era (plateau with Omega_r >= 0.9 before any recollapse), and separately for which >= 1 e-fold with Omega_r >= 0.9 and H > 0 occurs; width versus Y.
3. Fine tuning: |Lambda_res|/H0^2 and |d - d*_exact|/|d*| needed for (a) one radiation-dominated e-fold, (b) a radiation era until BBN (T = 1 MeV) for reheating temperatures 10 MeV ... 1e16 GeV, with the exact series slope dLambda/dd = delta/54 + O(delta^2) for general delta, compared with the observed vacuum (~1e-122 in Planck units, `LITERATURE_2022_2026.md`). Stated as a restatement of the cosmological-constant problem, not a solution.
4. Why Omega_vac < 0 at d* and recollapse estimate tau_rec = pi/(2 sqrt|Lambda_res|) (closed form for H^2 = Lambda + C/a^4, exact algebra), plus the reduced-model numbers.

## 4. Controls

- C1 (calibration): EFT hilltop growth rate for the registered model (d = 0, Y = 0, t -> 0) against the known 5D value 1.65719 (delta = 1e-3, Chat 14) and the closed form mu^2 = -7.7198; tolerance 2e-3.
- C2: the plus-branch solver at d = 0 reproduces the M2 value H_+^2/H0^2 = 0.406168, phi = 0.991435 (tolerance 1e-5).
- C3 (wrong formula): Gamma_E = Y (missing f^(-3/2)) and the Jordan Friedmann equation without the f-dot term; they must change T2 by more than its tolerance, otherwise the T1/T2 comparison is declared insensitive.
- C4 (informative, after the verdict): Lambda_res(d) and phi_s(d) compared with the A1 late values.

---

### Dated note, 1 October 2026 (after the first validation run; model, targets, tolerances and rule unchanged)

1. **Seed for T5-T8.** The A1 detuned controls (1.01 d\*, 1.05 d\*) exist only for the dc = 1e-2 seed, and section 2 says "seed dc = 1e-2 unless stated". The first driver also compared the dc = 1e-4 reduced runs with these dc = 1e-2 targets, which is a like-for-like error (the seed shifts all times by about 2.2 H0^-1). T5-T8 are now evaluated with dc = 1e-2 only; T1-T4 still use both seeds. This changes T5/T7 from FAIL to PASS but not the verdict, which is fixed by T1-T4 (FAIL in both versions).
2. **Post-registration variants and diagnostics (do not change the verdict).** After the registered model failed T1-T4, I added (a) an Einstein-frame energy-dump matching, (b) a radiation-only late sector (Weyl set to zero at the match), (c) the registered matching recipe applied to the A1 5D trajectory itself, and (d) a comparison of the phase-I EFT trajectory with the A1 5D trajectory at equal phi_b. All are labelled post-registration in `B4_VALIDATION.json` and the README.
3. **Static-branch fit range.** The exact '+1' solver does not converge for d/d\* > 0.98 (the vertex offset underflows below 1e-19), so the cubic fit uses the converged points d/d\* in [0.90, 0.98] (9 points) instead of the planned set up to 0.995. Quadratic and cubic fits differ by < 5e-7 H0^2 at d\*.

Correction to item 3 of the note above (same day): the fit uses 10 converged points (d/d\* = 0.90, 0.91, ..., 0.97, 0.975, 0.98), and the quadratic-cubic difference at d\* is 5.4e-7 H0^2 (4.5e-6 H0^2 at 1.05 d\*), not "< 5e-7". Values from `B4_STATIC_LAMBDA.json`.
