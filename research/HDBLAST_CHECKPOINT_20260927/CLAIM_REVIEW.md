# Claim-by-claim review — HDBLAST checkpoint, 27 September 2026

Prepared with AI assistance for Ricardo Maldonado's HDBLAST program. This is an internal review, not external peer review. The table is generated from `RESULTS_SUMMARY.json` by `synthesis/build_claim_review.py`; the critic corrections (last section) are then applied by `critic_corrections.py`.

## How to read this table

- **Producer label**: what the workstream itself called the result.
- **Verifier verdict**: the independent auditor's verdict (each workstream folder has a `VERIFICATION.md` and a `verify/` folder with the auditor's own scripts). *not-independently-audited* means only the workstream's own checks exist (literature, program map, and the synthesis checks done for this report).
- **Status here**: the label this checkpoint uses after applying the verifier's corrections. *corrected* means the original claim was refuted and the corrected statement is given instead.
- **Established** (✔) only when the verifier confirmed the claim and it is not inconclusive or corrected.

Status vocabulary: exact-verified (symbolic identity), numerical (floating point with measured tolerances, never interval-certified here), conditional (valid under stated assumptions), negative (the idea fails the test), inconclusive, corrected, context, procedural-negative.

Totals: 66 claims. Verdicts: partially-confirmed 13, confirmed 39, refuted 3, not-independently-audited 11. Status: numerical 27, exact-verified 11, conditional 10, negative 6, corrected 3, context 4, inconclusive 4, procedural-negative 1.

## Refuted or corrected claims (read these first)

- **AS-4** — original: *The delta^9 shell coefficients depend on cone regularity and need global matching.* → **Correction:** The delta^9 coefficients are local and exactly computable (resonance constant drops out): at registered c, eta_b: 3.249147254269e-6, H^2: -2.759694789e-7; they match the workstream's fits and an independent BVP.
- **PH-8** — original: *Perturbative decay chi->psi psibar gives a radiation/vacuum ratio <=1.3e-4.* → **Correction:** Stale README number: the JSON maximum is 2.34e-4 at (phi*=0.5, G=300, y=1). Conclusion (<<1) unchanged.
- **MX-10** — original: *The registered linear tension cannot tune away the vacuum: the initial-shell family ends near c->0+ and the shell cannot be built at c*.* → **Correction:** The stop was a parametrisation artefact (phi_h=-1+10^x). With the sign allowed, a static shell exists all the way to c*=-0.99307 at delta=0.1 (phi_h=-1.00695, phi_b=-1.7287), and the c* pilot (single grid, one seed sign, Y=0, delta=0.1 only) rolls from phi_b=-1.73 across phi=-1 to 0.80. It is a qualitatively different initial state; its endpoint is open.

## Other corrections from the audits (partially confirmed)

- **GI-1**: Independent-polish agreement in rho_b is <=8.5e-14 relative, not 6e-15 as first quoted.
- **GI-8**: Label split: mu^2>=0 (Sturm-Liouville) is exact; absence of bound states in (0,9/4) is numerical.
- **GI-9**: C1 pass criterion was chosen after seeing an O(delta) deviation; C2's root (-28555) lies outside the main scan, so it is a weak demonstration of detection power.
- **TD-2**: True for the fully discrete RK4 system and for shell observables. The semi-discrete operator also has grid-scale shell-localised modes with Re(lambda)~+0.245 that do not shrink with refinement (numerical artefacts damped by RK4).
- **TD-4**: Supported precision is -1.500 +/- 1e-3 (Richardson depends on assumed order). The 'shell-localised' second profile is actually a bulk bump (near-shell weight 2.4e-10); the time-domain test barely excites the shell (~1e-8), so the spectral evidence carries the no-bound-state claim.
- **PH-2**: G x (relative deviation) is constant to 0.4%/0.1% (not '4 digits'); D at phi*=0.25 and 0.9 are closed-form values (measured 10.6476 and 1.3415).
- **EC-2**: The geometric gain applies to the un-normalised H+2M; the background-normalised H has an extra factor ~360. The t=0.05 density diagnostic depends on the time step.
- **EC-6**: Measured factors are ~25-200 at t=0.5-1 and ~8 at t=0.25; at t=1 'A+B smallest' is a tie with B.
- **MX-2**: The N_eff input values come from search snippets and are unverifiable here.
- **MX-3**: Algebra exact, but the exclusion rests on a cited stability theorem and the Z2 topology: conditional, not exact. The linear Lambda_4 formula holds only near lambda=4/3.
- **MX-4**: Covers only ~1.05-1.10 e-folds of resolved contraction; negative only on that range.
- **MX-11**: The runs stopped at the solver's phi_b>0.995 cutoff, not the chart freeze. The H^2 budget has large cancellations (vacuum -117%, Weyl +121%, radiation +69%). With the cutoff lifted, Weyl/radiation is 1.23 at the end of the reliable window (about 12x the Delta N_eff-derived limit of 0.1). Values after the chart freeze are not usable: proper time advances by only ~0.2 H0 tau while the ratio falls through 0.6 and below zero and then swings back as H turns negative (gauge artefact or reversal); the earlier '~0.6 after the freeze' figure is withdrawn (CRITIC_CHECKS.json, post_freeze_dstar).

## Stability of the +1 branch, Method A (gauge-invariant perturbation theory) (`stability_gauge_invariant/`)

| ID | Claim | Producer label | Verifier verdict | Status here | Est. | Evidence | Correction / note |
|---|---|---|---|---|---|---|---|
| GI-1 | The +1 static background is reproduced to all stated digits (phi_b=0.9999159473169134, rho_b=129.92476284968, H^2=5.924014794329e-5, H/H0=0.606721732). | numerical | partially-confirmed | numerical |  | `stability_gauge_invariant/BACKGROUNDS.json`<br>`stability_gauge_invariant/VERIFICATION.md` | Independent-polish agreement in rho_b is <=8.5e-14 relative, not 6e-15 as first quoted. |
| GI-2 | Gauge-invariant linearized system (X, Z=Psi/phi') and shell condition M=BX+3(mu^2+4)Z/rho^2=0 with B=phi''/phi'+sigma''/2; no brane bending in longitudinal gauge; tensor h''+4Hh'+mu^2h/rho^2=0, h'(y_b)=0. | exact-verified | confirmed | exact-verified | ✔ | `stability_gauge_invariant/DERIVATION_RESULTS.json`<br>`stability_gauge_invariant/verify/indep_derivation.json` | 55/55 symbolic checks, 4/4 wrong-formula controls; auditor's independent derivation 26/26. |
| GI-3 | Calibration on the original unstable shell: one bound state mu^2=-7.717871625260, growth 1.657193631259 (recorded Chat 9: -7.717871625176). | numerical | confirmed | numerical | ✔ | `stability_gauge_invariant/SPECTRUM_RESULTS.json`<br>`stability_gauge_invariant/WINDING_RESULTS.json` | Validation against a known internal target, not a blind test. |
| GI-4 | On the +1 branch at delta in {0.0003,0.001,0.003,0.01,0.03,0.1} there is no scalar eigenvalue with mu^2<9/4 (min normalized mismatch >=0.99274). | numerical | confirmed | numerical | ✔ | `stability_gauge_invariant/SPECTRUM_RESULTS.json`<br>`stability_gauge_invariant/verify` | Real mu^2 scanned on a 360-point grid in [-400, 2.2499]; mu^2<-400 and complex mu^2 outside the winding contours are excluded only by the conditional B>0 theorem (GI-6). Regularity (normalisability) at the cone is an assumed criterion. Auditor's margin test: an instability needs B in [-0.659, 0.0081], actual B=3.479-3.555. [critic 28 Sept 2026] |
| GI-5 | No complex scalar eigenvalues inside [-60,2]x[-30,30] and [-400,2]x[-200,200] (winding 0; 1 on calibration). | numerical | partially-confirmed | numerical |  | `stability_gauge_invariant/WINDING_RESULTS.json` | Auditor re-ran only a subset (small contour, 3 cases). |
| GI-6 | If B>0 at the shell, all scalar eigenvalues are real and mu^2>-4; B=3.479-3.555 on the +1 branch and -2.68e-4 on the original shell. | conditional | confirmed | conditional | ✔ | `stability_gauge_invariant/DERIVATION_RESULTS.json`<br>`stability_gauge_invariant/verify/identities.json` |  |
| GI-7 | No physical l=1 mode (B phi'_b != 0) and no l=0 static zero mode (junction Jacobian condition number 15-17). | numerical | confirmed | numerical | ✔ | `stability_gauge_invariant/SPECTRUM_RESULTS.json` |  |
| GI-8 | Tensor sector: only the normalizable massless graviton; no KK bound state in (0,9/4); no tachyons. | exact-verified | partially-confirmed | numerical |  | `stability_gauge_invariant/SPECTRUM_RESULTS.json`<br>`stability_gauge_invariant/VERIFICATION.md` | Label split: mu^2>=0 (Sturm-Liouville) is exact; absence of bound states in (0,9/4) is numerical. |
| GI-9 | Controls C1-C5 behave as expected (Gegenbauer limit, flipped sigma'' creates a tachyon, dropping gravity removes the calibration root, 4D limit, perturbed c). | numerical | partially-confirmed | numerical |  | `stability_gauge_invariant/CONTROLS_RESULTS.json` | C1 pass criterion was chosen after seeing an O(delta) deviation; C2's root (-28555) lies outside the main scan, so it is a weak demonstration of detection power. |
| GI-10 | The static +1 branch is linearly stable in the scalar (including brane bending and both junctions) and tensor sectors at all six tested delta; the only marginal mode is the massless graviton. | numerical | confirmed | numerical | ✔ | `stability_gauge_invariant/README.md`<br>`stability_gauge_invariant/VERIFICATION.md` | Floating point, not interval-certified; vector sector, nonlinear stability, tunnelling and dynamical attraction not studied. |

## Stability of the +1 branch, Method B (time domain and discrete operator) (`stability_time_domain/`)

| ID | Claim | Producer label | Verifier verdict | Status here | Est. | Evidence | Correction / note |
|---|---|---|---|---|---|---|---|
| TD-1 | Time-domain calibration on the original shell: rate 1.657193368 (linear) / 1.657193408 (nonlinear) vs 1.6571936312. | numerical | confirmed | numerical | ✔ | `stability_time_domain/results/KEY_RESULTS.json`<br>`stability_time_domain/VERIFICATION.md` |  |
| TD-2 | Only gauge (lambda->1) and far-boundary (~0.6/L) modes have positive rates on the +1 branch. | numerical | partially-confirmed | numerical |  | `stability_time_domain/results/KEY_RESULTS.json`<br>`stability_time_domain/VERIFICATION.md` | True for the fully discrete RK4 system and for shell observables. The semi-discrete operator also has grid-scale shell-localised modes with Re(lambda)~+0.245 that do not shrink with refinement (numerical artefacts damped by RK4). |
| TD-3 | Every shell-supported mode lies on Re(lambda)=-3/2 (dS continuum mu^2>9/4) to ~1e-8; no scalar bound state with mu^2<9/4 in the homogeneous sector. | numerical | confirmed | numerical | ✔ | `stability_time_domain/results/KEY_RESULTS.json`<br>`stability_time_domain/verify` | Homogeneous (FRW-symmetric) sector only, at delta=0.001 with 0.003 and 0.01 as perturbed-parameter checks. The time-domain signal barely reaches the shell (peak /f_b/~8.5e-9), so this claim rests on the spectral evidence. [critic 28 Sept 2026] |
| TD-4 | Leading physical rate on the +1 branch converges to -3/2 (in H_+ units). | numerical | partially-confirmed | numerical |  | `stability_time_domain/results/KEY_RESULTS.json` | Supported precision is -1.500 +/- 1e-3 (Richardson depends on assumed order). The 'shell-localised' second profile is actually a bulk bump (near-shell weight 2.4e-10); the time-domain test barely excites the shell (~1e-8), so the spectral evidence carries the no-bound-state claim. |
| TD-5 | Flipping sigma'' in the scalar junction produces lambda=167.48925 (partner sums to -3). | numerical | confirmed | numerical | ✔ | `stability_time_domain/results/KEY_RESULTS.json` |  |
| TD-6 | Whole-domain constraints do not converge in the far stretched region (under-resolved outgoing radiation); near-shell constraints converge. | negative | confirmed | negative | ✔ | `stability_time_domain/results/KEY_RESULTS.json` |  |
| TD-7 | Model equations, junctions and normalisations re-derived (32/32 sympy checks). | exact-verified | confirmed | exact-verified | ✔ | `stability_time_domain/verify` |  |

## Exact analytic structure of the +1 branch (`analytic_structure/`)

| ID | Claim | Producer label | Verifier verdict | Status here | Est. | Evidence | Correction / note |
|---|---|---|---|---|---|---|---|
| AS-1 | The 22 Sept O(delta) and O(delta^2) coefficients are reproduced exactly. | exact-verified | confirmed | exact-verified | ✔ | `analytic_structure/SERIES_COEFFICIENTS.json`<br>`analytic_structure/verify/V1_SERIES_INDEPENDENT.json` |  |
| AS-2 | phi_b, H^2, rho_b, y_b and eta_h expansions exact in c through O(delta^8); e.g. h3=-c^2(11767c+11362)/2826240. | exact-verified | confirmed | exact-verified | ✔ | `analytic_structure/SERIES_COEFFICIENTS.json`<br>`analytic_structure/verify/V1_SERIES_INDEPENDENT.json` |  |
| AS-3 | 30-65 digit solutions of the nonlinear boundary-value problem confirm the new coefficients (e.g. H^2(delta=1e-3)=5.9240147943288795996237552780e-5). | numerical | confirmed | numerical | ✔ | `analytic_structure/SERIES_VS_BVP.json`<br>`analytic_structure/runs/BVP_SCAN_reg.json`<br>`analytic_structure/verify/V4_COMPARE.json` |  |
| AS-4 | The delta^9 shell coefficients depend on cone regularity and need global matching. | numerical | refuted | corrected |  | `analytic_structure/VERIFICATION.md`<br>`analytic_structure/verify/V1_SERIES_INDEPENDENT.json` | The delta^9 coefficients are local and exactly computable (resonance constant drops out): at registered c, eta_b: 3.249147254269e-6, H^2: -2.759694789e-7; they match the workstream's fits and an independent BVP. |
| AS-5 | Exact resonance at (m,j)=(2,7) with S=-24576/17; no log in shell observables at delta^9; the cone value eta_h gets a delta^16 ln(delta) term. | numerical | confirmed | exact-verified | ✔ | `analytic_structure/RESONANCE_PROBE.json`<br>`analytic_structure/verify/V1_SERIES_INDEPENDENT.json` | Verifier upgraded 'no log at delta^9' from fitted to exact. |
| AS-6 | Leading cone displacement eta_h=-(111537/131072)c(1+c)^7 delta^8[1+O(delta)]. | exact-verified | confirmed | exact-verified | ✔ | `analytic_structure/SERIES_COEFFICIENTS.json`<br>`analytic_structure/verify/V4_COMPARE.json` |  |
| AS-7 | Degree-14 Gegenbauer structure: n(n+4)=252, n=Delta-4 with Delta=3W''/W=18 at phi=+1; closed forms; no polynomial at phi=-1. | exact-verified | confirmed | exact-verified | ✔ | `analytic_structure/GEGENBAUER_STRUCTURE.json`<br>`analytic_structure/verify/V3_GEGENBAUER_SPECTRUM.json` | Wording: at phi=-1, s=3W''/W=-18/5 equals 4-Delta, not Delta. |
| AS-8 | Leading-order tensor quantization (mu-3/2)P^{-mu}_{1/2}(cosh u_b)=0: massless graviton plus continuum from m=3H/2 (the gap itself is known literature). | exact-verified | confirmed | exact-verified | ✔ | `analytic_structure/SPECTRUM_LEADING_ORDER.json`<br>`analytic_structure/verify/V3_GEGENBAUER_SPECTRUM.json` | Leading order in delta only. |
| AS-9 | Decoupled bulk scalar: no zero mode or tachyon for any H; no mode below 9H^2/4 for H/k<16.65. | conditional | confirmed | conditional | ✔ | `analytic_structure/SPECTRUM_LEADING_ORDER.json` | Neglects scalar-metric mixing; superseded for the actual branch by GI-4/GI-10, which include mixing. |
| AS-10 | Stability of the +1 branch is not established by this workstream. | negative | confirmed | context | ✔ | `analytic_structure/README.md` | Resolved within scope by the two stability workstreams (GI-10, TD-3). |

## Instant preheating along the archived roll-off (`preheating/`)

| ID | Claim | Producer label | Verifier verdict | Status here | Est. | Evidence | Correction / note |
|---|---|---|---|---|---|---|---|
| PH-1 | Flat linear-crossing calibration reproduces the instant-preheating spectrum to <=8.4e-7. | numerical | confirmed | numerical | ✔ | `preheating/CONTROLS.json`<br>`preheating/verify/INDEP_MODES.json` | The two wrong-formula controls are fixed by Gaussian algebra; weak as solver tests. |
| PH-2 | On the archived trajectory N=q^{3/2}/(8pi^3)(1+D/q+...), with D=10.648, 7.0034, 3.0032, 1.3416 at phi*=0.25,0.5,0.75,0.9. | numerical | partially-confirmed | numerical |  | `preheating/CONTROLS.json` | G x (relative deviation) is constant to 0.4%/0.1% (not '4 digits'); D at phi*=0.25 and 0.9 are closed-form values (measured 10.6476 and 1.3415). |
| PH-3 | Closed-form O(1/q) coefficient D for a curved, expanding crossing matches toy backgrounds to ~3e-5. | numerical | confirmed | numerical | ✔ | `preheating/CORRECTION_FIT.json`<br>`preheating/verify/DDP_EXPONENT.json` | Derivation note corrected: exponent-only h^2 term is 15/(2pi), so c_hh and c_hA follow from the mass shift plus exponent; only the +/-pi/8 pieces remain fitted. No novelty claimed. |
| PH-4 | At the M5 cutoff the produced energy is <=1.7e-4 of the residual-vacuum radiation threshold for screen-passing couplings (<=2.8e-3 over the whole scan). | conditional | confirmed | conditional | ✔ | `preheating/SUMMARY.json`<br>`preheating/verify/CHECK_CLAIMS.json` | Evaluated at the data end; the maximum over the post-window profile is 5.9e-4. Conditional on the fixed archived trajectory, the M5 cutoff and instantaneous full conversion. |
| PH-5 | Reaching the threshold needs chi masses >=18 M5 (>=6.3 M5 if only post-crossing masses are constrained); stronger coupling makes it worse (R ~ K/sqrt(G)). | conditional | confirmed | conditional | ✔ | `preheating/SUMMARY.json` |  |
| PH-6 | Wherever the channel could matter (R~1), backreaction on the scalar junction is 10-79%, so the fixed-trajectory calculation stops being self-consistent. | conditional | confirmed | conditional | ✔ | `preheating/SUMMARY.json` |  |
| PH-7 | Unsubtracted one-loop tension shift is 1.1-3.7 x /sigma'/ at the favoured points: keeping the registered sigma requires tuned counterterms. | conditional | confirmed | conditional | ✔ | `preheating/PREHEATING_RESULTS.json` |  |
| PH-8 | Perturbative decay chi->psi psibar gives a radiation/vacuum ratio <=1.3e-4. | conditional | refuted | corrected |  | `preheating/SUMMARY.json`<br>`preheating/verify/CHECK_CLAIMS.json` | Stale README number: the JSON maximum is 2.34e-4 at (phi*=0.5, G=300, y=1). Conclusion (<<1) unchanged. |
| PH-9 | Instant preheating on the archived +1 roll-off does not produce a radiation era in any scanned coupling range. | conditional | confirmed | negative | ✔ | `preheating/README.md`<br>`preheating/VERIFICATION.md` | Conditional negative; the bulk-coupled feedback (all three junctions) is not solved. |

## Constraint control in 5D evolutions (numerical method) (`evolution_constraints/`)

| ID | Claim | Producer label | Verifier verdict | Status here | Est. | Evidence | Correction / note |
|---|---|---|---|---|---|---|---|
| EC-1 | The 22 Sept constraint stall is reproduced (finest-pair order 0.24 at t=0.5, 0.27 at t=1). | numerical | confirmed | numerical | ✔ | `evolution_constraints/ANALYSIS.json`<br>`evolution_constraints/verify/compare_reruns.json` |  |
| EC-2 | The error front is seeded within ~0.02 of the shell for t<~0.05 and carried outward by the exact transport law of e^{3A}(H+2M). | numerical | partially-confirmed | numerical |  | `evolution_constraints/ANALYSIS.json` | The geometric gain applies to the un-normalised H+2M; the background-normalised H has an extra factor ~360. The t=0.05 density diagnostic depends on the time step. |
| EC-3 | The stall comes from an initial-data defect floor near the shell that does not decrease with h. | numerical | confirmed | numerical | ✔ | `evolution_constraints/ANALYSIS.json`<br>`evolution_constraints/verify/recompute_saved.json` |  |
| EC-4 | Discrete-constraint projection of the initial data restores convergence (orders ~2.9-3.1 at the finest pair; H_max(t=0.5) 6.38e-3 -> 9.35e-4). | numerical | confirmed | numerical | ✔ | `evolution_constraints/ANALYSIS.json` | The projection correction itself does not converge below h=1e-4 (~1e-9). |
| EC-5 | Damped transport law (d_t-d_z)[e^{3A}C+]=-kappa e^{3A}C+ for the outgoing constraint source. | exact-verified | confirmed | exact-verified | ✔ | `evolution_constraints/CONTROLS.json`<br>`evolution_constraints/verify/sympy_independent.json` |  |
| EC-6 | Outgoing-characteristic damping (kappa=10, off near the shell) lowers H by 13-100x. | numerical | partially-confirmed | numerical |  | `evolution_constraints/ANALYSIS.json` | Measured factors are ~25-200 at t=0.5-1 and ~8 at t=0.25; at t=1 'A+B smallest' is a tie with B. |
| EC-7 | Asymptotic fourth order for the eps=0.01 seed is not demonstrated (phi_t and B_t converge at ~2). | inconclusive | confirmed | inconclusive |  | `evolution_constraints/ANALYSIS.json`<br>`evolution_constraints/verify/recompute_saved.json` | Verifier localised the deficit on the shell-seeded front z~-t; no constraint remedy touches it. |
| EC-8 | Uniform damping, H-only damping, KO dissipation, tighter inversion and order 6 without projection do not help or are unstable. | negative | confirmed | negative | ✔ | `evolution_constraints/ANALYSIS.json` |  |

## Screen of mechanisms for a radiation era (`mechanisms/`)

| ID | Claim | Producer label | Verifier verdict | Status here | Est. | Evidence | Correction / note |
|---|---|---|---|---|---|---|---|
| MX-1 | The blast produces a positive bulk Weyl (dark-radiation-like) term, peaking at 0.0806-0.0809 H0^2 (~15% of H^2) at H0 tau~6.24; late value not converged. | numerical | confirmed | numerical | ✔ | `mechanisms/M1_WEYL_ALONG_CHAT14.json` | 3 independent grids (not 4). |
| MX-2 | Bulk Weyl radiation cannot be the hot Big Bang and is capped by N_eff at a few percent of the radiation. | negative | partially-confirmed | negative |  | `mechanisms/M5_OBSERVATIONAL_SCREEN.json` | The N_eff input values come from search snippets and are unverifiable here. |
| MX-3 | Dark-bubble cosmology is not available in the registered model. | negative | partially-confirmed | conditional |  | `mechanisms/M3_EXACT_CHECKS.json`<br>`mechanisms/verify/V1_INDEPENDENT_DERIVATIONS.json` | Algebra exact, but the exclusion rests on a cited stability theorem and the Z2 topology: conditional, not exact. The linear Lambda_4 formula holds only near lambda=4/3. |
| MX-4 | The collapse fate is not ekpyrotic (w_eff=0.76-0.85). | negative | partially-confirmed | inconclusive |  | `mechanisms/M5_OBSERVATIONAL_SCREEN.json` | Covers only ~1.05-1.10 e-folds of resolved contraction; negative only on that range. |
| MX-5 | 4D energy budget with zero-mode Planck masses: radiation can dominate the residual RS vacuum for at most 0.59 e-folds (0.47 at delta=0.1); ~21.8 are needed. | conditional | confirmed | conditional | ✔ | `mechanisms/M2_ENERGY_BUDGET_AND_SCALES.json`<br>`mechanisms/verify/V4_ARITHMETIC_CHECKS.json` | H_vac^2/H0^2=0.368 only as delta->0 (0.372 at 0.01, 0.406 at 0.1). Conditional on the 4D EFT and NEC. |
| MX-6 | If H_vac is today's dark-energy rate, the blast e-folds in 6.4 Gyr and delta~1.1e-62..7.6e-68. | conditional | confirmed | conditional | ✔ | `mechanisms/M2_ENERGY_BUDGET_AND_SCALES.json` |  |
| MX-7 | Sudden tension-to-radiation conversion is gravitationally screened (the Weyl term must cancel 80-99.9% of the radiation terms). | exact-verified | confirmed | exact-verified | ✔ | `mechanisms/M3_EXACT_CHECKS.json` |  |
| MX-8 | A dissipative shell coupling at registered c captures 0.37-2.2% of the tension drop; radiation never reaches the residual-vacuum threshold (R/R_crit<=0.038). | negative | confirmed | negative | ✔ | `mechanisms/M6_PILOT_COUPLED_5D.json`<br>`mechanisms/verify/V2_PILOT_CONSISTENCY.json` | 5D pilots run at delta=0.1 (100x the registered detuning). [critic 28 Sept 2026] |
| MX-9 | Gravitational particle production is short by 85-122 orders of magnitude. | negative | confirmed | negative | ✔ | `mechanisms/M5_OBSERVATIONAL_SCREEN.json` | Order of magnitude, assumed coefficient. |
| MX-10 | The registered linear tension cannot tune away the vacuum: the initial-shell family ends near c->0+ and the shell cannot be built at c*. | numerical | refuted | corrected |  | `mechanisms/verify/V3B_CONTINUE_TO_CSTAR.json`<br>`mechanisms/verify/V5_CSTAR_PILOT_SIGNED.json` | The stop was a parametrisation artefact (phi_h=-1+10^x). With the sign allowed, a static shell exists all the way to c*=-0.99307 at delta=0.1 (phi_h=-1.00695, phi_b=-1.7287), and the c* pilot (single grid, one seed sign, Y=0, delta=0.1 only) rolls from phi_b=-1.73 across phi=-1 to 0.80. It is a qualitatively different initial state; its endpoint is open. |
| MX-11 | With a quadratic tension term tuned to remove the vacuum (d*~-3.1, a model change), shell radiation reaches 69% of H^2, but Weyl radiation is ~1.8x the radiation when the run stops. | inconclusive | partially-confirmed | inconclusive |  | `mechanisms/M8_QUADRATIC_TENSION_TUNING.json`<br>`mechanisms/verify/V2_PILOT_CONSISTENCY.json`<br>`mechanisms/verify/pilot/ext_dstar_Y1_dzf2e-3_phimax1.3_summary.json` | The runs stopped at the solver's phi_b>0.995 cutoff, not the chart freeze. The H^2 budget has large cancellations (vacuum -117%, Weyl +121%, radiation +69%). With the cutoff lifted, Weyl/radiation is 1.23 at the end of the reliable window (about 12x the Delta N_eff-derived limit of 0.1). Values after the chart freeze are not usable: proper time advances by only ~0.2 H0 tau while the ratio falls through 0.6 and below zero and then swings back as H turns negative (gauge artefact or reversal); the earlier '~0.6 after the freeze' figure is withdrawn (CRITIC_CHECKS.json, post_freeze_dstar). Pilots at delta=0.1 only; d* convergence rests on one grid pair (dz_fine 1e-3 and 2e-3). [critic 28 Sept 2026] |

## Literature sweeps (theory, observations, methods) (`literature/`)

| ID | Claim | Producer label | Verifier verdict | Status here | Est. | Evidence | Correction / note |
|---|---|---|---|---|---|---|---|
| LT-1 | 'Our universe as a wall born in a 5D event' has substantial published prior art (1998-2026); HDBLAST cannot claim the concept. | context | not-independently-audited | context |  | `literature/theory.md`<br>`literature/theory_checks/theory_sources.json` | Snippet-level only. |
| LT-2 | Registered model sits inside standard RS/KR/DFGK formulas: balanced tension = critical tension; thin-brane H^2 differs from the branch by exactly -c^2 delta^2/384 (7.869 ppm in H). | exact-verified | not-independently-audited | exact-verified |  | `literature/theory_checks/rs_thin_brane_crosscheck.json` | Self-checked with 4 wrong-formula controls; no separate auditor. The O(delta^2) term predicts -7.849 ppm of the measured -7.869 ppm shift in H; higher orders supply the rest. [critic 28 Sept 2026] |
| LT-3 | The registered PTA knee is at 0.9972 x (1/yr), on the pulsar position/proper-motion blind spot; no 2025-2026 release tests it at likelihood level. | numerical | not-independently-audited | inconclusive |  | `literature/observations_checks/knee_and_scale_checks.json`<br>`literature/observations.md` | Knee position numerical; the blind-spot statement is from snippets. |
| LT-4 | Tabletop gravity bounds put the static branch at >=TeV vacuum scale and sub-mHz horizon-scale GW frequency, ~1600-2600x the knee frequency. | conditional | not-independently-audited | conditional |  | `literature/observations_checks/knee_and_scale_checks.json` |  |
| LT-5 | Methods sweep ran no new searches (budget exhausted); its list is second-hand. | procedural-negative | not-independently-audited | procedural-negative |  | `literature/methods.md` |  |
| LT-6 | Certifying the +1 branch needs the eta=phi-1 formulation: the regular cone germ is amplified by 6.288e18 at delta=1e-3 (a direct phi formulation needs ~39 digits). | numerical | not-independently-audited | numerical |  | `literature/methods_checks/certification_pilot.json` |  |

## Program map (`program_map/`)

| ID | Claim | Producer label | Verifier verdict | Status here | Est. | Evidence | Correction / note |
|---|---|---|---|---|---|---|---|
| PM-1 | Program map of HDBLAST Aug 2025 - 23 Sept 2026 with claims ledger; 34/34 self-checks pass. | context | not-independently-audited | context |  | `program_map/PROGRAM_MAP.md`<br>`program_map/program_map_checks.json` | Zenodo record 22922928 metadata (version, concept, files) unverified: zenodo.org blocked. |

## Synthesis checks done for this report (`synthesis/`)

| ID | Claim | Producer label | Verifier verdict | Status here | Est. | Evidence | Correction / note |
|---|---|---|---|---|---|---|---|
| SY-1 | Methods A and B agree where they overlap: calibration growth to 1.6e-8 (eigenvalue) and 1.6e-7 (time domain); the artificial flipped-sigma'' tachyon to 3.7e-10; both find no scalar mode below the 9/4 continuum; B's leading rate matches the continuum edge -3/2 (A's criterion). | numerical | not-independently-audited | numerical |  | `synthesis/RECONCILIATION.json`<br>`synthesis/reconcile_stability.py` | Comparison of saved outputs; wrong mu^2 maps are rejected as controls. Slowest-decay agreement is within the supported +/-1e-3 of TD-4. [critic 28 Sept 2026] |
| SY-2 | Method A's shell coefficient B extrapolates to 3.5555554 as delta->0, matching the exact 32/9 = 14/9 (regular Delta=18 growth rate) + 2 (sigma''/2) to 3e-8. | numerical | not-independently-audited | numerical |  | `synthesis/RECONCILIATION.json` |  |
| SY-3 | Cross-workstream reconciliation: the theory sweep's 'no dark radiation' statement applies to static pure-tension configurations; the dynamical roll-off does produce a positive Weyl term (MX-1). | context | not-independently-audited | context |  | `literature/theory.md`<br>`mechanisms/M1_WEYL_ALONG_CHAT14.json` |  |
| SY-4 | All 56 key numbers quoted in this checkpoint's report are found in the workstreams' JSON outputs within stated tolerances; 3 negative controls pass. | numerical | not-independently-audited | numerical |  | `synthesis/QUOTED_NUMBER_CHECKS.json`<br>`synthesis/check_quoted_numbers.py` |  |

## Reconciliation of the two stability methods

| Question | Method A (gauge-invariant) | Method B (time domain) | Agree? | Why / scope |
|---|---|---|---|---|
| Calibration growth of the original shell | 1.657193631259 | 1.657193368 (pencil), 1.6571936573 (eigenvalue) | yes, 1.6e-7 / 1.6e-8 | same physical tachyon; B's residual is time-step/grid error |
| Artificial tachyon (sign of σ″ flipped), δ=0.001 | μ²=−28555.115 → rate 167.4892463 | rate 167.4892463 (shift-invert) | yes, 3.7e-10 | independent codes see the same wrong-model instability, so both can detect one |
| Scalar bound states below the 9/4 continuum on the +1 branch | none at 6 δ | none (homogeneous sector) at δ=0.001, 0.003, 0.01 | yes | B covers only the homogeneous (FRW-symmetric) sector; every dS₄ harmonic μ² has a homogeneous representative, so both test the same eigenvalue condition |
| Slowest decay | continuum edge μ²=9/4 ⇒ e^{−3H₊τ/2} | −1.500 ± 0.001 in H₊ units | yes | = −0.910 H₀ |
| Tensor sector | only the massless graviton, gap 3H/2 | not tested | — | analytic_structure's leading-order spectrum {0} ∪ [9/4, ∞) agrees with A |
| Positive discrete rates | none physical | gauge mode (λ→1), far-boundary artefacts, grid-scale modes | no conflict | B's positive rates have no shell support and change with L or h; they are properties of the truncated grid |
| The λ≈1 mode | ℓ=1 harmonic μ²=−4 is pure gauge (no physical mode since Bφ′_b≠0) | residual gauge mode, λ=0.99999993 at L=8 | consistent | λ=1 ⇔ μ²=−λ(λ+3)=−4; mapped value agrees with −4 to <1e-6 |
| δ→0 limit | B-coefficient → 3.5555554 (fit) | — | matches exact 32/9 to 3e-8 | new cross-check (`synthesis/RECONCILIATION.json`) |

Overall: the two methods agree wherever they overlap, and their scopes are complementary. Neither is interval-certified; the vector sector, nonlinear stability, quantum tunnelling and dynamical attraction remain open.

## Citations

Paper sites were blocked for downloads and the shared web-search budget ran out during the audits, so no auditor could re-verify any citation. All citations in this round are search-snippet-only or taken from earlier programme records (`LITERATURE_2022_2026.md`).

A cross-check by `critic_corrections.py` (`CRITIC_CHECKS.json`, key `citations`) compares every arXiv identifier and link in the top-level files, workstream READMEs and audit reports with the literature files (`literature/`, `LITERATURE_2022_2026.md`, `synthesis/LITERATURE_MERGED.json`). The following have **no record** there. Their provenance is only the workstream's own statement, so treat them as unverified:

- `analytic_structure`: arXiv gr-qc/0010035, hep-th/0011156, hep-th/0012044.
- `evolution_constraints`: arXiv 1503.03436, gr-qc/0407110, gr-qc/0504114.
- `mechanisms`: arXiv 2312.09042, 2603.13226, astro-ph/0208133, gr-qc/0111013, hep-th/0111279, hep-th/0206146, hep-th/9906064; 3 non-arXiv link(s): <https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevD.64.123522/fulltext>, <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.102.063514>, <https://link.springer.com/article/10.1007/JHEP02(2024)102>.
- `preheating`: arXiv 0905.2284, 1402.5669, hep-ph/0205240, hep-th/0405016; 1 non-arXiv link(s): <https://journals.aps.org/prd/abstract/10.1103/PhysRevD.9.341>.
- `stability_gauge_invariant`: arXiv astro-ph/9702174; 2 non-arXiv link(s): <https://www.osti.gov/etdeweb/biblio/5746769>, <https://www.sciencedirect.com/science/article/abs/pii/0550321384903948>.
- `top level`: arXiv 2603.13226.

Corrections to the earlier version of this section: arXiv:2609.21421 *does* have a search-snippet record (theory E11, methods C9). Not every unrecorded citation is context only: the ΔN_eff input 2.990 ± 0.070 (arXiv:2603.13226) feeds the dark-radiation limit used in MX-2 and MX-11, and the ACT DR6 N_eff value is quoted differently in two sweeps (2.86 ± 0.13 and 2.89 ± 0.11). The negative conclusions of MX-2 and MX-11 hold for either value, but the thresholds must be checked against the papers before any external use. No computed (non-observational) result depends on these citations.

## Process notes

- The analytic_structure auditor accidentally regenerated that folder's `MANIFEST.sha256.json`, adding 25 `verify/` entries. `synthesis/QUOTED_NUMBER_CHECKS.json` confirms that all 29 original entries still match the current files. The file was left as found; the root `MANIFEST.sha256.json` covers every file.
- No claim in this checkpoint is a discovery, a proof, a hot Big Bang mechanism, or observational support for the hypothesis. Novelty is claimed only relative to this project.

## Critic review (28 September 2026)

A final completeness and overclaim review of the top-level documents (`00_READ_FIRST.md`, `README.md`,
`PUBLIC_SUMMARY.md`, `NEXT_TESTS.md`, `REPRODUCE.md`, this file and `RESULTS_SUMMARY.json`) against every workstream
README, audit report and JSON output. Workstream outputs were not modified. Checks: `critic_corrections.py` →
`CRITIC_CHECKS.json`.

Corrections made:

1. **d\* pilot after the chart freeze (MX-11).** The report quoted Weyl/radiation "≈ 0.6 after the chart freeze,
   ≥ 6× the N_eff limit". No JSON contained that value, and the mechanisms workstream itself rules post-freeze values
   out as gauge artefacts. Recomputed from the saved time series: after the freeze proper time advances by only
   ≈ 0.2 H₀τ while the ratio falls through 0.6 and below zero and then swings back as H turns negative, so there is
   no plateau. The figure is withdrawn; the supported statement is 1.23 at the end of the reliable window (≈ 12× the
   ΔN_eff-derived limit of 0.1).
2. **Stability scope.** "Two completely different methods" and "stable at δ = 0.0003–0.1" were stronger than the
   evidence: Method B covers only the homogeneous sector at three δ; Method A scans real μ² on [−400, 2.2499] and
   relies on the conditional B > 0 theorem beyond it and on cone regularity as the normalisability criterion. The
   planted-tachyon controls test only a strong instability; the auditor's B\* margin is the relevant detection-power
   evidence. Method B's time-domain runs barely excite the shell (≈ 8.5×10⁻⁹), so its conclusion rests on the spectra.
3. **Decay-rate precision.** The reconciliation row "−1.5001, agree within 1×10⁻⁴" contradicted the audit
   (−1.500 ± 0.001, order-sensitive Richardson); corrected.
4. **Pilots at δ = 0.1.** The dissipative-coupling, d\* and c\* pilots all ran at δ = 0.1 (100× the registered
   detuning); the c\* run is a single-grid, single-seed, Y = 0 exploration; d\* convergence rests on one grid pair.
   These caveats were missing from the summary documents.
5. **Preheating bound.** The fixed-trajectory condition and the "post" cutoff maximum (≤ 4.0×10⁻³) were added next
   to the headline ≤ 1.7×10⁻⁴.
6. **Smaller numerical wording.** δ–Λ lock range is 1.1×10⁻⁶² to 7.6×10⁻⁶⁸ (not "10⁻⁶²–10⁻⁶⁷"); the O(δ²)
   thin-brane term explains −7.849 of the −7.869 ppm shift (not all of it); H_vac²/H₀² is 0.372 at δ = 0.01 and
   0.406 at δ = 0.1; the multiprecision solvers' agreement is 10⁻³⁹–10⁻⁴¹, their working precision 40–65 digits.
7. **Contradictions between documents.** README and 00_READ_FIRST said every workstream had an independent audit,
   while the literature sweeps and program map are self-checked only; REPRODUCE said every audited script was re-run
   bit-identically, while several were re-run only in subsets, some agree to 10⁻⁹–10⁻⁸, and m7 was not re-run.
8. **Reproduction.** Commands for the audit scripts that produce quoted corrections (mechanisms v1/v2/v4 and the
   extended d\* pilot, preheating `check_claims.py`, evolution `recompute_saved.py`, time-domain
   `v4_compare_and_extract.py`) and the order of the post-build critic step were added to `REPRODUCE.md`.
9. **Citations.** The unverified-citation list was completed (section above); CMB-S4, which has no literature record,
   was removed from `NEXT_TESTS.md`; the N_eff normalisation mismatch (ρ_DR/ρ_SM at production vs ρ_dr/ρ_γ at
   recombination) is now flagged in the pre-registered decision rule of next test A1.
10. **AI product names** were removed from the top-level files.

Remaining gaps the critic could not fix (outside the top level): workstream READMEs keep their pre-audit wording
(for example the stale preheating decay column and the "13–100×" damping range); `synthesis/build_results_summary.py`
and `synthesis/build_claim_review.py` still emit the pre-critic text, which this script re-corrects;
`analytic_structure/MANIFEST.sha256.json` still carries the 25 accidental `verify/` entries; `__pycache__` folders
are present in some workstream folders.
