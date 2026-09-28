#!/usr/bin/env python3
"""Build RESULTS_SUMMARY.json for the 27 Sept 2026 HDBLAST checkpoint.

Every claim carries: the producer's label, the independent verifier's verdict,
the status used in this checkpoint after applying the verifier's corrections,
and evidence paths (relative to the checkpoint root) that are checked to exist.

Rules enforced (the script refuses to write if violated):
  - status in the allowed vocabulary;
  - "established" (the boolean) only when the verifier verdict is "confirmed";
  - a "refuted" verdict must carry a correction and final status "corrected";
  - every evidence path exists.
A negative control feeds one deliberately bad record through the validator and
requires it to be rejected.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

STATUSES = {"exact-verified", "numerical", "conditional", "negative", "inconclusive", "corrected",
            "procedural-negative", "context"}
VERDICTS = {"confirmed", "partially-confirmed", "refuted", "unverifiable", "not-independently-audited"}

C = []


def claim(cid, ws, text, producer, verdict, status, evidence, correction="", note=""):
    C.append({"id": cid, "workstream": ws, "claim": text, "producer_label": producer,
              "verifier_verdict": verdict, "status": status, "established": verdict == "confirmed"
              and status not in ("corrected", "inconclusive"),
              "evidence": evidence, "correction": correction, "note": note})


GI = "stability_gauge_invariant"
TD = "stability_time_domain"
AS = "analytic_structure"
PH = "preheating"
EC = "evolution_constraints"
MX = "mechanisms"
LT = "literature"
PM = "program_map"
SY = "synthesis"

# ---------------------------------------------------------------- stability A
claim("GI-1", GI, "The +1 static background is reproduced to all stated digits (phi_b=0.9999159473169134, "
      "rho_b=129.92476284968, H^2=5.924014794329e-5, H/H0=0.606721732).", "numerical", "partially-confirmed",
      "numerical", [f"{GI}/BACKGROUNDS.json", f"{GI}/VERIFICATION.md"],
      "Independent-polish agreement in rho_b is <=8.5e-14 relative, not 6e-15 as first quoted.")
claim("GI-2", GI, "Gauge-invariant linearized system (X, Z=Psi/phi') and shell condition M=BX+3(mu^2+4)Z/rho^2=0 "
      "with B=phi''/phi'+sigma''/2; no brane bending in longitudinal gauge; tensor h''+4Hh'+mu^2h/rho^2=0, h'(y_b)=0.",
      "exact-verified", "confirmed", "exact-verified",
      [f"{GI}/DERIVATION_RESULTS.json", f"{GI}/verify/indep_derivation.json"],
      note="55/55 symbolic checks, 4/4 wrong-formula controls; auditor's independent derivation 26/26.")
claim("GI-3", GI, "Calibration on the original unstable shell: one bound state mu^2=-7.717871625260, growth "
      "1.657193631259 (recorded Chat 9: -7.717871625176).", "numerical", "confirmed", "numerical",
      [f"{GI}/SPECTRUM_RESULTS.json", f"{GI}/WINDING_RESULTS.json"],
      note="Validation against a known internal target, not a blind test.")
claim("GI-4", GI, "On the +1 branch at delta in {0.0003,0.001,0.003,0.01,0.03,0.1} there is no scalar eigenvalue "
      "with mu^2<9/4 (min normalized mismatch >=0.99274).", "numerical", "confirmed", "numerical",
      [f"{GI}/SPECTRUM_RESULTS.json", f"{GI}/verify"])
claim("GI-5", GI, "No complex scalar eigenvalues inside [-60,2]x[-30,30] and [-400,2]x[-200,200] (winding 0; 1 on "
      "calibration).", "numerical", "partially-confirmed", "numerical", [f"{GI}/WINDING_RESULTS.json"],
      note="Auditor re-ran only a subset (small contour, 3 cases).")
claim("GI-6", GI, "If B>0 at the shell, all scalar eigenvalues are real and mu^2>-4; B=3.479-3.555 on the +1 "
      "branch and -2.68e-4 on the original shell.", "conditional", "confirmed", "conditional",
      [f"{GI}/DERIVATION_RESULTS.json", f"{GI}/verify/identities.json"])
claim("GI-7", GI, "No physical l=1 mode (B phi'_b != 0) and no l=0 static zero mode (junction Jacobian "
      "condition number 15-17).", "numerical", "confirmed", "numerical", [f"{GI}/SPECTRUM_RESULTS.json"])
claim("GI-8", GI, "Tensor sector: only the normalizable massless graviton; no KK bound state in (0,9/4); no "
      "tachyons.", "exact-verified", "partially-confirmed", "numerical",
      [f"{GI}/SPECTRUM_RESULTS.json", f"{GI}/VERIFICATION.md"],
      "Label split: mu^2>=0 (Sturm-Liouville) is exact; absence of bound states in (0,9/4) is numerical.")
claim("GI-9", GI, "Controls C1-C5 behave as expected (Gegenbauer limit, flipped sigma'' creates a tachyon, dropping "
      "gravity removes the calibration root, 4D limit, perturbed c).", "numerical", "partially-confirmed",
      "numerical", [f"{GI}/CONTROLS_RESULTS.json"],
      "C1 pass criterion was chosen after seeing an O(delta) deviation; C2's root (-28555) lies outside the main "
      "scan, so it is a weak demonstration of detection power.")
claim("GI-10", GI, "The static +1 branch is linearly stable in the scalar (including brane bending and both "
      "junctions) and tensor sectors at all six tested delta; the only marginal mode is the massless graviton.",
      "numerical", "confirmed", "numerical", [f"{GI}/README.md", f"{GI}/VERIFICATION.md"],
      note="Floating point, not interval-certified; vector sector, nonlinear stability, tunnelling and "
           "dynamical attraction not studied.")
# ---------------------------------------------------------------- stability B
claim("TD-1", TD, "Time-domain calibration on the original shell: rate 1.657193368 (linear) / 1.657193408 "
      "(nonlinear) vs 1.6571936312.", "numerical", "confirmed", "numerical",
      [f"{TD}/results/KEY_RESULTS.json", f"{TD}/VERIFICATION.md"])
claim("TD-2", TD, "Only gauge (lambda->1) and far-boundary (~0.6/L) modes have positive rates on the +1 branch.",
      "numerical", "partially-confirmed", "numerical", [f"{TD}/results/KEY_RESULTS.json", f"{TD}/VERIFICATION.md"],
      "True for the fully discrete RK4 system and for shell observables. The semi-discrete operator also has "
      "grid-scale shell-localised modes with Re(lambda)~+0.245 that do not shrink with refinement (numerical "
      "artefacts damped by RK4).")
claim("TD-3", TD, "Every shell-supported mode lies on Re(lambda)=-3/2 (dS continuum mu^2>9/4) to ~1e-8; no scalar "
      "bound state with mu^2<9/4 in the homogeneous sector.", "numerical", "confirmed", "numerical",
      [f"{TD}/results/KEY_RESULTS.json", f"{TD}/verify"])
claim("TD-4", TD, "Leading physical rate on the +1 branch converges to -3/2 (in H_+ units).", "numerical",
      "partially-confirmed", "numerical", [f"{TD}/results/KEY_RESULTS.json"],
      "Supported precision is -1.500 +/- 1e-3 (Richardson depends on assumed order). The 'shell-localised' second "
      "profile is actually a bulk bump (near-shell weight 2.4e-10); the time-domain test barely excites the shell "
      "(~1e-8), so the spectral evidence carries the no-bound-state claim.")
claim("TD-5", TD, "Flipping sigma'' in the scalar junction produces lambda=167.48925 (partner sums to -3).",
      "numerical", "confirmed", "numerical", [f"{TD}/results/KEY_RESULTS.json"])
claim("TD-6", TD, "Whole-domain constraints do not converge in the far stretched region (under-resolved outgoing "
      "radiation); near-shell constraints converge.", "negative", "confirmed", "negative",
      [f"{TD}/results/KEY_RESULTS.json"])
claim("TD-7", TD, "Model equations, junctions and normalisations re-derived (32/32 sympy checks).",
      "exact-verified", "confirmed", "exact-verified", [f"{TD}/verify"])
# ---------------------------------------------------------------- analytic
claim("AS-1", AS, "The 22 Sept O(delta) and O(delta^2) coefficients are reproduced exactly.", "exact-verified",
      "confirmed", "exact-verified", [f"{AS}/SERIES_COEFFICIENTS.json", f"{AS}/verify/V1_SERIES_INDEPENDENT.json"])
claim("AS-2", AS, "phi_b, H^2, rho_b, y_b and eta_h expansions exact in c through O(delta^8); e.g. "
      "h3=-c^2(11767c+11362)/2826240.", "exact-verified", "confirmed", "exact-verified",
      [f"{AS}/SERIES_COEFFICIENTS.json", f"{AS}/verify/V1_SERIES_INDEPENDENT.json"])
claim("AS-3", AS, "30-65 digit solutions of the nonlinear boundary-value problem confirm the new coefficients "
      "(e.g. H^2(delta=1e-3)=5.9240147943288795996237552780e-5).", "numerical", "confirmed", "numerical",
      [f"{AS}/SERIES_VS_BVP.json", f"{AS}/runs/BVP_SCAN_reg.json", f"{AS}/verify/V4_COMPARE.json"])
claim("AS-4", AS, "The delta^9 shell coefficients depend on cone regularity and need global matching.",
      "numerical", "refuted", "corrected", [f"{AS}/VERIFICATION.md", f"{AS}/verify/V1_SERIES_INDEPENDENT.json"],
      "The delta^9 coefficients are local and exactly computable (resonance constant drops out): at registered c, "
      "eta_b: 3.249147254269e-6, H^2: -2.759694789e-7; they match the workstream's fits and an independent BVP.")
claim("AS-5", AS, "Exact resonance at (m,j)=(2,7) with S=-24576/17; no log in shell observables at delta^9; the "
      "cone value eta_h gets a delta^16 ln(delta) term.", "numerical", "confirmed", "exact-verified",
      [f"{AS}/RESONANCE_PROBE.json", f"{AS}/verify/V1_SERIES_INDEPENDENT.json"],
      note="Verifier upgraded 'no log at delta^9' from fitted to exact.")
claim("AS-6", AS, "Leading cone displacement eta_h=-(111537/131072)c(1+c)^7 delta^8[1+O(delta)].",
      "exact-verified", "confirmed", "exact-verified", [f"{AS}/SERIES_COEFFICIENTS.json", f"{AS}/verify/V4_COMPARE.json"])
claim("AS-7", AS, "Degree-14 Gegenbauer structure: n(n+4)=252, n=Delta-4 with Delta=3W''/W=18 at phi=+1; "
      "closed forms; no polynomial at phi=-1.", "exact-verified", "confirmed", "exact-verified",
      [f"{AS}/GEGENBAUER_STRUCTURE.json", f"{AS}/verify/V3_GEGENBAUER_SPECTRUM.json"],
      "Wording: at phi=-1, s=3W''/W=-18/5 equals 4-Delta, not Delta.")
claim("AS-8", AS, "Leading-order tensor quantization (mu-3/2)P^{-mu}_{1/2}(cosh u_b)=0: massless graviton plus "
      "continuum from m=3H/2 (the gap itself is known literature).", "exact-verified", "confirmed",
      "exact-verified", [f"{AS}/SPECTRUM_LEADING_ORDER.json", f"{AS}/verify/V3_GEGENBAUER_SPECTRUM.json"],
      note="Leading order in delta only.")
claim("AS-9", AS, "Decoupled bulk scalar: no zero mode or tachyon for any H; no mode below 9H^2/4 for H/k<16.65.",
      "conditional", "confirmed", "conditional", [f"{AS}/SPECTRUM_LEADING_ORDER.json"],
      note="Neglects scalar-metric mixing; superseded for the actual branch by GI-4/GI-10, which include mixing.")
claim("AS-10", AS, "Stability of the +1 branch is not established by this workstream.", "negative", "confirmed",
      "context", [f"{AS}/README.md"], note="Resolved within scope by the two stability workstreams (GI-10, TD-3).")
# ---------------------------------------------------------------- preheating
claim("PH-1", PH, "Flat linear-crossing calibration reproduces the instant-preheating spectrum to <=8.4e-7.",
      "numerical", "confirmed", "numerical", [f"{PH}/CONTROLS.json", f"{PH}/verify/INDEP_MODES.json"],
      note="The two wrong-formula controls are fixed by Gaussian algebra; weak as solver tests.")
claim("PH-2", PH, "On the archived trajectory N=q^{3/2}/(8pi^3)(1+D/q+...), with D=10.648, 7.0034, 3.0032, 1.3416 "
      "at phi*=0.25,0.5,0.75,0.9.", "numerical", "partially-confirmed", "numerical", [f"{PH}/CONTROLS.json"],
      "G x (relative deviation) is constant to 0.4%/0.1% (not '4 digits'); D at phi*=0.25 and 0.9 are closed-form "
      "values (measured 10.6476 and 1.3415).")
claim("PH-3", PH, "Closed-form O(1/q) coefficient D for a curved, expanding crossing matches toy backgrounds to "
      "~3e-5.", "numerical", "confirmed", "numerical", [f"{PH}/CORRECTION_FIT.json", f"{PH}/verify/DDP_EXPONENT.json"],
      "Derivation note corrected: exponent-only h^2 term is 15/(2pi), so c_hh and c_hA follow from the mass shift "
      "plus exponent; only the +/-pi/8 pieces remain fitted. No novelty claimed.")
claim("PH-4", PH, "At the M5 cutoff the produced energy is <=1.7e-4 of the residual-vacuum radiation threshold for "
      "screen-passing couplings (<=2.8e-3 over the whole scan).", "conditional", "confirmed", "conditional",
      [f"{PH}/SUMMARY.json", f"{PH}/verify/CHECK_CLAIMS.json"],
      note="Evaluated at the data end; the maximum over the post-window profile is 5.9e-4. Conditional on the "
           "fixed archived trajectory, the M5 cutoff and instantaneous full conversion.")
claim("PH-5", PH, "Reaching the threshold needs chi masses >=18 M5 (>=6.3 M5 if only post-crossing masses are "
      "constrained); stronger coupling makes it worse (R ~ K/sqrt(G)).", "conditional", "confirmed", "conditional",
      [f"{PH}/SUMMARY.json"])
claim("PH-6", PH, "Wherever the channel could matter (R~1), backreaction on the scalar junction is 10-79%, so the "
      "fixed-trajectory calculation stops being self-consistent.", "conditional", "confirmed", "conditional",
      [f"{PH}/SUMMARY.json"])
claim("PH-7", PH, "Unsubtracted one-loop tension shift is 1.1-3.7 x |sigma'| at the favoured points: keeping the "
      "registered sigma requires tuned counterterms.", "conditional", "confirmed", "conditional",
      [f"{PH}/PREHEATING_RESULTS.json"])
claim("PH-8", PH, "Perturbative decay chi->psi psibar gives a radiation/vacuum ratio <=1.3e-4.", "conditional",
      "refuted", "corrected", [f"{PH}/SUMMARY.json", f"{PH}/verify/CHECK_CLAIMS.json"],
      "Stale README number: the JSON maximum is 2.34e-4 at (phi*=0.5, G=300, y=1). Conclusion (<<1) unchanged.")
claim("PH-9", PH, "Instant preheating on the archived +1 roll-off does not produce a radiation era in any scanned "
      "coupling range.", "conditional", "confirmed", "negative", [f"{PH}/README.md", f"{PH}/VERIFICATION.md"],
      note="Conditional negative; the bulk-coupled feedback (all three junctions) is not solved.")
# ---------------------------------------------------------------- evolution constraints
claim("EC-1", EC, "The 22 Sept constraint stall is reproduced (finest-pair order 0.24 at t=0.5, 0.27 at t=1).",
      "numerical", "confirmed", "numerical", [f"{EC}/ANALYSIS.json", f"{EC}/verify/compare_reruns.json"])
claim("EC-2", EC, "The error front is seeded within ~0.02 of the shell for t<~0.05 and carried outward by the exact "
      "transport law of e^{3A}(H+2M).", "numerical", "partially-confirmed", "numerical", [f"{EC}/ANALYSIS.json"],
      "The geometric gain applies to the un-normalised H+2M; the background-normalised H has an extra factor ~360. "
      "The t=0.05 density diagnostic depends on the time step.")
claim("EC-3", EC, "The stall comes from an initial-data defect floor near the shell that does not decrease with h.",
      "numerical", "confirmed", "numerical", [f"{EC}/ANALYSIS.json", f"{EC}/verify/recompute_saved.json"])
claim("EC-4", EC, "Discrete-constraint projection of the initial data restores convergence (orders ~2.9-3.1 at the "
      "finest pair; H_max(t=0.5) 6.38e-3 -> 9.35e-4).", "numerical", "confirmed", "numerical",
      [f"{EC}/ANALYSIS.json"], note="The projection correction itself does not converge below h=1e-4 (~1e-9).")
claim("EC-5", EC, "Damped transport law (d_t-d_z)[e^{3A}C+]=-kappa e^{3A}C+ for the outgoing constraint source.",
      "exact-verified", "confirmed", "exact-verified", [f"{EC}/CONTROLS.json", f"{EC}/verify/sympy_independent.json"])
claim("EC-6", EC, "Outgoing-characteristic damping (kappa=10, off near the shell) lowers H by 13-100x.",
      "numerical", "partially-confirmed", "numerical", [f"{EC}/ANALYSIS.json"],
      "Measured factors are ~25-200 at t=0.5-1 and ~8 at t=0.25; at t=1 'A+B smallest' is a tie with B.")
claim("EC-7", EC, "Asymptotic fourth order for the eps=0.01 seed is not demonstrated (phi_t and B_t converge at "
      "~2).", "inconclusive", "confirmed", "inconclusive", [f"{EC}/ANALYSIS.json", f"{EC}/verify/recompute_saved.json"],
      note="Verifier localised the deficit on the shell-seeded front z~-t; no constraint remedy touches it.")
claim("EC-8", EC, "Uniform damping, H-only damping, KO dissipation, tighter inversion and order 6 without projection "
      "do not help or are unstable.", "negative", "confirmed", "negative", [f"{EC}/ANALYSIS.json"])
# ---------------------------------------------------------------- mechanisms
claim("MX-1", MX, "The blast produces a positive bulk Weyl (dark-radiation-like) term, peaking at 0.0806-0.0809 "
      "H0^2 (~15% of H^2) at H0 tau~6.24; late value not converged.", "numerical", "confirmed", "numerical",
      [f"{MX}/M1_WEYL_ALONG_CHAT14.json"], note="3 independent grids (not 4).")
claim("MX-2", MX, "Bulk Weyl radiation cannot be the hot Big Bang and is capped by N_eff at a few percent of the "
      "radiation.", "negative", "partially-confirmed", "negative", [f"{MX}/M5_OBSERVATIONAL_SCREEN.json"],
      "The N_eff input values come from search snippets and are unverifiable here.")
claim("MX-3", MX, "Dark-bubble cosmology is not available in the registered model.", "negative",
      "partially-confirmed", "conditional", [f"{MX}/M3_EXACT_CHECKS.json", f"{MX}/verify/V1_INDEPENDENT_DERIVATIONS.json"],
      "Algebra exact, but the exclusion rests on a cited stability theorem and the Z2 topology: conditional, not "
      "exact. The linear Lambda_4 formula holds only near lambda=4/3.")
claim("MX-4", MX, "The collapse fate is not ekpyrotic (w_eff=0.76-0.85).", "negative", "partially-confirmed",
      "inconclusive", [f"{MX}/M5_OBSERVATIONAL_SCREEN.json"],
      "Covers only ~1.05-1.10 e-folds of resolved contraction; negative only on that range.")
claim("MX-5", MX, "4D energy budget with zero-mode Planck masses: radiation can dominate the residual RS vacuum for "
      "at most 0.59 e-folds (0.47 at delta=0.1); ~21.8 are needed.", "conditional", "confirmed", "conditional",
      [f"{MX}/M2_ENERGY_BUDGET_AND_SCALES.json", f"{MX}/verify/V4_ARITHMETIC_CHECKS.json"],
      note="H_vac^2/H0^2=0.368 only as delta->0 (0.372 at 0.01, 0.406 at 0.1). Conditional on the 4D EFT and NEC.")
claim("MX-6", MX, "If H_vac is today's dark-energy rate, the blast e-folds in 6.4 Gyr and delta~1e-62..1e-67.",
      "conditional", "confirmed", "conditional", [f"{MX}/M2_ENERGY_BUDGET_AND_SCALES.json"])
claim("MX-7", MX, "Sudden tension-to-radiation conversion is gravitationally screened (the Weyl term must cancel "
      "80-99.9% of the radiation terms).", "exact-verified", "confirmed", "exact-verified", [f"{MX}/M3_EXACT_CHECKS.json"])
claim("MX-8", MX, "A dissipative shell coupling at registered c captures 0.37-2.2% of the tension drop; radiation "
      "never reaches the residual-vacuum threshold (R/R_crit<=0.038).", "negative", "confirmed", "negative",
      [f"{MX}/M6_PILOT_COUPLED_5D.json", f"{MX}/verify/V2_PILOT_CONSISTENCY.json"])
claim("MX-9", MX, "Gravitational particle production is short by 85-122 orders of magnitude.", "negative",
      "confirmed", "negative", [f"{MX}/M5_OBSERVATIONAL_SCREEN.json"], note="Order of magnitude, assumed coefficient.")
claim("MX-10", MX, "The registered linear tension cannot tune away the vacuum: the initial-shell family ends near "
      "c->0+ and the shell cannot be built at c*.", "numerical", "refuted", "corrected",
      [f"{MX}/verify/V3B_CONTINUE_TO_CSTAR.json", f"{MX}/verify/V5_CSTAR_PILOT_SIGNED.json"],
      "The stop was a parametrisation artefact (phi_h=-1+10^x). With the sign allowed, a static shell exists all "
      "the way to c*=-0.99307 at delta=0.1 (phi_h=-1.00695, phi_b=-1.7287), and the c* pilot rolls from "
      "phi_b=-1.73 across phi=-1 to 0.80. It is a qualitatively different initial state; its endpoint is open.")
claim("MX-11", MX, "With a quadratic tension term tuned to remove the vacuum (d*~-3.1, a model change), shell "
      "radiation reaches 69% of H^2, but Weyl radiation is ~1.8x the radiation when the run stops.", "inconclusive",
      "partially-confirmed", "inconclusive",
      [f"{MX}/M8_QUADRATIC_TENSION_TUNING.json", f"{MX}/verify/V2_PILOT_CONSISTENCY.json",
       f"{MX}/verify/pilot/ext_dstar_Y1_dzf2e-3_phimax1.3_summary.json"],
      "The runs stopped at the solver's phi_b>0.995 cutoff, not the chart freeze. The H^2 budget has large "
      "cancellations (vacuum -117%, Weyl +121%, radiation +69%). With the cutoff lifted, Weyl/radiation is 1.23 at "
      "the end of the reliable window and ~0.6 after the chart freeze, still >=6x above the N_eff limit (<=0.1); "
      "then H turns negative (gauge artefact or reversal).")
# ---------------------------------------------------------------- literature
claim("LT-1", LT, "'Our universe as a wall born in a 5D event' has substantial published prior art (1998-2026); "
      "HDBLAST cannot claim the concept.", "context", "not-independently-audited", "context",
      [f"{LT}/theory.md", f"{LT}/theory_checks/theory_sources.json"], note="Snippet-level only.")
claim("LT-2", LT, "Registered model sits inside standard RS/KR/DFGK formulas: balanced tension = critical tension; "
      "thin-brane H^2 differs from the branch by exactly -c^2 delta^2/384 (7.869 ppm in H).", "exact-verified",
      "not-independently-audited", "exact-verified", [f"{LT}/theory_checks/rs_thin_brane_crosscheck.json"],
      note="Self-checked with 4 wrong-formula controls; no separate auditor.")
claim("LT-3", LT, "The registered PTA knee is at 0.9972 x (1/yr), on the pulsar position/proper-motion blind spot; "
      "no 2025-2026 release tests it at likelihood level.", "numerical", "not-independently-audited", "inconclusive",
      [f"{LT}/observations_checks/knee_and_scale_checks.json", f"{LT}/observations.md"],
      note="Knee position numerical; the blind-spot statement is from snippets.")
claim("LT-4", LT, "Tabletop gravity bounds put the static branch at >=TeV vacuum scale and sub-mHz horizon-scale GW "
      "frequency, ~1600-2600x the knee frequency.", "conditional", "not-independently-audited", "conditional",
      [f"{LT}/observations_checks/knee_and_scale_checks.json"])
claim("LT-5", LT, "Methods sweep ran no new searches (budget exhausted); its list is second-hand.",
      "procedural-negative", "not-independently-audited", "procedural-negative", [f"{LT}/methods.md"])
claim("LT-6", LT, "Certifying the +1 branch needs the eta=phi-1 formulation: the regular cone germ is amplified by "
      "6.288e18 at delta=1e-3 (a direct phi formulation needs ~39 digits).", "numerical",
      "not-independently-audited", "numerical", [f"{LT}/methods_checks/certification_pilot.json"])
# ---------------------------------------------------------------- program map
claim("PM-1", PM, "Program map of HDBLAST Aug 2025 - 23 Sept 2026 with claims ledger; 34/34 self-checks pass.",
      "context", "not-independently-audited", "context",
      [f"{PM}/PROGRAM_MAP.md", f"{PM}/program_map_checks.json"],
      note="Zenodo record 22922928 metadata (version, concept, files) unverified: zenodo.org blocked.")
# ---------------------------------------------------------------- synthesis (this folder)
claim("SY-1", SY, "Methods A and B agree where they overlap: calibration growth to 1.6e-8 (eigenvalue) and 1.6e-7 "
      "(time domain); the artificial flipped-sigma'' tachyon to 3.7e-10; both find no scalar mode below the 9/4 "
      "continuum; B's leading rate matches the continuum edge -3/2 (A's criterion).", "numerical",
      "not-independently-audited", "numerical", [f"{SY}/RECONCILIATION.json", f"{SY}/reconcile_stability.py"],
      note="Comparison of saved outputs; wrong mu^2 maps are rejected as controls.")
claim("SY-2", SY, "Method A's shell coefficient B extrapolates to 3.5555554 as delta->0, matching the exact "
      "32/9 = 14/9 (regular Delta=18 growth rate) + 2 (sigma''/2) to 3e-8.", "numerical",
      "not-independently-audited", "numerical", [f"{SY}/RECONCILIATION.json"])
claim("SY-3", SY, "Cross-workstream reconciliation: the theory sweep's 'no dark radiation' statement applies to "
      "static pure-tension configurations; the dynamical roll-off does produce a positive Weyl term (MX-1).",
      "context", "not-independently-audited", "context", [f"{LT}/theory.md", f"{MX}/M1_WEYL_ALONG_CHAT14.json"])
claim("SY-4", SY, "All 56 key numbers quoted in this checkpoint's report are found in the workstreams' JSON "
      "outputs within stated tolerances; 3 negative controls pass.", "numerical", "not-independently-audited",
      "numerical", [f"{SY}/QUOTED_NUMBER_CHECKS.json", f"{SY}/check_quoted_numbers.py"])


def validate(c):
    errs = []
    if c["status"] not in STATUSES:
        errs.append("bad status")
    if c["verifier_verdict"] not in VERDICTS:
        errs.append("bad verdict")
    if c["established"] and c["verifier_verdict"] != "confirmed":
        errs.append("established without confirmation")
    if c["verifier_verdict"] == "refuted" and (c["status"] != "corrected" or not c["correction"]):
        errs.append("refuted claim without correction")
    for p in c["evidence"]:
        if not os.path.exists(os.path.join(ROOT, p)):
            errs.append("missing evidence " + p)
    return errs


bad = {"id": "CTRL", "status": "numerical", "verifier_verdict": "refuted", "established": True,
       "evidence": ["does/not/exist.json"], "correction": ""}
control_rejected = len(validate(bad)) >= 3

problems = {c["id"]: validate(c) for c in C if validate(c)}
if problems or not control_rejected:
    print("REFUSING TO WRITE:", problems, "control_rejected:", control_rejected)
    sys.exit(1)

counts = {}
for key in ("status", "verifier_verdict"):
    counts[key] = {}
    for c in C:
        counts[key][c[key]] = counts[key].get(c[key], 0) + 1
out = {
    "checkpoint": "HDBLAST research checkpoint, 27 September 2026 (calculations and audits completed 28 September 2026)",
    "author": "Ricardo Maldonado (HDBLAST program); prepared with Claude (Anthropic)",
    "hypothesis": "A higher-dimensional 'blast' (a five-dimensional gravitational event) sparked our Big Bang.",
    "registered_model": {"W": "1-phi+phi^3/3", "U": "(1/2)W_phi^2-(2/3)W^2", "sigma": "2W+delta(1+c phi)",
                         "delta": 0.001, "c": 0.5975949350280132, "bulk": "doubled (Z2), single shell"},
    "headline": [
        "The static +1 branch (the corrected late-time endpoint of 22 Sept) is linearly stable in the scalar and tensor "
        "sectors at all six tested detunings; two independent methods agree and both reproduce the original shell's "
        "tachyon. Status: numerical, independently verified; not interval-certified.",
        "No screened mechanism turns the registered 5D event into a hot, radiation-dominated Big Bang: instant "
        "preheating delivers <=1.7e-4 of the residual-vacuum threshold at the M5 cutoff (conditional), and a 4D "
        "energy budget allows at most ~0.59 e-folds of radiation domination vs ~22 needed (conditional).",
        "Two model-changing leads remain open and unresolved: a quadratic tension term that removes the residual "
        "vacuum (radiation reaches 69% of H^2 but bulk Weyl radiation stays above N_eff limits so far), and a "
        "sign-continued initial-shell family reaching the vacuum-tuned c* (found by the auditor).",
    ],
    "allowed_status": sorted(STATUSES), "allowed_verdicts": sorted(VERDICTS),
    "established_rule": "established=true only if the independent verifier confirmed the claim",
    "validator_negative_control_rejected": control_rejected,
    "counts": counts, "n_claims": len(C), "claims": C,
}
with open(os.path.join(ROOT, "RESULTS_SUMMARY.json"), "w") as f:
    json.dump(out, f, indent=1, ensure_ascii=False)
print(json.dumps(counts, indent=1), "n_claims", len(C))
