# Independent review of fresh general-family benchmark

8 October 2026. Read-only source/results audit; no full sweep rerun. Internal AI-assisted review, not external peer review or an interval proof. The recorded script SHA256 matches the reviewed file. See independent-family-review-summary.json for recalculated coefficient diagnostics.

**Verdict:** the implementation and recorded finite-detuning trends support the conditional second-order scalar-feedback law. No equation or cone-series coefficient error was found. This is stronger than residual-only checking, but three finite positive detunings do not prove an asymptotic remainder theorem or existence.

The dimensionless variables are t=k*y and R=k*rho. Accordingly R_t=rho_y, eta_t=phi_y/k, the second-order geometric equation is R_tt=-R*(eta_t²/4+U/(6k²)), and both junctions retain their stated k factors. U_phi=W_phi*(W_phiphi-4W/3), U_phiphi=W_phiphi*(W_phiphi-4W/3)+W_phi*(W_phiphiphi-4W_phi/3) are correctly differentiated. The regular-cone R cubic/quintic and eta quadratic/quartic coefficients are the proper scaled forms. The Riccati linear seed has the correct positive-cone leading logarithmic derivative. It initializes shooting and does not impose the nonlinear response coefficient on the integrated solution.

All six parameter sets have k>0, w>4k and f0>0; all tested detunings are positive. They therefore probe the assumed growing, large-radius positive-curvature regime. The code solves both second-order Einstein/scalar equations, rather than building the radial constraint into the integrator. The constraint and boundary-curvature identity supply distinct useful checks, although they share the model functions and floating-point implementation. Constraint samples (257) are not continuum enclosures; adaptive tolerances and setting agreement are not certified global errors.

For each nonzero-coupling family, define Q(delta)=(H²-H²_metric)/delta² and q=-k*f1²/[24(w-2k)]. The recalculated |Q-q| halving ratios range **1.99956–2.00043** over .002→.001→.0005, closely supporting O(delta) coefficient error. At delta=.0005, relative coefficient discrepancies are:

| Family | Relative discrepancy |
|---|---:|
| Registered HDBLAST | 0.12496% |
| Changed cubic/tension | 0.09025% |
| Negative scalar coupling | 0.02643% |
| Small bulk curvature | 0.19090% |
| Stronger tension curvature | 0.04400% |

The corresponding scalar-shift errors also halve (ratios 1.99867–2.00206). Third-order normalized H² remainders remain approximately constant, consistent with the stated truncation. The largest standard/refined or identity difference divided by the nonzero scalar-feedback signal is below 2.2e-7 across these cases. Thus the observed feedback is resolved well above these measured numerical differences. This ratio is evidence from the compared settings, not a complete error budget.

The f1=0 control retains eta=0 exactly in the implementation; f2 need not vanish for this constant solution. Its H² difference from the exact metric formula is about 1e-19, a floating-point floor. Dividing that floor by delta²/delta³ magnifies it; do not demand an O(delta³) trend or quote relative errors for a zero coefficient. The control is a useful geometry/unit check, but its scalar value is deliberately fixed and it is not an independently fitted scalar root.

Actionable qualifications:

1. Call these **36 radial integrations / 30 fitted shooting integrations**, or 18 configurations at two settings. The six zero-coupling runs use the analytic shell position and eta=0 rather than a fitted root. All recorded success flags and residual checks pass, so this counting distinction changes no result.
2. The numerical_checks_pass threshold tests finite residuals, not agreement with the susceptibility law. Keep the trend diagnostics above as a separate scientific check; a pass flag alone does not validate the law.
3. No matched pair currently holds k,w,f0,f1 fixed while changing W'''=v or f''=f2. The changed-cubic case changes other parameters too. State that exact algebra establishes their absence from second order, while these data demonstrate consistency across varied models. To claim a direct numerical independence test, add matched variations and compare their delta³-scaled differences at multiple detunings.
4. Retain the qualified status: no negative-detuning branch was solved; no convergence radius, interval existence/uniqueness, stability, dynamical attraction, external novelty or physical confirmation follows. A final residual pass can be valid even if a root solver reports failure, but such cases must be disclosed; none is recorded here.

The remaining numerical limitations do not invalidate the present benchmark. A focused research note can accurately report the exact local derivation plus these finite-detuning tests, with the above distinctions explicit.
