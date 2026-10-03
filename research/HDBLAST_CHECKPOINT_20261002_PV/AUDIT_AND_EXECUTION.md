# Internal audit and execution

Internal AI-assisted analytical and numerical checks are not external peer review.

The original protocol was committed before mode integration at 32f1fc617b0691bbc10783e7e12b337afb4a3af9. The pre-execution addendum and initial source were saved at 06b92943aad8a112879180d11efd6369f394320e. Both protocol files are hashed and checked by code/verify_protocol.py. No physical numerical outcomes were known when those choices were fixed.

The first executed scientific source commit is 9829992e2389f4195e7f062b34f731bbf56ecfaa. [Run 36949788139](https://github.com/maldonado-research/HDblast-archive/actions/runs/36949788139), job 110659940363, completed successfully. Compilation, protocol hashes, independent 70-digit matching/state tests, all eight dynamic cases, the algebraic static null and artifact upload passed. Source SHA-256 is recorded in outputs/pv_actual.json. The final public replay status is linked in the publication PR.

Independent mathematical reviewers checked the common-action variation, curvature pressure, UV coefficients, exact Jost series and out projection, analytic phase, stable static-tail formulas, matched trace including finite-K correction, and the direct positive-reference representation. No blocking sign/factor error was identified. Code was also read independently.

The high-precision checks use mpmath at 70 decimal digits, independent analytic zero-point integrals, differentiation of the analytic tail, hypergeometric functions, and direct phase quadrature. The largest recorded absolute discrepancy is 1.14e-13 in the phase; static functions agree within 3.91e-14. The test covers 54 static combinations and additional state/phase points.

Registered outcomes: maximum exact-occupation discrepancy 5.74e-10; maximum Wronskian drift 1.37e-9; maximum normalized work residual 5.59e-14; maximum nine-point pressure-trace residual 2.34e-11 in units of m_infinity^4. The trace uses independent finite differences and reports stencil/coarsening sensitivity, rather than differentiating via the mode equations. Work is an auxiliary quadrature sharing the ODE integrator; it is not a separate solver.

The work normalizer includes the natural floor m_infinity^4=16. Raw absolute residuals divided by that scale and source magnitudes are also provided. A natural-scale tolerance is not a relative error on a small source.

The estimator uses the exact normalization relation to reconstruct bilinears. Conservative integrated differences from raw-mode reconstruction are recorded. Primary bounds are 1.69e-7 for energy, 5.37e-8 for pressure and 9.46e-12 for S; at K=192 the energy bound is 1.61e-6. These are diagnostic bounds on the reconstruction discrepancy, not rigorous bounds on every integration error. Small quadrature differences cannot establish physical source accuracy below these scales.

The finite regulator sequence passes its declared natural-scale target, but its Lambda=8 to 16 changes are still roughly 5–7% of the crossing stress and about 4.5% of the crossing current. Regulator removal is not mathematically certified. Direct-limit cutoff and trace comparisons are retained separately in outputs/analysis_summary.json.

Figures and further comparisons are computed from the saved eight-run arrays. This postprocessing does not rerun or retune the experiment. No failed variant was omitted and no threshold was relaxed after the results.

Saved-data analysis and three Matplotlib SVG figures were executed successfully in [run 36950532920](https://github.com/maldonado-research/HDblast-archive/actions/runs/36950532920), job 110662258303, at commit 4ee128225620adb3c0d7b2eaa131100ed25175ff. This reused the first run's arrays. The curated analysis_summary.json preserves that output.

Two reviewers independently derived the exact finite-K direct trace factor and leading omitted UV tails after these results. Root and adversarial reviewer each executed a nine-point V8 finite-difference check on the saved primary curves, agreeing on maximum absolute residual 1.6729e-9 when the factor is included. outputs/posthoc_tail_checks.json records that check and the crossing addbacks. The updated analysis script repeats the exact finite-K trace check for all eight cases during the final replay. These are labeled post hoc; no prospective thresholds were changed.

Full regenerated per-case NPZ arrays are available as Actions artifacts with finite retention, and are reproducible from the public source. The permanent curated package preserves all case summaries and spectra, the 2001-row primary source history, analytic checks, analysis summaries and figures.
