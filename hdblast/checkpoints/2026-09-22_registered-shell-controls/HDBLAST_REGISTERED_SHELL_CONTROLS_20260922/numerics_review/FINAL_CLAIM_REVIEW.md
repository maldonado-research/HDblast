# Final scientific-claim review

22 September 2026. Fresh review of the root `00_READ_FIRST.md` and `CONTROL_ANALYSIS.md`, compared with `CONTROL_ANALYSIS.json`, the independent solver/constraint/seed reviews, saved balanced-candidate results, the primary-literature report, and public Zenodo metadata. This review did not edit the root reports, rerun an evolution, recompute spectra or shooting roots, or certify physical endpoints.

**Verdict: no blocking scientific-scope or novelty overstatement found.** The reports make a concrete but bounded claim: numerical recovery of a previously known registered scalar growth rate, short small-amplitude checks of the nonlinear equations, exact diagnostic identities, and candidate initial data. They clearly leave registered nonlinear fate, a successful coordinate repair, visible matter production and observational validation open.

## Checks of the presentation

- The quoted spectral and fitted rates agree with the machine-readable control analysis. The target rate was known in advance and is explicitly labeled a validation target. The spectral difference is not presented as an interval error bound or proof of full spectral stability.
- “Six spectra” and “six short evolutions” match the saved analysis. The two shifted-tension runs are preserved as negative initial-data controls because their constraints worsen with refinement. They are not silently grouped into a successful nonlinear-convergence claim.
- The constraint-transport signs, factor `exp(3A)`, outgoing direction and approximate gain agree with the independent audit. The README includes the necessary continuum/no-dissipation qualification and avoids attributing every failed earlier evolution to this mechanism.
- The initial-data defect identity and second-time junction obstruction agree with the independent symbolic derivation. Exact local formal compatibility is conditional on a genuine compensated root and an open bump-free brane neighborhood. Numerical roots and rounded residuals are not existence certificates.
- Constraint normalization is disclosed as a background scale, not a relative perturbation-error bound. The independent solver review separately records the larger one-sided field-derivative residuals and raw eigenvector velocity defects.
- The proper-clock result is described as a local coordinate identity, with positive Jacobian required; no PDE implementation or additional spacetime coverage is claimed.
- The current model remains the homogeneous, matter-free Einstein–scalar system with a spacelike extra dimension. Adding visible matter, modified gravity or extra fields is correctly identified as a new model branch.

## Small correction sent to the root author

The reviewed README stated that floating direct `D_b` and higher-corner residuals “remain nonzero.” Two saved `independent_mass_transport` entries instead round to `D_b=0.0` and signed-zero corner residuals. Suggested wording:

> In floating arithmetic these quantities are not certified zero; they are generally nonzero, and a rounded zero is not an exact certificate.

This correction does not change the mathematical conclusion or the candidate status. The independent seed review's eight saved-candidate checks concern the four signed amplitudes at two tolerances; they should not be described as an independent rerun of all shooting formulations. Its exact compatibility derivation is broader than those numerical checks.

## Primary-literature and publication scope

The main report's prior-art comparison is supported by [BraneCode](https://arxiv.org/abs/hep-th/0309001), while retaining the two-brane versus single-shell distinction. Its statement that reheating need not require oscillations is supported by [Instant Preheating](https://arxiv.org/abs/hep-ph/9812289), and is properly conditional on introducing the missing matter sector and coupling.

The recent papers in the detailed review require assumptions absent from the registered model. They are not evidence that its collapse becomes a bounce or that it produces a thermal bath. September journal dates are distinguished from earlier preprint dates. The public [Zenodo record 22347452](https://zenodo.org/records/22347452), verified through its API, is the 5 September v24 conditional linear-response checkpoint, not publication of the later shell calculation.

No assertion of external mathematical novelty, a verified Big Bang mechanism, independent external peer review or a prize-level discovery appears as an established conclusion. The source-v2 snapshot caveat and bounded file/web audit are appropriate; a changing external source must remain distinguishable from frozen inputs.

## Remaining acceptance gates

The next claim should depend on demonstrated constraints for the initial data actually sampled on the evolution grid; finite-amplitude and resolution controls; adequate resolution of the departing wall; and independently benchmarked late-time coordinates. A heating claim additionally requires specified and normalized matter couplings, energy exchange and backreaction, a demonstrated thermal distribution or justified thermalization treatment, radiation domination and controlled residual vacuum/bulk energy. None of these gates is closed by the current growth-rate agreement.

Release resolution: the root report adopted the rounded-zero correction. Final replay also extended the saved-candidate consistency checks from eight to ten after addition of the zero reference; the coverage remains nominal/refined data, not an independent shooting replay.
