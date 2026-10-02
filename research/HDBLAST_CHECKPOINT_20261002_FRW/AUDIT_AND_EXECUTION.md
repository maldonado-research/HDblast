# Execution and independent internal audit

## Prospective records and actual runs

Registration commit: 34a2be98b5e972ce6fa0263cda201bae56488720, SHA-256 939e7f6ad178c32cdcf638d01b8694d2194e805f0aad38964416518e3b72a2b5.
Pre-execution addendum commit: 2dee11d88753dabe1b07e1ab4a7cc500075d676c, SHA-256 e7a98bea92852ac47a57a720e16efadd48444283ae71bbe72c5daa99575337c4.
Neither document was rewritten after results.

| Run | Source | Actual outcome |
|---|---|---|
| [36957517031](https://github.com/maldonado-research/HDblast-archive/actions/runs/36957517031) | 20a60ab48c8b28eb106696866ba980ecfff45d35 | Compilation failed; no physical calculation executed |
| [36957629852](https://github.com/maldonado-research/HDblast-archive/actions/runs/36957629852) | 0ee14880ffa04b938bda62cb54fad26ba2c70b93 | Algebra and oscillator controls passed; four converted higher-degree continuity checks failed |
| [36959579410](https://github.com/maldonado-research/HDblast-archive/actions/runs/36959579410) | 2bad8db9b1531732b917a7783e45e8257b0f7478 | All eight native-reconstruction checks and all controls passed; comparisons, smoothing limit and figures executed |
| [36959791142](https://github.com/maldonado-research/HDblast-archive/actions/runs/36959791142) | 84d33b785edf351426adc083510c4c4c923445cd | Same controls passed; added explicit nodal fidelity and native/de Boor comparison |

The curated primary geometry/figure report comes from the first successful repair run. Fidelity comes from the supplementary run. The sources differ only by the added diagnostic/workflow; the geometry algorithm is unchanged. Later packaging replays verify the checked-in ZIP and regenerate controls on the PR head; their links are recorded in the PR description to avoid embedding a commit's own hash inside itself.

## Representation repair

See NUMERICAL_REPAIR.md and the complete original failed geometry report. Native B-spline coefficients/knots now govern every field evaluation; piecewise polynomial conversion is used only to suggest roots, which are refined and checked natively. One-sided knot limits use cached 80-digit derivative coefficients and de Boor recurrence. True higher-derivative jumps remain. Finite-value checks cover fields, side limits, conformal-time quadrature and tail coefficients. No theoretical zero is substituted for an evaluated continuity residual and no threshold is relaxed.

## Implementation differences and post hoc additions

- Registered stored-value interpolation fidelity was completed explicitly by the supplemental diagnostic; derivative fidelity was already present. No new pass threshold was introduced.
- Conformal-time quadrature partitions at each reconstruction's own breakpoints; for Hermite these are the Hermite intervals. This avoids integrating across a piece boundary of the chosen fit.
- Initial fixed-band interference output is a separate two-jump toy. The later 299-knot fine-archive check uses inverse minimum conformal spacing with explicitly post hoc scale factors. Both evaluate leading asymptotic amplitudes, not archived exact modes.
- Reconstruction discrepancy ratios in geometry_audit.json use the **maximum absolute value of the primary Hermite curve**, not peak-to-peak range. Posthoc C4/C6 ratios use the maximum absolute C4 sample. Stored-derivative ratios use the full reliable history's maximum absolute stored derivative, not only the audit window. These denominators are retained and disclosed rather than relabelled as range normalizations.
- Full-grid C4/C6 comparisons, the smooth-width energy quadrature, nodal/native diagnostics and the final archive-interference implementation were added after primary outcomes. They are descriptive/supplementary; they do not alter registered criteria.
- Finite output JSON uses null only for archival unavailable metadata. Physical field, jump and quadrature nonfinites cause errors. Original input summary JSON preserves late-failure NaN bytes and is read by the original loader.

## Independent internal review

Separate agents derived/reviewed the FRW counterterms, wrote a real Radau oscillator implementation with an analytic Jacobian, screened other projects, checked the repair and audited actual output logs. A separate JavaScript Hermite formula check predates the Python audit. The reviewer found no blocking mathematical/implementation defect after repair, requested explicit nodal fidelity, and emphasized that the largest C4/C6 high-derivative differences occur at the initial boundary. This is AI-assisted internal review, not external peer review or an interval-certified proof.

## Privacy and scope

The four input files are scientific inputs already included in the prior public quantum checkpoint; their original bytes and provenance are retained. Only newly prepared scientific reports, executable code, computed outputs, source metadata and these curated inputs are published. No private new-files raw archive, chats, correspondence or personal documents are included. Public sister-project methods are cited, not copied wholesale. Drive handoffs informed context only.

The verdict is numerical-control success with physical readiness blocked. Neither an absolute FRW source, an admissible global quantum state, coupled backreaction, a thermal radiation bath, a hot Big Bang nor external mathematical novelty is claimed.

The reviewer independently executed the 299-knot archive interference sum with compensated accumulation, compared three- and six-pair Ci expansions, and derived/executed the phase-independent inequality. The ratio agrees with the primary JavaScript result within floating roundoff. The bound concerns the leading-amplitude model with supplied floating coefficients; no EFT validity at the displayed large momenta is assumed.

Supplementary Python check code/verify_archive_interference.py uses scipy.special.sici, compares the independently executed JavaScript ratio, and verifies the phase-independent relative bound. Its executed outcome is recorded in the final PR CI logs; the frozen JavaScript report remains in the package.
