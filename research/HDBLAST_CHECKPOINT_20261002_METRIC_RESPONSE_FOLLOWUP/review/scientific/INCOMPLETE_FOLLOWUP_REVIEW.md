# Read-only scientific review of the incomplete continuum followup

The separately registered continuum followup at public freeze
`5ff571cdb7729398834ff3c025b50191c4ee29f2` **failed before scientific
validation**. The replay attempted 28 of its 32 commands, with 27 successful
commands. Both frozen physical producer commands were attempted; only the
primary producer completed. Source verification passed before and after the
replay. The failed registration, partial archives and logs must be retained
with this outcome. This reviewer read saved evidence only, ran no physical
evaluation, and changed no numerical source, configuration or gate.

The direct-mode producer saved the positive_B coarse and fine runs, then its
post-run resource check rejected peak RSS of **306088 KiB** against the frozen
**262144 KiB** budget. Its internal elapsed time was
153.80774093999935 seconds; the wrapper recorded 153.93690977899678 seconds and
exit code 1. The exception was `RuntimeError: Frozen 256-MiB peak-memory budget
exceeded`. Neither signed_uB run was attempted. The producer failed before
combining the four runs, checking its complete refinements/internal gates, or
entering either normal or optimized cross-route/raw-archive validator.

The primary producer completed all twelve source/time rows in
37.52687140800117 seconds. Its saved result contains 432 normalized quantity
error estimates across continuum and three finite cutoffs. The maximum saved
estimate is 1.9193656658868163e-12 for positive_B, eta=-2.5, K=64,
q_second; the maximum saved density-derivative estimate is
3.8640698661404397e-13 at the same source/time/cutoff. These are empirical
producer diagnostics, including fixed 50/70-digit continuum comparisons and
quadrature estimates. No completed independent validation establishes their
accuracy. Completing the primary calculation at the original failure point
does not make the full experiment pass.

There is a separate warning in the incomplete saved direct-mode metadata:
the declared fine positive_B Ward endpoint maximum is
**1.9552444892255006e-05**, at eta=-1.5, K=256, against the unchanged 2e-06
fine endpoint gate. The coarse declared maximum is 0.02107360235728546, at
eta=-2.5, K=256. These are existing producer declarations, not raw-audited
measurements or a newly executed official gate result. The maximum declared
canonical Wronskian residuals are 1.0579640034553224e-17 (coarse) and
2.477612729749456e-17 (fine); small Wronskian diagnostics alone do not establish
accurate stresses or Ward histories. An archive-only resource repair should
preserve all these computed values and must carry no expectation of a
scientific PASS.

The reviewer compared the original and continuum-followup registrations.
Geometry, sources, initial/final/observation times, cutoffs, normalization,
independent configuration, runtime, all scientific gates, Ward acceptance,
and continuum comparison policy are exactly equal. Eleven checked source
files, including the direct-mode producer, generic WKB, stable baselines,
metric contact coefficients, source jets, tail module/certificate and both
scientific validators, are byte-identical. Only the primary continuum
algorithm/metadata and followup description change. The machine-readable
review records all comparisons and pins the actual saved evidence.

The original public freeze
`57668b8fadd75df8738565e0bbd1eb852c1ebae8` remains a distinct **FAIL**:
its primary stopped at positive_B, eta=-4.5 when SciPy's weighted-log
quadrature raised the roundoff warning that the original registration made
fatal. It attempted 27 commands, with 26 successes; no independent physical
run began. The continuum followup addresses that integration implementation,
while its own later memory failure is preserved separately. A future expected
failure control must not relabel either scientific result as a PASS.

The action/operator evidence remains a useful analytic prerequisite. At the
special point H=1, x=r=2, xi=0, the physical homogeneous conformal perturbation
has canonical forcing `g=4 L^2 h-2 L h'-h''=-a0^2 deltaR/6`. The measured
variance/stresses require the full metric operator, measure and W2/W4
subtraction contacts; current=q for fixed scalar source. The fixed-K physical
Ward identity retains `3 h' a0^4 (rho0,K+p0,K)` and is independently tested
through direct history integration, rather than defining either stress. The
producer and raw auditor retain unprojected complex modes, actual finite-K
baselines and distinct algebraic groupings. None of these implementation
checks substitutes for the unexecuted scientific comparisons.

The directed tail module encloses only the omitted k>K band. Its complete
N2 bound uses g through its fifth derivative, hence h through its seventh;
the uniform bounds are deliberately loose. The failed replay never reached
its post-run continuum comparisons. Even a completed comparison would not
certify finite quadrature, mode/time errors or precise continuum pressure/sign.

The actual shifted-root propagators and full metric contact kernels remain
to be calculated. A general lapse plus scale response or explicit coordinate
pullback is needed for compact gauge controls: a conformal h alone is a
physical curvature perturbation. Bulk/boundary completion, state and
constrained initial data, coupled evolution, stability and particle/heating
analysis remain independent prerequisites. This failed prescribed-background
calibration establishes no Big Bang origin mechanism or discovery.

`INCOMPLETE_FOLLOWUP_REVIEW.json` contains exact evidence hashes, the source
invariance table, incomplete-run declarations and the uncompleted checks.
