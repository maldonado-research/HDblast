# Independent review of the separately registered K384 follow-up

The separate K384 experiment passes all 13 registered gates. The original
K192-limited experiment remains FAIL because its curved pressure cutoff change
exceeded its original allowance. Neither its record nor its acceptance threshold
was rewritten. The original review remains in [ORIGINAL_REVIEW.md](ORIGINAL_REVIEW.md).

## Prospective separation and fixed inputs

The follow-up was committed before K384 execution at
dd1f00fc62a5b2d257f47b8490b79de1514bd2f3. Its protocol SHA256 is
b067afda194895255081553e7bd9a3833b781f5277d3011b619829b4b3033d51;
wrapper SHA256 is
334e0f347df0782765edbbaba891f79f8d8bc48ce9d58275804bd843e8e4c430.
It imports the unchanged repaired core, SHA256
32f6632bb84cc10dd6cef22c2311941acdbe223d846eb7157a2f1e9459fb160a.

Independent pre-execution review verified the imported core, fixed settings,
unchanged gates, frozen baseline hashes, and rejection of mismatched time grids.
The compact baseline's original summary is byte-identical to the preserved
original summary. Every retained t,K192_rho,K192_p,K192_Q array element matches
the original full archives exactly. This lossless export introduces no new
physical calculation. No K768 extension or extrapolated acceptance was allowed
or executed.

## Result and cutoff evidence

| Observable | Maximum K192→384 change | Registered allowance | Result |
| --- | ---: | ---: | --- |
| rho | 0.0000167988 | 0.0007216120 | PASS |
| p | 0.0050869144 | 0.0111218967 | PASS |
| Q | 0.0000152253 | 0.0003376823 | PASS |

The pressure change is about 0.2291% of its comparison scale, below the unchanged
0.5% relative allowance plus absolute floor. The eight/twelve-node quadrature
comparison, half-step physical solver comparison, both historical K192 overlap
comparisons, time ledger, exchange, finite-K trace, future spectral energy,
Wronskians, and new K384 static null all pass.

This supports finite-cutoff numerical convergence evidence for the specified
smooth prescribed background and explicit finite-action convention. It is not
a proof of an infinite-cutoff value or a certified bound on the omitted tail.
The earlier shell shrink ratios were not used as an acceptance extrapolation.

## Numerical floor and remaining qualifications

The larger cutoff increases cancellation sensitivity. The new static null
maxima are 1.745e-6 in rho, 6.00e-7 in p and 1.77e-11 in Q, below the unchanged
2e-5 absolute gate. Halving the physical solver phase step changes rho by at
most 1.425e-6 and p by 5.21e-7. These checks target the new high-momentum band;
the independent Radau comparison remains the original selected k<=48 check.
There was no independent Radau evolution at k=384.

For the selected twelve-node K384 curves, the conformal-energy ledger residual
is 5.29e-6 and the future full-energy/occupation-energy mismatch is 2.38e-6.
The fine local exchange residual is 6.19e-4, compared with 3.08e-4 on the coarse
grid. Both pass, but this is not monotonic improvement under temporal refinement;
the larger cutoff exposes a larger floating-point/solver floor. The independently
sampled pressure trace residual is 1.46e-5. Its direct-equation counterpart is
1.26e-11 and primarily tests algebra/subtraction consistency.

The calculation retains raw numerical modes without post hoc normalization.
The maximum raw Wronskian error across the four new cases is 8.77e-14.
Future coherent pressure and Q remain retained; maximum differences from
occupation-only estimates are approximately 0.05055 and 0.02348. A particle
energy curve is therefore insufficient to replace the full paired source.

The action-matching qualifications from the original review remain unchanged:
finite local counterterms define a choice, FRW does not determine the independent
Weyl/Euler extension, and the heavy-field argument has no supplied certified
uniform remainder constant. No coupled shell, backreaction, energy-consumption
model, decay or thermalization was solved.

## Independent post-execution and CI review

The read-only artifact audit passes for all four new cases, in addition to the
13 original cases. It loads only numeric NPZ with allow_pickle=False, checks
matching physical/observable grids, finite arrays, saved-mode reintegration,
raw Wronskians and recomputed future occupation energy. Fifteen cutoff, solver,
quadrature and historical overlap comparison records were separately recomputed
from saved curves; their differences, scales, limits and verdicts agree exactly.

The curated original-outcome validator preserves scientific_status FAIL while
returning EXPECTED_SCIENTIFIC_FAILURE_CONFIRMED for the documented negative
result. Independent synthetic tests show it rejects rewriting the original
status as PASS and rejects an additional failed final acceptance gate. It does
not promote report-only coarse diagnostics to new scientific acceptance gates.
The replay wrapper records the actual original command exit, validates that
negative outcome, and evaluates the K384 result separately using the pinned
portable historical baseline. Package/replay validation success must remain
distinct from the original scientific outcome.

Evidence: [follow-up artifact audit](outputs/followup_integrity.json),
[follow-up summary](outputs/followup/summary.json),
[preserved original summary](outputs/original_matrix/summary.json), and
[unchanged independent mode evidence](outputs/independent_radau.json). This is independent internal
implementation review, not external peer review or support for a discovery claim.
