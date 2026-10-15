# Hash-specific Sol review of the final waveform checkpoint

15 October 2026 UTC. **Accepted within the stated conditional scope.**

The reviewed main note is `UNIFORM_FIRST_VARIATION_WAVEFORM.md`, SHA-256
`6eada42c742aec6141925c244df21dc263e3555ae2a064e1e5d3fb7dcc9a82ad`.
`FINAL_WAVEFORM_REVIEW.json` pins this note, the root evaluator, generator,
checker, computed CSV/JSON, and the supporting reviewed inputs. This
acceptance does not automatically apply to changed hashes.

The independent Sol characteristic and positive boundary-kernel proof
controls the actual derivative observation. Its exact global error, tail
error and support-aware observation-norm bound equal those used by the root
calculation. The master is the retained coupled scalar master with both
prescribed initial profiles; the dynamic shell law, the factor
`gamma=C_b/2`, and the negative signed physical gain are retained.
The prefront zero interval and the bounded post-pulse tail are properly
distinguished. The main note preserves the strict causal/evolution margins
and the conditional physical gauge and exact-family hypotheses.

I reviewed the continuous-cell range argument independently of its plots.
On each support cell, reciprocal bump coordinates have exact extrema;
both `exp(4-v)` and `v^2 exp(4-v)` decrease for `v>=4`. Interval products
therefore enclose the signed derivative, including cells touching the flat
support endpoints. Outward negative exponential evaluation uses a reduced
Taylor remainder, repeated interval squaring and a conservative large-input
tail. Signed multiplication preserves the negative gain. The clock-factor
bound follows from `log cosh(t)<=t^2/2` and
`exp(-3t^2/4)>=1-3t^2/4>=1-t^2`; the exported interval is outward.
Real continuous-cell sums enclose the comparator norm. Square-root
enclosures, the observation-norm triangle inequality, positive prefactor
and `1/ell` Jacobian then give the stated leading coefficient interval.

The peak statement is also consistent with a uniform observation error:
the independently replayed exterior comparator gaps exceed twice the Sol
published error plus `ell`. The main note's outward height intervals hold
using that same Sol bound. This supports exclusion of actual extrema
outside the stated open boxes, without promoting comparator critical-point
uniqueness to actual-waveform uniqueness or identifying the absolute-peak
lobe.

Normal and optimized checks were regenerated after changing the Sol input
path to capsule-relative `inputs/previous`; scientific values are unchanged.
The root finite checker, independent CSV/norm auditor and peak checker
pass in both modes with identical parsed results. Fresh normal/optimized
4096-cell generations reproduce both committed computed files byte for
byte. Five explicit Sol failure injections each fail at their named gate
in both modes; a failed injection writes no passing artifact. All replay
results and exact hash pins are recorded in `FINAL_WAVEFORM_REVIEW.json`.

The acceptance covers the narrow specified first variation and leading
response coefficient. It does not provide an effective nonlinear epsilon,
a finite-amplitude waveform, observed time calibration, external peer
review, novelty determination, or formal verification of the retained
analytic PDE/evolution hypotheses.
