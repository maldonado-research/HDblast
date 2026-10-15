# Astra final cross-review of the uniform waveform checkpoint

15 October 2026 UTC. **Conditional acceptance; no blocking finding.**

Reviewed main note SHA-256:
`6eada42c742aec6141925c244df21dc263e3555ae2a064e1e5d3fb7dcc9a82ad`.
The companion `FINAL_ROOT_WAVEFORM_REVIEW.json` records exact pins for all
reviewed production, computed, graphic, gain, and auxiliary proof files.

The characteristic proof closes. Positive shell inertia gives positive
boundary damping, and integrating the known bump derivative by parts bounds
its boundary-velocity contribution by 1/mu. The displacement/velocity matrix
has positive inverse. The incoming characteristic then controls the observed
normal derivative directly. The stated global 1.0000004e-7 bound, 4.752e-14
tail bound, and 1.0000005e-7 observation-norm bound follow from the exact
rational majorants. The independent Astra error-field derivation agrees and
is slightly sharper. Neither proof infers pointwise control from bulk energy.

The source region and boundary lift stay within the retained collar. The
zero prefront claim concerns the linear first variation, and the tail bound
retains the earlier boundary history. The nonlinear correction margin,
compatible smooth data, reconstruction/gauge identification, exact-family
hypotheses, and missing effective nonlinear amplitude radius remain explicit.

The signed physical gain is negative on the actual retained branch. Its
correlated shell-constraint audit and the retained output formula imply the
claimed physical sign reversal. The proper-time cosh factors cancel exactly
in the response integral, with the rescaled-time Jacobian 1/ell included.
No observed dimensional time or finite nonlinear response is inferred.

The continuous-cell production formulas have been inspected, including
outward signed integer arithmetic, degree-24 exponential Taylor bounds and
repeated squaring, the large-argument tail, monotonic reciprocal bump ranges,
and flat support endpoints. The separate cross-review script imports no
producer: it checks the exported grid, all 4098 error additions and physical
interval products, squared-cell integral sums, outward square roots, norm
propagation, and K/ell response scaling with exact arithmetic. It authenticates
the CSV, gain input, and producer hashes. The production checker was also
replayed in normal and optimized Python, with identical passing receipts.
This finite arithmetic audit is distinct from the analytic cell-range proof.

The peak statements are appropriately limited to the normalized trace H.
Exact comparator geometry and value-gap bounds place every actual global
maximum in (1.3033,1.3034) and every actual global minimum in
(1.6966,1.6967). The stated rounded heights safely contain the independent
certified intervals. No actual zero/extremum uniqueness or identification of
the absolute-maximum lobe is claimed. High precision auxiliary critical
points are explicitly comparator quantities.

The PNG was visually inspected. Its negative-then-positive physical pulse
matches the negative gain. Its band is generated from the certified cells;
the centerline is labeled illustrative, the axis is model time, and the
tail panel explicitly displays an unresolved shape bounded by a small band.
The scientific-limit caption is legible and consistent with the main note.
The approximate centerline is display-only and is not used in any certificate.

The rounded response claim 6.06e15 < C < 7.44e15 encloses the exact computed
endpoints. This is a tighter conditional result for the same state and
normalization, not an estimate for a wider pulse or a calibrated physical
event. Internal AI cross-review and executable rational checks remain
distinct from external peer review or formal verification of the analytic
PDE and retained exact-family construction.

Reproduce the finite cross-review after the two root checker receipts exist:

```sh
python -I -B astra/review_root_waveform.py
python -I -O -B astra/review_root_waveform.py
```

Root files were not edited during this review.
