# Independent review of the unchanged-candidate initial refinement

30 October 2026. This review concerns the initial-data error certificate for
the existing `W1_quadratic` reconstruction. It does not change its initial
data, the actual target, the full original quadratic source, the shell law,
the outer quadratic continuation, or either convolution width.

## Why subdivision can improve a valid certificate

Write the registered paired comparison as three telescoping channels
`actual-nominal=sum_ch d_ch`. On any closed leaf `I`, interval evaluation
of the order-five registered algebra supplies component caps for scaled
`u_rr` and `v_r`. The square root of the positive weighted sum of their
squares bounds the H norm of that channel at every point of the leaf.
Therefore, for a finite exact cover `I_j`,

```
sup_r ||sum_ch d_ch(r)||_H
  <= sum_ch max_j channel_H_cap(ch,I_j).
```

This is a statement about a maximum over a cover, not an assumption that
the independently chosen channel maxima occur at the same point. Triangle
inequality permits the sum of maxima. The older aggregation took component
maxima before packing the H norm; the displayed aggregation is also valid
and can be sharper. Replacing one coarse interval with smaller exact
covering intervals leaves the analytic comparison unchanged.

The fixed nominal displacement is C1, so its distributional second
derivative is its piecewise classical second derivative without point-mass
terms. A positive even convolution kernel contributes the retained global
charge

```
M3_u * (35 delta_u/128) + (maximum isolated curvature-jump cap)/2
```

componentwise. The outer quadratic continuation matches curvature, so it
adds no outer jump. The first inner join at `65/64` has zero nominal
curvature jump in every component. The second inner join is retained.
The separation of joins exceeds twice `delta_u=2^-30`; a kernel footprint
contains at most one join. No derivative of an unknown Taylor remainder
is taken: each needed background derivative has its separately enclosed
order-five remainder in the authenticated source core.

The velocity remains the same rounded-node piecewise linear profile,
smoothed at `delta_v=1/2048`. Its secant, endpoint-radius, and convolution
charges are retained. Improving the paired actual-minus-nominal comparison
does not remove those charges.

## Numerical audit design

The independent checker imports no refinement producer code. It hashes the
unchanged retained `initial_second_core.py`, authenticates its original
source closure, recomputes each new leaf's component caps using the original
registered algebra, and checks the report's rational H caps directly
against the positive weighted squares. It verifies exact active-domain
coverage, that only the original cutoff cell is subdivided, and equality
of all retained per-cell caps to the authenticated original report.

It then reconstructs the same inner-jump, third-derivative convolution,
rounded-node secant and velocity convolution charges and the final uniform
H bounds. Structural and arithmetic negative controls must reject omitted
leaves, changed retained cells, inconsistent scales, incorrect Gram weights
and omitted smoothing charges. Authentication is repeated after the audit.

This is an independent aggregation and implementation audit built on the
same previously verified analytic core, not an independent derivation of
every upstream Einstein-scalar or initial-data premise. Numerical checks
and exact rational arithmetic do not establish physical relative accuracy,
small covariant residuals, a finite-amplitude nonlinear family, or a
higher-dimensional origin of the Big Bang.

## Result

Both normal and optimized isolated Python executions passed, with
byte-identical reports. Each recomputes all ten new leaves, verifies all
47 retained cells, and verifies the 57-cell exact active-domain cover.
There are 374 explicit check entries and seven rejected negative controls
per mode. These are internal reproducibility and audit counts.

The independently reconstructed uniform initial displacement-curvature
error cap is approximately `1.1845751478843721e28`, strictly below
`1.185e28`. The old upper cap was approximately `6.218904973049681e33`;
their ratio exceeds 500000. This is a ratio of two upper certificates for
the same analytic candidate and target, so no reduction of the actual
error is inferred. The smoothing charge remains approximately
`4.131484717143961e27`. The full velocity-gradient cap is retained
exactly at its published value, below `9.672e24`; no small paired velocity
term is substituted for the full interpolation-plus-convolution ledger.

The producer's proof note and source were also read directly. Its exact
leaf coverage, ordinary derivative scaling, isolated-jump charge and
quadratic exterior treatment agree with this review. No unresolved
mathematical or implementation blocker was found within this audit scope.
