# First minimum passage in the coupled closed universe

On the fitted potential with `phi_m = −0.5`, all four time-symmetric quadrupoles develop a local minimum crossing while the total proper volume is expanding. The sign of the initial quadrupole determines which part of the sphere crosses. Positive seeds first cross at the two poles and later cross back on those same geodesics. Negative seeds first cross at the equator after the poles have begun contracting. On the `phi_m = +0.5` potential, the flank rolling away from the minimum reaches the edge of the fitted field range first for all four seeds.

All times below use the pinned solver's dimensionless Gaussian proper time. A dash denotes an event absent by the domain guard. The local contraction time is the first zero of the minimum of the normal expansion `H_local = theta/3` after the initial time-symmetric slice. The volume turnaround is the positive-to-negative crossing of `H_V = d(log V)/(3 dx)`. A local expansion sign concerns the trace of the expansion tensor.

| phi_m | Signed A | First minimum | First local contraction | First return through minimum on entry geodesic | Volume turnaround | Domain guard | H_V at guard |
|---:|---:|---:|---:|---:|---:|---:|---:|
| +0.5 | +0.0010 | — | 0.341599 | — | — | 0.370300 | +0.291900 |
| +0.5 | −0.0010 | — | 0.426248 | — | 0.449844 | 0.450842 | −0.024716 |
| +0.5 | +0.0005 | — | 0.395123 | — | — | 0.421088 | +0.335717 |
| +0.5 | −0.0005 | — | 0.479288 | — | — | 0.501847 | +0.017980 |
| −0.5 | +0.0010 | 0.488395 | 0.611053 | 0.642846 | 0.672399 | 0.751516 | −0.541783 |
| −0.5 | −0.0010 | 0.589134 | 0.500840 | — | — | 0.648388 | +0.107950 |
| −0.5 | +0.0005 | 0.551851 | 0.680019 | 0.706671 | 0.743866 | 0.816778 | −0.495716 |
| −0.5 | −0.0005 | 0.653080 | 0.570531 | — | — | 0.713378 | +0.151375 |

## What the turning geodesic does

The scalar offset is `q = Theta − Theta_top`. The fitted minima are `q_m = −0.385805180` on the positive branch and `+0.296939358` on the negative branch. On the positive branch, the deepest approach to the minimum by the guard is `−0.337189` for the larger negative seed and `−0.336128` for the smaller one. Positive seeds reach only about `−0.0745` on that side. The upper domain limit intervenes elsewhere in each case.

On the negative branch with `A = +0.001`, the first polar crossing has `Pi = +2.18304`, `H_local = +0.24530`, and `H_V = +0.41816`. At `x = 0.564100` the field turns at `q = 0.403051`. The fitted potential gives `−46.8439` in the local Klein–Gordon acceleration, the spatial Laplacian contributes `−0.0421`, and the expansion term vanishes with the velocity. At the return crossing, `Pi = −1.97888`, `H_local = +0.17262`, and `H_V = +0.16212`. The corresponding turning time for the half-amplitude seed is `0.627588`, with `q = 0.401715`, potential acceleration `−46.0118`, and Laplacian `−0.0367`. The same-worldline minimum-to-minimum intervals are `0.154451` and `0.154821`. Both returns occur after the equator starts to contract and before the volume turnaround.

For negative seeds the crossing equator is still expanding, with `H_local = +0.32420` and `+0.37016` for the larger and smaller seeds. The polar values of `H_local` are already `−1.25387` and `−1.18984`. The equatorial field is still increasing at the global domain stop. Across sampled fixed worldlines, the maximum number of minimum passages is two for the positive seeds and one for the negative seeds. A third passage in the original direction would complete a crossing-to-crossing period. None occurs in the available evolution. A partial cycle at the poles gives no particle number or radiation density.

## Domain, numerical control, and replay

The calculation reuses the prior nonlinear constraint initial data and six-field Einstein–scalar evolution over the whole reflection-symmetric three-sphere. The fitted spline covers `−1.318340064 < q < +0.538544129` to the precision shown. The primary DOP853 integration stops when a monitored field extremum comes within `1e−4` of either endpoint. Prior to each RHS evaluation a polynomial interpolation bound checks the full angular domain and recursively bisects any uncertain intervals. A trial stage at risk of leaving the hard domain is discarded and retried with a shorter step. The angular grid includes poles and equator. Its global extreme candidates are polished in angle when storing field extrema. The upper endpoint limits all positive-branch runs, and the lower endpoint limits all negative-branch runs. Reducing the guard to `1e−6` changes the times by less than `1.25e−5` and preserves the event sequence. At both margins the limiting field is moving outward. These finite-margin times mark the limit of the calibrated calculation. Gaussian chart extension and physics beyond the fit require additional information.

Every accepted time step is checked for the Hamiltonian and physical radial momentum residuals on an angular monitor distinct from the evolution grid. The numerical stop is `max(|Hamiltonian|/(6h²), |momentum|/h²) = 1e−4`. Volume and volume expansion use independent higher-order Gauss quadrature weighted by the evolving spatial volume element. Simultaneous expansion and contraction occur within one smooth chart. The smallest spatial Jacobian relative to the round reference at any of the eight primary endpoints exceeds 0.45.

Primary degree 192 and check degree 256 reach the domain guard for every seed. Their physical and stop times agree within `1e−10`. The largest sampled primary momentum residual is `1.37e−5`. At degree 256 the largest is below `2.2e−8`. Degrees 64 and 112 instead hit the numerical constraint guard for the two positive seeds on the negative branch. For `A=+0.001`, those stops occur at `0.646501` and `0.702985`, respectively. For `A=+0.0005`, they occur at `0.710362` and `0.766783`. These cutoffs move with angular resolution. Each positive-seed polar return precedes even the degree-64 numerical stop. A separate degree-192 run halves the maximum time step from 0.003 to 0.0015 and tightens the relative tolerance by a factor of ten. Its event times differ by less than `1e−10`. At the exact saved times of the earlier cap, all eight primary trajectories reproduce the peak field, proper mean, volume expansion and proper volume within `2e−13` in the solver's units. Machine-dependent small differences in residuals are recorded alongside the events.

From the study root, run

```sh
OPENBLAS_NUM_THREADS=1 python rollout/rollout.py --gravity gravity --output rollout/ROLLOUT.json
```

The JSON contains full primary event records, compact angular, temporal and domain-guard checks, and the maximum observed differences used in replay validation. An optional short run accepts `--levels coarse --amplitudes .001 -.001`.
