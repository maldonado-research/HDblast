# Separate exact Cauchy packet response

`W1_exact_Cauchy_packet` keeps the exact-packet displacement and replaces the
smoothed P1 velocity with the complete nominal analytic velocity. This is a
separate function and candidate. The earlier exact-packet function and its
same-function bound refinement remain intact.

The input adapter explicitly takes the refined absolute displacement caps
from the absolute-refinement report, the unchanged displacement-error caps
from PR109, and both fresh absolute velocity caps and fresh paired velocity
errors from the new native velocity proof. The new absolute velocity caps
are larger than the old smoothed P1 allowances; they are nevertheless used in
full. No old velocity allowance is inherited for this candidate.

Both native and independent whole-domain initial and fresh local velocity
gates must pass first. The adapter checks the full 23-component norm, exact
completed velocity value/r affine outer traces, the independent C1 velocity
regularity audit, parity isometry, flat origin collar, and separated nominal
denominators. Displacement is globally C2; velocity is globally C1. The
physical-side C2 response retains the published source/coefficient regularity,
compatible corner and one-sided shell-trace premises.

The response uses exactly the same original target `U`, `F`, unfiltered `F_x`,
`G`, `G_x`, `M`, operators, domain and full 23-component `H` norm. All 35 PR109
ledger nodes and constants are reproduced before substitution. Both initial
acceleration charges `f00`/`r00`, full Robin forcing, derivative commutator,
feedback inverses and `r1(0)=0` cancellation are retained. Separate absolute
`W1` prescribed-forcing caps are computed as well.

| Conditional upper certificate | PR109 | New Cauchy packet | Old/new |
| --- | ---: | ---: | ---: |
| C1 X | 1.60214883e23 | 1.30333362e23 | 1.22927 |
| xx and xr H | 2.13407313e25 | 3.62261174e24 | 5.89098 |
| rr H | 2.14322475e25 | 3.69988619e24 | 5.79268 |

Exact rationals and the required gates, every second-error upper certificate
at most half PR109 and C1 no larger, are in `reference/RESPONSE.json`.
These are global coordinate upper-certificate comparisons. The separate
fresh initial velocity witness supplies a local true-error result. No global
true-error factor, evolved point/quadrature table, natural or constraint
residual, invariant observable, nonlinear shell-contact solution, novelty or
cosmological result follows from this calculation.

The arithmetic helper is the already frozen local Fraction-only recipe; no
numerical producer is imported. All selected original and proof sources are
authenticated before helper import and read back afterward. Normal and
optimized Python canonical reports must be byte identical. The independent
ledger audit reconstructs all numerical nodes separately. All 60 fresh subprocess controls are rejected in normal and optimized Python.
The controls include stale old velocity charges, force/Robin/commutator/initial
charge omissions, table/source changes, and false report caps, constants,
roles or scope flags. Complete canonical report verification rejects every
advertised-field mutation.

```sh
python3 checks/run_replay.py \
  --repo-root /path/to/published/repository \
  --same-response-dir /path/to/same-function-response \
  --absolute-dir /path/to/absolute-refinement \
  --velocity-dir /path/to/new-native-velocity \
  --independent-velocity-dir /path/to/new-independent-velocity \
  --output-dir /path/to/fresh/replay
```

The exact minimal input inventory is `checks/FREEZE.json`. All roots are
explicit CLI inputs; canonical reports contain role-relative source hashes.
