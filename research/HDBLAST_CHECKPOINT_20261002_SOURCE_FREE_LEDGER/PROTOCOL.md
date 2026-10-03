# Fixed prospective diagnostic protocol

The parent public science commit is
`71d00cc423e9049c8166b7ee7afbefbac28d8a18`. Its three metric experiments remain
scientific **FAIL**, including 29 endpoint and 30 ledger-refinement failures in
the complete memory-only run. The current test was motivated by those known
failures. Its new numerical quantities remain uncomputed before registration.

## Cases and immutable data

Use `positive_B` and `signed_uB`, each at saved `coarse` and `fine` resolution,
and `K=64,128,256`. The only diagnostic interval is `[-2.5,-1.5]`. Compare every
saved density and pressure sample in this inclusive interval: 129 coarse and
257 fine samples. Earlier samples between the end of support and the anchor
are not part of this profile comparison. No interval, cutoff, state, tolerance
or algorithm is selected after seeing the diagnostic result.

The capsules copy 17 required NPY members verbatim from each original archive.
Their manifest binds the full original archive size, SHA256 and Git blob at
the parent commit, each selected member's bytes and headers, and each capsule.
Extraction and checksum checks do not interpret physical array values. Complete
stored histories remain available for the unchanged Simpson prefix arithmetic.
All four original archives and their original results remain preserved.

Input long doubles are converted through exact integer ratios. Complex inputs
are split into their long-double real and imaginary components; conversion
through Python `float` or `complex` is forbidden. Numeric roundtrips check the
actual value, including separately recorded zero signs where relevant, rather
than unused long-double padding bytes.

Match the producer's numerical normalization constants exactly:

| Constant | Exact binary-rational value |
| --- | --- |
| `np.longdouble('0.0001')` | `3777893186295716171 / 37778931862957161709568` |
| `np.arccos(np.longdouble(-1))` | `14488038916154245685 / 4611686018427387904` |

Momentum nodes, weights and time values retain their exact stored ratios. The
analytic background uses `L=-1/eta` evaluated from those exact times at the
declared precision. The geometric rounding in the old producer is consequently
included in the comparison with its stored history; it is not projected away.

## Exact free flow and ledger

For each momentum, start from the actual saved `u_4,w_4` at `eta_a=-2.5`:

```
Omega = 2*k
d = w_a/(i*Omega)
c = u_a-d
E(t) = exp(i*Omega*(t-eta_a))
u(t) = c+d*E(t)
w(t) = w_a*E(t)
```

Retain `Re(c)` without normalization or Wronskian projection. The exact
antiderivative of the independently expanded free Ward integrand is

```
I_k = 3*Re(c)*(L_b^2-L_a^2)/(2*k*epsilon)
      +Re(d*(3*[L^2*E]_a^b-i*Omega*[L*E]_a^b))/(2*k*epsilon).
```

The canonical `I_ab` evaluates the endpoint phase directly with MP exponential
at each declared precision. A recurrence-derived primitive is reported
separately and must agree in absolute value within `1e-12` at **each** precision;
agreement between 80 and 100 digits alone cannot hide a shared primitive bias.

Integrate with the original momentum weights and numerical measure
`k^2/(2*pi^2)`. Neither measured final density nor a conservation-defined stress
may be used to construct this ledger. The preserved proposal and symbolic proof
give the unrestricted density, pressure and signed mode-flow defect formulas.

## Stored Simpson arithmetic and signed comparisons

Reconstruct the original `simpson_to` expression in NumPy long double on the
full stored ledger integrand. Compute global prefixes `S_b` and `S_a`, subtract
them **in long double**, then convert that result exactly for subsequent
high-precision arithmetic. This quantity is `S_ab`. Also report direct reset
Simpson and its difference from the prefix subtraction; it does not replace
the frozen definition.

For the fine same-trajectory control, use the full global history `F_fine[::2]`
and double its time spacing, preserving the corresponding even prefix indices.
Again subtract the two prefixes in long double. Only Simpson sampling/weights
change; no modes or stresses are reevaluated to create this coarsened history.

From the exact stored endpoint values and independently integrated ledger form

```
DeltaR = R_b-R_a
D_S = DeltaR-S_ab
D_cont = DeltaR-I_ab
E_Q = I_ab-S_ab
D_S = D_cont+E_Q.
```

Keep signed quantities. Separately report the weighted signed final-density
projection of the measured endpoint mode defects and its triangle bound. The
bound applies to those measured defects, not to the inherited initial-state
error or total physical error.

## Fixed criteria and interpretation

All twelve cases must satisfy the profile agreement scale `2e-7`, continuous
ledger residual `abs(D_cont)<=2e-7`, and signed mode-flow density projection
`abs(E_flow)<=2e-7`. The arithmetic precision-gap, closure and serialization
threshold is `1e-12`. Implementations use fixed decimal precision levels 80 and
100; their exact arithmetic algorithms, rounding and reported quantity sets
are specified in the separately frozen implementation notes.
Every common scalar and profile in `OUTPUT_SCHEMA.json`, including the fine-only
fields, must agree across the two implementations within an absolute `1e-12`.
Within each implementation every reported non-metadata scalar/profile is in
the precision-gap universe. The independent per-node streamed endpoint defects
also participate; their weighted prefixes must reproduce the reported ledgers
and flow projections. Normal and optimized validators must give identical
substantive results. These are internal arithmetic checks, not physical error
certificates or independently derived physical models.

Only call a post-support ledger error demonstrated at the original gate scale
when every consistency criterion passes and at least one case has
`abs(D_S)>2e-6`, `abs(D_cont)<=0.1*abs(D_S)` and
`abs(E_Q)>=0.9*abs(D_S)`. Otherwise identify each failed criterion or report that
no gate-scale attribution was demonstrated. A positive result cannot change
the earlier calibration status. No factor-15 Richardson estimate is treated as
a certified error bound, and high-precision agreement is an empirical arithmetic
check rather than a certificate of trajectory accuracy.

Resource accounting and the exact reproduction interface are frozen with the
implementations before execution. Each separately executed route has one
900-second and 262144-KiB budget for **all twelve cases and both precisions**,
starting before provenance checks/member decoding and ending after output
serialization. Outer process receipts are authoritative. The independent route
is a separately labeled complete execution with the same budget; no combined
runtime claim follows from two individual passes. The physical driver runs
normal Python 3.12.14 with the pinned dependency versions and x87 binary80 ABI.

A completed diagnostic returns exactly `LEDGER_ERROR_DEMONSTRATED`,
`NO_GATE_SCALE_ATTRIBUTION`, or `CONSISTENCY_FAILURE`. A scientific consistency
failure is a retained negative result with successful *execution*. Integrity,
resource, arithmetic, schema, serialization, or cross-route failures make the
execution fail and prevent a verified classification. Every outcome preserves
`old_metric_status="FAIL"`. The public freeze precedes all physical diagnostic
evaluation. Results will appear in `reports/RESULTS.md`; no outcome is inserted
into these frozen prospective files afterward.
