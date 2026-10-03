# Internal critical review of an active-source saved-trajectory diagnostic

Prepared 2 October 2026 in the user's Pacific time zone. This is an internal,
AI-assisted design review. No physical array payload, source value, mode flow,
or response integral was evaluated in preparing it. Prior sources, formulas,
protocols and published result summaries were read as text.

The earlier source-free diagnostic demonstrated an integration-ledger error on
the fixed interval `[-2.5,-1.5]`. It did not validate active-source propagation
or change the three earlier metric experiments' scientific FAIL status. A new
active-source diagnostic is motivated by those known results and must say so.

## Minimal identifiable question

Use both existing sources, both saved resolutions and all three existing
cutoffs. A fixed active interval `[-4.5,-3.5]` has retained raw modes at both
ends and contains nonzero source jets. Reconstruct forced flow from the actual,
unprojected initial modes, independently evaluate direct density and pressure,
and independently integrate their full finite-band Ward integrand. Compare the
stored stress increment, unchanged Simpson ledger, continuous-flow ledger and
measured mode-flow defect with signed quantities.

This is a local, anchored propagation/operator/integration test. A passing
result cannot remove an error inherited before the initial anchor. The whole
support interval `[-5,-3]` has no raw mode record at its final endpoint; it does
not permit the same measured endpoint-flow attribution without new data.

## Scientific requirements before a public freeze

1. **Independent physical operators.** Density and pressure must come from the
   minimal-scalar operator and complete retained fixed-reference subtraction.
   Neither may be obtained by solving the Ward equation. An analytic primitive
   reduced to an endpoint expression is legitimate only after its derivative
   is independently proved equal to the separately expanded integrand, with
   unrestricted initial complex amplitudes and all local contacts retained.
2. **Finite-band source contact.** For scaled response quantities the integrated
   Ward expression is `F=L*(R-3*P)-3*h_prime*(R0K+P0K)`. The final term generally
   does not vanish at fixed comoving cutoff. Using a continuum anomaly or a
   physical-momentum moving cutoff changes the calculation's target.
3. **Explicit continuous approximation target.** If a surrogate is used,
   canonical forcing and every source/contact derivative must derive from the
   same surrogate, or its force mismatch must be separately reported. For
   `u'=w`, `w'=i*2*k*w-epsilon*f`, direct metric stresses using `g[h]` satisfy
   `R'-F=L*(f-g)/(2*k)` before momentum weighting. Exact analytic bump values at
   independently fixed quadrature nodes avoid defining a new physical source,
   but source-flow quadrature and dense-output errors still need controls.
4. **Separate approximations.** Arithmetic precision agreement, source/forcing
   interpolation or quadrature refinement, time-ledger refinement, momentum
   refinement, and comparison with saved trajectories answer different
   questions. Eighty versus one hundred digits cannot detect a shared source
   approximation or shared phase bias. If only final reductions use MP while
   source, phases and modes use long double, the output must state that scope.
5. **No hidden accuracy claim.** A signed projection and triangle envelope of
   measured endpoint mode differences concern those differences only. They
   do not bound inherited state error, total trajectory error, momentum
   quadrature error or physical model error. Exact-ratio conversion preserves
   the saved represented values; it does not recover unrounded data.
6. **All registered cases and intervals.** No source, resolution, cutoff,
   observation, subinterval or control may be selected after outcome inspection.
   Keep global-prefix Simpson subtraction in its original long-double
   arithmetic. A same-fine-history doubled-step Simpson control changes only
   sampled integration and must not silently reevolve modes or stresses.
7. **Signed attribution and failure semantics.** Report `D_S=DeltaR-S`,
   `D_cont=DeltaR-I`, `E_Q=I-S` and their signed closure. A failed science gate
   is a retained negative result, distinct from integrity/schema/resource
   failure. Classification cannot retroactively replace an old FAIL receipt.
8. **Independent canonical-integral control.** If the primary uses an endpoint
   phase primitive, a second route should numerically integrate independently
   expanded direct `F` on dense fixed nodes, with a fixed higher-order or
   subdivided control. Two copies of the same endpoint expression check
   implementation, not ledger quadrature. Compare each control at each
   declared precision, not only the cross-precision difference.
9. **Actual resource scope.** Freeze one wall/RSS budget for every source,
   resolution, cutoff, approximation control and precision of each route,
   including provenance/member decoding and serialization. A budget cannot be
   reset per case. Preserve partial progress and exact resource failures.
10. **Authentic serialization and provenance.** Reject binary64 loss of raw
   long-double complex components, nonfinite values, missing cases, duplicates,
   changed histories, changed source/subtraction inventories and changed
   contact signs. Decimal strings must roundtrip within a fixed arithmetic
   scale. Physical decoding must remain inaccessible without the verified
   public freeze and explicit registration byte hash.

## Limits of a successful result

A complete consistent result would locate an error in an existing numerical
pipeline during a preselected active interval. It would not establish a new
mathematical law, a successful full metric calibration, shifted-root stability,
the lapse/bulk/state response, coupled Einstein evolution, heating, or support
for a higher-dimensional origin of the hot Big Bang. All prior scientific FAIL
receipts remain unchanged regardless of the new classification.

## Textual sources inspected

- `research/HDBLAST_CHECKPOINT_20261002_SOURCE_FREE_LEDGER/PROTOCOL.md`
- `research/HDBLAST_CHECKPOINT_20261002_SOURCE_FREE_LEDGER/PRIMARY_IMPLEMENTATION.md`
- `research/HDBLAST_CHECKPOINT_20261002_SOURCE_FREE_LEDGER/reports/RESULTS.md`
- `research/HDBLAST_CHECKPOINT_20261002_SOURCE_FREE_LEDGER/theory/NEXT_REGISTERED_LEDGER_DIAGNOSTIC.md`
- `research/HDBLAST_CHECKPOINT_20261002_SOURCE_FREE_LEDGER/theory/verify_discrete_ward.py`
- `research/HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP/theory/primary/METRIC_ACTION_AND_CONTACT_DERIVATION.md`
- `research/HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP/independent/forced_metric.py`
- `research/HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP/code/metric_primary.py`

This early review is provisional. A GO/NO-GO assessment requires reading the
actual staged protocol, implementation, output schema, mutation tests and
resource plan before physical evaluation.
