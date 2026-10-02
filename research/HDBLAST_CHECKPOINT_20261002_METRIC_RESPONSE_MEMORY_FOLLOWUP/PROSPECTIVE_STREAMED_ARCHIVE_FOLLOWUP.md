# Prospective storage-only memory followup

Proposed checkpoint:
`HDBLAST_CHECKPOINT_20261002_METRIC_RESPONSE_MEMORY_FOLLOWUP`.
This document explains a bounded numerical-resource repair. It contains no
new physical evaluation and is not itself a public registration.

The high-precision followup completed the primary and the independent
positive_B coarse/fine runs, then stopped at peak RSS=306088 KiB above its
262144-KiB gate. The next separate registration can change archive retention,
serialization and file hashing while preserving the entire numerical
calculation, physical model, sources, precision schedule, gates and budgets.
The previous two failed registrations remain immutable.

## Restrict the change to storage

Retain the same H=1, x=r=2, xi=0 background, fixed BD state, fixed physical
mass/scalar/reference, two metric pulses, amplitude, source jets, canonical
forcing, observation times, cutoffs, quadratures, 50/70-digit primary
continuum arithmetic, independent longdouble/clongdouble arithmetic,
coarse/fine step and momentum settings, minimal stress/subtraction formulas,
Ward ledger, controls, normalizations and all acceptance gates.

Retain the **900-second producer limit and 262144-KiB peak-memory limit**.
Do not fit the budget by dropping source profiles, observations, modes,
history points, raw integrands, subtraction terms or metadata, or by reducing
precision. The independently scheduled four evolutions and both normal and
optimized validators remain required.

The storage repair should stream SHA-256 in bounded chunks rather than
reading an entire NPZ file into memory. If archive members are streamed or
written incrementally, retain every existing array, member name, shape,
dtype, complex phase, longdouble precision and historical node. Streaming
storage must not replace raw arrays with summaries or binary64 copies.
Preserve all per-run diagnostics and failure captures.

Source comparison should demonstrate that numerical source, propagator,
observable, cutoff integration, baseline, ledger and refinement function
bodies remain unchanged. Any source edits needed for archive lifetime,
incremental serialization, bounded hash reads, output finalization or
provenance must be listed explicitly as storage changes. Updating the
independent manifest and registration is required because its executable
source changes, even when the numerical bodies do not.

## Verify storage without a new physical solve

Before freezing, pure or synthetic checks can establish streamed hash
agreement, full archive-member coverage, exact shape/dtype preservation,
round-trip equality for real/complex arrays, finalization order and failure
behavior. Existing saved arrays may be inspected or reserialized as
historical storage fixtures without re-evolving their modes; such checks
must be labeled storage tests, not a new response experiment.

The archive should be finalized and its streaming hash recorded before its
run is marked complete. Failure should retain the completed members and
available active snapshots with honest completion status. Resource checks
must continue to measure peak process RSS, not only current live memory;
freeing an allocation after an excessive peak cannot turn that run into a
PASS. No allocation guard may be moved or weakened to conceal a transient
peak.

The actual numerical replay remains the necessary resource test. Lower
memory is a prediction of the storage repair until the full registered
producer finishes inside the unchanged budget. A static explanation of the
earlier allocation cannot establish that result.

## Freeze, execute and classify separately

Publish the distinct memory-followup source, complete configuration,
unchanged gates/budgets, independent manifest, retained numerical-body
comparison and exact storage-test evidence before its first physical call,
then verify the remote frozen tree. Preserve and reference the original
roundoff failure and the later resource failure with their exact freezes,
source pins, partial outputs and raw archives.

Execute the complete scheduled replay with fresh output paths. Report
primary completion, all four independent runs, actual peak RSS and time,
cross-route validation, refinement, raw reconstruction, Ward/trace checks,
continuum-tail comparisons and controls according to their actual statuses.
Any warning, memory/time overrun or failed gate stops this new round with
all attempted output preserved. No unregistered setting changes follow a
failure.

A successful memory followup could establish the bounded homogeneous
special-point metric calibration on prescribed geometry. It would leave
actual shifted-root numerical response, the full retarded metric/state
matrix, gravitational/scalar bulk constraints, moving-shell matching,
physical EFT validity, coupled stability and heating as separate
prerequisites. The existing UV bounds remain conservative; finite-K accuracy
would not automatically certify a precise continuum pressure or its sign.
