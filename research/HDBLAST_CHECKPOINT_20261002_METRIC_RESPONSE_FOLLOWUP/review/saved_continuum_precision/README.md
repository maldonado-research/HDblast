# Optional saved-continuum precision audit

This post-run audit is outside the prospectively frozen 32-command experiment.
It evaluates saved JSON decimals, floats, arithmetic identities and hashes. It
never imports or calls a physical source evaluator, mode solver, stress
producer, or quadrature function. It changes no frozen source or gate.

The final source is `audit_saved_continuum_precision.py`, SHA256
`a2027962cee12cde7b1be2164b1efd17178b6a5d47a32cbc4430675213843b40`.
The actual primary results audited have SHA256
`3ebb2fe671404b3e314b65541d9d7313a014dd58ab565d6348efefa285620dc2`.
Final captures are `normal-002/` and `optimized-001/`. Both pass and produce
identical counts and row diagnostics. The initial passing audit source is
preserved under `development/initial-saved-data-audit/`, with its actual
capture under `normal-001/`. No audit failures occurred; any future failure
is captured in a fresh output directory and never overwritten.

The audit verifies all 12 original rows and both precision runs per row:

- 72 qjet/contact reconstructions from saved memories and endpoint jets.
- 72 raw quadrature-error normalization checks.
- 36 reported precision gaps against saved-decimal rounding ranges.
- 36 final qjet float conversions and 36 raw-memory float conversions.
- 48 saved double forcing recurrences from already-saved h jets, and
  descriptive comparisons of all four endpoint jets across both precisions.
- 36 effective-allowance float bindings and upward normalization checks.
- 120 propagated error fields, 36 metric variance sums, 24 stress sums,
  12 fixed-phi current equalities and 12 direct trace arithmetic checks.
- Original result/provenance hashes and all 331 prospectively frozen files,
  unchanged before and after the audit.
- Five rejected mutations of saved metadata, without any physical evaluation.

The largest declared qjet precision gap is 8.8327068e-51. The largest saved
qjet float-conversion diagnostic is 1.1798583e-16, and the largest normalized
serialized allowance is the same at reported precision. These are saved
numerical diagnostics, not certified physical errors.

An important limit is that printing a low-precision qjet with 50 significant
decimal digits does not retain every bit of its internal 50-dps binary
mantissa. The audit therefore checks that each reported 50/70 gap lies in the
rounding range compatible with both printed qjets; it cannot reconstruct the
original gap exactly to 70 digits. Likewise, underlying B/uB endpoint values,
actual quadrature nodes, degree reached, and the mathematical accuracy of the
quadrature estimates are source-bound declarations rather than new independent
physical calculations. Decimal/context-rounding consistency tolerances in this
audit concern saved arithmetic only. They are not additions to the frozen
physical acceptance gates or interval enclosures of the original integrals.

Run with the pinned environment, using a fresh output outside the checkpoint:

```bash
source /workspace/hdblast-cloud-setup/activate-frw.sh
python -B audit_saved_continuum_precision.py \
  --checkpoint CHECKPOINT_DIRECTORY \
  --results PRIMARY_RUN_DIRECTORY/results.json \
  --output FRESH_EXTERNAL_AUDIT_DIRECTORY
```

The saved results directory must also contain its original `started.json`
and `EXECUTION.json`. Use `python -B -O` for the optimized audit. All checks
remain active under optimization. `SUMMARY.json` records the final counts and
scope, and each capture includes its actual command and source hash.
