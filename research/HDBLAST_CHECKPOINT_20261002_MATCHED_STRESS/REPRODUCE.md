# Reproduce the registered stress calibration

Use a Linux Python 3.12 runtime whose NumPy long double has at least 63
fraction bits. Install the exact wheels in `requirements-replay.txt`:

```bash
python -m pip install --only-binary=:all: -r requirements-replay.txt
python code/replay_stress.py --output /an/unused/external/directory
```

The driver sets numerical thread counts to one, verifies the complete
payload manifest and prospective inputs before and after execution, and
uses fresh output paths. It requires no private folder or full repository.
Two included immutable support files are used solely to reconstruct the
original source-pin paths for the frozen symbolic proof in a temporary
directory. The driver changes no frozen source or physical setting.

The 22-command replay runs primary algebra, the source-polynomial bridge,
stress-tail algebra, complete pulse-root proofs, finite-cutoff Ward algebra,
and validator guard tests in normal and optimized Python. It then runs the
primary physical calculation, all four independent mode resolutions, and
both actual validators. Post-run saved-mode audits and the separate
analytic second-order proposal are also checked normally and optimized,
followed by the regenerated summary and nine figures. Real command logs, return codes, timing and failures
are retained. Each physical producer has its registered 900-second budget;
the hosted workflow allows 45 minutes for the full replay.

The independent archive preserves extended-precision arrays and direct
stress histories. Long-double storage padding may vary between processes;
scientific array values and each run's actual archive hash are authoritative.
Byte equality between distinct fresh NPZ files is not a reproducibility gate.
Acceptance gates test the frozen experiment, not merely self-consistency
within one producer. No numerical grid or tolerance is adjusted after seeing
results. Original evidence remains under `outputs/`.

`code/build_package.py` checks all frozen hashes, generates the outer manifest,
constructs the deterministic ZIP twice, and verifies every member against
the source. The outer manifest excludes only itself, the ZIP and its hash
sidecar. Published local receipts and exact-head hosted CI establish what
was actually executed; merely preparing this driver does not establish a
successful replay. The omitted-band interval bounds do not certify the
complete numerical solution.

The [first complete local replay](evidence/replay_local_18/VALIDATION.json)
passed its original 18 commands. Its 444 fresh mode-array fields, containing
4,938,184 elements, equal the original numerical values exactly; NPZ storage
padding is excluded from that comparison. The final extracted-ZIP replay
receipt is [published beside the checkpoint](../HDBLAST_MATCHED_STRESS_FINAL_ZIP_REPLAY_20261002.json)
to avoid embedding a receipt containing the archive hash inside that archive.
The hosted job replays the complete final 22-command program.
