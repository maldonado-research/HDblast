# Saved metric diagnostics, with science FAIL preserved

The V2 CLI reads existing primary JSON, independent JSON and mode archives. It imports the source-pinned validator's pure schema/provenance helpers and a separate raw auditor. It never imports or runs physical producers, evolves modes, samples a new physical grid or computes new response histories. It reconstructs weighted observation quantities and Simpson Ward ledgers from saved modes and direct histories.

The original strict prototype and failures remain in `positive-partial-001/`, `positive-partial-002/` and `development/`. The first attempt stopped on an overly broad public-registration selector; it never entered numerical reconstruction. The second verified the actual prospective public record, then the frozen raw auditor rejected coarse baseline histories at gap `1.283694347524729309e-12` against its unchanged allowance `1.0000555111512312578e-12`. The fine run stopped on C_1 metadata reconstruction at gap `0.0120210368144252441775`. These original failures are retained, with zero accepted raw comparison coverage.

The background has `C=a''/a=2a²`, `a' = a²`. Therefore `C_n=2*(n+1)!*a^(n+2)`. The frozen raw auditor used `(n+2)!*a^(n+2)`, which agrees only at n=0. V2 changes that metadata expression alone; the physical producer, modes, arrays, grids, normalization, subtraction kernel and scientific gates remain unchanged. `verify_c_metadata_v2.py` proves the all-n recurrence, eight direct derivatives and wrong-formula witnesses. All 13 checks pass in normal and optimized Python. The original frozen auditor and its hash remain preserved. This correction explains a validator metadata rejection; it does not cure any science or numerical-gate failure.

Use the pinned environment and actual immutable checkpoint pins:

```bash
/workspace/hdblast-cloud-setup/venv-frw/bin/python \
  /workspace/hdblast-research-work/metric-saved-diagnostics/diagnose_saved_metric_v2.py \
  --checkpoint CHECKPOINT_DIRECTORY \
  --primary EXTERNAL_RESULTS_DIRECTORY/primary/results.json \
  --modes EXTERNAL_RESULTS_DIRECTORY/independent/results.json \
  --registration-sha256 ACTUAL_REGISTRATION_SHA256 \
  --public-freeze-commit ACTUAL_EARLIER_PUBLIC_FREEZE_COMMIT \
  --public-freeze-evidence ACTUAL_SAVED_PUBLIC_VERIFICATION_RECEIPT.json \
  --failure-evidence ACTUAL_SAVED_NONZERO_EXIT_RECEIPT.json \
  --output FRESH_DIRECTORY_OUTSIDE_CHECKPOINT \
  --record-numeric-mismatches --plots
```

`--record-numeric-mismatches` permits finite arrays of the same shape to continue after their original raw or JSON comparator threshold fails. Every original gap, allowance and failure is recorded; no threshold is raised. Dtype, structural, nonfinite and other original guards remain fatal. Archive consistency is explicitly UNVERIFIED, and the independently reconstructed quantities may be used only as diagnostics. Exit zero means the diagnostic completed; it never means scientific PASS or cross-validation. Exit 2 indicates that at least one archive remained blocked by a fatal audit error. A fresh directory is required, and source/data hashes are checked before and after the diagnostic.

Complete four-run coverage contains 72 reconstructed observation/cutoff points and 72 direct-history Ward endpoints, 180 core cross comparisons, 144 baseline/current comparisons, 180 core refinements, 36 current refinements and 36 Ward-ledger refinements. The earlier two positive-source archives provide only half that coverage. Missing runs and failed raw audits are always reported. Coarse Ward endpoints and baseline refinement differences have no extra invented gate. The 48 primary pointwise Ward checks use saved derivatives and stresses from the same primary producer; they provide a consistency diagnostic, not an independent derivative measurement.

The optional figure plots gaps divided by the unchanged gates and is titled with underlying science FAIL and actual partial/full saved coverage. It cannot promote a resource failure, failed ledger, incomplete archive set or raw mismatch to a calibration PASS. It does not claim a certified continuum limit, stability, heating, general lapse/actual-root/bulk/state closure or a cosmological discovery.
