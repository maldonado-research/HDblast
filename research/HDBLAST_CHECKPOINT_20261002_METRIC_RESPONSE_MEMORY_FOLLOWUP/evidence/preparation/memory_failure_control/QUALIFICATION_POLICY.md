# Memory followup expected scientific failure control

The control invokes the unchanged frozen 32-command replay once. It succeeds only when the registered experiment fails its internal scientific gates after completing all 12 primary rows and all four independent coarse/fine runs. A control PASS is recorded as `EXPECTED_SCIENTIFIC_FAILURE_CONFIRMED`; the underlying scientific result remains `FAIL`, and `registered_metric_calibration_passed` remains false.

Invoke the helper **from its registered checkpoint copy**, normally under the pinned Python 3.12 runtime:

```text
python CHECKPOINT/code/reproduce_metric_expected_scientific_failure.py \
  --checkpoint CHECKPOINT \
  --output FRESH_EXTERNAL_DIRECTORY \
  --public-freeze-commit PUBLIC_COMMIT_40_HEX \
  --registration-sha256 FULL_REGISTRATION_SHA256 \
  --independent-manifest-sha256 INDEPENDENT_MANIFEST_SHA256
```

The supplied bindings must describe the independently published prospective freeze. The helper requires them to match `FREEZE_RECEIPT.json`, authenticates every full-registration input before importing the replay's pure verifier, and requires its own running bytes to match `code/reproduce_metric_expected_scientific_failure.py` in the registration. Sources are reauthenticated after execution. Public Git publication is established by the calling workflow's provenance; this helper performs no network publication check.

The exact replay plan and executed command arguments are checked. Expected execution is 28 attempted commands, 27 successful commands, and an ordinary nonzero exit from `independent_modes`; the remaining validation/summary/figure commands are not reached. Failure before that stage, timeout, signal, source change, unsupported producer exception, unexpected success, or an output hash discrepancy fails the control.

All four independent archives must match their recorded SHA256 and the complete registered `independent/PRESERVED_ARCHIVE_SCHEMA.json`: exactly 631 unique members each, correct coarse/fine shapes, endian dtypes, C order, payload lengths, and ZIP CRCs. This is saved-byte inspection and performs no physical evaluation. A complete archive alone does not establish completion because the producer saves an archive before its final resource check.

Actual primary budget metadata uses `primary.producer_wall_seconds=900`. The independent final elapsed time must be at most 900 seconds and peak RSS at most 262144 KiB. All scientific failures are reconstructed from saved run JSON and the unchanged frozen gates. The known required witness is the positive_B fine Ward endpoint at eta=-1.5 and K=256 exceeding 2e-6. Declared maxima and combined/fine values must agree with saved runs. Unexpected failure kinds are fatal.

The replay runs in a new process group. Outer timeout or interruption kills and reaps the group before output hashing, so descendant producers cannot continue mutating failed output.

Preparation checks use explicitly fabricated metadata and zero-valued streamed NPZ fixtures, with the configuration and schema read from the actual third checkpoint. Final normal and optimized receipts each report 48 passing checks, zero physical evaluations, zero replay subprocesses, and `production_control_executed=false`. The saved prior memory-budget failure is separately rejected as an unsupported scientific-failure control outcome.

Final receipts:

- `proof/FINAL_ACTUAL_NORMAL/CHECKS.json`
- `proof/FINAL_ACTUAL_OPTIMIZED/CHECKS.json`

The actual post-freeze control has not been executed during this preparation. Cross-route calibration, conservation acceptance, precise continuum stress, coupled evolution, heating, stability, and actual-root/bulk/state results are not established by a successful expected-failure control.
