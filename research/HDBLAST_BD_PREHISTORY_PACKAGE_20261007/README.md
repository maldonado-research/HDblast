# Portable additive BD prehistory replay

Place this directory at `research/HDBLAST_BD_PREHISTORY_PACKAGE_20261007`, beside
the immutable `HDBLAST_CHECKPOINT_20261007_BD_PREHISTORY` directory. It introduces
no Python source inside that frozen checkpoint. A self-contained archive must
retain both sibling directories.

The custodian supplies `PUBLIC_GO.json`, `REPLAY_PINS.json` and
`PAYLOAD_MANIFEST.json`. Expected SHA256 values for the pin file and full payload
manifest must come from the independently reviewed release readback and be
provided explicitly on the CLI. Local self-produced pins do not establish
public authentication. The checkpoint's frozen guard independently checks the
full registration, root public-GO receipt, immutable commit/blob identities and
complete review before any numerical route import.

The payload manifest has exactly `schema_version:1` and `files`, mapping every
checkpoint-relative leaf to `{bytes,sha256}`. It includes the full registration,
registered preparation, sealed original numerical results and continuous-model
certificate. It excludes this sibling directory, the manifest itself and new
replay outputs. Unlisted leaves, changed bytes and symlinks fail.

The pin file has exactly these keys:

```
schema_version: 1
payload_manifest_sha256: <complete checkpoint manifest SHA256>
registration_sha256: <frozen FULL_REGISTRATION.json SHA256>
public_go_sha256: <authenticated root PUBLIC_GO.json SHA256>
reference_data:
  primary: evidence/actual/primary-normal/DATA.json
  independent: evidence/actual/independent-normal/DATA.json
reference_continuous_certificate: evidence/actual/CONTINUOUS_MODEL_CERTIFICATE.json
```

Use Linux with `libseccomp.so.2`, an existing nonroot user, exact Python3.12.14,
and the checkpoint's `requirements-replay.txt` (`python-flint0.9.0`,
`sympy1.14.0`, `mpmath1.3.0`). Hosted Ubuntu CI is supported. Install requirements
into an explicit virtual environment outside the repository; preserve the pins
and TLS/artifact verification. The helper itself uses standard-library code.

To authenticate saved evidence without source evaluation:

```
python -B research/HDBLAST_BD_PREHISTORY_PACKAGE_20261007/replay_checkpoint.py \
  --pins-sha256 a839dc8ccd4efb2be41711093987eabd2e7ee39bdf4fc2fd9be646db82c72243 \
  --payload-manifest-sha256 1bb65ecadc6ea11b9ffc58c4e88021b82728bb3696f50a7941468b4428a1b83a \
  --verify-only
```

Verification checks all payload bytes, full registration and public GO, both
original bounded entry receipts and source journals, all36 output rectangles
and72 component intersections, and the exact continuous-model mathematical
core with its original budget/entry provenance. It evaluates no new source.

To repeat the authorized registered target in four fresh bounded processes:

```
python -B research/HDBLAST_BD_PREHISTORY_PACKAGE_20261007/replay_checkpoint.py \
  --pins-sha256 a839dc8ccd4efb2be41711093987eabd2e7ee39bdf4fc2fd9be646db82c72243 \
  --payload-manifest-sha256 1bb65ecadc6ea11b9ffc58c4e88021b82728bb3696f50a7941468b4428a1b83a \
  --python /absolute/path/to/venv/bin/python \
  --output-directory /absolute/fresh/external/replay-directory
```

This calls the original frozen custodian for primary/independent routes in
normal and optimized Python, with the original900-second/512-MiB/20-MiB caps.
Each worker authenticates before source callbacks and denies network/process
creation. Every route's DATA must match both the sealed original and its other
Python mode exactly. All scientific BUDGET fields must also match the sealed
original, including source models, serialized traces and error budgets; only
the registered measured `wall_seconds` and `peak_rss_kib` fields may vary.
Each mode must pass all36 independent rectangle and72
component intersections. Frozen theorem, fabricated-target, outer-guard and
continuous-model self-tests must match retained preparation receipts. The
independent cap verifier matches the added
`evidence/actual/INDEPENDENT_CAP_CURRENT_{NORMAL,OPTIMIZED}.json` receipts
regenerated from the unchanged frozen verifier. Its old preparation receipts
remain preserved: they recorded31 checks before the verifier gained the32nd
check, `twice integrated nested drift source identity`. The first replay stopped on that receipt
mismatch after both normal-mode numerical routes passed; that failed run is
retained separately. No numerical source, gate or theorem is changed.
The continuous coefficient/source certificate is reproduced from each fresh
independent budget/entry; its complete mathematical core must match, while
execution-specific provenance is independently authenticated and retained.

The helper writes a final `REPLAY_RECEIPT.json` and preserves failed step logs
and underlying exit codes. Outputs remain outside both immutable sibling
directories. A full replay evaluates88 authorized real source jets; no saved
quantum arrays are decoded. It certifies the prescribed nine-probe target and
the separately scoped continuous source/coefficient representation. Original
binary80 state error stays NOT_ENCLOSED, full twelve-case certificate UNRESOLVED,
metric calibration FAIL and higher-dimensional Big Bang cause NOT_ESTABLISHED.
