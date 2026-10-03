# Exact source-free ledger reproduction

This tooling authenticates registered bytes, runs two separate complete
twelve-case diagnostics at fixed 80 and 100 decimal digits, and consumes the
authoritative scientific classifier. It never changes the earlier metric FAILs.
Each entire route has one 900 second / 262144 KiB allowance, including both
precisions. A negative scientific classification is a completed calculation;
nonzero exits, missing outputs, exceeded resources or arithmetic/integrity
failures remain execution failures with their actual receipts preserved.

At preparation time all new diagnostic quantities are UNCOMPUTED. Do not run
the physical commands until the frozen sources and FULL_REGISTRATION.json have
been published, their remote commit verified, and FREEZE_RECEIPT.json recorded.
Checksumming inputs and the synthetic tests perform zero physical evaluations.

Use Python 3.12.14 and install the five exact versions in requirements-replay.txt.
The driver checks installed metadata before launching any route, runs physics
under normal Python, fixes all numerical thread settings to one, and runs the
postrun validator in normal and optimized Python. Its output directory must be
new and outside the checkpoint. The expected inputs are:

```json
{
  "reproduction": {
    "routes": {
      "primary": "code/diagnostic_primary.py",
      "independent": "independent/diagnostic_independent.py"
    },
    "validator": "code/validate_ledger.py",
    "result_filename": "diagnostic.json",
    "validation_filename": "VALIDATION.json",
    "postfreeze_output_directories": ["outputs", "reports", "figures", "evidence"]
  }
}
```

FULL_REGISTRATION.json has schema_version 1 and files mapping every frozen
relative file path to its SHA256. The registration itself and subsequent
FREEZE_RECEIPT.json cannot appear in that mapping. The post-freeze receipt
contains public_freeze_commit (full 40-character lower-case commit) and
registration_sha256 (SHA256 of the complete FULL_REGISTRATION.json bytes).
Registration may additionally bind experiment_sha256, input_manifest_sha256
and frozen_configuration; the integrity helper checks them when present.
Unregistered payload files are forbidden outside the explicitly frozen output
directories. Extras in those directories must have known data/report suffixes;
unregistered modules, archives and executable files are rejected.

Both physical routes receive:

```text
--checkpoint-root CPP --registration-sha256 HEX --freeze-commit FULL40 --output-dir FRESH
```

They write diagnostic.json, progress.jsonl and any route-specific detailed
outputs under their fresh directory. The validator receives the same pin
arguments plus --primary FILE, --independent FILE, --execution-receipts FILE
and --output-dir FRESH. Its VALIDATION.json contains classification in
LEDGER_ERROR_DEMONSTRATED, NO_GATE_SCALE_ATTRIBUTION or CONSISTENCY_FAILURE,
and old_metric_status equal to FAIL. Normal and optimized reports must agree
after removing only the optional python_optimization field. Successful replay
status COMPLETED_REPLAY reports execution completion, with the scientific
classification retained separately.

The public Git verifier is an extra read-only check when the complete Git
objects are available. It verifies the exact registered blob bytes at the
recorded freeze commit; its PASS_LOCAL_FREEZE_GIT_BYTES status describes local
Git bytes, not evidence of public publication or its timing. The root separately
verifies the remote public freeze before issuing permission to evaluate. The
portable replay also works from the standalone ZIP
without Git, using the explicit external commit/hash and protected receipt.

```bash
python CPP/code/test_replay_package_guards.py
python -O CPP/code/test_replay_package_guards.py
python CPP/code/verify_public_freeze.py --checkpoint-root CPP --repository-root REPO --registration-sha256 HEX --freeze-commit FULL40
python CPP/code/replay_ledger.py --checkpoint-root CPP --registration-sha256 HEX --freeze-commit FULL40 --output-dir FRESH
```

--plan-only authenticates and records the complete command plan while launching
zero children. Replay preserves stdout, stderr, signed child exit codes,
timeouts, peak resident memory, elapsed wall time and current fresh outputs.
Linux child-subreaper isolation rejects surviving descendants even if they
escape the original process group, kills them and preserves their exit/resource
receipts. The registered routes use one Python process and one numeric thread
each. The child receives RLIMIT_NPROC=0, prohibiting creation of descendants or
extra threads in the kernel; execution requires a nonroot Linux UID. Resident
memory accounting combines authoritative wait4 peak RSS and live process-group
samples. Resource enforcement is for these fixed single-process routes.
It runs the other independent route after a producer exit failure and skips
scientific validation if either route failed. A frozen-source change prevents
further execution. Source bytes and the entire registration/receipt are checked
before, between and after commands.

After recording the original computed outputs and reports inside CPP, build
the complete deterministic standalone package. MANIFEST.json hashes every
payload file including FREEZE_RECEIPT.json and recorded outputs. Only
MANIFEST.json itself, the fixed checkpoint .zip and .zip.sha256 bookkeeping
names, and Python bytecode caches are excluded. The input capsules stay inside
the payload. The ZIP uses a fixed timestamp, regular-file permissions, sorted
members and fixed compression settings; two constructions must match.

```bash
python CPP/code/build_package.py --checkpoint-root CPP --registration-sha256 HEX --freeze-commit FULL40 --output-dir FRESH_PACKAGE --extract-to FRESH_EXTRACTION
python FRESH_EXTRACTION/HDBLAST_CHECKPOINT_20261002_SOURCE_FREE_LEDGER/code/replay_ledger.py --checkpoint-root FRESH_EXTRACTION/HDBLAST_CHECKPOINT_20261002_SOURCE_FREE_LEDGER --registration-sha256 HEX --freeze-commit FULL40 --output-dir FRESH_ZIP_REPLAY
```

The final fresh ZIP replay receipt stays outside CPP and the ZIP payload. The
packaging receipt reports packaging integrity only. Retain a failed replay's
progress and actual error receipts without replacing them with inferred
success. An existing payload manifest or output directory is never overwritten.

The staged hdblast_source_free_ledger.yml must be installed into
.github/workflows only after the public freeze is verified.
Its literal __PUBLIC_FREEZE_COMMIT__ and __FULL_REGISTRATION_SHA256__ placeholders
must first be replaced with the verified original freeze pins; CI compares the
receipt against those reviewed values and requires the freeze to be an ancestor
of checked-out HEAD and the fetched public
origin/research/20261002-source-free-ledger branch. Main is not used as the
publication anchor. It runs on pushes
and pull requests touching this checkpoint, installs exact dependencies,
checks public freeze bytes, performs a complete new diagnostic and retains
all outputs even when replay fails. It has no scheduled or continuous service.
