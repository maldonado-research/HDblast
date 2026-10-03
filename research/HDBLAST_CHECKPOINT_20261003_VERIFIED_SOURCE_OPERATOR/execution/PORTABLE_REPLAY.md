# Portable completed-checkpoint replay

Copy `replay_checkpoint.py` to `execution/replay_checkpoint.py` before the
prospective source freeze. The eventual completed archive must contain the
unchanged helper, `FULL_REGISTRATION.json`, `FREEZE_RECEIPT.json`, all frozen
sources/proofs, and original actual executions under `results/primary` and
`results/independent`. The final `MANIFEST.json` covers every payload file except
itself, using exact schema1 `{files: {relative: {sha256,bytes}}}`.

Invocation requires explicit externally verified pins and interpreter:

```sh
python -B CHECKPOINT/execution/replay_checkpoint.py \
  --checkpoint CHECKPOINT \
  --expected-manifest-sha256 MANIFEST_SHA256 \
  --receipt-sha256 FREEZE_RECEIPT_SHA256 \
  --registration-sha256 FULL_REGISTRATION_SHA256 \
  --python EXACT_REGISTERED_ENVIRONMENT_PYTHON \
  --output-directory FRESH_EXTERNAL_DIRECTORY
```

`--receipt` may explicitly select another pinned receipt path; its only default
is `CHECKPOINT/FREEZE_RECEIPT.json`. `--optimized` propagates optimized mode to
both custodians and route drivers. There are no workspace or interpreter
defaults. The user must supply the registered environment interpreter whose
runtime the existing registered driver verifies.

Before loading guard/review modules, the helper verifies the mandatory manifest
hash, complete exact file membership, all member byte counts/hashes, canonical
paths, symlink rejection and its own executing-source pin against the archive.
It then uses the pinned source guard and `review_completed.py` to verify both
archived actual executions. No new source callback can be launched unless that
actual baseline review succeeds. Fabricated or missing archived execution
evidence is rejected by the real review function.

The helper invokes the existing bounded custodian twice, preserving its
explicit registration/readback pins, interpreter, route universe and nonroot
external resource enforcement. New outputs remain in the fresh external
directory. It reviews those actual receipts and compares exact canonical
numerical payloads and `stable_frame` values with the archived originals.
Runtime/resource metadata is omitted only by the existing reviewed
`stable_frame` definition. Every scientific interval, model/remainder field and
semantic validation remains part of the exact comparison.

`PASS_EXACT_FRESH_REPLAY` means the archived scientific state reproduced. A
negative or unresolved archived scientific state remains negative or unresolved
and is explicitly recorded; replay never upgrades it. Changed numerical values
produce a failure report while retaining both fresh outputs. Execution failure
retains the custodians' logs, partial source journals and failure evidence.
The whole frozen archive is reverified before each route and after review.

Preparation evidence uses only manufactured standard-library mock custodians,
guards and reviewers, not production scientific source callbacks. Eleven
controls pass in normal and optimized Python, covering exact replay, optimized
replay, preservation of negative science, changed-value failure, manifest/payload
mutation, extra files, wrong pins, pre-existing output and rejection before any
new route when actual archived evidence is absent. Actual source evaluation and
fresh scientific replay are root's later registered execution.
