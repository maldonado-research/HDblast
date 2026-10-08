# Prepared all-node comparison execution

The candidate directory is
`candidate/HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON`. Preparation
creates only a local fabricated registration. Root owns the independent final
review, final registration, public freeze and externally authenticated PUBLIC_GO.
No command below creates publication or a public authorization receipt.

The prepared environment uses
`/workspace/hdblast-research-work/environment/bin/python` (Python3.12.14,
python-flint0.9.0, sympy1.14.0, mpmath1.3.0), Linux and libseccomp.so.2. An
equivalent environment may use another interpreter path. Runtime byte hashes
are observation provenance, while exact versions and the native guard ABI are
execution requirements. Use the existing checkout; no new Git worktree is needed.

## Local preparation and fabricated integration

The following shell variables are illustrative non-secret paths. The current
registration SHA is printed by `prepare_candidate.py`; copy that literal hash
to `COMPARISON_REG_SHA` before launching. Every output directory must be new.

```bash
COMPARISON_WORK=/workspace/hdblast-research-work/continuation-state-comparison-20261007
COMPARISON_CANDIDATE="$COMPARISON_WORK/candidate/HDBLAST_CHECKPOINT_20261007_STORED_STATE_COMPARISON"
COMPARISON_PY=/workspace/hdblast-research-work/environment/bin/python
"$COMPARISON_PY" -B "$COMPARISON_WORK/execution/prepare_candidate.py" \
  --work "$COMPARISON_WORK" \
  --prior /workspace/HDblast/research/HDBLAST_CHECKPOINT_20261007_BD_PREHISTORY \
  --candidate "$COMPARISON_CANDIDATE"
COMPARISON_REG_SHA=REPLACE_WITH_PRINTED_LITERAL_REGISTRATION_SHA256
"$COMPARISON_PY" -B "$COMPARISON_CANDIDATE/execution/bounded_launcher.py" \
  --candidate "$COMPARISON_CANDIDATE" --registration-sha256 "$COMPARISON_REG_SHA" \
  --output "$COMPARISON_WORK/execution/fresh-fabricated-normal" --fabricated-only
"$COMPARISON_PY" -B -O "$COMPARISON_CANDIDATE/execution/bounded_launcher.py" \
  --candidate "$COMPARISON_CANDIDATE" --registration-sha256 "$COMPARISON_REG_SHA" \
  --output "$COMPARISON_WORK/execution/fresh-fabricated-optimized" --fabricated-only --optimized
"$COMPARISON_PY" -B "$COMPARISON_CANDIDATE/execution/readback.py" \
  --candidate "$COMPARISON_CANDIDATE" --registration-sha256 "$COMPARISON_REG_SHA" \
  --output "$COMPARISON_WORK/execution/fresh-fabricated-normal" --fabricated-only \
  --compare-output "$COMPARISON_WORK/execution/fresh-fabricated-optimized"
```

Full fabricated integration manufactures four deterministic19-member capsules
with the complete49,152-node geometry and uses the prior explicit finite
polynomial source fixture. It opens no historical NPY/NPZ, receives no actual
repository root and invokes no physical source callback. Manufactured arrays
and nodes have separate counters. The 18 fabricated diagnostics do not claim
intersection with old physical results.

The kernel/output security suite and registration/payload controls are distinct
from full integration:

```bash
HDBLAST_COMPARISON_TEST_OUTPUT="$COMPARISON_WORK/execution/kernel-controls-normal" \
  "$COMPARISON_PY" -B "$COMPARISON_CANDIDATE/execution/test_bounded_launcher.py"
HDBLAST_COMPARISON_TEST_OUTPUT="$COMPARISON_WORK/execution/kernel-controls-optimized" \
  "$COMPARISON_PY" -B -O "$COMPARISON_CANDIDATE/execution/test_bounded_launcher.py"
"$COMPARISON_PY" -B "$COMPARISON_CANDIDATE/execution/test_controls.py" \
  --candidate "$COMPARISON_CANDIDATE" --registration-sha256 "$COMPARISON_REG_SHA" \
  --evidence "$COMPARISON_WORK/execution/fresh-guard-payload-controls" \
  --output "$COMPARISON_WORK/execution/fresh-fabricated-normal"
```

## Actual execution after root's separate public GO

These commands remain unavailable until root supplies the external immutable
registration and PUBLIC_GO pins. The receipt must live outside the candidate.
The actual repository root must contain all exact pinned InputSpec paths.

```bash
"$COMPARISON_PY" -B "$COMPARISON_CANDIDATE/execution/bounded_launcher.py" \
  --candidate "$COMPARISON_CANDIDATE" --registration-sha256 "$COMPARISON_REG_SHA" \
  --public-go "$COMPARISON_PUBLIC_GO" --public-go-sha256 "$COMPARISON_PUBLIC_GO_SHA" \
  --repository-root "$COMPARISON_REPOSITORY_ROOT" --output "$COMPARISON_ACTUAL_NORMAL"
"$COMPARISON_PY" -B -O "$COMPARISON_CANDIDATE/execution/bounded_launcher.py" \
  --candidate "$COMPARISON_CANDIDATE" --registration-sha256 "$COMPARISON_REG_SHA" \
  --public-go "$COMPARISON_PUBLIC_GO" --public-go-sha256 "$COMPARISON_PUBLIC_GO_SHA" \
  --repository-root "$COMPARISON_REPOSITORY_ROOT" --output "$COMPARISON_ACTUAL_OPTIMIZED" --optimized
"$COMPARISON_PY" -B "$COMPARISON_CANDIDATE/execution/readback.py" \
  --candidate "$COMPARISON_CANDIDATE" --registration-sha256 "$COMPARISON_REG_SHA" \
  --public-go "$COMPARISON_PUBLIC_GO" --public-go-sha256 "$COMPARISON_PUBLIC_GO_SHA" \
  --output "$COMPARISON_ACTUAL_NORMAL" --compare-output "$COMPARISON_ACTUAL_OPTIMIZED"
```

For an offline standalone package, retain all four smallcapsule NPZ inputs, all
four pinned original archives, the producer and manifest files, and the prior
static-audit file at their exact repository-relative paths. Those paths are in
`provenance/INPUT_SPEC.json`; the combined old input bytes are streamed rather
than simultaneously resident. The new candidate includes its exact source
builder copies, prior successful source budget, both old diagnostic DATA files,
layout evidence and inventories. Third-party PDFs are unnecessary for execution.

## Results and limits

`ENTRY_RECEIPT.json` pins all nine worker payload files.
`EXECUTION.json` records the custodian's authoritative wait4 outcome, wall/RSS/CPU
and output acceptance evidence. The two empty actual-attempt journals are
expected only in fabricated mode. Actual mode requires22 source attempts and20
selected-array attempts before construction/decoding. `NODE_TARGETS.jsonl.gz`,
`SOURCE_CERTIFICATE.json`, `DIAGNOSTIC_CROSSCHECK.json`, `DATA.json` and
`SCIENCE_SUMMARY.json` retain exact scientific outputs.

The fixed source112/Arb512/momentum832/dyadic96 complete-L1<=1e-18 parameters
are not adaptive. Worker resources are900seconds wall/CPU,512MiB address space,
128MiB aggregate accepted output and per-file FSIZE, nonroot, no child processes,
no sockets/exec/io_uring and one thread. Aggregate output is monitored and
checked finally; FSIZE is per-file. Failed output directories and diagnostics
are preserved and never reused for a later attempt.

Finite weighted incoming-state comparison does not establish full continuous
pressure/contact certification, metric success, higher-dimensional origin or
novelty. Their unchanged statuses are recorded in every scientific summary.
