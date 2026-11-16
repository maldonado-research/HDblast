# Reproduce the common-observer replay

The capsule is `research/HDBLAST_OBSERVER_CLOCK_20261116`. The workflow
is `.github/workflows/observer-clock-20261116.yml`, named
**Common observer and proper-clock replay**. Use a trusted reviewed
revision of that workflow as the external manifest-pin source. The
capsule cannot provide its own trust anchor by hashing its downloaded
`MANIFEST.json` and declaring that hash trusted.

## Dependencies and fresh outputs

Use Python 3.12.14 and the exact binary-wheel hashes in
`requirements-checks.txt` (SymPy 1.14.0 and mpmath 1.3.0). For example,
from the repository root:

```bash
python3.12 -m venv /tmp/hdblast-observer-env
/tmp/hdblast-observer-env/bin/python -m pip install --require-hashes --only-binary=:all: -r research/HDBLAST_OBSERVER_CLOCK_20261116/requirements-checks.txt
```

Use that environment's Python for the commands below.

Set the capsule path and copy the 64-character lowercase hexadecimal
manifest pin from the trusted workflow revision. Choose an output
parent outside the capsule. The commands below create fresh report
names; the checker refuses existing outputs and symlink parents.

```bash
OBSERVER_CAPSULE="$PWD/research/HDBLAST_OBSERVER_CLOCK_20261116"
read -r OBSERVER_PIN  # paste the trusted 64-character workflow digest
OBSERVER_RUNS=$(mktemp -d /tmp/hdblast-observer-replay.XXXXXX)

python -I -B "$OBSERVER_CAPSULE/checks/check_checkpoint.py" \
  --manifest-sha256 "$OBSERVER_PIN" \
  --output "$OBSERVER_RUNS/normal.json"
python -I -B -O "$OBSERVER_CAPSULE/checks/check_checkpoint.py" \
  --manifest-sha256 "$OBSERVER_PIN" \
  --output "$OBSERVER_RUNS/optimized.json"
cmp "$OBSERVER_RUNS/normal.json" "$OBSERVER_RUNS/optimized.json"
```

The executed report is authoritative for the actual output and count.
The externally pinned campaign and GitHub workflow expose the same code.

## Authentication and mathematical rejection controls

```bash
python -I -B "$OBSERVER_CAPSULE/checks/run_controls.py" \
  --manifest-sha256 "$OBSERVER_PIN" \
  --output-dir "$OBSERVER_RUNS/authentication"
python -I -B "$OBSERVER_CAPSULE/checks/run_math_controls.py" \
  --manifest-sha256 "$OBSERVER_PIN" \
  --output-dir "$OBSERVER_RUNS/mathematics-normal"
python -I -B -O "$OBSERVER_CAPSULE/checks/run_math_controls.py" \
  --manifest-sha256 "$OBSERVER_PIN" \
  --output-dir "$OBSERVER_RUNS/mathematics-optimized"
```

The authentication runner first checks a positive baseline and then
tests independent source-byte, extra-file, missing-file, source-symlink
and stale-external-pin faults. Each rejected child must leave no success
report.

The mathematical runner authenticates a positive baseline before
importing the own pinned mutation roster. It creates one fresh full
capsule fixture and manifest for each single-claim mutation and invokes
a fresh full checker subprocess. Rejection must reach a substantive
named own-oracle identity or interface guard. Schema, parser, type and
engine failures are excluded from that evidence. The optimized run
explicitly preserves optimization in its children. Control receipts may
therefore record different optimization metadata even when their
mathematical case records agree.

## What is authenticated and replayed

The external pin authenticates the entire local closed capsule roster
except the manifest's own entry. Every declared member's size and
SHA-256 are recomputed. JSON parsing rejects duplicate keys and
nonfinite values. Member paths must be canonical relative paths;
all source files and directories must be regular and free of symlinks.
The report path is checked before evaluator import and created
exclusively only after success.

The evaluator authenticates the exact selected immutable source roster,
Git blob identities, qualified direct originating-member joins
and input provenance as data. It checks the retained scientific
conventions and typed new status before importing
`checks/observer_oracles.py`, whose exact bytes have their own guard.
No ancestor producer or checker is executed. The selected source counts are 25 snapshots and 18 direct joins.
The complete evaluator has 234 named checks, including 62 own oracle
checks, and the substantive wrong-claim roster contains 50 cases.

The own finite witnesses cover their stated clock, proper normalization,
ODE variation, full moving scalar/tensor/metric derivatives, determinant
and frame identities. They do not evaluate the actual sixteen-field
coefficient tensors or prove the continuum observer construction.
The written proofs retain that responsibility.
