# Reproduce the effective observer checkpoint

The capsule is `research/HDBLAST_EFFECTIVE_OBSERVER_20261116`. Its
workflow is `.github/workflows/effective-observer-20261116.yml`, named
**Effective observer domain replay**. Obtain the manifest SHA-256 from a
trusted reviewed revision of that workflow, outside the capsule. Hashing
the downloaded manifest does not by itself establish its trust.

## Runtime and fresh outputs

The recorded runtime is Python 3.12.14. The own oracle and four wrappers
use only the Python standard library. `requirements-checks.txt` is a
comment-only dependency record: no external wheels, symbolic package or
pip installation is required for this checkpoint. Prior capsules keep
their original dependency records unchanged.

From the repository root, set the capsule path, copy the trusted
64-character lowercase hexadecimal pin, and create an outside output
parent. Replace the pin placeholder before running:

```bash
EFFECTIVE_CAPSULE="$PWD/research/HDBLAST_EFFECTIVE_OBSERVER_20261116"
EFFECTIVE_PIN='<copy the trusted workflow manifest pin>'
EFFECTIVE_RUNS=$(mktemp -d /tmp/hdblast-effective-replay.XXXXXX)

python -I -B "$EFFECTIVE_CAPSULE/checks/check_checkpoint.py" \
  --manifest-sha256 "$EFFECTIVE_PIN" \
  --output "$EFFECTIVE_RUNS/normal.json"
python -I -B -O "$EFFECTIVE_CAPSULE/checks/check_checkpoint.py" \
  --manifest-sha256 "$EFFECTIVE_PIN" \
  --output "$EFFECTIVE_RUNS/optimized.json"
cmp "$EFFECTIVE_RUNS/normal.json" "$EFFECTIVE_RUNS/optimized.json"
```

The checker requires fresh output files outside the capsule, with
existing regular parents and no symlinks. It creates a report
exclusively only after successful evaluation. Existing reports are not
overwritten.

## Rejection controls

```bash
python -I -B "$EFFECTIVE_CAPSULE/checks/run_controls.py" \
  --manifest-sha256 "$EFFECTIVE_PIN" \
  --output-dir "$EFFECTIVE_RUNS/authentication"
python -I -B "$EFFECTIVE_CAPSULE/checks/run_math_controls.py" \
  --manifest-sha256 "$EFFECTIVE_PIN" \
  --output-dir "$EFFECTIVE_RUNS/mathematics-normal"
python -I -B -O "$EFFECTIVE_CAPSULE/checks/run_math_controls.py" \
  --manifest-sha256 "$EFFECTIVE_PIN" \
  --output-dir "$EFFECTIVE_RUNS/mathematics-optimized"
```

The authentication controls test independent source-byte, extra-file,
missing-file, source-symlink and stale-pin faults after a positive
baseline. Mathematical controls use the authenticated own mutation
roster. Each changes one typed claim, creates a fresh complete fixture
manifest and starts a fresh full checker subprocess. A substantive
named own-oracle identity or interface rejection is required; schema,
parser, type or engine failures do not count. Optimized parents must
preserve optimization in their children. Control receipts may therefore
contain different optimization metadata while agreeing on their cases.

These commands describe the protocol. [RESULTS.md](RESULTS.md) records
the frozen check/control rosters and acceptance conditions. Actual
final-capsule run evidence comes from the workflow and separate root
receipts, not an anticipated mutation count or an earlier capsule's run.

## Authentication and finite scope

Before evaluator import, the checker authenticates the external pin,
canonical relative member paths, exact closed regular-file roster,
member sizes and SHA-256 values, and strict JSON. Duplicate keys,
nonfinite JSON constants, unsafe paths and symlinks reject. The output
path is checked before import.

`checks/evaluate_effective_observer.py` authenticates the selected
immutable source roster, Git blob identities, direct originating-member
joins and provenance, then checks retained values and typed new claims.
Only afterward may it import the separately pinned own module
`checks/effective_oracles.py`. No ancestor scientific producer, checker
or oracle is imported or executed. Selected-member authentication is
distinct from re-running an entire ancestor capsule closure.

The own exact finite witnesses test the declared exponent, metric,
phase, curvature/Jacobi, causal-delay, determinant and scalar-envelope
implications. They do not evaluate the actual full coefficient tensors
or algorithmically prove the continuous physical field theorem,
Sobolev trace reserve, ODE continuation or sphere degree argument.
The mathematical proofs and the finite replay have separate roles.

Final manifest pins and execution receipts are kept outside the
capsule trust loop. [PROVENANCE.md](PROVENANCE.md) identifies the
immutable base, selected-source projection and original proof bytes.
