# Closed-capsule replay

The executable checks verify finite identities and numerical ledgers used
by the written analytic proofs. They do not replace the continuum proofs,
establish the retained background existence premises, prove novelty or
test the hypothesis against astronomical data.

Use Python 3.12.14 and the exact hash-pinned wheels in
`requirements-checks.txt`. The SymPy checks use exact rational expressions;
seed bounds use outward-rounded rational Taylor majorants, without bump
sampling. The seed producer is an explicitly authenticated prior program,
not an independent reimplementation. Symbolic oracle components and
finite weight/constant ledgers are independently recomputed in the new
modules.

## Authentication before import

The checker requires an external SHA256 pin for MANIFEST.json. The reviewed
pin is in the repository workflow
`.github/workflows/nonlinear-constraint-ball-20261113.yml`, outside this
closed capsule. Obtain that pin from a trusted reviewed commit or handoff.
Recomputing it from an untrusted local manifest does not authenticate it.

The parser rejects duplicate keys and nonfinite JSON values. The checker
authenticates the manifest and exact roster, rejects symlinks and extra
files, verifies each file's bytes, then loads the mathematical modules.
Source snapshots include original commit, path, Git blob ID, size and
SHA256 in SOURCE_PINS.json and INPUT_PROVENANCE.json. All project inputs
come from base `4c1c47ce4cf7f0106e3a462bfd48618b04ae3bdf`. The source checker
also tests meaningful model conventions; authentication alone is not a
mathematical test.

## Run

From the repository root, set `CAPSULE` to
`research/HDBLAST_NONLINEAR_CONSTRAINT_BALL_20261113` and set
`MANIFEST_SHA256` to the independently trusted external pin. Use a fresh
temporary directory outside the capsule for all output:

```bash
python3.12 -m venv /tmp/hdblast-nonlinear-check-env
/tmp/hdblast-nonlinear-check-env/bin/python -m pip install --require-hashes --only-binary=:all: -r "$CAPSULE/requirements-checks.txt"
/tmp/hdblast-nonlinear-check-env/bin/python -I -B "$CAPSULE/checks/check_checkpoint.py" --manifest-sha256 "$MANIFEST_SHA256" --output /tmp/hdblast-nonlinear-normal.json
/tmp/hdblast-nonlinear-check-env/bin/python -O -I -B "$CAPSULE/checks/check_checkpoint.py" --manifest-sha256 "$MANIFEST_SHA256" --output /tmp/hdblast-nonlinear-optimized.json
cmp /tmp/hdblast-nonlinear-normal.json /tmp/hdblast-nonlinear-optimized.json
/tmp/hdblast-nonlinear-check-env/bin/python -I -B "$CAPSULE/checks/run_controls.py" --manifest-sha256 "$MANIFEST_SHA256" --output-dir /tmp/hdblast-nonlinear-auth-fresh
/tmp/hdblast-nonlinear-check-env/bin/python -I -B "$CAPSULE/checks/run_math_controls.py" --manifest-sha256 "$MANIFEST_SHA256" --output-dir /tmp/hdblast-nonlinear-math-fresh
```

The output paths and control directories must be fresh. Do not write
reports or bytecode caches inside the capsule. For a new run, use new
temporary paths. All required checks use explicit conditions that remain
enabled under Python `-O`.

## What is recomputed

- A curved positive four-dimensional metric germ exercises the full
  inverse, determinant, scalar curvature and canonical matter Hessians,
  including a nonzero state connection relative to a fixed chart.
- Fixed-density and coordinate-density momentum formulas are compared,
  including both `h DPi` and `Pi Dh` and scalar momentum/gradient terms.
- Every finite leading derivative allocation through the H12 source
  order is enumerated, including derivatives on w^-2. The remaining
  weight and derivative budgets are checked exactly.
- The actual four-sector seed's radial jets are regenerated from the
  pinned ODE/bump producer and charged to the full tensor and graded
  envelopes.
- Embedding, product, normalization, background, correction, whole-ball
  Hessian and contraction inequalities are evaluated with exact integers
  and rational numbers. Large constants are stored as exponents, without
  rounding a tiny radius to zero in floating point.

Authentication faults deliberately corrupt the protocol. Mathematical
controls instead alter substantive coefficients or inequalities, refresh
the entire manifest of a separate full capsule, and invoke the closed
checker with that fresh external manifest pin. Their rejection must come
from a named mathematical test, not the original capsule's hash mismatch.
The records are negative fixtures, not alternate model data.

The CI workflow repeats normal and optimized evaluation and both sets of
fault controls. Final counts, head/merge commits, preservation, CI and
live-site receipts are saved in the external handoff after publication;
this premerge document does not assert those future outcomes.
