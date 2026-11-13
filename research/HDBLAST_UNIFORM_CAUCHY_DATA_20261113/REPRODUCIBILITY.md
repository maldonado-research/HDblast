# Closed source and exact replay protocol

Use Python 3.12.14 and the exact hash-pinned wheels in `requirements-checks.txt`. The new algebra uses exact rational polynomial/matrix witnesses. The numerical ledger uses integer and rational comparisons with formal decimal exponents; it does not round tiny amplitudes to floating-point zero.

The executable checks audit the written analytic calculation's finite identities and numerical allocations. They do not prove the continuum Sobolev/composition theorems, the retained background premises, mathematical novelty or the cosmological hypothesis.

## Authenticate before importing

The closed checker requires a trusted external SHA256 pin for MANIFEST.json. The reviewed pin is in `.github/workflows/uniform-cauchy-data-20261113.yml`, outside this capsule. Obtain it from a trusted reviewed commit or publication handoff. Computing a pin from an untrusted local manifest does not authenticate that manifest.

The parser rejects duplicate JSON keys and nonfinite constants. The checker authenticates every required code, proof, data and source byte; rejects symlinks and extra or missing files; and imports the evaluator only afterward. All project snapshots are fixed at base `a41940ed60dbda4808ea7286890443a931a44a74`. SOURCE_PINS.json and INPUT_PROVENANCE.json record path, commit, byte count, Git blob ID and SHA256. Snapshot creation executes no source producer. Literal producer scales and the model polynomial are checked from its syntax tree, while the prior seed/constraint and flow estimates are separately authenticated inputs.

## Execute

From the repository root, set `CAPSULE=research/HDBLAST_UNIFORM_CAUCHY_DATA_20261113` and set `MANIFEST_SHA256` to the independently trusted external pin. Use fresh output paths and directories outside the capsule:

```bash
python3.12 -m venv /tmp/hdblast-cauchy-check-env
/tmp/hdblast-cauchy-check-env/bin/python -m pip install --require-hashes --only-binary=:all: -r "$CAPSULE/requirements-checks.txt"
/tmp/hdblast-cauchy-check-env/bin/python -I -B "$CAPSULE/checks/check_checkpoint.py" --manifest-sha256 "$MANIFEST_SHA256" --output /tmp/hdblast-cauchy-normal.json
/tmp/hdblast-cauchy-check-env/bin/python -O -I -B "$CAPSULE/checks/check_checkpoint.py" --manifest-sha256 "$MANIFEST_SHA256" --output /tmp/hdblast-cauchy-optimized.json
cmp /tmp/hdblast-cauchy-normal.json /tmp/hdblast-cauchy-optimized.json
/tmp/hdblast-cauchy-check-env/bin/python -I -B "$CAPSULE/checks/run_controls.py" --manifest-sha256 "$MANIFEST_SHA256" --output-dir /tmp/hdblast-cauchy-auth-fresh
/tmp/hdblast-cauchy-check-env/bin/python -I -B "$CAPSULE/checks/run_math_controls.py" --manifest-sha256 "$MANIFEST_SHA256" --output-dir /tmp/hdblast-cauchy-math-fresh
```

Choose new temporary paths for another run. Reports and bytecode caches must remain outside the capsule. Every required check uses an explicit condition that remains enabled under `-O`.

## Finite witnesses and fault controls

The canonical oracle checks a nonzero curved positive spatial metric, nonzero momenta, affine lapse and their two-direction jets against the forward half-action canonical definitions. An independent five-dimensional Christoffel calculation solves twenty rational equations: five unknown initial gauge velocities, each through the coefficients `1,e,f,ef`. It retains the metric dependence in `g^BC GammaBar^A_BC` and compares the complete background-subtracted mixed velocity.

The chart oracle combines the actual square-root operation with a non-diagonal nonlinear polynomial witness. Direct parameter expansion is compared with the full velocity Hessian and its polarized form, including the third state differential. This tests a universal chain identity on finite witnesses; it does not approximate or numerically integrate the actual flow.

The arithmetic suite checks the stated graded background multipliers, determinant compositions, monomial counts, raw norm conversion, root/flow derivative budget, Sobolev chain count, initial patch and final one-/two-copy joins. Source-to-ledger interfaces bind the original admitted family and the actual retained frame and flow certificates.

Five authentication faults are rejected before evaluator import. Mathematical controls instead change substantive claims, rehash every file in a fresh complete capsule, and invoke a new full checker with that fixture's new manifest pin. Rejection must reach a named mathematical check, rather than the original capsule hash mismatch. These fixtures are deliberately invalid data.

CI repeats the authenticated normal/optimized reports and both sets of fault controls. Final counts and actual head, merge, CI and live-page receipts belong in the external handoff after publication. This premerge document does not assert those later outcomes.
