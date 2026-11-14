# Reproduction and independently rehashed faults

Use Python3.12.14 and the exact hash-pinned binary wheels in
`requirements-checks.txt`. From the repository root:

```sh
python -m venv /tmp/hdblast-limit-env
/tmp/hdblast-limit-env/bin/python -m pip install --require-hashes --only-binary=:all: -r research/HDBLAST_NONLINEAR_PICARD_LIMIT_20261114/requirements-checks.txt
```

Take `PIN` from the four literal external digests in
`.github/workflows/nonlinear-picard-limit-20261114.yml` at the trusted
commit. A locally computed untrusted manifest digest is not an
authenticity anchor.

```sh
/tmp/hdblast-limit-env/bin/python -I -B research/HDBLAST_NONLINEAR_PICARD_LIMIT_20261114/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/hdblast-limit-normal.json
/tmp/hdblast-limit-env/bin/python -O -I -B research/HDBLAST_NONLINEAR_PICARD_LIMIT_20261114/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/hdblast-limit-optimized.json
cmp /tmp/hdblast-limit-normal.json /tmp/hdblast-limit-optimized.json
/tmp/hdblast-limit-env/bin/python -I -B research/HDBLAST_NONLINEAR_PICARD_LIMIT_20261114/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/hdblast-limit-auth-fresh
/tmp/hdblast-limit-env/bin/python -O -I -B research/HDBLAST_NONLINEAR_PICARD_LIMIT_20261114/checks/run_math_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/hdblast-limit-math-fresh
```

All outputs must be fresh and outside the capsule. The checker first
authenticates its closed byte roster and strict JSON before importing
the evaluator. The evaluator authenticates29 immutable source snapshots,
their provenance and ancestor joins before mathematical imports.
Snapshot producers are never imported or executed. Symlinks, missing
files and extra files are rejected.

The convergence oracle checks auxiliary nonlinear source derivatives,
joint mean-value integration, distinct Picard/actual-solution difference
identities, weighted Volterra/factorial formulas and finite tail models.
It checks finite dependency bookkeeping through the strong lower rows
and the negative-order fourteenth row. Finite tails support algebraic
identities; they are not an infinite convergence theorem by themselves.

The extension oracle checks first-order symmetrization and the full
normal flux/characteristic polynomial with all four spatial directions
and sixteen fibers. It checks independent coefficient/fiber parity,
fourfold closure and meaningful trace rosters. It does not certify
Sobolev traces, weak conormal commutation, Friedrichs residual convergence
or top continuity; those are written analytic arguments.

The arithmetic module verifies the exact rational/exponent allocations,
unchanged admission, kernel/tail/stability constants and causal support
fractions. Huge powers remain symbolic; no practical nonlinear
wave-amplitude simulation or large floating-point expansion is run.

Five authentication faults test changed bytes, roster, symlinks and pin.
Mathematical faults receive fresh complete capsules, altered claim
objects, new manifests and new external digests; rejection must reach
the exact mathematical evaluator. The harness propagates Python's
optimization level to checker children and records it. Direct normal
and optimized checkpoint reports must match.

The optional static support figure uses Matplotlib3.10.8 and is not
required by the algebra checker. To regenerate equivalent geometry
with an installed Matplotlib, writing to a fresh external file:

```sh
MPLCONFIGDIR=/tmp/hdblast-plot-cache python -I -B research/HDBLAST_NONLINEAR_PICARD_LIMIT_20261114/figures/generate_causal_support.py --output /tmp/hdblast-support.svg
```

Fonts/platforms can change SVG bytes. The supplied SVG is byte-pinned
by the capsule; its coordinate fractions and support interpretation
are independently reviewable in the integrated theorem and generator.
