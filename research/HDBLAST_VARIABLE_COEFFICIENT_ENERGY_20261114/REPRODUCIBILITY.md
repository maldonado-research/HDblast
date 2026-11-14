# Reproduction and fresh fault controls

Use Python3.12.14 and the exact hash-pinned binary wheels in
`requirements-checks.txt`. From the repository root:

```sh
python -m venv /tmp/hdblast-variable-env
/tmp/hdblast-variable-env/bin/python -m pip install --require-hashes --only-binary=:all: -r research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/requirements-checks.txt
```

Take `PIN` from the four literal external manifest digests in
`.github/workflows/variable-coefficient-energy-20261114.yml` at a trusted
commit. Computing an untrusted manifest's digest does not authenticate it.

```sh
/tmp/hdblast-variable-env/bin/python -I -B research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/hdblast-variable-normal.json
/tmp/hdblast-variable-env/bin/python -O -I -B research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/hdblast-variable-optimized.json
cmp /tmp/hdblast-variable-normal.json /tmp/hdblast-variable-optimized.json
/tmp/hdblast-variable-env/bin/python -I -B research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/hdblast-variable-auth-fresh
/tmp/hdblast-variable-env/bin/python -I -B research/HDBLAST_VARIABLE_COEFFICIENT_ENERGY_20261114/checks/run_math_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/hdblast-variable-math-fresh
```

Outputs and directories must be fresh and outside the capsule. The closed
checker authenticates its complete byte roster and strict JSON before
importing the evaluator. The evaluator authenticates26 immutable base
snapshots and their provenance before importing mathematical modules.
No source producer is imported or executed. Symlinks and extra files
are rejected.

The finite oracle differentiates an auxiliary covariant energy density
and flux with a nontrivial split-preserving connection. It leaves
complementary boundary jets free, checks differentiated unisolated
normal rows through order fourteen and nonzero interior mixed terms,
enumerates normal triangularity and ordered commutator subsets, and
expands the exact predecessor difference-source identity. These
manufactured fields are algebraic witnesses, not evolved physical data.

The arithmetic ledger checks exact rational principal margins, finite
rotation/product counts, the six energy/recovery constants, all thirteen
normal stages and the conditional smooth-iteration admission. Nested
powers in the amplitude proof are compared symbolically; enormous
decimal integers and physical amplitudes are not computed.

Five authentication cases change source bytes, roster, symlinks and the
manifest pin. Mathematical cases alter declared claims in fresh complete
capsules, regenerate each manifest and invoke the full checker with its
fresh external digest. Rejection must reach the intended mathematical
oracle. The mathematical harness propagates the calling Python
optimization level to child checkers and records it; the commands above
use ordinary children. Direct normal/-O checkpoint reports must match.

These finite runs verify identities, counts and arithmetic. Continuum
rotation/covariant norm bridges, Sobolev products, smooth frozen
solvability and limit regularity require the written analytic arguments
and their stated premises. No nonlinear-limit or physical confirmation
claim follows from passing a replay.
