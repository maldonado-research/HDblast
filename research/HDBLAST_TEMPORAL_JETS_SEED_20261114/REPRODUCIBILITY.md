# Reproduction and authentication

Use Python 3.12.14 and the exact hash-pinned binary wheels in
`requirements-checks.txt`. From the repository root:

```sh
python -m venv /tmp/hdblast-temporal-env
/tmp/hdblast-temporal-env/bin/python -m pip install --require-hashes --only-binary=:all: -r research/HDBLAST_TEMPORAL_JETS_SEED_20261114/requirements-checks.txt
```

Take `PIN` from the four literal `--manifest-sha256` values in
`.github/workflows/temporal-jets-seed-20261114.yml` at the trusted commit.
The pin lives outside this capsule. A computed digest of an untrusted
manifest is not an authenticity anchor.

```sh
/tmp/hdblast-temporal-env/bin/python -I -B research/HDBLAST_TEMPORAL_JETS_SEED_20261114/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/hdblast-temporal-normal.json
/tmp/hdblast-temporal-env/bin/python -O -I -B research/HDBLAST_TEMPORAL_JETS_SEED_20261114/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/hdblast-temporal-optimized.json
cmp /tmp/hdblast-temporal-normal.json /tmp/hdblast-temporal-optimized.json
/tmp/hdblast-temporal-env/bin/python -I -B research/HDBLAST_TEMPORAL_JETS_SEED_20261114/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/hdblast-temporal-auth-fresh
/tmp/hdblast-temporal-env/bin/python -I -B research/HDBLAST_TEMPORAL_JETS_SEED_20261114/checks/run_math_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/hdblast-temporal-math-fresh
```

Each output directory must be fresh. The checker authenticates the closed
file roster, byte digests and strict JSON before importing the evaluator.
The evaluator checks the immutable base snapshots and provenance before
importing mathematical modules. Symlinks and extra files are rejected.
Snapshot producers are neither imported nor executed.

The exact finite replay checks the nonlinear binomial recursion through
row twelve, isolation of the leading acceleration coefficient, a manufactured nonlinear
solution, bundle/base connection correction and the nonzero-state horizontal
coefficient derivative. It checks all 225 degree-bounded Gaussian derivative
polynomials, every matching trace through fourteen, spectral weight
identities and evaluated integer/exponent admission inequalities.

Authentication controls change source bytes, the roster, symlink conditions
and the manifest pin. Mathematical controls change declared claims in fresh complete
capsules, regenerate their manifests, and run the same full checker with
the fresh external digest. Their rejection therefore tests the evaluator,
not a stale-manifest mismatch. Safety checks remain active under `-O`.

Finite identities and arithmetic do not certify Sobolev multiplication,
operator domains, strong spectral convergence or nonlinear existence.
Those are written analytic arguments with the retained premises and
independent internal review described in `AUDIT.md`.
