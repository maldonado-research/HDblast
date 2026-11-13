# Reproduce the authenticated all-mode finite replay

The written analytic proofs are in the four companion notes. The runner
checks finite algebra and coefficient/rational accounting; it is not an
analytic existence, regularity or physical-validation engine.

Use Python 3.12.14 and the hash-pinned SymPy 1.14.0 / mpmath 1.3.0 wheels:

```
python -m venv /tmp/hdblast-all-mode-env
/tmp/hdblast-all-mode-env/bin/python -m pip install --require-hashes --only-binary=:all: -r research/HDBLAST_ALL_MODE_INVERSE_20261113/requirements-checks.txt
```

Obtain the independent literal `PIN` from the checked workflow
`.github/workflows/all-mode-inverse-20261113.yml`. A hash calculated from
an unauthenticated replacement capsule is not an independent trust pin.
From the repository root, with fresh output paths outside the capsule:

```
/tmp/hdblast-all-mode-env/bin/python -I -B research/HDBLAST_ALL_MODE_INVERSE_20261113/checks/check_checkpoint.py --manifest-sha256 PIN --output /tmp/hdblast-all-mode-normal.json
/tmp/hdblast-all-mode-env/bin/python -O -I -B research/HDBLAST_ALL_MODE_INVERSE_20261113/checks/check_checkpoint.py --manifest-sha256 PIN --output /tmp/hdblast-all-mode-optimized.json
cmp /tmp/hdblast-all-mode-normal.json /tmp/hdblast-all-mode-optimized.json
/tmp/hdblast-all-mode-env/bin/python -I -B research/HDBLAST_ALL_MODE_INVERSE_20261113/checks/run_controls.py --manifest-sha256 PIN --output-dir /tmp/hdblast-all-mode-controls-fresh
```

The protocol authenticates its closed file roster, manifest and all bytes
before importing evaluator code. It rejects symlinks, missing/unlisted or
changed files, stale pins, duplicate JSON keys and nonfinite JSON. Exact
dependency wheels are authenticated by the hash-required install, and
runtime versions are checked.

The five mathematical modules cover joint angular/radial and endpoint
symbols; the generic scalar square-plus-remainder identity; actual
polynomial coefficient envelopes; exact scalar/vector graph and inverse
allocations; and complete momentum normal operators in every sector.
The actual potential coefficients are also tied to the authenticated
source snapshot. All calculations are exact rational, polynomial or
integer operations. Huge exponential ceilings are not numerically sampled.

Failure controls alter norm powers, principal factors, operator terms,
coefficient budgets, squared-gap accounting and graph/bridge/inverse
allocations. Mathematical cases have refreshed whole-capsule manifests
and run in fresh child processes, so their rejection must occur in a
recomputed mathematical suite rather than at a stale hash. Authentication
faults are tested separately. Normal and optimized canonical reports must
be byte-identical.

Local algebra does not alone prove a global spectral gap or an inverse
norm. The all-mode quantifiers, positive-operator estimates, spectral
absorption, tensor norm bridge, density convention and dotted-domain
passage are the written analytic arguments accepted in the internal
review. Neither those proofs nor these finite checks supply a nonlinear
amplitude radius or physical confirmation.
