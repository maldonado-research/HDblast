# Reproduce the finite initial-output calculation

Use Python 3.12 or later. The public replay uses only the standard
library, exact rational arithmetic and finite polynomial operations.
It does not run the historical PDE campaigns, solve the normal
operator numerically or propagate a nonlinear Einstein–scalar family.

Obtain this entire directory at the reviewed commit. Use the external
MANIFEST.json SHA256 pin recorded in
`.github/workflows/weighted-initial-output-20261106.yml` at that same
commit; a manifest recomputed from untrusted inputs is not a trust anchor.
From the repository root, set PIN to that literal 64-character value
and choose three fresh output locations outside the capsule:

```sh
python -I -S -B research/HDBLAST_WEIGHTED_INITIAL_OUTPUT_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/initial-output-normal.json
python -O -I -S -B research/HDBLAST_WEIGHTED_INITIAL_OUTPUT_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/initial-output-optimized.json
cmp /tmp/initial-output-normal.json /tmp/initial-output-optimized.json
python -I -S -B research/HDBLAST_WEIGHTED_INITIAL_OUTPUT_20261106/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/initial-output-controls-fresh
```

The controls directory must not already exist. Reports must stay outside
the closed capsule. The canonical reports should match byte for byte;
the control report should record its passing baseline and six intended
rejections. Source/code mutations that also replace the trusted pin
are outside the authenticity guarantee; the arithmetic control instead
checks that an internally rehashed but mathematically inconsistent
fixture is rejected by a specific calculation.

Read SCALAR_WEIGHTED_OUTPUT.md and ACTUAL_INITIAL_SOURCE.md for the
continuous arguments and assumptions. The finite tests authenticate
source provenance and arithmetic, not existence, external originality,
physical evolution, an amplitude radius or a Big Bang mechanism.
