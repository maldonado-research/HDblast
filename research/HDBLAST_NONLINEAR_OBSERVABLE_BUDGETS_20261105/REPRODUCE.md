# Reproduce the algebra and source checks

Use Python 3.12 with its standard library. No numerical PDE solver, downloaded
package, external field/source oracle, or prior campaign is required.

From the repository root, obtain the trusted manifest SHA-256 from the pinned
`nonlinear-observable-budgets-20261105.yml` workflow at the commit being checked.
Replace `TRUSTED_MANIFEST_SHA256` below with that value. Use fresh absolute
output directories outside the capsule:

```
python -I -S -B research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/checks/run_checkpoint.py --manifest-sha256 TRUSTED_MANIFEST_SHA256 --output-dir /tmp/hdblast-budget-replay-new
python -I -S -B research/HDBLAST_NONLINEAR_OBSERVABLE_BUDGETS_20261105/checks/run_controls.py --manifest-sha256 TRUSTED_MANIFEST_SHA256 --output-dir /tmp/hdblast-budget-controls-new
```

The first command authenticates the entire closed inventory, runs geometry,
odd-remainder and observer-flow checkers in normal and optimized Python modes,
and compares all six produced reports byte-for-byte with their pinned
references. Producer stdout must equal the actual report bytes. It rechecks
the source inventory afterward and writes `VERIFICATION.json` externally.

The second command runs two positive authentication-only controls and six
actual rejections: wrong manifest pin in both Python modes, changed selected
analytic input, changed reference, changed producer, and occupied output.
These are source/CLI transport controls, not mathematical proof mutants.
Rejected controls execute no scientific producer.

The finite algebra counts per mode are geometry **115**, odd remainder **63**,
and observer flow **176**. These accompany the continuous written arguments;
they do not evaluate nonlinear field budgets, prove retained PDE premises,
validate a common observer tube, or constitute external peer review.

An archive of this directory includes only this working checkpoint. Selected
older documents in `inputs/` are retained analytic premises with immutable
source URLs and byte/Git-blob pins in `INPUT_PROVENANCE.json`. Authenticating
them does not rerun or re-prove their original results. Private review logs,
editorial history, old complete capsules and publication-account credentials
are excluded.
