# Reproduce the high-order first-wave checkpoint

Python 3.12 or later with only the standard library is required. Obtain the
literal external manifest pin from `.github/workflows/high-order-first-wave-20261112.yml`
at the same reviewed commit and set PIN to that value. A pin recomputed from
untrusted inputs supplies no authenticity anchor. From the repository root,
use fresh paths outside the closed capsule:

```sh
python -I -S -B research/HDBLAST_HIGH_ORDER_FIRST_WAVE_20261112/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/high-order-normal.json
python -O -I -S -B research/HDBLAST_HIGH_ORDER_FIRST_WAVE_20261112/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output /tmp/high-order-optimized.json
cmp /tmp/high-order-normal.json /tmp/high-order-optimized.json
python -I -S -B research/HDBLAST_HIGH_ORDER_FIRST_WAVE_20261112/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir /tmp/high-order-controls-fresh
```

There are 221 named checks: 24 original source size/SHA256/Git-blob checks,
95 exact rotation/temporal polynomial checks and 102 finite norm/recurrence
arithmetic and derivative-allocation checks. Both modes must produce identical
canonical JSON. Authentication rejects changed/missing/extra files, symlinks,
stale pins, duplicate JSON keys and nonfinite constants before importing code.

The formal engine constructs the nine-dimensional scalar quadrupole module,
its six exact generator matrices and its rational sphere-moment Gram matrix.
It checks skew-adjointness and all ordered-word norm identities through
fourteen, the Killing tight frame, natural projected-gradient/metric identities,
and the vanishing translation-spin cross term at the north pole. A separate
indeterminate temporal polynomial verifies the required Leibniz cross factor.
These finite tests support the general natural tensor action proved in
SYMMETRY.md and HIGH_ORDER.md; they are not a componentwise simulation of the
sixteen-field PDE. Finite integer allocations check every recurrence index
range without enumerating billions of rotation words.

The seven corruption controls include five source/inventory/pin failures and
two internally rehashed mathematical faults. After source authentication passes,
Casimir 8 changed to 7 must fail at
`Arithmetic failed: rotations:selected_quadrupole_casimir`; time cross factor
2 changed to 1 must fail at
`Arithmetic failed: rotations:temporal_Leibniz_cross_factor`.

No giant exponential is evaluated in floating point; its logarithm is bounded
with exact integers and Fractions. The geometric/natural-action proofs and the
retained analytic background/linear existence premises remain mathematical
arguments, not conclusions of a finite test count. A passing report supplies
no nonlinear-family amplitude interval or proper-time actual-area remainder.
