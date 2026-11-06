# Reproduce the finite quadratic arithmetic

Python 3.12 standard library only. Obtain the external manifest SHA256 from
the immutable commit's `.github/workflows/quadratic-proper-area-20261106.yml`
or the actual publication receipt. Assign it to `PIN`. Merely hashing a
local manifest without an independent comparison does not authenticate it.

From the repository root, choose a fresh existing output directory outside
the capsule, then run the closed-source entry point:

```sh
OUT=$(mktemp -d)
python -I -S -B research/HDBLAST_QUADRATIC_PROPER_AREA_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output "$OUT/normal.json"
python -O -I -S -B research/HDBLAST_QUADRATIC_PROPER_AREA_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output "$OUT/optimized.json"
cmp "$OUT/normal.json" "$OUT/optimized.json"
python -I -S -B research/HDBLAST_QUADRATIC_PROPER_AREA_20261106/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir "$OUT/controls"
```

The normal and optimized reports must match byte for byte. The evaluator
checks named finite-arithmetic statements about the known quadratic
determinant and shear/expansion-square sources for the proper cut at
amplitude zero. Historical inputs are authenticated and consumed; their
master-field campaign is not rerun. The evaluator module supplies the
mathematical component; use this closed-source entry point to authenticate it.

The entry point authenticates the complete closed file inventory, rejects
duplicate JSON keys and symlinks, and prevents output inside the capsule.
Existing outputs are refused. Controls alter fresh copies to test the
wrong manifest pin, changed source bytes, an unlisted file, a rehashed
quadratic-arithmetic error, a source symlink and an alternate lexical
script path with an inside-capsule output. Each child must reach its
intended rejection, create no report and exit unsuccessfully. Their
fixture directory intentionally contains a symlink test; it is not the
public capsule.

These executions authenticate source bytes and evaluate finite arithmetic.
They do not establish the continuous proper-cut or nonlinear PDE argument.
The unknown second Ricci variation and uniform bounds for full amplitudes
remain unresolved, so no full second-order area conclusion or nonzero
amplitude radius follows. Continuous derivation and implementation review
are separate; internal review is not external peer review.
