# Reproduce the finite arithmetic

Python 3.12 standard library only. Obtain the external manifest SHA256 from
the immutable commit's `.github/workflows/first-proper-geometry-20261106.yml`
or the actual publication receipt. Assign it to `PIN`. Merely hashing a
local manifest without an independent comparison does not authenticate it.

From the repository root, choose a fresh existing output directory outside
the capsule, then run the closed-source entry point:

```sh
OUT=$(mktemp -d)
python -I -S -B research/HDBLAST_FIRST_PROPER_GEOMETRY_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output "$OUT/normal.json"
python -O -I -S -B research/HDBLAST_FIRST_PROPER_GEOMETRY_20261106/checks/check_checkpoint.py --manifest-sha256 "$PIN" --output "$OUT/optimized.json"
cmp "$OUT/normal.json" "$OUT/optimized.json"
python -I -S -B research/HDBLAST_FIRST_PROPER_GEOMETRY_20261106/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir "$OUT/controls"
```

The normal and optimized reports must match. The evaluator checks the
pinned historical report identities and registered state, then computes
the future/past trace and first proper-geometry bounds using exact rational
arithmetic. The closed-source entry point also checks the complete capsule
inventory and prevents report output inside it. The evaluator module's
standalone CLI is a mathematical component, not this source-authentication
entry point.

Controls alter fresh copies to test rejection, including an authenticated
arithmetic mutant. Their fixture directory intentionally includes a symlink
test; it is not the public capsule. Existing outputs are refused.

These executions do not rerun the historical master-field campaign, solve
the nonlinear PDE, verify continuous gauge/geodesic arguments, supply
uniform nonlinear parameter bounds or establish external originality.
Continuous first-variation derivation and implementation review are
separate, and internal review is not external peer review.
