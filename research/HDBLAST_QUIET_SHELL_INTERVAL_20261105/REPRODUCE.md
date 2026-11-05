# Reproduce the finite arithmetic and source checks

Python 3.12 standard library only. Obtain the external manifest SHA256 from
the immutable commit's `.github/workflows/quiet-shell-interval-20261105.yml`
or its actual publication receipt. Assign that value to `PIN`; hashing an
untrusted local manifest without an independent comparison is not source
authentication.

From the repository root, choose a fresh existing output directory outside
the capsule and run:

```sh
OUT=$(mktemp -d)
python -I -S -B research/HDBLAST_QUIET_SHELL_INTERVAL_20261105/checks/check_quiet.py --manifest-sha256 "$PIN" --output "$OUT/normal.json"
python -O -I -S -B research/HDBLAST_QUIET_SHELL_INTERVAL_20261105/checks/check_quiet.py --manifest-sha256 "$PIN" --output "$OUT/optimized.json"
cmp "$OUT/normal.json" "$OUT/optimized.json"
python -I -S -B research/HDBLAST_QUIET_SHELL_INTERVAL_20261105/checks/run_controls.py --manifest-sha256 "$PIN" --output-dir "$OUT/controls"
```

The checker authenticates the closed capsule, then checks the registered
rational support distance, delay, active duration, majorant factors, clock
threshold and unit-sphere moment normalization. Reports must be outside the
capsule, and existing outputs are refused. Controls deliberately modify
fresh copies, including an arithmetic input with a rebuilt manifest, to
check rejection. The control fixture directory intentionally contains a
symlink test and must not be confused with the public capsule.

The normal and optimized reports must match exactly. These finite checks
do not run an Einstein–scalar solver or verify continuous propagation,
localized uniqueness, speed bounds, physical identification, coefficient
suprema, amplitude radius or external originality. Read the continuous
derivation and retained hypotheses separately. Primary source pins provide
provenance; they do not assert that a qualitative theorem's constants have
been evaluated. The continuous argument and implementation receive distinct
internal reviews, which are not external peer review.
