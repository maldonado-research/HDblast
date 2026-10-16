# Actual full quadratic bulk forcing in the first relative-harmonic gauge

16 October 2026. The transformed first metric and scalar are evaluated
together on the retained actual background and actual solved gauge vector.
All mixed input jets are authenticated from the published first-gauge
capsule; the global collar is used throughout.

On ell=10^-6, 0<=x<=3, r>=0, x+r<=8, times all of S3:

| Absolute adapted-frame component bound | Strict cap |
| --- | ---: |
| Full natural Einstein Hessian QE | 9.444e20 |
| Registered reduced Einstein Hessian QEred | 9.444e20 |
| Full scalar Hessian Qphi | 3.855e15 |
| Quadratic gauge vector QC | 5971 |
| First covariant derivative of QC | 2.684e9 |

The Einstein bounds cover all 15 independent symmetric components and all
generated angular sectors. The norm is the maximum of absolute individual
components in the stated adapted background frame, including derivative
slots; it is not an energy or boost-invariant norm. These are second
parameter derivatives along the affine first-pair path. Second-equation
forcing has the opposite sign and the same absolute caps. The coefficient
of epsilon squared is one half the displayed Hessian.

Read `COMMON_GAUGE_BULK_BOUND_PROOF.md` for the field transformations,
global static reconstruction, ordered baseline derivatives, the complete
rank-three connection correction, and the universal contraction counts.
`COMMON_GAUGE_BULK_NORMAL.json` retains the exact rational bounds,
intermediate field caps, per-sector mixed jets and complete component
coverage. The optimized report has identical bytes.

After creating a new external output directory, run:

```sh
python -I -B check_common_gauge_bulk.py --output /tmp/new-common-gauge-bulk/normal.json
python -I -B -O check_common_gauge_bulk.py --output /tmp/new-common-gauge-bulk/optimized.json
python -I -B run_bulk_controls.py --output /tmp/new-common-gauge-bulk/controls.json
```

Only Python's standard library is required. Outputs cannot overwrite files
or be written inside this source/input tree. The explicit replay roster
and verification receipt record both modes and all reason-specific controls.

This is a continuous full *bulk forcing* certificate. It is not a complete
physical/gauge/constraint boundary-forcing certificate, a second evolution,
a nonlinear amplitude/remainder bound, or an observational claim.
