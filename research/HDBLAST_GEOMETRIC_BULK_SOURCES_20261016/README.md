# Quadratic geometric sources and an intrinsic shell diagnostic

`HDBLAST_GEOMETRIC_BULK_SOURCES_20261016` continues the
[numerical-input checkpoint](https://github.com/maldonado-research/HDblast/tree/main/research/HDBLAST_NUMERICAL_RESPONSE_INPUTS_20261016).
It concerns one registered conditional Einstein–scalar material-shell model,
its actual first harmonic field and a fixed second initial prescription.
Earlier frozen research and manuscripts are unchanged.

## Results and scope

- **Complete geometric source:** an indexed Ricci Hessian, natural and reduced
  Einstein source, scalar source, gauge covector and unit-principal forcing.
  The exact 23-amplitude maps retain both transverse-traceless orientations.
  All-point polynomial reconstruction is required after coefficient extraction;
  finite angular witnesses alone do not establish the identities.
- **Actual source enclosures:** authenticated first-field C3 derivative errors
  and background jets produce ordinary frame metric/scalar two-jet intervals.
  Eighty-one clipped cells cover the causal triangle, with a 64-cell evolved
  radial profile and six shell time cells. The bulk forcing norm satisfies
  `F_H <= 17387513156573`, improving the retained upper bound by a factor
  greater than 9.8. The full second-response residual is not updated by this bound.
- **Intrinsic shell diagnostic:** for the jointly pulled-back material scalar,
  `I=gamma^{-1}(d phi_shell,d phi_shell)` has vanishing background and first
  variations. Its second variation is invariant under the specified joint
  field/embedding gauge transformations. Its mean on the chosen background
  time sections has exactly two simple zeros in each scaled-time bracket
  `[5/4,11/8]` and `[13/8,7/4]`. Four distinct zeros are certified; additional
  zeros elsewhere are not excluded. This is a formal quadratic coefficient,
  not a finite-amplitude sign, physical force, or new physical clock.
- **An evolved second-mean corollary:** the retained global derivative bound
  proves opposite scalar-mean signs, with magnitude greater than `1.85e30`,
  on two whole radius cells for scaled times `[1/1000,1/200]`. At every such
  time there is a radial mean zero between the cells. This is quantitative
  persistence of initial signs, not an independent waveform computation or
  gauge-invariant second-field observable.

The universal curvature and gauge-variation formulas are established results.
The contribution here is the model-specific exact projection, actual interval
source evaluation, intrinsic diagnostic and checked corollaries. No external
novelty or priority assessment is claimed. Internal AI-assisted mathematical
and computational reviews are not external specialist peer review.

## Files

- `manuscript/GEOMETRIC_BULK_SOURCES_NOTE.pdf` and `.tex`, with the separate
  table source and checked build records, give the technical working note.
- `geometry/` implements the full indexed/frame source calculation and records
  the affine-path and general/off-shell versus actual harmonic conventions.
- `projection/` contains the coefficient-generic angular algebra, universal
  maps, exact reconstruction checks and public coefficient evaluator.
- `numeric/` contains actual two-jet enclosures, the exterior-only refinement,
  equivalent outward integer/Fraction arithmetic, bulk/shell certificates,
  complete cell cover and authenticated premise snapshots.
- `shell-invariant/` proves the intrinsic scalar formulas and continuous
  monotonicity certificate for the two pairs of simple mean zeros.
- `response-signs/` contains the evolved-mean persistence proof and certificate.
- `review/`, `numeric-review/` and `manuscript-review/` record independent
  internal checks, exact accepted hashes and clearly stated coverage limits.
- `literature/` compares verified primary sources and states the limits of
  transferring existing curvature, shell and numerical results.
- `GAUGE_COVECTOR_BOOKKEEPING.md` distinguishes the varying-metric lowered
  gauge Hessian from background lowering of its upper-vector Hessian.

Only the explicit public file selection is included. Private editorial history,
raw discovery logs and downloaded copyrighted primary bodies are excluded.

## Replay

The scientific replay needs Python's standard library and no network access.
All runtime inputs are local and pinned. Obtain the manifest SHA256 from the
separate trusted publication/freeze record, then run from the extracted root:

```sh
python -I -S -B reproduce.py \
  --manifest-sha256 "$verified_manifest_sha256" \
  --output-dir /tmp/hdblast-geometric-replay
```

The output directory must not exist and must be outside the immutable source.
The roster runs selected scientific classes in normal and optimized modes,
compares complete typed JSON and raw bytes, checks specific rejection reasons,
and requires rejection processes to leave no output. It distinguishes direct,
nested and embedded in-process controls; counts are processes, not theorems.
Historical reports that are authenticated but not replayed are identified.
The final archive was also checked after extraction to a fresh temporary path.

Rebuilding the optional PDF requires a LaTeX installation. From `manuscript/`,
run `pdflatex` twice with an external output directory; both `.tex` files are
needed. The scientific replay does not require LaTeX.

## Remaining work

The full numerical second waveform and a continuous response residual
certificate remain unfinished. Next is a continuous selected-space candidate
with certified source/input, causal quadrature, initial, shell, reflection and
reconstruction errors. Available first-field C3 data do not automatically
certify a derivative of the quadratic source. An exact finite-amplitude
nonlinear family, invariant force, observational fit or cosmological origin
is not established here. External novelty assessment and peer review remain
outstanding.
No new Zenodo version or journal submission accompanies this checkpoint.
