# Conditional shell and anchored-volume sign bounds

25 October 2026. A model-specific mathematical checkpoint.

A new three-term coupled-mean approximation, together with a complete
positive component error bound, places the conditional exact shell
coefficients `N0` and `J0` strictly below zero at scaled times `x=9/8` and
`x=5/4`. At `x=9/8`, bounding the full proper-time clock and quadratic
first-order terms also places the **anchored volume second parameter
derivative** `v2=V2/V0` below zero. This is a conditional kinematical jet
result; no full angular second-order Einstein solution is constructed.

The final exact rational intervals are in
[results/JOINED_NORMAL.json](results/JOINED_NORMAL.json). Rounded outward:

| Quantity | Scaled shell time | Conditional enclosure |
| --- | --- | --- |
| Total coordinate coefficient J0 | 9/8 | [−1.14562550e16, −1.00037137e16] |
| Total coordinate coefficient J0 | 5/4 | [−2.81184636e17, −2.68928911e17] |
| Anchored volume coefficient V2/V0 | 9/8 | [−1.71843851e16, −1.50055681e16] |

These bounds include candidate quadrature/arithmetic uncertainty and the
complete response error. The simpler generic ten-state comparison also
proves the two J0 signs and the volume sign, with wider intervals; the
conclusions do not depend on the sharper integrated bound alone.

## What changed

The old positive J0 candidate sample was a difference of terms about
`1.75e21`; its published global error bound did not certify a sign.
The exact registered displacement obeys `u0_w=C u0_J,r`. Preserving that
correlation by integration by parts removes a large boundary interpolation
defect. A boundary-only repair flips the candidate sign but still has an
inconclusive total-response interval. It is retained as a sensitivity
experiment, not used as the final sign proof.

The final candidate is `B0+S B0+S^2 B0`. Both feedback terms use the exact
registered rational displacement and its derivative; B0 retains the
previous continuous P1 velocity. S is the full frozen constant-matrix
bulk/Robin operator. One-dimensional moment identities evaluate the
second feedback at the shell while retaining all interior metric/scalar
couplings. A fixed 256-subinterval quadrature evaluates the moments.

The error comparison uses the actual variable-coefficient operator T,
actual initial uncertainty, omitted bulk/shell source, all coefficient
defects, and the **complete third-and-higher feedback tail**. Frozen P1
displacement/derivative envelopes bound norms only; they are absent from
the new candidate. Earlier source capsules and numerical claims are
preserved unchanged.

## Read the proofs

- [Finite component comparison and complete error identity](numeric/finiteproof.md).
- [Second-feedback moment identities and fixed quadrature](balanced/b2/B2_DERIVATION.md).
- [Boundary cancellation experiment](balanced/EXACT_BOUNDARY_PROOF.md).
- [Anchored observable bound, clock and normalization](analysis/OBSERVABLE_BOUND.md).
- [Finite proper-time tube and kinematical completion](analysis/FINITE_TUBE_OBSERVABLE.md).
- [Focused primary-source comparison and its search limits](analysis/LITERATURE_SCOPE.md).

All earlier regular-branch, forcing, projection, subsidiary and qualitative
coefficient-existence premises remain conditional. The volume coefficient
assumes smooth exact mean histories and the known smooth first-order
induced metric. A smooth kinematical completion supplies its second jet;
it need not satisfy the missing angular field equations. The fixed
physical initial anchor and intrinsic normal proper-time clock are part
of the definition. Coordinate-time volume is a different quantity.

## Replay

Use Python 3.12 and a repository checkout containing this capsule and the
seven frozen dependencies named in `numeric/finiteproof.md`. No package
installation is required. Retain the manifest SHA256 outside the capsule;
the dedicated repository workflow provides its publication pin.

```
python -I -S -B reproduce.py --manifest-sha256 EXTERNAL_PIN \
  --repo-root /path/to/HDblast --output-dir /tmp/fresh-shell-precision
```

The default replay authenticates the entire source closure before any
project execution, rejects seven altered-source controls, reruns the
bounded component/feedback/boundary calculations in normal and optimized
Python, and compares every generated result with the frozen reference
bytes. Independent polynomial, noncommuting-matrix and interval-assembly
checkers are also run. `--certificates-only` skips the new numerical
evaluations while retaining source authentication and certificate controls.
When an archive is supplied, add `--archive PATH --archive-sha256 PIN` to
verify its complete member list, CRC and byte equality with the source.
Outputs must be fresh and outside the source capsule.

The independent reviews are internal AI-assisted mathematical/software
reviews. Their source and result pins are retained. They are distinct
from external specialist peer review and proof-assistant formalization.
Here, a certified bound means an explicit conditional inequality with
outward arithmetic and a reproducible error ledger; it does not mean
approval by a journal or scientific institution.

## Limits and next scientific task

The negative signs concern two coordinate samples and one anchored
parameter second derivative. No whole-history sign, finite-amplitude
volume change, sign of the expansion-rate correction, or physical force
is established. The numerical magnitudes depend on the registered
perturbation normalization and do not measure the size of our universe.

The remaining scientific target is the compatible nonmean angular
second-order construction with its complete shell conditions, followed
by a full constraint/residual check. An exact finite-amplitude nonlinear
family and a distinguishable observational prediction remain separate
requirements. The proposed higher-dimensional origin of the Big Bang
has not been established. No new Zenodo version or journal submission
accompanies this checkpoint.
