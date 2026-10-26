# Full angular response approximation and a nonzero tensor coefficient

26 October 2026. Conditional mathematics for one retained quadrupole state.

A continuous approximation now covers **all 23 retained angular amplitudes**
on the stated causal domain. Its complete positive-norm error, including
scaled time and radial first derivatives, is below **1.1787044291e28**,
meeting the declared **1e29** target. A separate calculation bounds both
equal transverse-traceless (TT) coefficients of the actual total formal
response strictly below zero at one exterior point.

These are formal second parameter coefficients in the retained
Einstein–scalar shell model. They do not provide an exact finite-amplitude
nonlinear spacetime or evidence for the proposed Big Bang origin.

## Results

| Result | Conditional certificate |
| --- | --- |
| Full retained angular space | 23 amplitudes; orthogonal blocks 5+8+8+1+1 |
| Selected symmetry | T+=T− exactly; 22 active amplitudes with both TT slots audited |
| Continuous full-space approximation error | X error <1.1787044291e28 |
| Complete integral-equation residual norm | <1.1790566331e28 |
| Error from the residual/resolvent route | <1.1795592728e28 |
| Actual total T+=T− at (x,r)=(1/4,5/4) | [−0.001311761536, −0.000001509027] |
| Boundary structure | Five Dirichlet and eighteen Robin modal rows; no boundary row discarded |

The full-space norm uses the probability-normalized positive angular
Gram/Frobenius inner product. It is a comparison norm, not physical energy.
The error target is an absolute coefficient accuracy; it is not a relative
pointwise accuracy, a finite-amplitude error or an observational precision.
All displayed inequalities are rounded outward from exact rational data.

The global candidate has zero TT entries and its global error is too broad
to resolve the small actual TT coefficient. The separate source-resolved
TT calculation supplies that nonzero enclosure. The global approximation
also does not replace the sharper frozen mean and anchored-volume sign
certificates from 25 October.

## Why the tensor result matters within this model

Both TT initial displacements and velocities vanish. Their linear equations
decouple from every other amplitude, while their complete quadratic bulk
sources are equal. The evaluated backward cone does not reach the shell.
Thus the nonzero total TT response there is driven by the bulk quadratic
source; it is not a nonzero initial tensor inserted by hand or a shell term
omitted from the calculation. A second coefficient with both TT entries
identically zero cannot represent this registered solution.

The coefficients use the registered gauge and background angular basis.
A physical gravitational-wave observable, detector strain, spectrum or
cosmological signal has not been computed. Scalar-induced tensor effects
and their gauge issues are established subjects; no new general wave
mechanism, priority claim or external novelty assessment is asserted.

## Construction and proof

The candidate is `V=B0+B1`. B0 uses the exact expression-defined registered
rational displacement and one continuous P1 velocity. B1 retains every
constant-matrix bulk and Robin coupling, with independent P1 displacement
and radial-derivative approximants used only inside its forcing. All 23
components are reconstructed in the same fixed angular basis.

The fixed denominator-256 mesh has 190 shared candidate nodes and 192
actual-minus-nominal error cells. The rational cutoff is C1 after all
registered derivatives; higher interpolation bounds are one-sided inside
aligned analytic cells. The full retained source, actual input errors,
velocity reconstruction, coefficient defects and complete remaining
feedback tail are charged. The error bound holds on the whole cone, not
only at the six illustrative sample evaluations.

Read [DERIVATION.md](DERIVATION.md), then the detailed certificates:

- `initial/INITIAL_DERIVATION.txt` and `initial/results/INITIAL_REPORT.json`.
- `compute/response-normal/RESPONSE_REPORT.json` and the three-field sample
  table `EVALUATED_CANDIDATE.json` in the same directory.
- `structure/ANGULAR_STRUCTURE_PROOF.md` and the exact 22↔23 maps.
- `shell/FULL23_SHELL_AND_SUBSIDIARY.md` for all boundary rows, the rank-eight
  shell-source packing and the retained subsidiary premises.
- `tensor/TOTAL_TT_CERTIFICATE.txt` and `TT_TOTAL_RESPONSE_NORMAL.json`.
- `structure/response-review/INTEGRAL_RESIDUAL_LEMMA.md` and
  `ANGULAR_DERIVATIVE_LEMMA.md` for the residual and intrinsic-sphere derivative
  corollaries. The latter does not replace the physical warp/frame conversion.

The finite candidate generally has a nonzero Robin defect. It is bounded,
not silently set to zero. The integral residual is distinct from a
pointwise second-derivative Einstein residual. The inherited exact formal
construction and subsidiary theorem remain conditional premises; this
checkpoint supplies a numerical approximation and error certificate, not
a new existence theorem for nonlinear families.

## Replay

Python 3.12 and its standard library suffice. Supply a checkout containing
this capsule and the eight frozen upstream capsules identified in the
replay program. No network access is needed during scientific replay.
Keep the new manifest SHA256 outside the capsule; its publication pin is
in the dedicated GitHub workflow.

```
python -I -S -B reproduce.py --manifest-sha256 EXTERNAL_PIN \
  --repo-root /path/to/HDblast --output-dir /tmp/fresh-angular-replay
```

The default replay authenticates the closed source, reruns the bounded
new calculations in normal and optimized Python, compares reference bytes,
and runs independent algebra, integration, interval and corruption controls.
Use fresh external output paths. Archive verification additionally accepts
an archive path and its externally retained SHA256. No older production
campaign is rerun as a hidden dependency.

Internal AI-assisted reviews are recorded separately from primary
calculations. They are not external specialist peer review or
proof-assistant formalization. The focused primary-source literature note
explains the gauge and interpretation limits; it is not exhaustive.

## Remaining work

The next numerical target is sufficient higher-jet control for a pointwise
second-derivative residual or a well-defined tensor observable. Finite
angular closure supplies angular derivative control, but does not supply
missing second time/radial derivative error bounds. A finite-amplitude
nonlinear realization, an observational prediction and a physical Big Bang
origin mechanism remain open. Existing manuscripts and research are
preserved. No new Zenodo version or journal submission accompanies this
checkpoint. Private editorial history is excluded.
