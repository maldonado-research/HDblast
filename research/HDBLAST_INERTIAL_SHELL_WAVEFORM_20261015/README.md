# HDBLAST inertial shell waveform — 15 October 2026

The newly specified width ell=10^-6 has a certified **conditional linear
first-variation** waveform. Retaining shell inertia resolves its actual
post-pulse trace continuously:

- 3.2e-8<H(x)<5.3e-8 for every 2<=x<=3. The physical quadrupole coefficient
  has negative sign because its background gain is negative.
- At least one zero lies in each of (1.49998,1.49999) and (1.96,1.97).
  Uniqueness and a complete zero count are not asserted.
- A fresh full-domain integral gives 6.06e11<C_(3ell)<7.44e11.

Read [the main manuscript](INERTIAL_SHELL_WAVEFORM.md). The comparison wave
solves a free bulk problem with the exact constant boundary inertia;
analytic potential and boundary-coefficient error estimates enclose the
actual coupled linear solution. The new coefficient is not transferred
from the previous ell=10^-10 state.

## Scientific files

| File or directory | Role |
| --- | --- |
| [computed/INERTIAL_WAVEFORM_ENCLOSURE.json](computed/INERTIAL_WAVEFORM_ENCLOSURE.json) | Exact endpoints, strict tail and sign gates, fresh response integral |
| [computed/INERTIAL_WAVEFORM_CELLS.csv](computed/INERTIAL_WAVEFORM_CELLS.csv) | 4353 continuous-cell comparison, actual and physical enclosures |
| [compute_inertial_waveform.py](compute_inertial_waveform.py) | Production certificate using outward intervals and rational sums |
| [inertial_error_bounds.py](inertial_error_bounds.py) | Global, post-pulse and cell-specific residual bounds |
| [sol/FREE_INERTIAL_COMPARATOR.md](sol/FREE_INERTIAL_COMPARATOR.md) | Exact boundary kernel and its validated evaluator |
| [referee/INDEPENDENT_INERTIAL_BOUND.md](referee/INDEPENDENT_INERTIAL_BOUND.md) | Independent characteristic cancellation and complete memory error proof |
| [astra/INERTIAL_COMPARATOR_AND_POSITIVE_TAIL.md](astra/INERTIAL_COMPARATOR_AND_POSITIVE_TAIL.md) | Separate energy-kernel proof and conservative positive-tail certificate |
| [archive/BACKGROUND_DERIVATIVE_AUDIT.md](archive/BACKGROUND_DERIVATIVE_AUDIT.md) | Wider collar and independent coefficient derivative bounds |
| [millennium/HDBLAST_METHODS_TRANSFER_PUBLIC_SAFE.md](millennium/HDBLAST_METHODS_TRANSFER_PUBLIC_SAFE.md) | Independently derived scalar residual-energy method and its limits |
| [literature/INERTIAL_WAVEFORM_PRIMARY_LITERATURE.md](literature/INERTIAL_WAVEFORM_PRIMARY_LITERATURE.md) | Bounded primary-source review with applicability limits |
| [inputs/PRIOR_SOURCE_PINS.json](inputs/PRIOR_SOURCE_PINS.json) | Authenticated retained mathematical inputs |

The CSV time endpoints are integers divided by 4096. Field endpoints use
the value denominator in the JSON report (2^112). Each interval covers its
whole closed time cell. Decimal summaries and graphics are displays;
the exact integer/rational endpoints decide the certificate.

## Reproduction

The production waveform and rational estimators use Python's standard
library. The complete replay also requires the pinned SymPy dependency in
`requirements.txt` for the independent formal coefficient identities. The
replay entrypoint uses an externally supplied source-manifest digest. Replace
`SOURCE_MANIFEST_SHA256` with the digest in the published
`SOURCE_FREEZE_RECEIPT.json` and use a fresh output directory. That external
receipt accompanies the capsule in the dated continuation directory.

```sh
python -I -B reproduce.py --manifest-sha256 SOURCE_MANIFEST_SHA256 --output-dir /tmp/hdblast-inertial-fresh-replay
```

Individual mathematical components can also be replayed from the capsule:

```sh
python -I -B compute_inertial_waveform.py --output-dir /tmp/hdblast-inertial-fresh-normal
python -I -O -B compute_inertial_waveform.py --output-dir /tmp/hdblast-inertial-fresh-optimized
python -I -B astra/check_inertial_error.py --source-root .
python -I -O -B astra/check_inertial_error.py --source-root .
python -I -B sol/check_comparator.py --source-root .
```

Standalone notes describe their independent checker interfaces and fault
controls. The frozen `REPLAY_ROSTER.json` defines exactly which components
and controls the complete replay executes. Publication and byte-verification
receipts accompany the capsule in the dated continuation directory; scientific
source files alone do not establish that a remote publication has occurred.

The retained analytic background and PDE reconstruction are inputs. An
exact-family interpretation still requires the localized inverse,
reflection, gauge, whole-correction support and finite-regularity evolution
hypotheses for this wider tangent. There is no effective nonlinear epsilon,
finite-amplitude waveform, observed time scale, cosmological fit or Big Bang
origin claim. Internal independent review and exact finite arithmetic are
distinct from external peer review and formal proof verification.
