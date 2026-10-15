# Exact algebra reproducibility supplement

This supplement replays four exact symbolic checkers. It verifies finite algebra and the field, shell, radial, and support identities for a specified exact background. It does not establish a general PDE existence, range, stability, or shell/interface sufficiency theorem, certify the written proof in arbitrary dimension, or establish novelty.

Use Python **3.12.14**, SymPy **1.14.0**, and mpmath **1.3.0**. The runner checks these versions. Dependencies must already be installed (the versions are listed in `requirements.txt`); replay uses no network or Git.

From any working directory, run:

```sh
python /path/to/capsule/supplement/reproduce.py --output-dir /path/to/new-replay-output
```

The output directory must not exist and must be outside `supplement/`. The runner refuses to overwrite previous results. It writes full reports, stdout/stderr logs, and `REPLAY_RECEIPT.json`. Exit **0** and `PASS_EXACT_AGREEMENT` require all expected agreements and unchanged package hashes. A failed or interrupted run is not a successful receipt.

The runner starts eight positive subprocesses: each checker under ordinary Python and `python -O`. It compares each complete JSON report with its packaged expected report. It also compares the complete reports between modes, removing only the finite-matrix field `python_optimization`, whose actual values must be 0 and 1. No other field is excluded. The positive inventories are:

| Checker | Positive checks per mode | Embedded wrong-formula controls per mode |
| --- | ---: | ---: |
| Finite commuting geometry | 132 | 12 |
| Rotational shell algebra | 46 | 9 |
| Quadrupole admissibility algebra | 29 | 4 |
| Explicit background identities | 39 | 10 |
| Total | 246 | 35 |

The matrix program's **144 generated records already include its 12 controls**; those controls are not 12 additional checks. Its examples cover the five joint-multiplicity strata in dimension four, exact ranks and inverse formulas, velocity obstructions, and constructed commuting curves. The rotational and quadrupole programs check the exact expressions and angular integrals in their sources; together with the matrix program they retain 207 positive checks and 25 controls per mode. The background program independently computes the connection and Ricci tensor from the stated metric, then checks the scalar equation, shell normals and junction values, radial coefficients, and exact support bounds. These numbers are the inventories of actual captured runs, not a historical count inherited from an earlier note.

Six additional subprocesses exercise `--inject-failure` in both Python modes on the two earlier checkers and the background checker. They are separate from the 35 embedded controls per mode:

| Injection | Required result in each mode |
| --- | --- |
| Rotational shell checker | Exit 1; failure JSON with error `explicit_injected_failure: nonzero residual 1` |
| Quadrupole checker | Exit 1; no result JSON; final exception `RuntimeError: Exact identity failed: rotation_12_polynomial_covariance` |
| Explicit background checker | Exit 1; failure JSON with error `RuntimeError: injected background mismatch detected` |

The failure contracts compare the full failure JSON where present, exact exit code, result presence/absence, stdout, and the exact failure reason. The complete quadrupole traceback is preserved in stderr logs; its location-dependent traceback paths are not an equality criterion. The background injection mutates the Ricci radial coefficient from 24 to 23 and fails the metric-derived Ricci comparison. The finite-matrix checker has no injected-failure option, so no such process is claimed for it.

`INPUT_PINS.json` records SHA-256 and byte lengths for the five copied source files and nine expected-report files. Before execution the runner verifies those pins. It also hashes every file in the supplement before and after replay. The final capsule's top-level manifest supplies the independent package inventory.

The two earlier programs and the rotational companion note were copied byte-for-byte from commit `a121e78c7eac5c6c060cabdda4e9b5a53d011ee0`; the exact original paths are in `INPUT_PINS.json`. The new matrix and background checkers were copied byte-for-byte from `math/sol/finite_commuting_geometry.py` and `math/background/check_explicit_background.py` in the current 15 October 2026 agent draft and are not attributed to that earlier commit. The expected matrix reports agree in full with that draft's original ordinary/optimized receipts. The historical companion note is included because the rotational checker hashes it; the checks do not certify every statement in that note.

The explicit background uses its displayed potential and tension and demonstrates a specified example in the paper's background class. It is a separate model from the registered HDBLAST branch. The checks do not certify that branch or supply nonlinear continuation through the obstructed tangent or shell contact.

Expected reports are captured reference outputs, not independent scientific proofs. Fresh replay and exact report agreement demonstrate reproducibility of the implemented finite checks. The mathematical argument and its limitations must be assessed in the accompanying manuscript.
