# Executed full saved-trajectory audit

Status: **executed successfully** on GitHub Actions with Python 3.12 and NumPy 2.2.6.

Run [36931967040](https://github.com/maldonado-research/HDblast-archive/actions/runs/36931967040),
source commit 5d37bbf18ddc24b5218dd89d9e22be7a26c7c567, job 110603219023.
The fetched complete report has nine requested runs and no errors. Its script SHA-256
9e37f073f0027b98994f1600903cf3a294dd07c96f1940443da5bf49dc3c6605
matches the saved source independently hashed by root.

This reads existing NPZ arrays with pickle disabled. It does not rerun the five-dimensional equations.

The saved producer labels in the table are not the final audited reliability verdicts: Y=3 remains conditional and Y=5 remains inconclusive/unreliable. The reliable Y=2 result leads this checkpoint.

## New sampled result

The new retrospective gate is F_abs>=0.9, |Weyl|/linear_radiation<=0.03, H>0,
positive linear radiation, finite diagnostics and the original near-shell constraint prefix.
F_abs counts absolute vacuum/scalar/source contributions instead of cancelling them.

| Run | Longest full-gate DeltaN | Passing samples | Legacy point-gate DeltaN | Saved producer label (not audited verdict) |
|---|---:|---:|---:|---|
| fine_Y2_dc1e-2_dzf2.5e-4 | 0.314463407 | 67 | 0.926199364 | PASS-combined |
| fine_Y2_dc1e-4_dzf2.5e-4 | 0.313042599 | 67 | 0.838236008 | PASS-combined |
| fine_Y3_dc1e-2_dzf1.25e-4 | 0.349426998 | 63 | 0.936128315 | PASS-combined |
| fine_Y5_dc1e-2_dzf1.25e-4_cfl0.25 | 0.326283033 | 31 | 0.832621888 | PASS-combined |
| main_Y0.5_dc1e-2_dzf5e-4 | 0.000000000 | 0 | 0.000000000 | FAIL-Weyl |
| main_Y0.7_dc1e-2_dzf5e-4 | 0.000000000 | 0 | 0.000000000 | FAIL-Weyl |
| main_Y1.5_dc1e-2_dzf5e-4 | 0.000000000 | 0 | 0.000000000 | PASS-conservative |
| main_Y2_dc1e-2_dzf5e-4 | 0.314463817 | 67 | 0.851679195 | PASS-combined |
| main_Y2_dc1e-4_dzf5e-4 | 0.313043023 | 67 | 0.838341739 | PASS-combined |

For Y=2, dc=0.01, the finer history has 67 passing samples spanning **0.31446340736608747 e-fold**,
from H0 tau=7.397523901922649 to 11.220520380640787. The coarser history gives
0.3144638170146514, a difference about 4.1e-7. For dc=0.0001 the finer duration is
0.31304259944474744 and the coarser duration 0.31304302304987.

This earlier interval is before the later comparison plateau. At the saved dc=0.01 plateau,
F_abs=0.7654491546246199 and the scalar-free upper bound F0=0.7654491662636436;
that point fails the new criterion. The earlier maximum F_abs=0.9416382931240115
alone would not establish the simultaneous Weyl gate; the report tests them together.

The roughly 0.3195 e-fold **from plateau to turnaround** is a separate measurement and
must not be described as accepted radiation domination or added to the accepted interval.
The legacy point gate in this report also differs from B3's complete registered plateau rule.
Original B3 labels are preserved.

## Numerical limits and ledgers

These are contiguous **sampled** intervals, not rigorous lower bounds or proofs that the gate
holds between records. Finer continuations inherit parent histories before restart; the accepted
intervals occur after restart, but agreement is not independent end-to-end convergence.
Y3 remains conditional and Y5 unreliable under the old audit even where pointwise gates pass.
The illustrative half-e-fold flag in JSON is a retrospective comparison, not a new cosmological
requirement or registered discovery test.

Across the audited reliable records the reconstructed Friedmann relative closure is at most
4.15057625e-13.
That is largely an algebraic consistency check because Weyl was reconstructed from the identity;
it is not independent validation of the field equations.
For the finer Y2 histories the trapezoidal weighted radiation-ledger endpoint residual is about
6.1e-6. The report separately records differential residuals, intermittent constraint checks,
source file hashes and restart ancestry. Those diagnostics are numerical, not certificates.

This demonstrates a brief sampled conventional-radiation interval **in the phenomenological
friction model**. It establishes neither physical quantum production, thermalization, a
realistic temperature history, a sustained hot Big Bang, nor observational support.

The new simultaneous point gate does not implement the original B3 stationarity and half-e-fold plateau rule. Its measured duration is a different, retrospective sampled diagnostic. An illustrative half-e-fold comparison in the report is not a new preregistered or observational requirement.

## Public replay inputs

The curated package now supplies the eleven NPZ timeseries and matching summaries needed for these nine histories, including restart parents. Each binary and summary was SHA-256 matched to this executed report's provenance before copying. Run the unchanged audit with --scan-root data/B3_Y_scan; only location and source-ref metadata differ. See data/PROVENANCE.md. Numeric histories are preserved, not freshly evolved.
