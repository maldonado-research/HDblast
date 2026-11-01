# Conditional proper-time shell volume for explicitly prepared mean data

For the retained Einstein–scalar/tension model, this checkpoint specifies compact
second-order mean initial data and bounds their actual coupled mean evolution.
At the physical shell and the fixed background-calibrated proper time
`tau_star=rho_b*ell*(5/4)`, with `ell=10^-6`, it gives:

| Second parameter derivative | Safe outward enclosure |
|---|---|
| Properly averaged shell scalar, `M2` | `[1.1296e30,1.9779e30]` |
| Relative shell volume, `v2=V2/V0` | `[-5.4520e16,-2.6521e16]` |
| Proper-time volume expansion, `Theta2` | `[-9.0208e21,-4.4458e21]` |

These are **changes relative to the background**, not negative total volume or
an assertion of contraction. The background expansion is positive. For each
admitted exact family, the coefficient signs imply, for sufficiently small
nonzero parameters of either sign, a positive volume below the background volume
and a positive expansion rate below the background rate at that proper time.
The usable amplitude interval is uncomputed and can depend on the family.
Large coefficients therefore do not establish a large physical effect.

All numbers use model units without an SI calibration. These results do not
establish radiation production, a force, an observational signal or a
higher-dimensional origin of the Big Bang.

## What changed

The earlier weighted-only construction did not numerically specify its second
monopole. Here the mean data are deliberately prepared using the registered
packet bump and explicit radial quadratures. The retained constraint inverse
then supplies the remaining exact family without changing that prescribed
second mean. This is a **different nonlinear completion** from the earlier
uninstantiated weighted-only family. Its preparation is a mathematical choice,
not a derived blast-generation mechanism.

Direct evolution in the original coordinates suffers from large cancellations.
A checked pure second-order gauge change removes the radial Lie-derivative
metric term and retains the corresponding scalar term. The conditioned second
jet describes the same prepared physical response. It does not claim that a
background-linear gauge transformation fixes the full nonlinear harmonic gauge.

The calculation keeps all five mean couplings, the complete actual bulk/shell
source history, and the proper-clock and volume-measure corrections. It derives
a new total-response bound from the prepared data rather than importing the
older total cap or initial error allocation. The scalar mean then drives the
intrinsic Gauss trace equation, avoiding a second large cancellation in a
direct metric-volume evaluation.

Read [the assembled statement](DERIVATION.md), the preparation and gauge proofs,
the initial quadratures, and the response/independent reconstruction components.
The full quantitative result is at one registered proper time. Bounds on input
histories do not constitute a signed volume history over a numerical interval.

## Reproduce

From a checkout containing the pinned dependencies, use Python 3.12.14:

```sh
python -m pip install -r research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/requirements.txt
python -I -S -B research/HDBLAST_PREPARED_SHELL_VOLUME_20261101/checks/run_checkpoint.py --repo-root . --output-dir /tmp/hdblast-prepared-volume-fresh
```

Choose a fresh output directory outside the checkout. Closed-source
authentication runs before scientific imports and afterward. Canonical reports
are reproduced in normal and optimized Python. Real rejected subprocesses,
in-process formula checks and successful relocated controls are counted
separately. Initial-background evaluations are distinguished from evolved
first-field oracle calls; neither count is an independent observation.

The result retains the analytic background, localized-constraint, evolution,
regularity and actual-first-field premises. Byte authentication and finite exact
checks do not re-prove those PDE theorems. Independently written rational
aggregation shares recorded input boxes and does not independently solve their
PDE. Internal AI-assisted review is not external peer review or formal
proof-assistant verification.

A working manuscript is supplied for author and specialist review, not as a
submission-ready or accepted paper. Earlier capsules and manuscripts remain
unchanged. No journal submission or Zenodo edition is made in this checkpoint.
[Next targets](NEXT_STEP.md) distinguish further mathematical certification from
the physical preparation, matter and observational questions still open.

Working paper: [PDF](manuscript/PREPARED_SHELL_VOLUME.pdf) · [TeX source](manuscript/PREPARED_SHELL_VOLUME.tex). The retained PDF is an exposition artifact; the checkpoint does not claim a canonical PDF rebuild.
