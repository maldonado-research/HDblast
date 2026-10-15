# HDBLAST: a conditional quantitative shell-response enclosure

15 October 2026 UTC. GitHub research checkpoint.

For the retained near-phi=+1 regular background at exact delta=1/1000, and a
newly specified aligned quadrupole with ell=10^-10 in model coordinates,

\[
\mathcal R_\varepsilon(3\rho_b\ell)
=\varepsilon^2\mathcal C_{3\ell}+o(\varepsilon^2),\qquad
10^{15}<\mathcal C_{3\ell}<4\times10^{17}.
\]

This is a coarse analytic interval in a specified normalization. It follows
from the conditional exact-family construction, an explicit free-wave
comparison, an energy-plus-characteristic derivative-trace bound, and exact
rational arithmetic. It is not a computed waveform, an effective nonlinear
amplitude interval, an observational scale or evidence for a cosmological
origin. A narrow unit-normalized bump can have a large derivative coefficient;
the useful nonlinear parameter range remains unquantified.

Read [the primary proof](SHELL_RESPONSE_ERROR_CERTIFICATE.md) for the state,
physical proper-time definition, background, support and evolution assumptions.
The whole nonlinear correction requires a separate strict causal margin.
The result is internally reviewed conditional mathematics; no external peer
review, certified novelty or formal PDE verification is claimed.

The coefficient simplifies exactly to

\[
\mathcal C_T=\frac{3\pi^2}{\rho_b^2p_b^2}
\|w_z-(C_b/2)w\|_{L^2(0,T;z=0)}^2.
\]

The comparison controls this physical derivative observation, including
positive shell inertia. It does not infer it from displacement energy alone.

## Reproduce

Use Python 3.12 and the pinned symbolic dependency in `requirements.txt`.
Obtain the external manifest SHA256 from the dated publication receipt or
the repository workflow. Keep the extracted capsule's directory name.

```sh
python -I -S -B verify_manifest.py --expected EXTERNAL_MANIFEST_SHA256
python -I -B -m pip install --only-binary=:all: -r requirements.txt
python -I -B reproduce.py --manifest-sha256 EXTERNAL_MANIFEST_SHA256 --output-dir /tmp/hdblast-response-fresh
```

Choose a new output directory outside the capsule. The runner authenticates
the complete manifest before executing local scientific code, then copies
only declared dependencies into separate workspaces for normal and optimized
processes. Positive JSON must match its golden data. Wrong-formula controls
must fail for their named reasons. Logs, actual process counts and finite
identity counts are recorded in `REPRODUCTION_RECEIPT.json`; none are counts
of analytic theorems proved by a computer.

`certify_response.py` implements the exact-rational arithmetic consequence
of declared coefficient/residual bounds. It is not a waveform solver and does
not establish those analytic bounds by accepting numeric inputs. The included
actual-background bound is derived separately, using a pinned analytic branch
theorem and a strict first-exit collar argument. Square roots, exponential
growth and pi are enclosed with rational inequalities; floats are rejected.

## Independent work and source scope

- [Actual background and rigorous collar](archive/ACTUAL_BACKGROUND_AND_RIGOROUS_COLLAR_BOUNDS.md): same model and branch, with exact source pins. Historical floating profiles and alternative-branch data are reference material only; their archived programs are not run.
- [Astra derivation](astra/RESIDUAL_TO_PHYSICAL_RESPONSE.md) and [acceptance](astra/ASTRA_QUANTITATIVE_RESPONSE_ACCEPTANCE.md).
- [Sol characteristic derivation](sol/RESIDUAL_TO_PHYSICAL_RESPONSE.md) and [independent quantitative acceptance](sol/crossreview/SOL_QUANTITATIVE_ACCEPTANCE.md).
- [Critical review](referee/FINAL_QUANTITATIVE_RESPONSE_ACCEPTANCE.md) and its independent rational coarse calculation.
- [Focused primary-literature comparison](literature/VALIDATED_RESPONSE_PRIMARY_LITERATURE.md): bounded search, not an exhaustive literature or novelty review. Raw retrieved paper bodies are not distributed.

The two independent coarse arithmetic checks do not import the primary
certificate implementation. Both reproduce the stated interval by simpler
rational bounds. The Sol positive-majorant evaluator is an additional generic
tool; its manufactured checks are not HDBLAST waveform calculations.

Review receipts can mention private replay workspaces: those executable
scratch copies are not part of this capsule. Their result receipts and all
required source dependencies are included and covered by the manifest.

The next scientific target is a resolved waveform on a useful chosen scale,
with sharper validated response bounds, followed by an effective nonlinear
remainder estimate. No new Zenodo version or APS submission accompanies this
checkpoint.
