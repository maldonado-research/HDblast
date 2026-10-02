# HDBLAST — Could a higher-dimensional "blast" have sparked our Big Bang?

[Readable project overview](https://maldonado-research.github.io/projects/hdblast/) · [All research projects](https://maldonado-research.github.io/)


**Ricardo Maldonado · independent researcher · research program 2025–2026**

**Website:** <https://maldonado-research.github.io/HDblast/> ·
**Zenodo (all versions):** [10.5281/zenodo.17088132](https://doi.org/10.5281/zenodo.17088132) ·
**Static-branch companion:** [Zenodo 22922928](https://zenodo.org/records/22922928), version 1.0 (23 Sept 2026),
separate concept [10.5281/zenodo.22922927](https://doi.org/10.5281/zenodo.22922927).
Version and concept verified from the author-provided record screenshot; file list not independently checked.

> HDBLAST asks whether a violent gravitational event in a fifth dimension, the "blast", could have
> transferred energy into our four-dimensional universe and started the hot Big Bang.

This is a **speculative, testable research hypothesis. It is not a discovery.** It has not
proven that extra dimensions exist or what caused the Big Bang, and no observation supports it
yet. It has not been peer reviewed. What the program does have is a dated public record:
explicit equations, reproducible code, certified and numerical calculations, negative results,
and published corrections. Every result below carries a label saying exactly how well it is
established.

---

## The idea in one picture

Our universe may behave like a thin **shell (a "brane")** inside a larger five-dimensional
space. Gravity and a field called a *scalar* fill that larger space. If the shell sits in an
unstable balance, like a ball on a hilltop, it can roll off, and that motion releases energy.
HDBLAST asks whether such a higher-dimensional event could have produced the hot, dense early
universe that the standard Big Bang model starts from. The hypothesis would **extend** the
standard Big Bang story, not replace it.

![The registered shell sits on a hilltop; its 5D spectrum matches a closed-form effective theory](hdblast/figures/shell_hilltop_and_spectrum.png)

*Left: the registered shell sits on top of a potential "hill", so it is unstable. Right: the
instability rate from full five-dimensional theory agrees with a closed-form four-dimensional
formula. Floating-point calculations in the model; not observations.*

## Where the research stands (2 October 2026)

| Result | Status |
|---|---|
| A specific five-dimensional model is fixed in advance ("registered"): Einstein gravity plus a scalar field, with W(φ)=1−φ+φ³/3 and shell tension σ=2W+δ(1+cφ), δ=0.001 | Definition |
| Within the model, the registered shell solution exists (computer-assisted proof with local uniqueness; a statement about the equations, not about nature) | **Certified**, conditional on its listed premises |
| The shell has an unstable (tachyonic) scalar mode, m² ≈ −7.7179 H², growth e^{1.657 Hτ}; numerically it is the only one. No unstable mode exists in the tensor, vector or special-harmonic sectors (mode stability, not full linear stability). | **Certified**, conditional (the scalar mode exists) / **Numerical** (uniqueness) / **Exact** (other sectors) |
| A closed-form 4D effective theory reproduces the 5D instability to about 10⁻⁶ | **Exact + numerical** |
| Full nonlinear 5D evolution: the shell rolls off toward one of **two fates**, relaxation toward an empty de Sitter brane (endpoint not yet reached in the simulations) or reversal and collapse | **Numerical** |
| Neither fate produces a radiation-filled (hot Big Bang) universe **in the homogeneous, classical roll-off of this model with no added matter** | **Negative result** |
| The earlier late-time endpoint (constant φ=1) fails a boundary condition; a corrected static solution was found | **Exact** (the failure) / **Numerical** (the new solution) |
| A proposed extension adds a new matter field to the shell, with an exact energy-exchange law; later checkpoints test particle production under stated approximations | **Exact**, as a proposal |
| An earlier (July 2026) gravitational-radiation claim was **withdrawn** after an audit | **Published correction** |
| A pulsar-timing "knee" signature was proposed and screened against public data | **Screening only**; not derived from the 5D model; the strict registered test failed on pilot data |
| **New (27 Sept):** the corrected end state (the "+1 branch") is linearly stable in the sectors tested; two independent methods agree, each first calibrated on the known instability | **Numerical**, independently audited (not a proof) |
| **New (27 Sept):** exact small-parameter series for the end state through eighth order, confirmed to 30–40 digits | **Exact**, independently re-derived |
| **New (27 Sept):** particle production and six other heating routes screened; none gives a radiation era in the registered model. The obstacle is leftover vacuum energy (about 0.6 of the ~22 e-folds needed) | **Negative**, under stated assumptions |
| **New (27 Sept):** an added tension term cancels the leftover vacuum, but "dark radiation" from the fifth dimension stays above observational limits so far | **Inconclusive** (requires a model change) |
| **New (28 Sept):** with the tension tuned to cancel the leftover vacuum (model change) and a stand-in (friction) energy-transfer formula, a radiation-dominated phase appears; "dark radiation" is 7.3% of the radiation at transfer rate Y = 1, within the conservative limit. A 1% change in the tuning removes it | **Numerical**, conditional, fine-tuned ([verdict](research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum/A1_VERDICT.md)) |
| **New (30 Sept):** at Y = 2 the dark radiation falls to 2.1%, within the strict current limit (3%). The phase lasts under half an expansion e-fold before leftover negative vacuum energy recollapses the shell, needs the tuning to about 1 part in 10⁴, and still uses the stand-in formula | **Numerical**, model-internal, conditional |
| **New (30 Sept):** replacing the stand-in with particle production derived from the matter extension: with particle masses below the 5D gravity scale, dark radiation stays far larger than the radiation; no steady state in 32 runs | **Inconclusive**, leaning negative |
| **New (30 Sept):** the same test at the registered δ = 0.001 is not completed; the numerical breakdown is diagnosed but not fixed. A simple 4D model does not reproduce the 5D dark-radiation numbers | **Inconclusive** / **Negative** (4D shortcut) |
| **New (1 Oct):** a single-cohort energy envelope includes radiation and undecayed particles. The two locally screened archived crossing cases remain far below the radiation target at the tested endpoint, even with optimistically timed decay. This uses a fixed history and a cutoff screen, not a complete quantum calculation | **Analytical bound + exploratory numerical application**, conditional |
| **New (1 Oct):** an executed audit of saved trajectories finds an earlier radiation-dominated sampled interval of about **0.314 e-fold** in the finer Y=2 run under a new cancellation-resistant diagnostic. Its later plateau fails this diagnostic. The original registered verdict is unchanged | **Post hoc saved-data diagnostic**; stand-in source, sampled duration |
| **New (1 Oct):** the particle-mode solver passes 90 exact pulse occupation checks, including reflectionless cases, and 30 fine-resolution Wronskian checks | **Registered numerical control**; not a coupled 5D matter result |
| **New (1 Oct, quantum modes):** actual scalar modes on the archived source-free tuned shell give a particle count 0.0503% below the earlier corrected estimate for the declared finite-band state. Coherence adds 10.96% to the particle-only scalar current and reverses the relative pressure's sign. Ten registered variants and independent RK4 checks pass | **Numerical**, prescribed background and relative stress; not a radiation era |
| **New (1 Oct, quantum matching):** the proposed scalar loop requires additional quartic-potential and intrinsic-curvature terms in the shell action | **Analytical EFT requirement**; curved-shell matching coefficients and absolute source unresolved |
| **New (2 Oct):** a flat-background quantum source through a smooth mass-zero crossing passes eight registered variants, exact scattering, energy/work and independent pressure-trace checks. Potential and curvature matching are fixed at a stated reference | **Numerical benchmark**, prescribed Minkowski background; not coupled shell evolution |
| **New (2 Oct):** omitting curvature matching makes pressure drift with the regulator even though energy conservation passes. A direct physical-mode representation has the expected cutoff convergence and an exact finite-cutoff trace identity | **Analytical and numerical checks**; finite-PV sources retain several-percent regulator effects |

![Two fates of the unstable shell at the registered parameters](hdblast/figures/two_fates_registered_detuning.png)

*Full five-dimensional numerical relativity at the registered parameters. The shell heads
toward an empty de Sitter brane (blue) or reverses and collapses (orange). Neither branch
contains a hot radiation era.*

The registered homogeneous model has not produced a hot Big Bang. Tuned variants with a stand-in transfer formula meet the earlier radiation-fraction criterion briefly and then recollapse. The first 1 October checkpoint bounds the energy budget and diagnoses a sampled 0.314-e-fold interval. The new quantum-mode calculation tests the proposed matter field on the source-free tuned geometry: its integrated count agrees closely with the corrected crossing estimate, but coherent stress/current matter and the endpoint excitations are not radiation. The 2 October checkpoint now provides an absolute source in a fixed matching convention for a flat-space pulse, including the curvature contribution to pressure. Curved-shell matching, coupled evolution, decay and thermalization remain unresolved.

Latest: [matched quantum source and pressure checks](research/HDBLAST_CHECKPOINT_20261002_PV/00_READ_FIRST.md), [reproduction](research/HDBLAST_CHECKPOINT_20261002_PV/REPRODUCE.md), [scientific ZIP](research/HDBLAST_CHECKPOINT_20261002_PV/HDBLAST_CHECKPOINT_20261002_PV.zip), and [targeted literature review, including 2026 papers](research/HDBLAST_CHECKPOINT_20261002_PV/LITERATURE_20261002.md).

Previous: [quantum modes, results and limitations](research/HDBLAST_CHECKPOINT_20261001_QFT/00_READ_FIRST.md), [reproduction](research/HDBLAST_CHECKPOINT_20261001_QFT/REPRODUCE.md), and [downloadable scientific package](research/HDBLAST_CHECKPOINT_20261001_QFT/HDBLAST_CHECKPOINT_20261001_QFT.zip).

Earlier: [energy bounds and exact pulse controls](research/HDBLAST_CHECKPOINT_20261001/00_READ_FIRST.md), including a [targeted review of five recent primary papers](research/HDBLAST_CHECKPOINT_20261001/LITERATURE_REVIEW_20261001.md). The new checkpoint adds [five foundational sources on stress/current renormalization and states](research/HDBLAST_CHECKPOINT_20261001_QFT/LITERATURE_AND_RENORMALIZATION.md).

## What would count for or against the hypothesis

For the hypothesis to become credible, all of the following must happen, and it may fail at any step:

1. a consistent matter sector in which the blast produces particles, with their gravity fed back into all equations;
2. a sustained radiation-dominated era with realistic temperature and duration;
3. stable, convergent five-dimensional simulations, with energy fully accounted for;
4. a quantitative prediction that competing explanations do not make;
5. independent reproduction, peer review and observational tests.

## Explore the work

| Where | What |
|---|---|
| [`hdblast/`](hdblast/) | Start here: the plain-language guide, program map, claims ledger and Zenodo records |
| [`hdblast/checkpoints/`](hdblast/checkpoints/) | Dated research checkpoints (Sept 2026) with reports, code, data and internal AI referee reviews (not external peer review) |
| [`research/HDBLAST_CHECKPOINT_20260927/`](research/HDBLAST_CHECKPOINT_20260927/) | The 27 Sept 2026 checkpoint: stability, exact series, particle production, simulation-error control, new mechanisms and a 2022–2026 literature review, each with an internal audit report. Start with `00_READ_FIRST.md` or `PUBLIC_SUMMARY.md`. Raw simulation arrays (`.npz`) are not included. |
| [`research/HDBLAST_CHECKPOINT_20260928/`](research/HDBLAST_CHECKPOINT_20260928/) | The 28 Sept 2026 tuned-vacuum test (A1): pre-registration, solver, analysis and verdict. |
| [`research/HDBLAST_CHECKPOINT_20260930/`](research/HDBLAST_CHECKPOINT_20260930/) | The 30 Sept 2026 checkpoint: derived particle production, the δ = 0.001 diagnosis, the transfer-rate scan, a 4D cross-check and a fine-tuning measure, each with an independent audit and a final critic pass. Start with `00_READ_FIRST.md`. Raw arrays (`.npz`, `.pkl`) are not included. |
| [`research/HDBLAST_CHECKPOINT_20261001/`](research/HDBLAST_CHECKPOINT_20261001/) | The 1 October checkpoint: cohort bounds, radiation diagnostics, exact mode controls, executable code and internal audits. |
| [`research/HDBLAST_CHECKPOINT_20261001_QFT/`](research/HDBLAST_CHECKPOINT_20261001_QFT/) | Quantum modes on the archived source-free shell, coherent stress/current, EFT matching, raw inputs, figures and independent replay. |
| [`research/HDBLAST_CHECKPOINT_20261002_PV/`](research/HDBLAST_CHECKPOINT_20261002_PV/) | Matched flat-background quantum stress/current, regulator and cutoff tests, exact finite-cutoff trace, figures and reproducible package. |
| [`site/`](site/) | Source of the website |

Most checkpoints can be re-run with Python 3 (`numpy`, `scipy`, `sympy`, `mpmath`). Each
checkpoint has its own `REPRODUCE`/`REPLAY` instructions. The author's raw working archive
(drafts, chats and earlier packages) is kept private; the curated, citable material is here.

## How to cite

Maldonado, R. *HDBLAST: higher-dimensional blast research program.* Zenodo,
[doi:10.5281/zenodo.17088132](https://doi.org/10.5281/zenodo.17088132) (all versions). See
[`CITATION.cff`](CITATION.cff). Each Zenodo record states its own licence.

**License for this repository:** text, figures and data are licensed under
[CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); code is licensed under the MIT License.
Reuse is welcome with credit to Ricardo Maldonado. See [`LICENSE`](LICENSE).

## A note on honesty

This project publishes its failures and corrections alongside its successes. A result is
called established only when a saved script or proof supports it. The reviews inside the
checkpoints are internal and were run with AI assistants; **no external peer review has taken
place yet.** Specialists in general relativity, numerical relativity and cosmology are warmly
invited to check the work and report errors by opening a GitHub issue.

*Research assistance: parts of the 2025–2026 calculations and documents were prepared with AI
assistants under the author's direction: Claude (Anthropic) and OpenAI ChatGPT/Codex. The October checkpoints were prepared with AI
assistance and independent internal computational checks. Verification steps are recorded in each checkpoint.*
