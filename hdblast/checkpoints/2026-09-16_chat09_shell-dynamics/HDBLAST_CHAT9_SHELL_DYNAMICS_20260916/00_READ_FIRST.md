# HDBLAST Chat 9 — the candidate universe is unstable, stable ones exist, and the 5D problem now has a 4D formula

16 September 2026 · Ricardo Maldonado's higher-dimensional blast research · prepared with Claude (Anthropic)

**This continuation asked a question the project had not yet asked: what does the registered de Sitter shell — the
candidate "universe" of the model — actually do?** The answer is a certified mathematical result with a clear physical
reading. The registered shell is unstable: it carries exactly one growing mode, with `m² = −7.7179 H²`, so it falls
off its equilibrium in well under one Hubble time. The same analysis produced a closed-form four-dimensional
description of the shell that reproduces the five-dimensional numbers to one part in a million, and it identified
brane tensions for which the shell is **stable**.

These results advance the registered model. **They do not show that a higher-dimensional blast caused the Big Bang.**
In fact the follow-up calculation is negative on that point for the shell as currently registered: its roll-off does
not produce a hot, radiation-filled universe. That is stated plainly in Section 5, and it tells us what to change.

![The registered shell sits on a hilltop; the 5D spectrum converges to the closed form](figures/SHELL_HILLTOP_AND_SPECTRUM_v3.png)

Both panels are calculations inside the registered dimensionless model, not observations.
[Vector version](figures/SHELL_HILLTOP_AND_SPECTRUM.svg).

## 1. What the registered model is

Reading the M462 shell problem and the forward-experiment contract together identifies the registered bulk as a
five-dimensional Einstein–scalar model built from a superpotential,

    W = 1 − φ + φ³/3,     U = W'²/2 − (2/3)W²,

with two anti-de Sitter vacua (φ = −1 and φ = +1). With the detuning `t` switched off the background is the exact
flat wall `φ = −tanh y`. The registered solution is a mirror-symmetric de Sitter 3-brane ("shell") whose tension is
detuned from its balanced value, `σ_t = 2W + t(1 + c★φ)`, `t = 10⁻³`. A new global solver
([background/](background/hdblast_background.py)) integrates from the first cone
to the shell and reproduces the certified M462 root: δ to 1×10⁻⁹, η to 4×10⁻⁸. The v24 certified response lives in
a tiny corner of this spacetime (σ ≤ 4×10⁻⁷); the shell sits at σ ≈ 12.45.

## 2. Main result — a certified tachyon on the registered shell

For scalar perturbations of the shell with four-dimensional mass `m² = μ²H²`:

    exactly one bound state below the 9H²/4 continuum, and it is tachyonic:   μ² = −7.7178716.

Disturbances grow as `e^{1.657 Hτ}`: an e-folding time of 0.60 Hubble times.

| Check | Who / how | Result |
|---|---|---|
| Longitudinal-gauge derivation + shooting | orchestrator | −7.7178716187 |
| Gaussian-normal-gauge derivation, blind to the first | independent agent | −7.7178716252 |
| Published Frolov–Kofman equations, transcribed | literature agent | −7.7178716236 |
| Own integrator + brute-force linearisation of all 15 Einstein components (residual ≤ 1.4×10⁻⁹) | adversarial verifier | −7.717871623 |
| **Interval-arithmetic certificate** (74 gates, replayable) | certificate agent | **−7.71788 < μ² < −7.71786** |
| Zero count in the complex μ² plane (Re ∈ [−400, 2], \|Im\| ≤ 300) | orchestrator, verifier | exactly one |
| Sign of the norm | two agents, two methods | positive: a tachyon, not a ghost |

The certificate ([certified_instability/](certified_instability/)) is conditional
on the archived M462 root enclosure and on the perturbation equations; those equations were derived in two gauges and
checked numerically against the full field equations. Uniqueness of the mode and positivity of its norm are supported
numerically and by analytic identities, not interval-certified. The tensor (graviton) sector is healthy: one massless
graviton, `M₄² = 2.0677`, no other bound state ([tensor_and_calibration/](tensor_and_calibration/TENSOR_SECTOR_AND_CALIBRATION.md)).

## 3. A closed-form four-dimensional description

Expanding the radial Hamilton–Jacobi equation of the bulk gives an effective action for the shell's position `φ_b`:

    S₄ = ∫√−g [ f R/2 − ½ Z (∂φ_b)² − V ],    f = 2I(φ_b),   Z = −W f'/W',   V = σ_t − 2W,

where `I` solves `W' I' − (2/3) W I = −1`. Consequences, all checked
([derivation](effective_theory_orchestrator/HJ_EFFECTIVE_THEORY_DERIVATION.md),
[exact-arithmetic checks](exact_checks/verify_exact_identities.py)):

* **It explains the registered constant.** The condition for a static shell at `φ_b = 0` is exactly
  `c★ = 2/I₊ − 4/3`.
* **It predicts the tachyon in closed form:** `μ²(t→0) = −4(3c² − 4c + 8)/(c(3c + 4)) = −7.719796`, equivalently
  `μ² = −4 − Δ` with `Δ = 3.72`.
* **The next order was derived analytically** ([slope_analytic/](slope_analytic/SLOPE_DERIVATION.md)):
  `μ² = −7.719796 + 1.9243896 t + O(t²)`; five-dimensional shooting gives 1.9243896(5). The four-dimensional formula
  and the five-dimensional calculation agree to about 10⁻⁶ in five different configurations.
* The five-dimensional response of the shell to a source reproduces the predicted kinetic normalisation to 10⁻⁵.

The method is established (de Boer–Verlinde–Verlinde; Brax–van de Bruck–Davis–Rhodes). The explicit results for this
model appear to be new to the project and were not found in a bounded literature search; no external novelty is claimed.

## 4. Stable shells exist

The instability belongs to the registered *linear* detuning. With a curved tension,
`σ_t = 2W + t(1 + c★φ + dφ²/2)`, the full five-dimensional problem at `t = 10⁻³` gives:

| d | scalar bound states below 9H²/4 | 4D prediction (t → 0) |
|---|---|---|
| 0 (registered) | −7.7178716 | −7.7198 |
| 1.0 | −0.7895812 | −0.7869 |
| 1.3 | +1.2949033 | +1.2930 |
| 1.4 | +1.9871330 | +1.9863 |
| 1.6, 2.0 | none | above the continuum |

A finer scan at `t = 10⁻³` ([STABILITY_WINDOW_SCAN.json](stability/orchestrator_derivation/STABILITY_WINDOW_SCAN.json))
locates the window: the mode is tachyonic up to `d = 1.11` (μ² = −0.077), positive from `d = 1.12` (+0.086), and it
merges into the continuum between `d = 1.43` (μ² = 2.195) and `d = 1.44` (no bound state) — the 4D formula predicts the
merger at `d₀ + ¾Z_E = 1.438`. Seventeen values of `d` follow the 4D prediction to better than 1% except within ±0.02
of `d₀`, where the finite-`t` shift of the equilibrium matters.

So for `d > d₀ = 1.1135` the shell is a stable de Sitter braneworld. The referee found this; the orchestrator
reproduced it independently ([stable_shell_recheck.py](stability/orchestrator_derivation/stable_shell_recheck.py)).
This is the constructive outcome of the checkpoint: it says exactly which shells are worth registering next.

## 5. What the roll-off does — and why it is not yet a Big Bang

Within the four-dimensional description ([blast_dynamics/](blast_dynamics/BLAST_DYNAMICS_REPORT.md)):

* The shell leaves the hilltop after about 3.4 e-folds (2.6 if it starts at the de Sitter bounce); the roll itself
  lasts about one Hubble time.
* **Toward φ = +1** it settles, without oscillating, into another de Sitter brane with `H = 0.606 H_top`. Acceleration
  never ends, there is no radiation era, and particle production is tiny: at most 2.4×10⁻⁸ of the energy for a
  minimally coupled field, 5×10⁻⁵ with a strong explicit coupling.
* **Toward the AdS throat** the effective potential eventually turns negative and the brane-frame expansion reverses.
  The five-dimensional fate there is unknown.
* Tuning the tension to make the hilltop flat does not rescue inflation: the potential has a large cubic term, so
  `n_s ≤ 1 − 4/N ≈ 0.93`, excluded by Planck, and inflation never ends.
* A modified tension with an interior minimum at zero energy does reheat by modulus decay
  (`T_rh ≈ 3 (H/M_Pl)^{1.5} M_Pl`), but with at most ~10 e-folds the horizon and flatness problems are unsolved and a
  closed shell recollapses.

**Plain conclusion:** the registered shell's instability is a real, calculable higher-dimensional event, but as it
stands it is a transition between two vacuum states, not the origin of a hot universe. A viable version needs
(i) a long quasi-de Sitter phase from some other ingredient, and (ii) an exit into radiation. Neither is present yet.

## 6. Other results

* **Moving-shell theorem** ([moving_shell/](moving_shell/MOVING_SHELL_FRIEDMANN_THEOREM.md)):
  for any shell trajectory, `H² + 1/a² = (λ+ρ)²/36 + U(φ_b)/6 − (∇φ)²_b/12` (the Maeda–Wands equation with vanishing Weyl
  term); the static case is exactly M462 (7.1). A moving shell with a prescribed tension is over-determined in the frozen
  bulk — the bulk must respond — and all vacuum shells are time-shifted hyperboloids with a closed-form tension.
* **Entropy identity** ([instanton_orchestrator/](instanton_orchestrator/SHELL_INSTANTON_ENTROPY_IDENTITY.md)):
  the Euclidean action of the shell solution is exactly `−16π²∫ρ² dy = −8π²M₄²/H²`, minus the de Sitter entropy and
  minus the five-dimensional horizon area over 4G₅. Because `μ² < −4` the saddle has 1 + 5 negative modes in the
  four-dimensional counting, so no tunnelling-probability reading is offered.
* **Physical scales** for an assumed shell Hubble rate (the model has no preferred scale; these are conversions, not
  predictions): `H = 10¹³ GeV` gives `M₅ = 1.3×10¹⁷ GeV`, `L₀ = 2.5×10⁻³¹ m`, growth time 4×10⁻³⁸ s; every case satisfies
  the laboratory bound on the AdS length by more than 17 orders of magnitude.

## 7. Corrections made during review

The referee found errors in the orchestrator's first notes; all are fixed in the files and recorded here.

1. "A tuned tension gives n_s ≈ 0.965" was wrong (cubic term ignored; finite-t floor `μ² ≤ −2.31√t` ignored).
2. "Field space ends at the throat" was wrong: `Z(φ_b = −1) = 9/23` is finite and the description continues.
3. The φ = +1 end value is `(1+c)t/81 = 0.019723 t`, not 0.01976 t.
4. `REPLAY.py` had a syntax error introduced by a late edit; repaired and re-run (all five gates pass).

## 8. Prior art

The instability of inflating branes with bulk scalars is known (Frolov–Kofman hep-th/0309002; BraneCode
hep-th/0309001); the perturbation equations here are theirs, equation for equation, in the registered chart. The
"−4" is Gen–Sasaki / Garriga–Vilenkin. The effective action is the holographic-RG / moduli-space method. What a
bounded search (~30 papers, 2000–2026, including dark-bubble cosmology) did not find: the soft-boundary regime
`μ² < −4`, a single shell with a regular cone on a back-reacted wall, the closed forms, and the O(t) test. Citations
gathered by agents through web summaries must be re-checked by a person before publication
([literature/](literature/PRIOR_ART_SHELL_STABILITY.md)).

## 9. Status table

| Question | Status after this checkpoint |
|---|---|
| Does the registered shell have a tachyonic scalar mode? | **Certified** (interval arithmetic; conditional on M462 enclosure and the linearised equations) |
| Is it the only bound state, with positive norm? | Supported numerically and by analytic identities; not certified |
| Closed-form 4D description valid through O(t)? | **Derived and confirmed to ~10⁻⁶** |
| Stable de Sitter shells in the model class? | **Yes for d > 1.1135** (floating point, two independent runs) |
| Non-linear five-dimensional fate of the roll-off? | **Unknown** |
| Hot radiation universe from the roll-off? | **No, for the registered tension** |
| Viable inflation from the shell modulus? | **No** (n_s ≤ 0.93, no exit) |
| Higher-dimensional origin of the Big Bang? | **Not established** |
| External novelty? | Not claimed |

The v24 certified response and the PRD manuscript are untouched by this work.

## 10. What to do next

1. **Register and certify a stable shell** (for example `d = 1.6`), so the programme has a candidate universe that does
   not fall apart. The certificate machinery now exists.
2. **Five-dimensional non-linear evolution of the roll-off** (an SO(4)-symmetric 1+1 code with the shell as boundary).
   This is the decisive calculation for both end states and for the energy radiated into the bulk.
3. **Four-derivative Hamilton–Jacobi terms**, to decide where the 4D description fails during the fast roll.
4. A person should open the seven key papers listed in the referee review and confirm the equation correspondences.
5. Only after 1–3: look for an exit into radiation. If none survives, report the cosmological programme for this shell
   as negative and move to a different mechanism.

For Zenodo, this supports a clearly labelled **stability and effective-theory supplement** with its negative
cosmological conclusions included. It does not support any origin-of-the-Big-Bang claim. Nothing was published, sent or
uploaded; no original research file was modified.

## Reproduction

In ``: `python3 REPLAY.py` (system Python 3.9 + numpy, ~4 minutes, five gates);
`python3 certified_instability/REPLAY.py` replays the interval certificate (~6 minutes, standard library only).
Each track folder has its own scripts and JSON results. The full referee report is
[verification/REFEREE_REVIEW.md](verification/REFEREE_REVIEW.md).
Agents were independent AI reviewers, not external peer review. A few stale log files and two duplicate figure PNGs
could not be deleted by the sandbox; they are harmless and can be removed by hand.
