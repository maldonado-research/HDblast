# HDBLAST Chat 14 — the roll-off at the registered detuning: same two fates

22 September 2026 · Ricardo Maldonado's higher-dimensional blast research · prepared with Claude (Anthropic)

**Chat 13 computed the full five-dimensional evolution of the unstable shell at tension detunings of 0.1 and 0.03,
because the code could not run the registered value 0.001. That caveat is now removed.** With a stretched grid and a
fourth-order boundary treatment the same code runs at `t = 0.001`, reproduces the certified linear growth rate to five
digits, and finds the same two fates: relaxation toward an empty Randall–Sundrum de Sitter brane on one side, reversal
and collapse on the other. No branch produces a radiation-filled universe.

This removes the detuning caveat from the Chat 9–13 sequence at the parameters you actually registered. It remains a
negative result for the shell's position mode as the sole "blast" mechanism; it does not bear on higher-dimensional
origins with ingredients the model lacks. Two referee agents reviewed the code and this document; their corrections are
incorporated and listed in Section 4.

![The two fates at the registered detuning](TWO_FATES_REGISTERED_v7.png)

Solid curves: `t = 0.001` (the `+1` curve is the run with the refined far region, shown to the end of its record; the
throat curve is the default-grid run). Faint curves: the Chat 13 result at `t = 0.1`, shifted in brane time by −0.04
(`+1` side) and −0.09 (throat side) so that the two roll-offs cross `|φ_b| = 0.05` together.
[Vector version](TWO_FATES_REGISTERED.svg).

## 1. What was fixed

Chat 13's referee attributed the failure at `t = 0.001` to under-resolution of the wall (0.013 wide in the conformal
coordinate at `ρ_b = 79`) and suggested a fourth-order one-sided closure as one possible remedy. Both were needed
([rolloff5d_v2.py](rolloff5d_v2.py), which imports the Chat 13 code
for the model and the static background):

* **A stretched grid.** A smooth coordinate map puts 75–134 grid points across the wall and coarse spacing far from
  it (1,400–5,200 points in total instead of the 60,000 a uniform grid would need).
* **A fourth-order one-sided boundary closure at the shell.** The symmetric ghost-point closure of Chat 13 makes the
  second derivative at the boundary only first-order accurate, and at `t = 0.001` — where the physics is a `10⁻⁴`
  residual of two large, nearly cancelling terms — that error dominates: it gives a growth rate of 1.90 at 75 points per
  wall, 1.71 at 134, and 1.6565 (still rising) at 222 in the referee's own run. The one-sided closure (first derivative
  fourth-order, second derivative third-order at the boundary) gives 1.6572 at 75 points and at 134.

Also new: the seed is tapered smoothly at the far boundary (Chat 13's clamp launched a weak inward front), and the
time-derivative fields receive consistent boundary slopes, so the shell-layer constraint monitor is meaningful (the
momentum constraint vanishes identically at the shell under the junction data; it stays below `10⁻³` in all runs).
The numerical dissipation in this version is applied per grid index and is therefore about 2,500 times weaker than in
Chat 13 — effectively off; the runs are stable and show no grid noise regardless, but the code's docstring now says so.

**Validation at `t = 0.001`.**

| Test | Result |
|---|---|
| Linear growth rate, 75 points/wall (seed `Δc = 10⁻⁸`, sliding windows `t = 1–5`) | 1.65714–1.65717 |
| Same, 134 points/wall | 1.65718–1.65719 |
| Same, seed `Δc = 3×10⁻⁸` (referee) | 1.65714 — seed-independent |
| Same, 40 points/wall (referee) | 1.65640 — error 27× larger than at 75 points, consistent with `h⁵` convergence |
| Certified linear theory (Chats 9–11): `μ² = −7.7178716` | 1.65719 |
| Reproduction of Chat 13 at `t = 0.1` with the new grid | `φ_b(t = 5) = 0.20386` (Chat 13: 0.20386); zero-seed rate 1.62697 (Chat 13: 1.62694) |
| `+1` side at `H₀τ = 6.90`, wall resolution 75 vs 134 points, default far grid | `φ_b = 0.9946` both; `H_J/H₀ = 0.6111 / 0.6118` — but see the next row |
| `+1` side, far-grid spacing `0.008` vs `0.002` (same wall resolution) | agree to `5×10⁻⁵` at `H₀τ = 6.0`, 0.2% at 6.5, then separate: at 6.9 `0.611` vs `0.639`; the refined run keeps the bulk constraint below `10⁻²` while the default grid's rises to 0.65 |
| Throat side, wall resolution 75 vs 134 points | `φ = −1` crossing at `H₀τ = 5.6085 / 5.6083`; turnaround at `5.8947 / 5.8944` with `φ_b = −1.9417 / −1.9412` (interpolated); `H_J = −5H₀` at `6.200 / 6.199`; agreement 1–2% down to `−5H₀`, 8% at `−10H₀`, diverging beyond |
| Throat side, far-grid spacing `0.008` vs `0.002` | `φ = −1` crossing `5.6085 / 5.6085` (`H_J = 0.7978` both); turnaround `5.8947 / 5.8947`, `φ_b = −1.9417 / −1.9416`; `H_J = −5H₀` at `6.1999 / 6.1995`, `−10H₀` at `6.283 / 6.285` — the far region is irrelevant here |

The lesson of the second `+1` row: the domain wall recedes into the bulk behind the shell, and once it leaves the finely
gridded zone it must still be resolved. The default far spacing (`0.008`, under two points per wall width) is not enough
there; a fourfold refinement brings the constraint down by two orders of magnitude and changes the late `H_J` by 5%.
The refined run is one refinement step short of demonstrated convergence for that late value, which is therefore quoted
to one digit. Earlier in the roll-off (before `H₀τ ≈ 6.5`) and on the throat side the far region does not matter: the
throat-side run with the refined far grid reproduces the default-grid crossing, turnaround and collapse to four digits.

## 2. The two fates at `t = 0.001`

**Toward `φ = +1`.** The shell crosses the wall in about one Hubble time, reaches `φ_b = 0.99` by `H₀τ ≈ 6.8`, and its
expansion rate falls monotonically, without oscillation. When the chart's lapse begins to collapse (`H₀τ ≈ 6.9–7.0`;
the brane's proper time stops advancing at about 7.05) the rate is `H_J/H₀ = 0.63–0.64` and still decreasing at about
`−0.07 H₀²`. The exact Randall–Sundrum de Sitter brane of this tension in the `φ = +1` vacuum has `H/H₀ = 0.6067`,
which is also Chat 9's four-dimensional prediction (0.606). The 5D value is 4–5% above it, consistent with relaxation
toward it but not demonstrating it; the asymptotic value beyond the chart freeze is not determined. The same situation
was found at the larger detunings:

| detuning `t` | 5D value of `H_J/H₀` at the end of the reliable window | exact RS brane in AdS₊ |
|---|---|---|
| 0.1 | 0.649 (still falling) | 0.638 |
| 0.03 | 0.63 (still falling) | 0.616 |
| **0.001** | **0.63–0.64 (still falling)** | **0.6067** |

No radiation, no reheating: the shell becomes another empty inflating universe with about 62–64% of its original
expansion rate.

**Toward the throat.** The shell crosses `φ_b = −1` at `H₀τ = 5.61` with `H_J = 0.798 H₀` (Chat 9's four-dimensional
theory: 0.797). Beyond `φ_b = −1/c = −1.67` the detuning term of the tension, `t(1 + cφ_b)`, turns negative — the brane's
tension falls below its balanced value — and at `H₀τ = 5.89` the brane's expansion reverses, with `φ_b = −1.942` at both
wall resolutions (four-dimensional prediction, derived for `t → 0`: `−1.942`). The contraction is kinetic-dominated and
the five-dimensional curvature at the brane grows without bound: `H_J` passes `−5H₀` at `H₀τ = 6.20` (two resolutions
within 1.4%) and `−30H₀` before `6.32`, beyond which the resolutions diverge; the classical evolution ends in a
singularity of unresolved form within about 0.43 Hubble times of the turnaround, with the brane's proper time saturating
near `H₀τ ≈ 6.33`.

For comparison, at `t = 0.1` the turnaround was at `H₀τ = 5.99` with `φ_b = −1.942` when interpolated the same way (Chat
13 quoted the first recorded sample, `−1.96`). The registered detuning behaves the same way, slightly earlier in proper
time.

**On the constraint monitor.** For `t < 3.5` the dominant bulk violation is an outgoing pulse launched from the shell by
the seed's initial junction mismatch, travelling at light speed (`z = −t`). Later the maximum sits at `z ≈ −0.9` and grows
with the unstable mode; on the default far grid it reaches 0.65 on the `+1` side by `t = 8`, on the refined far grid it
stays below `10⁻²`. This feature is not fully explained, but it lives behind the shell and the refined-grid results show
it does not alter the fates.

## 3. What this does and does not establish

* **Numerically supported (floating point, two wall resolutions, two far-region resolutions on the `+1` side, referee
  re-derivations and independent runs):** at the registered parameters, the homogeneous, classical, single-scalar,
  matter-free roll-off of the unstable shell ends either in an empty de Sitter brane or in a collapse. There is no
  radiation-dominated branch.
* **Not established:** the asymptotic value of the `+1` plateau beyond the chart freeze (0.63–0.64 at the freeze, one
  refinement step short of converged, 4–5% above the exact RS value); the final approach to the singularity; the origin
  of the late constraint feature at `z ≈ −0.9`; anything inhomogeneous, quantum, or with brane matter.
* **No novelty claimed.** Stretched grids and one-sided closures are standard numerical relativity; the outcomes are of
  the BraneCode type (Martin, Felder, Frolov, Peloso, Kofman, hep-th/0309001).

## 4. Referee corrections

Two agents (numerics; claims) re-derived the closure algebra, reproduced every quoted number from the saved data, ran
four further evolutions and read this document line by line ([numerics report](referee_numerics/REFEREE_NUMERICS.md),
[claims report](referee_claims/REFEREE_CLAIMS.md)). They confirmed
the code, the growth rate and the two fates, and corrected the first draft in these places, all fixed above:

1. The `+1` plateau had been quoted as 0.61 from runs whose far region under-resolves the receding wall; the refined run
   gives 0.63–0.64 with a hundredfold smaller constraint, and the "two resolutions" quoted shared the same far grid.
2. The turnaround `φ_b = −1.953` was a sampling artefact (first record after the crossing); interpolated it is `−1.942`,
   matching the four-dimensional prediction to three digits — better than claimed. Likewise the crossing value is 0.798.
3. "Converged to `−10H₀`" was too strong (8% there); "`−300H₀` by 6.32" was single-resolution.
4. The outgoing-pulse description of the constraint monitor holds only for `t < 3.5`.
5. Closure orders and the "150 points" figure were misstated; the dissipation strength in v2 was misdescribed.
6. Chat 13's referee had proposed the one-sided closure as a remedy and is credited; the earlier "BPS cancellation"
   wording is not revived.
7. Phrases such as "closes the last open item" were softened; this checkpoint removes a caveat, it does not close the
   programme.

## 5. Where the programme stands

Chats 9–14 give a refereed account of the shell in the registered model: the registered shell is unstable
(interval-certified), stable shells exist and are certified, a closed-form four-dimensional theory links them, the
equations underneath are verified symbolically, and the non-linear fate of the instability is known at the registered
parameters. The remaining options were listed at the end of Chat 13: a consolidated Zenodo supplement of Chats 9–14
(recommended next), or new physics — brane matter with a derived coupling, or a tension with a trapping minimum — to look
for a reheating branch.

## Reproduction

In `` with the system `python3` (numpy), single-threaded
(`export OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1`; several processes with default threading oversubscribe the machine):

    python3 rolloff5d_v2.py --tdet 1e-3 --dc 1e-8 --closure 4 --L 10 --tf 14 --tag runs/rate                        # ~15 min: linear growth rate
    python3 rolloff5d_v2.py --tdet 1e-3 --dc 1e-4 --closure 4 --dzc 2e-3 --L 10 --tf 11 --phimax 0.9999 --tag runs/plus   # ~25 min: +1 side, refined far grid
    python3 rolloff5d_v2.py --tdet 1e-3 --dc=-1e-4 --closure 4 --L 10 --tf 8 --phimin -40 --tag runs/minus              # ~8 min: throat side
    python3 analyse_v2.py

`runs/ANALYSIS_V2.json` holds the fitted numbers (its turnaround entries are first-record samples; the referee's
`referee_numerics/ref_analyse.py` interpolates). Logs of runs interrupted by session restarts are kept with the suffix
`_KILLED`; their trajectories agree with the reruns. Nothing was published or sent; no original research file was modified.
