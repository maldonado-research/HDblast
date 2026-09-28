# Referee report (physics and claims) — Chat 13, nonlinear 5D roll-off

Referee: `referee_physics` (Claude, skeptical-referee role). Date: 21 September 2026.
Scope: `00_READ_FIRST.md` (one level up), `rolloff5d.py`, `derive_evolution_equations.py`, `analyse_runs.py`, `runs/*`.
Everything I computed is in this folder (`check_numbers.py` reproduces the numbers; `runs/` holds three short
check runs, 2e-3 resolution, ~11 min total). Nothing outside `referee_physics/` was modified.

Legend: **[EST]** established (derivation or converged numerics), **[NUM]** numerically supported but not converged /
short window, **[NOT]** not established, **[ERR]** error in the document.

---

## 0. Verdict in one paragraph

The calculation is sound and the two qualitative fates are real: at `t = 0.1` the homogeneous roll-off toward `φ = +1`
relaxes toward a slower de Sitter brane, and the roll toward the throat turns the brane's expansion around and sends it
into accelerating contraction. The linear regime and the turnaround are cleanly reproduced and match Chat 9's 4D EFT
(turnaround at `φ_b = -1.957` vs `-1.94`; `H_J/H₀ = 0.793` vs `0.797` when the brane crosses `φ_b = -1`). But the
document overstates what the runs show in four places: (i) the "new de Sitter" end state is observed for only
**0.32 e-folds / 0.49 Hubble times of brane proper time** before the chart freezes (at `H₀τ ≈ 7.39`, not 8.3), and
`H_J` is still falling; (ii) the quoted "simple estimate 0.603" is wrong — the exact de Sitter-brane-in-AdS₊ value at
`t = 0.1` is **0.638** (the analysis script dropped the `O(t²)` term of `σ²`), so the residual is 1.3 %, not 7 %, and is
consistent with incomplete relaxation, not with a discrepancy; (iii) the late collapse is **not converged** between
resolutions (`H_J` differs by a factor 2–3 at equal proper time after `φ_b ≈ -5`), the power-law `H_J ≈ -0.65/(τ*-τ)`
is a fit artefact (the local exponent wanders between 1.2 and 2.6 and the fitted `α` drifts 0.49–0.70 with the window),
and the 5D curvature at the brane is three orders of magnitude above the model's own AdS₋ curvature well before the
quoted `τ*`; (iv) the sentence "the scalar runs to `−∞` because `U` is unbounded below" is not the mechanism: the
turnaround happens where `U(φ_b) = +3.9 > 0` and `W > 0`; what drives it is the detuning term `t(1 + cφ_b)` turning
negative (tension below the BPS value) beyond `φ_b = -1/c = -1.673`, and the runaway of `φ_b` in the contracting phase
is the ordinary kinetic-dominated behaviour of a negative-energy collapse (Felder–Frolov–Kofman–Linde, hep-th/0202017),
which does not need `U` to be unbounded below. The bottom line "no branch of the homogeneous classical roll-off yields a
hot radiation universe at `t = 0.1, 0.03`" stands; "at the full nonlinear 5D level" must carry "homogeneous, classical,
one bulk scalar, no brane matter, `t ≥ 0.03`".

---

## 1. Is `H_J` the brane Hubble rate and `H₀τ` the brane proper time? — **[EST] yes**

Induced metric on the shell (`z = 0` in `ds² = e^{2B}(-dt²+dz²) + e^{2A}dx₃²`):
`h = -e^{2B(t,0)}dt² + e^{2A(t,0)}dx₃²`. Hence `dτ = e^{B_b}dt`, `a(τ) = e^{A_b}`, and
`H = (1/a)(da/dτ) = A_t e^{-B}` at `z = 0`. The code uses `A = ln ρ + t + a`, `B = ln ρ + b`, so `A_t = 1 + ∂_t a`,
`e^{B_b} = ρ_b e^{b_b}`; `HJ = (1 + pa[-1])/eB` (line 199) is exactly `A_t e^{-B}` at the last grid point, which *is*
the shell (Neumann data enter through ghost points, so no half-cell offset). The static solution has `A_t = 1`,
`e^{B_b} = ρ_b`, so `H₀ = 1/ρ_b` and the stored `HJ·ρ_b` is `H_J/H₀`. The `tau` accumulator (lines 211/213) adds
`½dt·e^{b_b}` before and after each RK4 step, i.e. a trapezoidal `∫e^{b_b}dt = H₀∫e^{B_b}dt = H₀τ`; second-order in
`dt = 5e-4`, negligible error. Residual gauge freedom of the conformal chart with the shell pinned at `z = 0` is a single
reparametrisation `t → f(t)`, under which both `H_J` and `τ` are invariant. The junction conditions
`A_z = B_z = (σ/6)e^B`, `φ_z = -(σ'/2)e^B` are the Z₂ Israel conditions for `K^i_j = -(σ/6)δ^i_j` with unit normal
`e^{-B}∂_z`; I re-derived them (Gauss–Codazzi, below) and they agree with the static check in
`EVOLUTION_EQUATIONS_CHECK.json`. "J" is a fair label: this is the frame in which brane-localised matter, if any, would
be minimally coupled (Chat 9's Jordan frame).

One consequence the document under-reports: because `e^{b_b} ∝ e^{-0.95t}` after the roll (measured slope
`db_b/dt = -0.951` at `t = 0.1`, `-0.863` at `t = 0.03`), the chart covers brane proper time only up to
`H₀τ_∞ = τ_end + e^{b_end}/0.951 = 7.385` (`t = 0.1`) and `7.23` (`t = 0.03`). **Not 8.3** (no run in `runs/` reaches
`H₀τ > 7.38`; `8.3` appears nowhere in the logs) **[ERR]**. In the same window the deviation front moves out at
`dz/dt = 1` (half-width `z₅₀ = -0.26, -1.41, -3.36, -5.36` at `t = 6, 8, 10, 12`) while `e^{2B}` collapses in the whole
region between brane and front: the chart is degenerating onto a null surface, i.e. `t = ∞` is a Cauchy horizon of the
coordinates through the brane event `τ_∞`, not the brane's future infinity. Everything after `H₀τ ≈ 7.39` on the
`+1` branch is uncomputed.

## 2. End state on the `φ = +1` side — **[NUM] consistent with the RS de Sitter brane; the quoted estimate is wrong**

Exact de Sitter brane in AdS₊ (Z₂, `κ₅ = 1`): `H² = σ_t(1)²/36 − 1/l₊²`, `1/l₊ = √(−U(1)/6) = 1/9`,
`σ_t(1) = 2/3 + t(1 + c)`. Expanding, `H² = (1+c)t/27 + (1+c)²t²/36`. Chat 9's `0.6064 = √((1+c)t/27)/H_top` is the
`O(t)` term only; `analyse_runs.py` (line 26) also uses only the `O(t)` term. With the actual `H₀ = 1/ρ_b`:

| `t` | `H₀ = 1/ρ_b` | `H₀²/t` | exact RS `H/H₀` | `O(t)` only | 5D run (last value) | window after `φ_b > 0.98` |
|---|---|---|---|---|---|---|
| 0.1 | 0.127612 | 0.1629 | **0.6379** | 0.6028 (doc: "0.603") | 0.6458 at `H₀τ = 7.378`, still falling | 0.49 Hubble times, 0.32 e-folds |
| 0.03 | 0.069604 | 0.1615 | **0.6161** | 0.6053 | 0.6372 at `H₀τ = 7.077` (run stopped by `φ_b > 0.995`) | 0.33 Hubble times, 0.22 e-folds |
| 1e-3 | 0.012686 | 0.1609 | 0.6067 | 0.6064 (Chat 9) | not run | — |

So the "5D value 0.646 vs estimate 0.603" gap in Section 2 of the document is an artefact of the estimate **[ERR]**; the
true gap is 0.646 vs 0.638 (1.3 %) at `t = 0.1` and 0.637 vs 0.616 (3 %) at `t = 0.03`, in both cases with `H_J` still
decreasing when the record ends. The `t = 0.03` run should not be quoted as "`H_J/H₀ ≈ 0.64`, behaves identically": it was
stopped by the `--phimax 0.995` guard while still rolling (`φ_b` had reached 0.995 and `H_J` was falling at
`−0.03 H₀` per Hubble time).

*Does `φ_b = 0.991 ≠ 1` matter?* No. Using the `nn` Gauss equation at the brane for a quasi-static state
(`Ḣ = φ̇ = 0`): `6H² = σ²/6 − σ'²/8 + U(φ_b)`, i.e. `H² = σ²/36 − σ'²/48 + U(φ_b)/6` (reduces to the RS formula at
`φ_b = 1`). At `φ_b = 0.9908`, `t = 0.1`: `σ' = 0.023`, `H/H₀ = 0.6375` — a 0.06 % shift, because
`σ²/36 + U/6` is stationary at `φ = 1` to `O(t)`. The bulk at the end is not AdS₊ but a nearly homogeneous
`φ ≈ 0.991` region (the whole first 0.6 in `z`), sourced by the brane through `φ_z = −σ'(φ_b)e^B/2` with `σ'(0.991) > 0`;
`φ_b = 0.991` is where `σ_t' ` is small enough to be balanced by the massive (`m² = U''(1) = 3.1`, `m l₊ ≈ 16`)
scalar's tail, so it is plausibly the true quasi-static position, not a transient.

*Is the residual relaxation understood?* Fitting `H_J² = H_∞² + C·a^{−n}` over `t ∈ [8, 11]` (`N` = e-folds since
`t = 8`) gives `H_∞/H₀ = 0.637` for `n = 4` (rms `2e-4`) and `0.642` for `n = 6` (rms `4e-4`). The `n = 4` fit lands on
the RS value to 0.1 %, with `C/H₀² = 0.032`, i.e. a ≈ 7 % radiation-like component when `φ_b` first exceeds 0.98,
redshifting away. In a homogeneous 5D evolution this is exactly the kind of term one expects (bulk Weyl / "dark
radiation" and scalar-gradient energy sourced during the roll). I record this as **[NUM]**: the window is half a
Hubble time and the last few points (`t > 11.5`, where `e^{b_b} < 0.01`) drift faster than any power law, which I
attribute to the degenerate chart (`t = L = 12` is also when the far boundary first reaches the brane; my `L = 16` check
run is in §7). The honest statement: the `+1` branch is heading toward the RS de Sitter brane with `H = 0.638 H₀`
(`t = 0.1`), it has been followed to within 1.3 % of it, and nothing in the record suggests any other attractor.

## 3. Throat side: turnaround **[EST]**, crunch **[NUM]**, stated cause **[ERR]**

*Turnaround.* `H_J = 0` at `H₀τ = 5.9912` (`dz = 1e-3`) / `5.9909` (`2e-3`), `φ_b = −1.9566 / −1.9562`: converged to four
digits. Chat 9's EFT (t-independent in Hubble units): `−1.94`. Crossing of `φ_b = −1`: `H_J/H₀ = 0.793` (5D) vs
`0.797` (EFT). Crossing of the EFT's `V_E = 0` point `φ_b = −1/c = −1.673`: `H_J/H₀ = 0.263`, i.e. the brane is already
decelerating hard there. The 4D EFT is therefore quantitatively right up to and through the turnaround at `t = 0.1`,
which is a stronger validation of Chat 9 than the document claims.

*Mechanism.* Check the numbers at the turnaround: `U(−1.96) = +3.90`, `W(−1.96) = +0.45`, `σ_t(−1.96) = +0.883`.
The potential is positive, the BPS tension is positive, and the shell has not reached `W = 0` (`φ = −2.104`), let alone
the sextic tail. What has changed sign is the detuning: `σ_t − 2W = t(1 + cφ_b) < 0` for `φ_b < −1.673` — the brane is
now *below* BPS tension. Chat 9's effective potential is `V_E = t(1 + cφ_b)/f²`, negative exactly there, and has a
**negative minimum** (`−0.0895 t` near `φ_b ≈ −2.9`), so along the EFT it is not unbounded below at all. The `O(t)`
quasi-static Gauss relation `H² = (t/36)(4W(1+cφ) − 3cW')` changes sign near `φ_b ≈ −1.3` (with the caveat, already in
Chat 9, that the frozen-bulk junction with `W' ≠ 0` is not the right description once `φ̇_b` matters). Once the effective
energy has passed through zero, a homogeneous universe with a canonical scalar contracts, the kinetic term grows as
`a^{−6}`, and the field runs logarithmically in `a` regardless of the potential's shape — this is the FFKL result
(hep-th/0202017, "practically independent of V(φ)", true for a bounded negative minimum as well as for an unbounded
potential). The 5D run shows exactly that signature: `dφ_b/d ln a = −2.0` (constant to ±5 %) through the contraction.
So: **the turnaround is caused by the brane's tension deficit relative to BPS (the detuning), and the runaway of `φ_b` is a
consequence of the contraction, not its cause.** The sentence "the scalar field runs off to `−∞` because the registered
potential is unbounded below" and the status-table line "scalar runs away because `U` is unbounded below" should be
replaced (wording in §6). The sextic tail becomes relevant only after the collapse is under way (`W < 0`, negative
tension `σ_t(−6.84) = −198`).

*Is it a big crunch?* What is established: after the turnaround the brane contracts by `Δln a = −1.9` (`dz = 1e-3`) /
`−2.2` (`2e-3`) with `|H_J|` reaching `40–84 H₀` in `ΔH₀τ = 0.43`, and the two resolutions agree to 3 % down to
`φ_b ≈ −5` (`H_J/H₀ = −7.45` vs `−7.70`, `Δln a = −0.82`). Not established:
(a) the asymptotics — at equal proper time `H₀τ = 6.415` the fine run has `H_J = −29 H₀`, the coarse `−84 H₀`; the
fitted `α` in `H_J ≈ −α/(τ*−τ)` is 0.70/0.63/0.49 (`dz = 1e-3`) for thresholds `H_J < −3/−10/−30`, 0.65/0.57/0.47 (`2e-3`),
and the local exponent `−d(1/H_J)/dτ` runs 2.6 → 1.2 → 2.1 along the record; neither the 4D kination value 1/3 nor the
`ρ²`-dominated brane value 1/6 nor any constant is seen. `τ*` itself is robust only because little proper time is left
(`H₀τ* = 6.43–6.44`, but the record ends at 6.415–6.419, so `τ*` is a 0.3 % extrapolation whatever the exponent).
(b) the nature of the singularity — the 5D Ricci scalar at the brane, `R = φ_n² − φ̇² + (10/3)U(φ_b)` with
`φ_n = −σ_t'/2`, grows from `0.4` (static) through `21` (turnaround) to `−3.4e3` (`φ_b = −5.6`) and `−1.6e4` at the end
of the record, i.e. `2.6e3` times the AdS₋ curvature `|R| = 6.2` that sets the model's own scale, and the proper
thickness of the region carrying the field excursion shrinks from 1.5 to 0.09 (5D units). This is a 5D curvature
singularity forming at the brane (compare BraneCode's Kasner-like strong-gravity regime), not merely a 4D `a → 0`;
either way it is beyond any classical control long before `τ*`. Constraint violation: 1.9 % (`1e-3`), 7.8 % (`2e-3`) —
the document's "`≲2×10⁻²`" holds for the fine run only.

Safe statement: *the throat-side brane reaches a turnaround at `φ_b = −1.96` and then contracts by at least two e-folds
with rapidly growing curvature; the classical 5D evolution ends in a singularity within `ΔH₀τ ≈ 0.44` of the turnaround
(extrapolated `H₀τ* ≈ 6.43`), whose detailed asymptotics are not resolved.* "Genuine big crunch with `H_J ≈ −0.65/(τ*−τ)`"
is not supported.

## 4. Is "no branch gives a hot radiation universe" justified? — **yes, with the qualifiers made explicit**

What the runs establish: for the *homogeneous, classical* 5D Einstein–scalar system with a pure-tension shell, at
`t = 0.1` (and, more weakly, `0.03`), starting from an `O(10⁻⁴)` seed on the unstable mode, neither sign of the seed leads to
an oscillating or radiation-dominated brane: the `+1` side is overdamped (Chat 9: `m² = 2H²` at the end point) and relaxes
to de Sitter with a redshifting radiation-like admixture of at most a few per cent; the throat side collapses. A homogeneous
classical evolution with one bulk scalar simply has no channel for producing a thermal bath; the only "radiation" it can
produce is the bulk-Weyl/`a^{−4}` term, and that is small and decaying here. The conclusion is justified *for that system*.

What is outside the runs and could change the answer for the model (not for these runs):
1. **Inhomogeneity.** `μ² = −7.7` means every mode with `k/a ≲ 2.8H` grows and the seed sign is random per Hubble patch
   (Chat 9's own caveat). The real outcome is a domain structure of `+`-patches and crunching `−`-patches separated by
   walls in `φ_b`; wall collisions, gravitational and scalar radiation from the walls, and the fate of crunching regions
   embedded in de Sitter (black holes on the brane? bulk singularities eating the brane?) are all unaddressed and are
   precisely where a "hot" component could come from. The present code (one bulk coordinate, flat homogeneous 3-space)
   cannot address any of it.
2. **Brane matter.** No brane-localised fields exist in the model, so nothing can be reheated. Gravitational particle
   production during the ~1-Hubble-time roll (Chat 9: `ρ_r ~ 10⁻³–10⁻² H⁴`) is the only generic channel and is negligible.
3. **The potential / tension function.** The result depends on the cubic `W` and the linear detuning `1 + cφ`; a tension
   function with a minimum of `V_E` at positive value, or a quadratic detuning, changes the fate (Chat 9 §5 already lists
   this). That is a different model.
4. **Quantum effects.** Seed statistics, tunnelling out of the `+1` de Sitter state, and the quantum fate of the crunch
   (Hertog–Horowitz hep-th/0406134 argue a crunch of this AdS-type is a genuine end in the dual description) are open.
5. **`t = 10⁻³`.** Not run. The linear theory gives `μ²(t)` almost flat and the two detunings agree, so the qualitative
   fates are expected to persist; but "expected" is the right word, and the boundary closure's failure there (growth rates
   4–12 instead of 1.657) is a numerical problem that has not been diagnosed, only sidestepped.
6. **KK / bulk-wave content.** The 1+1 setting includes only the homogeneous sector of bulk gravitons; brane-frame "dark
   radiation" from inhomogeneous bulk modes is excluded by construction.

## 5. Literature

* **BraneCode** — J. Martin, G. N. Felder, A. V. Frolov, M. Peloso, L. Kofman, "Braneworld dynamics with the BraneCode",
  hep-th/0309001, Phys. Rev. D 69, 084017 (2004). Abstract (verified today): static warped configurations with de Sitter
  branes are often unstable because of a tachyonic radion mass during inflation; unstable systems restructure toward
  lower-curvature configurations (flattened branes, lower effective 4D cosmological constant); brane collisions are
  common; with a bulk scalar, colliding branes approach a universal homogeneous anisotropic Kasner-like strong-gravity
  asymptotic, with equal Kasner indices along the brane and a different one in the extra dimension. The document's
  one-line summary is accurate. The present setting (single Z₂ shell, regular cone/horizon, no second brane) is different;
  the "lower-curvature relaxation" is the analogue of the `+1` branch, and the Kasner-like strong-gravity regime is the
  relevant comparison for the throat-side end, which is a bulk curvature singularity at the brane, not a collision.
* **Negative-potential collapse** — G. Felder, A. Frolov, L. Kofman, A. Linde, "Cosmology with negative potentials",
  hep-th/0202017, Phys. Rev. D 66, 023507 (2002): a universe whose scalar reaches a negative-energy region does not enter
  AdS but reaches a turning point where the energy density vanishes and then contracts to a singularity "practically
  independent of V(φ)". This is the correct reference for the throat side (Chat 9 already cited it).
* **AdS crunches** — T. Hertog, G. T. Horowitz, "Towards a big crunch dual", hep-th/0406134, JHEP 0407:073 (2004) (smooth
  AdS-invariant initial data evolving to a crunch; scalar solitons imply negative-energy configurations). Relevant as lore
  for what a negative-tension/negative-energy brane state means, not as a direct comparison.
* **Single moving dilatonic wall** — H. A. Chamblin, H. S. Reall, "Dynamic dilatonic domain walls", hep-th/9903225,
  Nucl. Phys. B 562, 133 (1999): single wall in a bulk with a scalar (Liouville potentials), static bulk, singular bulks with
  horizons, time-dependent wall motion and worldvolume inflation. Closest exact-solution precedent for a single shell
  moving in a scalar bulk.
* **Brane Friedmann equation with a bulk scalar** — K. Maeda, D. Wands, "Dilaton-gravity on the brane", hep-th/0008188,
  Phys. Rev. D 62, 124009 (2000): the projected equations used implicitly in §2–3 (the `nn` Gauss relation and the
  non-closure through the bulk Weyl term).
* I found no paper that follows a single dS brane numerically through a roll-off to a crunch in a 5D bulk with a scalar;
  the two-brane BraneCode runs and the 4D negative-potential literature are the nearest. That is a fair basis for the
  document's "no novelty claimed / the specific setting is the project's own".

## 6. Line-by-line: overclaims and imprecisions in `00_READ_FIRST.md`

| Where | Statement | Problem | Fix |
|---|---|---|---|
| title, §0, §2, §5 | "settles into another empty de Sitter universe … stays there" | observed for 0.32 e-folds / 0.49 Hubble times; `H_J` still falling; chart freezes at `H₀τ = 7.39` | "relaxes toward a slower de Sitter brane (followed to within 1.3 % of the RS value 0.638 H₀ before the chart freezes at `H₀τ ≈ 7.4`)" |
| §0, §2, §5 | "expanding at 0.65 of its original rate", "`H = 0.65 H₀`" | last value 0.646, still decreasing; RS attractor 0.638 | "≈ 0.64 H₀, consistent with the exact de Sitter brane in AdS₊ (0.638 H₀)" |
| §2 | "the simple de Sitter-brane-in-AdS₊ estimate at this `t` gives 0.603; the 5D value is 0.646" | **[ERR]** estimate omits the `t²` term of `σ²`; exact value 0.638 | replace 0.603 by 0.638 (and 0.605 by 0.616 at `t = 0.03`) and drop the implied 7 % discrepancy |
| §2 | "Chat 9's four-dimensional effective theory predicted 0.606 for `t → 0`" | correct, but it is the same RS formula at `O(t)`; say so | "Chat 9's EFT value 0.606 is the `t → 0` limit of the same formula" |
| §2 | "the brane's proper time saturates near `H₀τ ≈ 8.3`" | **[ERR]** it saturates at 7.39 (`t = 0.1`), 7.23 (`t = 0.03`) | correct the number; add that this bounds the observed window to ~0.5 Hubble times of the new state |
| §2 | "The `t = 0.03` run behaves identically (settling at `φ_b ≈ 0.995`, `H_J/H₀ ≈ 0.64`)" | run was stopped by the `φ_b > 0.995` guard while still rolling (0.22 e-folds after `φ_b > 0.98`) | "the `t = 0.03` run follows the same path and was stopped at `φ_b = 0.995`, `H_J = 0.637 H₀`, still relaxing (RS value 0.616)" |
| §0, §2, §5 | "collapses to a big crunch at a finite proper time `H₀τ* = 6.43`" / "genuine big crunch: `H_J ≈ −0.65/(τ*−τ)`" | power law not converged (`α` = 0.47–0.70 depending on window and resolution; local exponent 1.2–2.6); late stage differs by 2–3× between resolutions; `τ*` is a short extrapolation | "the expansion reverses and the brane contracts by ≥ 2 e-folds with rapidly growing curvature; the classical evolution ends in a singularity within `ΔH₀τ ≈ 0.44` (extrapolated `H₀τ* ≈ 6.43`), whose asymptotic form is not resolved" |
| §0, §2, §5 | "the scalar runs to `−∞` because the registered potential is unbounded below", "nothing stops it" | **[ERR]** mechanism: turnaround occurs at `U = +3.9`, `W > 0`; cause is the sub-BPS tension `t(1+cφ_b) < 0` beyond `φ_b = −1.673`; runaway is kination in the contracting phase (FFKL); the EFT potential has a negative *minimum* | "beyond `φ_b = −1/c = −1.67` the tension falls below its BPS value, the brane's effective energy turns negative, and the brane contracts; in the contraction the scalar's kinetic energy dominates and `φ_b` runs away (`dφ_b/d ln a ≈ −2`), as in any negative-energy collapse" |
| §1 table | "The linear regime is reproduced to five digits" | control run 1.62694 vs 1.62702 (5e-5 relative, four–five digits); seeded runs 1.625–1.628 (1e-3) | "to `1e-3` (seeded runs) and `5e-5` (zero-seed control)" |
| §1 table | "Momentum constraint … `≲2×10⁻²` in the final collapse" | 1.9 % at `dz = 1e-3`, 7.8 % at `2e-3` | quote both |
| §1 | "the non-linear results move by under 1 % between resolutions" | true for turnaround, end `φ_b`, `τ*`; false for the late collapse (`H_J` at equal `τ`: −29 vs −84 `H₀`) | "… under 1 % for the turnaround, the end-state `φ_b` and the extrapolated `τ*`; the last ~0.03 Hubble times of the collapse are not converged" |
| §0 | "computed without approximation beyond finite resolution" | also: homogeneous, classical, one bulk scalar, no brane matter, seed of a specific form and size, fixed `t = 0.1` | list them |
| §0, §5 | "No branch … now holds at the full nonlinear 5D level" | correct for the homogeneous classical system; the model-level statement needs the §4 qualifiers | "… for the homogeneous classical 5D evolution; inhomogeneous fragmentation and brane matter are not covered" |
| §3 | "found relaxation to a lower-Hubble state or a brane collision with Kasner-like asymptotics" | accurate | — |
| §4 | "second-order boundary closure … corrupts the unstable mode there" | supported by `formulation_lab*.json` (spurious modes with growth 4–12 are localised at the shell and violate the constraint by factors 10²–10⁵), but no higher-order closure was tried, so "second-order" is a hypothesis | "spurious constraint-violating modes localised at the shell appear at `t = 10⁻³` at the resolutions tried; presumably the boundary closure; not yet fixed" |
| §1 | "the shell (kept at `z = 0` …) supplies Neumann data" | fine; add that pinning the shell exhausts the residual conformal freedom up to `t → f(t)`, so `H_J` and `τ` are gauge-invariant | optional |

## 7. Independent check runs (this folder, `dz = 2e-3`)

Three runs of the unmodified `rolloff5d.py`, all `t = 0.1`, `dz = 2e-3` (logs and `*_timeseries.npz` in `runs/`; total
≈ 12 min):

| run | change vs. the folder's run | result |
|---|---|---|
| `plus_L16` | far boundary at `L = 16` instead of 12, `tf = 12` | identical to `t01_plus_coarse` to all printed digits (`φ_b = 0.99038`, `H_J/H₀ = 0.63923`, `H₀τ = 7.3810` at `t = 12`). The far boundary plays no role; the late `dz`-dependence of `H_J` (0.639 coarse vs 0.646 fine at `t = 12`) is truncation error in the degenerating chart, not boundary contamination. |
| `plus_dc1e-3` | seed ten times larger (`Δc = 10⁻³`) | same end state: `φ_b = 0.9905`, `H_J/H₀ = 0.643` at `H₀τ = 5.96`, i.e. the whole history shifted earlier by `ln 10 / 1.627 = 1.41` Hubble times as linear theory predicts. The `+1` attractor does not depend on the seed amplitude. |
| `minus_ko0.1` | Kreiss–Oliger coefficient doubled (0.10 instead of 0.05), `tf = 9` | turnaround `H₀τ = 5.9909`, `φ_b = −1.9562` (identical); end `H_J/H₀ = −83.3`, `φ_b = −7.188` vs `−83.6`, `−7.190` for the folder's coarse run. The late-stage collapse is insensitive to dissipation, so the factor-2 difference between `dz = 2e-3` and `1e-3` in that stage is genuine under-resolution, not a dissipation artefact. |

These support the document's turnaround and end-state numbers and sharpen the two caveats above (chart freeze; unresolved
late collapse). I did not attempt `t = 10⁻³`.

## 8. Safe abstract-style paragraph

> We follow the roll-off of the unstable de Sitter shell of the detuned-BPS model numerically in the full five-dimensional
> Einstein–scalar theory, in a 1+1 conformal chart with the Z₂ shell as a Neumann boundary, for a homogeneous brane and a
> single bulk scalar, at tension detunings `t = 0.1` and `0.03` (the registered `10⁻³` was not run). The linear growth rate
> agrees with the 4D perturbation theory to `10⁻³`. Seeded toward `φ = +1`, the brane's expansion rate falls monotonically
> and, when the chart freezes after about 7.4 Hubble times of brane proper time, is within 1.3 % of the exact
> Randall–Sundrum de Sitter brane in AdS₊ (`H = 0.638 H₀` at `t = 0.1`), still decreasing, with no oscillation and a
> redshifting radiation-like admixture of at most a few per cent. Seeded toward the throat, the brane crosses `φ_b = −1`
> with `H = 0.79 H₀`, its tension falls below the BPS value beyond `φ_b = −1.67`, its expansion reverses at
> `φ_b = −1.96` (4D EFT: `−1.94`), and it contracts by at least two e-folds with the scalar running away kinetically and
> the five-dimensional curvature at the brane growing by more than three orders of magnitude; the classical evolution ends
> in a singularity about 0.44 Hubble times after the turnaround, whose asymptotic form our resolutions do not fix. Within
> the homogeneous classical evolution neither branch produces a radiation-dominated brane. Inhomogeneous fragmentation of
> the shell, brane-localised matter and quantum effects are outside this calculation.

---

Files in this folder: `REFEREE_PHYSICS.md` (this), `check_numbers.py` (all numbers in §1–3 and §7 from `../runs/*.npz`),
`runs/` (three check runs and their logs).
