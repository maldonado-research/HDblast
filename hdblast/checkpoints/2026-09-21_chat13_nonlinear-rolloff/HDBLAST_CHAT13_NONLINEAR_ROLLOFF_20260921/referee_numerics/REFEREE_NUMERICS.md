# Referee report (numerics) on HDBLAST Chat 13 — non-linear 5D roll-off

Referee: numerics referee (Claude, adversarial role). Date: 21 September 2026.
Scope: `rolloff5d.py`, `derive_evolution_equations.py`, `formulation_lab.py`, `analyse_runs.py`, `runs/`, and the claims in `00_READ_FIRST.md`.
Everything I wrote is under `referee_numerics/` (this file, `referee_derivation.py`, `rolloff5d_ref.py`, `analyse_referee_runs.py`,
`REFEREE_ANALYSIS.json`, `derivation_output.txt`, `runs/`). Nothing else was touched. Compute used: about 25 minutes wall clock,
runs in parallel, all with the system `python3` (numpy 2.0.2); SymPy only through the Chat 11 venv.

## 0. Verdict in one paragraph

The equations, the junction conditions and the deviation formulation are correct (independently re-derived, Section 1); the code
reproduces bit-for-bit, its linear regime matches the Chat 9 growth rate to 1e-3, and the results through the roll-off are robust to
seed amplitude, Kreiss–Oliger strength, boundary-closure variant, domain size and resolution (Section 3). **Two of the headline
numbers are, however, overstated.** (i) The `φ = +1` branch does **not** "settle at `H_J/H₀ = 0.646`": the conformal chart's lapse
collapses (`e^{B_b} ∝ e^{−t}`), the brane's proper time saturates at `H₀τ ≈ 7.38`, the domain wall is still receding into the
bulk at the speed of light, and at the last numerically clean time (`t ≈ 10`, `H₀τ ≈ 7.33`) the brane Hubble rate is
`0.649 H₀` **and still falling at `dH/dτ ≈ −0.035 H₀²`**. The asymptotic value is not determined by these runs; `0.65` is an upper
bound and the Chat 9 / Randall–Sundrum estimate `0.60` is entirely compatible with the trend. The stated "plateau 0.639–0.646" is
the interval over which the code's late-time output (`t > 10`) is being corrupted by the collapsing lapse (the constraint at the
shell reaches 3 % at `t = 10`, 18 % at `t = 12`, 100 % at `t = 15`; the author's monitor excludes the six shell points and reports
≲1e-3). (ii) The throat branch: the turnaround (`H₀τ = 5.991`, `φ_b = −1.956`) is solid; the runaway of `φ_b` and `H_J` is
solid up to `H₀τ ≈ 6.30` (`H_J ≈ −5 H₀`, two resolutions agree to 1–2 %); beyond that the lapse collapses again, the two resolutions
differ by a factor 2 in `H_J` at equal `t`, and `H_J ≈ −α/(τ*−τ)` is an extrapolation whose coefficient `α` is not converged
(0.48–0.70 depending on threshold and resolution). `τ*` is bracketed in `[6.42, 6.45]` only because every run reaches
`H₀τ = 6.415` with `H_J < −40 H₀` and still accelerating. "Big crunch" is the physically natural reading; it is not what the code
has computed. The qualitative conclusion — no branch produces a radiation-dominated, decelerating universe — survives; the
sentence "now holds at the full nonlinear 5D level" should carry the caveats above and the `t = 10⁻³` caveat, whose stated cause
(BPS cancellation in the boundary closure) my tests do not support (Section 4).

## 1. Equations and junction conditions (independent re-derivation)

`referee_derivation.py` (SymPy, Chat 11 venv; output in `derivation_output.txt`) goes through the Einstein tensor and the scalar
stress tensor, `G_AB − T_AB` with `κ₅² = 1`, instead of the trace-reversed form:

* `G_AB − T_AB` equals the trace-reverse of `R_AB − ∂_Aφ∂_Bφ − (2/3)U g_AB` in 5D (so `U = V`, no factor lost).
* `G_tz − T_tz` is **exactly** the code's momentum constraint `M` (ratio 1). `G_tt − T_tt` is exactly `½` of the code's `Hc`.
* Substituting the three closed-form evolution equations (`A_tt`, `B_tt`, `φ_tt` as in the code header) into **all** components of
  `G_AB − T_AB` leaves only the `(tt)`, `(tz)`, `(zz)` components, which are `Hc/2`, `M`, `Hc/2`: the evolution system is complete
  and the constraints are the only remaining content. Klein–Gordon residual 0.
* Junction: for `z = const` with unit normal `e^{−B}∂_z` (out of the bulk `z < 0`), `K^t_t = e^{−B}B_z`, `K^x_x = e^{−B}A_z`.
  The sign is anchored on Randall–Sundrum rather than on a memorised convention: kept side `y > 0` with `e^{−2k|y|}`, outward normal
  `−∂_y`, gives `K^μ_ν = +k δ^μ_ν = (σ/6)δ^μ_ν` for positive tension `σ = 6k`. So on the kept side, with outward normal,
  `K^μ_ν = (σ/6) δ^μ_ν`, i.e. `A_z = B_z = (σ_t/6)e^B`, as in the code and in `solve_shell()` (`ρ_y/ρ = σ/6`, warp factor maximal at
  the brane). Scalar: `□φ = U' + σ_t'(φ)δ`, Z₂ ⇒ `φ_y(0⁻) = −σ_t'/2` ⇒ `φ_z = −(σ_t'/2)e^B`. Matches the code.
* **Consistency identity** (new check): with `A_z = B_z = σe^B/6` and `φ_z = −σ'e^B/2` imposed for all `t`, the momentum constraint
  at the shell vanishes identically, `M|_shell = −3∂_t(σe^B/6) + 3(σe^B/6)B_t + φ_t σ'e^B/2 = 0`. This fixes the relative sign of
  the gravitational and scalar junction conditions independently of the static solution, and it shows that the Neumann data are
  compatible with the constraint (the pair is fixed up to an overall sign, which the RS check and the static background fix).
  It also gives `A_tz|_shell = ∂_t(σe^B/6) ≠ 0`: the code's constraint monitor uses `pa_z = 0` there, which is why it must skip the
  last points; my copy evaluates `pa_z` one-sidedly and monitors the shell layer separately.

Re-running the author's `derive_evolution_equations.py` (from my folder, so its JSON lands here) reproduces all six `True`.

**Deviation formulation** (read line by line): with `A_s = ln ρ + t`, `B_s = ln ρ`, the static parts subtracted in `ra`, `rbb`, `rf`
are exactly `A_s,zz − 3 + 3H_c² + (2/3)ρ²U_s`, `B_s,zz + 3 − 3H_c² + φ_s,z²/2 − ρ²U_s/3`, `φ_s,zz + 3H_cφ_s,z − ρ²U_s'`, and
`SU = e^{2B}U(φ) − ρ²U(φ_s) = ρ²[(e^{2b}−1)U(φ) + (U(φ)−U(φ_s))]` is what `dU_exact` returns (I checked the sextic's coefficients
`_Uc` by hand: `−1/6, 4/3, −5/3, −4/9, 17/18, 0, −2/27`). Consistent in all three equations. `H_J = A_t e^{−B}` and
`H₀τ = ∫e^{b_b}dt` (trapezoidal) are right. The ghost closure `f_{+j} = f_{−j} + 2jh g` gives `f_z(0) = g` exactly with the
4th-order stencil and an `O(h f''')` error in `f_zz` at the last point (so the scheme is globally second-order at best; the
measured constraint violation indeed scales as `dz²`: the author's `dz = 2e-3` and `1e-3` runs have `M` in the ratio 4.0–4.1 at
every time). The KO operator sign is now correct (damping). One inconsistency: the KO ghosts of the time-derivative fields
(`pa, pbv, pf`) use slope 0, while the true slope is `∂_t g`; this is an `O(ε·∂_t g)` forcing on the last three points (an
`O(h)` boundary error). I tested the consistent slope (`--pslope 1` in my copy): it changes `H_J` by `3e-7` and `φ_b` by `3e-8`
through `t = 12` — negligible, and it does **not** remove the late-time shell-layer constraint growth (Section 3.3), so that
growth is not caused by this inconsistency.

## 2. Code reproducibility

`ref_control` (`t = 0.1`, `dc = 0`, `dz = 2e-3`, `L = 8`, `t_f = 8`): `φ_b(8)` identical to the author's `t01_control` to all
printed digits (`2.335299995e-3`), round-off growth rate `1.62694` (linear theory `1.62702`). My additional monitors on the static
solution: shell-layer momentum constraint `≤ 9e-12`, Hamiltonian `≤ 4e-10` (interior) and `1.4e-11` (shell). The static solution
is static to round-off, as claimed.

## 3. Attempts to refute the numerical results

All referee runs at `dz = 2e-3` unless stated; comparisons are against the author's `dz = 2e-3` runs at equal conformal time.
Numbers in `REFEREE_ANALYSIS.json`.

### 3.1 Seed sign flipped, half amplitude (`dc = −5e-5`)
Expected time shift `ln 2 / s = 0.4260`. Measured: `φ_b` crosses `−0.5, −1, −1.5` later by `0.420, 0.430, 0.430`; turnaround later
by `0.430` (`t = 6.43`, `H₀τ = 6.420`, `φ_b = −1.971` vs `−1.956`). In the `(φ_b, H_J)` plane the trajectories agree to 1–2 %
(`H_J` at `φ_b = −1.5`: 0.438 vs 0.448; at `−3`: −1.673 vs −1.640; at `−5`: −7.64 vs −7.70). The evolution is the growth of the
single unstable mode followed by a seed-independent non-linear trajectory. **Passed.**

### 3.2 Kreiss–Oliger strength (`ε = 0.02, 0.05, 0.10`), `φ = +1` branch
Through `t = 10`: `|Δφ_b| ≤ 3e-6`, `|ΔH_J| ≤ 3e-5` between the three values. Same for the throat branch with `ε = 0.10`
(`H_J(t=9) = −83.3` vs `−83.6`). End states do not depend on dissipation. **Passed.** But the shell-layer momentum constraint at
`t = 12` scales with `ε`: `0.08, 0.18, 0.30` — the late-time boundary error is dissipation-driven (see 3.3).

### 3.3 Domain size `L = 16` (`t_f = 16`) — far boundary and the late-time drift
Bit-identical to `L = 12` through `t = 12` (`|ΔH_J| ≤ 4e-8`): the far boundary is causally irrelevant, as claimed. **But** the
`L = 16` run shows the same late drop as the author's discarded `t01_plus` run (`H_J/H₀`: 0.639 at `t = 12`, 0.625 at 13, 0.586 at 14,
0.480 at 15, 0.192 at 16), while `H₀τ` moves only from 7.381 to 7.388. So the drop is **not** a far-boundary reflection. It is
accompanied by the shell-layer constraint going to `O(1)`:

| `t` | 8 | 10 | 11 | 12 | 13 | 14 | 15 |
|---|---|---|---|---|---|---|---|
| `M` shell layer (last 6 points, one-sided `pa_z`) | 5e-3 | 2.9e-2 | 7.5e-2 | 0.18 | 0.37 | 0.62 | 0.81 |
| `M` interior (author's slice) | 5e-4 | 5e-4 | 1.4e-3 | 3.7e-3 | 1e-2 | 2.8e-2 | 8e-2 |
| Hamiltonian, shell layer | 3.5e-3 | 7.8e-3 | 1.1e-2 | 1.6e-2 | 3.2e-2 | 7.5e-2 | 0.20 |
| `1 + A_t` at shell | 0.19 | 0.032 | 0.012 | 0.0044 | — | — | — |

`H_J = (1 + pa)e^{−b}/ρ_b` is the ratio of two quantities that decay like `e^{−t}`; a fixed absolute error becomes an exponentially
growing relative one. The resolution difference in `H_J` between the author's two runs grows the same way: `1.9e-4` (t=8), `4.5e-4`
(9), `1.0e-3` (10), `2.5e-3` (11), `6.6e-3` (12). **Conclusion: the code's output at the shell is trustworthy to ≲0.3 % up to
`t ≈ 10` (`H₀τ ≈ 7.33`) and not beyond.** The "plateau" values quoted at `t = 12` are already 1 % contaminated at `dz = 2e-3`.

### 3.4 What the `φ = +1` branch actually shows
From the author's `dz = 1e-3` run (`t01_plus_v2`), `H_J/H₀` and its proper-time derivative:

| `t` | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|
| `H₀τ` | 7.055 | 7.254 | 7.335 | 7.366 | 7.378 |
| `H_J/H₀` (dz=1e-3 / 2e-3) | 0.6611 / 0.6609 | 0.6517 / 0.6513 | 0.6489 / 0.6478 | 0.6475 / 0.6449 | 0.6458 / 0.6392 |
| `dH_J/dτ / H₀²` | −0.059 | −0.038 | −0.036 | −0.07 (*) | −0.27 (*) |
| `φ̇_b/H₀` | +0.027 | +0.008 | +0.004 | +0.001 | −0.006 (*) |

(*) already in the corrupted window. The bulk snapshots (`snap_f_8/10/12`) show the domain wall (`φ` from −1 to +1) detached from
the brane and moving toward the throat at `dz/dt ≈ 1.1` (at `z ≈ −1.3, −3.3, −5.5`); between it and the brane `φ ≈ 0.99`. So the
"new de Sitter brane" is a brane in a growing AdS₊ region bounded by a light-like receding wall — a sensible picture (the wall
falls through the bulk horizon of the new state) — but the brane's `H` is still relaxing when the chart freezes. The number `0.646`
is `H_J` at `H₀τ ≈ 7.35`, about 1.3 Hubble times after the shell reached `φ_b = 0.9`, not an asymptote. The Chat 9 EFT/RS value at
this detuning is `0.603`; the gap `0.045` is about one Hubble time of the observed `dH/dτ`. Neither value is established; the run
gives `H_∞/H₀ ≤ 0.65`. The claim "`H = 0.65 H₀`" should become "`H_J` had fallen to `0.65 H₀` and was still decreasing when the
chart froze; consistent with relaxation toward the RS value 0.60".

Also: the `t = 0.03` run (`t003_plus`) did not "settle at `φ_b ≈ 0.995`" — it was **stopped by the `--phimax 0.995` criterion**
(`stop_reason: "phi_b left [-1.6,0.995]"`, `φ_b,max = 0.9950030`) while `φ_b` was still rising (`φ̇_b = +0.016 H₀`) and `H_J` was
still falling at `−0.058 H₀²` (`0.637` at the stop, `H₀τ = 7.08`). The sentence in `00_READ_FIRST.md` §2 should be corrected.

### 3.5 Throat branch and the crunch fit
Turnaround: `H₀τ = 5.9909` (dz=2e-3) / `5.9912` (1e-3) / `5.9909` (ε=0.1) / `5.9909` (pslope), `φ_b = −1.956` in all; half seed
gives `φ_b = −1.971` (1 % seed dependence). **Solid.** Comparison of the two resolutions at equal `t`: `H_J` agrees to 1 % up
to `t = 6.5` (`H₀τ = 6.264`, `H_J = −4.0`), 6 % at `t = 7` (`−8.6` vs `−9.2`), 40 % at `t = 8` (`−21` vs `−29.5`), factor 2 at `t = 9`
(`−40` vs `−84`), while `H₀τ` only moves from 6.36 to 6.42. This is the lapse collapse again (`b_b` from −0.2 at turnaround to −5.3
at `t = 9`). The fit `−1/H_J` linear in `τ` gives (fine / coarse): threshold `H_J < −3`: `τ* = 6.4396 / 6.4260`, `α = 0.70 / 0.65`;
`< −10`: `6.4364 / 6.4231`, `α = 0.63 / 0.57`; `< −20`: `6.4324 / 6.4210`, `α = 0.52 / 0.48`. Free-exponent fits give `p = 1.01–1.07`.
So: `α` is not a physical number (it drifts with threshold and resolution and is measured entirely inside the non-converged window);
`τ*` is robust only in the weak sense that all runs pass `H₀τ = 6.415` with `H_J < −40` and `|H_J|` accelerating, so
if the trend continues the crunch is within `0.03` of `6.42`. The brane scale factor from my `a_b` record: `A_b − ln ρ_b = 5.82,
5.42, 4.85, 4.14, 3.65` at `t = 6, 6.5, 7, 8, 9` — decreasing, but with `A_t = 1 + pa` shrinking (`−0.91, −0.61, −0.44, −0.32, −0.23`),
i.e. in this chart `ln a_b` has not reached `−∞`; the chart runs out (constant-`t` slices pile up at `τ ≈ 6.42`) before the code
sees `a_b → 0`. Physically the collapse with `φ̇_b ≈ −66 H₀` at `t = 9` and `U → −∞` leaves little room for anything but a crunch, and
I have no counter-scenario to offer; but the statement should be "the brane collapses with `H_J ≤ −40 H₀` by `H₀τ = 6.42`; a
`1/(τ*−τ)` extrapolation places the crunch at `H₀τ* ≈ 6.43`; the code cannot follow it further because the lapse collapses", and
the "`H_J ≈ −0.65/(τ*−τ)`" coefficient should not be quoted as a result.

## 4. The failure at `t = 10⁻³`

The author attributes it to a second-order ghost closure spoiled by the BPS cancellation at the shell. I ran the seeded problem
(`dc = 1e-4`, `L = 3`, `t_f = 3`) at three resolutions and a zero-seed control:

| `dz` | 1e-3 | 5e-4 | 2.5e-4 |
|---|---|---|---|
| growth rate of `φ_b` (windows 0.5–3) | ≈ 4–6, blow-up (NaN) at `t = 1.47` | ≈ 4.0, 2.9, 3.8; blow-up at `t = 2.61` | 2.04, 2.32, 2.05, 2.15, 2.29 |
| linear theory | 1.657 | | |
| `M` shell layer at `t = 1` | 3e-4 | 7e-5 | 4e-6 |
| `M` interior (author's slice) | 0.8 | 0.6 | 0.4–0.5 |

Points across the domain wall (proper width ≈ 1, i.e. `Δz ≈ 1/ρ_b = 0.0127` at `t = 10⁻³` versus `0.128` at `t = 0.1`): 6, 13, 25, 51
for `dz = 2e-3, 1e-3, 5e-4, 2.5e-4`; the `t = 0.1` production runs have 64–128. In the wall the linearised scalar equation balances
`f_zz ~ f/(0.0127)² ≈ 6·10³ f` against `−ρ_b²U''(φ_s) f ≈ +2·10⁴ f` (tachyonic in the wall core); under-resolving that balance
produces spurious growing modes with rates of order the ones seen. **The spectrum of measured rates converging downward with
resolution (5 → 4 → 2.2) is what a resolution problem looks like, not a boundary-closure defect**: at `dz = 2.5e-4` the shell layer
and the wall are constraint-clean (`M ≤ 4e-3` on the last 6 points, `≤ 1e-6` in the wall). What is *not* clean at `t = 10⁻³` is the
outgoing front deep in the AdS₋/cone region: at `t = 2.5` the normalised momentum constraint is `0.54` at `z ≈ −2.5` (`ρ = 0.29`),
with `|M| ≈ 7` carried by `∂_z A_t` across a front whose amplitude (`~10⁻³`) already exceeds the seed by 100. The residual 30 %
excess of the growth rate at `dz = 2.5e-4` and this front are the actual open problem; whether they are one problem, I could not
settle in the budget. (A zero-seed control at `dz = 2.5e-4`, `t_f = 2`, turned out to be uninformative: with `dc = 0` the two
shell solves returned identical parameters, the deviation fields were exactly `0.0` and stayed `0.0` — the deviation formulation
has no round-off seed of its own; the author's `t = 0.1` control was seeded by a `1e-14` difference between the two solves.) The author's own
`diag_dz0.001` (zero seed, `t_det = 10⁻³`) reached a normalised `M = 3.7` within `t = 1` from round-off while `diag_dz0.002` stayed
at `1.5e-4`, i.e. a numerical instability with rate `> 30` that appears at *finer* resolution — that pattern also points to the
under-resolved wall/stiff-coefficient regime (`ρ_b² ≈ 6200` multiplies every potential term) rather than to a boundary closure.

Suggested fixes, cheapest first: (a) a non-uniform grid or a coordinate `ζ` with `dz = ρ⁻¹ dζ`-type stretching so that the wall
gets ≥ 60 points while the cone region is not over-resolved (a uniform `dz = 2.5e-4` on `L = 12` costs ~50k points × 100k steps,
too slow in pure numpy but feasible with a stretched grid of ~8k points); (b) evolve in the proper-distance-like coordinate `y` of
the static solution (gauge `ds² = e^{2B}(−dt² ) + dy²...` is not conformal, but a 1+1 code with lapse and non-unit shift is standard);
(c) if the shell closure is to be improved anyway: 4th-order one-sided (5-point) stencils for `f_z` and `f_zz` at the last two
points with the Neumann value substituted, which removes the `O(h f''')` term; (d) a lapse condition (e.g. harmonic slicing or a
periodic re-slicing to the static gauge of the current configuration) to escape the lapse collapse that ends every run at
`H₀τ ≈ 7.4` (plus branch) and `6.42` (throat branch) — this is the single change that would turn "still falling at 0.65" into a
number. Imposing the junction through the Hamiltonian constraint (`formulation_lab.py` "hamB") is not needed for correctness
(Section 1 shows the Neumann data are constraint-compatible) and I would not pursue it.

## 5. Smaller points

* `analyse_runs.py` computes the growth-rate fit on `5·dφ(0) < dφ < 5e-3`, `t > 0.5`; the summary fit in `rolloff5d.py` uses
  `30·seed < dφ < 2e-3`. Both give 1.625–1.628 at `t = 0.1`; fine. The `t01_plus` entry in `ANALYSIS.json` ("plateau 0.534–0.645")
  is the corrupted `t_f = 16` run and should be deleted from the table or labelled.
* The constraint monitor slice `10:-6` excludes exactly the layer where the late-time violation lives; the "≲1e-4 through the
  roll" statement is correct for the interior, but the shell layer should be reported too (it is 3e-2 at `t = 10`, `dz = 2e-3`).
* `stop_phi = (−1.6, 0.995)` defaults: the upper stop truncated the `t = 0.03` run (Section 3.4); use `--phimax 1.5` for plus runs.
* `formulation_lab.py` has the leftover `damp*M/3.0*0` placeholder; harmless.
* Prior art: hep-th/0309001 (Martin, Felder, Frolov, Peloso, Kofman, "Braneworld dynamics with the BraneCode") — title and
  authors verified online today; their finding (unstable dS-brane warped geometries restructure violently toward lower-curvature
  configurations) matches the plus branch here qualitatively, as the author says.

## 6. Status table (referee's version)

| Claim | Status |
|---|---|
| Evolution equations, constraints, junction signs | **Established** (two independent SymPy routes; RS-anchored sign; constraint-compatibility identity) |
| Static solution static to round-off; linear rate 1.627 reproduced | **Established** (bit-for-bit reproduction) |
| Trajectory through the roll independent of seed, KO, closure, domain, resolution (to ≲1 %) | **Numerically supported** (this report) |
| `φ = +1` branch: brane in AdS₊ region, no oscillation, no radiation, `H_J` falls monotonically to `≤ 0.65 H₀` | **Numerically supported** up to `H₀τ ≈ 7.33` |
| "settles at `H_J = 0.646 H₀`" (and "0.64" at `t = 0.03`, "`φ_b ≈ 0.995`") | **Not established**: still relaxing (`dH/dτ ≈ −0.035 H₀²`) when the chart freezes; `t = 0.03` run was cut by the stop criterion |
| Throat branch: turnaround at `H₀τ = 5.99`, `φ_b = −1.96`; `H_J < −5 H₀` by `H₀τ = 6.30` | **Numerically supported** (two resolutions ≤ 2 %) |
| Crunch at `H₀τ* = 6.43`, `H_J ≈ −0.65/(τ*−τ)` | **Extrapolation**: `τ*` plausible within ±0.03; `α` not converged; code stops at `τ = 6.42` because the lapse collapses |
| `t = 10⁻³` failure = boundary closure / BPS cancellation | **Not supported**: rates converge downward with `dz` (5 → 4 → 2.2 vs 1.657); shell layer clean at `dz = 2.5e-4`; the remaining violation is a bulk front deep in the cone |
| "No branch gives a hot radiation universe, now at the non-linear 5D level" | **Supported for `t = 0.1, 0.03`, homogeneous, up to the times above**; the `t = 10⁻³` case remains uncomputed |

## Appendix — runs performed (all in `referee_numerics/runs/`)

`ref_control` (t=0.1, dc=0), `ref_minus_half` (dc=−5e-5), `ref_plus_ko002`, `ref_plus_ko01`, `ref_plus_L16` (t_f=16),
`ref_plus_pslope`, `ref_minus_ko01`, `ref_minus_pslope`, `ref_t1e3_dz1e-3`, `ref_t1e3_dz5e-4`, `ref_t1e3_dz5e-4_pslope`,
`ref_t1e3_dz2.5e-4`, `ref_t1e3_control_dz2.5e-4`. Record columns in `*_timeseries.npz` are the author's twelve plus
`[M_shell, H_interior, H_shell, a_b]`. The `*_checkpoint.npz` files (1 MB each) were deleted after use; the front diagnosis
in Section 4 was made from the `t = 2.5` checkpoint of `ref_t1e3_dz2.5e-4` with the constraint evaluated on the full grid
(code fragment reproducible from `rolloff5d_ref.py`'s `constraints()`). All floating point; nothing here is a certificate.
