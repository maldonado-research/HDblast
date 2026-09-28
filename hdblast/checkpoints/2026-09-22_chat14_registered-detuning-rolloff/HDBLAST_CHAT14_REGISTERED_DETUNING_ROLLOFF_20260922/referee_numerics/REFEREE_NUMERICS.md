# Referee report (numerics) on HDBLAST Chat 14 — roll-off at the registered detuning t = 10⁻³

Referee: numerics referee (Claude, adversarial role). Date: 22 September 2026.
Scope: `rolloff5d_v2.py` against `rolloff5d_v1.py`, `analyse_v2.py`, `make_figure_v2.py`, `runs/` (including the two runs that finished
during this review, `v2_t1e3_minus_c4_fine` and `v2_t1e3_plus_c4_finecoarse`), the claims in `00_READ_FIRST.md`, and the Chat 13
referee reports. Everything I wrote is under `referee_numerics/` (this file; `ref_closure_check.py` + `closure_check_output.txt`;
`ref_analyse.py` + `REFEREE_ANALYSIS.json` (author's runs) + `REFEREE_ANALYSIS_refruns.json` (mine) + `analysis_output_authors_runs.txt`;
`runs/`). Nothing else was touched. Compute: four single-threaded runs (`OMP_NUM_THREADS=1`), never more than two at a time,
about 31 minutes of run time in total.

## 0. Verdict in one paragraph

The stretched grid, chain rule, CFL condition, 4th-order Neumann closure, time-derivative slopes and taper are implemented correctly
(Section 1: the closure algebra is written out and checked to 10⁻¹⁵; a constant-coefficient eigenvalue test shows the closure is
stable). The linear growth rate is genuinely reproduced: 1.65714–1.65719 against 1.65719, seed-independent to 2·10⁻⁵ (my Δc = 3·10⁻⁸
run), and 75 points per wall is converged to 3·10⁻⁵ in the rate (my 40-point run gives 1.6564, an error 27× larger for a 1.9× coarser
grid — fourth-order behaviour). The symmetric closure at 222 points per wall (my run, 6.5 Hubble times) gives 1.6565 and is still creeping up toward 1.657
(1.90 at 75, 1.71 at 134 points), so "closure order, not only resolution" is supported (Section 5). The throat-side numbers reproduce from the saved data, with two small corrections of the
quoted values that are sampling artefacts (turnaround φ_b = −1.942, not −1.953; H_J = 0.798 H₀ at the φ = −1 crossing, not 0.793 —
both move toward the 4D prediction). **The +1-side plateau number is not established as stated.** The author's own `finecoarse` run
(coarse spacing 2·10⁻³ instead of 8·10⁻³, finished during this review) gives `H_J/H₀ = 0.639` at `H₀τ = 6.9` instead of 0.611, and a
plateau that flattens at ≈ 0.63 (`0.632` at τ = 7.0, `0.629` at τ = 7.05, slope −0.06 H₀²) instead of falling through 0.61 toward
0.6067. The two "resolutions" of the README (75 and 134 points per wall) share the same coarse grid, so their 0.1 % agreement did
not test this. The two grids agree to 5·10⁻⁵ at τ = 6.0 and 0.17 % at τ = 6.5, then separate: 1.5 % at 6.7, 3.8 % at 6.8, 4.5 % at
6.9, 16 % at 7.0. The default grid's bulk constraint monitor reaches 65 % before τ = 6.9 (a second light-like pulse launched from
the shell at t ≈ 7.0, not the seed pulse the README describes); on the refined coarse grid the same monitor stays below 1 %.
Consequently: (a) the clean window of the default grid ends at τ ≈ 6.6 (H_J ≈ 0.67), not 6.9; (b) the best available value at
the chart freeze is `H_J/H₀ ≈ 0.63`, still slowly falling, 3.5 % above the RS value, with an unknown residual coarse-grid error
(only two coarse resolutions exist, no Richardson estimate is possible); (c) the sentence "0.61, still falling, approaching 0.6067
from above" and the last row of the trend table should be replaced by "0.63–0.64 at the freeze (two coarse resolutions: 0.611 /
0.639), consistent with relaxation toward 0.6067 but not converged". The physical conclusion — empty de Sitter brane, no
radiation — survives; the quantitative plateau does not. Two further things the README should say: the Kreiss–Oliger
dissipation of v2 is effectively off (applied with `h = 1` in the index coordinate, 2500× weaker than v1; harmless as far as I can
tell, Section 1), and the Chat 13 referee proposed the closure fix that Chat 14 implemented (Section 5).

## 1. Code review: v2 against v1

**Coordinate map.** `build_grid` integrates `dz/dξ = −s(z)`, `s = dz_c + (dz_f − dz_c)·½(1 + tanh((z + z_f)/w))`, by RK4 at unit ξ
steps from the shell outward and reverses the node list; `z' = s(z_i)` and `z'' = s'(z_i)s(z_i)` are exact at the nodes
(`d²z/dξ² = (ds/dz)(dz/dξ)`). `dz1 = D1ξ/z'` and `dz2 = D2ξ/z'² − D1ξ z''/z'³` are the correct chain-rule forms. Checked on the actual
1392-node grid with `f = sin 3z`: max error 2·10⁻³ in `f_z`, 0.16 in `f_zz`, both in the transition zone `|z + 0.06| < 0.05`, where the
map stretches by **17.7 % per cell** (`dz_{i+1}/dz_i = 1.215` at `z = −0.061`). That is steep by the usual standard (≤ 5–10 % per cell):
consistent, but of reduced order there, and a plausible partial reflector. With `dz_c = 2·10⁻³` the maximum stretch is 4 % per cell.
The actual fine spacing for `--dzf 1.5e-4` is `s(0) = 1.694·10⁻⁴` (because `sig(0) = 0.9975`), hence 75 points per wall; "150" in the
README's closure-2 sentence should read 134 (`tiny_fine.log` header). Points across a wall-width feature: 75 at the shell, 1.6 in the
coarse region — but the receding wall broadens as `1/ρ(z)`, so at `z = −0.3` (ρ = 5.5) it spans ~25 coarse cells (Section 4).

**CFL.** `dt = 0.5·min(z')`: CFL 0.5 at the shell, smaller elsewhere; characteristic speed 1 in the conformal chart. Correct. Record
cadence and the trapezoidal `H₀τ = ∫e^{b_b}dt` are as in v1.

**Kreiss–Oliger.** `ko` applies `ε·δ⁶F/(64h)` with `h = 1` (index spacing); v1 used `h = dz`. The grid mode (`δ⁶F = −64F`) is damped at
`ε/h = 0.05` per unit `t` in v2 versus `ε/dz = 125` (dz = 4·10⁻⁴) in v1; per time step (`dt = 8.5·10⁻⁵`) that is `4·10⁻⁶` against v1's
`0.025`. **v2 runs with essentially no artificial dissipation.** Evidence that this is harmless: the constant-coefficient RK4 step map
has spectral radius ≤ 1 for both closures at ε = 0 (below); the end-state `f` profiles have a second-difference-to-field ratio of
3–6·10⁻⁵ in the fine region for all closure-4 runs (no grid noise; RK4 at CFL 0.5 itself damps the highest modes by ~1.4 % per step).
A stretched-grid KO would scale as `ε/z'_i` node-wise. This must be stated; the docstring's "Kreiss–Oliger dissipation" is not what
the v2 runs used, and the Chat 13 referee's "insensitive to KO" check does not carry over.

**4th-order closure — algebra.** Nodes `N−3 … N` (shell at `N`), ghosts `N+1, N+2`, `h = 1`, `g_ξ = g_z z'(0)`. Quartic extrapolation
through `N−3 … N+1` (vanishing fifth difference): `F_{N+2} = 5F_{N+1} − 10F_N + 10F_{N−1} − 5F_{N−2} + F_{N−3}` (code `e2` ✓). Centred
5-point `D1` at `N`: `(F_{N−2} − 8F_{N−1} + 8F_{N+1} − F_{N+2})/12 = g_ξ`; substituting, `3F_{N+1} + 6F_{N−2} − 18F_{N−1} + 10F_N − F_{N−3}
= 12g_ξ`, i.e. `F_{N+1} = 4g_ξ − 2F_{N−2} + 6F_{N−1} − (10/3)F_N + F_{N−3}/3` (code `e1` ✓). Equivalently the shell derivative is the
standard one-sided formula `f'_N = (−F_{N−3} + 6F_{N−2} − 18F_{N−1} + 10F_N + 3F_{N+1})/12`, verified to 10⁻¹⁵ (`closure_check_output.txt`, A).
Closure-4 ghosts reproduce polynomials of degree ≤ 4 to round-off; closure 2 fails at degree 3. Measured truncation orders
(`sin(2z + 0.3)`): closure 2 — `D1` exact at the shell, `O(h²)` at `N−1`, **`D2` first-order (`2hf'''/9`) at the shell and at `N−1`**; closure 4 —
`D1` `O(h⁴)` at `N−1`, `D2` `O(h³)` at both. So the README's "second-order at the boundary" for closure 2 is generous (its `D2` is
first-order) and "fourth-order" for closure 4 holds for `D1` (its `D2` is third-order). The eigenvalue errors behave accordingly:
closure 2 +0.24 (75 pts), +0.06 (134); closure 4 −8·10⁻⁴ (40), −3·10⁻⁵ (75), ≲10⁻⁵ (134).

**Stability.** For `u_tt = u_zz` (Dirichlet left, Neumann right, 60 nodes) the semi-discrete `D2` with either closure has purely real,
negative eigenvalues (largest `−6.87·10⁻⁴ = −(π/2/60.5)²`); the RK4 step map at CFL 0.5 has spectral radius `1.0000000000` (ε = 0) and
`0.99999999` (ε = 0.05) for both closures. No growing or complex modes in the constant-coefficient problem; nonlinear stability rests on
the 11–14-Hubble-time runs.

**Time-derivative slopes.** `neumann_t` gives `ġ_A = (σ'φ_t + σB_t)e^B/6`, `ġ_F = −(σ''φ_t + σ'B_t)e^B/2`, `σ'' = 4φ + t·d`: correct
(the static parts are time-independent). At the shell `M = −(σ'φ_t + σB_t)e^B/2 + σB_te^B/2 + σ'φ_te^B/2 = 0` identically, so the
shell-layer monitor is meaningful. Since the momenta's spatial derivatives enter only the (effectively absent) KO term, this change
affects the monitor, not the evolution.

**Taper.** `w = ½(1 − cos(π(z + L)/(L·taper)))` on `z < −L(1 − taper)`; `C¹`, `w → 0` at `−L`. It violates the bulk constraints by `O(Δc)`
in `z < −8.5`, causally disconnected from the shell for `t < 8.5`; all quoted +1-side numbers are at `t ≤ 8.4`, throat numbers at
`t ≤ 8.0`. Harmless for everything quoted; the tiny-seed and `ext` runs beyond `t = 8.5` are used only qualitatively.

**Static background on the stretched grid.** `hz = min(dz_f/4, 2.5·10⁻⁵)`, cubic-Hermite onto the nodes; junction residuals on the grid
`3.6·10⁻¹¹ / 4.3·10⁻¹⁰` (75) and `3.2·10⁻¹⁰ / 7.6·10⁻⁹` (134). `φ_b,static` differs in the 7th decimal between the grids (2.2780 vs
2.2628·10⁻⁵, from the landing step); irrelevant.

**Regression against Chat 13.** `v2_t01_control` (t = 0.1, closure 2, 128 pts/wall): late-window rate 1.62697 vs 1.62702 theory and
Chat 13's 1.62694; `v2_t01_plus` reproduces the Chat 13 picture (`H_J/H₀ = 0.627` at `H₀τ = 7.38`).

## 2. Reproduction of the key numbers from the saved data (`ref_analyse.py`, independent of `analyse_v2.py`)

Sliding 1-Hubble-time windows of `ln|φ_b − φ_b,static|` restricted to `|Δφ_b| < 10⁻³`:

| run | pts/wall | windows 2–3, 2.5–3.5, 3–4, 3.5–4.5, 4–5 | theory 1.65719 |
|---|---|---|---|
| `tiny_c4` (Δc = 10⁻⁸) | 75 | 1.65717, 1.65716, 1.65716, 1.65715, 1.65714 | −3·10⁻⁵ |
| `tiny_c4_fine` | 134 | 1.65719, 1.65719, 1.65719, 1.65718, 1.65718 | ≲10⁻⁵ |
| `ref_c4_dc3e-8` (mine, Δc = 3·10⁻⁸) | 75 | 1.65714, 1.65714, 1.65713, 1.65713, 1.65711 | seed-independent to 2·10⁻⁵ |
| `ref_c4_dzf3e-4` (mine) | 40 | 1.65642, 1.65641, 1.65640, 1.65639, 1.65638 | −8·10⁻⁴ |
| closure 2 (`tiny_c2.log`, t = 4–6) | 75 | 1.897 | +0.24 |
| closure 2 (`tiny_fine.log`, t = 2–3) | 134 | 1.71 | +0.06 |
| closure 2, `ref_c2_dzf3.75e-5_long` (mine) | 222 | 1.6535, 1.6544, 1.6551, 1.6556, 1.6559, 1.6562, 1.6563, 1.6565 (windows 2–3 … 5.5–6.5, still rising, increments shrinking) | −7·10⁻⁴ at t = 6, → ≈ −3·10⁻⁴ extrapolated |

The README's "1.6571–1.6572" is confirmed (at 134 points the agreement in windows 2–4 is 10⁻⁶). The rate drifts down after t ≈ 5
(1.6570 at 5–6, 1.6564 at 7–8) because `Δφ_b` approaches the linear cut, not for numerical reasons. `b_b` grows at 1.51–1.64 in the
same windows: it carries a non-modal (gauge/pulse) component and is rightly not used as the rate estimator.

**+1 side, at fixed brane proper time** (linear interpolation in `H₀τ`; `dz_c` is the coarse spacing, `finecoarse` is the author's run
with `dz_c = 2·10⁻³` that finished during this review):

| `H₀τ` | 75 pts, `dz_c = 8·10⁻³`: `H_J/H₀`, `φ_b`, `b_b`, `t` | 134 pts, `dz_c = 8·10⁻³` | 82 pts, `dz_c = 2·10⁻³` |
|---|---|---|---|
| 6.0 | 0.78886, 0.7535, −0.269, 6.155 | 0.78886, 0.7535, −0.257, 6.148 | 0.78882, 0.7535, −0.298, 6.172 |
| 6.5 | 0.68458, 0.9517, −0.630, 6.927 | 0.68457, 0.9517, −0.598, 6.904 | 0.68342, 0.9517, −0.753, 6.995 |
| 6.7 | 0.64711, 0.9819, −0.887, 7.349 | 0.64954, 0.9818, −0.830, 7.309 | 0.65693, 0.9813, −1.123, 7.501 |
| 6.8 | 0.62302, 0.9899, −1.155, 7.623 | 0.62438, 0.9899, −1.059, 7.563 | 0.64688, 0.9891, −1.410, 7.854 |
| 6.9 | 0.61106, 0.9946, −1.760, 8.040 | 0.61181, 0.9946, −1.571, 7.926 | **0.63864**, 0.9938, −1.853, 8.360 |
| 7.0 | 0.54374 (`ext`), 0.9972, −3.93, 9.54 | — | 0.63192, 0.9966, −2.753, 9.323 |
| 7.05 | (ext: τ saturates at 7.016) | — | 0.6287 at τ = 7.053 (t = 11.0, end of run) |
| `dH_J/dτ`, 6.7–6.9 | −0.167 | −0.185 | −0.090 |
| max bulk monitor, τ ≤ 6.9 | 0.65 | 0.63 | 0.0063 |
| max shell monitor, τ ≤ 6.9 | 3.8·10⁻⁵ | 1.9·10⁻⁵ | 5.8·10⁻⁵ |

The two fine resolutions (same coarse grid) agree to 0.1 % in `H_J(τ)`; the refined coarse grid agrees with them to 5·10⁻⁵ at τ = 6.0
and 0.17 % at 6.5, then departs: 1.5 % (6.7), 3.8 % (6.8), 4.5 % (6.9), 16 % (7.0). The gauge lapse `b_b` is unconverged in every
comparison (5–20 % between fine resolutions, 10–20 % between coarse resolutions), so the coordinate time of a given τ-slice differs by
up to 0.4 between grids. The tiny-seed runs reach the same end state as the Δc = 10⁻⁴ runs on the same grid (0.6092/0.6104 at
φ_b = 0.995, shifted by 5.58 Hubble times; linear theory 5.56): the attractor is seed-independent, but that is a statement about the
grid's attractor, and the grid matters at the 4 % level from τ ≈ 6.7 on.

**Throat side** (interpolated to the exact crossings; both runs on the default coarse grid):

| quantity | 75 pts | 134 pts | README | 4D EFT (Chat 9) |
|---|---|---|---|---|
| turnaround `H₀τ` | 5.8947 | 5.8944 | 5.90 | |
| turnaround `φ_b` | −1.9417 | −1.9412 | −1.953 | −1.94 |
| `φ_b = −1` crossing `H₀τ` | 5.6085 | 5.6083 | 5.61 | |
| `H_J/H₀` there | 0.7978 | 0.7978 | 0.793 | 0.797 |
| `φ_b = −1/c = −1.673` at `H₀τ` | 5.8335 | 5.8333 | | |
| `H₀τ` at `H_J = −5, −10, −30 H₀` | 6.1999, 6.2829, 6.3139 | 6.1980, 6.2763, 6.3107 | 6.20 | |
| `H_J = −100, −300 H₀` | 6.3217, 6.3226 | run ended at −39 (τ = 6.3130) | −300 by 6.32 | |
| max shell monitor / bulk monitor, whole run | 8·10⁻⁴ / 0.53 | 3.5·10⁻⁶ / 0.12 | | |

The README's −1.953 and 0.793 are the first *recorded samples* after the crossings (`analyse_v2.py` uses `argmax`; `φ_b` moves by 0.011
per 0.01 t there). The interpolated values −1.942 and 0.798 agree with the 4D theory to 0.1 %; please quote those. The two fine
resolutions agree to 3·10⁻⁴ in τ down to `H_J = −5`, 0.1 % at −10, 0.05 % at −30, and the fine run stops at −39; "−300 H₀ by 6.32"
is a single-resolution number. **No throat-side run with the refined coarse grid exists**; given that the +1 side agrees across coarse
grids to 0.17 % up to τ = 6.5 and departs afterwards, the turnaround (τ = 5.89) is probably safe and the collapse beyond τ ≈ 6.2 is
not checked against the coarse grid. The closure-2 run at 75 points (`v2_t1e3_minus`) has the turnaround at τ = 5.185, φ_b = −1.826
(12 % off): closure 2 at this resolution is not usable at t = 10⁻³.

## 3. My runs (all `--tdet 1e-3 --L 10 --dc 1e-8` unless stated)

| run | purpose | result |
|---|---|---|
| `ref_c4_dc3e-8` (closure 4, 75 pts, Δc = 3·10⁻⁸, t_f = 6.5; 374 s) | seed independence | 1.65714 (windows 2–4) vs 1.65717 for Δc = 10⁻⁸: **seed-independent to 2·10⁻⁵**; the residual is the earlier onset of the nonlinear correction for the 3× larger seed |
| `ref_c4_dzf3e-4` (closure 4, 40 pts, t_f = 7; 220 s) | wall-resolution convergence | 1.65640 ± 0.00002, error −7.9·10⁻⁴ vs −3·10⁻⁵ at 75 and ≲10⁻⁵ at 134: ratio 27 for a 1.9× coarser grid (≈ h⁵). **75 points per wall is converged to 3·10⁻⁵ in the rate** |
| `ref_c2_dzf3.75e-5` (closure 2, 222 pts, t_f = 3.5; 448 s) | does closure 2 converge to 1.657? | windows 1–2 … 2.5–3.5: 1.6490, 1.6519, 1.6535, 1.6544 — still drifting upward, so rerun longer: |
| `ref_c2_dzf3.75e-5_long` (same, t_f = 6.5; 740 s) | | 1.6551, 1.6556, 1.6559, 1.6562, 1.6563, 1.6565 (windows 3–4 … 5.5–6.5): monotone approach to 1.6572 **from below**, increments 5, 3, 3, 2, 1·10⁻⁴; the +0.24 (75) and +0.06 (134) excesses are gone. Closure 2 does converge toward the certified rate, slowly, and its approach is not a clean power law (sign change between 134 and 222 points), which is what one expects when a first-order boundary error competes with the interior fourth-order error |

The closure-2 sequence 1.90 (75) → 1.71 (134) → 1.6565 (222, still rising) versus closure 4's 1.6572 at 75 confirms that at any resolution the
author could afford, the symmetric closure's first-order `D2` at the shell dominates the rate error, and that the error does go
away with resolution — both halves of "closure, not only resolution".

## 4. The "clean window" on the +1 side

What supports a cut at `H₀τ = 6.9` on the default grid: the shell-layer monitor is ≤ 4·10⁻⁵ up to there (10⁻⁴ at t = 8.5, 10⁻³ at
t = 10, 9·10⁻³ at t = 12); the two fine resolutions agree to 0.1 %; the lapse has fallen to `e^{b_b} = 0.17–0.21` and proper time
saturates near 7.02 — the chart is ending. What does not: **the refined coarse grid gives a different `H_J` from τ ≈ 6.6 on** (Section 2
table), and the default grid's own bulk monitor tells the story. Its maximum sits (i) at `z = −t` for t ≲ 3 (the seed pulse the README
describes; ≤ 5·10⁻⁵), (ii) at a fixed `z ≈ −0.92` (ρ ≈ 1.7) for 3.5 ≲ t ≲ 7, growing at the unstable mode's rate (×2.3 per 0.5 Hubble
time, 8·10⁻⁵ → 7 %) — present at all fine resolutions and on the refined coarse grid at a 2.4× smaller level, so a property of the
discretised mode's tail in the cone region and not yet understood, and (iii) from t ≈ 7.2 on, a second light-like pulse launched
from the shell at t ≈ 7.0 (τ ≈ 6.54): `z = −(t − 7.0)`, 32 % at t = 7.5, 65 % at t = 8.0 (τ = 6.89), 88 % at t = 10. On the refined
coarse grid the same pulse stays ≤ 0.9 % throughout. The 65 % is inside the quoted window, and the divergence of `H_J` between the
coarse grids begins exactly when this pulse appears. Whether the pulse is created by the coarse cells themselves or by the 17.7 %
per-cell transition at z ≈ −0.06 (both change with `dz_c`) I cannot separate; either way it belongs to the default grid.

Recommended statement: on the default grid the clean window ends at τ ≈ 6.6 (`H_J/H₀ ≈ 0.67`, still falling); on the refined grid the
window extends to at least τ ≈ 7.0 with `H_J/H₀ = 0.632`, and to τ = 7.05 with 0.629 (shell monitor 8·10⁻⁴ at the end); a third coarse
resolution (e.g. `dz_c = 4·10⁻³` or `10⁻³`) is needed before any plateau digit is quoted. The README's trend table row for t = 10⁻³
("0.61, still falling") should read "0.63–0.64 (two coarse resolutions 0.611 / 0.639 at τ = 6.9; 0.63 flattening on the finer one)
versus RS 0.6067". The statement "bounded above by the quoted value, approaching the RS value from above" should become "the
better-resolved run sits 3.5 % above RS and is flattening; convergence to RS is consistent with, not shown by, these runs".
Note that this is the same situation the Chat 13 referee found at t = 0.1 (0.649 at the freeze, RS 0.638), so the trend across
detunings is "2–4 % above RS at the freeze in every case", not the tightening sequence the README presents.

The late negative `H_J` in `plus_c4_ext`: the first negative sample is at t = 11.89 (τ = 7.0160) with `1 + a_t = −8.4·10⁻⁶` and
`e^{b_b} = 1.5·10⁻³`. Since `H_J/H₀ = (1 + a_t)/e^{b_b}`, a physical `0.6 H₀` at that lapse would need `1 + a_t ≈ 9·10⁻⁴`; an absolute
error of 10⁻³ in `A_t` (static value 1) is amplified 700× by the collapsed lapse. Plainly numerical, a 0/0 of the lapse collapse — and
the README should say that rather than "not trustworthy". On the refined grid `H_J` never goes negative (0.629 at the end, t = 11).

The receding wall: at the end of the 75-point run (t = 8.12) the φ = 0 crossing is at z = −0.29 (ρ = 5.5), i.e. the wall moved 0.3 in
z in the 2.5 Hubble times after the roll-off; between t = 8 and 12 it moves from −0.29 to −3.84, i.e. close to the speed of light. Its
conformal width at z = −0.3 is ~1/ρ ≈ 0.2, ~25 coarse cells, so "leaves the finely resolved region" is true but "unresolved" would be
an overstatement; the coarse-region second-difference indicator at the wall is 3·10⁻³ on the default grid and 7·10⁻⁶ on the refined one.

## 5. Is "the closure, not only resolution, was Chat 13's problem" fair to the Chat 13 referee?

The Chat 13 referee measured closure-2 rates 5 → 4 → 2.2 at 13, 25, 51 points per wall, called the downward convergence "what a
resolution problem looks like, not a boundary-closure defect", listed the closure hypothesis as **not supported**, and then
recommended as fix (c) exactly "4th-order one-sided (5-point) stencils … with the Neumann value substituted, which removes the
O(h f''') term". Chat 14's data continue the sequence (1.90 at 75, 1.72 at 134, 1.6565 and rising at 222 in my run) and show that at fixed
resolution the closure order alone moves the rate from 1.90 to 1.6572. So the Chat 13 referee was right that closure 2 converges and
right about the cure, but wrong to call the closure hypothesis unsupported: the boundary's first-order second derivative is the
dominant error at every affordable resolution. Chat 14's sentence is fair in substance. The README should (i) add that the Chat 13
report proposed the closure fix, (ii) not call closure 2 "second-order at the boundary" (its `D2` is first-order), and (iii) not
revive the Chat 13 author's "BPS cancellation" wording — what is demonstrated is an ordinary boundary truncation error whose size is
large because `ρ_b² ≈ 6200` multiplies the potential terms, which is what the Chat 13 referee said.

## 6. Summary table

| claim (README) | status |
|---|---|
| Coordinate map, chain rule, CFL, taper, `ġ` slopes correct | **Established** (Section 1) |
| Closure-4 algebra as described; stable | **Established** (algebra to 10⁻¹⁵; constant-coefficient spectral radius ≤ 1) |
| "Kreiss–Oliger dissipation" | **Not what the code does**: 2500× weaker than v1 (h = 1 in ξ); runs stable regardless; must be documented |
| Linear rate 1.6571–1.6572 at 75 and 134 pts/wall | **Established** (1.65714–1.65719 vs 1.65719); seed-independent to 2·10⁻⁵; 40 pts give 1.6564; 75 pts converged to 3·10⁻⁵ |
| Closure 2 gives 1.90 / 1.72 and converges slowly | **Numerically supported** (logs; my 222-point run: 1.6565 and still rising after 6.5 Hubble times) |
| "closure, not only resolution" | **Fair**, with the credit noted in Section 5 |
| +1 side: `H_J/H₀ = 0.612/0.610` at `H₀τ = 6.9`, relaxing toward 0.6067 | **Not established**: 0.611 on the default coarse grid, 0.639 on the 4× refined one; the latter flattens at 0.63 (τ = 7.0–7.05); the fine-resolution pair shares the coarse grid and did not test this |
| Clean window ends at 6.9 | **Too late for the default grid** (coarse grids separate from τ ≈ 6.6; 65 % bulk violation from a second pulse inside the window); the refined grid is clean to τ ≈ 7.0 |
| Late negative `H_J` numerical | **Established** (0/0 of the collapsed lapse; absent on the refined grid) |
| Throat: turnaround 5.90, φ_b = −1.953; φ = −1 at 5.61 with 0.793 | **Established as 5.895, −1.942; 5.608, 0.798** (interpolated; two fine resolutions; default coarse grid only) |
| "−300 H₀ by 6.32", "everything up to −10 H₀ converged" | −10: yes (0.1 %, two fine resolutions); −300: single resolution; no coarse-grid check on the throat side |
| Dominant bulk violation = seed pulse at z = −t | **Incomplete**: true for t ≲ 3; then a stationary maximum at z ≈ −0.92 growing at the mode rate; then a second pulse from the shell at t ≈ 7 that the refined grid does not have |
| "Same two fates, no radiation branch" at t = 10⁻³ | **Numerically supported** (homogeneous, classical, matter-free; both grids give an inflating brane with `H_J/H₀ ≈ 0.61–0.64` and a throat-side collapse) |

Suggested next runs (cheapest first): (1) `--dzc 4e-3` and `--dzc 1e-3` on the +1 side to bracket the plateau (the 2·10⁻³ run cost
24 min); (2) a throat-side run with `--dzc 2e-3` to t = 7 (`≈ 15 min`); (3) fix the KO scaling (`ε/z'_i`) and confirm nothing
changes; (4) a gentler map (`width = 0.05`) to separate the transition-zone reflection from coarse-cell dispersion.

## Appendix — files

`ref_closure_check.py` / `closure_check_output.txt` (Section 1); `ref_analyse.py`, `REFEREE_ANALYSIS.json` (author's runs),
`REFEREE_ANALYSIS_refruns.json` (my runs), `analysis_output_authors_runs_v2.txt` (complete text output including the `finecoarse`
and `minus_c4_fine` runs; `analysis_output_authors_runs.txt` is the earlier pass without them) (Sections 2, 4); `runs/ref_c4_dc3e-8*`,
`runs/ref_c4_dzf3e-4*`, `runs/ref_c2_dzf3.75e-5*`, `runs/ref_c2_dzf3.75e-5_long*` (Section 3).
