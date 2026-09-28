# Referee report (claims) — Chat 14, roll-off at the registered detuning t = 10⁻³

Referee: claims referee (Claude, adversarial role). Date: 22 September 2026, 20:45 PDT.
Scope: `00_READ_FIRST.md` (it sits in `untitled folder 151/`, one level above the run folder, with links relative to that
location), `runs/*.log`, `runs/*_summary.json`, `runs/ANALYSIS_V2.json`, `runs/*_timeseries.npz`, `analyse_v2.py`,
`make_figure_v2.py`, `rolloff5d_v2.py`, and Chat 13's refereed README and two referee reports. No runs were launched; all
numbers below were recomputed from the stored time series with numpy (interpolation between records where the README uses
"first record past the threshold"). Nothing outside `referee_claims/` was written.

Legend: **[EST]** established by the stored data · **[NUM]** numerically supported with the stated caveat · **[NOT]** not
established by the evidence in the folder · **[ERR]** wrong number or wrong statement.

## 0. Verdict in one paragraph

The central numbers are in the logs and most are quoted correctly: the closure-4 linear growth rate (1.6572 in windows
t = 1–4.5 at 75 and 134 points per wall, against 1.65719), the closure-2 rate 1.90 at 75 points, the blow-up of the
closure-2 dc = 10⁻⁴ run at t = 4.71, φ_b = 0.995 at H₀τ = 6.91, the throat crossing at H₀τ = 5.61, the turnaround at
H₀τ = 5.89, H_J = −5H₀ at 6.20, proper time saturating near 6.32–6.33 and 7.02, the exact RS value 0.6067 and the Chat 9
EFT values (0.606, 0.797, −1.94). Three things are stronger than the evidence and one is contradicted by the folder's
own logs. (i) The README says the dominant bulk constraint violation is the outgoing seed pulse at z = −t; the logs show
that is true only for t ≲ 3. From t ≈ 3.5 the maximum sits at z ≈ −0.8 to −1.0 (behind the receding wall, in the
coarse-grid region, dz = 8×10⁻³ in *both* "resolutions"), grows with the mode, and reaches 0.32 at t = 7.5 and
0.65–0.69 at t = 8 — exactly where the headline value H_J/H₀ = 0.61 is read off. At the corresponding phase the t = 0.1
run had < 10⁻². By the standard Chat 13's referees applied (plateau quoted only while the monitor is small), the last
constraint-clean point at t = 10⁻³ is H_J/H₀ ≈ 0.78 (monitor 10⁻², H₀τ = 6.0) or ≈ 0.66 (monitor 0.1, H₀τ = 6.63); the
0.61 at 6.9 must be quoted with its monitor value, and the "two resolutions" test does not cover it because both runs share
the same coarse spacing where the violation lives. **The `v2_t1e3_plus_c4_finecoarse` run (dz_c = 2×10⁻³, 82 points per
wall), which reached t = 8.97 while this report was being written, settles the point against the README**: its bulk
monitor stays ≤ 9×10⁻³ through H₀τ = 6.97 (140× smaller than 0.65), but its brane Hubble rate is H_J/H₀ = 0.644 at
H₀τ = 6.83, 0.637 at 6.92 and 0.634 at 6.97, against 0.62 / 0.609 / 0.58 in the dz_c = 8×10⁻³ runs at the same proper
times; the two agree to ≤ 0.007 up to H₀τ ≈ 6.5 and then separate. The headline "0.61, still relaxing toward 0.607" is
therefore a coarse-region truncation error of about 0.03 (5 %); the cleaner number is ≈ 0.63–0.64 H₀ at the freeze, 4–5 %
above the RS value and still falling at dH/dτ ≈ −0.07 H₀², i.e. the same situation as Chat 13's t = 0.1 run (0.649 vs
0.638), not a "clean trend" toward the exact value. Whether 0.634 is itself converged in dz_c is untested (one refinement
step only). (ii) The resolution comparison in the validation
table compares the two runs at different proper times (H₀τ = 6.892 vs 6.911); at equal H₀τ = 6.90 the numbers are
φ_b = 0.9946 (both) and H_J/H₀ = 0.6111 (75) / 0.6118 (134) — closer than stated, and the finer grid is *higher*, so
"0.612 (0.610 at the finer grid)" in Section 2 is an artefact of the time mismatch. (iii) "Everything up to H_J ≈ −10H₀ is
converged" is not: at equal H₀τ the two resolutions differ by 1.4 % at −5H₀, 8 % at −10H₀, 15 % at −14H₀ and 23 % at
−22H₀; Chat 13's refereed wording ("1–2 % until −5H₀, not beyond") still applies. (iv) The framing words — "closes the last
open item", "Chats 9–14 now form a closed, refereed account", "converged", "clean" — outrun a folder whose new numerics
(stretched grid, one-sided closure, d/dt Neumann slopes) have not yet been refereed and whose decisive resolution test is
still running. The qualitative conclusion (same two fates as Chat 13, no radiation branch, within the homogeneous
classical single-scalar matter-free system) is supported.

## 1. Number-by-number check of `00_READ_FIRST.md`

| README statement | What the folder shows | Status |
|---|---|---|
| Growth rate 1.6571–1.6572 (75 pts, windows t = 1–5, dc = 10⁻⁸) | ANALYSIS_V2: 1.6571, 1.6572 ×5, 1.6571 ×3 for windows 1.0–5.0 | [EST] |
| Same, 134 pts: 1.6572 | 1.6572 for windows 1.0–4.5; 1.6571 at 5.0 | [EST] |
| "reproduces the certified linear growth rate to five digits" | early windows 1.6572 vs 1.65719 (6×10⁻⁶); the automatic fit over the whole linear window (t = 2.2–7.6) gives 1.65704 / 1.65707 (1×10⁻⁴); windows drift to 1.6564 by t = 7 (deviation 7×10⁻⁴, onset of nonlinearity or numerical) | [NUM] say "to 10⁻⁵ in the early windows, 10⁻⁴ over the full linear window" |
| Certified theory: μ² = −7.7178716 → 1.65719 | √7.7178716 = 2.7781 … the quoted 5D rate 1.6571907 is the Chat 9 conversion (S_LIN in analyse_v2.py); consistent with Chats 9–13 | [EST] (inherited) |
| Closure 2: 1.90 at 75 pts/wall | `v2_t1e3_tiny_c2.log`, t = 4–6: ln-slope 1.897 | [EST] |
| Closure 2: 1.72 "at 150" pts/wall | the run is `v2_t1e3_tiny_fine` at **134** points, stopped at t = 3; slope 1.72 (t = 2.0–2.5), 1.69 (2.5–3.0); a t ≤ 3 window, no sliding windows 1–5 | [ERR] number of points; [NUM] rate (short window) |
| dc = 10⁻⁴ closure-2 run blew up at t = 4.7 | `stop_reason: non-finite`, t_end = 4.708 | [EST] |
| Reproduction of Chat 13 at t = 0.1: φ_b(t = 5) = 0.20386 (Chat 13: 0.20386), rate 1.627 | v2: 0.2038595; Chat 13 `t01_plus_v2`: 0.2038589; rate 1.6271 from the zero-seed control (round-off seeded) | [EST] — but both v2 t = 0.1 runs use closure 2, so this validates the grid map, not closure 4 |
| "+1 side, two resolutions at H₀τ = 6.9: φ_b = 0.9943 / 0.9950, H_J/H₀ = 0.6122 / 0.6104" | these are the last log lines, at H₀τ = 6.892 (75) and 6.911 (134). At equal H₀τ = 6.90: φ_b = 0.99462 / 0.99461, H_J/H₀ = 0.61106 / 0.61181 | [ERR] mismatched times; corrected numbers agree better |
| Section 2: "0.612 (0.610 at the finer grid) at H₀τ = 6.9" | see above: 0.611 / 0.612 at equal time; the fine grid is higher, not lower | [ERR] |
| "reaches φ_b = 0.995 by H₀τ ≈ 6.9" | stop criterion φ_b > 0.995 fired at H₀τ = 6.912 / 6.913 | [EST] (note: the window ends because of the stop criterion, not a physics judgement) |
| "falls monotonically, without oscillation" | H_J non-increasing after φ_b > 0.5 in both runs to t = 8 | [EST] |
| "still relaxing slowly" | dH_J/dτ ≈ −0.10 to −0.12 H₀² over H₀τ = 6.85–6.91 in the dz_c = 8×10⁻³ runs (3× steeper than Chat 13's −0.035 at t = 0.1); the extended run crosses the RS value 0.6067 at H₀τ = 6.925, i.e. 0.03 Hubble times after the quoted point. In the finecoarse run the slope is −0.07 H₀² and H_J stays ≥ 0.634 to H₀τ = 6.97 | [ERR] the quoted value and slope are coarse-region artefacts; see §2.1 |
| Table row t = 0.001: "0.61 (still falling)" vs RS 0.6067; "the trend across detunings is now clean" | finecoarse: 0.644 (H₀τ = 6.83), 0.637 (6.92), 0.634 (6.97), monitor ≤ 9×10⁻³, still falling | [ERR] at the 0.03 level; the row should read "≈ 0.63–0.64 (still falling)", 4–5 % above RS, and the "clean trend" sentence goes |
| "lapse collapses … proper time saturates near H₀τ ≈ 7.0" | ext run: H₀τ = 7.016 at t = 12, b_b = −6.6 | [EST] |
| "the domain wall, receding into the bulk at light speed, leaves the finely resolved region" | at t = 8 the φ = 0 point of the wall is at z = −0.28 (75) / −0.24 (134), the fine region is |z| ≲ 0.08; the wall left it near φ_b ≈ 0.9 (H₀τ ≈ 6.3), well inside the quoted window; it moves at coordinate light speed only after t ≈ 8 (z = −3.8 at t = 12) | [ERR] as timing; the wall is on the coarse grid for the whole quoted plateau |
| "numbers after H₀τ ≈ 6.9 are not trustworthy" (H_J < 0 by t = 12) | M_shell rises from 4×10⁻⁵ to 9×10⁻³, M_bulk to 0.88, b_b to −6.6 — plausible, but the positive test (finecoarse) is pending | [NUM] |
| Exact RS value 0.6067 | √(σ₁²/36 − 1/81)·ρ_b with σ₁ = 2/3 + t(1 + c), ρ_b = 78.828: 0.60673 | [EST] |
| "also Chat 9's four-dimensional prediction (0.606)" | O(t) form √((1+c)t/27)·ρ_b = 0.6064 (BLAST_RESULTS.json: 0.60640) | [EST] |
| Table row t = 0.1: 0.649 (still falling) / 0.638 | Chat 13 refereed values; RS recomputed 0.63787 | [EST] |
| Table row t = 0.03: 0.63 / 0.616 | Chat 13: 0.637 at the φ_b = 0.995 stop; RS recomputed 0.61609 | [EST]; column header "at the freeze" is wrong for 0.03 and 0.001 (both ended at the stop criterion) |
| "61 % of its original expansion rate" | mixes the measured 0.61 (still falling) with the RS 0.6067 | [NUM] fine if phrased "≈ 0.61 H₀, consistent with 0.607" |
| Throat: φ = −1 crossing at H₀τ = 5.61, H_J = 0.793 H₀ (4D: 0.797) | first-record-past-threshold: 5.613 / 0.7927 (75), 5.612 / 0.7940 (134). Interpolated: 5.6085 / **0.7978** at both resolutions; EFT 0.7973 | [EST]; the agreement with the EFT is 6×10⁻⁴, better than quoted |
| Turnaround at H₀τ = 5.90, φ_b = −1.953 (4D: −1.94) | first-record: −1.9529 (75) but **−1.9412 (134)**; φ_b changes by 0.03–0.04 per record here. Interpolated: H₀τ = 5.8947 / 5.8944, φ_b = **−1.9417 / −1.9412**; EFT −1.9418 | [EST] as H₀τ; [ERR] −1.953 is a sampling artefact (same artefact gave Chat 13 its −1.96; Chat 13's run interpolates to −1.942) |
| H_J = −5H₀ at H₀τ = 6.20 | 6.1999 / 6.1980; resolutions differ 1.4 % at equal τ | [EST] |
| "−300H₀ by 6.32" | 75-pt run only (−300 at 6.3226 = last record); 134-pt run ended at t = 7.49 with −39 at 6.313; at equal H₀τ = 6.31 the two give −21.6 vs −28.1 | [NUM] single resolution, non-converged regime |
| "proper time saturating near H₀τ ≈ 6.33" | τ(t = 8) = 6.3226 with dτ/dt → 0; the run ended there (t_f = 8) | [NUM] short extrapolation |
| "everything up to H_J ≈ −10H₀ is [converged]" | at equal H₀τ: 0.2 % at −0.8, 1.4 % at −5, 3.6 % at −7.5, **8 % at −10**, 15 % at −14, 23 % at −22 | [ERR] as stated; Chat 13's "1–2 % to −5H₀" is the honest version |
| "at t = 0.1 the turnaround was at H₀τ = 5.99, φ_b = −1.96, collapse 0.44 Hubble times later" | 5.988; φ_b interpolated −1.942 (−1.956 first-record); 0.44 is Chat 13's *extrapolated* τ* − τ_turn | [NUM] label −1.96 → −1.94 and 0.44 as extrapolated |
| "The registered detuning behaves the same way, slightly earlier in proper time" | turnaround 5.89 vs 5.99; −5H₀ at 6.20 vs 6.30; τ* − τ_turn ≈ 0.43 vs 0.44 | [EST] |
| "the momentum constraint … now reads ≲ 10⁻⁴ there throughout" (shell layer) | max M_shell: 4.4×10⁻⁵ / 2.2×10⁻⁵ (+1), 3.5×10⁻⁶ (throat, 134), **8.2×10⁻⁴** (throat, 75, at t = 8 in the non-converged collapse), 9×10⁻³ (discarded ext window) | [NUM] "≲ 10⁻⁴ within the quoted windows" |
| "The dominant bulk constraint violation is an outgoing pulse … z = −t … never returns" | true for t ≤ 3 (logs: @z = −0.5, −1.0, …, −3.0). From t = 3.5 the maximum is at z ≈ −0.92 (fixed), then −0.84 … −1.02, growing ∝ the mode: 7.7×10⁻⁵ (t = 3.5), 9.5×10⁻⁴ (5), 7.3×10⁻³ (6), 6.7×10⁻² (7), 0.32 (7.5), 0.65 (8). max_bulk_constraint = 0.688 (75) / 0.677 (134): **does not decrease with the fine-region resolution** because it lives in the shared coarse region. Same in the dc = 10⁻⁸ runs (0.69 / 0.63) and on the throat side (0.53 / 0.12 — there it does improve) | [ERR] as a description of the runs; the README does not report M_bulk at all |
| "1,400–1,950 points instead of the 60,000 a uniform grid would need" | 1392 / 1458 / 1954; L/dz_f = 10/1.5×10⁻⁴ ≈ 67,000 | [EST] |
| "crosses the wall in about one Hubble time" (copied from Chat 13) | φ_b from 0.1 to 0.9 takes H₀τ 4.45 → 6.30 = 1.85 Hubble times | imprecise; "about two" |
| "Established (floating point, two resolutions, referee-checked code lineage)" | v1 was refereed; the v2 additions (grid map with z′, z″ chain rule, closure 4, `neumann_t`, taper) are not; the two resolutions differ only inside |z| ≲ 0.08 | [NOT] as worded |
| "Not established: … asymptotic value … (bounded above by the quoted value, approaching the RS value from above)" | "bounded above by 0.61" needs the monitor caveat (§2.1); "approaching from above" is an expectation — the run crosses 0.6067 0.03 Hubble times later and keeps falling | [NOT] |
| "This closes the last open item"; "Chats 9–14 now form a closed, refereed account" | Chat 14 is under review; its numerics referee has not reported; finecoarse pending; Chat 13's referees' open items (lapse/slicing to get the asymptotic H, inhomogeneity, brane matter) stand | [NOT] |

## 2. What a hostile reader will attack, in order

### 2.1 The bulk constraint monitor at the quoted +1 plateau, and what the coarse-region refinement shows (most damaging)

Facts from the logs (75-pt run; the 134-pt run is the same within 10 %):

| t | H₀τ | φ_b | H_J/H₀ | M_bulk (relative) | location |
|---|---|---|---|---|---|
| 5.10 | 5.07 | 0.258 | 0.957 | 10⁻³ | z ≈ −0.9 |
| 6.20 | 6.03 | 0.771 | 0.781 | 10⁻² | z ≈ −0.8 |
| 7.20 | 6.63 | 0.974 | 0.663 | 0.1 | z ≈ −0.8 |
| 7.50 | 6.76 | 0.987 | 0.632 | 0.32 | z ≈ −0.6 |
| 8.00 | 6.89 | 0.994 | 0.612 | 0.65 | z ≈ −1.0 |

For comparison the t = 0.1 run (`v2_t01_plus`) first exceeds 10⁻³ at H₀τ = 7.11 (H_J = 0.658) and 10⁻² at 7.35
(0.644): Chat 13's quoted 0.649 was read where the monitor was < 10⁻². Chat 13's referees excluded the frozen window
precisely because the constraint became order one there. At t = 10⁻³ the same criterion excludes the quoted 0.61. Two
mitigating facts: the shell-layer monitor is clean (≤ 4×10⁻⁵), and the H_J(τ) curves of the two runs agree to 10⁻³ — but
both runs use dz = 8×10⁻³ at z ≈ −1, so that agreement does not test the region where the violation lives. I could not
decide whether the relative normalisation makes the monitor pessimistic there (the time-derivative fields are not stored),
so the monitor alone would not have been a verdict. The finecoarse run is (log lines at 20:59 PDT; the run was still going,
its `_summary.json` and time series not yet written; `phimax` was evidently raised, since it passed φ_b = 0.995 without
stopping):

| run | dz_c | t | H₀τ | φ_b | H_J/H₀ | M_bulk |
|---|---|---|---|---|---|---|
| finecoarse | 2×10⁻³ | 5.98 | 5.854 | 0.668 | 0.825 | 1.5×10⁻³ |
| 75-pt (ext) | 8×10⁻³ | 6.00 | 5.876 | 0.681 | 0.820 | 7.3×10⁻³ |
| finecoarse | 2×10⁻³ | 6.98 | 6.492 | 0.950 | 0.685 | 2.5×10⁻³ |
| 75-pt (ext) | 8×10⁻³ | 7.00 | 6.537 | 0.959 | 0.679 | 6.7×10⁻² |
| finecoarse | 2×10⁻³ | 7.48 | 6.692 | 0.981 | 0.658 | 2.6×10⁻³ |
| 75-pt (ext) | 8×10⁻³ | 7.50 | 6.757 | 0.987 | 0.632 | 0.32 |
| finecoarse | 2×10⁻³ | 7.98 | 6.828 | 0.991 | **0.644** | 4.5×10⁻³ |
| 75-pt (ext) | 8×10⁻³ | 8.00 | 6.892 | 0.994 | **0.612** | 0.65 |
| finecoarse | 2×10⁻³ | 8.47 | 6.917 | 0.994 | **0.637** | 7.2×10⁻³ |
| 75-pt (ext) | 8×10⁻³ | 8.50 | 6.956 | 0.996 | 0.597 | 0.79 |
| finecoarse | 2×10⁻³ | 8.97 | 6.974 | 0.996 | **0.634** | 8.6×10⁻³ |
| 75-pt (ext) | 8×10⁻³ | 9.00 | 6.985 | 0.997 | 0.577 | 0.85 |

Reading: the two agree to ≤ 0.007 in H_J at equal proper time up to H₀τ ≈ 6.5, separate from H₀τ ≈ 6.7 on, and differ by
0.03 at H₀τ = 6.9 — exactly where the README reads off its headline. The refined run's monitor is two orders of magnitude
smaller, its proper time freezes later (H₀τ = 6.97 at t = 9 and still advancing 0.1 per unit t), and its H_J is 0.634 and
falling at −0.07 H₀². So: (a) the README's 0.61 is a coarse-region artefact, (b) the "two resolutions agree" row was
testing the wrong region, (c) the right statement is "≈ 0.63–0.64 H₀ at the freeze, still falling, 4–5 % above the RS
value" — and (d) 0.634 is one refinement step, not a converged number; a dz_c = 10⁻³ run (≈ 4× the cost, ~80 min
single-threaded) is needed before any plateau value is quoted to two digits.

### 2.2 The "two resolutions" claim on the +1 side tests the shell layer only

`dz_fine` was halved; `dz_coarse` = 8×10⁻³ and the transition at |z| ≈ 0.06–0.08 are identical. The receding wall is at
z ≈ −0.25 to −0.3 by the end of the window and the constraint feature at z ≈ −0.6 to −1.0. The convergence statement must
say "two shell-layer resolutions"; the coarse-region resolution test is the finecoarse run.

### 2.3 Chat 13's diagnosis is contradicted without saying so

Chat 13 (refereed README §5, numerics report §4): the t = 10⁻³ failure "was misdiagnosed as a boundary-closure problem";
cause = under-resolution. Chat 14: "That was half the story." The evidence for the closure half is the 75-pt closure-2
rate 1.90 vs closure-4 1.657 — sound — but the 134-pt closure-2 number comes from a 3-unit run in the initial transient.
Write it as: "the resolution diagnosis was right (uniform-grid rates 5 → 4 → 2.2 → 1.90 at 6–75 points per wall); the
residual 15 % excess at 75 points is the second-order closure (1.90 → 1.657 on switching to the one-sided fourth-order
closure at the same grid)", and either rerun closure 2 at 134 points to t ≥ 6 or drop the "1.72" number.

### 2.4 The late collapse

"−300H₀ by 6.32", "proper time saturating near 6.33", "up to −10H₀ converged" — see the table. Keep Chat 13's refereed
formula: converged to 1–2 % until H_J ≈ −5H₀ (here H₀τ = 6.20), diverging beyond; both resolutions pass −30H₀ before
H₀τ = 6.315; the final approach and its 4D/5D nature are open. Note also that the 134-pt throat run was stopped at
t_f = 7.5 (H_J = −39), so nothing beyond −39H₀ has a second resolution.

### 2.5 Sampling artefacts in the turnaround numbers

`analyse_v2.py` takes the first record with H_J < 0; records are 0.01 apart in t and φ_b moves 0.03–0.04 per record at the
turnaround. The interpolated values −1.9417 (75) / −1.9412 (134) agree with the EFT −1.9418 to 10⁻⁴, and the interpolated
H_J at the φ = −1 crossing (0.7978) agrees with 0.7973 to 6×10⁻⁴. This is a *better* result than the README claims; fix
the script (linear interpolation) and quote it. The same artefact is in Chat 13's −1.96 (interpolates to −1.942).

### 2.6 Framing

"closes the last open item", "closed, refereed account", "converged", "settled" (the word does not appear, good), "the
trend across detunings is now clean" — replace by "computed", "numerically supported", "consistent". The Big Bang
sentence in §0 ("No branch produces a radiation-filled universe") must carry the Chat 13 qualifiers in the same sentence:
homogeneous, classical, one bulk scalar, no brane matter, within the computed window. §3 has them; §0 does not.

### 2.7 Small points

* "Validation at t = 0.001" table: the t = 0.1 reproduction row is closure 2; a closure-4 t = 0.1 run (≈ 2 min) would
  close the obvious question "does the new closure change the refereed t = 0.1 result?".
* `H0tau_end` for the 75-pt +1 run is 6.912, for the 134-pt 6.913; the log line quoted in the table is at t = 8.00 (6.892).
* "1.6571–1.6572 … vs 1.65719 theory" — also give the whole-window fit (1.65704 / 1.65707) so nobody accuses you of
  picking windows.
* The zero-seed controls at t = 10⁻³ stay at 0.0 exactly (as Chat 13's numerics referee predicted); the t = 0.1 control
  grows from a 10⁻¹⁴ solve difference. Say which is which.
* `_KILLED` logs: the claim "their trajectories agree with the reruns" — `v2_t1e3_minus_c4_KILLED.log` vs the rerun agree
  to all printed digits (checked t = 0.5–5.5); fine.

## 3. `make_figure_v2.py` and the caption

* **Time shift.** The faint t = 0.1 curves are shifted so that |φ_b| crosses 0.05 at the same H₀τ as the t = 10⁻³ curves;
  the shifts are −0.043 (+1) and −0.091 (throat) Hubble times, computed separately per branch. The README caption says
  only "time-shifted to align the roll-offs"; the shift values and the criterion are in SVG `<title>` tooltips, which the
  PNG loses. Put both in the caption ("shifted by −0.04 and −0.09 H₀τ so that |φ_b| = 0.05 coincides"). The small size of
  the shifts is itself a nice fact — say it.
* **Truncation.** `T_CLEAN = 8.0` cuts the solid +1 curve at H₀τ = 6.89; the only disclosure is the in-panel note "0.61 at
  H₀τ = 6.9, chart freezes beyond". The README caption does not say the curve is cut, nor that the faint t = 0.1 reference
  is *not* cut (it runs to H₀τ = 7.38, into the window Chat 13 declared unclean after 7.33). Either cut the reference at
  Chat 13's clean time or say in the caption that it is shown in full. State which runs are plotted (defaults: the 75-pt
  runs) and that the throat curves leave the panel at H_J/H₀ = −3 (y-range), not at a physics cutoff.
* The "0.61 at H₀τ = 6.9" label should not be the only number on the panel if §2.1 stands; consider marking where the
  bulk monitor passes 0.1 (H₀τ ≈ 6.6).
* Footer "Floating-point … Not observations; not a certificate" — good; keep.

## 4. Consistency with Chat 13's refereed wording

Chat 13 (after referee corrections) says: "relaxes toward the RS de Sitter brane … H ≤ 0.65 H₀ … when the chart freezes";
"the turnaround is converged; the final approach is not; 4D/5D nature open"; "no radiation branch for the homogeneous,
classical, single-scalar, matter-free system"; "registered detuning: qualitative behaviour expected to be the same; not
directly computed". Chat 14 should read as the same sentences with t = 10⁻³ substituted and the numbers updated — and
with one honest addition Chat 13 did not need: the +1 plateau at t = 10⁻³ is quoted with a bulk-constraint monitor of
0.65, pending the coarse-region refinement. Where Chat 14 currently departs from Chat 13's wording ("converged to
−10H₀", "closes", "closed, refereed", the closure diagnosis), Chat 13's wording is the safer one.

## 5. Safe abstract paragraph

> We repeat the Chat 13 five-dimensional evolution of the unstable shell at the registered tension detuning t = 10⁻³,
> which the uniform-grid code could not resolve. A smooth stretched grid (75 and 134 points across the wall) and a
> fourth-order one-sided Neumann closure at the shell reproduce the certified linear growth rate, 1.6572 against 1.65719,
> to 10⁻⁵ in the early windows and 10⁻⁴ over the full linear window; with the symmetric second-order closure the same grid
> gives 1.90, so both resolution and closure order were needed. Seeded toward φ = +1, the brane's expansion rate falls
> monotonically and without oscillation; when the shell reaches φ_b ≈ 0.995 and the conformal chart's lapse begins to
> collapse (H₀τ ≈ 6.9–7.0) it is 0.63–0.64 H₀ and still decreasing at about −0.07 H₀² per Hubble time, 4–5 % above the
> exact Randall–Sundrum de Sitter brane in AdS₊ (0.607 H₀), as at t = 0.1 (0.649 vs 0.638). This late value depends on
> the resolution of the bulk region behind the receding wall (0.61 with the far-grid spacing 8×10⁻³, 0.64 with 2×10⁻³,
> where the constraint monitor is also 100× smaller); one refinement step is available, so the plateau is quoted to one
> digit and its asymptotic value beyond the chart freeze is not determined. Seeded toward the
> throat, the brane crosses φ_b = −1 at H₀τ = 5.61 with H_J = 0.798 H₀ (4D EFT: 0.797), its tension falls below the BPS
> value beyond φ_b = −1.67, its expansion reverses at H₀τ = 5.89 with φ_b = −1.941 (EFT: −1.942), and it contracts with
> kinetic-dominated scalar runaway; H_J passes −5H₀ at H₀τ = 6.20 (two resolutions within 1.4 %) and −30H₀ before 6.32,
> beyond which the resolutions diverge and the classical evolution ends in a singularity of unresolved form within about
> 0.43 Hubble times of the turnaround. Within the homogeneous, classical, single-scalar, matter-free evolution neither
> branch produces a radiation-dominated brane, at the registered parameters as at t = 0.1 and 0.03. Inhomogeneous
> fragmentation, brane matter, quantum effects and the asymptotic value of the +1 plateau beyond the chart freeze are outside
> this calculation.

## 6. Corrected statements (drop-in replacements)

1. §0: "That caveat is now removed" → "That caveat is now addressed: the code runs at t = 0.001 and reproduces the linear
   theory; the two fates recur, with the caveats below."
2. §0: "reproduces the certified linear growth rate to five digits" → "to 10⁻⁵ in the early windows (1.6572 vs 1.65719)
   and 10⁻⁴ over the full linear window".
3. §0: "No branch produces a radiation-filled universe" → "Within the homogeneous, classical, single-scalar, matter-free
   evolution, no branch produces a radiation-filled universe."
4. §0: "This closes the last open item in the Chat 9–13 sequence" → "This completes the t = 0.001 computation left open
   by Chat 13; the numerics of the new grid and closure await their own referee, and the late +1 plateau awaits a
   coarse-region refinement."
5. §1: "That was half the story" → "That diagnosis was right and accounts for most of the error (rates 5 → 4 → 2.2 → 1.90
   as the wall goes from 6 to 75 points); the remaining 15 % is the second-order closure."
6. §1 table: "1.72 at 150" → "1.7 at 134 points per wall (a t ≤ 3 window only)" or rerun.
7. §1 table, +1 row → "at H₀τ = 6.90: φ_b = 0.9946 / 0.9946, H_J/H₀ = 0.6111 / 0.6118 (75 / 134 points per wall; the
   two runs share the coarse spacing outside |z| ≈ 0.08)".
8. §1: "≲ 10⁻⁴ there throughout" → "≲ 10⁻⁴ at the shell within the quoted windows (8×10⁻⁴ at the end of the coarse throat
   run, 9×10⁻³ in the discarded extended window)".
9. §1: the outgoing-pulse sentence → "For t ≲ 3 the largest bulk constraint violation is the outgoing pulse from the seed's
   junction mismatch, at z = −t. From t ≈ 3.5 it is a feature at z ≈ −0.8 to −1.0, behind the receding wall in the
   coarse-grid region, which grows with the mode and reaches 0.1 at H₀τ = 6.6 and 0.65 at H₀τ = 6.9 (at both shell-layer
   resolutions, which share the coarse spacing). Refining the coarse spacing to 2×10⁻³ reduces it to < 10⁻² and changes
   H_J/H₀ at H₀τ = 6.9 from 0.61 to 0.64."
10. §2: "reaches φ_b = 0.995 by H₀τ ≈ 6.9, and its expansion rate falls … to H_J/H₀ = 0.612 (0.610 at the finer grid) at
    H₀τ = 6.9, still relaxing slowly" → "reaches φ_b = 0.995 near H₀τ = 6.9. Its expansion rate falls monotonically to
    0.63–0.64 H₀ at H₀τ = 6.9–7.0 (dz_c = 2×10⁻³; 0.61 with dz_c = 8×10⁻³, where the bulk constraint monitor is 0.65), still
    decreasing at dH/dτ ≈ −0.07 H₀² when the chart freezes." Delete "The trend across detunings is now clean" and change
    the table row to "≈ 0.63–0.64 (still falling) | 0.6067".
11. §2: "the domain wall, receding into the bulk at light speed, leaves the finely resolved region" → "the domain wall,
    which left the finely resolved layer near φ_b ≈ 0.9 and by t = 8 sits at z ≈ −0.3, recedes toward light speed after
    t ≈ 8".
12. §2 table header "5D value of H_J/H₀ at the freeze" → "5D value at the end of the clean window (t = 0.1: chart
    freeze; t = 0.03 and 0.001: φ_b = 0.995 stop)".
13. §2: "61 % of its original expansion rate" → "about 0.61 of its original expansion rate, consistent with the exact
    value 0.607".
14. §2 throat: "H₀τ = 5.90 with φ_b = −1.953 (… −1.94; the φ = −1 crossing happens at H₀τ = 5.61 with H_J = 0.793 H₀
    against 0.797)" → "H₀τ = 5.89 with φ_b = −1.941 (both resolutions, interpolated between records; 4D EFT: −1.942);
    the φ = −1 crossing is at H₀τ = 5.61 with H_J = 0.798 H₀ (EFT 0.797)".
15. §2 throat: "H_J passes −5H₀ at H₀τ = 6.20 and −300H₀ by 6.32, with the brane's proper time saturating near
    H₀τ ≈ 6.33. … everything up to H_J ≈ −10H₀ is [converged]" → "H_J passes −5H₀ at H₀τ = 6.20 (two resolutions within
    1.4 %) and −30H₀ before 6.32 (the resolutions then differ by 8 % at −10H₀ and 23 % at −22H₀); the coarse run reaches
    −300H₀ at 6.3226 with proper time saturating there. As in Chat 13 the final approach is not converged beyond
    H_J ≈ −5H₀ and its 4D/5D nature is open."
16. §2: "at t = 0.1 the turnaround was at H₀τ = 5.99, φ_b = −1.96, and the collapse followed 0.44 Hubble times later" →
    "… φ_b = −1.94 (interpolated; −1.96 in Chat 13 was the first record past zero), with the extrapolated singularity
    0.44 Hubble times later; here 0.43".
17. §3 "Established (floating point, two resolutions, referee-checked code lineage)" → "Numerically supported (floating
    point; two shell-layer resolutions; Chat 13 code lineage refereed, the v2 grid and closure not yet)".
18. §3 "(bounded above by the quoted value, approaching the RS value from above)" → "(0.63–0.64 H₀ and still falling at
    the chart freeze, one coarse-region refinement step from converged; relaxation to the RS value 0.607 is the expectation,
    which the runs neither confirm nor exclude)".
19. §4 "Chats 9–14 now form a closed, refereed account" → "Chats 9–13 are refereed; Chat 14 extends them to the registered
    detuning and is under review."
20. Figure caption → "Solid: t = 0.001 (75 points per wall). Faint: Chat 13's t = 0.1 runs, shifted by −0.04 (+1) and
    −0.09 (throat) Hubble times so that |φ_b| = 0.05 coincides. The solid +1 curve is cut at H₀τ = 6.9 (φ_b = 0.995 stop,
    lapse collapse, constraint caveat); the faint +1 curve is shown to its own chart freeze at 7.38. Throat curves leave the
    panel at H_J/H₀ = −3."

## 7. Status table (referee's version)

| Claim | Status |
|---|---|
| Code runs at t = 10⁻³; closure-4 linear rate 1.6572 vs 1.65719 at 75 and 134 points per wall | **Established** (logs, ANALYSIS_V2) |
| Closure 2 at the same 75-pt grid gives 1.90; dc = 10⁻⁴ closure-2 run blows up at t = 4.7 | **Established**; the "1.72 at 150" number is a 134-pt, t ≤ 3 transient |
| Both resolution and closure were needed (Chat 13's diagnosis "half the story") | **Numerically supported**, phrase as a refinement of Chat 13's correct diagnosis |
| +1 side: monotone fall of H_J, no oscillation, φ_b → 0.995 at H₀τ = 6.91 | **Established** within the run |
| +1 side: H_J/H₀ = 0.61 at H₀τ = 6.9, "0.610 at the finer grid", "still relaxing slowly", "trend across detunings clean" | **Not established / superseded**: coarse-region truncation error; with dz_c = 2×10⁻³ the value is 0.64 at 6.9 and 0.634 at 6.97 (monitor < 10⁻²), still falling at −0.07 H₀²; one refinement step only |
| +1 side approaches the RS value 0.607 from above | **Not established** (expectation; the extended run crosses it 0.03 Hubble times later and keeps falling; chart degenerating) |
| Throat: φ = −1 crossing (5.61, 0.798), turnaround (5.89, −1.941), −5H₀ at 6.20, all vs 4D EFT | **Established** (two resolutions ≤ 1.4 %; EFT agreement 10⁻³–10⁻⁴ once interpolated) |
| Throat: "converged up to −10H₀", "−300H₀ by 6.32" | **Not established** as stated: 8 % at −10H₀, 23 % at −22H₀; −300 is single-resolution |
| "Same two fates as Chat 13 at the registered parameters; no radiation branch" (homogeneous, classical, one scalar, no matter) | **Numerically supported** |
| "Closes the last open item"; "closed, refereed account" | **Not established**: v2 numerics unrefereed; finecoarse pending; Chat 13's open items stand |

## 8. What I did not do

No runs of my own (the folder's own data sufficed for a claims check). I did not verify the v2 discretisation (chain-rule
metric terms, the closure-4 stencil, `neumann_t`), the meaning of the relative constraint normalisation at z ≈ −1, or the
SVG rendering; those are for the numerics referee. The finecoarse comparison in §2.1 is from log lines up to t = 8.97
(20:59 PDT); when its `_summary.json` and `_timeseries.npz` appear, `analyse_v2.py` will add it to ANALYSIS_V2.json and the
interpolated H_J at H₀τ = 6.90 (expected ≈ 0.638) should replace my log-line reading. Before any plateau value is quoted to
two digits, a second coarse-region refinement (dz_c = 10⁻³, same seed and closure) is required; the README's figure and
tables should then be regenerated from the refined runs, and the 75/134-point runs demoted to a shell-layer convergence
check.
