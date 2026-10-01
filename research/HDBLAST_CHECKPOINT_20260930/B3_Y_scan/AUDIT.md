# Independent audit of B3 (r(Y) scan, tuned model d = d*(0.1), δ = 0.1, friction closure)

Audit date: 1 October 2026. Auditor wrote only this file and `audit/` (scripts, their JSON outputs, a copy of the solver and one re-run).
Nothing else in the producer folder was modified. Cores used: 1.

## What was checked and how

| Check | Script / output | Method |
|---|---|---|
| A. Registration timing and integrity | file birth/modify times (`stat`), text of REGISTRATION.md | REGISTRATION.md born 10:11:20; first scan job list (`jobs_main.txt`, `jobs_pre.txt`) born 10:11:30; first pre-run finished 10:12:26; first main run finished 10:18:47. Only the D1 diagnosis (10:04), the dev run (10:03–10:09) and `diag/` (10:01–10:10) predate it, all disclosed in §0. |
| B. Code provenance | `diff` against A1 originals, `sha256sum` | `a1_analyze.py`, `static_w.py` identical to A1. `evolve_a1.py` differs only by default-off options (`--save_T`, `--restart`, `--stop_recollapse`, `--order`); the hash matches `code_sha256` stored in the run summaries. |
| C. Independent re-classification | `audit/audit_b3.py` → `audit/AUDIT_CHECKS.json` | Own loader, own hybrid splice, own plateau finder (written from the rule text, not the producer's code), run on every B3 and A1 run. |
| D. Cross-resolution r at the decisive plateau times | same | r of every resolution interpolated to the finest run's plateau time. |
| E. Common-ancestor error of hybrid runs | same | 1e-3 vs 5e-4 difference at the restart time (the hybrids inherit the 5e-4 history before it). |
| F. Fit and trend | same | Re-fit of the power law; local log-log crossing of r = 0.03. |
| G. Ledger identity R a⁴ = Y ∫v²a⁴dτ | same | Recomputed by trapezoid from saved series. |
| H. C10, C11, D2 | `audit/audit_c10_c11_d2.py` → `audit/AUDIT_C10_C11_D2.json` | C10 with an independent settled reference; C11 extended to H₀τ = 30; D2 eigenproblem re-run and the full RK4 amplification factor evaluated on every eigenvalue (not only the real-axis limit). |
| I. Reproducibility of a decisive PDE run | `audit/rerun/`, `audit/compare_rerun.py` → `audit/AUDIT_RERUN.json` | Fresh re-run of `main_Y2_dc1e-2_dzf5e-4` (parent of the decisive Y = 2 hybrid) with the same code and arguments (without `--save_T`). |

No new resolutions were run (one core; a 2.5e-4 or 1.25e-4 run costs 3–11 min of wall time but needs a parent state; the decisive
comparison was instead made across the four existing resolution levels).

## Findings by claim

1. **Y = 3 diagnosis (A1 no-plateau = 4th-order drift, lapse-amplified, ~0.6 e-folds window): partially confirmed.**
   - Confirmed: at the Y = 3 plateau time (H₀τ = 23.55) the relative W a⁴ offsets from the 1.25e-4 run are 1.29 (1e-3), 0.087 (5e-4), 0.0051 (2.5e-4); successive ratios 14.8 and 16.9, i.e. order 3.9–4.1 (audit, independent of D1). D1's 14–23 / 3.8–4.5 reproduced in spirit. Settled-to-turnaround 0.57–0.64 e-folds and late-window Weyl identity p95 0.32–0.36 are in `D1_Y3_DIAGNOSIS.json`.
   - Not demonstrated: "amplified by the bounded chart's growing shell lapse" is a correlation (drift grows while the lapse grows); no run with a different chart was made to show causation.
   - Inconsistencies in the text: the `conclusions` strings in `D1_Y3_DIAGNOSIS.json` say "~0.75 e-folds" between friction-off and turnaround and "lapse 10 -> 50", while the computed fields give 0.995–1.016 e-folds (friction-off → turnaround) and 0.57–0.64 (settling → turnaround). REGISTRATION §0 quotes Ω_r ≈ 1.05, Ω_vac ≈ −0.06 at settling; the D1 JSON gives Ω_r 1.09–1.12, Ω_vac −0.10 to −0.13. These do not affect any class.

2. **Y = 5 stiffness (D2): confirmed, with a qualification.** Re-run reproduces |λ|dz = 2.31 (Y = 2), 5.03 (Y = 3), 9.76 (Y = 5). Evaluating the RK4 amplification factor on all eigenvalues: max |amp| = 12.3 at Y = 5, cfl 0.5 (unstable), 1.0 at cfl 0.25 and at Y = 2, 3. Qualification: a large real (friction) eigenvalue exists only for Y ≥ 3; at Y ≤ 2 the largest |λ| is the interior wave mode (≈ 2.3–2.9, imaginary), so "|λ|dz grows ~ Y" holds only for Y ≳ 3.

3. **Y = 0.5, 0.7 FAIL-Weyl; Y = 1.5 PASS-conservative: confirmed.** Audit plateau finder gives identical plateau times and r (0.2237, 0.1326, 0.0360); pairs agree to < 3%; near-shell constraints < 1e-4.

4. **Y = 2 PASS-combined, both seeds (headline): confirmed (numerical).**
   - Audit plateau finder: 5e-4 run r = 0.02148 (H₀τ 18.20), hybrid 2.5e-4 r = 0.02130 (same time); dc = 1e-4: 0.02147 / 0.02130 at H₀τ 21.66. Same class in both runs of each pair.
   - Robustness beyond the rule: at the plateau time the 1e-3 run (no plateau) has r = 0.0242, the 5e-4 run 0.02148, the 2.5e-4 hybrid 0.02130; differences shrink by 15.5×, i.e. 4th-order convergence; every level is below 0.03. Richardson 0.02128–0.02129 confirmed from JSON.
   - Common-ancestor error (hybrid shares the 5e-4 history up to H₀τ 5.69 / 9.17): 1e-3 vs 5e-4 differ there by ≤ 2e-4 in W and ≤ 5e-5 in R (estimated 5e-4 error ~1e-5), so the shared history cannot affect the class. C10 (restart from the coarse run reproduces the full 5e-4 drift; audit recomputation 0.5–4% agreement) also supports the restart method.
   - Fresh re-run of the 5e-4 parent: see "Re-run" below.
   - Caveat (correctly disclosed by the producer, but material): the PASS sits at H/H₀ ≈ 0.03, with Ω_r = 1.36 only because Ω_vac = −0.39, and the shell turns around 0.32 e-folds later. Ω_r ≥ 0.9 is met because the negative vacuum lowers H², not because radiation dominates a long expansion. Also note that the A1 r-agreement tolerance max(20%, 0.02 absolute) is 100% of r at r ≈ 0.02; it is the "same class" requirement, not the r tolerance, that does the work near the 0.03 threshold.

5. **Y = 3 PASS-combined (tertiary pair), labelled conditional: label appropriate; value confirmed; class depends on a post-hoc level.**
   - Audit finder: 2.5e-4 hybrid r = 0.01017, 1.25e-4 hybrid r = 0.01012, both plateau at H₀τ 23.55. The value r ≈ 0.010 is robust (even the 1e-3 run gives r = 0.023 < 0.03 at that time); the PASS class is not: the 5e-4 run's minimum 0.5-e-fold spread is 6.2% (needs 5%).
   - The tertiary level was added after the 5e-4 / 2.5e-4 outcome was seen, and the Y = 3 r value was known before registration (both disclosed).
   - `TERTIARY_ORDER.json`: the quoted order 3.8–4.0 is measured at H₀τ = 43.8 and 52.4, near/after the turnaround (H = 0 at 47.7). At the plateau-relevant times (H₀τ 20–26) the 1.25e-4 drift is −0.09 to −0.15% (opposite sign to the coarser drifts) and the apparent order is undefined. The magnitude (0.15%) is far below the 5% tolerance, so this does not threaten the class, but "order 3.8–4.0" does not describe the plateau epoch. The dc = 1e-4 repeat was not run; the YES therefore rests on Y = 2, as the producer states.
   - Both tertiary runs restart from the same 5e-4 state and keep dz_coarse = 1e-3, so the tertiary pair does not test the coarse-region or pre-restart error; E shows the latter is ~1e-5.

6. **Y = 5 UNRELIABLE / inconclusive: confirmed.** Only the 1.25e-4 hybrid has a plateau (r = 0.0044); the 2.5e-4 hybrid's best spread is 10.6%; r at the plateau time ranges 0.0044 (1.25e-4), 0.0054 (2.5e-4), 0.020 (5e-4); W a⁴ differs by 23% between the two finest. Not converged.

7. **Trend r ~ 0.069 Y^−1.67, Y(r = 0.03) ≈ 1.64: confirmed (descriptive).** Re-fit: p = 1.670, a = 0.0685, max deviation 10.6%; local log-log interpolation between Y = 1.5 and 2 gives Y = 1.66; fit without the conditional Y = 3 point: p = 1.63. Local slopes 1.42 → 1.84 confirmed, i.e. not a pure power law. The ledger identity R a⁴ = Y ∫v²a⁴dτ/ρ_b is reproduced to ≤ 1.1e-4 on every plateau run. The "mechanism" (Weyl yield per ∫v²a⁴ falls as Y^−0.67) is a re-expression of the fitted r(Y) through that identity, not an independent explanation; the in-words physical story is an interpretation.

8. **Plateau later / fewer e-folds as Y grows: partially confirmed.** Measured e-folds plateau → turnaround: 0.383 (1.5), 0.319 (2), 0.227 (3), 0.093 (5) — match the producer and the ¼ ln((Ω_r+Ω_W)/|Ω_vac|) estimate to ≤ 1%. For Y = 0.5, 0.7 and Y = 2 at dc = 1e-4 the runs stopped before the turnaround (H stays > 0.012–0.03; Y = 0.5/0.7 end "non-finite" while still expanding), so the 0.50 / 0.48 / 0.32 values quoted there are the estimate only, not measurements (the audit's within-run lower bounds are 0.17, 0.26, 0.23).

9. **Strengthened C5: confirmed numbers; discriminating, with a note.** Windowed ratios (R dropped 40–675×; 3R / −R 75–1317×) and 4H→3H 124–2023× match the JSON. Correction: the summary says the all-record R-dropped control fails "for Y ≤ 1.5"; it also fails at Y = 2, dc = 1e-4 (8.6×) — the README states this correctly. The window (R-term weight ≥ 1e-3) lies in the friction-on epoch (H₀τ ≈ 2–12), so it tests the R coefficient during the roll, not in the late plateau epoch; it was chosen after viewing A1 values (disclosed).

10. **C10 restart validation and C11 L-independence: confirmed.** C10 recomputed with an independent reference: restart vs full 5e-4 drift 0.0051/0.0053, 0.0191/0.0194, 0.0468/0.0470, 0.0957/0.0960 at H₀τ 18–24, coarse 0.15–1.48 (registered criteria met). C11 extended to H₀τ = 30: max relative ΔW a⁴ ≤ 1.9e-7, ΔH ≤ 6e-11 (L = 10 vs 16 effectively identical, as expected from causality at T < 20).

11. **C6 ledger: confirmed.** Medians 2.7e-6 to 1.07e-4 on the finest runs; Y = 0.5 fails marginally (1.07e-4; 1.13e-4 on its 1e-3 run), disclosed.

## Registration

- Written before any scan run (file times above); disclosures in §0 match the files that predate it (D1, dev run, diag/). Whether the pre-note text was ever edited cannot be proven (no version history was consulted; the file's modification time reflects the appended notes); nothing in the text looks rewritten, and the dated notes are ordered consistently with the run file times.
- Deviations, all disclosed in dated notes: cfl 0.25 for Y = 5; trigger evaluated at the last common reliable time; speculative dc = 1e-4 scheduling; same-grid extension for Y = 5; tertiary level for Y = 3 and Y = 5 (post hoc).
- Minor undisclosed items: `secondary_trigger.py` and `tertiary_order.py` were modified at 12:15 (after all runs) and their JSON regenerated at 12:32; the trigger outcomes (orders 3.78–4.15) are far from the 3/5 bounds, so this cannot have changed any decision. REGISTRATION §0 mentions a "dev check of `b3_analyze.py`" before 10:11, but the present `b3_analyze.py` was created at 12:15 (an earlier version is not preserved).
- The registered per-run rule and reliability criterion are applied unchanged (identical code; audit re-implementation agrees on every run).

## Re-run

Fresh re-run of `main_Y2_dc1e-2_dzf5e-4` (436 s wall, same code hash): all 743 records are bit-identical to the saved run (max |Δ| = 0 in H₀τ, H, 𝒲, R, φ_b), same stop reason, and the audit's plateau finder gives the same plateau (H₀τ 18.20, r = 0.021475, Ω_r = 1.357, Ω_vac = −0.386). This shows the saved outputs come from the stated code and arguments, and that the result is deterministic. It does not test convergence; that is covered by D above.

## Overall

The headline "a Y with r ≤ 0.03 exists within the friction closure at δ = 0.1, d = d*" is supported at Y = 2 for both seeds by the registered secondary pair, and the r value (0.0213) is robust across four resolution levels with 4th-order convergence. The Y = 3 PASS is correctly labelled conditional. Status labels are not stronger than the evidence. Remaining corrections are textual (D1 conclusion strings, Ω values in §0, the C5 "Y ≤ 1.5" wording, e-folds for Y = 0.5/0.7/2(dc 1e-4) being estimates, the lapse causation being correlational, tertiary order measured near turnaround). The physical significance is limited exactly as the producer says: the "radiation era" is a ≤ 0.5-e-fold interval ending in recollapse, Ω_r > 1 reflects the negative residual vacuum, Y is phenomenological, and d is fine-tuned.

## Reproduction

```
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260930/B3_Y_scan/audit
export OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
python3 audit_b3.py              # AUDIT_CHECKS.json
python3 audit_c10_c11_d2.py      # AUDIT_C10_C11_D2.json
python3 -W ignore evolve_a1.py --out $PWD/rerun --tag rerun_Y2_dc1e-2_dzf5e-4 --dstar --Y 2 --dc 1e-2 --dzf 5e-4 --dzc 2e-3 --L 10 --Tf 20 --kappa 10 --project 9.9 --stop_recollapse --xc 4.6
python3 compare_rerun.py         # AUDIT_RERUN.json
```
