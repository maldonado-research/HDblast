# Independent audit of the "mechanisms" workstream

Audit date: 28 September 2026. This is an adversarial check of the mechanism screen in this folder (README.md, m1–m8, pilot_5d). The auditor wrote only this file and the `verify/` subfolder. No workstream file was modified. Every audit number comes from a script in `verify/` that writes JSON.

## Verdict in one paragraph

Most of the screen reproduces. Every script re-ran bit-identically. The key identities hold under an independent symbolic derivation (the H² identity with shell radiation, the Weyl transport identity, the Codazzi ledger and the dark-bubble junction). A separate data test also passes: the Weyl transport identity with matter holds on the pilot time series, and it fails when the friction term is dropped.

Two conclusions do not hold as written.

1. **Claim (i), "the c-family ends near c→0⁺; the initial shell cannot be built at c\*", is refuted as stated.**
   - The continuation stops at c=0 only because the solver writes the cone value as φ_h=−1+10^x, which cannot represent φ_h<−1. At c=0 the family crosses the trivial solution φ≡−1; this is exact, since U′(−1)=0 and σ′(−1)=δc.
   - With the sign flipped, the same shell equations and junction residual give a static shell all the way to c\*=−0.99307 at δ=0.1: φ_h=−1.00695, φ_b=−1.7287, H₀²/δ=0.487, residual 3×10⁻¹⁵.
   - Along the continued family f_H falls smoothly to 0, so "f_H≥0.21 along the family" is false.
   - The c\* pilot then runs. Its shell rolls from φ_b=−1.73 across φ=−1 to φ_b=0.80 by the end of the reliable window (H/H₀=0.14).
   - The initial state differs qualitatively from the registered one: the shell starts beyond the φ=−1 vacuum instead of near φ≈0. Whether this roll ends on a vacuum-free brane is unresolved.
2. **The d\* pilot runs were stopped by the solver's φ_b>0.995 cutoff, not by the chart freeze.**
   - At the stop, d(H₀τ)/dt was still 0.31–0.66, above the 0.2 freeze threshold. The tuned "+1" endpoint lies at φ_b≈1.036 (1.0346 at d=−3.0).
   - The "69% radiation share of H²" sits inside a budget with large cancellations. The vacuum part is −117% of H², the Weyl term +121%, radiation +69% and the friction cross terms +23%.
   - With the cutoff lifted, the reliable window extends to H₀τ=3.75. There 𝒲/(radiation terms)=1.23. After the freeze it falls to about 0.6, still ≥6× above the ΔN_eff limit. Then H turns negative, which is a gauge artefact or a reversal.
   - The workstream's "inconclusive" label stands, but for the stated reason it is mis-described.

## Verdict table

| # | Workstream claim | Verdict | Main evidence |
|---|---|---|---|
| 1 | The blast makes a positive Weyl term, peaking at 0.0806–0.0809 H₀² at H₀τ≈6.24 (≈15% of H²); the late value is not converged | **confirmed** (with caveats) | M1 re-run bit-identical. H² identity and Weyl transport identity derived independently (V1 D1, D2). Caveats: the transport residual is 5×10⁻⁸–7×10⁻⁷ at the median, but its 95th percentile is 3–4% (124 in the `ext` run). `ext` continues `plus_c4` (identical peak), so there are 3 independent grids, not 4. |
| 2 | Dark radiation cannot be the hot Big Bang; ΔN_eff≤0.3 (0.09) gives ρ_DR/ρ_SM≤0.105 (0.031) at g\*=106.75 | **partially confirmed** | Mapping recomputed exactly (V4: 0.1050, 0.0315). The N_eff input values (ACT DR6 2.89±0.11; 2.990±0.070) could not be checked, because the web-search budget was exhausted: **unverifiable**. |
| 3 | Dark-bubble cosmology is excluded in the registered model (BPS ΔW=4/3 is critical; brane ρ coefficient −5/54) | **partially confirmed** | All algebra reproduced (V1 D4, D5). ΔW=3(1/ℓ₋−1/ℓ₊) holds for any W with ℓ=3/W, which is the standard statement that BPS walls are critical, not a special feature of this model. The linear Λ₄=10/81−5λ/54 holds only near λ=4/3: at λ=1/2 the exact H² is 0.637, against 0.077 from the linear formula. The exclusion rests on the cited stability theorem and on the Z2 topology, so "exact-verified" overstates it; it is conditional. |
| 4 | The collapse branch is not ekpyrotic (w_eff=0.76–0.85; w>1 for 19–21%) | **partially confirmed** | Numbers reproduce (M5 re-run), but they cover only about 1.05–1.10 e-folds of resolved contraction (−10<H/H₀<−0.5). M1's own median over −30<H<−0.5 is −0.13 to 0.22, with spikes to w≈21. "Negative" holds only on the resolved range. |
| 5 | 4D budget: f_E=0.0858, N_max=0.59 (0.59→0.47 over δ); ≈21.8 e-folds needed | **confirmed** (conditional) | M2 bit-identical. f_E=f_H·M_i²/M_f² rederived from V_E∝H²/M² (V4). The inverted-ratio control gives 1.58. Target 8.78×10³⁷ gives N=21.84. Remains conditional on the 4D EFT and the NEC. Headline wording: f_H is 0.368 only as δ→0; it is 0.372 at δ=0.01 and 0.406 at δ=0.1, so it is not "fixed for every δ". |
| 6 | δ–Λ lock: e-folding time 6.4 Gyr; δ≈1.1×10⁻⁶²…7.6×10⁻⁶⁸ | **confirmed** (conditional) | V4 recomputes 6.417 Gyr, 1.129×10⁻⁶² and 7.58×10⁻⁶⁸. The 38.6 μm AdS-radius proxy was not re-derived. |
| 7 | Sudden tension→radiation conversion is gravitationally screened (𝒲/rad = −0.80…−0.999) | **confirmed** | E5 algebra reproduced. The V1 derivation shows H² depends on the junction only through k_s=(σ+R)/6, so H is continuous at fixed σ+R. |
| 8 | A dissipative coupling at registered c gives negligible radiation (0.37/0.98/1.7/2.2% capture; R/R_crit≤0.038; Y=2 slows the rate 1.605→0.366) | **confirmed** | M4 and M6 re-runs bit-identical. The Y=0 pilot re-run reproduces the archived Chat14 run to 7.6×10⁻¹². Codazzi ledger derived (V1 D3). **New data check (V2):** the Weyl transport identity with matter holds on all pilot series at median 2×10⁻⁵–1.8×10⁻⁴. The controls fail: friction dropped gives 3×10⁻³–5×10⁻², 4H→3H gives 5×10⁻⁴–8×10⁻³. |
| 9 | Gravitational particle production is short by 85–122 orders of magnitude | **confirmed** (order-of-magnitude) | Arithmetic reproduced (5.5×10⁻¹²², 4.1×10⁻⁸⁵). This is essentially the H²/M_P² suppression, with an assumed coefficient. |
| 10 | The registered linear tension cannot tune away the vacuum: the family ends near c→0⁺, f_H≥0.21, and the c\* pilot fails because the shell cannot be built | **refuted** (as stated) | V1 D7, V3, V3b and V5; see above. What survives: the c\* initial shell is a different configuration (φ_b≈−1.73), and its endpoint is unresolved. |
| 11 | A quadratic tension term d\* removes the vacuum and keeps the initial shell; the pilot reaches 69% radiation share with 𝒲/rad=1.76; 𝒲 still falling "when the chart froze"; inconclusive | **partially confirmed** | d\* and initial-shell survival re-run bit-identically (M8). The pilot re-run (dzf=2×10⁻³) is identical. Corrections: (a) the stop was the φ_b cutoff, not the chart freeze; (b) the budget cancellation above; (c) lifting the cutoff reduces 𝒲/rad only to 1.23 inside the reliable window. The M8 series check never reaches d\* itself or δ=0.1; V3 finds the δ=0.1 endpoint consistent with the series to 3% at d=−3.0. |

## Checks run

| Script (in `verify/`) | What it does | Output |
|---|---|---|
| `rerun/` + `compare_json.py` | Copies m1–m8 and re-runs m1, m2, m3, m4, m5, m6 (via a symlink to the saved runs) and m8, then compares every leaf of the JSON. Max relative difference 0 for all seven. m7 was not re-run; its rows are reproduced by V3's control (φ_b=−0.994864 at c=0.002657). m8 took 7.9 min, not the ~2 min stated in the README. | `COMPARE_*.json` |
| `pilot/` | Re-runs `pilot_dstar_Y1…dzf2e-3` (identical summary) and `pilot_Y0` (identical; 7.6×10⁻¹² from the archive). An extension run uses `--phimax 1.3`. | `pilot/*_summary.json` |
| `v1_independent_derivations.py` | Independent SymPy derivation from the SMS projected equations: (D1) H² identity with R; a junction-factor control fails. (D2) Weyl transport with matter; a wrong-ledger-sign control fails. (D3) Codazzi ledger. (D4) BPS tautology and U″(−1)=76/9. (D5) Exact vs linear dark-bubble H². (D6) Leading H_vac². (D7) Trivial solution at c=0. | `V1_INDEPENDENT_DERIVATIONS.json` |
| `v2_pilot_consistency.py` | Weyl-with-matter identity on the pilot data, with three wrong-formula controls. Also the H² budget decomposition at the window end and the a⁴-scaled trends. | `V2_PILOT_CONSISTENCY.json` |
| `v3_c_family_through_zero.py`, `v3b_continue_to_cstar.py` | Signed continuation of the initial-shell family through c=0 to c\* at δ=0.1, plus the tuned "+1" endpoint at δ=0.1. The control reproduces M7's last row to machine precision. | `V3_*.json`, `V3B_*.json` |
| `v4_arithmetic_checks.py` | Recomputes f_E, N_max, the target N, the ΔN_eff mapping, the δ–Λ lock, the bounce H and the Weyl-identity tails. | `V4_ARITHMETIC_CHECKS.json` |
| `pilot/v5_cstar_pilot_signed.py`, `v5_analysis.py` | The c\* pilot with the unmodified pilot code and a signed shell root. Stop window widened to [−3, 1.3]. Single grid, single seed sign, Y=0: exploratory. | `V5_CSTAR_PILOT_SIGNED.json` |

## Other caveats found

- The R_crit used in M6 for the d\* runs belongs to the registered-c vacuum, so R/R_crit=0.15 there has no meaning. The README does not use it.
- Pilot convergence for d\* rests on one grid pair (dz_fine 10⁻³ and 2×10⁻³, which agree to 0.1%). L, dz_coarse and the seed were not varied for d\*.
- The literature list includes two items the workstream says it saw only as search results. The auditor could not open any paper.
- Novelty statements in the README are appropriately limited to "new to this project".

## Reproduce the audit

```
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260927/mechanisms/verify
export OMP_NUM_THREADS=1
(cd rerun && python3 m3_exact_checks.py && python3 m1_weyl_along_chat14.py && python3 m2_energy_budget_and_scales.py \
  && python3 m5_observational_screen.py && python3 m4_dissipative_junction_toy.py && python3 m8_quadratic_tension_tuning.py \
  && python3 m6_pilot_analysis.py)        # rerun/pilot_5d holds a copy of rolloff5d_v1.py and a symlink runs -> ../../../pilot_5d/runs
python3 compare_json.py M1_WEYL_ALONG_CHAT14.json M2_ENERGY_BUDGET_AND_SCALES.json M3_EXACT_CHECKS.json M5_OBSERVATIONAL_SCREEN.json
python3 compare_json.py M4_DISSIPATIVE_JUNCTION_TOY.json; python3 compare_json.py M6_PILOT_COUPLED_5D.json; python3 compare_json.py M8_QUADRATIC_TENSION_TUNING.json
(cd pilot && python3 rolloff5d_matter.py --Y 1 --dc 1e-2 --L 17 --tf 15 --dzf 2e-3 --d -3.106933495673783 --tag rerun_dstar_Y1_dzf2e-3 \
  && python3 rolloff5d_matter.py --Y 0 --tag rerun_pilot_Y0 \
  && python3 rolloff5d_matter.py --Y 1 --dc 1e-2 --L 17 --tf 15 --dzf 2e-3 --d -3.106933495673783 --phimax 1.3 --tag ext_dstar_Y1_dzf2e-3_phimax1.3)
python3 v1_independent_derivations.py; python3 v2_pilot_consistency.py
python3 v3_c_family_through_zero.py; python3 v3b_continue_to_cstar.py; python3 v4_arithmetic_checks.py
(cd pilot && python3 v5_cstar_pilot_signed.py 12); python3 v5_analysis.py
```
