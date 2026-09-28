# HDBLAST mechanism screen: what could turn the five-dimensional event into a hot Big Bang?

Checkpoint workstream, 27–28 September 2026. This is new calculation and literature screening for Ricardo Maldonado's HDBLAST program. It has not been externally reviewed or published. Earlier folders were only read; nothing under `new-files/` was modified.

## Question

Chats 9–14 found two fates for the registered shell. It either relaxes toward an empty Randall–Sundrum (RS) de Sitter brane or reverses and collapses. Neither fate has a radiation era. This workstream asks which additional ingredient could produce a hot, radiation-dominated Big Bang while remaining within Einstein's equations. Each candidate is screened with numbers, not just described.

Registered model (κ₅²=1, model length unit L): W=1−φ+φ³/3, U=½W_φ²−(2/3)W², σ=2W+δ(1+cφ), δ=0.001, c=0.5975949350280132, with a doubled (Z2) bulk and one shell. The bulk vacua are φ=−1 (AdS radius ℓ₋=9/5) and φ=+1 (ℓ₊=9).

## Bottom line

1. **In the registered model, no screened mechanism gives a hot Big Bang.** The decisive obstruction is not the missing matter but the residual RS vacuum.
   - The endpoint keeps H_vac²/H₀²=0.368 of the initial expansion rate squared. This fraction does not depend on δ.
   - With the graviton zero-mode Planck masses, a 4D effective energy budget allows radiation to dominate for at most **0.59 e-folds**. Our universe needs about 22 e-folds.
   - A coupled 5D pilot with shell radiation stays far inside that bound: R/R_crit≤0.04. It shows that the large 5D tension drop is not available as heat on the Hubble time.
2. **The residual vacuum cannot be removed with the registered linear tension.** Setting c≈−1 would do it algebraically, but the initial unstable shell ceases to exist long before that: φ_b→−1 as c→0⁺.
3. **Adding one quadratic tension term does remove it** (d≈−3.19, a model change). The initial unstable shell is essentially unchanged. The pilot then shows, for the first time in this program, shell radiation growing to **69% of H²** by the end of the reliable window.
   - The same run carries a Weyl ("dark radiation") term 1.76× larger than the radiation terms. ΔN_eff allows at most 0.03–0.1.
   - That Weyl term is still decaying (𝒲a⁴ fell 31% in the last 0.06 e-fold) when the coordinate chart freezes. Whether it falls below the ΔN_eff limit is the open question.
   - This is a **conditional lead, not a result**. It needs the vacuum fine-tuned (the cosmological-constant problem restated), and it has so far been run only at δ=0.1, with a phenomenological coupling, up to the chart freeze.

## Method

Every numerical claim comes from a script in this folder. Each script writes a JSON result.

| Script | What it does | Output |
|---|---|---|
| `m1_weyl_along_chat14.py` | Reads 7 archived Chat14 trajectories at δ=10⁻³ without rerunning them. Uses the exact shell identity H²=σ²/36+v²/12−σ′²/48+U/6+𝒲 to extract the Weyl scalar 𝒲 produced by the blast. It also tests the Weyl transport identity on the data, tracks tension flow and measures w_eff during collapse. | `M1_WEYL_ALONG_CHAT14.json` |
| `m2_energy_budget_and_scales.py` | Computes the initial unstable shell and the final "+1" static branch for δ=3×10⁻⁴…0.1, at two tolerances. For each it adds the graviton zero-mode normalisation 1/κ₄²=(2/κ₅²)∫(ρ/ρ_b)²dy. From these it derives the 4D effective energy budget, physical scales, sudden-conversion algebra and a temperature window. | `M2_ENERGY_BUDGET_AND_SCALES.json` |
| `m3_exact_checks.py` | Runs 15 SymPy identities: the BPS structure, critical wall tension, the dark-bubble junction expansion, Weyl conservation, the sudden-conversion identity and the short-wave energy partition. | `M3_EXACT_CHECKS.json` |
| `m4_dissipative_junction_toy.py` | Flat 1D numerical calibration of a dissipative scalar junction, with and without a bulk mass. | `M4_DISSIPATIVE_JUNCTION_TOY.json` |
| `m5_observational_screen.py` | Converts ΔN_eff into a Weyl/Standard-Model radiation ratio. Also runs the collapse-branch ekpyrosis test, the bounce requirement and a gravitational-particle-production estimate. | `M5_OBSERVATIONAL_SCREEN.json` |
| `pilot_5d/rolloff5d_matter.py`, `m6_pilot_analysis.py` | **Coupled 5D pilot.** Chat14's nonlinear solver is extended with shell radiation R=κ₅²ρ (p=R/3) fed by a friction-type scalar source κ₅²j=Yv. All three modified junctions of the 22 Sept matter extension are used: nA=(σ+R)/6, nB=(σ−3R)/6, nφ=−(σ′+Yv)/2. Optional c and d (quadratic tension) settings are supported. Runs are at δ=0.1, the Chat13 detuning, because it is about 100× cheaper. | `M6_PILOT_COUPLED_5D.json`, `pilot_5d/runs/` |
| `m7_initial_shell_vs_c.py` | Follows the initial unstable shell by continuation in c toward c\*≈−1, at δ=0.1, 0.01 and 0.001. | `M7_INITIAL_SHELL_VS_C.json` |
| `m8_quadratic_tension_tuning.py` | Computes d\*(δ) for σ=2W+δ(1+cφ+dφ²/2). Checks the generalised H_vac² series numerically and tests whether the initial shell survives the tuning. | `M8_QUADRATIC_TENSION_TUNING.json` |

## Results

Status labels follow earlier checkpoints: **exact-verified** means symbolic identity; **numerical** means floating-point with measured tolerances; **conditional** means valid under stated assumptions; **negative** means the mechanism fails the screen; **inconclusive** means unresolved.

| # | Candidate | Key numbers | Verdict |
|---|---|---|---|
| (a) | **Bulk Weyl / dark radiation 𝒲∝a⁻⁴ created by the blast** | The blast produces a positive 𝒲. It peaks at 0.0806–0.0809 H₀² at H₀τ≈6.24 in four +1 runs, which is ≈15% of H² at that time. The late value is not converged: 0.026 H₀² on the best grid, versus 0.0004, 0.0018 and −0.37 on the others. ΔN_eff≤0.3 (or 0.09) allows ρ_DR/ρ_SM≤0.105 (or 0.031) at production (g*=106.75). | **Negative** as the source of the hot Big Bang: it is geometry, not a Standard-Model plasma, and it is capped at a few percent of the radiation. **Numerical** that the blast produces it. It is therefore a constraint on any reheating. |
| (b) | **Dark-bubble cosmology** | Exact fake-SUSY structure: φ′=W′ and A′=−W/3 solve the model's static equations. The flat BPS wall tension ΔW=4/3 equals the critical thin-wall value 3(1/ℓ₋−1/ℓ₊)=4/3. In the thin-wall dark-bubble orientation: brane-localised ρ coefficient −5/54, exterior black hole +5/(4a⁴), Λ₄=10/81−5λ/54. | **Negative within the registered model.** The algebra is exact; the no-decay conclusion depends on the cited fake-supergravity stability theorem. No subcritical wall exists for this W, and the mechanism also needs a non-Z2 topology. |
| (c) | **Brane collision / ekpyrotic bounce** (collapse fate) | Resolved contraction (three grids): ln a-weighted mean w_eff=0.76–0.85, with w_eff>1 for only 19–21% of it. A BBN-compatible bounce needs \|H_b\|≥9.4×10¹⁸ H_vac. | **Negative/inconclusive.** The collapse is not ekpyrotic. A bounce or collision partner is outside Einstein's equations plus the NEC. |
| (d) | **Bulk scalar energy absorbed by the shell / dissipative shell coupling** | Exact local capture fraction = Yv/\|σ′\|, which reduces to Y/(2+Y) in the short-wave limit (toy agrees to ≤2.4×10⁻⁴). **5D pilot, registered c:** Y=0.1/0.3/0.6/1 capture 0.37/0.98/1.7/2.2% of the tension drop, and R/R_crit≤0.038 before the chart freezes. Y=2 cuts the instability rate from 1.605 to 0.366 H₀. | **Negative** at registered c: radiation never even momentarily dominates the residual vacuum. |
| (e) | **Gravitational particle production** | ρ~C_g g H⁴ gives ≤5.5×10⁻¹²² (δ fixed by the observed Λ) or ≤4×10⁻⁸⁵ (vacuum tuned, T_max=5 MeV) of the released energy. | **Negative** by 85–122 orders of magnitude. |
| (f) | *New:* **4D energy budget with zero-mode Planck masses** | f_H=H_vac²/H₀²=0.3681, M_i²/M_f²=0.2330 (local RS guess 0.3340), f_E=0.0858. Radiation/vacuum ≤10.7, so N_max=0.59, and 0.59→0.47 across δ=3×10⁻⁴…0.1. H₀²/δ→0.1609 and H_vac²/δ→0.0592=(1+c)/27. | **Conditional (4D EFT, NEC); structural negative for the registered model.** The target is ≈21.8 e-folds (ρ_r/ρ_Λ=8.8×10³⁷ at 5 MeV). |
| (g) | *New:* **δ–Λ–timescale lock** | If H_vac is the observed dark-energy rate, the blast grows with e-folding time 6.4 Gyr, and δ=1.1×10⁻⁶² (ℓ₊=38.6 μm)…7.6×10⁻⁶⁸ (ℓ₊=0.1 μm). | **Conditional.** At registered c the "blast" happens at the dark-energy scale. |
| (h) | *New:* **sudden tension→radiation conversion** | nA is continuous, so H is continuous. The Weyl term must then be 𝒲/(radiation terms)=−0.80…−0.999 (ε=1%…100%). | **Exact-verified.** The 5D tension drop (Δσ=1.33, ~840 R_crit) is gravitationally screened, so tension energy is not free 4D energy. |
| (i) | *New:* **tuning the vacuum with c** | The series gives c\*=−0.99993 at δ=10⁻³. But the initial-shell family runs φ_b→−0.995 as c→0.0027 at δ=0.1, and follows the same curve at δ=0.01 and 10⁻³ (continuation time-limited at c=0.032 and 0.198). Along the family, f_H≥0.21. Direct pilot runs at c\* fail because the initial shell cannot be built. | **Numerical negative.** The registered linear tension cannot remove the vacuum while keeping the blast's initial state. |
| (j) | *New:* **tuning the vacuum with a quadratic tension term** (model change) | d\*=−3.107 (δ=0.1), −3.186 (0.01), −3.194 (10⁻³), close to −2(1+c). The generalised series is confirmed to ≤0.007 δ³ along the path. The initial shell survives: H₀²/δ goes 0.16285→0.16276 at δ=0.1 and is unchanged at 0.01. **Pilot at δ=0.1:** Y=1 reaches a radiation share of H² of 0.690, with 𝒲/(radiation terms)=1.756; two shell spacings agree to 0.1%. Over the last 0.06 e-fold, 𝒲a⁴ fell 31% and Ra⁴ grew 9%. Y=0.3 reaches 0.52 (𝒲/rad=3.46). | **Conditional lead, inconclusive.** Radiation becomes a large share of H², but Weyl radiation is still ~2× the radiation when the chart freezes, while ΔN_eff needs ≤0.03–0.1. Whether 𝒲 keeps decaying is the key open question. It requires the vacuum fine-tuned: \|A\|=\|1+c+d/2\| must match the observed Λ. |

### Details and controls

- **Weyl extraction (a).** The Weyl transport identity is satisfied by the archived data to a median normalised residual of 5×10⁻⁸–7×10⁻⁷. A deliberately wrong factor (4H→3H) gives 1.7–7.9×10⁻⁵. Dropping the σ′²/48 term gives |𝒲₀|>1 instead of ≈5×10⁻⁵. In the frozen-scalar phase 𝒲a⁴ is exactly conserved (M3 E4), so a Weyl term left after the roll behaves as dark radiation.
- **Energy budget (f).**
  - Known limit: the "+1" zero-mode integral matches the analytic dS-in-AdS result [sinh(2ky)/4k−y/2]/sinh²(ky_b) to 1.1×10⁻⁹.
  - Wrong-formula control: an unsquared warp factor gives 8.40 instead of 4.44.
  - Archive checks: ρ_b and H² reproduce the archived values.
  - Tolerance check: tightening changes f_E by ≤2×10⁻⁹.
  - The bound itself uses that total density is non-increasing under the NEC in an expanding FRW phase, so ρ_r,end≤V_i−V_f in the Einstein frame. The ratio to the final vacuum is frame-independent.
- **Dissipative junction (d).**
  - Toy: energy is conserved to ≤1.5×10⁻³, and halving the grid spacing changes the capture fraction by 1.5×10⁻⁴.
  - With a bulk mass m=3, the shell stops at the analytic equilibrium 1−φ²=mφ (0.30278), and the capture fraction saturates near 0.52 instead of 0.8, so the partition depends on the regime.
  - In the registered slow roll (rate ≈1.66 H₀ ≪ 1/ℓ) only a few percent is captured, as the pilot confirms.
- **5D pilot controls.**
  - With Y=0 the extended code reproduces the archived Chat14 run `v2_t01_plus`: |Δφ_b|=7.6×10⁻¹² and |ΔH/H₀|=1.7×10⁻⁹.
  - The radiation ledger dR/dτ+4HR=Yv² holds to a median relative 7×10⁻⁷–1.6×10⁻⁵.
  - Seed size: Y=0.3 with dc=10⁻⁴ and with dc=10⁻² agree at the end of the reliable window to 2.4%.
  - For d\* with Y=1, dz_fine=10⁻³ and 2×10⁻³ agree to 0.1%.
  - All pilot quantities are taken before the conformal chart freezes, where d(H₀τ)/dt<0.2. Values after the freeze (e.g. the summary value H/H₀=0.294 in the dc=10⁻² Y=0 run) are gauge artefacts; they are kept in the files but not used.
- **Failed runs kept.** Both c\* pilots (`pilot_cstar_*.log`) stop at the initial-shell root. The four dc=10⁻² registered-c runs are valid only inside the window.

## Observational constraints used

- **ΔN_eff.** ACT DR6 gives N_eff=2.89±0.11 ([arXiv:2503.14454](https://arxiv.org/pdf/2503.14454)). A 2026 combined BBN+CMB+BAO analysis gives 2.990±0.070 ([arXiv:2603.13226](https://arxiv.org/html/2603.13226v2)). Braneworld dark radiation at BBN is limited to −12.1%…+6.2% of the photon background ([Ichiki et al.](https://arxiv.org/abs/astro-ph/0203272); also [astro-ph/0208133](https://arxiv.org/pdf/astro-ph/0208133)). We use ΔN_eff≤0.3 (BBN-level) and ≤0.09 (combined). Only the expansion-rate part transfers directly to a Weyl fluid.
- **AdS radius.** No Yukawa deviation from Newtonian gravity at strength α=1 is seen above 38.6 μm ([Eöt-Wash, arXiv:2002.11761](https://arxiv.org/pdf/2002.11761)). This is used as a proxy upper bound on ℓ₊; the RS-specific bound was not re-derived.
- **Hot Big Bang target.** T_RH≥5 MeV, and ρ_r/ρ_Λ at 5 MeV is 8.8×10³⁷, equivalent to N≈21.8 e-folds. Inputs: h=0.674, Ω_Λ=0.685, g*=10.75.

## Literature (URLs obtained in this session, or cited by the read-only packages)

- Dark bubble: [Banerjee et al. 2018](https://arxiv.org/abs/1807.01570); [The dark bubbleography (JHEP 2024)](https://link.springer.com/article/10.1007/JHEP02(2024)102); [Shedding light on dark bubble cosmology](https://arxiv.org/pdf/2310.15032); [Experimental tests of dark bubble cosmology (PRD 109 026003)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.109.026003); [Weak gravity at micron scales from dark bubble cosmology (2025)](https://arxiv.org/html/2511.21362); [Dark bubbles, dark dimensions and fat gravitons (2026)](https://arxiv.org/pdf/2606.20942). The last two were seen only as search results.
- Brane–bulk exchange and dark radiation: [Kiritsis et al., hep-th/0207060](https://arxiv.org/abs/hep-th/0207060); [Langlois, Sorbo, Rodríguez-Martínez, hep-th/0206146](https://arxiv.org/pdf/hep-th/0206146).
- Fake supergravity and domain-wall stability: [Freedman, Núñez, Schnabl, Skenderis, PRD 69 104027](https://ui.adsabs.harvard.edu/abs/2004PhRvD..69j4027F/abstract).
- Collisions and ekpyrosis: [Langlois, Maeda, Wands, gr-qc/0111013](https://arxiv.org/pdf/gr-qc/0111013); [hep-th/0111279](https://arxiv.org/pdf/hep-th/0111279); [Khoury, Ovrut, Steinhardt, Turok, PRD 64 123522](https://harvest.aps.org/v2/journals/articles/10.1103/PhysRevD.64.123522/fulltext); [Reheating after S-brane ekpyrosis, PRD 102 063514](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.102.063514).
- Gravitational particle production: [Kolb and Long, RMP 96 045005 / arXiv:2312.09042](https://arxiv.org/abs/2312.09042). The C_g≈10⁻² coefficient is our order-of-magnitude assumption.
- RS2: [hep-th/9906064](https://arxiv.org/abs/hep-th/9906064v1). Instant preheating and BraneCode, as cited in the earlier packages: [hep-ph/9812289](https://arxiv.org/abs/hep-ph/9812289), [hep-th/0309001](https://arxiv.org/abs/hep-th/0309001).

Dark radiation, dark bubbles, ekpyrosis, gravitational particle production and fine-tuning of the brane tension are established ideas. The following are new to this project only, with no claim of novelty beyond it: the Weyl extraction from the HDBLAST blast; the zero-mode energy budget and its δ-independence; the δ–Λ lock; the c-family obstruction; the quadratic-tension tuning that preserves the initial shell; and the coupled pilot.

## Limitations

- The 4D budget (f) assumes that a 4D scalar-tensor description with zero-mode Planck masses captures the energetics, and that the NEC holds. Chat 9's 4D theory matched the 5D crossing and turnaround to 3 digits, which supports this but does not prove it.
- The pilot runs, like Chat14, lose the conformal chart shortly after φ_b≈0.99. Asymptotic outcomes, including the late Weyl amplitude, are therefore not established.
- The pilot matter sector is a phenomenological friction closure, not a derived quantum particle-production calculation.
  - It was run only at δ=0.1.
  - Convergence evidence is one two-grid pair (d\*, Y=1) plus the Y=0 archival match and the seed-size check.
  - `neumann_t`, which is used only by the weak Kreiss–Oliger boundary ghosts, omits the Y·dv/dt and dR/dt pieces.
- c\* and d\* use the O(δ²) series. For d\*, the series is checked to 0.007 δ³ up to 90% of d\*. The residual vacuum at the pilot's d\* is O(δ³), not exactly zero.
- The c-continuation at δ=0.01 and 10⁻³ stopped on its time limit, not on a failure. The inference that the family ends near c→0⁺ rests on the δ=0.1 continuation and on φ_b(c) agreeing across δ to 3 decimals.
- The ΔN_eff mapping uses entropy-conserving g* scaling. Gravitational-particle-production and bounce numbers are order-of-magnitude.
- The literature search was targeted, not exhaustive.

## The single most promising next calculation

**Follow the vacuum-tuned (d\*) coupled roll-off past the conformal-chart freeze and measure the asymptotic Weyl-to-radiation ratio.**

1. **Model.** σ=2W+δ(1+cφ+dφ²/2) with registered c and d=d\*(δ) from `m8_quadratic_tension_tuning.py`. This is an explicit model change.
2. **Code.**
   - Start from `pilot_5d/rolloff5d_matter.py`, which already has the three modified junctions, the radiation ODE, the ledger checks and the `--d` option.
   - Replace the conformal time slicing near the shell with the proper-clock gauge proved in folder 152 (`literature/PROPER_CLOCK_PROOF.md`), so the evolution continues past H₀τ≈3.6.
   - Add dR/dt and the source's time derivative to `neumann_t`.
3. **Runs.**
   - (i) δ=0.1, Y∈{0.3, 1, 3}, dc=10⁻², L=17, two shell spacings.
   - (ii) The same at δ=10⁻³ with the Chat14 refined grid: `--tdet 1e-3 --closure 4 --dzf 1.5e-4 --dzc 2e-3 --zfine 0.06 --L 10`.
4. **Measure.**
   - Radiation share Ω_r(τ).
   - 𝒲/(radiation terms), together with 𝒲a⁴ and Ra⁴, until both plateau.
   - Constraint residuals.
5. **Decision rule.**
   - If 𝒲a⁴ plateaus with 𝒲/(radiation terms)≤0.1 (ΔN_eff≤0.3; ≤0.03 for the combined bound) while Ω_r→1, the tuned model yields a radiation era compatible with the dark-radiation constraint. Then replace the friction with the derived instant-preheating χ source of the 22 Sept extension and compute T_max(δ, ℓ₊). With the vacuum tuned, T≥5 MeV needs δ≳10⁻²⁶ at ℓ₊=38.6 μm.
   - If 𝒲 plateaus above that ratio, the tuned model is excluded by N_eff, and no screened mechanism survives.

Two things remain even if the test passes. Perturbations, meaning the primordial spectrum, are untouched. And the fine-tuning of 1+c+d/2 against the observed Λ is an assumption, not an explanation.

## Reproduction

The system python3 has numpy 2.3.5, scipy 1.16.3 and sympy 1.14. Run from this folder with `export OMP_NUM_THREADS=1`:

```
python3 m1_weyl_along_chat14.py            # ~5 s (reads the archived Chat14 zip, read-only)
python3 m2_energy_budget_and_scales.py     # ~10 s
python3 m3_exact_checks.py                 # ~5 s, asserts 15 identities
python3 m4_dissipative_junction_toy.py     # ~3 min
python3 m5_observational_screen.py         # ~2 s (needs M1, M2)
python3 m7_initial_shell_vs_c.py           # ~12 min (240 s continuation budget per delta)
python3 m8_quadratic_tension_tuning.py     # ~2 min
cd pilot_5d                                 # each run ~5-10 min on one core
python3 rolloff5d_matter.py --Y 0   --tag runs/pilot_Y0
python3 rolloff5d_matter.py --Y 2   --tag runs/pilot_Y2
python3 rolloff5d_matter.py --Y 0.1 --L 21 --tf 18 --tag runs/pilot_Y0.1_L21
python3 rolloff5d_matter.py --Y 0.3 --L 21 --tf 18 --tag runs/pilot_Y0.3_L21
python3 rolloff5d_matter.py --Y 0   --dc 1e-2 --tag runs/pilot_Y0_dc1e-2
python3 rolloff5d_matter.py --Y 0.3 --dc 1e-2 --tag runs/pilot_Y0.3_dc1e-2
python3 rolloff5d_matter.py --Y 0.6 --dc 1e-2 --L 17 --tf 15 --tag runs/pilot_Y0.6_dc1e-2_L17
python3 rolloff5d_matter.py --Y 1   --dc 1e-2 --L 17 --tf 15 --tag runs/pilot_Y1_dc1e-2_L17
python3 rolloff5d_matter.py --Y 0   --dc 1e-2 --L 17 --tf 15 --c -0.9930694751920546 --tag runs/pilot_cstar_Y0_dc1e-2_L17    # fails (documented)
python3 rolloff5d_matter.py --Y 0.3 --dc 1e-2 --L 17 --tf 15 --c -0.9930694751920546 --tag runs/pilot_cstar_Y0.3_dc1e-2_L17  # fails (documented)
python3 rolloff5d_matter.py --Y 0   --dc 1e-2 --L 17 --tf 15 --d -3.106933495673783 --tag runs/pilot_dstar_Y0_dc1e-2_L17
python3 rolloff5d_matter.py --Y 0.3 --dc 1e-2 --L 17 --tf 15 --d -3.106933495673783 --tag runs/pilot_dstar_Y0.3_dc1e-2_L17
python3 rolloff5d_matter.py --Y 1   --dc 1e-2 --L 17 --tf 15 --d -3.106933495673783 --tag runs/pilot_dstar_Y1_dc1e-2_L17
python3 rolloff5d_matter.py --Y 1   --dc 1e-2 --L 17 --tf 15 --dzf 2e-3 --d -3.106933495673783 --tag runs/pilot_dstar_Y1_dc1e-2_L17_dzf2e-3
cd .. && python3 m6_pilot_analysis.py
```

The inputs are read from the 22 Sept checkpoint: `frozen/registered_solver.py`, `static_branch/solve_plus_branch.py` (imported, not modified) and `source_audit/inputs/CHAT14_COMPLETED_REFERENCE.zip`. `pilot_5d/rolloff5d_v1.py` is an unmodified copy of the Chat14 model module (sha256 9ea4d8d6…). The extended solver derives from Chat14's `rolloff5d_v2.py` (sha256 7dc0cd58…), with the changes listed in its docstring. The pilot default z_fine=0.4 reproduces the archived δ=0.1 grid (1833 points).
