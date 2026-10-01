# A1 pre-registration: tuned-vacuum roll-off past the chart freeze

Registered 28 September 2026, before any tuned-model (d = d\*) evolution with the new solver was run.
Ricardo Maldonado's HDBLAST program; prepared with AI assistance. This file is not edited after registration; changes are appended below as dated notes with reasons.

Work done before this registration (none of it is a tuned-model run): exact symbolic checks (`S0_SYMBOLIC_CHECKS.json`), the static-shell/reference controls (`t1_static_controls.py`), the diagnosis of the chart freeze from archived pilot data (`D0_FREEZE_DIAGNOSIS.json`), and two short development runs of the new solver in the old chart with the registered tension (`dev/smoke_*`, compared with the archived pilot; they are development checks, not the formal calibration below).

## 1. Model (MODEL CHANGE, labelled)

- Registered 5D Einstein-scalar model, κ₅² = 1: W = 1 − φ + φ³/3, U = ½W_φ² − (2/3)W², Z2-doubled bulk, one shell.
- **Model change for A1:** σ(φ) = 2W + δ(1 + cφ + dφ²/2), with the registered c = 0.5975949350280132 and d = d\*(δ) from `HDBLAST_CHECKPOINT_20260927/mechanisms/M8_QUADRATIC_TENSION_TUNING.json`: d\*(0.1) = −3.106933495673783, d\*(10⁻³) = −3.1942416958680835. (d = 0 is the registered model and is used only for calibration.)
- Shell matter (22 Sept matter extension, as in the pilot): radiation R = κ₅²ρ, p = R/3, fed by the phenomenological friction closure κ₅²j = Y v (v = dφ_b/dτ). Junctions: n·∂A = (σ+R)/6, n·∂B = (σ−3R)/6, n·∂φ = −(σ′+Yv)/2; ledger dR/dτ + 4HR = Yv².
- Seed: the pilot's seed (static shell with c + dc, same d), dc ∈ {10⁻², 10⁻⁴}, positive (the +1 side).

## 2. Dark-radiation normalisation (fixed here)

**Convention.** All quantities are evaluated on the shell from the exact shell identity (as in the pilot and in the V2 audit):

H² = (σ+R)²/36 + v²/12 − (σ′+Yv)²/48 + U/6 + 𝒲.

- "Radiation terms": rad ≡ σR/18 + R²/36 (the whole contribution of R to H²).
- Weyl (dark-radiation) term: 𝒲 as defined by the identity.
- **Ratio: r ≡ 𝒲 / rad**, evaluated at the time the radiation era would begin, i.e. on the plateau defined in §3 (after the friction source has switched off, where 𝒲a⁴ and Ra⁴ are both constant, so r is constant thereafter).
- Radiation share Ω_r ≡ rad/H².

**Mapping to ΔN_eff (exact algebra, `S0_SYMBOLIC_CHECKS.json`).** In the low-energy regime (R ≪ σ, where rad → σR/18 = (8πG₄/3)ρ_r) both 𝒲 and rad scale as a⁻⁴, so r = ρ_DR/ρ_r at production. If the brane radiation is the Standard-Model plasma with g\* = g\*_s = 106.75 at production and entropy is conserved afterwards, ρ_SM a⁴ ∝ g\*^{−1/3}, so

- ρ_DR/ρ_SM(BBN, T ≈ 1 MeV, g\* = 10.75) = r × (10.75/106.75)^{1/3} = 0.46524 r, and ΔN_eff = ρ_DR/ρ_SM(BBN)/(7/43) = r/0.349904;
- after e⁺e⁻ annihilation: ρ_DR/ρ_γ = (7/8)(4/11)^{4/3} ΔN_eff = 0.227107 ΔN_eff.

**Thresholds adopted (from NEXT_TESTS A1):**

| r (at production) | ΔN_eff | ρ_DR/ρ_SM at BBN | ρ_DR/ρ_γ at CMB |
|---|---|---|---|
| 0.1 (conservative, BBN-level) | 0.286 | 0.0465 | 0.0649 |
| 0.03 (combined CMB+BBN+BAO) | 0.0857 | 0.0140 | 0.0195 |

The ΔN_eff input values (≈0.3 and ≈0.09) are snippet-level (`LITERATURE_2022_2026.md`); they were not verified against the papers here (paper sites are blocked). The ACT DR6 snippet ρ_dr/ρ_γ < 0.016 corresponds to r ≈ 0.025, between the two thresholds. The mapping counts only the expansion-rate effect of the Weyl fluid; perturbation-level differences from free-streaming radiation are ignored. **Sign:** NEXT_TESTS states the rule one-sided (r ≤ 0.1). A negative Weyl term lowers N_eff and is also constrained, so this registration uses |r|. If a run ends with r < −0.1, it is reported under both readings.

The R²/36 correction is included in rad. R/σ is reported at the plateau; if R/σ > 0.1 the low-energy identification r = ρ_DR/ρ_r is flagged as approximate.

## 3. Pass/fail rule (fixed now)

**Plateau.** Let ln a = A_b − A_b(0) (shell scale factor). A run reaches a plateau at time τ_p if, over the last Δln a = 0.5 before τ_p, H > 0 throughout and both |Δ(𝒲a⁴)|/|𝒲a⁴| ≤ 0.05 and |Δ(Ra⁴)|/(Ra⁴) ≤ 0.05. The run's r and Ω_r are those at τ_p (the first time the plateau criterion holds).

**Per run:**
- **PASS-conservative:** plateau reached, Ω_r ≥ 0.9 at τ_p, and |r| ≤ 0.1. **PASS-combined** additionally requires |r| ≤ 0.03.
- **FAIL-Weyl:** plateau reached with |r| > 0.1.
- **FAIL-no-radiation-era:** (a) the shell recollapses (H/H₀ < −0.05) before any time with Ω_r ≥ 0.9 and |r| ≤ 0.1, or (b) plateau reached with Ω_r < 0.5 (vacuum or Weyl dominated).
- **INCONCLUSIVE:** none of the above within the reliable evolution (end time, wall-clock limit of 25 min per run, or loss of reliability); or plateau with 0.5 ≤ Ω_r < 0.9 and |r| ≤ 0.1.

**Reliability of a run's classification.** A classification counts only if (i) the two shell spacings (dz_fine = 10⁻³ and 5×10⁻⁴, coarse spacing refined by the same factor) give the same class and agree on r at the classification time to 20% relative or 0.02 absolute, whichever is larger (for a recollapse: agree on the proper time of H = 0 to 0.05 H₀⁻¹); (ii) the relative Hamiltonian and momentum residuals within Z > −1 of the shell stay below 0.05 up to the classification time on the finer grid.

**Aggregate verdict for A1 at δ = 0.1** (Y ∈ {0.3, 1, 3}, dc ∈ {10⁻², 10⁻⁴}, two spacings):
- **PASS** (tuned model has a candidate radiation era at δ = 0.1) if, for at least one Y, both seeds give a reliable PASS-conservative. Report PASS-combined separately.
- **FAIL** (d\* model excluded at δ = 0.1 within the friction closure for the scanned Y) if every Y gives a reliable FAIL (either kind) for both seeds.
- **INCONCLUSIVE** otherwise.
- δ = 10⁻³ runs, if completed, are reported with the same per-run rule; if they disagree with δ = 0.1, the δ = 10⁻³ result governs statements about the registered-scale detuning and the aggregate becomes INCONCLUSIVE unless both are FAIL or both PASS.

**Meaning of the outcomes.** PASS: the d\* model has a candidate radiation era compatible with the N_eff bound in this closure; the next step is to replace the friction closure with the derived χ source (and the tuning 1 + c + d/2 ≈ 0 remains an unexplained fine tuning, the cosmological-constant problem restated). FAIL: the d\* model is excluded by N_eff or has no radiation era; together with the 27 Sept screen, no screened mechanism survives. INCONCLUSIVE: neither can be claimed.

## 4. Method (fixed now; details in README)

- Chart: conformal gauge with the shell fixed at Z = 0, in null labels U = F(u), V = F(v) with one common F (PROPER_CLOCK_PROOF.md setting), F(x) = ln(1+e^{−x_c}) − ln(e^{−x} + e^{−x_c}). Reason (`D0_FREEZE_DIAGNOSIS.json`): in the old chart d B_b/dt → −1.000, i.e. the shell reaches the future light cone of the static vertex at finite proper time; F cancels the e^{−t} decay and is regular across that light cone. x_c is a gauge choice, set from a coarse old-chart pre-run as b_b + t at the first time d b_b/dt < −0.9 (rounded to 0.1; if never reached, x_c = ∞).
- Initial data and far boundary: exact static seed in the new chart (`static_w.reference`), including beyond the vertex light cone; exact far data while V < 0.
- Deviation form (exact static seed to round-off) until T_switch = F_∞ − 1, then full fields. 4th-order differences, degree-5 Hermite ghost at the shell, RK4 with Δt = 0.5 min ΔZ.
- Remedy A: discrete Hamiltonian projection of the initial warp on Z ≥ −0.2. Remedy B: κ = 10, switched off for Z > −0.02, full for Z < −0.05. Time derivatives of the Neumann data include B_T, R_T and Yφ_TT.
- L = 16 (shell exact for T < 32), T_final = 14 (runs stop earlier at the 25-min wall limit, H/H₀ outside [−20, 20], or shell lapse outside [10⁻⁴, 10⁴]).

## 5. Calibration and controls (must pass before the tuned runs are trusted)

- **C1 (old-chart reproduction):** registered c, δ = 0.1, Y = 0, dc = 10⁻², x_c = ∞: φ_b and H/H₀ at H₀τ = 1, 2, 3, 4 within 10⁻³ of the archived pilot `pilot_Y0_dc1e-2`.
- **C2 (growth rates):** registered c, x_c = ∞, deviation mode: δ = 0.1 rate within 10⁻³ of Chat 13/14's 1.6269–1.6270; δ = 10⁻³ rate within 2×10⁻⁴ of 1.65719 (skipped if the δ = 10⁻³ run does not fit in the wall-clock budget; then reported as not done).
- **C3 (chart independence):** the same registered-c run with x_c = ∞ and two finite x_c: agreement of φ_b, H/H₀, 𝒲 versus proper time to ≤ 10⁻³ before the old freeze; the two finite-x_c runs agree to ≤ 10⁻² after it.
- **C4 (known end state, informative):** registered c, Y = 0, δ = 0.1 continued in the new chart; report whether H/H₀ approaches the +1-branch value 0.63731 (M2: H₊²/H₀² = 0.406168) and φ_b → 0.991435. If it clearly goes elsewhere without explanation, tuned-run verdicts are downgraded to INCONCLUSIVE.
- **C5:** Weyl transport identity with matter (V2 form) on each tuned run: median relative residual ≤ 10⁻³; the wrong-factor controls (4H → 3H; R dropped from n·A) must be ≥ 10× worse.
- **C6:** radiation ledger median relative residual ≤ 10⁻⁴.
- **C7 (Remedies):** one tuned run repeated with κ = 0 and one without projection; the classification must not change.
- **C8 (perturbed parameter):** Y = 1 with d = 0.95 d\* and 1.05 d\* (sensitivity of r and of the recollapse to the tuning).
- **C9:** the tuned model with Y = 0 (no radiation) as a baseline.

---

### Dated note, 28 September 2026 (after the first calibration runs; still before any tuned-model run)

**Change of the chart map (method only; model, thresholds and decision rule unchanged).** The registered-c calibration run in the first chart map, F(x) = ln(1+e^{−x_c}) − ln(e^{−x} + e^{−x_c}), continued through the old freeze (H₀τ ≈ 4.56) and settled on the +1 branch (φ_b = 0.991435, H/H₀ = 0.63739 against the static values 0.991435 and 0.637313). It then showed a second lapse decay ∝ e^{−T} with the proper time saturating at H₀τ ≈ 6.59. This is a coordinate singularity of that map: as a function of the regular label U_K it is singular at U_K = e^{−x_c−c₀} > 0, beyond the vertex light cone. The map is replaced by F(x) = C − arcsinh(e^{−x}/(2e^{−x_c})), C = arcsinh(e^{x_c}/2), which has the same behaviour at early times (F ≈ x) and at the light cone (F′ ∝ e^{−x}), but is an entire function of U_K, so it has no singularity beyond the light cone. The x_c rule of §4 is kept (with this map the late shell lapse tends to about twice the value it had with the first map). Runs made with the first map are kept in `runs/cal/` (tags without `_asinh`) and reported as such.

### Second dated note, 28 September 2026 (before the main tuned-model batch)

**Method changes (model, thresholds, pass/fail rule and reliability criteria unchanged):**

1. **Chart map.** The asinh map of the first note is regular across the vertex light cone but, being unbounded in U_K, drives the far part of the grid into the crunch of the light-cone interior (the static seed's Milne region recollapses; its crunch is at |w| ≈ 6×10³ w_b). The default map is now `bounded` (`static_w.Chart`): s = e^{−u} = h(σ), σ = 2e^{−x_c} sinh(C − U), h(σ) = s_max σ/(s_max + ε ln(1+e^{−σ/ε})), s_max = 20, ε = e^{−x_c}. It equals the asinh map up to and across the light cone and saturates smoothly (U_K ≤ 20e^{−c₀}), keeping every static grid point at |w| ≤ 20 w_b. A side effect: when the shell's arrival label U_K tends to a finite limit (as for a de Sitter end state, where the brane's infinite future lies at finite U_K), the shell lapse in this chart diverges at a finite T. The run then stops when |B_b| > 40 or when the constraint criterion of §3 fails; the reliable part is what counts. The first map (`softplus`), whose coordinate singularity lies at U_K = e^{−x_c−c₀}, is kept as a cross-check chart (C3).
2. **Remedy A window.** Projection on Z ≥ −0.2 moved the initial Hamiltonian defect to the window edge and made it larger (3.3×10⁻³ against 1.1×10⁻³ at the shell corner; registered-c seed). The window is now Z ≥ −(L − 0.1), as in lab.py (−0.85L there), which puts the edge defect next to the far boundary, causally disconnected from the shell for T < 2L − 0.1.
3. **Disclosure.** To choose between the two chart maps, two coarse tuned-model development runs were made (d\*, Y = 1, dc = 10⁻², dz_fine = 10⁻³, `dev/x_dstar_Y1_bnd`, `dev/x_dstar_Y1_sp`). Their outcome was therefore seen before the main batch: both charts agree where they overlap (e.g. r = 0.0888 and 0.0889 at H₀τ = 5.0), and the bounded chart continues to a Weyl/radiation plateau with r ≈ 0.073 and a small negative residual vacuum (≈ −3.4×10⁻⁴ H₀²). These development runs are not used for the verdict; the main batch (§3 grid of runs, both spacings, both seeds, Y ∈ {0.3, 1, 3}) and the calibration runs are rerun with the final code.
