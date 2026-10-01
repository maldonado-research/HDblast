# A1 verdict: tuned-vacuum roll-off (28 September 2026)

Ricardo Maldonado's HDBLAST program; prepared with AI assistance. The rule applied here is the one fixed in `REGISTRATION.md` before the main batch. The numbers come from `A1_RESULTS.json`, produced by `a1_analyze.py` from the run summaries in `runs/`. Status: **numerical**, internal, not peer reviewed.

## Registered verdict at δ = 0.1: PASS (conservative), within the friction closure

| Y | Seeds dc = 10⁻², 10⁻⁴ (two spacings each) | r = 𝒲/rad at plateau | Ω_r | Reliable class |
|---|---|---|---|---|
| 0.3 | both | 0.463 | 0.785 | FAIL-Weyl |
| **1** | **both** | **0.073–0.074** | **1.12** | **PASS-conservative** |
| 3 | both | none (no plateau before the runs ended, H₀τ ≈ 31–45) | none | UNRELIABLE / INCONCLUSIVE |
| 0 (baseline C9) | both | none (no plateau, no radiation by construction) | 0 | INCONCLUSIVE |

- At Y = 1, both seeds agree at both spacings (r within 0.0005). The near-shell constraint residuals are ≤ 8×10⁻⁵.
- r = 0.073 corresponds to ΔN_eff ≈ 0.21. It passes the conservative threshold (r ≤ 0.1). It **fails the combined CMB+BBN+BAO threshold** (r ≤ 0.03) and the ACT DR6 snippet level (r ≈ 0.025).
- PASS-combined: no run.

## Controls

- **C1–C3 (calibration, chart independence):** pass.
- **C2 at δ = 10⁻³:** **fail**. The fitted rate is 1.698 against the target 1.65719. The sliding rate is still drifting (1.659 in the last window), so the run was too short to settle.
- **C4 (known end state):** the bounded-chart run ends at H/H₀ = 0.6322 and φ_b = 0.99268, against the static branch values 0.6373 and 0.99143. That is within 0.5%, and still settling when the run stopped.
- **C5 (Weyl identity, Y = 1):**
  - The median residual is 3.8×10⁻⁵, which passes.
  - The 4H → 3H wrong-factor control is 360× worse, which passes.
  - **The "R dropped" control is only about 3× worse, which fails the registered 10× requirement.** This control is too weak to discriminate at this level, so it cannot confirm the identity; it is not evidence that the identity fails.
- **C6 (radiation ledger):** median 7.7×10⁻⁵, which passes (at Y = 0.3 it is 1.26×10⁻⁴, just above 10⁻⁴).
- **C7 (remedies):** with κ = 0, and without the projection, the class is still PASS, with r = 0.0738 in both cases. The finest spacing (dz = 2.5×10⁻⁴) also gives PASS, with r = 0.0733. The softplus cross-check chart fails at its known coordinate singularity before any plateau.
- **C8 (tuning sensitivity, Y = 1):**

  | d / d\* | Outcome |
  |---|---|
  | 0.95 | FAIL: vacuum dominated (Ω_r = 0.05) |
  | 0.99 | FAIL: vacuum dominated (Ω_r = 0.32) |
  | **1.00** | **PASS** |
  | 1.01 | H crosses zero at H₀τ = 13.5 and recollapses at 18.9; classed INCONCLUSIVE |
  | 1.05 | FAIL: recollapse at H₀τ = 8.2 |

## δ = 10⁻³ (the registered-scale detuning): not completed

- The pre-run never reached the light-cone criterion, so x_c was left at ∞ (the old chart). Both δ = 10⁻³ runs stopped when the shell lapse diverged:
  - dz = 2×10⁻⁴ at H₀τ = 11.0;
  - dz = 1.5×10⁻⁴ at H₀τ = 3.1.
- The two spacings do not agree on how long they last, and neither gives a classification. The growth-rate calibration there also failed (C2).
- Per the registration, the δ = 0.1 verdict therefore stands on its own and says nothing yet about δ = 10⁻³.

## What this does and does not show

- **Shown (numerical, within the model and closure):** with the quadratic tension tuned to d\* and the friction closure at Y = 1, the roll-off produces a radiation-dominated phase with a dark-radiation fraction below the conservative N_eff bound. This is the first HDBLAST run that meets a registered radiation-era criterion.
- **Not shown:**
  - **Fine tuning.** A 1% change in d removes the radiation era: the residual vacuum energy either dominates or turns the expansion around. The PASS requires the residual vacuum energy to be tuned to nearly zero, which is the cosmological-constant problem restated, not solved.
  - **A transient era.** Even at d\*, the residual vacuum is slightly negative (Ω_vac ≈ −0.21 at the plateau and growing as the radiation dilutes). The radiation phase therefore ends in recollapse unless something else changes.
  - **A stand-in source.** The energy transfer is the phenomenological friction closure (κ₅²j = Yv), not the derived χ-field source. Only one Y value of the three passes.
  - **Only δ = 0.1.** The registered detuning δ = 10⁻³ is not completed.
  - **The strict bound.** r = 0.073 fails the stricter combined N_eff bound.
- **Disclosure:** two development runs at d\* and Y = 1 were seen before the main batch (second dated note in `REGISTRATION.md`).

## Next steps

1. Replace the friction closure with the derived χ source.
2. Finish δ = 10⁻³: fix the x_c rule for runs that never reach the light-cone criterion, and run a longer growth calibration.
3. Extend the Y = 3 runs to a plateau.
4. Scan Y between 0.3 and 3 to see whether any Y reaches r ≤ 0.03.
