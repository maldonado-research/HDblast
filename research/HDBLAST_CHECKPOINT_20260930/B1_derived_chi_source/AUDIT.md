# Independent audit of workstream B1 (derived χ source in place of the friction closure)

Audit date: 1 October 2026. Auditor scope: read-only on the producer's files. The only files written are this AUDIT.md and the `audit/` subfolder. Cores used: 1.

## Overall finding

**The registered verdict, INCONCLUSIVE at λ_c = 1, is correct and reproducible.** The headline descriptive numbers at b_cut reproduce exactly:
- r = 94, 674, 3.56×10³ and 8.47×10³ for y = 1, with Ω_r ≤ 0.012;
- spacing agreement better than 0.1%.

I found three things to correct, none of which changes the registered verdict:

1. **Production-model controls are reported against the wrong baseline.** README §6.5 and the producer summary say the controls "change the end r from 0.17 to −338, −1149 and 22.5", and they label −338, −1149 and 22.5 as "r at run end".
   - Those three values are actually r at the *reliable end* (H₀τ = 1.85, 1.57 and 2.62). The baseline 0.17 is r at the *run end* (H₀τ = 5.19).
   - Compared at the run end on both sides, the controls give r = 0.316 (ramp 1.5), 0.082 (instant) and 0.169 (no D), against a baseline of 0.170. That is +86%, −52% and +0.3%.
   - So production-model sensitivity at b_need is real for the ramp and the instant insertion, but it is about 2× in magnitude, not orders of magnitude. The O(1/q) D correction has no effect.
2. **The cell φ\* = 0.9, G = 100, y = 1 at b_cut does not meet the "χ fully decayed, 𝒲a⁴ and Ra⁴ flat to ≤ 1%" description.**
   - At the end, Ω_χ/Ω_r ≈ 0.095, and Ra⁴ rises 16% over the last Δln a ≈ 0.12.
   - Its r = 3.56×10³ is also inflated by about 1.9× by the field-ramp underproduction (the producer discloses this).
   - The conclusion r ≫ 0.1 is unaffected.
3. **The y = 0.1 b_cut values (r = 306 to 1.1×10⁵) are snapshots, not lower bounds.**
   - At the end, undecayed χ carries 0.6 to 44 times the radiation energy.
   - Under a linear "all χ → R" estimate, r would fall to about 52 (φ\*=0.5, G=100), 194, 2.5×10³ and 4.4×10³. Since massive χ redshifts more slowly than radiation, the true values could be lower still.
   - These are still ≫ 0.1. The ≥ 94 statement holds only for y = 1, which is how the producer phrases it.

Two smaller points:
- The "constraints ≤ 1.2×10⁻⁴" figure holds only on the fine grid. On the coarse grid the maximum is 4.2×10⁻⁴.
- The b_cut cells appear as "UNRELIABLE" in `B1_RESULTS.json` → `aggregate.b_cut_classes`. That is because no plateau r exists to compare. Their constraints are fine; "INCONCLUSIVE at both spacings" is the more accurate label.

## Registration integrity

- **Birth time.** `REGISTRATION.md` was created at 2026-09-30 07:17:11 UTC (stat Birth). This is before:
  - `B1_REDUCED.json` (07:19:14);
  - the first calibration runs (07:21);
  - the superseded time-ramp grid cell (07:24–07:26);
  - every final-build grid run (1 Oct, 08:36–09:50).
- **The body text reads as unedited.**
  - It still describes the time ramp (§1.3).
  - It mentions "20 exact identities", while the current symbolic file has 22 (C4c and C4d were added for the field ramp, as the second note says).
  - It discloses the pre-registration dev runs `smoke_chi` and `dev2`, whose files predate the birth time (06:59 and 07:13).
- **What I could not verify (git not used).**
  - The current modification time (1 Oct 09:29:15) equals the time of `jobs_explore.txt`, which fits the third note being appended then.
  - I cannot verify from timestamps that the second note was appended before the grid batch started at 08:36. Its content is consistent with that, but this is unverifiable.
- **Post-registration changes, all disclosed with reasons:**
  - y = 0.1 cells moved to the y = 1 b_need, made after seeing the reduced estimate;
  - time ramp replaced by the field-space ramp, made after seeing the time-ramp cell stop non-finite and the dev3/dev4 outcomes.
  - Neither change could flip the aggregate, because every b_cut cell is INCONCLUSIVE either way.
- **Superseded and dev runs.** Their code hashes match the disclosure: 888cfc92… for the superseded runs, a0ae79e1… for the final build and dev4.

## Checks performed (scripts and JSON in `audit/`)

| Check | Script → output | Result |
|---|---|---|
| Independent sympy re-derivation: j = n∂ω/∂φ; Ward identity with time-dilated decay; Q cancelling in the total; shell budget with n·φ = −(σ′+κ²J)/2; preheating integrals; rad decomposition. Controls: no time dilation; + sign | `audit_algebra.py` → `AUDIT_ALGEBRA.json` | all 0; both controls nonzero |
| Rerun of a copy of `b1_symbolic.py` | `b1_symbolic_copy.py` → `AUDIT_SYMBOLIC_RERUN.json` | identical: 22/22 identities and 8/8 controls |
| Rerun of a copy of `b1_reduced.py` | `b1_reduced_copy.py` → `AUDIT_REDUCED_RERUN.json` | 0 differences from `B1_REDUCED.json` (relative tolerance 1e-9) |
| K1 and K2 compared directly with the A1 npz archives, without the producer's analysis code | `audit_calibration.py` → `AUDIT_CALIBRATION.json` | K1: 0.0 difference on 305 records. K2: φ_b, H, 𝒲 within 1.5e-11; same end time |
| Recomputation of r, Ω_r, Ω_χ, 𝒲a⁴, Ra⁴, reliable end, R_j and constraints for every run; brute-force longest A1 plateau window (5%, H > 0) | `audit_timeseries.py` → `AUDIT_TIMESERIES.json` | see below |
| 5D rerun: exact rerun of `main_ps0.5_G100_y1_bcut_dzf1e-3`, a third spacing dz_f = 7.5e-4, and a wrong-sign (j → −j) run | `rerun.sh`, `audit_reruns.py` → `AUDIT_RERUNS.json` | bit-identical rerun. Third spacing: r = 94.044 (vs 94.054 and 94.041). Wrong sign: r = 94.19 (+0.14%), confirming the test-field regime at b_cut |

**Recomputed values, confirming the producer's tables:**
- b_cut y = 1: r = 94.05/94.04, 673.9/673.9, 3562.7/3562.2 and 8474/8473, with Ω_r = 0.0116, 0.0016, 3.1e-4 and 1.3e-4.
- Explore cell: r = 0.2820/0.2819 with Ω_r = 0.819 and R_j = 0.083.
- b_need run ends: 0.17, 0.91, −0.09, 0.058 (y = 1).

**Longest window satisfying the A1 plateau tolerance, over all runs:** at most Δln a = 0.394 (b_cut) and 0.251 (explore). Both are below 0.5, so no run reaches a plateau and INCONCLUSIVE is independently confirmed.

## Per-claim verdicts

| # | Claim | Verdict |
|---|---|---|
| 1 | Closed shell system exact-verified | confirmed |
| 2 | Calibration K1, K2 | confirmed |
| 3 | Registered aggregate INCONCLUSIVE | confirmed |
| 4 | r ≥ 94 (y = 1) at λ_c = 1, Ω_r ≤ 0.012, descriptive | partially confirmed (numbers exact; the "fully decayed / flat" description fails for φ\*=0.9, G=100; the constraint bound holds for the fine grid only; y = 0.1 values are not bounds) |
| 5 | Reduced estimate FAIL-Weyl / recollapse (conditional) | confirmed (exact rerun). The 5D agreement of 3–4% is in end values; during the decay R_5D/R_red is 0.74–0.92 because of the ramp versus instant insertion |
| 6 | λ_c needed 11–50 (17–74), unreliable | confirmed as conditional and unreliable (arithmetic checked, e.g. 53.61 × 0.01^{1/3} = 11.55) |
| 7 | Exploratory r = 0.282 at λ_c = 5.36 | confirmed (not registered; r ∝ b^−1.15 arithmetic checked) |
| 8 | "Negative": channel cannot replace the friction closure below M₅ | partially confirmed. Well supported for y = 1, but the label should read "negative (conditional/descriptive)": the registered outcome is INCONCLUSIVE |
| 9 | b_need unreliable; production sensitivity "large" | partially confirmed. Unreliability and the φ = −1 crossing at Z ≈ −0.4, T = 1.85 are confirmed. The control r values are mis-paired (see Overall finding, item 1) |
| 10 | K3 consistency (ledger, Weyl identity, controls) | confirmed as reported, including the disclosed exceptions. The wrong-sign control is non-discriminating at b_cut, as expected (it is registered for b_need only) |

## Reproduction of this audit

```bash
cd /home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source
export OMP_NUM_THREADS=1
python3 audit/audit_algebra.py
python3 audit/b1_symbolic_copy.py
python3 -W ignore audit/b1_reduced_copy.py
python3 -W ignore audit/audit_calibration.py
python3 -W ignore audit/audit_timeseries.py
audit/rerun.sh                              # about 6 min on 1 core
python3 -W ignore audit/audit_reruns.py
```
