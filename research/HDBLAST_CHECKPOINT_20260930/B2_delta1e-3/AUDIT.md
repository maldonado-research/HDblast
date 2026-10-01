# AUDIT of workstream B2_delta1e-3 (independent, skeptical review)

Audit date: 1 October 2026. Auditor wrote only this file and `audit/` (scripts, JSON outputs, one rerun). Nothing else in the
folder was modified. One core used.

## Overall

The registered verdict **INCONCLUSIVE at δ = 10⁻³ is supported**. I reproduced one main run bit for bit and recomputed the key
numbers from the saved time series. An estimator I chose myself confirms the C2 seed-offset explanation independently. Three
problems remain:

- the headline overgeneralises the breakdown mechanism to "every Y = 1 run";
- "x_c ruled out as a cause" is contradicted by the producer's own controls;
- the README reinterprets the pre-stated acceptance criteria of Amendment 2 after the fact.

## Checks performed (scripts and outputs in `audit/`)

| Check | Script / output | Result |
|---|---|---|
| Rerun of `main_Y1_dc1e-2_S1` (same command, no checkpoint, 369 s) | `audit/rerun/`, `audit/compare_rerun.py` → `audit/AUDIT_RERUN.json` | **Bit-identical** to the producer's run: T_end 11.5800, H₀τ_end 12.473571386862394; max diff 0 in T, H₀τ, φ_b, H, lapse, R, 𝒲, rad over all 1164 records |
| C2 registered estimator recomputed independently (np.polyfit) | `audit/audit_c2_derivative.py` → `audit/AUDIT_C2.json` | Matches: C2a/C2a8 rate 1.656301, spread 2.44×10⁻³; C2b 1.656318, spread 2.39×10⁻³ → registered C2 FAIL confirmed |
| **Independent offset test**: log-slope of d(dev)/dτ, which removes a constant offset B exactly with no fit, on the *registered* windows 2–5 | same | C2a8 **1.657184** (spread 2.6×10⁻⁵); C2b **1.657187** (spread 1.6×10⁻⁵); C2c 1.657195 (4.5×10⁻⁵); C2d 1.657048 (6.0×10⁻⁵). Target 1.65719. The seed-offset explanation is confirmed by a method the producer did not use. (This audit estimator is post hoc; it does not change the registered FAIL.) |
| Run end times, φ_b, Ω_r, proper time at F_inf, recomputed from the time series | `audit/audit_runs.py` → `audit/AUDIT_RUNS.json` | All quoted end times, reliable ends and φ_b maxima reproduce (see corrections for two definitional mismatches) |
| D6 old-chart freeze | same | dB_b/dT = −0.65, −0.83, −0.93, −0.97, −0.984 at T = 12, 13, 14, 15, 15.5; b_b + T → 11.364; H₀τ → 11.35; φ_b 0.094; H/H₀ 1.002 → confirmed |
| d\*(10⁻³) and Yρ_b | same | d\* = −3.1942416958680835 matches M8 JSON; 1 + c + d\*/2 = 4.7×10⁻⁴; Yρ_b = 78.8 |
| Location and convergence of the breakdown | inline read of time series (numbers below) | At T ≈ 11.43–11.5 the near-shell momentum residual is 0.38 (S1), 0.039 (S2) and ~0.001–0.01 (S3′), at z ≈ −0.1 to −0.2. At fixed T it converges by about 10× per doubling. The breakdown is a steepening feature that each resolution resolves for a little longer, not an O(1) continuum blow-up at a fixed time |
| Chart variants | `static_w.Chart` source | `bounded` equals `asinh` exactly (h(σ) = σ) until the light cone, and s_max acts only beyond it. E1 and E3 are therefore the same near-shell gauge as the default, which explains their 10-digit identity. They test far-region saturation only, not the near-shell gauge |
| a1_analyze.py unchanged | sha256 | Identical to the A1 file (5c3000ef…) |
| Registration timing | `stat` birth times | See below |

## Registration

- **Original text.** `REGISTRATION.md` birth time is 2026-09-30 07:32:05. That is before `jobs_stage1.txt` (07:32:29) and before
  the first Y = 1 pre-run output (07:32:57). It is after the dev runs it discloses (D3 at 07:30:27). This is consistent with
  "registered before any Y = 1 run". Whether the original body is unedited cannot be checked: there is no snapshot, and the
  mtime (1 Oct 09:50:02) reflects the appended notes. The body's content is consistent with what was disclosed.
- **Note 1 (Amendments 1–2)** came after the κ = 10 development pre-run had been seen. The README discloses this.
- **Note 2** ("before any main run"): XC_SELECTION was written at 08:41:10 and the main queue started at 08:41:34. This is
  plausible but cannot be verified to the minute.
- **Note 3** (09:20) follows the D6 continuation, which was launched at 09:14 without a prior note. D6 is a diagnostic, so this
  is minor.
- **Note 4** (dated 09:48) "adds" S3′, but the S3′ log was created at 09:47:11, so the run was launched about a minute before
  the note. This is minor: S3′ is diagnostic only.
- **Note 5** (09:50) precedes E4 (09:51:12): OK.
- **Amendment 2 post-hoc reinterpretation.** Amendment 2 accepts the seed-transient interpretation only if three criteria hold:
  (i) the late estimator passes; (ii) the window rates increase monotonically with geometrically shrinking gaps; (iii) C2d shows
  the same behaviour.
  - (ii) fails as written: `B2_RESULTS.json` reports `monotone: false`, because the late windows turn down.
  - For (iii), C2d *fails* the late estimator (diff −2.2×10⁻⁴).
  - The README says "(ii) holds up to the 5–6 window" and "(iii) holds: same gap ratios". That reinterprets the pre-stated
    criteria after the data were seen.
  - Practical effect: none, since no tuned run was classified. My derivative estimator independently supports the
    interpretation anyway. But by the amendment's own rule the interpretation is **not accepted**, and the README should say so.
- README says "four dated notes"; there are five.

## Claim-by-claim verdicts

1. **κ = 0 pre-run unstable; κ = 10 gives x_c = 11.3** — *confirmed (numerical).*
   - XC_SELECTION shows lapse-rise failures at κ = 0 and x_c = 11.2/11.3 at 32/42.7 points per wall.
   - The D6 freeze at b_b + T → 11.364 independently corroborates 11.3.
   - The dc = 10⁻⁴ value (κ = 10: 18.7) rests on one resolution.
2. **A1 spacings stopped for numerical reasons** — *partially confirmed.*
   - D2 supports both diagnoses: the projection diverged (7×10⁻⁷ → 1.6×10⁻² → 2×10⁻⁶, warp shift 9×10⁻³), and the residual
     passed 0.05 at T ≈ 7.0 near z ≈ −0.05.
   - However, the dz = 2×10⁻⁴ run ended at H₀τ 10.97, close to the physical light-cone crossing (11.36) that B2 itself found.
     Its *end* is therefore not cleanly separable from the crossing. Only its loss of reliability (T ≈ 7) is clearly numerical.
3. **Registered C2 fails; the cause is a neutral seed offset; the A1 drift was the offset plus the table kick** — *confirmed.*
   - Recomputed the registered rates exactly.
   - The offset-free derivative estimator gives 1.65718–1.65719 with spread ≤ 3×10⁻⁵ on the registered windows.
   - The gap ratios of 0.42–0.45 against e^{−λ/2} = 0.437 are in the JSON.
   - The 35× table-kick ratio is in D2/D1.
4. **The supplementary late estimator passes; the mode fit gives 1.65720 ± 3×10⁻⁵** — *partially confirmed.*
   - The numbers reproduce from D4 and B2_RESULTS.
   - The passes are marginal: spreads of 1.83×10⁻⁴ and 1.93×10⁻⁴ against a limit of 2×10⁻⁴, with late windows drifting downward.
   - The estimator and fit are post-registration, and Amendment 2's own acceptance criteria are not met as written (see above).
5. **Registered tuned test INCONCLUSIVE** — *confirmed.*
   - Rerun bit-identical.
   - None of the four main Y = 1 runs reaches a plateau or recollapse.
   - Max Ω_r ≤ 0.0032 within the reliable part.
   - C5/C6 medians pass as registered.
   - Caveat: the 4H → 3H control discriminates only at the median. At p95 the identity residual (e.g. 4.3×10⁻⁵, S1) is ≥ the 3H
     control (3.7×10⁻⁵), and "R dropped" is indistinguishable. C5 is a weak check in this regime.
6. **Breakdown mechanism (crossing before roll; x_c ruled out)** — *partially confirmed / partly refuted.*
   - Confirmed:
     - the D6 crossing at H₀τ ≈ 11.36 with φ_b ≈ 0.094 and H/H₀ ≈ 1.002;
     - lapse growth of 3.55/4.00/4.36 per unit T (S1/S2/S3′) against 1.24 at δ = 0.1;
     - E4 ends within 0.0014 in T of the default.
   - Refuted:
     - **"Ruled out: the x_c value."** The end tracks x_c. With dc = 10⁻², moving x_c from 11.3 to 12.3 moves the run end from
       H₀τ 12.47 to 13.22 and the reliable end from 11.44 to 12.70, a larger gain than any resolution doubling. With dc = 10⁻⁴,
       moving x_c from 11.8 to 18.7 moves the end from 12.9 to 19.8. What is robust is that every run dies within about 0.3–0.9
       chart units of F_inf.
     - **"Every Y = 1 run breaks down 1–2.5 H₀⁻¹ after the crossing"** is false for the main dc = 10⁻⁴ pair. Those runs end
       *before* F_inf (T_end − F_inf = −0.20 and −0.10), at H₀τ 12.9/13.6. That is about 5–6 H₀⁻¹ *before* the crossing their own
       κ = 10 pre-run and sensitivity set locate (H₀τ ≈ 18.7–18.8). They fail because the registered x_c = 11.8 (a κ = 0
       artifact) puts the chart's F_inf early. The claim is also false for C7 (κ = 0, ends at T = 6.79).
   - E1/E3 do not test the near-shell gauge, because the charts coincide there.
7. **Resolution scaling; plain refinement cannot finish** — *partially confirmed.*
   - Confirmed: end H₀τ 12.47/13.10/13.77 and reliable end 11.44/12.12/12.55.
   - "Cannot fix", "~10 more doublings" and "10⁶× cost" are a linear extrapolation from three points of one seed, at one x_c, on
     one chart family. They should be labelled conditional/extrapolated, not numerical.
   - The fixed-T residual converges at roughly 10× per doubling. This supports under-resolution of a steepening feature but
     does not prove that refinement is hopeless.
8. **Y = 0 baseline works for a fast roll** — *partially confirmed.*
   - Confirmed: φ_b → 1.000, and the runs go 1.27–1.56 chart units past F_inf.
   - H/H₀ is 0.14 at F_inf and 0.057 at the end; the quoted 0.08–0.13 is approximate.
   - However, both S1 Y = 0 runs lose reliability (residual > 0.05) at H₀τ 2.46 and 3.41, i.e. *at or before* the crossing (2.59,
     3.69). Only S2 stays reliable to H₀τ 7.5–8.2. "The solver works for a fast roll" holds at S2 only.
9. **C7 κ = 0 unusable; class unchanged** — *confirmed.*
   - The near-shell H residual reaches 0.33 and M ≈ 1.0 at the stop, with the lapse at 9. The summary's "residual 1e-2"
     understates this.

## Corrections requested

- Restrict the headline mechanism to the dc = 10⁻² runs, the sensitivity set and the S3′/C3′/E runs. State that the registered
  dc = 10⁻⁴ pair ends before F_inf, about 5–6 H₀⁻¹ before its physical crossing, because of the registered x_c = 11.8 artifact.
- Remove "x_c value" from the "ruled out" list. x_c shifts the reach past the crossing by as much as or more than a resolution
  doubling; only the qualitative failure is x_c-independent.
- State that, by Amendment 2's own criteria (ii) and (iii) as written, the seed-transient interpretation is not formally
  accepted; the agreement rests on post-hoc estimators. Optionally cite the offset-free derivative estimator (audit) as
  independent support.
- Relabel "plain grid refinement cannot fix this / ~10⁶ cost" as conditional (an extrapolation).
- README: say five dated notes, not four; note that S3′ was launched about a minute before note 4.
- Minor quantities:
  - δ = 0.1 crossing H/H₀: D7 gives 0.20 at F_inf, not 0.26.
  - B2_RESULTS `H0tau_end_run` (last record, e.g. 12.37) differs from the quoted summary H₀τ_end (12.47); state which is meant.
  - Max Ω_r ≤ 0.003 holds only within the reliable part; over all records it reaches 0.013 (C3′).
- Y = 0 baseline: note that the S1 runs fail the residual criterion near the crossing.
