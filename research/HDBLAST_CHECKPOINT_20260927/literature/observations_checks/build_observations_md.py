#!/usr/bin/env python3
"""Validate the observation-sweep records and render ../observations.md.

Inputs : sources_observations.py (records), knee_and_scale_checks.json (numbers)
Outputs: ../observations.md, observations_sources.json
Run    : python3 build_observations_md.py            (build)
         python3 build_observations_md.py --selftest (negative controls: the
                                                      validator must reject corrupted records)
"""
from __future__ import annotations

import copy
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import sources_observations as SO  # noqa: E402

IMPLICATIONS = {"supports", "constrains", "method", "context", "challenges"}
ACCESSES = {"search-snippet-only", "full-text-read"}
SECTIONS = {"PTA": "Pulsar timing arrays (PTA) and the nanohertz background",
            "CMB": "CMB: dark radiation, spectral tilt and primordial tensors",
            "BBN": "Big Bang nucleosynthesis and braneworld dark radiation",
            "DESI": "DESI and evolving dark energy",
            "LVK": "LIGO/Virgo/KAGRA, GW170817 and gravitational-wave propagation",
            "ISL": "Tabletop tests of Newton's inverse-square law",
            "H0": "The Hubble tension",
            "JWST": "JWST early galaxies"}


def validate(records, queries):
    errs = []
    seen = set()
    for r in records:
        rid = r.get("id", "?")
        for k in ["id", "title", "url", "year", "venue", "access", "finding", "relevance", "implication", "queries", "section"]:
            if not r.get(k) and r.get(k) != 0:
                errs.append(f"{rid}: missing {k}")
        if not str(r.get("url", "")).startswith("https://"):
            errs.append(f"{rid}: url not https")
        if r.get("url") in seen:
            errs.append(f"{rid}: duplicate url")
        seen.add(r.get("url"))
        if r.get("implication") not in IMPLICATIONS:
            errs.append(f"{rid}: bad implication {r.get('implication')}")
        if r.get("access") not in ACCESSES:
            errs.append(f"{rid}: bad access {r.get('access')}")
        if r.get("section") not in SECTIONS:
            errs.append(f"{rid}: bad section")
        y = str(r.get("year"))
        if not (re.fullmatch(r"(19|20)\d\d", y) or y in {"n/a", "unknown"}):
            errs.append(f"{rid}: bad year {y}")
        for q in r.get("queries", []):
            if not (1 <= q <= len(queries)):
                errs.append(f"{rid}: query index {q} out of range")
        m = re.search(r"arxiv\.org/(?:abs|html|pdf)/(\d{4})\.(\d{4,5})", r.get("url", ""))
        if m and re.fullmatch(r"20\d\d", y):
            yy, mm = int(m.group(1)[:2]), int(m.group(1)[2:])
            if not (1 <= mm <= 12):
                errs.append(f"{rid}: arXiv id month invalid")
            # journal year may be later than the arXiv year, never earlier
            if int(y) < 2000 + yy:
                errs.append(f"{rid}: year {y} earlier than arXiv id year 20{yy:02d}")
    if len(queries) < 15:
        errs.append("fewer than 15 queries")
    return errs


def number_checks(J):
    """Numbers quoted in the prose of the records/markdown must match the JSON."""
    A, B, C = J["A_knee_benchmark"], J["B_braneworld_scale_setting"], J["C_dark_radiation"]
    pp = A["amplitude_at_1_per_yr"]["per_PTA"]
    ts = B["tabletop_scale_setting"]
    ew, lr = ts["EotWash2020_Yukawa_range_38.6um_proxy"], ts["LivingReview_0.1mm"]
    de = B["if_brane_dS_is_todays_dark_energy"]
    pairs = [
        ("hc at 1/yr 2.40e-15", A["amplitude_at_1_per_yr"]["hc_SBPL_at_f_yr_h0674"], 2.40e-15, 0.01),
        ("hc at 1/yr h=0.735 2.62e-15", A["amplitude_at_1_per_yr"]["hc_SBPL_at_f_yr_h0735"], 2.62e-15, 0.01),
        ("ratio vs NG15 CNM 1.31", pp["NG15_CNM_A_fixed"]["SBPL_over_PTA_Omega_at_f_yr"], 1.31, 0.01),
        ("strain ratio vs NG15 CNM 1.14", pp["NG15_CNM_A_fixed"]["SBPL_hc_over_A"], 1.14, 0.01),
        ("ratio vs PPTA 1.44", pp["PPTA_DR3_A_13_3"]["SBPL_over_PTA_Omega_at_f_yr"], 1.44, 0.01),
        ("ratio vs EPTA DR2new 0.92", pp["EPTA_DR2new_A_13_3"]["SBPL_over_PTA_Omega_at_f_yr"], 0.92, 0.01),
        ("ratio vs MPTA 13/3 0.25", pp["MPTA_4p5yr_A_13_3"]["SBPL_over_PTA_Omega_at_f_yr"], 0.25, 0.02),
        ("ratio vs MPTA free 0.10", pp["MPTA_4p5yr_A_free"]["SBPL_over_PTA_Omega_at_f_yr"], 0.10, 0.04),
        ("f_k/f_yr 0.9972", A["knee_location"]["f_k_over_f_yr"], 0.9972, 1e-4),
        ("NG15 bins from 1/yr 0.045", abs(A["knee_location"]["per_dataset"]["NANOGrav_15yr"]["knee_minus_1_per_yr_in_bins"]), 0.045, 0.02),
        ("MPTA bins from 1/yr 0.013", abs(A["knee_location"]["per_dataset"]["MPTA_4p5yr"]["knee_minus_1_per_yr_in_bins"]), 0.013, 0.05),
        ("EPTA bins from 1/yr 0.03", abs(A["knee_location"]["per_dataset"]["EPTA_DR2new"]["knee_minus_1_per_yr_in_bins"]), 0.03, 0.1),
        ("MPTA knee bin 4.49", A["knee_location"]["per_dataset"]["MPTA_4p5yr"]["knee_in_units_of_1_over_T"], 4.49, 0.01),
        ("CPTA knee bin 2.99", A["knee_location"]["per_dataset"]["CPTA_DR1"]["knee_in_units_of_1_over_T"], 2.99, 0.01),
        ("MPTA resolution 7.04 nHz", A["knee_location"]["per_dataset"]["MPTA_4p5yr"]["frequency_resolution_1_over_T_nHz"], 7.04, 0.01),
        ("EPTA resolution 3.08 nHz", A["knee_location"]["per_dataset"]["EPTA_DR2new"]["frequency_resolution_1_over_T_nHz"], 3.08, 0.01),
        ("Omega 25 Hz 2.28e-35", A["LVK_25Hz"]["Omega_SBPL_25Hz"], 2.28e-35, 0.01),
        ("26 orders below LVK", A["LVK_25Hz"]["orders_of_magnitude_below_limit"], 26.0, 0.01),
        ("DeltaNeff ACT upper 0.071", C["ACT_DR6_gaussian_95_upper_DeltaNeff"], 0.071, 0.01),
        ("rho_dr/rho_gamma ACT 0.016", C["ACT_DR6_rho_dr_over_rho_gamma_upper"], 0.016, 0.02),
        ("DeltaNeff BBN upper 0.46", C["BBN_2024_gaussian_95_upper_DeltaNeff"], 0.46, 0.01),
        ("rho_dr/rho_gamma BBN 0.10", C["BBN_2024_rho_dr_over_rho_gamma_upper"], 0.104, 0.01),
        ("H ell 0.069", B["registered_H_ell"]["H_ell_from_rho_b"], 0.0693, 0.01),
        ("rho/lambda 0.0024", B["registered_H_ell"]["rho_over_lambda_estimate=(H ell)^2/2"], 0.0024, 0.01),
        ("ell for dark energy 372 Mpc", de["ell_needed_Mpc"], 372.0, 0.01),
        ("ratio to 38.6 um 3e29", de["ratio_to_38.6um"], 3.0e29, 0.02),
        ("delta for dark energy ~1e-62", de["delta_needed_for_H_dS_at_38.6um_leading_order"], 1.13e-62, 0.01),
        ("M5 at 38.6um 3.1e8 GeV", ew["M5_GeV"], 3.1e8, 0.02),
        ("lambda^1/4 at 38.6um 5.5 TeV", ew["lambda_quarter_GeV"], 5.5e3, 0.01),
        ("rho^1/4 at 38.6um 1.2 TeV", ew["rho_vac_quarter_GeV"], 1.22e3, 0.01),
        ("M5 at 0.1mm 2.3e8 GeV", lr["M5_GeV"], 2.27e8, 0.01),
        ("f0 at 38.6um 8.4e-5", ew["f_today_horizon_scale_Hz(f*/H*=1)"], 8.36e-5, 0.01),
        ("f0 at 0.1mm 5.2e-5", lr["f_today_horizon_scale_Hz(f*/H*=1)"], 5.19e-5, 0.01),
        ("T for knee g*=106.75 0.19 GeV", B["temperature_for_knee"]["106.75"], 0.19, 0.01),
        ("T for knee g*=10.75 0.28 GeV", B["temperature_for_knee"]["10.75"], 0.278, 0.01),
        ("bin1 ratio 0.035", A["NG15_bins"]["rows"][0]["ratio_SBPL_over_NG15_A_13_3"], 0.035, 0.05),
        ("bin3 ratio 0.15", A["NG15_bins"]["rows"][2]["ratio_SBPL_over_NG15_A_13_3"], 0.153, 0.01),
        ("20-yr resolution 1.58 nHz", A["hypothetical_20yr_span"]["resolution_nHz"], 1.58, 0.01),
        ("20-yr 1/yr bin ~20", A["hypothetical_20yr_span"]["one_per_yr_bin"], 20.0, 0.01),
        ("NG15 hc reproduces 2.4e-15 within 0.1%", A["amplitude_at_1_per_yr"]["hc_SBPL_at_f_yr_h0674"], 2.4e-15, 0.001),
        ("DeltaNeff exact 7.09e-4", A["delta_neff"]["DeltaNeff_exact"], 7.09e-4, 0.001),
        ("T for 25 Hz f*/H*=1 ~1.5e8 GeV", B["temperature_for_band_GeV"]["LVK_25Hz"]["fH1"], 1.5e8, 0.01),
        ("T for 25 Hz f*/H*=100 ~1.5e6 GeV", B["temperature_for_band_GeV"]["LVK_25Hz"]["fH100"], 1.5e6, 0.01),
    ]
    bad = []
    for name, got, want, rtol in pairs:
        if abs(got / want - 1) > rtol:
            bad.append(f"{name}: got {got}, quoted {want}")
    return bad, len(pairs)


def fmt(x, n=3):
    return f"{x:.{n}g}"


SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def sci(x, n=2):
    """2.6e+03 -> 2.6×10³ (n significant digits)."""
    m, e = f"{x:.{n-1}e}".split("e")
    return f"{m}×10{str(int(e)).translate(SUP)}"


PRETTY = {"EotWash2020_Yukawa_range_38.6um_proxy": "ℓ ≤ 38.6 µm (Eöt-Wash 2020 Yukawa range, proxy)",
          "LivingReview_0.1mm": "ℓ ≤ 0.1 mm (review statement)"}


def render(records, J):
    A, B, C = J["A_knee_benchmark"], J["B_braneworld_scale_setting"], J["C_dark_radiation"]
    kl = A["knee_location"]["per_dataset"]
    ng = A["NG15_bins"]["rows"]
    amp = A["amplitude_at_1_per_yr"]["per_PTA"]
    ts = B["tabletop_scale_setting"]
    ew, lr = ts["EotWash2020_Yukawa_range_38.6um_proxy"], ts["LivingReview_0.1mm"]
    de = B["if_brane_dS_is_todays_dark_energy"]
    bk = A["project_BK_family_windows"]
    L = []
    w = L.append
    w("# Literature sweep 2: observations that could test or constrain HDBLAST (2023–2026)")
    w("")
    w("**Work in progress, 28 September 2026.** Nothing here is a claim of the programme until it appears in a dated checkpoint or a Zenodo record.")
    w("")
    w("**How this was done.** There were 48 distinct web searches, listed at the end; three more were refused when the shared search budget ran out. "
      "**No paper could be opened.** arxiv.org, zenodo.org and journal sites are blocked for downloads here. "
      "Every entry below therefore rests on **search-result snippets only**: titles, URLs and short summaries. Numbers quoted from papers are transcribed from those snippets and must be checked against the papers before reuse. "
      "Where a snippet summary did not make clear which page a statement came from, the entry says so. "
      "Every number *computed* here comes from `observations_checks/knee_and_scale_checks.py`, which writes `knee_and_scale_checks.json` and passes 12 built-in checks and controls. "
      f"The builder cross-checks {J['_n_number_checks']} numbers quoted in this text against that JSON.")
    w("")
    w("Labels: **exact-verified**, **numerical**, **conditional**, **negative**, **inconclusive**; for sources, *supports / constrains / challenges / method / context*.")
    w("")
    # --------------------------------------------------------------------------------------------
    w("## 1. Bottom line")
    w("")
    w("1. **The registered PTA knee is still untested at likelihood level. No 2025–2026 release performs the test it needs** (**inconclusive**). "
      "The new PTA results that bear on it do so only indirectly: amplitudes, index trends, running and piecewise reconstructions. "
      "The NANOGrav 20-yr and IPTA DR3 background analyses, which could test it directly, had not appeared in the searches.")
    w(f"2. **The knee sits exactly on a known blind spot** (**numerical** position; the blind spot comes from a snippet). "
      f"f_k = 3.16×10⁻⁸ Hz is {fmt(A['knee_location']['f_k_over_f_yr'],4)} × (1/yr), within {fmt(abs(kl['NANOGrav_15yr']['knee_minus_1_per_yr_in_bins']),2)} frequency bins of 1/yr for NANOGrav 15-yr "
      f"and {fmt(abs(kl['MPTA_4p5yr']['knee_minus_1_per_yr_in_bins']),2)} bins for MeerKAT. "
      "Fitting each pulsar's sky position and proper motion removes sensitivity around 1/yr. "
      "A bend there is therefore maximally degenerate with the timing model. Unless the analysis propagates that absorption, the knee's location is not identifiable. "
      "This is the most important methodological finding of this sweep. "
      f"A scan of {J['_scan']['files_scanned']:,} project text files found **no** discussion of this degeneracy "
      "(`observations_checks/project_scan_1yr.json`; the positive control passed). It is standard in the PTA literature but new to this project.")
    w(f"3. **Inside the NANOGrav band the frozen curve is effectively a γ = 3 power law.** Here γ is the timing-residual index, with Ω ∝ f^(5−γ). "
      f"Its unweighted effective index over the 14 standard bins is γ ≈ {fmt(A['NG15_bins']['effective_gamma_all14'],3)} (**numerical**). "
      "Curvature appears only in bins 13–14. "
      "The 2025–2026 noise-model reanalyses (EPTA; NANOGrav customized chromatic noise) move the inferred background toward γ = 13/3 and lower amplitude. "
      "That trend is **unfavourable** to the frozen low-frequency slope, but it is not a rejection. MeerKAT's free index, γ = 3.60 (+1.31/−0.89), still allows γ = 3.")
    w(f"4. **Amplitude.** The frozen height reproduces the 2023 NANOGrav amplitude at 1/yr (h_c = {sci(A['amplitude_at_1_per_yr']['hc_SBPL_at_f_yr_h0674'],3)}; **numerical**, a calibration). "
      f"Relative to the 2025–2026 values it has {fmt(amp['NG15_CNM_A_fixed']['SBPL_over_PTA_Omega_at_f_yr'],3)}× the power of NANOGrav's customized-noise amplitude "
      f"and {fmt(amp['PPTA_DR3_A_13_3']['SBPL_over_PTA_Omega_at_f_yr'],3)}× that of PPTA. "
      f"Against MeerKAT it has only {fmt(amp['MPTA_4p5yr_A_13_3']['SBPL_over_PTA_Omega_at_f_yr'],2)}× (γ = 13/3) to {fmt(amp['MPTA_4p5yr_A_free']['SBPL_over_PTA_Omega_at_f_yr'],2)}× (free index). "
      "The spread between teams already exceeds a factor of 4 in Ω, so one frozen height cannot match every team.")
    w(f"5. **Scale-setting sharply challenges the idea that the knee is the blast's signature** (**conditional**, order of magnitude). "
      f"In the registered model H·ℓ = {fmt(B['registered_H_ell']['H_ell_from_rho_b'],4)} (dimensionless; ℓ = 9 is **exact-verified**). "
      f"Torsion-balance tests bound ℓ ≲ 38.6–100 µm. That puts the static branch's brane vacuum energy at ρ^(1/4) ≳ {fmt(lr['rho_vac_quarter_GeV']/1e3,2)}–{fmt(ew['rho_vac_quarter_GeV']/1e3,2)} TeV. "
      f"Its horizon-scale GW frequency today would be f ≳ {sci(lr['f_today_horizon_scale_Hz(f*/H*=1)'])}–{sci(ew['f_today_horizon_scale_Hz(f*/H*=1)'])} Hz. "
      f"That is about {round(lr['f_today_over_knee'],-2):,.0f}–{round(ew['f_today_over_knee'],-2):,.0f} times the knee frequency, in the sub-mHz band, not the nHz band. "
      f"Putting the knee there instead would need δ ≈ {sci(ew['delta_required_for_knee_at_this_ell_leading_order'])}–{sci(lr['delta_required_for_knee_at_this_ell_leading_order'])} rather than 10⁻³, "
      f"or ℓ ≈ {fmt(ew['ell_required_for_knee_at_registered_delta_m'],3)} m, which is excluded.")
    w(f"6. **Dark energy (DESI).** The relaxing fate is an empty RS de Sitter brane, with w = −1. Identifying it with today's acceleration needs ℓ ≈ {fmt(de['ell_needed_Mpc'],3)} Mpc, "
      f"{sci(de['ratio_to_38.6um'],1)} times the tabletop bound, or δ ≈ {sci(de['delta_needed_for_H_dS_at_38.6um_leading_order'])} (**conditional**). "
      "This is the Randall–Sundrum form of the cosmological-constant problem. HDBLAST currently makes **no** DESI prediction.")
    w(f"7. **Dark radiation is the cleanest direct constraint on a 5D blast remnant.** A bulk black hole or black brane (for example from the collapse fate) appears on the brane as C/a⁴. "
      f"ACT DR6 (N_eff = 2.86 ± 0.13) implies ρ_dr/ρ_γ ≲ {fmt(C['ACT_DR6_rho_dr_over_rho_gamma_upper'],2)}, and BBN 2024 implies ≲ {fmt(C['BBN_2024_rho_dr_over_rho_gamma_upper'],2)} "
      "(Gaussian approximations, **numerical**). This applies once a radiation era exists, which the registered model has not yet produced.")
    w("8. **Unchanged guardrails** (**numerical**, reproduced). The knee spectrum's ΔN_eff = "
      f"{sci(A['delta_neff']['DeltaNeff_exact'],3)}, which is exact via a Beta-function closed form and matches the archived ≈7×10⁻⁴. "
      f"Its Ω(25 Hz) = {sci(A['LVK_25Hz']['Omega_SBPL_25Hz'],3)} is {fmt(A['LVK_25Hz']['orders_of_magnitude_below_limit'],3)} orders below the LVK O1–O4a limit. Neither observation tests the knee.")
    w("9. **Not yet testable.** The CMB tilt (n_s ≈ 0.968–0.974), tensors (r < 0.034), JWST early galaxies and H0 all require a primordial perturbation spectrum or an expansion history after a radiation era. "
      "The registered model has neither yet. These are **open requirements**, not passes or failures.")
    w("")
    # --------------------------------------------------------------------------------------------
    w("## 2. What each observation would constrain in a braneworld / higher-dimensional-origin model")
    w("")
    w("| Observation | What it constrains in a braneworld origin model | Current value (snippet) | Status for HDBLAST |")
    w("|---|---|---|---|")
    w("| PTA background (NANOGrav, EPTA/InPTA, PPTA, CPTA, MPTA) | Spectral shape (knees, breaks, infrared tail), amplitude, isotropy and correlation pattern of any relic GW. Also any radion/KK phase transition at T ~ 0.1–1 GeV | NG15 A = 2.4e-15 (2023), 2.1e-15 with customized noise (2026); EPTA 2.5e-15; PPTA 2.0e-15; MPTA 4.8e-15 | Knee untested at likelihood level; the knee sits on the 1/yr blind spot; index trend unfavourable (conditional) |")
    w("| CMB N_eff (Planck, ACT DR6, SPT-3G) | Bulk Weyl dark radiation C/a⁴, KK/radion relics, GW energy density before recombination | N_eff = 2.86 ± 0.13 (ACT DR6) | ρ_dr/ρ_γ ≲ 0.016 (Gaussian approx.); binding once a radiation era exists |")
    w("| CMB tensors (BICEP/Keck, SPT-3G) | Energy scale of any inflation-like phase; RS high-energy enhancement of tensors | r < 0.034 (2025 combination) | No HDBLAST tensor prediction yet |")
    w("| CMB scalar tilt | Mechanism for adiabatic, nearly scale-invariant perturbations | n_s = 0.968 ± 0.003 (CMB), 0.973–0.974 with BAO | Open requirement |")
    w("| BBN | Expansion rate at ~1 MeV: ρ²/2λ term (needs λ^(1/4) ≫ MeV), dark radiation of either sign | ΔN_eff = −0.09 ± 0.28 (one 2024 configuration); braneworld DR −12.1% to +6.2% at 10 MeV (2017) | Automatically satisfied for tabletop-scale λ; DR bound open |")
    w("| DESI BAO + SNe | Late-time w(z); a braneworld gives w = −1 on an RS dS brane, phantom-like behaviour only with induced gravity (DGP-type) | w0waCDM preferred at 2.8–4.2σ (DR2); 3.2σ after DES-Dovekie; Lyα DR2 2026 consistent with ΛCDM | Outside the registered model unless a late sector is added |")
    w(f"| LVK stochastic background | High-frequency relic GWs; 25 Hz corresponds to horizon-scale production at T ≈ {sci(B['temperature_for_band_GeV']['LVK_25Hz']['fH1'],2)} GeV (f*/H* = 1) or {sci(B['temperature_for_band_GeV']['LVK_25Hz']['fH100'],2)} GeV (f*/H* = 100) | Ω(25 Hz) ≤ 2.0e-9 (2/3), ≤ 2.8e-9 (flat) | Knee 26 orders below; a TeV-scale blast would land in the sub-mHz band instead |")
    w("| GW propagation (GW170817, GWTC-3/4 sirens, graviton mass) | Leakage/damping into non-compact dimensions, bulk shortcuts, KK dispersion | D = 3.95 (+0.09/−0.07) (GWTC-3); D = 4.38 (+1.91/−1.01) (GWTC-4 dark sirens); m_g < 1.92e-23 eV | Bites only if ℓ or a crossover scale is cosmological, which tabletop tests already exclude |")
    w("| Tabletop inverse-square law | The AdS radius ℓ, hence M5, λ and the physical unit of the registered model | Yukawa range < 38.6 µm (2020); ℓ < 0.1 mm, M5 > 1e8 GeV (review) | Sets ρ^(1/4) ≳ 0.8–1.2 TeV for the static branch (conditional) |")
    w("| JWST early galaxies | Small-scale primordial power, early growth | z ≈ 14.4 galaxies; debated overabundance | Untestable until HDBLAST predicts perturbations |")
    w("| Hubble tension | Early expansion (dark radiation), late expansion | Local 73.50 ± 0.81 vs early 67.24 ± 0.35 (7.1σ) | No HDBLAST H0 prediction; ACT DR6 disfavours dark-radiation fixes |")
    w("")
    # --------------------------------------------------------------------------------------------
    w("## 3. The project's registered PTA knee test and the 2025–2026 data")
    w("")
    w("**The registered test** (project page `GPD-site-deploy/library/answers/checked/what-could-prove-hd-blast-wrong.json` and Zenodo record 17968738, 17 Dec 2025) is a smooth broken power law. "
      "The formula is copied from the archived `NEXTLEVEL_DeltaNeff_Audit_PHYSICAL.py`:")
    w("")
    w("    Omega(f) = Omega_k [((f/f_k)^(a1 D) + (f/f_k)^(a2 D))/2]^(-1/D),  f_k = 3.16e-8 Hz, Omega_k = 8.00e-9, a1 = +3, a2 = -2, D = 2")
    w("")
    w("Failure criteria stated by the project: a straight power law through the knee region; the bend fading under a full noise-aware analysis; teams disagreeing about where it sits; or SMBHBs explaining everything. "
      "A later 'BK' (Bessel-K, ν = 2) fingerprint gate, from the 26 Mar 2026 handoff, adds 'LowK' (≈2.96 nHz) and 'HighK' (≈33.1 nHz) knee families. "
      "The same handoff records that the IPTA DR2 common-process free spectrum and the archived NANOGrav HD-like export **fail** the frozen BK gate (project-internal; not re-verified here).")
    w("")
    w("### 3.1 Shape facts (exact/numerical, script)")
    w("")
    sh = A["shape"]
    w(f"- Local slope d ln Ω/d ln f → +{fmt(sh['slope_far_below_knee'],3)} far below the knee, {fmt(sh['slope_far_above_knee'],3)} far above, and {fmt(sh['slope_at_knee'],3)} at the knee. "
      f"The equivalent residual indices are γ = {fmt(sh['gamma_equiv_far_below'],3)}, {fmt(sh['gamma_equiv_far_above'],3)} and {fmt(sh['gamma_equiv_at_knee'],3)}. "
      "Wrong-formula control: reading a1 = +3 as the low-frequency slope, i.e. Ω ∝ f³, is **rejected**. The archived formula gives f² below the knee, matching the book's wording.")
    w("- A causal, horizon-limited source in the radiation era has a universal **f³** infrared tail (arXiv 2010.03568). "
      "The frozen **f²** therefore needs a stated mechanism if it is to be cosmological.")
    w("")
    w("### 3.2 Where the knee sits relative to each data set (numerical)")
    w("")
    w("| Data set | Span T (yr) | Resolution 1/T (nHz) | Knee in units of 1/T | Knee − 1/yr (bins) |")
    w("|---|---|---|---|---|")
    for n, d in kl.items():
        w(f"| {n} | {d['T_yr']} | {fmt(d['frequency_resolution_1_over_T_nHz'],3)} | {fmt(d['knee_in_units_of_1_over_T'],4)} | {fmt(d['knee_minus_1_per_yr_in_bins'],2)} |")
    w("")
    w(f"- NANOGrav 15-yr's 14 standard bins end at {fmt(kl['NANOGrav_15yr']['highest_standard_bin_nHz'],3)} nHz, so the knee lies **above** the standard band.")
    w("- Spans: NANOGrav 16.03 yr is adopted. It is consistent with 'nearly 16 years' in the snippet and with the project page's '14 frequencies from about 2 to 28 nHz'. The EPTA DR2new, MeerKAT and CPTA spans are from snippets.")
    w(f"- Project BK windows: the LowK narrow window is {fmt(bk['LowK_narrow_nHz']['width_nHz'],3)} nHz wide, "
      f"{fmt(bk['LowK_narrow_nHz']['width_over_resolution']['NANOGrav_15yr'],2)} of NANOGrav's frequency resolution. No current data set can resolve a knee position that finely. "
      "Sampling a free spectrum at 0.025 nHz spacing, as the handoff requests, yields strongly correlated values, not independent information. This is a general Fourier-resolution point, not taken from a cited source. "
      f"The HighK window {'**contains** 1/yr' if bk['HighK_narrow_nHz']['contains_1_per_yr'] else 'does not contain 1/yr'}, so it has the same timing-model degeneracy as the SBPL knee.")
    w("- Perturbed-parameter control: moving f_k by −10% / +10% lowers Ω at 1/yr by "
      f"{fmt(100*(1-A['perturbed_f_k_control']['0.9']['Omega_at_f_yr_over_frozen']),2)}% / {fmt(100*(1-A['perturbed_f_k_control']['1.1']['Omega_at_f_yr_over_frozen']),2)}%. "
      "The power at 1/yr alone barely pins the knee position. Locating the knee needs the bins on both sides of 1/yr, which is exactly where the timing-model fit absorbs power.")
    w("")
    w("### 3.3 Frozen curve vs a γ = 13/3 power law in the NANOGrav bins (numerical, descriptive)")
    w("")
    w("| Bin | f (nHz) | Ω_SBPL | local γ | Ω_SBPL / Ω_PL(A = 2.4e-15) | Ω_SBPL / Ω_PL(A = 2.1e-15, 2026) |")
    w("|---|---|---|---|---|---|")
    for r in ng:
        if r["bin"] in (1, 2, 3, 5, 8, 10, 12, 13, 14):
            w(f"| {r['bin']} | {fmt(r['f_nHz'],3)} | {r['Omega_SBPL']:.2e} | {fmt(r['local_gamma'],3)} | {fmt(r['ratio_SBPL_over_NG15_A_13_3'],2)} | {fmt(r['ratio_SBPL_over_NG15_CNM_A_fixed'],2)} |")
    w("")
    w("Reading: *if* the true background were the published fixed-index (γ = 13/3) fit, the frozen curve would carry only 3.5–15% of its power in the three lowest bins. "
      "This compares a model with a model summary, not with the data. Free-index fits have broad errors: EPTA DR2 γ = 4.19 (+0.73/−0.63), MeerKAT γ = 3.60 (+1.31/−0.89). "
      "The curves cross near 25 nHz. PTA sensitivity to red processes is concentrated at the lowest frequencies (general PTA knowledge, not a cited snippet). "
      "The frozen curve therefore stands or falls mainly on the low-bin power and the index, not on the knee itself. "
      "**This is not a likelihood analysis.** A proper test fits the frozen SBPL, with no free parameters, to the corrected free-spectrum KDEs or the PPL posterior, together with the timing-model transmission.")
    w("")
    w("### 3.4 Do the 2025–2026 releases bear on the knee?")
    w("")
    w("| Result (year) | Bears on | Direction for the frozen knee |")
    w("|---|---|---|")
    w("| NANOGrav running of the spectral index (2025) | curvature within 2–28 nHz | no curvature required (β consistent with 0, BF 0.69): **inconclusive** |")
    w("| NANOGrav piecewise power-law reconstruction (2026) | broken/doubly-broken spectra, BMA | the right comparison target; numbers not visible in snippets: **inconclusive** |")
    w("| NANOGrav customized chromatic noise (2026) | amplitude, HD significance | amplitude down to 2.1e-15; frozen height now 1.31× in Ω, within errors: **mildly unfavourable** |")
    w("| NANOGrav 15-yr free-spectrum correction / erratum (2026) | the posteriors the project screened | project screens need a rerun on corrected KDEs: **method** |")
    w("| EPTA DR2 improved noise models (2025) | index, amplitude | toward γ = 13/3, lower amplitude: **unfavourable to the γ = 3 low branch** |")
    w("| MeerKAT PTA 4.5 yr (2024/2025) | amplitude, index at high cadence | amplitude 4–10× the frozen height in Ω; γ = 3 allowed: **cross-team spread** |")
    w("| Five-PTA combination (Dec 2025) | amplitude, power-law exponent | power-law only; could host a knee test: **method** |")
    w("| NANOGrav 20-yr, IPTA DR3 | direct test | not released in the searches: **pending** |")
    w("")
    w("### 3.5 Frequency–temperature map (numerical, conditional on f*/H* = 1)")
    w("")
    tk = B["temperature_for_knee"]
    w(f"A horizon-scale feature at f_k corresponds to production at T ≈ {fmt(tk['106.75'],2)}–{fmt(tk['10.75'],2)} GeV (g* from 106.75 down to 10.75), i.e. the QCD era. "
      f"The script's map reproduces the commonly quoted coefficient: f ≈ {fmt(B['frequency_map_control']['f0_at_T_1GeV_gstar100_Hz'],3)} Hz at T = 1 GeV, g* = 100. "
      "A five-dimensional origin of a feature at that temperature has published precedent: a radion/KK phase transition at MeV–GeV (Phys. Rev. D 108, 095017). "
      "That is a different mechanism from the registered HDBLAST shell.")
    w("")
    # --------------------------------------------------------------------------------------------
    w("## 4. Scale-setting the registered model with tabletop gravity (conditional)")
    w("")
    ef = B["exact_model_facts"]
    w(f"**Exact-verified** (sympy): W(1) = {ef['W(1)']}, U(1) = {ef['U(1)']}, 1/ℓ = √(−U(1)/6) = {ef['1/ell=sqrt(-U(1)/6)']}, so ℓ = {ef['ell']}. "
      f"The critical tension 2W(1) = {ef['sigma_crit=2W(1)']} satisfies σ_c/6 = 1/ℓ ({'verified' if ef['sigma_crit/6 == 1/ell'] else 'FAILED'}). This agrees with the theory sweep's cross-check. "
      f"**Numerical:** H·ℓ = 9/ρ_b = {fmt(B['registered_H_ell']['H_ell_from_rho_b'],6)}, and 9√(H²) agrees to {B['registered_H_ell']['rel_diff']:.1e}. "
      f"The O(δ²) expansion of H² matches the registered value to {B['registered_H_ell']['rel_diff_expansion2_vs_registered']:.1e} (relative). "
      f"ρ/λ ≈ (Hℓ)²/2 = {fmt(B['registered_H_ell']['rho_over_lambda_estimate=(H ell)^2/2'],2)}, so the static branch is in the low-energy RS regime.")
    w("")
    w("| ℓ bound used | M5 (GeV) | λ^(1/4) (TeV) | ρ_vac^(1/4) (TeV) | T if thermalised (TeV) | f today, f*/H* = 1 (Hz) | f / f_knee |")
    w("|---|---|---|---|---|---|---|")
    for lab, s in ts.items():
        w(f"| {PRETTY.get(lab, lab)} | {sci(s['M5_GeV'])} | {fmt(s['lambda_quarter_GeV']/1e3,3)} | {fmt(s['rho_vac_quarter_GeV']/1e3,3)} | {fmt(s['T_if_fully_thermalised_GeV(g*=106.75)']/1e3,3)} | {sci(s['f_today_horizon_scale_Hz(f*/H*=1)'])} | {round(s['f_today_over_knee'],-1):,.0f} |")
    w("")
    w("All energies scale as ℓ^(−1/2). A smaller ℓ only raises them, so the table gives **lower bounds** on the energy scales and the frequency. "
      "Assumptions: RS2 with a Z2-doubled bulk, M_Pl² = M5³ℓ, ρ_vac = 3M_Pl²H² on the brane, O((Hℓ)²) corrections ignored, full and instantaneous thermalisation, f*/H* = 1, standard radiation era afterwards. "
      "The Eöt-Wash 38.6 µm figure is a **Yukawa** range, used only as a proxy. The RS correction is a power law (1 + 2ℓ²/3r²), and the 2026 inverse-square-law review has the proper power-law limits. "
      "The review statement 'ℓ < 0.1 mm ⇒ M5 > 10⁸ GeV' is reproduced (control).")
    w("")
    w("**Interpretation (conditional).** If the registered δ = 10⁻³ static branch is the pre-radiation state, and ℓ obeys tabletop bounds, the natural relic-GW frequency lies in the sub-mHz band. "
      "A PTA knee would then have to come from later physics, not from the blast's horizon scale. "
      "This does not falsify HDBLAST. It shows that the PTA knee and the registered 5D model are currently **disconnected**, and any link needs an explicit mechanism.")
    w("")
    # --------------------------------------------------------------------------------------------
    w("## 5. Annotated bibliography")
    w("")
    w("Each entry gives the URL, the year, the finding as reported in the search snippet, its relevance to HDBLAST, and the implication label. **Access: every entry is snippet-only; no full text was read.**")
    w("")
    for sec, title in SECTIONS.items():
        rs = [r for r in records if r["section"] == sec]
        if not rs:
            continue
        w(f"### {title}")
        w("")
        for r in rs:
            w(f"**{r['title']}** ({r['year']}). {r['venue']}.  ")
            w(f"URL: {r['url']}  ")
            w(f"*Access:* {r['access']} · *Implication:* {r['implication']} · *Queries:* {', '.join('Q'+str(q) for q in r['queries'])}")
            w("")
            w(f"*Finding.* {r['finding']}")
            w("")
            w(f"*Relevance to HDBLAST.* {r['relevance']}")
            if r.get("attribution"):
                w("")
                w(f"*Attribution note.* {r['attribution']}")
            w("")
    # --------------------------------------------------------------------------------------------
    w("## 6. Limitations")
    w("")
    w("- Snippet-level evidence only. Snippets can garble numbers or merge several papers, and the entries flag the cases noticed. Nothing here substitutes for reading the papers.")
    w("- Search coverage is incomplete, and the search budget ran out: three planned queries (IPTA DR3 status, the 'What is the source of the PTA GW signal' paper, EPTA DR2 implications) were refused. Absence of a 2026 IPTA DR3 or NANOGrav 20-yr background result is 'not found', not 'does not exist'.")
    w("- The knee comparisons in §3 are descriptive, not likelihood analyses. They use published power-law summaries, not the free-spectrum posteriors.")
    w("- The scale-setting in §4 is order of magnitude and conditional on RS2 low-energy relations and on the thermalisation assumptions listed there. The registered model has not produced a radiation era.")
    w("- The Gaussian ΔN_eff upper limits are built from central values ± 1σ. They are not the collaborations' own 95% limits.")
    w("")
    # --------------------------------------------------------------------------------------------
    w("## 7. Search log")
    w("")
    for i, q in enumerate(SO.QUERIES, 1):
        w(f"- Q{i}: `{q}`")
    for q in SO.REFUSED_QUERIES:
        w(f"- refused (budget): `{q}`")
    w("")
    return "\n".join(L) + "\n"


def main():
    selftest = "--selftest" in sys.argv
    recs = SO.S
    if selftest:
        bad1 = copy.deepcopy(recs)
        bad1[0]["implication"] = "proves"
        bad2 = copy.deepcopy(recs)
        bad2[1]["url"] = "http://example.com"
        bad3 = copy.deepcopy(recs)
        bad3[2]["year"] = "2019"          # earlier than its arXiv id 2601 -> must fail
        bad4 = copy.deepcopy(recs)
        bad4.append(copy.deepcopy(recs[0]))  # duplicate URL
        results = {name: bool(validate(r, SO.QUERIES)) for name, r in
                   [("bad_label", bad1), ("bad_url", bad2), ("bad_year", bad3), ("duplicate", bad4)]}
        results["clean_passes"] = not validate(recs, SO.QUERIES)
        J = json.loads((HERE / "knee_and_scale_checks.json").read_text())
        J2 = copy.deepcopy(J)
        J2["A_knee_benchmark"]["LVK_25Hz"]["Omega_SBPL_25Hz"] *= 10
        results["number_check_detects_corruption"] = bool(number_checks(J2)[0])
        print(json.dumps(results, indent=1))
        ok = all(results.values())
        (HERE / "builder_selftest.json").write_text(json.dumps({"results": results, "all_pass": ok}, indent=1))
        print("selftest all_pass:", ok)
        sys.exit(0 if ok else 1)
    errs = validate(recs, SO.QUERIES)
    if errs:
        print("VALIDATION FAILED:\n" + "\n".join(errs))
        sys.exit(1)
    J = json.loads((HERE / "knee_and_scale_checks.json").read_text())
    if not J.get("all_checks_pass"):
        print("knee_and_scale_checks.json does not report all checks passing")
        sys.exit(1)
    bad, n = number_checks(J)
    if bad:
        print("NUMBER CHECKS FAILED:\n" + "\n".join(bad))
        sys.exit(1)
    J["_n_number_checks"] = n
    scan = json.loads((HERE / "project_scan_1yr.json").read_text())
    if not scan.get("positive_control_passed"):
        print("project scan positive control failed")
        sys.exit(1)
    J["_scan"] = scan
    md = render(recs, J)
    (HERE.parent / "observations.md").write_text(md)
    counts = {}
    for r in recs:
        counts[r["implication"]] = counts.get(r["implication"], 0) + 1
    (HERE / "observations_sources.json").write_text(json.dumps(
        {"n_sources": len(recs), "n_queries": len(SO.QUERIES), "refused_queries": SO.REFUSED_QUERIES,
         "implication_counts": counts, "sources": recs, "queries": SO.QUERIES, "number_checks": n}, indent=1, ensure_ascii=False))
    print(f"wrote observations.md ({len(md)} chars); {len(recs)} sources; {len(SO.QUERIES)} queries; {n} number checks; counts {counts}")


if __name__ == "__main__":
    main()
