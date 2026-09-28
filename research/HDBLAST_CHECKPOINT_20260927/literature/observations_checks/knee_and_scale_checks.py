#!/usr/bin/env python3
"""Observation-facing checks for the HDBLAST literature sweep 2 (observations).

Three blocks, all cheap (< 10 s, 1 core):

A. The registered PTA "knee" benchmark (Zenodo 17968738, 17 Dec 2025):
   smooth broken power law (SBPL), exactly as in the archived script
   NEXTLEVEL_DeltaNeff_Audit_PHYSICAL.py:
       Omega(f) = Omega_k * [((f/f_k)^(a1*D) + (f/f_k)^(a2*D))/2]^(-1/D)
   with f_k = 3.16e-8 Hz, Omega_k = 8.00e-9, a1 = +3, a2 = -2, D = 2.
   - local slopes, equivalent PTA timing-residual index gamma, strain slope
   - conversion to characteristic strain and comparison with published PTA
     amplitudes (numbers transcribed from search-result snippets; see
     observations.md for the URLs)
   - frequency-bin placement of the knee for several data spans, distance
     from 1/yr (where pulsar sky-position fits remove sensitivity)
   - Delta N_eff of the SBPL: archived trapezoid method reproduced, plus an
     exact closed form (Beta function) as an independent control
   - Omega at 25 Hz versus LVK O1-O4a upper limits

B. Braneworld scale-setting translation of the registered static "+1 branch"
   (CONDITIONAL, order of magnitude): model units have kappa_5^2 = 1 and
   AdS radius ell = 3/W(1) = 9 at phi = +1.  H*ell is scale free.  If the
   physical ell is bounded by torsion-balance tests of Newton's law, the
   brane vacuum energy, M5, brane tension, a thermalisation temperature and
   the present-day horizon-scale GW frequency follow.

C. Dark-radiation bookkeeping: Delta N_eff bounds -> rho_dr / rho_gamma.

Every literature number used here is listed in LIT_INPUTS with its source URL.
Outputs: knee_and_scale_checks.json (same folder).
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np
import sympy as sp
from scipy import integrate

OUT = Path(__file__).with_suffix(".json")

# ----------------------------------------------------------------------------
# Physical constants (CODATA-level values; precision far beyond what is needed)
# ----------------------------------------------------------------------------
HBAR_eVs = 6.582119569e-16          # eV s
HBARC_eVm = 1.973269804e-7          # eV m
KB_eV_per_K = 8.617333262e-5        # eV / K
MPC_m = 3.0856775814913673e22       # m
YEAR_s = 365.25 * 86400.0           # Julian year, s
F_YR = 1.0 / YEAR_s                 # Hz
M_PL_RED_eV = 2.435323e27           # reduced Planck mass, eV
T0_K = 2.7255
G_S0 = 43.0 / 11.0                  # entropy dof today

# ----------------------------------------------------------------------------
# Registered knee benchmark (project record; archived defaults)
# ----------------------------------------------------------------------------
KNEE = dict(f_k=3.16e-8, Omega_k=8.00e-9, a1=3.0, a2=-2.0, D=2.0, h=0.674)

# Project-internal "BK" knee families (HDBLAST handoff 26 Mar 2026, sec. 13);
# used only to compare window widths with frequency resolution 1/T.
BK_FAMILIES = {
    "LowK_narrow_nHz": [2.815, 3.048],
    "LowK_broad_nHz": [2.539, 3.275],
    "HighK_narrow_nHz": [28.504, 38.719],
    "HighK_broad_nHz": [25.148, 43.801],
}

# ----------------------------------------------------------------------------
# Literature inputs (transcribed from search-result snippets; URLs in md)
# ----------------------------------------------------------------------------
LIT_INPUTS = {
    "NG15_A_13_3": dict(A=2.4e-15, gamma=13 / 3, note="NANOGrav 15-yr, median, 90% CI +0.7/-0.6",
                        url="https://arxiv.org/abs/2306.16213"),
    "NG15_CNM_A_fixed": dict(A=2.1e-15, gamma=13 / 3, note="NANOGrav 15-yr customized chromatic noise, +0.6/-0.5 (2026)",
                             url="https://arxiv.org/abs/2606.28554"),
    "EPTA_DR2new_A_13_3": dict(A=2.5e-15, gamma=13 / 3, note="EPTA DR2new+InPTA, +-0.7",
                               url="https://arxiv.org/abs/2306.16214"),
    "EPTA_DR2full_HD_free": dict(A=10 ** -14.54, gamma=4.19, note="EPTA DR2 full, HD process, log10A=-14.54(+0.28/-0.41), gamma=4.19(+0.73/-0.63)",
                                 url="https://arxiv.org/abs/2306.16214"),
    "PPTA_DR3_A_13_3": dict(A=2.0e-15, gamma=13 / 3, note="PPTA DR3 common-spectrum process, +-0.2",
                            url="https://arxiv.org/abs/2306.16215"),
    "MPTA_4p5yr_A_13_3": dict(A=4.8e-15, gamma=13 / 3, note="MeerKAT PTA 4.5 yr, h_c,yr at alpha=-2/3, +0.8/-0.9",
                              url="https://arxiv.org/abs/2412.01153"),
    "MPTA_4p5yr_A_free": dict(A=7.5e-15, gamma=3.0 - 2 * (-0.26), note="MeerKAT PTA 4.5 yr, h_c,yr at alpha=-0.26, +0.8/-0.9",
                              url="https://arxiv.org/abs/2412.01153"),
}
DATA_SPANS_yr = {
    # NANOGrav 15-yr: 14 frequencies from ~2 to ~28 nHz (project page, and
    # 'nearly 16 years' in the NANOGrav release snippet); 16.03 yr adopted.
    "NANOGrav_15yr": 16.03,
    "EPTA_DR2new": 10.3,       # snippet: 'most recent 10.3 yr'
    "MPTA_4p5yr": 4.5,         # snippet: '4.5 years'
    "CPTA_DR1": 3.0,           # snippet: 'close to three years'
}
N_BINS_NG15 = 14

LVK_O4A_LIMITS_25Hz = {"alpha_2over3": 2.0e-9, "flat": 2.8e-9}   # arXiv 2508.20721 snippet
EOTWASH_YUKAWA_RANGE_m = 38.6e-6    # Lee et al. 2020, gravitational-strength Yukawa, 95% CL
LIVING_REVIEW_ELL_m = 1.0e-4        # Maartens & Koyama review statement 'ell < 0.1 mm'
ACT_DR6_NEFF = (2.86, 0.13)         # ACT DR6 extended models
BBN_2024_DNEFF = (-0.09, 0.28)      # 2024 BBN baryon abundance update (one configuration)
NEFF_SM = 3.044


# ----------------------------------------------------------------------------
# SBPL helpers
# ----------------------------------------------------------------------------
def omega_sbpl(f, f_k, Omega_k, a1, a2, D):
    x = np.asarray(f, dtype=float) / f_k
    return Omega_k * ((x ** (a1 * D) + x ** (a2 * D)) / 2.0) ** (-1.0 / D)


def slope_sbpl(f, f_k, a1, a2, D, **_):
    """Analytic d ln Omega / d ln f."""
    x = np.asarray(f, dtype=float) / f_k
    t1 = x ** (a1 * D)
    t2 = x ** (a2 * D)
    return -(a1 * t1 + a2 * t2) / (t1 + t2)


def H0_s(h):
    return h * 1.0e5 / MPC_m


def hc_from_omega(Omega, f, h):
    return math.sqrt(3.0 * H0_s(h) ** 2 * Omega / (2.0 * math.pi ** 2 * f ** 2))


def omega_from_powerlaw(A, gamma, f, h):
    hc = A * (f / F_YR) ** ((3.0 - gamma) / 2.0)
    return 2.0 * math.pi ** 2 * f ** 2 * hc ** 2 / (3.0 * H0_s(h) ** 2)


def omega_gamma(h, T=T0_K):
    return 2.472e-5 * (T / 2.7255) ** 4 / h ** 2


NEFF_CONV = (8.0 / 7.0) * (11.0 / 4.0) ** (4.0 / 3.0)


def block_A():
    k = KNEE
    out = {}
    fk, Ok, a1, a2, D, h = k["f_k"], k["Omega_k"], k["a1"], k["a2"], k["D"], k["h"]
    # A1. Shape facts ---------------------------------------------------------
    s_lo = float(slope_sbpl(fk * 1e-6, fk, a1, a2, D))
    s_hi = float(slope_sbpl(fk * 1e6, fk, a1, a2, D))
    s_k = float(slope_sbpl(fk, fk, a1, a2, D))
    # finite-difference control of the analytic slope
    fd = []
    for xf in [0.1, 0.5, 1.0, 2.0]:
        f0 = xf * fk
        eps = 1e-6
        num = (math.log(omega_sbpl(f0 * (1 + eps), fk, Ok, a1, a2, D)) - math.log(omega_sbpl(f0 * (1 - eps), fk, Ok, a1, a2, D))) / (math.log1p(eps) - math.log1p(-eps))
        fd.append(abs(num - float(slope_sbpl(f0, fk, a1, a2, D))))
    out["shape"] = {
        "Omega_at_fk_over_Omega_k": float(omega_sbpl(fk, fk, Ok, a1, a2, D) / Ok),
        "slope_far_below_knee": s_lo, "slope_far_above_knee": s_hi, "slope_at_knee": s_k,
        "gamma_equiv_far_below": 5 - s_lo, "gamma_equiv_far_above": 5 - s_hi, "gamma_equiv_at_knee": 5 - s_k,
        "strain_alpha_far_below": (s_lo - 2) / 2, "strain_alpha_far_above": (s_hi - 2) / 2,
        "analytic_vs_finite_difference_slope_max_abs_diff": max(fd),
        "wrong_formula_control": {
            "claim_tested": "naive reading 'Omega ~ f^{+a1} below the knee'",
            "naive_low_slope": a1, "actual_low_slope": s_lo,
            "naive_reading_rejected": abs(s_lo - a1) > 0.5,
        },
        "convention_note": "gamma is the timing-residual PSD index, Omega_gw ~ f^(5-gamma); SMBHB circular GW-driven gamma=13/3",
    }
    # A2. Knee location --------------------------------------------------------
    out["knee_location"] = {
        "f_k_Hz": fk, "f_yr_Hz": F_YR, "f_k_over_f_yr": fk / F_YR,
        "relative_offset_from_1_per_yr": (fk - F_YR) / F_YR,
        "per_dataset": {},
    }
    for name, T in DATA_SPANS_yr.items():
        Ts = T * YEAR_s
        d = {"T_yr": T, "frequency_resolution_1_over_T_nHz": 1e9 / Ts,
             "knee_in_units_of_1_over_T": fk * Ts, "one_per_yr_in_units_of_1_over_T": F_YR * Ts,
             "knee_minus_1_per_yr_in_bins": (fk - F_YR) * Ts}
        if name == "NANOGrav_15yr":
            d["highest_standard_bin_nHz"] = 1e9 * N_BINS_NG15 / Ts
            d["knee_above_highest_standard_bin"] = fk > N_BINS_NG15 / Ts
        out["knee_location"]["per_dataset"][name] = d
    # BK family windows vs resolution
    bk = {}
    for fam, (lo, hi) in BK_FAMILIES.items():
        bk[fam] = {"width_nHz": hi - lo,
                   "width_over_resolution": {n: (hi - lo) * 1e-9 * T * YEAR_s for n, T in DATA_SPANS_yr.items()},
                   "contains_1_per_yr": lo * 1e-9 <= F_YR <= hi * 1e-9}
    out["project_BK_family_windows"] = bk
    # hypothetical future span (NANOGrav ~20 yr data set, not yet analysed for the background)
    T20 = 20.0 * YEAR_s
    out["hypothetical_20yr_span"] = {"resolution_nHz": 1e9 / T20, "one_per_yr_bin": F_YR * T20,
                                     "knee_bin": fk * T20}
    # A3. Local slope across NG15 bins ---------------------------------------
    Ts = DATA_SPANS_yr["NANOGrav_15yr"] * YEAR_s
    fb = np.arange(1, N_BINS_NG15 + 1) / Ts
    sl = slope_sbpl(fb, fk, a1, a2, D)
    om_sb = omega_sbpl(fb, fk, Ok, a1, a2, D)
    coeff = np.polyfit(np.log(fb), np.log(om_sb), 1)
    coeff5 = np.polyfit(np.log(fb[:5]), np.log(om_sb[:5]), 1)
    rows = []
    for i, f in enumerate(fb):
        r = {"bin": i + 1, "f_nHz": f * 1e9, "Omega_SBPL": float(om_sb[i]), "local_slope": float(sl[i]),
             "local_gamma": float(5 - sl[i])}
        for key in ["NG15_A_13_3", "NG15_CNM_A_fixed"]:
            L = LIT_INPUTS[key]
            r[f"ratio_SBPL_over_{key}"] = float(om_sb[i] / omega_from_powerlaw(L["A"], L["gamma"], f, h))
        rows.append(r)
    out["NG15_bins"] = {
        "rows": rows,
        "effective_unweighted_loglog_slope_all14": float(coeff[0]), "effective_gamma_all14": float(5 - coeff[0]),
        "effective_unweighted_loglog_slope_bins1to5": float(coeff5[0]), "effective_gamma_bins1to5": float(5 - coeff5[0]),
        "caveat": "descriptive only; not a likelihood fit; PTA sensitivity is dominated by the lowest bins",
    }
    # A4. Amplitude comparison at 1/yr ---------------------------------------
    om_yr = float(omega_sbpl(F_YR, fk, Ok, a1, a2, D))
    hc_yr = hc_from_omega(om_yr, F_YR, h)
    amp = {"Omega_SBPL_at_f_yr": om_yr, "hc_SBPL_at_f_yr_h0674": hc_yr,
           "hc_SBPL_at_f_yr_h0735": hc_from_omega(om_yr, F_YR, 0.735),
           "note_h": "Omega_k is stated as Omega_gw with h=0.674 in the archived script; Omega_gw h^2 is the H0-independent quantity",
           "Omega_k_h2": Ok * h ** 2, "per_PTA": {}}
    for key, L in LIT_INPUTS.items():
        om_pl = omega_from_powerlaw(L["A"], L["gamma"], F_YR, h)
        amp["per_PTA"][key] = {"A": L["A"], "gamma": L["gamma"], "Omega_powerlaw_at_f_yr": om_pl,
                               "SBPL_over_PTA_Omega_at_f_yr": om_yr / om_pl, "SBPL_hc_over_A": hc_yr / L["A"],
                               "note": L["note"], "url": L["url"]}
    # calibration control: NANOGrav A=2.4e-15 at f_yr -> Omega ~ 8e-9 (the frozen Omega_k)
    amp["calibration_control"] = {
        "Omega_from_NG15_A_at_f_yr": omega_from_powerlaw(2.4e-15, 13 / 3, F_YR, h),
        "matches_Omega_k_within_5pct": abs(omega_from_powerlaw(2.4e-15, 13 / 3, F_YR, h) / Ok - 1) < 0.05,
    }
    out["amplitude_at_1_per_yr"] = amp
    # A5. Delta N_eff -----------------------------------------------------------
    s = -a2 / ((a1 - a2) * D)
    q = 1.0 / D - s
    mp.mp.dps = 30
    beta = mp.beta(s, q)
    I_exact = float(Ok * mp.power(2, 1 / D) * beta / ((a1 - a2) * D))
    # adaptive quadrature in u = ln x over a wide finite range
    g = lambda u: float(omega_sbpl(fk * math.exp(u), fk, Ok, a1, a2, D))
    I_quad, I_quad_err = integrate.quad(g, -60, 60, points=[0.0], limit=400, epsabs=0, epsrel=1e-12)
    trap = {}
    for n in [50_000, 100_000, 200_000, 400_000]:
        lnf = np.linspace(math.log(1e-12), math.log(1e4), n)
        trap[str(n)] = float(np.trapezoid(omega_sbpl(np.exp(lnf), fk, Ok, a1, a2, D), lnf))
    Og = omega_gamma(h)
    out["delta_neff"] = {
        "closed_form": "int Omega dln f = Omega_k 2^(1/D) B(s, 1/D - s) / ((a1-a2) D), s = -a2/((a1-a2) D)",
        "s": s, "q": q, "Beta": float(beta),
        "integral_exact": I_exact, "integral_quad": I_quad, "quad_error_estimate": I_quad_err,
        "integral_trapezoid_archived_band_by_n": trap,
        "rel_diff_quad_vs_exact": abs(I_quad / I_exact - 1),
        "rel_diff_trap200k_vs_exact": abs(trap["200000"] / I_exact - 1),
        "Omega_gamma": Og, "conversion_factor": NEFF_CONV,
        "DeltaNeff_exact": NEFF_CONV * I_exact / Og,
        "DeltaNeff_archived_method_n200k": NEFF_CONV * trap["200000"] / Og,
        "archived_value_quoted_in_project": "about 7e-4",
        "ACT_DR6_95pct_gaussian_upper_DeltaNeff": ACT_DR6_NEFF[0] + 1.96 * ACT_DR6_NEFF[1] - NEFF_SM,
        "caveat": "applies only if this spectrum existed as GW radiation before BBN/CMB; a late astrophysical origin would not count",
    }
    # A6. LVK band ----------------------------------------------------------------
    om25 = float(omega_sbpl(25.0, fk, Ok, a1, a2, D))
    out["LVK_25Hz"] = {"Omega_SBPL_25Hz": om25, "archived_value_quoted": "about 2.3e-35",
                       "limits_O1_O4a": LVK_O4A_LIMITS_25Hz,
                       "orders_of_magnitude_below_limit": math.log10(LVK_O4A_LIMITS_25Hz["alpha_2over3"] / om25)}
    # A7. Perturbed-parameter control: shift f_k by +-10% ------------------------
    pert = {}
    for fac in [0.9, 1.1]:
        pert[str(fac)] = {"Omega_at_f_yr_over_frozen": float(omega_sbpl(F_YR, fk * fac, Ok, a1, a2, D) / om_yr),
                          "knee_in_NG15_bins": fk * fac * DATA_SPANS_yr["NANOGrav_15yr"] * YEAR_s}
    out["perturbed_f_k_control"] = pert
    return out


# ----------------------------------------------------------------------------
# Block B: braneworld scale-setting translation (conditional)
# ----------------------------------------------------------------------------
REG = dict(delta=0.001, c=0.5975949350280132, phi_b=0.9999159473169134,
           rho_b=129.9247628497, H2=5.92401479433e-5)


def exact_rs_facts():
    p = sp.symbols("phi")
    W = 1 - p + p ** 3 / 3
    U = sp.Rational(1, 2) * sp.diff(W, p) ** 2 - sp.Rational(2, 3) * W ** 2
    W1 = sp.simplify(W.subs(p, 1))
    U1 = sp.simplify(U.subs(p, 1))
    inv_ell = sp.sqrt(-U1 / 6)
    return {"W(1)": str(W1), "U(1)": str(U1), "1/ell=sqrt(-U(1)/6)": str(sp.nsimplify(inv_ell)),
            "ell": str(sp.nsimplify(1 / inv_ell)), "sigma_crit=2W(1)": str(2 * W1),
            "sigma_crit/6 == 1/ell": bool(sp.simplify(2 * W1 / 6 - inv_ell) == 0),
            "W_phi(1)": str(sp.diff(W, p).subs(p, 1))}


def H_star_s(T_eV, gstar):
    return math.sqrt(math.pi ** 2 * gstar / 90.0) * T_eV ** 2 / M_PL_RED_eV / HBAR_eVs


def f_today_Hz(T_eV, gstar, gs=None, f_over_H=1.0):
    gs = gstar if gs is None else gs
    T0_eV = KB_eV_per_K * T0_K
    return f_over_H * H_star_s(T_eV, gstar) * (G_S0 / gs) ** (1 / 3) * T0_eV / T_eV


def T_for_f_today(f, gstar):
    # f is linear in T at fixed g*
    return f / f_today_Hz(1.0, gstar)


def block_B():
    out = {"exact_model_facts": exact_rs_facts()}
    ell_model = 9.0
    d, c = REG["delta"], REG["c"]
    Hl_rho = ell_model / REG["rho_b"]
    Hl_H2 = ell_model * math.sqrt(REG["H2"])
    H2_lo = d * (1 + c) / 27
    H2_2 = H2_lo + d ** 2 * ((1 + c) ** 2 / 36 - c ** 2 / 384)
    out["registered_H_ell"] = {
        "H_ell_from_rho_b": Hl_rho, "H_ell_from_H2": Hl_H2, "rel_diff": abs(Hl_rho / Hl_H2 - 1),
        "H2_expansion_O_delta": H2_lo, "H2_expansion_O_delta2": H2_2,
        "rel_diff_expansion2_vs_registered": abs(H2_2 / REG["H2"] - 1),
        "rho_over_lambda_estimate=(H ell)^2/2": Hl_rho ** 2 / 2,
        "note": "rel. diff of the O(delta^2) expansion should be O(delta^2) ~ 1e-6 or smaller (control)",
    }
    # perturbed-parameter control: c -> 1.01 c changes H ell at O(delta)
    out["perturbed_c_control"] = {
        "H_ell_leading_order_c": math.sqrt(3 * d * (1 + c)),
        "H_ell_leading_order_1.01c": math.sqrt(3 * d * (1 + 1.01 * c)),
    }
    scen = {}
    gstar = 106.75
    fk = KNEE["f_k"]
    for label, ell_m in [("EotWash2020_Yukawa_range_38.6um_proxy", EOTWASH_YUKAWA_RANGE_m),
                         ("LivingReview_0.1mm", LIVING_REVIEW_ELL_m)]:
        inv_ell_eV = HBARC_eVm / ell_m
        H_eV = Hl_rho * inv_ell_eV
        rho_vac = 3 * M_PL_RED_eV ** 2 * H_eV ** 2
        M5 = (M_PL_RED_eV ** 2 * inv_ell_eV) ** (1 / 3)
        lam = 6 * M_PL_RED_eV ** 2 * inv_ell_eV ** 2
        T_th = (30 * rho_vac / (math.pi ** 2 * gstar)) ** 0.25
        f0 = f_today_Hz(T_th, gstar)
        scen[label] = {
            "ell_m": ell_m, "1/ell_eV": inv_ell_eV, "H_brane_eV": H_eV,
            "H_brane_over_H0(h=0.674)": H_eV / (HBAR_eVs * H0_s(0.674)),
            "rho_vac_quarter_GeV": rho_vac ** 0.25 / 1e9, "M5_GeV": M5 / 1e9,
            "lambda_quarter_GeV": lam ** 0.25 / 1e9,
            "T_if_fully_thermalised_GeV(g*=106.75)": T_th / 1e9,
            "f_today_horizon_scale_Hz(f*/H*=1)": f0,
            "f_today_over_knee": f0 / fk,
            "delta_required_for_knee_at_this_ell_leading_order": d * (fk / f0) ** 4,
            "ell_required_for_knee_at_registered_delta_m": ell_m * (f0 / fk) ** 2,
        }
    out["tabletop_scale_setting"] = scen
    # if the relaxed RS de Sitter brane were today's dark energy
    h, OL = 0.674, 0.685
    H_dS = H0_s(h) * math.sqrt(OL)
    ell_needed = Hl_rho * 2.99792458e8 / H_dS
    out["if_brane_dS_is_todays_dark_energy"] = {
        "H_dS_s^-1": H_dS, "ell_needed_m": ell_needed, "ell_needed_Mpc": ell_needed / MPC_m,
        "ratio_to_38.6um": ell_needed / EOTWASH_YUKAWA_RANGE_m,
        "assumptions": "H_dS = H0 sqrt(Omega_L), h=0.674, Omega_L=0.685; registered H*ell",
        # leading order H ell = sqrt(3 delta (1+c)) -> delta needed for H = H_dS at ell = ell_max
        "delta_needed_for_H_dS_at_38.6um_leading_order": (H_dS * EOTWASH_YUKAWA_RANGE_m / 2.99792458e8) ** 2 / (3 * (1 + c)),
        "delta_needed_for_H_dS_at_0.1mm_leading_order": (H_dS * LIVING_REVIEW_ELL_m / 2.99792458e8) ** 2 / (3 * (1 + c)),
    }
    # control: textbook coefficient f0 ~ 1.65e-7 Hz at T=1 GeV, g*=g_s=100
    f_1GeV = f_today_Hz(1e9, 100.0)
    out["frequency_map_control"] = {"f0_at_T_1GeV_gstar100_Hz": f_1GeV,
                                    "agrees_with_commonly_quoted_1.65e-7_within_2pct": abs(f_1GeV / 1.65e-7 - 1) < 0.02}
    # which temperature puts f*/H*=1 at the knee
    out["temperature_for_knee"] = {str(gs): T_for_f_today(fk, gs) / 1e9 for gs in [10.75, 17.25, 61.75, 106.75]}
    out["temperature_for_knee_units"] = "GeV, g* = g_s given by key, f*/H* = 1"
    # other bands, f*/H* = 1 and f*/H* = 100 (sub-horizon source), g* = 106.75
    out["temperature_for_band_GeV"] = {
        band: {"fH1": T_for_f_today(fb, 106.75) / 1e9, "fH100": T_for_f_today(fb, 106.75) / 1e9 / 100}
        for band, fb in [("LVK_25Hz", 25.0), ("space_3mHz", 3e-3)]}
    # control check M5 against the review statement M5 > 1e8 GeV at ell = 0.1 mm
    out["M5_control_vs_review_1e8GeV"] = scen["LivingReview_0.1mm"]["M5_GeV"] > 1e8
    out["label"] = "CONDITIONAL, order of magnitude: RS2 low-energy relations M_Pl^2 = M5^3 ell (Z2-doubled bulk), rho_vac = 3 M_Pl^2 H^2, instantaneous thermalisation, f*/H*=1, standard radiation era afterwards; O((H ell)^2) corrections to M_Pl ignored"
    return out


def block_C():
    frac = (7.0 / 8.0) * (4.0 / 11.0) ** (4.0 / 3.0)
    act_up = ACT_DR6_NEFF[0] + 1.96 * ACT_DR6_NEFF[1] - NEFF_SM
    bbn_up = BBN_2024_DNEFF[0] + 1.96 * BBN_2024_DNEFF[1]
    return {"rho_dr_over_rho_gamma_per_unit_DeltaNeff": frac,
            "ACT_DR6_gaussian_95_upper_DeltaNeff": act_up,
            "ACT_DR6_rho_dr_over_rho_gamma_upper": act_up * frac,
            "BBN_2024_gaussian_95_upper_DeltaNeff": bbn_up,
            "BBN_2024_rho_dr_over_rho_gamma_upper": bbn_up * frac,
            "caveat": "Gaussian approximations built from quoted central values and 1-sigma errors; not the collaborations' own 95% limits"}


def main():
    res = {"script": Path(__file__).name, "numpy": np.__version__, "sympy": sp.__version__,
           "precision": "IEEE double; mpmath 30 digits for the Beta function; sympy exact for model facts",
           "A_knee_benchmark": block_A(), "B_braneworld_scale_setting": block_B(),
           "C_dark_radiation": block_C(), "literature_inputs": LIT_INPUTS,
           "data_spans_yr": DATA_SPANS_yr}
    A = res["A_knee_benchmark"]
    B = res["B_braneworld_scale_setting"]
    checks = {
        "Omega(f_k)=Omega_k": abs(A["shape"]["Omega_at_fk_over_Omega_k"] - 1) < 1e-14,
        "slopes_limit_+2_and_-3": abs(A["shape"]["slope_far_below_knee"] - 2) < 1e-6 and abs(A["shape"]["slope_far_above_knee"] + 3) < 1e-6,
        "slope_analytic_vs_fd": A["shape"]["analytic_vs_finite_difference_slope_max_abs_diff"] < 1e-6,
        "wrong_formula_rejected": A["shape"]["wrong_formula_control"]["naive_reading_rejected"],
        "calibration_NG15_amplitude_gives_Omega_k": A["amplitude_at_1_per_yr"]["calibration_control"]["matches_Omega_k_within_5pct"],
        "DeltaNeff_quad_vs_exact_1e-9": A["delta_neff"]["rel_diff_quad_vs_exact"] < 1e-9,
        "DeltaNeff_trap_vs_exact_1e-6": A["delta_neff"]["rel_diff_trap200k_vs_exact"] < 1e-6,
        "H_ell_two_routes_1e-9": B["registered_H_ell"]["rel_diff"] < 1e-9,
        "H2_expansion_O_delta2_within_1e-5": B["registered_H_ell"]["rel_diff_expansion2_vs_registered"] < 1e-5,
        "sigma_crit_matches_1_over_ell": B["exact_model_facts"]["sigma_crit/6 == 1/ell"],
        "frequency_map_textbook_coefficient": B["frequency_map_control"]["agrees_with_commonly_quoted_1.65e-7_within_2pct"],
        "M5_review_bound_reproduced": B["M5_control_vs_review_1e8GeV"],
    }
    res["checks"] = checks
    res["all_checks_pass"] = all(checks.values())
    OUT.write_text(json.dumps(res, indent=1, default=float))
    print(json.dumps(checks, indent=1))
    print("all_checks_pass:", res["all_checks_pass"])


if __name__ == "__main__":
    main()
