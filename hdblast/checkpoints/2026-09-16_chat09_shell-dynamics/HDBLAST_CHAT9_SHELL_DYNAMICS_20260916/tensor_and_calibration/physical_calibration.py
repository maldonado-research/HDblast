#!/usr/bin/env python3
"""HDBLAST Chat 9 - physical calibration of the registered (dimensionless) shell, t = 1e-3.
Units restored as  S = M_5^3 [ int d^5x sqrt(-g) (R/2 - (d phi)^2/2 - U/L_0^2) - int d^4x sqrt(-g_4) sigma_t/L_0 ],  phi dimensionless.
  M_4^2 = (2 I_kept) M_5^3 L_0,   H = Hhat/L_0,   k_- = 5/(9 L_0),  k_+ = 1/(9 L_0),  sigma = sigma_hat M_5^3/L_0,  delta sigma = t(1+c phi_b) M_5^3/L_0.
Given M_Pl (reduced) and H:  L_0 = Hhat/H,  M_5^3 = M_Pl^2 H/(2 I_kept Hhat)."""
import json, math, os
HERE = os.path.dirname(os.path.abspath(__file__))
T = json.load(open(os.path.join(HERE, "TENSOR_SPECTRUM.json")))
S = json.load(open(os.path.join(HERE, "..", "background", "REGISTERED_SHELL_FLOAT_SOLUTION.json")))
MPl = 2.435e18; hbar = 6.582119569e-25; hbarc = 1.973269804e-16     # GeV, GeV s, GeV m
Hhat = T["H"]; m4 = T["M4sq"]; t = S["t"]; c = S["c"]; phib = S["phi_b"]
sig_hat = 2*(1 - phib + phib**3/3) + t*(1 + c*phib); dsig_hat = t*(1 + c*phib)
s_growth = (-3 + math.sqrt(9 + 4*7.7178716))/2
LMW = lambda x: 1/(math.sqrt(1 + x*x) - x*x*math.asinh(1/x))
print("Hhat = %.8f, 2 I_kept = %.8f, sigma_hat = %.6f, delta sigma_hat = %.6e, 3 Hhat^2 (2 I_kept) = %.6e, s = %.5f (growth time %.4f/H)" % (Hhat, m4, sig_hat, dsig_hat, 3*Hhat**2*m4, s_growth, 1/s_growth))
print("H/k_- = %.5f, H/k_+ = %.5f, H L_0 = %.5f ; LMW F^2(H/k_-) = %.6f, F^2(H/k_+) = %.6f ; shell zero-mode F^2 = I_+/I_kept = %.6f" % (Hhat*9/5, Hhat*9, Hhat, LMW(Hhat*9/5), LMW(Hhat*9), T["ratio_Ip_over_Ikept"]))
rows = []
def case(label, H):
    L0 = Hhat/H; M5 = (MPl**2*H/(m4*Hhat))**(1/3)
    r = dict(label=label, H_GeV=H, L0_invGeV=L0, L0_m=L0*hbarc, M5_GeV=M5, M5L0=M5*L0, k_minus_GeV=5/(9*L0), k_plus_GeV=1/(9*L0),
             l_minus_m=1.8*L0*hbarc, l_plus_m=9*L0*hbarc, tension_GeV4=sig_hat*M5**3/L0, tension_quarter_GeV=(sig_hat*M5**3/L0)**0.25,
             detuning_GeV4=dsig_hat*M5**3/L0, detuning_quarter_GeV=(dsig_hat*M5**3/L0)**0.25, check_3H2MPl2_GeV4=3*H*H*MPl**2,
             KK_gap_GeV=1.5*H, tachyon_mass_GeV=math.sqrt(7.7178716)*H, growth_time_s=hbar/(s_growth*H), hubble_time_s=hbar/H,
             P_T=2*H*H/(math.pi**2*MPl**2)*T["ratio_Ip_over_Ikept"]*0 + 2*H*H/(math.pi**2*MPl**2))
    rows.append(r); return r
for lab, H in (("H=1e13 GeV", 1e13), ("H=1e10 GeV", 1e10), ("H=1e5 GeV", 1e5)): case(lab, H)
for lmax_um in (50.0, 30.0):
    L0 = lmax_um*1e-6/1.8/hbarc                                         # l_- = 1.8 L_0 = lmax
    case("sub-mm limit l_- = %.0f micron" % lmax_um, Hhat/L0)
hdr = ("H_GeV", "L0_invGeV", "L0_m", "M5_GeV", "M5L0", "k_minus_GeV", "k_plus_GeV", "l_minus_m", "l_plus_m", "tension_quarter_GeV", "detuning_quarter_GeV", "detuning_GeV4", "check_3H2MPl2_GeV4", "KK_gap_GeV", "growth_time_s")
for r in rows:
    print("\n[%s]" % r["label"]); print("   " + "  ".join("%s=%.4g" % (k, r[k]) for k in hdr))
json.dump(dict(Hhat=Hhat, M4sq_hat=m4, sigma_hat=sig_hat, dsigma_hat=dsig_hat, s_growth=s_growth, LMW_F2_kminus=LMW(Hhat*9/5), F2_shell=T["ratio_Ip_over_Ikept"], cases=rows),
          open(os.path.join(HERE, "PHYSICAL_CALIBRATION.json"), "w"), indent=1)
