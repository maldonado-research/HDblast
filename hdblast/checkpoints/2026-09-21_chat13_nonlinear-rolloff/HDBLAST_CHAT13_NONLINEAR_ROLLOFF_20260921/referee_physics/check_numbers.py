#!/usr/bin/env python3
"""Referee (physics) - reproduces every number quoted in REFEREE_PHYSICS.md from ../runs/*.npz and referee_physics/runs/*.npz.
Run from the HDBLAST_CHAT13 folder with the system python3 (numpy only)."""
import numpy as np, math, json, os, glob
np.seterr(all="ignore")
HERE = os.path.dirname(os.path.abspath(__file__)); TOP = os.path.dirname(HERE)
Ip = 1.0357712571566784; c = 2/Ip - 4/3
W = lambda p: 1 - p + p**3/3; W1 = lambda p: p*p - 1; U = lambda p: 0.5*W1(p)**2 - (2/3)*W(p)**2
P = json.load(open(os.path.join(TOP, "runs/linear_predictions.json")))
trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

print("== Section 2: de Sitter brane in AdS_+ (kappa5=1, Z2):  H^2 = sigma^2/36 - 1/l+^2,  1/l+ = %.4f" % math.sqrt(-U(1)/6))
for t, rho_b in [(0.1, P["0.1"]["rho_b"]), (0.03, P["0.03"]["rho_b"]), (1e-3, 78.82817923791542)]:
    H0 = 1/rho_b; sig = 2*W(1) + t*(1 + c)
    print("  t=%g H0=%.6f H0^2/t=%.4f exact %.4f | O(t) only %.4f" % (t, H0, H0**2/t, math.sqrt(sig**2/36 - 1/81)/H0, math.sqrt((1 + c)*t/27)/H0))
    for pb in [0.9908, 0.995, 1.0]:
        s = 2*W(pb) + t*(1 + c*pb); s1 = 2*W1(pb) + t*c
        print("     quasi-static (nn Gauss, Hdot=phidot=0) at phi_b=%.4f: sigma'=%+.4f H/H0=%.4f" % (pb, s1, math.sqrt(s*s/36 - s1*s1/48 + U(pb)/6)/H0))

def load(tag, base=os.path.join(TOP, "runs")):
    r = np.load(os.path.join(base, tag + "_timeseries.npz"))["rec"]
    return dict(zip("t phib HJ bb pa dphidtau M mf mb tau pmin pmax".split(), r.T))

print("\n== Section 1/2: plus-branch windows")
for tag, rho_b in [("t01_plus_v2", P["0.1"]["rho_b"]), ("t01_plus_coarse", P["0.1"]["rho_b"]), ("t003_plus", P["0.03"]["rho_b"])]:
    d = load(tag); t, phib, HJ, bb, tau = d["t"], d["phib"], d["HJ"], d["bb"], d["tau"]
    j = np.argmax(phib > 0.98); late = t > t[j] + 0.5; s = np.polyfit(t[late], bb[late], 1)[0]
    print("  %-16s end: t=%.2f H0tau=%.4f phi_b=%.5f HJ/H0=%.4f | after phi_b>0.98: dH0tau=%.3f efolds=%.3f | db/dt=%.3f -> H0 tau_inf=%.3f"
          % (tag, t[-1], tau[-1], phib[-1], HJ[-1], tau[-1] - tau[j], trapz(HJ[j:], tau[j:]), s, tau[-1] + math.exp(bb[-1])/(-s)))
    if tag == "t01_plus_v2":
        N = np.concatenate([[0], np.cumsum(0.5*(HJ[1:] + HJ[:-1])*np.diff(tau))])
        for (t1, t2) in [(8.0, 11.0), (8.5, 12.0)]:
            m = (t >= t1) & (t <= t2)
            for n in [4, 6]:
                X = np.exp(-n*(N[m] - N[m][0])); A = np.vstack([np.ones_like(X), X]).T; sol = np.linalg.lstsq(A, HJ[m]**2, rcond=None)[0]
                print("     fit H^2 = Hinf^2 + C a^-%d on t in [%g,%g]: Hinf/H0=%.4f C/H0^2=%.4f rms=%.1e" % (n, t1, t2, math.sqrt(sol[0]), sol[1], np.sqrt(np.mean((A@sol - HJ[m]**2)**2))))

print("\n== Section 3: throat branch")
print("  U(-1.96)=%.3f W(-1.96)=%.3f sigma_t(-1.96,t=0.1)=%.3f ; V_E=0 at -1/c=%.4f ; W=0 at -2.1038 ; U(-6.84)=%.0f sigma_t=%.0f" % (U(-1.96), W(-1.96), 2*W(-1.96) + 0.1*(1 + c*(-1.96)), -1/c, U(-6.84), 2*W(-6.84) + 0.1*(1 + c*(-6.84))))
ph = np.linspace(-2.2, -1.0, 12001); g = 4*W(ph)*(1 + c*ph) - 3*c*W1(ph); print("  O(t) quasi-static H^2 sign change at phi_b=%.3f" % ph[np.argmax(g > 0)])
rho_b = P["0.1"]["rho_b"]; Rst = (10/3)*U(-1)
for tag in ["t01_minus_ext", "t01_minus_coarse"]:
    d = load(tag); t, phib, HJ, tau, pd = d["t"], d["phib"], d["HJ"], d["tau"], d["dphidtau"]/rho_b
    j0 = np.argmax(HJ < 0); i1 = np.argmax(phib < -1); i2 = np.argmax(phib < -1/c)
    print("  %-16s cross phi_b=-1: H0tau=%.3f HJ/H0=%.3f | V_E=0 pt: HJ/H0=%.3f | turnaround H0tau=%.4f phi_b=%.4f | end H0tau=%.4f HJ/H0=%.1f phi_b=%.2f dln a=%.2f maxM=%.3f"
          % (tag, tau[i1], HJ[i1], HJ[i2], tau[j0], phib[j0], tau[-1], HJ[-1], phib[-1], trapz(HJ[j0:], tau[j0:]), np.nanmax(d["M"])))
    for p in [-5, -6]:
        i = np.argmax(phib < p); print("     phi_b=%g: H0tau=%.4f HJ/H0=%.2f dln a since turnaround %.3f" % (p, tau[i], HJ[i], trapz(HJ[j0:i], tau[j0:i])))
    for thr in [-3, -10, -30]:
        k = HJ < thr; A_ = np.polyfit(tau[k], -1/HJ[k], 1); print("     fit HJ=-alpha/(tau*-tau), HJ<%g: alpha=%.3f H0tau*=%.4f" % (thr, -1/A_[0], -A_[1]/A_[0]))
    k = HJ < -3; al = -np.gradient(-1/HJ[k], tau[k]); print("     local alpha min/max: %.2f %.2f" % (al.min(), al.max()))
    lna = np.concatenate([[0], np.cumsum(0.5*(HJ[1:] + HJ[:-1])*np.diff(tau))]); kk = HJ < -1
    print("     dphi_b/dln a in contraction: %.2f" % np.polyfit(lna[kk], phib[kk], 1)[0])
    i = np.argmin(abs(tau - 6.415)); print("     at H0tau=6.415: HJ/H0=%.1f" % HJ[i])
    for tt in [0, 6.0, 7.2, 9.0]:
        i = np.argmin(abs(t - tt)); p = phib[i]; s1 = 2*W1(p) + 0.1*c; R = (s1/2)**2 - pd[i]**2 + (10/3)*U(p)
        print("     t=%.1f phi_b=%+.2f R_5D=%.3g |R/R_AdS-|=%.3g" % (tt, p, R, abs(R/Rst)))

print("\n== Section 7: referee check runs")
for fn in sorted(glob.glob(os.path.join(HERE, "runs", "*_timeseries.npz"))):
    tag = os.path.basename(fn).replace("_timeseries.npz", ""); d = load(tag, os.path.join(HERE, "runs"))
    t, phib, HJ, tau, bb = d["t"], d["phib"], d["HJ"], d["tau"], d["bb"]
    print("  %-14s end t=%.2f H0tau=%.4f phi_b=%.4f HJ/H0=%.4f maxM=%.3f" % (tag, t[-1], tau[-1], phib[-1], HJ[-1], np.nanmax(d["M"])))
    for tt in [8, 9, 10, 11, 12]:
        i = np.argmin(abs(t - tt))
        if abs(t[i] - tt) < 0.02: print("     t=%d H0tau=%.4f phi_b=%.5f HJ/H0=%.5f" % (tt, tau[i], phib[i], HJ[i]))
    if (HJ < 0).any():
        j0 = np.argmax(HJ < 0); print("     turnaround H0tau=%.4f phi_b=%.4f ; end dln a=%.2f" % (tau[j0], phib[j0], trapz(HJ[j0:], tau[j0:])))
