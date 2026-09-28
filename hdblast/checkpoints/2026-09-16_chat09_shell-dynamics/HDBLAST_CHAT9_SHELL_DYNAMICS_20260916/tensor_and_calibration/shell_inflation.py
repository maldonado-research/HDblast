#!/usr/bin/env python3
"""HDBLAST Chat 9 - slow-roll predictions of the shell-modulus EFT (orchestrator result R3) with detuning
    sigma_t - 2W = t (1 + c phi + d phi^2/2),   V_E = t num/f^2,  canonical field Theta (dTheta/dy_b = +sqrt(Z_E,y); Theta>0 = AdS_- throat side).
Einstein-frame units M_Pl,E^2 = M_5^3 L_0 = 1.  Exact homogeneous dynamics in e-folds:  Th'' = -(3 - Th'^2/2)(Th' + G),  G = dlnV/dTheta.
Hubble-flow: e1 = Th'^2/2, e2 = dln e1/dN, e3 = dln e2/dN.  n_s-1 = -2e1-e2-2e1^2-(2C+3)e1e2-C e2e3, r = 16 e1 (1 + C e2), C=-0.7296,
alpha_s = -2e1e2 - e2e3,  A_s = H_E^2/(8 pi^2 e1) [x t/(M_5 L_0)^3 when units are restored]."""
import json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
D = np.load(os.path.join(HERE, "..", "effective_theory_orchestrator", "eft_landscape.npz"))
m = (D["y"] <= 9.0) & (D["y"] >= -16.0)
Y, PHI, I, F, K, TH = [D[k][m][::-1] for k in ("y", "phi", "I", "f", "ZEy", "Theta")]      # ascending y (and ascending Theta)
Wf = lambda p: 1 - p + p**3/3
Ip = 1.0357712571566782; c_star = 2/Ip - 4/3
FY = 2*(-1 + (2*Wf(PHI)/3)*I); PY = -(1 - PHI**2); SK = np.sqrt(K)
FYY = 2*((2*(PHI**2 - 1)*PY/3)*I + (2*Wf(PHI)/3)*FY/2)
NG = 400001; TG = np.linspace(TH[0], TH[-1], NG); DT = TG[1] - TG[0]
CE = -0.7296

def tables(c, d, e=0.0):
    num = 1 + c*PHI + d*PHI**2/2 + e*PHI**3/6
    lnV = np.log(num) - 2*np.log(F)
    G = ((c + d*PHI + e*PHI**2/2)*PY/num - 2*FY/F)/SK
    GT = np.gradient(G, Y)/SK
    return np.interp(TG, TH, lnV), np.interp(TG, TH, G), np.interp(TG, TH, GT)

def lin(tab, th):
    x = np.clip((th - TG[0])/DT, 0, NG - 1.000001); i = x.astype(int); w = x - i
    return tab[i]*(1 - w) + tab[i+1]*w

def evolve(c, d, th0, e=0.0, dN=5e-3, Nmax=4000.0):
    lnV, G, GT = tables(c, d, e)
    th = float(th0); p = -float(lin(G, np.array([th]))[0])
    rhs = lambda th, p: (p, -(3 - p*p/2)*(p + float(lin(G, np.array([th]))[0])))
    rec = []; N = 0.0; ended = False
    while N < Nmax:
        rec.append((N, th, p))
        if p*p/2 >= 1: ended = True; break
        if th <= TG[2] or th >= TG[-3]: break
        k1 = rhs(th, p); k2 = rhs(th + dN/2*k1[0], p + dN/2*k1[1]); k3 = rhs(th + dN/2*k2[0], p + dN/2*k2[1]); k4 = rhs(th + dN*k3[0], p + dN*k3[1])
        th += dN/6*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0]); p += dN/6*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1]); N += dN
    rec = np.array(rec); Ns, th, p = rec.T
    g = lin(G, th); gt = lin(GT, th); A = 3 - p*p/2; B = p + g
    p1 = -A*B; p2 = p*p1*B - A*(p1 + gt*p)
    e1 = p*p/2; e2 = 2*p1/p; e3 = (2*p2/p - 2*p1*p1/(p*p))/e2
    v = np.exp(lin(lnV, th))
    return dict(N=Ns, th=th, e1=e1, e2=e2, e3=e3, v=v, ended=ended)

def observables(tr, Nstar, iend=None):
    iend = len(tr["N"]) - 1 if iend is None else iend
    Nb = tr["N"][iend] - tr["N"]; i = int(np.argmin(np.abs(Nb - Nstar)))
    if tr["N"][iend] < Nstar + 5: return None
    e1, e2, e3, v, th = [tr[k][i] for k in ("e1", "e2", "e3", "v", "th")]
    ns = 1 - 2*e1 - e2 - 2*e1*e1 - (2*CE + 3)*e1*e2 - CE*e2*e3
    r = 16*e1*(1 + CE*e2); al = -2*e1*e2 - e2*e3
    HE2 = v/(3 - e1)                                            # x t / L0^2, Einstein frame
    t_over_N5cubed = 2.1e-9*8*math.pi**2*e1/HE2                 # t/(M_5 L_0)^3 from A_s = 2.1e-9
    phi = float(np.interp(th, TH, PHI)); f = float(np.interp(th, TH, F))
    return dict(Nstar=Nstar, Theta=float(th), phi_b=phi, f=f, e1=float(e1), e2=float(e2), e3=float(e3), n_s=float(ns), r=float(r), alpha_s=float(al),
                t_over_M5L0_cubed=float(t_over_N5cubed), H_over_MPl=float(math.sqrt(2.1e-9*8*math.pi**2*e1)))

def calibrate(o, t):
    MPl = 2.435e18; N5 = (t/o["t_over_M5L0_cubed"])**(1/3)      # M_5 L_0
    M5 = MPl/math.sqrt(o["f"]*N5)                               # M_Pl^2 = f M_5^3 L_0 (f at horizon exit; O(1) ambiguity, see text)
    L0 = N5/M5
    return dict(t=t, M5L0=N5, M5_GeV=M5, L0_invGeV=L0, L0_m=L0*1.973269804e-16, H_GeV=o["H_over_MPl"]*MPl, k_minus_GeV=5/(9*L0), H_over_kminus=o["H_over_MPl"]*MPl*9*L0/5)

if __name__ == "__main__":
    out = {}
    lnV, G, GT = tables(c_star, 0.0)
    i0 = int(np.argmin(np.abs(TG))); 
    # hilltop data vs d
    print("=== A. registered family c = c_star, hilltop at phi_b = 0 ===")
    for d in (0.0, 1.0, 1.105924, 1.1134966):
        lnV, G, GT = tables(c_star, d)
        GTT = np.gradient(GT, TG)
        print(" d=%.7f: eta_V(0)=V''/V=%+.6f (mu2=%+.5f), gamma=V'''/V=%+.4f ; max eps_V on phi>0 side=%.4f ; zeros of V' on throat side (Theta>1e-3): %s" % (
            d, GT[i0] + G[i0]**2, 3*(GT[i0] + G[i0]**2), GTT[i0] + 3*G[i0]*GT[i0], (G[:i0]**2/2).max(),
            [round(float(TG[j]), 5) for j in np.where(np.sign(G[i0+20:-1]) != np.sign(G[i0+21:]))[0] + i0 + 20]))
    out["A"] = []
    for d in (1.08, 1.10, 1.105924, 1.11, 1.1134):
        lnV, G, GT = tables(c_star, d); eta = GT[i0]
        tr = evolve(c_star, d, -2e-4 if eta < 0 else -2e-4)     # roll toward phi -> +1 (Theta < 0)
        imax = int(np.argmax(tr["e1"]))
        print(" d=%.6f eta0=%+.5f: roll to phi->+1: inflation ends (eps_H=1)? %s ; max eps_H=%.4f at Theta=%.4f after N=%.1f ; total N integrated %.1f" % (d, eta, tr["ended"], tr["e1"][imax], tr["th"][imax], tr["N"][imax], tr["N"][-1]))
        for Ns in (50, 55, 60):
            o = observables(tr, Ns, iend=imax)
            if o: print("     N*=%d before max-eps point: n_s=%.5f r=%.3e alpha_s=%+.2e Theta*=%+.5f | 1-4/N=%.5f, -2|eta|coth(|eta|N/2) formula: %.5f" % (Ns, o["n_s"], o["r"], o["alpha_s"], o["Theta"], 1 - 4/Ns, 1 - 2*abs(eta)/math.tanh(abs(eta)*Ns/2))); out["A"].append(dict(d=d, **o))
    json.dump(out, open(os.path.join(HERE, "SHELL_INFLATION_A.json"), "w"), indent=1)
