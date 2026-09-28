#!/usr/bin/env python3
"""FLOAT validation (guidance, NOT a proof) of OSCILLATION_THEOREM.md.
For each shell in shells.json: integrate the cone-regular solution psi ~ y^alpha for a grid of mu2 (lam = mu2 + 4) and
  Z     = number of zeros of psi (= of P = rho^2 psi) on (0, y_b)               [sign changes on the grid]
  F     = p P'/P at y_b = (psi' + 2H psi)/(phi'^2 rho^2 psi)
  G     = F - lam/beta,         mhat = m/psi_b,  m = B chi + 3 lam psi/(rho^2 phi')   (identity: m/(B phi' psi) = -3 rho^2 G)
  Theta = pi Z + arccot(F) - arccot(lam/beta)   (B > 0)
  N_thm(lam*) = [B<0] + Z + [G <= 0]   (Theorem 2/3)   versus   N_root(lam*) = number of sign changes of m-hat-normalised mismatch below lam*.
Also checks: dF/dlam = -int w P^2 / P_b^2 (Lemma 3), the identity int(pP'^2+qP^2) = lam (int wP^2 + P_b^2/beta) at the roots (Thm 1),
and the y^(2 nu) decay of the two cone boundary terms (Lemma 2 / item v)."""
import sys, math, json
import numpy as np
BG = "/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/background"
sys.path.insert(0, BG)
import hdblast_background as bg

def grid(y0, y_b, ratio=1.02, ymid=0.3, dy=2e-3):
    g = [y0]
    while g[-1]*ratio < min(ymid, y_b): g.append(g[-1]*ratio)
    n = max(1, int(math.ceil((y_b - g[-1])/dy))); g += list(np.linspace(g[-1], y_b, n + 1)[1:])
    return np.array(g)

def shoot(mu2, phi_h, y_b, y0=1e-3, record_at=None):
    """returns dict of arrays over mu2.  psi normalised to psi = y^alpha (1 + o(1))  (so C = 1)."""
    mu2 = np.asarray(mu2, float); M = mu2.size; lam = mu2 + 4
    alpha = 0.5 + np.sqrt(np.maximum(2.25 - mu2, 0.0))
    b = bg.series(phi_h, y0); s0 = b[2]
    Y = np.zeros((7, M)); Y[0] = b[0]; Y[1] = b[1]; Y[2] = b[2]
    Y[3] = y0**alpha; Y[4] = alpha*y0**(alpha - 1)
    two_nu = 2*alpha - 1
    with np.errstate(divide="ignore"):
        Y[5] = np.where(two_nu > 0, Y[3]**2*y0/(s0*s0*two_nu), np.inf)                       # int_0^y0 w P^2
        Y[6] = np.where(two_nu > 0, (alpha + 2)**2*Y[3]**2*y0/(s0*s0*two_nu), np.inf)        # int_0^y0 p P'^2 (+ q P^2 negligible)
    def rhs(Y):
        rho, phi, s, psi, dpsi = Y[:5]
        rp = np.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); H = rp/rho
        spp = bg.U1(phi) - 4*H*s; g = spp/s
        ddpsi = -2*(H - g)*dpsi - (-(4/3)*s*s - 4*H*g + (2 + mu2)/rho**2)*psi
        return np.array([rp, s, spp, dpsi, ddpsi, psi*psi/(s*s), (dpsi + 2*H*psi)**2*rho*rho/(s*s) + (2/3)*rho*rho*psi*psi])
    ys = grid(y0, y_b); Z = np.zeros(M, int); prev = np.sign(Y[3]); rec = {}
    for i in range(len(ys) - 1):
        h = ys[i+1] - ys[i]
        k1 = rhs(Y); k2 = rhs(Y + h/2*k1); k3 = rhs(Y + h/2*k2); k4 = rhs(Y + h*k3)
        Y = Y + h/6*(k1 + 2*k2 + 2*k3 + k4)
        sg = np.sign(Y[3]); Z += (sg*prev < 0); prev = np.where(sg != 0, sg, prev)
        if record_at is not None and any(abs(ys[i+1] - r) < 1e-12 for r in record_at): rec[float(ys[i+1])] = Y.copy()
    rho, phi, s, psi, dpsi, Iw, Idir = Y
    rp = np.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); H = rp/rho; spp = bg.U1(phi) - 4*H*s
    chi = -3*(dpsi + 2*H*psi)/s
    P = rho*rho*psi; pPp = (dpsi + 2*H*psi)/(s*s)
    return dict(lam=lam, Z=Z, P=P, pPp=pPp, F=pPp/P, psi=psi, chi=chi, Iw=Iw, Idir=Idir, rho=rho[0], s=s[0], phi=phi[0], spp=spp[0], ys=ys, rec=rec)

def arccot(x): return np.pi/2 - np.arctan(x)

def analyse(sh, t=1e-3):
    d = sh["d"]; phi_h = sh["phi_h"]; y_b = sh["y_b"]
    mu2 = np.concatenate([np.linspace(-60, -8.5, 104), np.linspace(-8.4, 2.24, 267), [2.2499, 2.25]])
    r = shoot(mu2, phi_h, y_b)
    sig2 = 2*bg.W2(r["phi"]) + t*d; B = r["spp"]/r["s"] + sig2/2; beta = r["rho"]**4*r["s"]**2*B
    lam = r["lam"]; G = r["F"] - lam/beta
    m = B*r["chi"] + 3*lam*r["psi"]/(r["rho"]**2*r["s"]); mhat = m/r["psi"]
    ident = np.max(np.abs(m/(B*r["s"]*r["psi"]) + 3*r["rho"]**2*G)/(np.abs(3*r["rho"]**2*G) + 1e-300))
    # roots of the mismatch D = beta pP' - lam P (smooth, no poles)
    Dm = beta*r["pPp"] - lam*r["P"]
    roots = []
    for i in range(len(mu2) - 1):
        if Dm[i]*Dm[i+1] < 0:
            a, b_ = mu2[i], mu2[i+1]
            for _ in range(5):
                g = np.linspace(a, b_, 17); rr = shoot(g, phi_h, y_b); v = beta*rr["pPp"] - (g + 4)*rr["P"]
                j = np.where(v[:-1]*v[1:] < 0)[0][0]; a, b_ = g[j], g[j+1]
            roots.append(0.5*(a + b_))
    roots = np.array(roots)
    N_root = np.array([(roots <= x).sum() for x in mu2])
    N_thm = (1 if B < 0 else 0)*np.ones(len(mu2), int) + r["Z"] + (G <= 0)
    if B < 0:   # for lam* <= 0 and B < 0 the count is [G >= 0] (Theorem 3(b))
        N_thm = np.where(lam <= 0, (G >= 0).astype(int), N_thm)
    mism = int((N_thm != N_root).sum())
    # monotonicity of F and of Theta
    finite = np.isfinite(r["Iw"])
    dF_num = np.gradient(r["F"], lam); dF_thm = -r["Iw"]/r["P"]**2
    sel = finite.copy(); sel[[0, -1, -2, -3]] = False; sel[103:105] = False
    dF_err = float(np.max(np.abs(dF_num[sel]/dF_thm[sel] - 1)))
    theta_b = np.pi*r["Z"] + arccot(r["F"])
    Theta = theta_b - arccot(lam/beta) if B > 0 else None
    out = dict(d=d, B=float(B), beta=float(beta), phi_prime_b=float(r["s"]), roots_mu2=[float(x) for x in roots],
               max_zeros_Z=int(r["Z"].max()), identity_m_vs_G_relerr=float(ident),
               count_mismatches_on_grid=mism, n_grid=len(mu2),
               F_decreasing=bool(np.all(np.diff(r["F"]) < 0)), theta_b_increasing=bool(np.all(np.diff(theta_b) > 0)),
               theta_b_range=[float(theta_b[0]), float(theta_b[-1])], dF_dlam_vs_minus_intwP2_over_Pb2_max_relerr=dF_err,
               threshold=dict(mu2=2.25, Z=int(r["Z"][-1]), F=float(r["F"][-1]), lam_over_beta=float(lam[-1]/beta), G=float(G[-1]),
                              mhat=float(mhat[-1]), mhat_times_sign_phi_prime=float(mhat[-1]*np.sign(r["s"])),
                              N_total_theorem=int((1 if B < 0 else 0) + r["Z"][-1] + (G[-1] < 0))),
               N_total_roots=len(roots))
    if Theta is not None:
        out["Theta_increasing"] = bool(np.all(np.diff(Theta) > 0)); out["Theta_over_pi_at_threshold"] = float(Theta[-1]/np.pi)
        out["N_total_from_Theta"] = int(max(0, math.ceil(Theta[-1]/np.pi)))
    # norm identity at the roots
    if len(roots):
        rr = shoot(roots, phi_h, y_b); lr = roots + 4
        brack = rr["Iw"] + rr["P"]**2/beta
        out["norm_identity"] = [dict(mu2=float(roots[k]), dirichlet=float(rr["Idir"][k]), lam_times_bracket=float(lr[k]*brack[k]),
                                     ratio=float(rr["Idir"][k]/(lr[k]*brack[k])), bracket_sign=float(np.sign(brack[k])),
                                     zeros=int(rr["Z"][k])) for k in range(len(roots))]
    return out

def cone_terms(sh):
    """boundary terms T1 = P pP' and T2 = p (P P_lam' - P' P_lam) at small y: both should scale as y^(2 nu)."""
    phi_h = sh["phi_h"]; res = []
    for mu2 in (-7.7, 0.0, 1.5):
        nu = math.sqrt(2.25 - mu2); e = 1e-5; pts = [0.01, 0.02, 0.04]
        ys = grid(1e-3, 0.05); pick = [ys[np.argmin(abs(ys - q))] for q in pts]
        r = shoot(np.array([mu2 - e, mu2, mu2 + e]), phi_h, 0.05, record_at=pick)
        rows = []
        for yk in sorted(r["rec"]):
            Y = r["rec"][yk]; rho, phi, s, psi, dpsi = Y[:5]
            rp = np.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); H = rp/rho
            P = rho*rho*psi; pPp = (dpsi + 2*H*psi)/(s*s)
            Pl = (P[2] - P[0])/(2*e); pPl = (pPp[2] - pPp[0])/(2*e)
            rows.append((yk, P[1]*pPp[1], P[1]*pPl - pPp[1]*Pl))
        sl1 = math.log(rows[-1][1]/rows[0][1])/math.log(rows[-1][0]/rows[0][0])
        sl2 = math.log(rows[-1][2]/rows[0][2])/math.log(rows[-1][0]/rows[0][0])
        a = bg.U1(phi_h)/5
        res.append(dict(mu2=mu2, two_nu=2*nu, slope_P_pPprime=sl1, slope_lambda_wronskian=sl2,
                        T1_over_prediction=[x[1]/((2.5 + nu)*x[0]**(2*nu)/a**2) for x in rows],
                        T2_over_prediction=[x[2]/(-x[0]**(2*nu)/(2*nu*a**2)) for x in rows]))
    return res

if __name__ == "__main__":
    shells = json.load(open("shells.json")); out = []
    for sh in shells:
        o = analyse(sh); out.append(o)
        print(json.dumps(o), flush=True)
    ct = cone_terms(shells[0]); print(json.dumps(ct))
    json.dump(dict(shells=out, cone_terms_d0=ct), open("osc_validate_output.json", "w"), indent=1)
