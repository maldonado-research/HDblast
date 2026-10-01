#!/usr/bin/env python3
"""B4 reduced model E (REGISTRATION.md section 1): Chat 9 two-derivative 4D EFT roll with friction (phase I),
matching at phi = phi_m, exact late Friedmann sector with Lambda_res(d) from the exact static branch (phase II).

Phase I (Einstein frame, kappa5^2 = 1 model units):  f = 2I(phi), W_phi I' - (2/3) W I = -1, I(0) = I+ ;
  Z = 2W(1 - 2WI/3)/W_phi^2,  Z_E = Z/f + 1.5 (f'/f)^2,  V = delta(1 + c phi + d phi^2/2),  V_E = V/f^2,
  3H_E^2 = Z_E pi^2/2 + V_E + X e^{-4N_E},   Z_E(pi' + 3H_E pi) + Z_E' pi^2/2 + V_E' = -Gamma_E pi,  Gamma_E = Y f^{-3/2},
  X = rho_E a_E^4,  X' = Gamma_E pi^2 e^{4N_E};  brane frame: dtau = dt_E/sqrt f, ln a = N_E - ln f/2, R = f^2 X e^{-4N_E},
  H = sqrt(f)(H_E - f' pi/(2f)),  v = dphi/dtau = sqrt(f) pi.
Matching variants:
  'registered' : brane-frame H, R continuous; Weyl_m = H_m^2 - Lambda_res - rad(R_m) (identity, scalar frozen at phi_s(d)).
  'einstein'   : (post-registration variant, dated note) scalar energy above the EFT endpoint, in the Einstein frame,
                 is given to the Weyl fluid: r = [Z_E pi^2/2 + V_E(phi_m) - V_E(1)]/rho_E, then the late brane state is
                 R a^4 (continuous), Weyl a^4 = r * rad a^4, and H from the identity (H jumps at the match).
Phase II: N' = h, h' = E'(N)/2, h^2 = E(N) = lam + K_r e^{-4N} + K_rr e^{-8N} + K_w e^{-4N}  (exact for v = 0).
All dimensionless quantities in units of each model's own H0 (EFT H0 for phase I, 5D H0 for Lambda_res)."""
import json, math, sys, hashlib, argparse
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
HERE = Path(__file__).resolve().parent
C = 2/1.0357712571566784 - 4/3
I_PLUS = 1.0357712571566784
DELTA = 0.1
DSTAR = -3.106933495673783
ST = json.loads((HERE/'B4_STATIC_LAMBDA.json').read_text())
FIT = ST['fits']['deg3']

def W(p): return 1 - p + p**3/3
def W1(p): return p*p - 1
def W2(p): return 2*p
def U(p): return 0.5*(p*p - 1)**2 - (2/3)*W(p)**2
def U1(p): return 2*p*(p*p - 1) - (4/3)*W(p)*(p*p - 1)

class EFT:
    def __init__(self, rtol=1e-12, pmax=0.999):
        rhs = lambda p, I: [(-1 + (2/3)*W(p)*I[0])/W1(p)]
        self.lo = solve_ivp(rhs, (0.0, -0.5), [I_PLUS], method='DOP853', rtol=rtol, atol=1e-14, dense_output=True)
        self.hi = solve_ivp(rhs, (0.0, pmax), [I_PLUS], method='DOP853', rtol=rtol, atol=1e-14, dense_output=True)
    def I(self, p): return float(self.hi.sol(p)[0]) if p >= 0 else float(self.lo.sol(p)[0])
    def funcs(self, p):
        I = self.I(p); w, w1, w2 = W(p), W1(p), W2(p)
        I1 = (-1 + 2*w*I/3)/w1
        I2 = ((2*w1*I/3 + 2*w*I1/3)*w1 - (-1 + 2*w*I/3)*w2)/w1**2
        f, f1, f2 = 2*I, 2*I1, 2*I2
        Z = -w*f1/w1
        # Z' = -(w1 f1 + w f2)/w1 + w f1 w2/w1^2
        Z1 = -(w1*f1 + w*f2)/w1 + w*f1*w2/w1**2
        ZE = Z/f + 1.5*(f1/f)**2
        ZE1 = Z1/f - Z*f1/f**2 + 3*(f1/f)*(f2/f - (f1/f)**2)
        return f, f1, f2, Z, ZE, ZE1

def tension(d, c=C, delta=DELTA):
    s = lambda p: 2*W(p) + delta*(1 + c*p + d*p*p/2)
    s1 = lambda p: 2*W1(p) + delta*(c + d*p)
    s2 = lambda p: 4*p + delta*d
    V = lambda p: delta*(1 + c*p + d*p*p/2)
    V1 = lambda p: delta*(c + d*p)
    return s, s1, s2, V, V1

def VE_and_slope(e, V, V1, p):
    f, f1, *_ = e.funcs(p)
    return V(p)/f**2, (V1(p)*f - 2*V(p)*f1)/f**3

def hilltop(e, d, c):
    """EFT stationary point of V_E near phi = 0 for tension (c, d)."""
    _, _, _, V, V1 = tension(d, c)
    g = lambda p: VE_and_slope(e, V, V1, p)[1]
    return brentq(g, -0.05, 0.05, xtol=1e-15)

def growth_rate(e, d=0.0, c=C, h=1e-4):
    """linear growth rate / H at the hilltop: lambda^2 + 3 lambda + mu^2 = 0, mu^2 = V_E''/(Z_E H_E^2) (canonical, at rest)."""
    _, _, _, V, V1 = tension(d, c)
    p0 = hilltop(e, d, c)
    VEpp = (VE_and_slope(e, V, V1, p0 + h)[1] - VE_and_slope(e, V, V1, p0 - h)[1])/(2*h)
    VE0 = VE_and_slope(e, V, V1, p0)[0]; ZE = e.funcs(p0)[4]
    mu2 = VEpp/(ZE*VE0/3)
    return dict(phi_h=p0, mu2=mu2, rate_over_H=(-3 + math.sqrt(9 - 4*mu2))/2)

def run(e, Y, d, dc=1e-2, phi_m=0.95, rtol=1e-10, variant='registered', wrong=None, tmax_E=400.0):
    s, s1, s2, V, V1 = tension(d)
    p_seed = hilltop(e, d, C + dc)                     # rest at the seed's stationary point, evolved with c
    f0 = e.funcs(p_seed)[0]
    H0 = math.sqrt(V(hilltop(e, d, C))/(3*e.funcs(hilltop(e, d, C))[0]))     # EFT initial-shell Hubble rate (brane frame)
    def HE_of(y):
        p, pi, N, X, tau = y
        f, f1, f2, Z, ZE, ZE1 = e.funcs(p)
        return math.sqrt(max(0.5*ZE*pi*pi + V(p)/f**2 + X*math.exp(-4*N), 0.0))
    def rhs(t, y):
        p, pi, N, X, tau = y
        f, f1, f2, Z, ZE, ZE1 = e.funcs(p)
        HE = math.sqrt(max(0.5*ZE*pi*pi + V(p)/f**2 + X*math.exp(-4*N), 0.0)/3)
        G = Y if wrong == 'gamma' else Y*f**-1.5
        VE1 = (V1(p)*f - 2*V(p)*f1)/f**3
        dpi = (-3*HE*ZE*pi - 0.5*ZE1*pi*pi - VE1 - G*pi)/ZE
        return [pi, dpi, HE, G*pi*pi*math.exp(4*N), 1/math.sqrt(f)]
    ev = lambda t, y: y[0] - phi_m
    ev.terminal = True; ev.direction = 1
    y0 = [p_seed, 0.0, 0.5*math.log(f0), 0.0, 0.0]
    sol = solve_ivp(rhs, (0, tmax_E), y0, method='DOP853', rtol=rtol, atol=1e-13, events=ev, dense_output=True, max_step=0.05)
    if sol.status != 1: return dict(ok=False, msg='phi_m not reached', t_end=float(sol.t[-1]), phi_end=float(sol.y[0, -1]))
    # brane-frame records of phase I
    ts = np.concatenate([np.linspace(0, sol.t[-1], 4000)])
    rec = []
    for t in ts:
        p, pi, N, X, tau = sol.sol(t)
        f, f1, f2, Z, ZE, ZE1 = e.funcs(p)
        HE = math.sqrt(max(0.5*ZE*pi*pi + V(p)/f**2 + X*math.exp(-4*N), 0.0)/3)
        H = math.sqrt(f)*(HE - f1*pi/(2*f)) if wrong != 'nofdot' else math.sqrt(f)*HE
        v = math.sqrt(f)*pi; R = f*f*X*math.exp(-4*N); lna = N - 0.5*math.log(f)
        # identity-defined Weyl along the EFT trajectory (local shell terms, physical units)
        loc = (s(p) + R)**2/36 + v*v/12 - (s1(p) + Y*v)**2/48 + U(p)/6
        rad = s(p)*R/18 + R*R/36
        rec.append([tau*H0, H/H0, p, v/H0, lna, R, (H*H - loc)/H0**2, rad/H0**2, (s(p)**2/36 - s1(p)**2/48 + U(p)/6)/H0**2])
    rec = np.array(rec)
    # phase II
    lam = float(np.polyval(FIT['coef_Lambda'], d))           # Lambda_res/H0_5D^2
    sig_s = float(np.polyval(FIT['coef_sigma'], d)); phi_s = float(np.polyval(FIT['coef_phi'], d))
    p, pi, N, X, tau = sol.y[:, -1]
    f, f1, f2, Z, ZE, ZE1 = e.funcs(p)
    HE = math.sqrt((0.5*ZE*pi*pi + V(p)/f**2 + X*math.exp(-4*N))/3)
    Hm = math.sqrt(f)*(HE - f1*pi/(2*f)) if wrong != 'nofdot' else math.sqrt(f)*HE
    Rm = f*f*X*math.exp(-4*N); lnam = N - 0.5*math.log(f); hm = Hm/H0
    radm = (sig_s*Rm/18 + Rm*Rm/36)/H0**2
    if variant == 'registered':
        wm = hm*hm - lam - radm
    elif variant == 'radonly':
        wm = 0.0                       # post-registration variant: no Weyl fluid after the match (H from the identity, jumps)
    else:
        _, _, _, VV, VV1 = tension(d)
        VE1end = VV(1.0)/81.0          # EFT endpoint vacuum V(1)/f(1)^2, f(1) = 2 I(1) = 9
        rW = (0.5*ZE*pi*pi + V(p)/f**2 - VE1end)/(X*math.exp(-4*N))
        wm = rW*radm
    E = lambda n: lam + (sig_s*Rm/18/H0**2)*math.exp(-4*(n - lnam)) + (Rm*Rm/36/H0**2)*math.exp(-8*(n - lnam)) + wm*math.exp(-4*(n - lnam))
    dE = lambda n: -4*(sig_s*Rm/18/H0**2 + wm)*math.exp(-4*(n - lnam)) - 8*(Rm*Rm/36/H0**2)*math.exp(-8*(n - lnam))
    h_start = math.sqrt(E(lnam)) if E(lnam) > 0 else float('nan')
    t0 = tau*H0
    if not np.isfinite(h_start): return dict(ok=False, msg='E<0 at match')
    ev2 = lambda t, y: y[1] + 0.5
    ev2.terminal = True
    s2_ = solve_ivp(lambda t, y: [y[1], 0.5*dE(y[0])], (t0, t0 + 400.0), [lnam, h_start], method='DOP853', rtol=rtol, atol=1e-14,
                    events=ev2, dense_output=True, max_step=0.05)
    tt = np.linspace(t0, s2_.t[-1], 6000)
    N2, h2 = s2_.sol(tt)
    R2 = Rm*np.exp(-4*(N2 - lnam)); rad2 = (sig_s*R2/18 + R2*R2/36)/H0**2; W2 = wm*np.exp(-4*(N2 - lnam))
    rec2 = np.column_stack([tt, h2, np.full_like(tt, phi_s), 0*tt, N2, R2, W2, rad2, np.full_like(tt, lam)])
    allr = np.vstack([rec, rec2[1:]])
    cols = ['H0tau', 'H_over_H0', 'phi_b', 'v_over_H0', 'ln_a', 'R', 'Wy', 'rad', 'vac']
    D = {k: allr[:, i] for i, k in enumerate(cols)}
    return dict(ok=True, D=D, match=dict(H0tau=t0, h=hm, h_after=h_start, R=Rm, lna=lnam, Weyl=wm, rad=radm, r=wm/radm, lam=lam,
                                       phi_s=phi_s, sigma_s=sig_s, H0_EFT=H0, Omega_r=radm/hm**2, Omega_W=wm/hm**2))

# ---------- A1 analysis rule (copied from a1_analyze.py: plateau_index and recollapse definitions)
def plateau_index(lna, H, Wa4, Ra4, dlna=0.5, tol=0.05, need_W=True):
    n = len(lna); j0 = 0
    for p in range(n):
        if lna[p] - lna[0] < dlna: continue
        while j0 < p and lna[j0] < lna[p] - dlna: j0 += 1
        if lna[p] < lna[j0] + dlna - 0.02: continue
        win = slice(j0, p + 1)
        if np.any(H[win] <= 0): continue
        w = Wa4[win]; okW = ((w.max() - w.min()) <= tol*abs(Wa4[p]) if Wa4[p] != 0 else False) if need_W else True
        rr = Ra4[win]; okR = Ra4[p] > 0 and (rr.max() - rr.min()) <= tol*Ra4[p]
        if okW and okR: return p, j0
    return None, None

def analyse(res):
    D = res['D']; H = D['H_over_H0']; tau = D['H0tau']; lna = D['ln_a']
    Wa4 = D['Wy']*np.exp(4*lna); Ra4 = D['R']*np.exp(4*lna)
    with np.errstate(divide='ignore', invalid='ignore'):
        Om_r = D['rad']/H**2; r = np.where(D['rad'] > 0, D['Wy']/D['rad'], np.nan); Om_v = D['vac']/H**2
    out = {}
    neg = np.flatnonzero(H < -0.05)
    out['recollapse_H0tau'] = float(tau[neg[0]]) if len(neg) else None
    i0 = np.flatnonzero(H <= 0)
    out['H_zero_H0tau'] = float(np.interp(0.0, -H[i0[0]-1:i0[0]+1], tau[i0[0]-1:i0[0]+1])) if len(i0) and i0[0] > 0 else None
    p, j0 = plateau_index(lna, H, Wa4, Ra4, need_W=bool(res['match']['Weyl'] != 0.0))
    good = (Om_r >= 0.9) & (np.abs(r) <= 0.1)
    if p is not None:
        out['plateau'] = dict(H0tau=float(tau[p]), r=float(r[p]), Omega_r=float(Om_r[p]), Omega_vac=float(Om_v[p]), H_over_H0=float(H[p]))
    else: out['plateau'] = None
    if len(neg) and (p is None or p > neg[0]): cls = 'FAIL-no-radiation-era (recollapse)'
    elif p is not None:
        om, rr = out['plateau']['Omega_r'], out['plateau']['r']
        cls = ('FAIL-Weyl' if abs(rr) > 0.1 else 'FAIL-no-radiation-era (Omega_r < 0.5)' if om < 0.5 else
               ('PASS-combined' if abs(rr) <= 0.03 else 'PASS-conservative') if om >= 0.9 else 'INCONCLUSIVE (0.5<=Omega_r<0.9)')
    else: cls = 'INCONCLUSIVE (no plateau)'
    out['classification'] = cls
    # radiation era existence (tuning window): plateau with Omega_r >= 0.9 before any recollapse (r not required)
    out['rad_era'] = bool(p is not None and out['plateau']['Omega_r'] >= 0.9 and (not len(neg) or p < neg[0]))
    # >= 1 e-fold with Omega_r >= 0.9 and H > 0
    m = (Om_r >= 0.9) & (H > 0)
    best = 0.0; start = None
    for i in range(len(m)):
        if m[i] and start is None: start = i
        if (not m[i] or i == len(m) - 1) and start is not None:
            end = i if m[i] else i - 1; best = max(best, float(lna[end] - lna[start])); start = None
    out['max_efolds_Omega_r_ge_0.9'] = best
    return out

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--Y', type=float, default=1.0); ap.add_argument('--q', type=float, default=1.0)
    ap.add_argument('--dc', type=float, default=1e-2); ap.add_argument('--phim', type=float, default=0.95); ap.add_argument('--variant', default='registered')
    a = ap.parse_args()
    e = EFT()
    print(growth_rate(e))
    res = run(e, a.Y, a.q*DSTAR, a.dc, a.phim, variant=a.variant)
    print(res.get('match'), res.get('msg'))
    if res['ok']: print(analyse(res))
