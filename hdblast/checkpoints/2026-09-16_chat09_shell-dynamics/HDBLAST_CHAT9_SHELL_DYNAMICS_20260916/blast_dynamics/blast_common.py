#!/usr/bin/env python3
"""HDBLAST Chat 9 - blast dynamics: common tables for the shell-modulus EFT (numpy only).

EFT (orchestrator result R3):  S = int sqrt(-g)[ f R/2 - Z (d phi_b)^2/2 - V ],  f = 2I, V = sigma_t - 2W.
Einstein frame (M_pl = 1): g_E = f g_J,  V_E = V/f^2,  Z_E,y = 2 W D/f + 6 D^2/f^2  (kinetic coefficient for y_b),
with D := 1 - 2 W I/3 = -dI/dy_b.   We integrate the cancellation-free system
      I' = -D ,   D' = (2W/3) D - (2/3) W_y I ,   W_y = W_phi phi_y = W_phi^2
along BOTH BPS branches that are regular in the phi=-1 throat:
   'tanh' branch  phi = -tanh y  (y from +12 to -100; phi_b in (-1,1))
   'coth' branch  phi = -coth y  (y from +12 down to ~0.4; phi_b < -1)   <- continuation through phi_b = -1
Canonical field Theta: dTheta = sqrt(Z_E,y) |dy|, Theta = 0 at phi_b = 0, Theta increases toward/through the throat.
"""
import numpy as np
W  = lambda p: 1 - p + p**3/3
W1 = lambda p: p*p - 1
Ip = 1.0357712571566782
c_star = 2/Ip - 4/3
KD = 9.6*9/46          # D -> KD e^{-4y} deep in the throat (both branches)

def _flow(ys, phifun):
    n = len(ys); I = np.empty(n); D = np.empty(n)
    D[0] = KD*np.exp(-4*ys[0]); I[0] = 3*(1-D[0])/(2*W(phifun(ys[0])))
    def rhs(y, I_, D_):
        p = phifun(y); return -D_, (2*W(p)/3)*D_ - (2/3)*W1(p)**2*I_
    for i in range(n-1):
        y = ys[i]; h = ys[i+1]-ys[i]
        a1,b1 = rhs(y, I[i], D[i]); a2,b2 = rhs(y+h/2, I[i]+h/2*a1, D[i]+h/2*b1)
        a3,b3 = rhs(y+h/2, I[i]+h/2*a2, D[i]+h/2*b2); a4,b4 = rhs(y+h, I[i]+h*a3, D[i]+h*b3)
        I[i+1] = I[i] + h/6*(a1+2*a2+2*a3+a4); D[i+1] = D[i] + h/6*(b1+2*b2+2*b3+b4)
    return I, D

def build_tables(ymax=12.0, ymin=-100.0, dy=1e-3, ycoth_min=0.40, dyc=2.5e-4):
    # tanh branch
    ys = np.linspace(ymax, ymin, int(round((ymax-ymin)/dy))+1)
    I, D = _flow(ys, lambda y: -np.tanh(y))
    phi = -np.tanh(ys); f = 2*I; ZEy = 2*W(phi)*D/f + 6*D**2/f**2
    s = np.sqrt(ZEy)
    cum = np.concatenate([[0], np.cumsum(0.5*(s[1:]+s[:-1])*np.abs(np.diff(ys)))])   # from y=ymax downward
    # distance from y=ymax to +infinity: s ~ s0 e^{-2(y-ymax)} -> s0/2
    tail_throat = s[0]/2
    i0 = np.argmin(np.abs(ys)); Theta = -(cum - cum[i0])
    Theta_end_throat = Theta[0] + tail_throat
    T = dict(y=ys, phi=phi, I=I, D=D, f=f, ZEy=ZEy, Theta=Theta, dThdy=s)       # dTheta/dy = +s on this branch
    # coth branch
    yc = np.linspace(ymax, ycoth_min, int(round((ymax-ycoth_min)/dyc))+1)
    Ic, Dc = _flow(yc, lambda y: -1/np.tanh(y))
    phic = -1/np.tanh(yc); fc = 2*Ic; ZEc = 2*W(phic)*Dc/fc + 6*Dc**2/fc**2
    ok = np.ones(len(yc), bool)
    bad = np.where((ZEc <= 0) | (fc <= 0))[0]
    if len(bad): ok[bad[0]:] = False
    sc = np.sqrt(np.where(ok, ZEc, 0))
    cumc = np.concatenate([[0], np.cumsum(0.5*(sc[1:]+sc[:-1])*np.abs(np.diff(yc)))])
    Thc = Theta_end_throat + sc[0]/2 + cumc
    C = dict(y=yc, phi=phic, I=Ic, D=Dc, f=fc, ZEy=ZEc, Theta=Thc, dThdy=-sc, ok=ok)
    return T, C, Theta_end_throat

class Landscape:
    """Merged monotone table in Theta with a user potential V(phi)/t (and dV/dphi)."""
    def __init__(self, Vfun, dVfun, T=None, C=None, Th_end=None):
        if T is None: T, C, Th_end = build_tables()
        self.T, self.C, self.Th_throat = T, C, Th_end
        m = C['ok']
        # order by increasing Theta: tanh branch reversed (y from -100 up to 12), then coth branch (y from 12 down)
        def cat(key, sign=1):
            return np.concatenate([T[key][::-1], C[key][m]])
        self.Th = cat('Theta'); self.phi = cat('phi'); self.f = cat('f'); self.D = cat('D'); self.y = cat('y')
        self.dThdy = cat('dThdy'); self.branch = np.concatenate([np.zeros(len(T['y'])), np.ones(m.sum())])
        phi_y = W1(self.phi)                       # BPS: phi_y = W_phi on both branches
        f_y = -2*self.D
        self.v = Vfun(self.phi)/self.f**2
        dvdy = dVfun(self.phi)*phi_y/self.f**2 - 2*Vfun(self.phi)*f_y/self.f**3
        self.dv = dvdy/self.dThdy
        self.fTh = f_y/self.dThdy                  # df/dTheta
        self.phiTh = phi_y/self.dThdy              # dphi_b/dTheta
        assert np.all(np.diff(self.Th) > 0), "Theta table not monotone"
        self.Th_min = self.Th[0]
        # analytic tail toward phi_b -> +1 :  f = 9 - (3/2) dTh^2, Theta_end = Th_min - dTh_min
        self.dTh_min = np.sqrt((9 - self.f[0])/1.5); self.Th_plus_end = self.Th_min - self.dTh_min
        self.i0 = np.argmin(np.abs(self.Th)); self.v0 = np.interp(0.0, self.Th, self.v)
    def interp(self, th, arr): return np.interp(th, self.Th, arr)

V_reg  = lambda p: 1 + c_star*p
dV_reg = lambda p: c_star + 0*p

if __name__ == "__main__":
    import json, os
    T, C, Th_end = build_tables()
    L = Landscape(V_reg, dV_reg, T, C, Th_end)
    out = {}
    # validation against orchestrator npz
    d = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "effective_theory_orchestrator", "eft_landscape.npz"))
    for yy in [8.0, 3.0, 1.0, 0.0, -1.0, -5.0, -20.0, -30.0]:
        i = np.argmin(np.abs(d['y']-yy)); j = np.argmin(np.abs(T['y']-yy))
        print("y=%+6.1f  Theta orch %.8f mine %.8f | f orch %.10f mine %.10f | V orch %.8f mine %.8f" %
              (yy, d['Theta'][i], T['Theta'][j], d['f'][i], T['f'][j], d['V'][i], V_reg(T['phi'][j])/T['f'][j]**2))
    out['I_plus_check'] = float(T['I'][np.argmin(np.abs(T['y']))])
    out['Theta_throat_end(phi_b=-1)'] = float(Th_end)
    out['Theta_plus_end(phi_b=+1)'] = float(L.Th_plus_end)
    # Z(phi=-1): zero-mode norm 2*int e^{-(10/9)s} e^{-4s} ds = 9/23
    j = np.argmin(np.abs(T['y']-8.0)); Zphi = 2*W(T['phi'][j])*T['D'][j]/W1(T['phi'][j])**2
    out['Z(phi_b->-1) numeric'] = float(Zphi); out['Z(phi_b=-1) exact 9/23'] = 9/23
    # smoothness across phi_b=-1: dv/dTheta on both sides
    jt = np.argmin(np.abs(T['y']-6.0)); jc = np.argmin(np.abs(C['y']-6.0))
    it = np.argmin(np.abs(L.Th-T['Theta'][jt])); ic = np.argmin(np.abs(L.Th-C['Theta'][jc]))
    out['dlnv/dTheta just before/after phi_b=-1'] = [float(L.dv[it]/L.v[it]), float(L.dv[ic]/L.v[ic])]
    out['dphi/dTheta just before/after'] = [float(L.phiTh[it]), float(L.phiTh[ic])]
    # +1 end: m^2/H^2 -> 2
    k = 5000; th = L.Th[:k]; d2 = np.gradient(L.dv[:k], th)
    out['3 v_ThTh / v at phi_b->+1 (expect 2)'] = float(3*d2[100]/L.v[100])
    out['v_inf=(1+c)/81'] = (1+c_star)/81; out['v at table end'] = float(L.v[0]); out['v0 hilltop'] = float(L.v0)
    out['mu2 hilltop from table 3 v''/v'] = float(3*np.gradient(L.dv, L.Th)[L.i0]/L.v0)
    # coth branch landmarks
    m = C['ok']
    out['coth branch: phi_b range covered'] = [float(C['phi'][m][0]), float(C['phi'][m][-1])]
    out['coth branch: first failure (ZEy<=0 or f<=0) at phi_b'] = (float(C['phi'][~m][0]) if (~m).any() else None)
    for ph in [-1.2, -1.5, -1/c_star, -2.0, -2.104, -2.3]:
        j = np.argmin(np.abs(C['phi']-ph))
        if m[j]: print("coth: phi_b=%.4f y=%.4f f=%.5f D=%.5f ZEy=%.5f Theta=%.5f V_E/t=%.5f" % (C['phi'][j], C['y'][j], C['f'][j], C['D'][j], C['ZEy'][j], C['Theta'][j], V_reg(C['phi'][j])/C['f'][j]**2))
    out['Theta at V=0 (phi_b=-1/c)'] = float(C['Theta'][np.argmin(np.abs(C['phi']+1/c_star))])
    print(json.dumps(out, indent=1))
    json.dump(out, open("BLAST_TABLE_CHECKS.json", "w"), indent=1)
