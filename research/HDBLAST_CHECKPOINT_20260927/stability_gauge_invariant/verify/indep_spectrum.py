"""Independent auditor re-computation of the scalar/tensor shell spectrum (does NOT import gi_core.py).

Differences from the audited code:
  * background: (rho, e, e') with rho' = sqrt(1 + rho^2(e'^2/12 - U/6)) (first integral, derived in indep_derivation.py).
    NOTE: a first attempt integrating rho'' as a free 4th variable was abandoned: the constraint
    C = H^2 - 1/rho^2 - ... is conserved, so O(1e-16) round-off at the cone (where H^2 ~ 1/y0^2) becomes a
    relative error ~1e-6 at the shell; the constraint form (as in the audited code) is required;
  * static shell re-solved with my own integrator (scipy root on both junctions), starting from the audited (e_h, y_b);
  * scalar sector integrated as a single Pruefer angle theta for (X, c0 Z):
        theta' = -2(H+g) sin cos - (c0/3) cos^2 - ((3 lam/rho^2 - 2 phi'^2)/c0) sin^2 ,   g = phi''/phi'
    (no overflow for very negative mu2, no normalisation), with regular cone data tan(theta0) = -c0 y0/(3(s+4));
    shell condition  B cos(theta) + (3 lam/(c0 rho^2)) sin(theta) = 0 ;
  * tensor sector as a Pruefer angle for (h, rho h'):  theta' = -3H s c - (mu2/rho) c^2 - s^2/rho, Neumann <=> sin(theta_b) = 0;
  * integrators: Radau (implicit) and DOP853 as the two precision variants.
Output: indep_spectrum.json
"""
import sys, json, math, time, hashlib
from pathlib import Path
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import root, brentq

HERE = Path(__file__).resolve().parent
T0 = time.time()
C_REG = 0.5975949350280132
AUD = json.loads((HERE.parent/'BACKGROUNDS.json').read_text())

def upoly(v):
    e, f = sp.symbols('e f')
    W = 1 - f + f**3/3
    U = sp.expand((sp.diff(W, f)**2/2 - sp.Rational(2, 3)*W**2).subs(f, v + e))
    return [float(sp.Poly(U, e).coeff_monomial(e**k)) for k in range(7)]

class BG:
    def __init__(self, v, delta, c=C_REG):
        self.v, self.d, self.c = v, delta, c
        a = upoly(v); self.a = np.array(a)
        self.a1 = np.array([k*a[k] for k in range(1, 7)])
    def U(self, e): return np.polyval(self.a[::-1], e)
    def U1(self, e): return np.polyval(self.a1[::-1], e)
    def sig(self, e):
        f = self.v + e; return 2*(1 - f + f**3/3) + self.d*(1 + self.c*f)
    def sig1(self, e): return 2*(2*self.v*e + e*e) + self.d*self.c
    def sig2(self, e): return 2*2*(self.v + e)
    def start(self, eh, y0):
        u, u1 = self.U(eh), self.U1(eh)
        # leading Frobenius terms (derived: rho = y - U y^3/36, e = e_h + U1 y^2/10)
        return [y0 - u*y0**3/36, eh + u1*y0**2/10, u1*y0/5]
    def rp(self, r, e, ep):
        return math.sqrt(1 + r*r*(ep*ep/12 - self.U(e)/6))
    def rhs(self, y, s):
        r, e, ep = s; rp = self.rp(r, e, ep)
        return [rp, ep, self.U1(e) - 4*rp*ep/r]
    def run(self, eh, yb, y0=1e-4, rtol=1e-12, method='DOP853', dense=False):
        atol = [1e-15, abs(eh)*1e-14 + 1e-300, abs(eh)*1e-14 + 1e-300]
        sol = solve_ivp(self.rhs, (y0, yb), self.start(eh, y0), method=method, rtol=rtol, atol=atol, dense_output=dense)
        assert sol.success, sol.message
        return sol
    def junction(self, eh, yb, **kw):
        r, e, ep = self.run(eh, yb, **kw).y[:, -1]; rp = self.rp(r, e, ep)
        return np.array([rp/r - self.sig(e)/6, ep + self.sig1(e)/2])
    def solve_static(self, eh0, yb0):
        sg = math.copysign(1, eh0)
        F = lambda x: self.junction(sg*10**x[0], x[1])/max(self.d, 1e-3)
        fit = root(F, [math.log10(abs(eh0)), yb0], tol=1e-14, options={'eps': 1e-9})
        self.eh, self.yb = sg*10**fit.x[0], fit.x[1]
        r, e, ep = self.run(self.eh, self.yb).y[:, -1]; rp = self.rp(r, e, ep)
        U = self.U(e)
        self.shell = dict(e_h=self.eh, y_b=self.yb, rho_b=r, phi_b=self.v + e, eta_b=e, phi1_b=ep,
                          H2=1/r**2, res=self.junction(self.eh, self.yb).tolist(),
                          constraint_rel=(rp*rp - 1 - r*r*(ep*ep/12 - U/6))/(rp*rp),
                          g_b=self.U1(e)/ep - 4*rp/r, B=self.U1(e)/ep - 4*rp/r + self.sig2(e)/2)
        # sign of phi' along the bulk (g = phi''/phi' must be regular)
        sol = self.run(self.eh, self.yb, dense=True)
        ys = np.linspace(1e-4, self.yb, 4001)
        epv = sol.sol(ys)[2]
        self.shell['phi1_sign_changes'] = int(np.sum(np.sign(epv[1:]) != np.sign(epv[:-1])))
        return self.shell

    # ---------------- scalar Pruefer angle
    def scalar_theta(self, mu2, c0=10.0, y0=1e-4, rtol=1e-12, method='DOP853'):
        lam = mu2 + 4; s = -1.5 + math.sqrt(2.25 - mu2)
        th0 = math.atan(-c0*y0/(3*(s + 4)))
        st = self.start(self.eh, y0) + [th0]
        U1 = self.U1
        def f(y, q):
            r, e, ep, th = q
            rp = self.rp(r, e, ep); H = rp/r; epp = U1(e) - 4*H*ep; g = epp/ep
            sn, cs = math.sin(th), math.cos(th)
            return [rp, ep, epp,
                    -2*(H + g)*sn*cs - (c0/3)*cs*cs - ((3*lam/(r*r) - 2*ep*ep)/c0)*sn*sn]
        atol = [1e-15, abs(self.eh)*1e-14 + 1e-300, abs(self.eh)*1e-14 + 1e-300, 1e-13]
        sol = solve_ivp(f, (y0, self.yb), st, method=method, rtol=rtol, atol=atol)
        assert sol.success, sol.message
        r, e, ep, th = sol.y[:, -1]; rp = self.rp(r, e, ep)
        B = U1(e)/ep - 4*rp/r + self.sig2(e)/2
        a, b = B*math.cos(th), 3*lam/(c0*r*r)*math.sin(th)
        # Bstar: the value the shell coefficient B would need for mu2 to be an eigenvalue (stability margin control)
        self._Bstar = -3*lam*math.tan(th)/(c0*r*r)
        return (a + b)/(abs(a) + abs(b)), th
    # ---------------- tensor Pruefer angle
    def tensor_theta(self, mu2, y0=1e-4, rtol=1e-12, method='DOP853'):
        s = -1.5 + math.sqrt(2.25 - mu2)
        st = self.start(self.eh, y0) + [math.atan(s)]
        def f(y, q):
            r, e, ep, th = q
            rp = self.rp(r, e, ep); H = rp/r; sn, cs = math.sin(th), math.cos(th)
            return [rp, ep, self.U1(e) - 4*H*ep, -3*H*sn*cs - (mu2/r)*cs*cs - sn*sn/r]
        atol = [1e-15, abs(self.eh)*1e-14 + 1e-300, abs(self.eh)*1e-14 + 1e-300, 1e-13]
        sol = solve_ivp(f, (y0, self.yb), st, method=method, rtol=rtol, atol=atol)
        return math.sin(sol.y[3, -1])

def scan(bg, grid, fun):
    vals = np.array([fun(m) for m in grid]); roots = []
    for i in range(len(grid) - 1):
        if vals[i]*vals[i + 1] < 0:
            roots.append(brentq(fun, grid[i], grid[i + 1], xtol=1e-12))
    return vals, roots

if __name__ == '__main__':
    quick = '--quick' in sys.argv
    GRID = np.concatenate([-np.logspace(math.log10(400), math.log10(4.05), 30), np.linspace(-4.0, 2.2499, 126 if quick else 251)])
    TGRID = np.linspace(-4.0, 2.2499, 64)
    out = dict(grid_range=[float(GRID[0]), float(GRID[-1])], n_grid=len(GRID), branches=[])
    cases = [('original', -1, d) for d in (0.001, 0.003, 0.01)] + [('plus', +1, d) for d in (0.0003, 0.001, 0.003, 0.01, 0.03, 0.1)]
    audp = {('plus', r['delta']): r['gi_core'] for r in AUD['plus']}
    audp.update({('original', r['delta']): r['gi_core'] for r in AUD['original']})
    for name, v, d in cases:
        t1 = time.time()
        a = audp[(name, d)]
        bg = BG(v, d); sh = bg.solve_static(a['e_h'], a['y_b'])
        rec = dict(branch=name, delta=d, shell=sh, rel_diff_rho_b_vs_audit=abs(sh['rho_b']/a['rho_b'] - 1),
                   diff_B_vs_audit=sh['B'] - a['B'])
        Bst = {}
        def M(m):
            v_ = bg.scalar_theta(m)[0]; Bst[float(m)] = bg._Bstar; return v_
        vals, roots = scan(bg, GRID, M)
        neg = [Bst[float(m)] for m in GRID if m < 0]
        rec['Bstar_control'] = dict(note='mu2<0 is an eigenvalue iff B = Bstar(mu2); unstable roots need B inside [min,max] of Bstar over mu2<0',
                                    Bstar_min_mu2_neg=float(min(neg)), Bstar_max_mu2_neg=float(max(neg)), B_actual=sh['B'],
                                    margin=float(sh['B'] - max(neg)))
        rec['scalar'] = dict(min_Mhat=float(vals.min()), max_Mhat=float(vals.max()), roots=roots,
                             growth=[-1.5 + math.sqrt(2.25 - r) for r in roots])
        # precision variants at a few points + for roots
        pts = [-400.0, -10.0, -4.0, -2.0, -1.0, 0.0, 2.0, 2.2499]
        var = {}
        for m in pts:
            vv = [bg.scalar_theta(m, rtol=rt, y0=y0, method=me, c0=c0)[0]
                  for (rt, y0, me, c0) in [(1e-10, 1e-3, 'DOP853', 10.0), (1e-12, 1e-4, 'DOP853', 10.0), (1e-13, 1e-5, 'DOP853', 3.0), (1e-10, 1e-4, 'Radau', 10.0)]]
            var[str(m)] = [float(min(vv)), float(max(vv))]
        rec['scalar']['variants_min_max'] = var
        if roots:
            rec['scalar']['root_variants'] = {repr(r): [brentq(lambda m: bg.scalar_theta(m, rtol=rt, y0=y0, method=me, c0=c0)[0], r - 0.05, r + 0.05, xtol=1e-12)
                                                        for (rt, y0, me, c0) in [(1e-10, 1e-3, 'DOP853', 10.0), (1e-13, 1e-5, 'DOP853', 3.0), (1e-10, 1e-4, 'Radau', 10.0)]]
                                              for r in roots}
        tv, troots = scan(bg, TGRID, lambda m: bg.tensor_theta(m))
        rec['tensor'] = dict(roots=troots, sin_theta_min=float(tv.min()), sin_theta_max=float(tv.max()),
                             tensor_at_mu2_0=bg.tensor_theta(0.0))
        rec['runtime_s'] = time.time() - t1
        print(json.dumps(dict(branch=name, delta=d, roots=roots, min=rec['scalar']['min_Mhat'], B=sh['B'], troots=troots,
                              res=sh['res'], rhodiff=rec['rel_diff_rho_b_vs_audit'], t=rec['runtime_s'])), flush=True)
        out['branches'].append(rec)
    o = [r for r in out['branches'] if r['branch'] == 'original' and r['delta'] == 0.001][0]
    out['calibration'] = dict(roots=o['scalar']['roots'], recorded_chat9=-7.717871625176294,
                              audited=-7.717871625260236,
                              diff_vs_recorded=(o['scalar']['roots'][0] + 7.717871625176294) if len(o['scalar']['roots']) == 1 else None)
    out['plus_any_root'] = any(r['scalar']['roots'] for r in out['branches'] if r['branch'] == 'plus')
    out['runtime_s'] = time.time() - T0
    out['sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(HERE/'indep_spectrum.json').write_text(json.dumps(out, indent=2, default=float) + '\n')
    print('calibration', out['calibration']); print('plus_any_root', out['plus_any_root'], 'runtime', out['runtime_s'])
