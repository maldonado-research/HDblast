#!/usr/bin/env python3
"""'Resummed' use of the second-order curvature expansion: solve the two junction conditions, truncated at O(X^2),
for (phi_b, X_b) at finite t WITHOUT expanding in t, then mu^2 = -4 - [T0 + X T1]/(3 [s0 + X s1(mu^2)]).
Error is O(t^2) with O(1) coefficients even when the modulus potential is tuned flat (where the plain t-series has a
small radius of convergence ~ mu0^2)."""
import sys, math, json, re
import numpy as np
import slope_formula as sf

class Table:
    def __init__(self, lo, hi, n=11):
        x = np.cos(np.pi*(np.arange(n) + 0.5)/n); self.lo, self.hi = lo, hi
        ph = [float(0.5*(lo + hi) + 0.5*(hi - lo)*xi) for xi in x]
        tab = sf.wall(ph)
        self.co = [np.polynomial.chebyshev.chebfit(x, [tab[p][k] for p in ph], n - 1) for k in range(3)]
    def __call__(self, p):
        x = (2*p - self.lo - self.hi)/(self.hi - self.lo)
        return tuple(float(np.polynomial.chebyshev.chebval(x, c)) for c in self.co)

def predict(t, Vf, tab, phi_guess):
    def junction(p):
        L = sf.local(p, *tab(p), Vf); V, Vp, _ = Vf(p); I, z2 = L['I'], L['z2']
        X = (-I + math.sqrt(I*I + 4*z2*t*V/6))/(2*z2)          # X I + X^2 z2 = t V/6
        return (L['p1'] + X*L['p2'])*V + 3*Vp*(I + X*z2), X, L   # (X p1 + X^2 p2) = -t V'/2 divided by first junction
    a = phi_guess; e = 1e-6
    for _ in range(60):
        fa = junction(a)[0]; d = (junction(a + e)[0] - junction(a - e)[0])/(2*e); step = -fa/d; a += step
        if abs(step) < 1e-15: break
    _, X, L = junction(a)
    mu = L['mu0']
    for _ in range(50):
        s = L['s0'] + X*(L['s1a'] + (3*mu + 12)*L['s1b'])
        mu = -4 - (L['T0'] + X*L['T1'])/(3*s)
    mu_first = -4 - L['T0']/(3*L['s0'])        # two-derivative EFT evaluated at the *true* second-order position
    return dict(phi_b=a, X_b=X, mu2=mu, mu2_lowest_order_at_phi_b=mu_first)

if __name__ == "__main__":
    out = []
    for tag, log in [("slowroll d=1.105924", "log_slowroll_first.txt"), ("slowroll small t", "log_slowroll.txt"), ("c_star", "log_cstar.txt"), ("d=1", "log_d1.txt")]:
        head = open(log).readline(); d = float(re.search(r" d=(\S+)", head).group(1)); d = 1.105924 if abs(d - 1.105924) < 1e-5 else d; phi0 = float(re.search(r"phi0=(\S+)", head).group(1))  # (header printed d with %g) phi0 = float(re.search(r"phi0=(\S+)", head).group(1))
        r = sf.slope(phi0, lambda p, I, I1: sf.linear_detuning(p, I, I1, d=d)); Vf = r['Vf']
        tab = Table(phi0 - 0.01, phi0 + 0.03)
        print(tag, " mu0=%.8f slope=%.6f" % (r['mu0'], r['slope']))
        for line in open(log):
            m = re.match(r"t=(\S+) phi_b=(\S+) .*X_b/t=(\S+) .* mu2=(\S+) \[", line)
            if not m: continue
            t, pb, xb, mu5 = [float(m.group(i)) for i in (1, 2, 3, 4)]
            p = predict(t, Vf, tab, phi0 + r['beta']*t/(1 + 300*t*abs(r['beta'])))
            print("  t=%.1e  5D: phi_b=%+.7e X/t=%.8f mu2=%.9f | resummed: phi_b=%+.7e X/t=%.8f mu2=%.9f (diff %.2e, diff/t^2=%.3f) | linear series %.9f (diff %.1e)" % (
                t, pb, xb, mu5, p['phi_b'], p['X_b']/t, p['mu2'], p['mu2'] - mu5, (p['mu2'] - mu5)/t**2, r['mu0'] + r['slope']*t, r['mu0'] + r['slope']*t - mu5))
            out.append(dict(case=tag, t=t, mu2_5D=mu5, mu2_resummed=p['mu2'], mu2_linear=r['mu0'] + r['slope']*t, phi_b_5D=pb, phi_b_resummed=p['phi_b']))
    json.dump(out, open("RESUMMED_VS_5D.json", "w"), indent=1)
