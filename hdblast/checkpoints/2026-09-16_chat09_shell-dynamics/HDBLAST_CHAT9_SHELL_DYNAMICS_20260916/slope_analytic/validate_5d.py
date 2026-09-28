#!/usr/bin/env python3
"""5D validation of the analytic slope: full background (cone -> shell) Newton solve + vectorised shooting of the master
equation with the scalar junction, for several t at fixed detuning profile; polynomial extrapolation in t.
usage: python3 validate_5d.py phi0 d tag    (phi0 = leading-order equilibrium position, d = quadratic detuning coefficient)"""
import sys, math, json, os
N1 = int(os.environ.get('NSTEP', '8000'))
import numpy as np
HERE = "/Users/ricardomaldonado/Documents/D-Blast 3/untitled folder 146/HDBLAST_CHAT9_SHELL_DYNAMICS_20260916/"
sys.path.insert(0, HERE + "background"); sys.path.insert(0, HERE + "slope_analytic")
import hdblast_background as bg
import slope_formula as sf

def residuals(phi_h, v_b, t, V, dv):
    rho, phi, s = bg.integrate(phi_h, v_b, dv)
    rp = math.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6))
    v, vp, vpp = V(phi)
    return np.array([rp/rho - (2*bg.W(phi) + t*v)/6, s + (2*bg.W1(phi) + t*vp)/2]), (rho, phi, s)

def solve(t, V, guess, dv=5e-4, tol=1e-13):
    x = np.array(guess, float)
    f = lambda xx: residuals(-1 + 10**xx[0], xx[1], t, V, dv)[0]
    for it in range(40):
        F = f(x); J = np.zeros((2, 2))
        for j in range(2):
            xp = x.copy(); xp[j] += 1e-7; xm = x.copy(); xm[j] -= 1e-7
            J[:, j] = (f(xp) - f(xm))/2e-7
        dx = np.linalg.solve(J, -F); lam = 1.0
        while lam > 1e-4 and np.linalg.norm(f(x + lam*dx)) > np.linalg.norm(F): lam /= 2
        x = x + lam*dx
        if np.linalg.norm(lam*dx) < tol: break
    res, (rho, phi, s) = residuals(-1 + 10**x[0], x[1], t, V, dv)
    return dict(phi_h=-1 + 10**x[0], y_b=x[1], rho_b=rho, phi_b=phi, s_b=s, resid=[float(res[0]), float(res[1])], x=[float(x[0]), float(x[1])])

def shoot(mu2, phi_h, y_b, t, V, y0=1e-3, n=12000, want_ratio=False):
    mu2 = np.asarray(mu2, float); M = mu2.size
    alpha = 0.5 + np.sqrt(2.25 - mu2)
    b = bg.series(phi_h, y0)
    Y = np.zeros((5, M)); Y[0] = b[0]; Y[1] = b[1]; Y[2] = b[2]; Y[3] = 1.0; Y[4] = alpha/y0
    def rhs(Y):
        rho, phi, s, psi, dpsi = Y
        rp = np.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); Hh = rp/rho
        spp = bg.U1(phi) - 4*Hh*s; g = spp/s
        return np.array([rp, s, spp, dpsi, -2*(Hh - g)*dpsi - (-(4/3)*s*s - 4*Hh*g + (2 + mu2)/rho**2)*psi])
    dv = (y_b - y0)/n
    for _ in range(n):
        k1 = rhs(Y); k2 = rhs(Y + dv/2*k1); k3 = rhs(Y + dv/2*k2); k4 = rhs(Y + dv*k3)
        Y = Y + dv/6*(k1 + 2*k2 + 2*k3 + k4)
    rho, phi, s, psi, dpsi = Y
    rp = np.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6)); Hh = rp/rho
    spp = bg.U1(phi) - 4*Hh*s
    chi = -3*(dpsi + 2*Hh*psi)/s
    vpp = V(phi)[2]
    # junction in the form (3 mu^2 + 12) psi / rho^2 + (phi'' + phi' sigma''/2) chi = 0   (equivalent to R1's BC)
    Bc = (3*mu2 + 12)*psi/rho**2 + (spp + s*(2*bg.W2(phi) + t*vpp)/2)*chi
    if want_ratio: return Bc, psi/chi
    return Bc

def root(f, a, b):
    for _ in range(7):
        g = np.linspace(a, b, 17); v = f(g)
        j = np.where(v[:-1]*v[1:] <= 0)[0]
        if len(j) == 0: raise RuntimeError("no bracket in [%g,%g]" % (a, b))
        a, b = g[j[0]], g[j[0]+1]
    fa, fb = f(np.array([a, b]))
    return a - fa*(b - a)/(fb - fa)

def throat_kappa():
    """probe scalar on AdS_- (dS slicing): delta'' + (4/l) coth(y/l) delta' = (76/9) delta, delta(0)=1 -> kappa e^{2y}."""
    l = 9/5; m2 = 76/9; y = 1e-3; d = np.array([1 + m2*y*y/10, m2*y/5]); hh = 1e-3
    f = lambda y, d: np.array([d[1], m2*d[0] - 4/(l*math.tanh(y/l))*d[1]])
    while y < 14:
        k1 = f(y, d); k2 = f(y + hh/2, d + hh/2*k1); k3 = f(y + hh/2, d + hh/2*k2); k4 = f(y + hh, d + hh*k3)
        d = d + hh/6*(k1 + 2*k2 + 2*k3 + k4); y += hh
    return d[0]*math.exp(-2*y)

def guess(t, phi0, x1, beta, kappa):
    l = 9/5; pb = phi0 + beta*t; yh = math.atanh(pb); rho_b = 1/math.sqrt(x1*t)
    A = (yh - (2/3)*math.log(math.cosh(yh)) + 1/(6*math.cosh(yh)**2))/3
    ysh = l*(math.log(2*rho_b/l) + (2/9)*math.log(2) - A)
    return (math.log10(2*math.exp(-2*ysh)/kappa), yh + ysh)

if __name__ == "__main__":
    phi0 = float(sys.argv[1]); d = float(sys.argv[2]); tag = sys.argv[3]
    ts = [float(a) for a in sys.argv[4:]] or [1e-3, 2e-3, 4e-3, 8e-3, 1.6e-2]
    r = sf.slope(phi0, lambda p, I, I1: sf.linear_detuning(p, I, I1, d=d))
    V = r['Vf']; c = V.c
    print("phi0=%g d=%g c=%.12f  mu0=%.10f  analytic slope=%.8f  beta=%.8f" % (phi0, d, c, r['mu0'], r['slope'], r['beta']), flush=True)
    kap = throat_kappa(); rows = []
    for t in ts:
        sol = solve(t, V, guess(t, phi0, r['x1'], r['beta'], kap))
        X = 1/sol['rho_b']**2
        cen = r['mu0'] + r['slope']*t
        mu = {}
        for n in (N1, 2*N1):
            f = lambda m: shoot(m, sol['phi_h'], sol['y_b'], t, V, n=n)
            mu[n] = root(f, cen - 0.05 - 30*t, cen + 0.05 + 30*t)
        mu_R = mu[2*N1] + (mu[2*N1] - mu[N1])/15
        _, sb = shoot(np.array([mu_R]), sol['phi_h'], sol['y_b'], t, V, n=2*N1, want_ratio=True)
        s_pred = r['s0'] + 0  # leading; full comparison below
        row = dict(t=t, phi_b=sol['phi_b'], X_b=X, resid=sol['resid'], nsteps=[N1, 2*N1], mu2_n1=mu[N1], mu2_n2=mu[2*N1], mu2=mu_R,
                   slope_est=(mu_R - r['mu0'])/t, phi_b_pred=phi0 + r['beta']*t, psi_over_chi=float(sb[0]))
        rows.append(row)
        print("t=%.1e phi_b=%+.6e (pred1 %+.6e) X_b/t=%.8f (x1=%.8f) resid=(%.0e,%.0e) mu2=%.10f [n-halving diff %.1e] (mu2-mu0)/t=%.7f  psi/chi=%.8f" % (
            t, sol['phi_b'], row['phi_b_pred'], X/t, r['x1'], sol['resid'][0], sol['resid'][1], mu_R, mu[2*N1]-mu[N1], row['slope_est'], sb[0]), flush=True)
    tt = np.array([q['t'] for q in rows]); ss = np.array([q['slope_est'] for q in rows])
    ext = {}
    for deg in range(1, len(tt)):
        co = np.polyfit(tt[:deg+1], ss[:deg+1], deg); ext[deg] = float(co[-1])
    print("t->0 extrapolations of (mu2-mu0)/t by degree (lowest t's):", ext, " analytic:", r['slope'], flush=True)
    json.dump(dict(phi0=phi0, d=d, c=c, mu0=r['mu0'], analytic_slope=r['slope'], beta=r['beta'], x1=r['x1'], s0=r['s0'],
                   s1_total=r['s1_total'], rows=rows, extrapolated_slope=ext), open("VALIDATE_%s.json" % tag, "w"), indent=1)
