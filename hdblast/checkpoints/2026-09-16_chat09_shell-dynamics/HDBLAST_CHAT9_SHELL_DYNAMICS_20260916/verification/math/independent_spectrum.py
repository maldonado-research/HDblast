#!/usr/bin/env python3
"""VERIFICATION (math agent) - independent numerics for the scalar-sector spectrum of the registered dS shell.

Differences from the orchestrator code:
  * own background solver (own cone series start y0 = 2e-3, graded step: uniform RK4 in x = ln y + y, own Newton);
  * first-order system in (psi, Zeta = chi/phi')  -- no phi''/phi' anywhere, no psi'':
        psi'  = -2 (rho'/rho) psi - phi'^2 Zeta/3
        Zeta' = 3 (mu2+4) psi/(rho^2 phi'^2) - 2 psi
    (validated against the full linearised Einstein equations by einstein_linear_check.py);
  * shell condition in the equivalent algebraic form  B = 3(mu2+4) psi/(rho^2 phi') + (phi'' + sigma''phi'/2) Zeta = 0;
  * vectorised over (complex) mu2; analytic normalisation B (y0/y_b)^alpha for the argument principle
    (the orchestrator divided by psi_b, which counts zeros MINUS zeros of psi_b).
Also computes the residue of the brane-to-brane response G = chi_b/j (source j(phi-phi_b0) on the shell => B = j/2)
and compares with the healthy-scalar EFT prediction  Res G = -rho_b^2/(f Z_E),  f = 2 I_plus, Z_E = c/2 + 3c^2/8.
A NEGATIVE residue = positive kinetic norm (tachyon, not ghost).
"""
import numpy as np, math, json, sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))

W  = lambda p: 1 - p + p**3/3
W1 = lambda p: p*p - 1
W2 = lambda p: 2*p
U  = lambda p: 0.5*W1(p)**2 - (2/3)*W(p)**2
U1 = lambda p: W1(p)*W2(p) - (4/3)*W(p)*W1(p)
U2 = lambda p: W2(p)**2 + 2*W1(p) - (4/3)*(W1(p)**2 + W(p)*W2(p))

def I_plus(n=200000, L=60.0):
    # own quadrature: composite Simpson on y in [0, L] of exp(2A), A = -(1/3) int_0^y W(-tanh)  (closed form integrated by hand)
    y = np.linspace(0, L, n + 1)
    # int_0^y W(-tanh u) du = y + ln cosh y - (1/3)[ ln cosh y - tanh^2 y /2 ]   (int tanh^3 = ln cosh - tanh^2/2)
    lc = np.logaddexp(y, -y) - math.log(2)
    intW = y + lc - (lc - np.tanh(y)**2/2)/3
    g = np.exp(-2*intW/3)
    h = L/n
    return h/3*(g[0] + g[-1] + 4*g[1:-1:2].sum() + 2*g[2:-1:2].sum())

Y0 = 2e-3
def cone(phi_h, y0=Y0):
    u, u1, u2 = U(phi_h), U1(phi_h), U2(phi_h)
    b = u1*(u2/280 + u/630)
    return y0 - u*y0**3/36, phi_h + u1*y0**2/10 + b*y0**4, u1*y0/5 + 4*b*y0**3

def integrate(phi_h, y_b, mu2=None, N=20000):
    """RK4 uniform in x = ln y + y. Returns background at y_b and (psi, Zeta) arrays if mu2 given."""
    rho, phi, s = cone(phi_h)
    y = Y0; x0 = math.log(Y0) + Y0; x1 = math.log(y_b) + y_b; h = (x1 - x0)/N
    pert = mu2 is not None
    if pert:
        mu2 = np.asarray(mu2)
        alpha = 0.5 + np.sqrt(2.25 - mu2)
        rp = math.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6))
        psi = np.ones_like(alpha); zeta = -3*(alpha/Y0 + 2*rp/rho)/s**2
    def f(y, rho, phi, s, psi, zeta):
        J = y/(1 + y)
        rp = math.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6)); Hh = rp/rho
        out = [J, J*rp, J*s, J*(U1(phi) - 4*Hh*s)]
        if pert:
            out += [J*(-2*Hh*psi - s*s*zeta/3), J*(3*(mu2 + 4)*psi/(rho*rho*s*s) - 2*psi)]
        else:
            out += [0.0, 0.0]
        return out
    st = [y, rho, phi, s, psi if pert else 0.0, zeta if pert else 0.0]
    for _ in range(N):
        k1 = f(*st)
        k2 = f(*[a + h/2*b for a, b in zip(st, k1)])
        k3 = f(*[a + h/2*b for a, b in zip(st, k2)])
        k4 = f(*[a + h*b for a, b in zip(st, k3)])
        st = [a + h/6*(b + 2*c + 2*d + e) for a, b, c, d, e in zip(st, k1, k2, k3, k4)]
    return st

def shell_res(x, t, c, N):
    phi_h = -1 + 10**x[0]
    y, rho, phi, s, _, _ = integrate(phi_h, x[1], None, N)
    rp = math.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6))
    sig = 2*W(phi) + t*(1 + c*phi); dsig = 2*W1(phi) + t*c
    return np.array([rp/rho - sig/6, s + dsig/2]), (rho, phi, s, rp)

def solve_bg(t, c, guess, N=20000):
    x = np.array(guess, float)
    for it in range(30):
        F, _ = shell_res(x, t, c, N)
        J = np.zeros((2, 2))
        for j in range(2):
            e = 1e-6
            xp = x.copy(); xp[j] += e; xm = x.copy(); xm[j] -= e
            J[:, j] = (shell_res(xp, t, c, N)[0] - shell_res(xm, t, c, N)[0])/(2*e)
        dx = np.linalg.solve(J, -F); lam = 1.0
        while lam > 1e-3 and np.linalg.norm(shell_res(x + lam*dx, t, c, N)[0]) > np.linalg.norm(F): lam /= 2
        x = x + lam*dx
        if np.linalg.norm(lam*dx) < 1e-12: break
    F, (rho, phi, s, rp) = shell_res(x, t, c, N)
    return dict(t=t, c=c, phi_h=-1 + 10**x[0], alpha=10**x[0], y_b=x[1], rho_b=rho, phi_b=phi, s_b=s, res=float(np.abs(F).max()))

def Bfun(mu2, bgd, N=20000, d2=0.0):
    y, rho, phi, s, psi, zeta = integrate(bgd['phi_h'], bgd['y_b'], mu2, N)
    rp = math.sqrt(1 + rho*rho*(s*s/12 - U(phi)/6)); spp = U1(phi) - 4*(rp/rho)*s
    sig2 = 2*W2(phi) + bgd['t']*d2
    B = 3*(np.asarray(mu2) + 4)*psi/(rho*rho*s) + (spp + sig2*s/2)*zeta
    chi_b = s*zeta
    alpha = 0.5 + np.sqrt(2.25 - np.asarray(mu2))
    norm = np.exp(-alpha*math.log(y/Y0))
    return B, chi_b, psi, norm

def find_roots(bgd, N, lo=-60.0, hi=2.2, m=400):
    grid = np.linspace(lo, hi, m); B = Bfun(grid, bgd, N)[0]
    roots = []
    for i in np.where(B[:-1]*B[1:] < 0)[0]:
        a, b = grid[i], grid[i + 1]
        for _ in range(5):
            g = np.linspace(a, b, 41); Bg = Bfun(g, bgd, N)[0]
            j = np.where(Bg[:-1]*Bg[1:] < 0)[0][0]
            # linear interpolation inside last bracket
            a, b = g[j], g[j + 1]
        Ba, Bb = Bfun(np.array([a, b]), bgd, N)[0]
        roots.append(a - Ba*(b - a)/(Bb - Ba))
    return roots, (grid, B)

if __name__ == "__main__":
    out = {}
    T0 = time.time()
    Ip = I_plus(); c = 2/Ip - 4/3
    print("own I_plus = %.13f   c_star = %.13f   (registered 1.0357712571567, 0.5975949350280)" % (Ip, c))
    ZE = c/2 + 3*c*c/8; mu2_eft = -4*(3*c*c - 4*c + 8)/(c*(3*c + 4))
    out.update(I_plus=Ip, c_star=c, mu2_eft=mu2_eft)
    # ---- registered shell
    t = 1e-3
    for N in (10000, 20000, 40000):
        bgd = solve_bg(t, c, (math.log10(0.2207*t**1.8), 8.2), N)
        roots, (grid, B) = find_roots(bgd, N)
        print("N=%6d  phi_h+1=%.8e  y_b=%.9f  rho_b=%.7f  phi_b=%.5e  s_b=%.8f  res=%.1e  | roots of B in [-60,2.2]: %s"
              % (N, bgd['alpha'], bgd['y_b'], bgd['rho_b'], bgd['phi_b'], bgd['s_b'], bgd['res'], ["%.9f" % r for r in roots]))
        out["registered_N%d" % N] = dict(bg=bgd, roots=roots)
    nsign = int((B[:-1]*B[1:] < 0).sum())
    print("   number of sign changes of B on 400-point grid [-60, 2.2]: %d" % nsign)
    mu0 = roots[0]
    # ---- residue of the brane-to-brane response  G = chi_b/(2B)
    d = 1e-4
    Bv, chib, psib, _ = Bfun(np.array([mu0 - d, mu0, mu0 + d]), bgd, 40000)
    dB = (Bv[2] - Bv[0])/(2*d)
    res_num = chib[1]/(2*dB); res_eft = -bgd['rho_b']**2/(2*Ip*ZE)
    print("residue of G=chi_b/j at the pole: numeric %.4f   EFT prediction -rho_b^2/(f Z_E) = %.4f   ratio %.6f  (negative = healthy norm)"
          % (res_num, res_eft, res_num/res_eft))
    out.update(residue_numeric=float(res_num), residue_eft=float(res_eft))
    # ---- argument principle with analytic normalisation
    for (x0, x1, yy) in [(-60.0, 2.0, 40.0), (-400.0, 2.0, 300.0)]:
        m = 500
        bottom = np.linspace(x0, x1, m) - 1j*yy; right = x1 + 1j*np.linspace(-yy, yy, m)
        top = np.linspace(x1, x0, m) + 1j*yy; left = x0 + 1j*np.linspace(yy, -yy, m)
        path = np.concatenate([bottom, right, top, left, bottom[:1]])
        Bv, chib, psib, norm = Bfun(path, bgd, 10000)
        for name, v in [("B*(y0/yb)^alpha", Bv*norm), ("psi_b*(y0/yb)^alpha", psib*norm), ("chi_b*(..)", chib*norm)]:
            ph = np.unwrap(np.angle(v)); wn = (ph[-1] - ph[0])/(2*math.pi)
            print("rect Re[%g,%g] |Im|<=%g : winding of %-22s = %+.5f   max phase jump between samples %.2f rad" % (x0, x1, yy, name, wn, np.abs(np.diff(ph)).max()))
            out["winding_%s_%g_%g" % (name, x0, yy)] = float(wn)
    # ---- t scan
    rows = []
    guess_yb = {1e-2: 6.16, 3e-3: 7.24, 1e-3: 8.23, 3e-4: 9.31, 1e-4: 10.30}
    for tt in (1e-2, 3e-3, 1e-3, 3e-4, 1e-4):
        b = solve_bg(tt, c, (math.log10(0.2207*tt**1.8), guess_yb[tt]), 20000)
        r, _ = find_roots(b, 20000, lo=-12.0, hi=2.2, m=72)
        rows.append((tt, r[0], len(r), b['rho_b'], b['res']))
        print("t=%.0e  mu2=%.9f  (#roots in [-12,2.2]: %d)  h=%.6e  bg residual %.1e   mu2 - mu2_EFT = %+.4e  ratio/t = %.4f"
              % (tt, r[0], len(r), 1/b['rho_b']**2, b['res'], r[0] - mu2_eft, (r[0] - mu2_eft)/tt))
    ts = np.array([r[0] for r in rows]); ms = np.array([r[1] for r in rows])
    A = np.vstack([np.ones_like(ts), ts, ts**2]).T
    coef = np.linalg.lstsq(A, ms, rcond=None)[0]
    print("quadratic fit mu2(t) = %.8f + %.4f t + %.2f t^2    (EFT limit %.8f; orchestrator slope 1.9243)" % (coef[0], coef[1], coef[2], mu2_eft))
    out.update(t_scan=[dict(t=r[0], mu2=r[1]) for r in rows], fit=[float(x) for x in coef])
    print("elapsed %.1f s" % (time.time() - T0))
    json.dump(out, open(os.path.join(HERE, "independent_spectrum_output.json"), "w"), indent=1, default=float)
