"""Auditor's independent high-precision shooting for the static +1 branch.

Independent of the audited mp_bvp.py: uses mpmath.odefun (mpmath's own Taylor integrator) on the FIRST-ORDER
system with the radial first integral ENFORCED (R_u = sqrt(1 + R^2 (eta_u^2/12 + V))), instead of the audited
second-order R equation with the first integral only monitored.  Cone start at small u0 from the linear AdS
solution eta = A C_14^(2)(cosh u)/680 plus the leading nonlinear correction; the shooting parameter is A = eta(0).
Shell: first zero of J1 = R_u/R - 3W - (3/2) delta (1 + c phi); A solved from J2 = eta_u + (9/2)(2W_phi + c delta) = 0.
usage: python3 v2_bvp_independent.py <dps> <u0> <c:'reg'|decimal> <delta1,delta2,...> <label>
Writes verify/runs/V2_BVP_<label>.json
"""
import sys, json, time
from pathlib import Path
import mpmath as mp

HERE = Path(__file__).resolve().parent
dps = int(sys.argv[1]); u0s = sys.argv[2]; cs = sys.argv[3]; dls = sys.argv[4].split(','); label = sys.argv[5]
mp.mp.dps = dps
c = (mp.mpf(2) / mp.mpf('1.0357712571566784') - mp.mpf(4) / 3) if cs == 'reg' else mp.mpf(cs)
u0 = mp.mpf(u0s)

def V(e): return 1 - 21 * e**2 - 25 * e**3 + mp.mpf(9) / 4 * e**4 + 6 * e**5 + e**6
def F(e): return 252 * e + 450 * e**2 - 54 * e**3 - 180 * e**4 - 36 * e**5
def W(e): return mp.mpf(1) / 3 + e**2 + e**3 / 3
def Wp(e): return 2 * e + e**2

# C_14^(2)(cosh u) = sum_j (j+1)(15-j) e^{(14-2j)u}  (standard Gegenbauer expansion; checked in v3)
def C14(u):
    return mp.fsum((j + 1) * (15 - j) * mp.exp((14 - 2 * j) * u) for j in range(15))
def C14p(u):
    return mp.fsum((j + 1) * (15 - j) * (14 - 2 * j) * mp.exp((14 - 2 * j) * u) for j in range(15))

def shoot(A, delta):
    # initial data at u0: linear solution + leading nonlinear correction (450 A^2 u^2 / 10); R from V(A)
    e0 = A * C14(u0) / 680 + 45 * A**2 * u0**2
    e0p = A * C14p(u0) / 680 + 90 * A**2 * u0
    kk = mp.sqrt(V(A))
    R0 = mp.sinh(kk * u0) / kk
    def rhs(u, y):
        R, e, ep = y
        Rp = mp.sqrt(1 + R**2 * (ep**2 / 12 + V(e)))
        return [Rp, ep, F(e) - 4 * Rp / R * ep]
    sol = mp.odefun(rhs, u0, [R0, e0, e0p])
    def J1(u):
        R, e, ep = sol(u)
        Rp = mp.sqrt(1 + R**2 * (ep**2 / 12 + V(e)))
        return Rp / R - 3 * W(e) - mp.mpf(3) / 2 * delta * (1 + c * (1 + e))
    # bracket the first zero of J1 on a coarse grid, then refine
    u = mp.mpf(1); h = mp.mpf('0.25'); g = J1(u)
    assert g > 0
    while True:
        g2 = J1(u + h)
        if g2 <= 0: break
        u += h; g = g2
        if u > 40: raise RuntimeError('no shell')
    ub = mp.findroot(J1, (u, u + h), solver='anderson')
    R, e, ep = sol(ub)
    J2 = ep + mp.mpf(9) / 2 * (2 * Wp(e) + c * delta)
    return ub, R, e, ep, J2

rows = []
for ds in dls:
    delta = mp.mpf(ds); t0 = time.time()
    # first-order guess for A: eta_h ~ -(111537/131072) c (1+c)^7 delta^8
    A0 = -mp.mpf(111537) / 131072 * c * (1 + c)**7 * delta**8
    if c == 0:
        A = mp.mpf(0)
    else:
        f = lambda a: shoot(a, delta)[4] / delta
        A = mp.findroot(f, (A0, A0 * (1 + mp.mpf('1e-2'))), solver='secant', tol=mp.mpf(10)**(-2 * dps + 8), maxsteps=40)
    ub, R, e, ep, J2 = shoot(A, delta)
    rows.append(dict(delta=ds, eta_h=mp.nstr(A, dps), u_b=mp.nstr(ub, dps), phi_b=mp.nstr(1 + e, dps), eta_b=mp.nstr(e, dps),
                     H2=mp.nstr(1 / (81 * R**2), dps), rho_b=mp.nstr(9 * R, dps), J2=mp.nstr(J2, 5), runtime_s=time.time() - t0))
    print(json.dumps(rows[-1]), flush=True)
(HERE / 'runs').mkdir(exist_ok=True)
(HERE / 'runs' / f'V2_BVP_{label}.json').write_text(json.dumps(dict(
    status='NUMERICAL (auditor-independent mpmath.odefun shooting, first integral enforced)', dps=dps, u0=u0s, c=mp.nstr(c, dps), rows=rows), indent=1) + '\n')
