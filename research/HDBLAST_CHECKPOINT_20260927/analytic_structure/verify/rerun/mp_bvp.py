"""High-precision (mpmath) solution of the full nonlinear static +1-branch boundary-value problem.

Independent of series_expansion.py.  Variables: u = y/9, R = rho/9, eta = phi - 1.
   R'' = R (V(eta) - eta'^2/4),          V = -U/(6k^2) = 1 - 21 e^2 - 25 e^3 + 9/4 e^4 + 6 e^5 + e^6
   eta'' = -4 (R'/R) eta' + F(eta),        F = U_eta/k^2 = 252 e + 450 e^2 - 54 e^3 - 180 e^4 - 36 e^5
Regular cone at u = 0 (R = u + ..., eta = eta_h + ...), started from a Frobenius power series.
Integration: fixed-step Taylor series method of order ORDER (automatic recurrences).
Shell: first u where J1: R'/R - 3W(phi) - (3/2) delta (1 + c phi) = 0; unknown eta_h solved from
J2: eta' + (9/2) sigma'(phi) = 0.  H^2 = 1/rho_b^2 = 1/(81 R_b^2).
Also reports the first-integral residual R'^2 - 1 - R^2 (eta'^2/12 + V) (not enforced).
"""
import mpmath as mp

def potV(e):
    e2 = e * e
    return 1 - 21 * e2 - 25 * e2 * e + mp.mpf(9) / 4 * e2 * e2 + 6 * e2 * e2 * e + e2 * e2 * e2

def Wphi(e):  # W at phi = 1 + e
    return mp.mpf(1) / 3 + e * e + e ** 3 / 3

class Series:
    """helpers for truncated power series stored as python lists of mpf"""
    @staticmethod
    def conv(a, b, n):
        return mp.fsum(a[i] * b[n - i] for i in range(n + 1))

def cone_series(eta_h, order):
    r = [mp.mpf(0)] * (order + 2); e = [mp.mpf(0)] * (order + 2)
    r[1] = mp.mpf(1); e[0] = eta_h
    # power series of powers of e, V, F, eta'
    for n in range(2, order + 1):
        # r_n from [R (V - eta'^2/4)]_{n-2}
        m = n - 2
        ser_e = e[:m + 1]
        pw = [ser_e]
        for k in range(5):
            pw.append([Series.conv(pw[-1], ser_e, i) for i in range(m + 1)])
        # pw[j] = e^(j+1)
        Vs = [(1 if i == 0 else 0) - 21 * pw[1][i] - 25 * pw[2][i] + mp.mpf(9) / 4 * pw[3][i] + 6 * pw[4][i] + pw[5][i] for i in range(m + 1)]
        ed = [(i + 1) * e[i + 1] for i in range(m + 1)]
        ed2 = [Series.conv(ed, ed, i) for i in range(m + 1)]
        G = [Vs[i] - ed2[i] / 4 for i in range(m + 1)]
        r[n] = Series.conv(r, G, m) / (n * (n - 1))
        if n % 2 == 0:
            # e_n from R eta'' + 4 R' eta' = R F  at u^{n-1}
            Fs = [252 * pw[0][i] + 450 * pw[1][i] - 54 * pw[2][i] - 180 * pw[3][i] - 36 * pw[4][i] for i in range(m + 1)]
            rhs = mp.fsum(r[a] * Fs[n - 1 - a] for a in range(1, n) if n - 1 - a <= m)
            other = mp.fsum(r[a] * e[n + 1 - a] * (n + 1 - a) * ((n - a) + 4 * a) for a in range(3, n + 1) if n + 1 - a >= 0)
            e[n] = (rhs - other) / (n * (n + 3))
    return r[:order + 1], e[:order + 1]

def taylor_coeffs(state, order):
    R0, R1, E0, E1 = state
    r = [mp.mpf(0)] * (order + 3); e = [mp.mpf(0)] * (order + 3)
    r[0], r[1], e[0], e[1] = R0, R1, E0, E1
    pw = [[] for _ in range(6)]          # e^1..e^6 coefficients
    q = []; ed2 = []
    for n in range(order + 1):
        pw[0].append(e[n])
        for k in range(1, 6):
            pw[k].append(Series.conv(pw[k - 1], e, n))
        Vn = (1 if n == 0 else 0) - 21 * pw[1][n] - 25 * pw[2][n] + mp.mpf(9) / 4 * pw[3][n] + 6 * pw[4][n] + pw[5][n]
        Fn = 252 * pw[0][n] + 450 * pw[1][n] - 54 * pw[2][n] - 180 * pw[3][n] - 36 * pw[4][n]
        edn = [(i + 1) * e[i + 1] for i in range(n + 1)]
        ed2.append(Series.conv(edn, edn, n))
        Rp_n = (n + 1) * r[n + 1]
        qn = (Rp_n - mp.fsum(r[k] * q[n - k] for k in range(1, n + 1))) / r[0]
        q.append(qn)
        Gn_list = None
        # [R (V - ed2/4)]_n
        Vs_n = Vn - ed2[n] / 4
        if n == 0:
            G = [Vs_n]
        else:
            G.append(Vs_n)
        r[n + 2] = Series.conv(r, G, n) / ((n + 2) * (n + 1))
        qe = mp.fsum(q[k] * edn[n - k] for k in range(n + 1))
        e[n + 2] = (-4 * qe + Fn) / ((n + 2) * (n + 1))
    return r[:order + 2], e[:order + 2]

def peval(c, t):
    s = mp.mpf(0)
    for x in reversed(c):
        s = s * t + x
    return s

def pderiv(c):
    return [(i) * c[i] for i in range(1, len(c))]

class Problem:
    def __init__(self, delta, cc, dps=50, order=50, h=mp.mpf('0.125'), u_cone=mp.mpf('0.25'), cone_order=90):
        self.delta = mp.mpf(delta); self.c = mp.mpf(cc)
        self.dps, self.order, self.h, self.u_cone, self.cone_order = dps, order, mp.mpf(h), mp.mpf(u_cone), cone_order

    def J1(self, R, Rp, e):
        return Rp / R - 3 * Wphi(e) - mp.mpf(3) / 2 * self.delta * (1 + self.c * (1 + e))

    def integrate(self, eta_h, full=False):
        r, e = cone_series(eta_h, self.cone_order)
        u = self.u_cone
        state = (peval(r, u), peval(pderiv(r), u), peval(e, u), peval(pderiv(e), u))
        g_prev = self.J1(state[0], state[1], state[2])
        if g_prev <= 0:
            raise RuntimeError('shell before cone start')
        while True:
            rc, ec = taylor_coeffs(state, self.order)
            h = self.h
            new = (peval(rc, h), peval(pderiv(rc), h), peval(ec, h), peval(pderiv(ec), h))
            g_new = self.J1(new[0], new[1], new[2])
            if g_new <= 0:
                drc, dec = pderiv(rc), pderiv(ec)
                f = lambda t: self.J1(peval(rc, t), peval(drc, t), peval(ec, t))
                t = mp.findroot(f, (mp.mpf(0), h), solver='anderson')
                if not (0 <= t <= h):
                    t = mp.findroot(f, (mp.mpf(0), h), solver='bisect')
                shell = (peval(rc, t), peval(drc, t), peval(ec, t), peval(dec, t))
                return u + t, shell, (r, e)
            state = new; u += h; g_prev = g_new
            if u > 60:
                raise RuntimeError('no shell found')

    def J2(self, shell):
        R, Rp, e, ep = shell
        return ep + mp.mpf(9) / 2 * (2 * e * (2 + e) + self.delta * self.c)

    def solve(self, eta_h_guess):
        with mp.workdps(self.dps):
            if self.c == 0:
                ub, shell, _ = self.integrate(mp.mpf(0))
                eta_h = mp.mpf(0)
            else:
                f = lambda x: self.J2(self.integrate(x)[1]) / self.delta
                x0 = mp.mpf(eta_h_guess)
                eta_h = mp.findroot(f, (x0, x0 * (1 + mp.mpf('1e-3'))), solver='secant', tol=mp.mpf(10) ** (-2 * self.dps + 10), maxsteps=60)
                ub, shell, _ = self.integrate(eta_h)
            R, Rp, e, ep = shell
            V = potV(e)
            fi = Rp ** 2 - 1 - R ** 2 * (ep ** 2 / 12 + V)
            return dict(delta=self.delta, c=self.c, eta_h=eta_h, u_b=ub, y_b=9 * ub, R_b=R, rho_b=9 * R,
                        eta_b=e, phi_b=1 + e, phi_y_b=ep / 9, H2=1 / (81 * R * R),
                        J1=self.J1(R, Rp, e), J2=self.J2(shell), first_integral_residual=fi,
                        first_integral_scale=1 + Rp ** 2 + R ** 2 * (ep ** 2 / 12 + abs(V)))
