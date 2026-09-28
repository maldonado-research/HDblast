"""Constraint-control laboratory for the balanced-disturbance evolutions.

Floating-point numerical experiments only; not an interval certificate.

The continuum model, the perturbation variables (a,b,f,pa,pb,pf), the
background subtraction and the shell junction data are those of the frozen
registered solver (copied verbatim into src/).  This module rebuilds the
spatial operators so that the following can be varied independently:

  order     centred interior differencing order p in {2,4,6} on the analytic
            stretched grid, with a Hermite brane ghost extension of degree p+1
            (values at p+1 nodes plus the prescribed junction slope);
  kappa     outgoing-characteristic constraint damping: a source
            s = -kappa*(H+2M)/(6*(A_t+A_z)) added to the B_tt equation only;
  ko        Kreiss-Oliger dissipation in the velocity equations only;
  xtol      absolute tolerance of the coordinate inversion used to sample the
            initial profiles (the original wrapper uses 2e-14).

With order=4, kappa=0, ko=0 and xtol=2e-14 the operators are identical to the
frozen solver's (checked by check_operator_identity in controls.py).
"""
from pathlib import Path
import importlib, math, sys
import numpy as np
from scipy import sparse
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'src'))
R = importlib.import_module('registered_solver')
E = importlib.import_module('evolve_balanced')      # copied wrapper (imports src copies)

D1W = {2: [-1/2, 0, 1/2],
       4: [1/12, -8/12, 0, 8/12, -1/12],
       6: [-1/60, 9/60, -45/60, 0, 45/60, -9/60, 1/60]}
D2W = {2: [1, -2, 1],
       4: [-1/12, 16/12, -30/12, 16/12, -1/12],
       6: [2/180, -27/180, 270/180, -490/180, 270/180, -27/180, 2/180]}


def hermite_ghost_weights(order):
    """Ghost values at x=1..r from a degree-(order+1) polynomial fitted to the
    values at x=0,-1,...,-order and the first derivative at x=0.
    Returns (r, order+2) array: weights on [u0,u_-1,...,u_-order, u'(0)]."""
    deg = order + 1
    nv = order + 1
    V = np.zeros((deg + 1, deg + 1))
    xx = -np.arange(nv, dtype=float)
    V[:nv] = xx[:, None] ** np.arange(deg + 1)[None, :]
    V[nv, 1] = 1.0
    r = order // 2
    X = np.array([(k + 1.0) ** np.arange(deg + 1) for k in range(r)])
    return X @ np.linalg.inv(V)


def build_operators(zp, zpp, order):
    n = len(zp); r = order // 2; gh = hermite_ghost_weights(order); nv = order + 1
    # Extension matrix built in COO form: assigning sparse.eye(n) into a LIL matrix
    # (as the frozen solver does) creates a dense n x n temporary (~10 GB at h=5e-5).
    rows = list(range(r, n + r)); cols = list(range(n)); vals = [1.0] * n
    gv = np.zeros(n + 2 * r)
    for k in range(r):
        for j in range(nv):
            rows.append(n + r + k); cols.append(n - 1 - j); vals.append(gh[k, j])
        gv[n + r + k] = gh[k, nv] * zp[-1]
    E_ = sparse.coo_matrix((vals, (rows, cols)), shape=(n + 2 * r, n)).tocsr()
    D1x = sum(c * E_[k:k + n] for k, c in enumerate(D1W[order]) if c != 0)
    D2x = sum(c * E_[k:k + n] for k, c in enumerate(D2W[order]) if c != 0)
    g1 = sum(c * gv[k:k + n] for k, c in enumerate(D1W[order]))
    g2 = sum(c * gv[k:k + n] for k, c in enumerate(D2W[order]))
    D1 = sparse.diags(1 / zp) @ D1x
    D2 = sparse.diags(1 / zp ** 2) @ D2x - sparse.diags(zpp / zp ** 3) @ D1x
    return D1.tocsr(), D2.tocsr(), g1 / zp, g2 / zp ** 2 - g1 * zpp / zp ** 3, gh


def make_profile_sampler(xtol):
    """Return a replacement for E.profile with a chosen inverse-map tolerance.
    The dense radial solutions are cached per epsilon within the process."""
    cache = {}
    rtol_inv = 4 * np.finfo(float).eps
    def profile(epsilon, rtol=8e-14):
        key = (epsilon, rtol)
        if key not in cache:
            cache[key] = E.B.run(epsilon, rtol=rtol, transport=True)
        sol, meta = cache[key]; zb = sol.y[4, -1]
        def sample(z):
            targets = np.asarray(z) + zb
            if np.min(targets) < sol.y[4, 0]:
                raise ValueError('Grid extends beyond cone start')
            ix = np.clip(np.searchsorted(sol.y[4], targets) - 1, 0, len(sol.t) - 2)
            ys = np.array([brentq(lambda y: sol.sol(y)[4] - tg, sol.t[i], sol.t[i + 1],
                                  xtol=xtol, rtol=rtol_inv) for tg, i in zip(targets, ix)])
            return sol.sol(ys)
        return sample, meta
    return profile


def frozen_matrices_lowmem(zp, zpp):
    """Drop-in replacement for the frozen R.matrices without the dense temporary;
    equality with R.matrices is checked in controls.py."""
    D1, D2, g1, g2, _ = build_operators(zp, zpp, 4)
    return D1, D2, g1, g2, R.hermite_ghosts()


_SAMPLERS = {}
def initial_data(epsilon, hmin, L, stretch, xtol=2e-14):
    """Balanced-seed initial data exactly as the copied wrapper builds it,
    optionally with a different inverse-map tolerance."""
    if xtol not in _SAMPLERS:
        _SAMPLERS[xtol] = make_profile_sampler(xtol)
    old, oldm = E.profile, R.matrices
    try:
        E.profile = _SAMPLERS[xtol]   # scipy's default brentq rtol is 4*eps, as here
        R.matrices = frozen_matrices_lowmem
        s, state, meta = E.setup(epsilon, hmin, L, stretch)
    finally:
        E.profile, R.matrices = old, oldm
    return s, state, meta


class Lab:
    def __init__(self, s, order=4, kappa=0.0, ko=0.0, damp_form='cplus', damp_off=None):
        # Copy the background from the wrapper-initialised frozen solver object.
        self.z, self.zp, self.zpp, self.n = s.z, s.zp, s.zpp, s.n
        self.L, self.hmin = s.L, s.hmin
        self.rho, self.phi, self.hc, self.phz = s.rho, s.phi, s.hc, s.phz
        self.rho2, self.rb, self.pb, self.pot = s.rho2, s.rb, s.pb, s.pot
        self.s0, self.s10, self.s20 = s.s0, s.s10, s.s20
        self.order, self.kappa, self.ko, self.damp_form = order, kappa, ko, damp_form
        self.r = order // 2
        # kappa(z): optional C2 switch-off of the damping near the shell,
        # kappa=0 for z>-damp_off[0], full kappa for z<-damp_off[1].
        if damp_off is None:
            self.kz = kappa * np.ones_like(s.z)
        else:
            x = np.clip((-s.z - damp_off[0]) / (damp_off[1] - damp_off[0]), 0, 1)
            self.kz = kappa * x ** 3 * (10 - 15 * x + 6 * x * x)
        self.nfix = max(2, self.r)
        self.D1, self.D2, self.g1, self.g2, self.gh = build_operators(self.zp, self.zpp, order)
        self.bc = s.bc            # identical nonlinear junction data
        self.scale = 1 + 6 * self.hc ** 2 + self.phz ** 2
        self.mask = (self.z > -.8 * self.L) & (np.arange(self.n) > 5)
        self.core = self.z > -.2
        # composite-trapezoid weights in z for L2 norms on the mask
        zz = self.z
        w = np.zeros(self.n); dz = np.diff(zz)
        w[:-1] += dz / 2; w[1:] += dz / 2
        self.wz = w

    # ---------------------------------------------------------------- pieces
    def _potential_terms(self, b, f):
        du = np.zeros_like(f); du1 = np.zeros_like(f)
        for j in range(6, 0, -1): du = (du + self.pot[j] / math.factorial(j)) * f
        for j in range(5, 0, -1): du1 = (du1 + self.pot[j + 1] / math.factorial(j)) * f
        eb = np.expm1(2 * b)
        su = self.rho2 * (eb * (self.pot[0] + du) + du)
        su1 = self.rho2 * (eb * (self.pot[1] + du1) + du1)
        return su, su1

    def _ko(self, v, g):
        """KO 2r-th difference (r=3 fixed, as in the frozen solver) with Hermite
        ghosts of the scheme; zero outer ghosts."""
        gh3 = R.hermite_ghosts()         # the frozen solver's three quintic ghosts
        gh = gh3[:, :5] @ v[-1:-6:-1] + gh3[:, 5] * self.zp[-1] * g
        ex = np.concatenate((np.zeros(3, dtype=v.dtype), v, gh))
        return self.ko * sum(c * ex[i:i + self.n] for i, c in enumerate([1, -6, 15, -20, 15, -6, 1])) / (64 * self.zp)

    def constraints(self, state, independent_edge=True):
        a, b, f, pa, pb, pf = state
        ga, gf, gpa, gpf = self.bc(b[-1], f[-1], pb[-1], pf[-1])
        d1 = self.D1 @ state[:3].T + self.g1[:, None] * np.array([ga, ga, gf])
        az, bz, fz = d1.T
        azz = self.D2 @ a + self.g2 * ga
        paz = self.D1 @ pa + self.g1 * gpa
        if independent_edge:
            one = np.array([25 / 12, -4, 3, -4 / 3, 1 / 4]) / self.zp[-1]
            paz = paz.copy(); paz[-1] = np.dot(one, pa[-1:-6:-1])
        su, _ = self._potential_terms(b, f)
        M = -3 * paz - 3 * (1 + pa) * (az - bz) + 3 * (self.hc + az) * pb - pf * (self.phz + fz)
        H = (-2 * su + 12 * pa + 6 * pa * pa + 6 * (1 + pa) * pb - 18 * self.hc * az + 6 * self.hc * bz
             - 12 * az * az + 6 * az * bz - 6 * azz - pf * pf - 2 * self.phz * fz - fz * fz)
        return H, M, az

    def rhs(self, state):
        a, b, f, pa, pb, pf = state
        ga, gf, gpa, gpf = self.bc(b[-1], f[-1], pb[-1], pf[-1])
        d1 = self.D1 @ state[:3].T + self.g1[:, None] * np.array([ga, ga, gf])
        d2 = self.D2 @ state[:3].T + self.g2[:, None] * np.array([ga, ga, gf])
        az, bz, fz = d1.T; azz, bzz, fzz = d2.T
        su, su1 = self._potential_terms(b, f)
        kin = 2 * pa + pa * pa; grad = 2 * self.hc * az + az * az
        out = np.empty_like(state); out[:3] = state[3:]
        out[3] = azz - 3 * kin + 3 * grad + (2 / 3) * su
        out[4] = bzz + 3 * kin - 3 * grad - .5 * pf * pf + self.phz * fz + .5 * fz * fz - su / 3
        out[5] = fzz - 3 * (1 + pa) * pf + 3 * (self.hc * fz + az * self.phz + az * fz) - su1
        if self.kappa:
            paz = self.D1 @ pa + self.g1 * gpa     # interior operator incl. exact shell slope
            M = -3 * paz - 3 * (1 + pa) * (az - bz) + 3 * (self.hc + az) * pb - pf * (self.phz + fz)
            H = (-2 * su + 12 * pa + 6 * pa * pa + 6 * (1 + pa) * pb - 18 * self.hc * az + 6 * self.hc * bz
                 - 12 * az * az + 6 * az * bz - 6 * azz - pf * pf - 2 * self.phz * fz - fz * fz)
            if self.damp_form == 'cplus':
                out[4] += -self.kz * (H + 2 * M) / (6 * (1 + pa + self.hc + az))
            elif self.damp_form == 'honly':   # deliberately naive variant (expected unstable where A_z>A_t)
                out[4] += -self.kz * H / (6 * (1 + pa))
            else:
                raise ValueError(self.damp_form)
        if self.ko:
            out[3] += self._ko(pa, gpa); out[4] += self._ko(pb, gpa); out[5] += self._ko(pf, gpf)
        out[:, :self.nfix] = 0
        return out

    def diagnostics(self, state, t):
        H, M, az = self.constraints(state)
        a = state[0]
        m, c = self.mask, self.core
        hn, mn = np.abs(H) / self.scale, np.abs(M) / self.scale
        wt = (self.rho / self.rb) ** 3 * np.exp(np.clip(3 * (t + a), -700, 700))
        cp, cm = wt * (H + 2 * M), wt * (H - 2 * M)
        ih = np.flatnonzero(m)[np.argmax(hn[m])]
        return dict(time=t, H_max=float(hn[m].max()), M_max=float(mn[m].max()),
                    H_L2=float(np.sqrt(np.sum(self.wz[m] * hn[m] ** 2))),
                    M_L2=float(np.sqrt(np.sum(self.wz[m] * mn[m] ** 2))),
                    H_peak_z=float(self.z[ih]),
                    core_H=float(hn[c].max()), core_M=float(mn[c].max()),
                    Cplus_w_max=float(np.abs(cp[m]).max()), Cminus_w_max=float(np.abs(cm[m]).max()),
                    shell_phi_dev=float(state[2, -1]), max_abs_f=float(np.abs(state[2]).max()),
                    finite=bool(np.isfinite(state).all())), H, M


def evolve(lab, state, tf=1.0, snap_times=(0.25, 0.5, 0.75, 1.0), record_dt=0.025, cfl=0.4,
           guard=0.5, progress=None):
    dt0 = cfl * lab.zp.min(); unit = round(1 / dt0)
    # choose dt so that every snapshot time (multiples of .25) is hit exactly
    unit = 4 * math.ceil(unit / 4); dt = 1.0 / unit
    steps = round(tf * unit)
    every = max(1, round(record_dt * unit))
    snap_steps = {round(ts * unit): ts for ts in snap_times if ts <= tf + 1e-12}
    v = state.copy(); rows = []; snaps = {}; reason = 'tf'
    for i in range(steps + 1):
        t = i * dt
        if i % every == 0 or i == steps or i in snap_steps:
            row, H, M = lab.diagnostics(v, t)
            if i % every == 0 or i == steps: rows.append(row)
            if i in snap_steps:
                snaps[snap_steps[i]] = dict(state=v.copy(), H=H.copy(), M=M.copy(), row=row)
            if progress and i % (every * 8) == 0: progress(row)
            if not row['finite'] or row['max_abs_f'] > guard or row['H_max'] > 1e6:
                reason = 'guard'; break
        if i == steps: break
        k1 = lab.rhs(v); k2 = lab.rhs(v + .5 * dt * k1); k3 = lab.rhs(v + .5 * dt * k2); k4 = lab.rhs(v + dt * k3)
        v += dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    return dict(dt=dt, steps=steps, reason=reason, rows=rows), snaps, v


# ---------------------------------------------------------------------------
# Remedy: discrete-constraint projection of the initial data.
def project_initial_hamiltonian(lab, state, z_left=None, tol=1e-13, maxit=12):
    """Replace a=b on z>=z_left by the solution of the *discrete* Hamiltonian
    constraint H_i=0 (the evolution's own operators and junction ghosts),
    holding f and all velocities fixed (pa=pb=pf=0 at t=0) and a fixed for
    z<z_left (taper region).  The momentum constraint stays exactly zero
    because a=b and the velocities vanish.  Newton iteration with an analytic
    sparse Jacobian.  Returns the new state and an iteration log."""
    if np.any(state[3:] != 0) or np.any(state[0] != state[1]):
        raise ValueError('projection assumes a=b and zero velocity perturbations')
    z_left = -.85 * lab.L if z_left is None else z_left
    idx = np.flatnonzero(lab.z >= z_left)
    if idx[0] < lab.nfix + lab.r: raise ValueError('projection window touches frozen nodes')
    v = state.copy(); log = []
    f = v[2]
    for it in range(maxit):
        a = v[0]
        H, M, az = lab.constraints(v, independent_edge=False)
        res = H[idx]
        log.append(dict(iteration=it, max_abs_H_window=float(np.abs(res).max()),
                        max_abs_M=float(np.abs(M[lab.mask]).max())))
        if np.abs(res).max() < tol * (1 + np.abs(lab.scale[idx]).max()) and it > 0: break
        # Jacobian of H wrt a (with b=a) on the window
        du = sum(lab.pot[j] * f ** j / math.factorial(j) for j in range(1, 7))
        dsu = 2 * lab.rho2 * np.exp(2 * a) * (lab.pot[0] + du)
        eb = np.exp(a[-1]); pbv = lab.pb; ff = f[-1]
        ds = lab.s10 * ff + 2 * pbv * ff * ff + (2 / 3) * ff ** 3
        dga = lab.rb * (lab.s0 * eb + eb * ds) / 6
        n = lab.n
        last = sparse.csr_matrix(([1.], ([0], [n - 1])), shape=(1, n))
        D1f = lab.D1 + sparse.csr_matrix((lab.g1 * dga)[:, None]) @ last
        D2f = lab.D2 + sparse.csr_matrix((lab.g2 * dga)[:, None]) @ last
        # the scalar junction slope also depends on b_b: gf=-rb(s10*expm1(b)+e^b*ds1)/2
        ds1 = lab.s20 * ff + 2 * ff * ff
        dgf = -lab.rb * (lab.s10 * eb + eb * ds1) / 2
        fz = lab.D1 @ f + lab.g1 * lab.bc(v[1, -1], ff, 0., 0.)[1]
        Jf = sparse.csr_matrix(((-2 * lab.phz - 2 * fz) * lab.g1 * dgf)[:, None]) @ last
        J = sparse.diags(-2 * dsu) + sparse.diags(-12 * lab.hc - 12 * az) @ D1f - 6 * D2f + Jf
        J = J.tocsr()[idx][:, idx]
        from scipy.sparse.linalg import spsolve
        # NATURAL ordering keeps the banded structure (COLAMD produced very large fill-in
        # at the finest grids; one order-6 setup reached ~10 GB before this change).
        da = spsolve(J.tocsc(), -res, permc_spec="NATURAL")
        v[0, idx] += da; v[1, idx] += da
    return v, log
