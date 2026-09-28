#!/usr/bin/env python3
"""Shared library: archived shell trajectories and chi mode evolution on the shell.

Units: every rate is in units of H0 = 1/rho_b(initial static shell) of the archived runs
(H0 = 1/78.828... in the registered length unit L).  s = H0*tau (shell proper time).
Shell FRW metric: ds^2 = -dtau^2 + a(tau)^2 dx^2.  Canonical shell scalar chi with
    m_chi^2(s)/H0^2 = mu0^2 + G^2 (phi_b(s) - phi_star)^2 ,   G = gbar/H0 .
Mode equation for X_k = a^{3/2} chi_k (comoving kappa = k/(a_* H0), a(s_*) = 1 at the crossing):
    X'' + Omega^2 X = 0,  Omega^2 = kappa^2/a^2 + m^2 - (9/4) h^2 - (3/2) h',  h = H/H0.
Nothing here modifies the archived data; the archive is read read-only from the checkpoint.
"""
import hashlib, io, json, math, zipfile
from pathlib import Path
import numpy as np
from scipy.interpolate import make_interp_spline
from scipy.integrate import solve_ivp

CKPT = Path('/home/user/unified-theory-maldonado/new-files/latest-work/'
            'HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922')
ARCHIVE = CKPT / 'source_audit/inputs/CHAT14_COMPLETED_REFERENCE.zip'
BRANCH = CKPT / 'static_branch/PLUS_BRANCH_RESULTS.json'
PREFIX = 'HDBLAST_CHAT14_REGISTERED_DETUNING_ROLLOFF_20260922/runs/'
T_CAUSAL = 8.5          # coordinate time of the earliest far-taper signal at the shell, L*(1-taper) as recorded in the archive summaries

DELTA = 1e-3
C_REG = 2/1.0357712571566784 - 4/3

def W(p): return 1 - p + p**3/3
def sigma(p): return 2*W(p) + DELTA*(1 + C_REG*p)
def sigma_p(p): return 2*(p*p - 1) + DELTA*C_REG

def sha256_file(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()

def load_branch():
    d = json.loads(BRANCH.read_text())
    row = next(r for r in d['rows'] if abs(r['delta'] - DELTA) < 1e-14)['refined']
    return row

class Trajectory:
    """Quintic-spline representation of phi_b(s), h(s) and ln a(s) from one archived run."""
    def __init__(self, tag, t_cut=T_CAUSAL, k=5):
        with zipfile.ZipFile(ARCHIVE) as z:
            self.meta = json.loads(z.read(PREFIX + tag + '_summary.json'))
            rec = np.load(io.BytesIO(z.read(PREFIX + tag + '_timeseries.npz')))['rec']
        self.tag = tag
        self.rec_all = rec
        causal = self.meta['L']*(1 - self.meta['taper'])
        self.t_cut = min(t_cut, causal)
        rec = rec[rec[:, 0] < self.t_cut]
        s = rec[:, 9]
        if np.any(np.diff(s) <= 0): raise RuntimeError('proper time not monotone')
        self.rec = rec
        self.s = s; self.s_min = float(s[0]); self.s_max = float(s[-1])
        self.phi = rec[:, 1]; self.h = rec[:, 2]; self.v_stored = rec[:, 5]
        self.sp = make_interp_spline(s, np.column_stack([self.phi, self.h]), k=k)
        # ln a(s) = int h ds (a normalised to 1 at s_min); cross-check with int (1 + pa) dt
        self.lna_sp = self.sp.antiderivative(1)
        t = rec[:, 0]
        lna_t = np.concatenate([[0.0], np.cumsum(0.5*np.diff(t)*((1 + rec[1:, 4]) + (1 + rec[:-1, 4])))])
        lna_s = np.array([self.lna(x) for x in s])
        self.lna_crosscheck = float(np.max(np.abs(lna_t - lna_s)))
        self.phi_fit_residual = float(np.max(np.abs(self.sp(s)[:, 0] - self.phi)))

    def bg(self, s, nu=0):
        return self.sp(s, nu)            # columns (phi, h)

    def lna(self, s):
        return float(self.lna_sp(s)[1] - self.lna_sp(self.s_min)[1])

    def crossing(self, phistar):
        y = self.phi - phistar
        idx = np.flatnonzero(np.sign(y[:-1]) != np.sign(y[1:]))
        if not len(idx): return None
        i = int(idx[0])
        from scipy.optimize import brentq
        return brentq(lambda x: float(self.bg(x)[0]) - phistar, self.s[i], self.s[i + 1], xtol=1e-14)


class LinearToy:
    """Calibration background: phi = phistar + v (s - s0), h = 0, a = 1 (flat, linear crossing)."""
    def __init__(self, v, s0=0.0, span=60.0):
        self.v = v; self.s0 = s0; self.s_min = s0 - span; self.s_max = s0 + span; self.tag = 'linear_toy'
    def bg(self, s, nu=0):
        s = np.asarray(s, dtype=float)
        if nu == 0: out = np.stack([self.v*(s - self.s0), np.zeros_like(s)], -1)
        elif nu == 1: out = np.stack([self.v*np.ones_like(s), np.zeros_like(s)], -1)
        else: out = np.zeros(s.shape + (2,))
        return out
    def lna(self, s): return 0.0
    def crossing(self, phistar): return self.s0


def background_terms(traj, s, G, phistar, mu0, lna_star):
    """Return m2, dm2, ddm2, a^-2, h, h' and the H-part of Omega^2 with its two derivatives."""
    p0, h0 = traj.bg(s, 0); p1, h1 = traj.bg(s, 1); p2, h2 = traj.bg(s, 2); p3, h3 = traj.bg(s, 3)
    d = p0 - phistar
    m2 = mu0*mu0 + G*G*d*d
    dm2 = 2*G*G*d*p1
    ddm2 = 2*G*G*(p1*p1 + d*p2)
    ia2 = math.exp(-2*(traj.lna(s) - lna_star))
    hterm = -2.25*h0*h0 - 1.5*h1
    dh = -4.5*h0*h1 - 1.5*h2
    ddh = -4.5*(h1*h1 + h0*h2) - 1.5*h3
    return m2, dm2, ddm2, ia2, h0, h1, hterm, dh, ddh


def omega2_fast(traj, s, kap, G, phistar, mu0, lna_star):
    p0, h0 = traj.bg(s, 0); _, h1 = traj.bg(s, 1)
    d = p0 - phistar
    return kap*kap*math.exp(-2*(traj.lna(s) - lna_star)) + mu0*mu0 + G*G*d*d - 2.25*h0*h0 - 1.5*h1


def omega2_and_derivs(traj, s, kap, G, phistar, mu0, lna_star):
    m2, dm2, ddm2, ia2, h0, h1, hterm, dh, ddh = background_terms(traj, s, G, phistar, mu0, lna_star)
    k2 = kap*kap*ia2
    O2 = k2 + m2 + hterm
    dO2 = -2*h0*k2 + dm2 + dh
    ddO2 = (4*h0*h0 - 2*h1)*k2 + ddm2 + ddh
    return O2, dO2, ddO2


def wkb_freq(O2, dO2, ddO2, order):
    """Adiabatic frequency W and dW/ds.  order 0: W = Omega; order 2: W^2 = Omega^2 - (1/2)[Omega''/Omega - (3/2)(Omega'/Omega)^2]."""
    O = np.sqrt(O2)
    dO = dO2/(2*O)
    if order == 0:
        return O, dO
    ddO = (ddO2 - 2*dO*dO)/(2*O)
    W2 = O2 - 0.5*(ddO/O - 1.5*(dO/O)**2)
    if np.any(W2 <= 0):
        raise FloatingPointError('second-order adiabatic frequency squared <= 0 (adiabatic expansion invalid)')
    W_ = np.sqrt(W2)
    return W_, dO   # derivative of W approximated by dOmega/ds (difference is adiabatic order 3)


def adiabaticity(traj, s, G, phistar, mu0, lna_star, kap=0.0):
    O2, dO2, ddO2 = omega2_and_derivs(traj, s, kap, G, phistar, mu0, lna_star)
    if O2 <= 0: return np.inf
    O = math.sqrt(O2)
    e1 = abs(dO2)/(2*O2*O)
    e2 = abs(ddO2)/(2*O2*O2)
    return max(e1, math.sqrt(e2))


def find_window(traj, sstar, G, phistar, mu0, lna_star, tol, n_scan=4000):
    """Largest-extent search: from s_* outward, first point beyond which the kappa=0 adiabaticity stays <= tol."""
    out = {}
    for side, lim in (('left', traj.s_min), ('right', traj.s_max)):
        grid = np.linspace(sstar, lim, n_scan)
        A = np.array([adiabaticity(traj, x, G, phistar, mu0, lna_star) for x in grid])
        bad = np.flatnonzero(A > tol)
        j = int(bad[-1]) + 1 if len(bad) else 1
        if j >= n_scan:
            out[side] = (float(lim), True, float(A[-1]))    # truncated at data edge
        else:
            out[side] = (float(grid[j]), False, float(A[j]))
    return out


def evolve_modes(traj, kaps, G, phistar, mu0, s_a, s_b, lna_star, order=2, rtol=1e-10, atol=1e-14,
                 s_eval=None):
    """Adiabatic vacuum (given WKB order) at s_a, integrate to s_b, return |beta|^2 per kappa in the same-order basis.
    Also returns n in the zeroth-order basis for comparison."""
    kaps = np.asarray(kaps, dtype=float); N = len(kaps)
    O2, dO2, ddO2 = omega2_and_derivs(traj, s_a, kaps, G, phistar, mu0, lna_star)
    if np.any(O2 <= 0): raise RuntimeError('Omega^2 <= 0 at window start')
    Wa, dWa = wkb_freq(O2, dO2, ddO2, order)
    X0 = 1/np.sqrt(2*Wa) + 0j
    P0 = (-1j*Wa - dWa/(2*Wa))*X0
    y0 = np.concatenate([X0.real, X0.imag, P0.real, P0.imag])

    def rhs(s, y):
        O2s = omega2_fast(traj, s, kaps, G, phistar, mu0, lna_star)
        Xr, Xi, Pr, Pi = y[:N], y[N:2*N], y[2*N:3*N], y[3*N:]
        return np.concatenate([Pr, Pi, -O2s*Xr, -O2s*Xi])

    sol = solve_ivp(rhs, (s_a, s_b), y0, method='DOP853', rtol=rtol, atol=atol,
                    t_eval=s_eval, dense_output=False)
    if not sol.success: raise RuntimeError(sol.message)

    def occupation(s, y, order_):
        X = y[:N] + 1j*y[N:2*N]; P = y[2*N:3*N] + 1j*y[3*N:]
        O2s, dO2s, ddO2s = omega2_and_derivs(traj, s, kaps, G, phistar, mu0, lna_star)
        Wb, dWb = wkb_freq(O2s, dO2s, ddO2s, order_)
        f = 1/np.sqrt(2*Wb); fp = (-1j*Wb - dWb/(2*Wb))*f
        beta = -1j*(f*P - fp*X)
        alpha = 1j*(np.conj(f)*P - np.conj(fp)*X)
        wr = np.abs(alpha)**2 - np.abs(beta)**2    # must equal 1 (Wronskian)
        return np.abs(beta)**2, wr
    yb = sol.y[:, -1]
    nb, wr = occupation(s_b, yb, order)
    nb0, _ = occupation(s_b, yb, 0)
    res = {'n': nb, 'n_order0_basis': nb0, 'wronskian_dev': float(np.max(np.abs(wr - 1))), 'nfev': int(sol.nfev)}
    if s_eval is not None:
        # instantaneous zeroth-order-basis occupation (basis dependent inside the nonadiabatic interval; illustration only)
        res['history'] = [(float(sol.t[i]), occupation(sol.t[i], sol.y[:, i], 0)[0]) for i in range(len(sol.t))]
    return res


def gl_nodes(kmax, n):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5*kmax*(x + 1), 0.5*kmax*w


def densities(traj, s, kaps, wts, nk, G, phistar, mu0, lna_star):
    """Particle-description (late-time adiabatic, interference terms dropped) densities at shell time s, in H0 units:
    number n, energy rho, pressure p, <chi^2>_part and source j = G^2 (phi - phi*) <chi^2>."""
    p0 = float(traj.bg(s)[0])
    a = math.exp(traj.lna(s) - lna_star)
    m2 = mu0*mu0 + G*G*(p0 - phistar)**2
    k2p = (kaps/a)**2
    om = np.sqrt(k2p + m2)
    pref = 1/(2*math.pi**2*a**3)
    n = pref*np.sum(wts*kaps**2*nk)
    rho = pref*np.sum(wts*kaps**2*nk*om)
    p = pref*np.sum(wts*kaps**2*nk*k2p/(3*om))
    chi2 = pref*np.sum(wts*kaps**2*nk/om)
    j = G*G*(p0 - phistar)*chi2
    return dict(s=float(s), a=a, phi=p0, m=math.sqrt(m2), n=float(n), rho=float(rho), p=float(p), chi2=float(chi2), j=float(j))


def run_case(traj, G, phistar, mu0=0.0, tol=1e-3, order=2, rtol=1e-10, nk=48, kfac=40.0,
             history_kappa=None, n_hist=0):
    """Full production calculation for one (trajectory, G, phi*, mu0).  Returns a dict of numbers (H0 units)."""
    sstar = traj.crossing(phistar)
    if sstar is None: return {'crossing': False}
    lna_star = traj.lna(sstar)
    v_star = float(traj.bg(sstar, 1)[0]); h_star = float(traj.bg(sstar)[1])
    q = G*abs(v_star)
    win = find_window(traj, sstar, G, phistar, mu0, lna_star, tol)
    s_a, s_b = win['left'][0], win['right'][0]
    # looser criterion (A <= 0.05) marking where a particle description is approximately usable after the crossing
    s_part = find_window(traj, sstar, G, phistar, mu0, lna_star, 0.05)['right'][0]
    # momentum grid: Gaussian estimate, extended until the tail is negligible
    kmax = math.sqrt(kfac*max(q, 1.0)/math.pi)
    for attempt in range(6):
        kaps, wts = gl_nodes(kmax, nk)
        try:
            r = evolve_modes(traj, kaps, G, phistar, mu0, s_a, s_b, lna_star, order=order, rtol=rtol)
            order_used = order
        except FloatingPointError:
            r = evolve_modes(traj, kaps, G, phistar, mu0, s_a, s_b, lna_star, order=0, rtol=rtol)
            order_used = 0
        tail = float(np.max(r['n'][-3:])); peak = float(np.max(r['n']))
        if tail <= 1e-9*max(peak, 1e-300) or peak == 0: break
        kmax *= 1.5
    out = {'crossing': True, 'tag': traj.tag, 'G': G, 'phistar': phistar, 'mu0': mu0, 'tol': tol, 'order': order, 'order_used': order_used,
           'rtol': rtol, 'nk': nk, 's_star': sstar, 'v_star': v_star, 'h_star': h_star, 'q': q,
           'window': {'s_a': s_a, 's_b': s_b, 'left_truncated': win['left'][1], 'right_truncated': win['right'][1],
                      'A_left_edge': win['left'][2], 'A_right_edge': win['right'][2], 's_particle_A0.05': s_part},
           'kmax': kmax, 'tail_over_peak': tail/max(peak, 1e-300), 'wronskian_dev': r['wronskian_dev'], 'nfev': r['nfev'],
           'kappa': kaps.tolist(), 'weights': wts.tolist(), 'n_k': r['n'].tolist(), 'n_k_order0_basis': r['n_order0_basis'].tolist()}
    # analytic instant-preheating comparison (a_* = 1): n_k = exp(-pi (kappa^2 + mu0^2)/q)
    na = np.exp(-math.pi*(kaps**2 + mu0**2)/q)
    I_num = float(np.sum(wts*kaps**2*r['n'])); I_ana = float(np.sum(wts*kaps**2*na))
    out['analytic'] = {'N_comoving_numeric': I_num/(2*math.pi**2), 'N_comoving_analytic': I_ana/(2*math.pi**2),
                       'N_exact_formula': q**1.5/(8*math.pi**3)*math.exp(-math.pi*mu0**2/q),
                       'rel_diff_number': I_num/I_ana - 1 if I_ana > 0 else None,
                       'max_abs_diff_nk': float(np.max(np.abs(r['n'] - na)))}
    out['_kaps'] = kaps; out['_wts'] = wts; out['_n'] = r['n']; out['_lna_star'] = lna_star
    if history_kappa is not None and n_hist:
        se = np.linspace(s_a, s_b, n_hist)
        rh = evolve_modes(traj, np.array(history_kappa, float), G, phistar, mu0, s_a, s_b, lna_star, order=order,
                          rtol=rtol, s_eval=se)
        out['history'] = [{'s': t, 'n': nn.tolist()} for t, nn in rh['history']]
    return out
