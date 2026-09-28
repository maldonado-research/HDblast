#!/usr/bin/env python3
"""HDBLAST Chat 9 - verification of the moving-shell Friedmann theorem.

Run from this folder:   python3 verify_moving_shell.py
Writes:                 MOVING_SHELL_CHECKS.json and MOVING_SHELL_VERIFY_LOG.txt (copy of the printed log)
Only numpy + stdlib.  PART A is exact rational arithmetic (python fractions): sampled exact
controls of the algebra, not a symbolic proof (the proof is in MOVING_SHELL_FRIEDMANN_THEOREM.md).
PART B is floating point on the real registered background (RK4 + cubic Hermite interpolation).

Conventions (kappa_5 = 1, signature -++++):
  y-chart     ds^2 = dy^2 + rho(y)^2 [ -d eta^2 + cosh^2 eta dOmega_3^2 ],  a = rho cosh eta
  shell       (y(tau), eta(tau)),  ydot = sinh psi, rho etadot = cosh psi   (psi = rapidity w.r.t. the static leaves)
  n_out       = (cosh psi, sinh psi / rho), pointing AWAY from the kept bulk (kept side: y < y_shell, contains the cone)
  alpha = rho'/rho, beta = tanh(eta)/rho
  X = n_out.grad ln a = alpha cosh psi + beta sinh psi        [= K^theta_theta on the kept side]
  H = u.grad ln a     = alpha sinh psi + beta cosh psi
  junctions   X = (lam + rho_m)/6 ;  n_out.a = psidot + alpha cosh psi = (lam - 2 rho_m - 3 p_m)/6 ;
              n_out.grad phi = cosh psi * phi' = -(lam' + J)/2
"""
import json, math, os, random, sys, time
from fractions import Fraction as Fr
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "background"))
import hdblast_background as bg

RESULTS = {"exact": {}, "float": {}}
random.seed(20260916)

# ----------------------------------------------------------------------------------------------
# PART A : exact rational controls
# ----------------------------------------------------------------------------------------------
def rq(lo=-3, hi=3, den=97):
    while True:
        x = Fr(random.randint(lo*den, hi*den), random.randint(1, den))
        if x != 0: return x

def rpos(hi=3, den=97):
    return abs(rq(-hi, hi, den))

def runit():           # rational in (-1,1), nonzero
    while True:
        x = Fr(random.randint(-96, 96), 97)
        if x != 0: return x

def Wx(p):  return 1 - p + p**3/3
def W1x(p): return p*p - 1
def W2x(p): return 2*p
def Ux(p):  return Fr(1, 2)*W1x(p)**2 - Fr(2, 3)*Wx(p)**2
def U1x(p): return W1x(p)*W2x(p) - Fr(4, 3)*Wx(p)*W1x(p)

def exact_checks(N=400):
    fails = {k: 0 for k in ["E1_conformal_K2_minus_H2", "E1b_conformal_friedmann", "E2_constraint_substitution",
                            "E3_static_reduction_M462_7p1", "E3b_M462_7p2_expanded", "E4_ychart_friedmann",
                            "E5_chart_equivalence", "E6_codazzi_constraint_propagation", "E7_BPS_no_force",
                            "E8_junction_trace_algebra", "E9_maeda_wands_Lambda4", "E10_kinematic_gradphi",
                            "E11_sum_rule_integrand", "E12_rapidity_law"]}
    n1b = 0
    for _ in range(N):
        # E1: conformal chart, arbitrary (T,R,Tdot,Rdot,w)
        T, R, Td, Rd, w = rq(), rq(), rq(), rq(), rq()
        sig = R*R - T*T
        X = Td/R + 2*w*(R*Td - T*Rd)
        Y = Rd/R + 2*w*(R*Rd - T*Td)
        Nn = Td*Td - Rd*Rd
        if X*X - Y*Y != Nn*(1/(R*R) + 4*w + 4*sig*w*w): fails["E1_conformal_K2_minus_H2"] += 1
        if Nn > 0:
            n1b += 1
            Om2 = 1/Nn                       # normalisation Omega^2 (Td^2 - Rd^2) = 1
            a2 = Om2*R*R
            if Y*Y + 1/a2 != X*X - 4*w*(1 + sig*w)/Om2: fails["E1b_conformal_friedmann"] += 1
            # E10: (n phi)^2 - phidot^2 = 4 sigma Hs^2 / Omega^2   (eps-independent)
            Hs = rq()
            nphi = 2*Hs*(R*Td - T*Rd); phid = 2*Hs*(R*Rd - T*Td)
            if nphi**2 - phid**2 != 4*sig*Hs*Hs/Om2: fails["E10_kinematic_gradphi"] += 1
        # E2: constraint substitution
        sg, w, Hs, Om2 = rq(), rq(), rq(), rpos()
        w_sig = w*w - Hs*Hs/3
        Uval = -6*(4*w + sg*w_sig + 3*sg*w*w)/Om2          # registered constraint solved for U
        if -4*w*(1 + sg*w)/Om2 != Uval/6 - (4*sg*Hs*Hs/Om2)/12: fails["E2_constraint_substitution"] += 1
        # E3: static reduction to M462 (7.1), (7.2)
        ph, t, c = rq(), rq(), rq()
        st = 2*Wx(ph) + t*(1 + c*ph); dst = 2*W1x(ph) + t*c
        lhs = st*st/36 + Ux(ph)/6 - (dst/2)**2/12          # phi' = -sigma_t'/2
        r71 = t*(Wx(ph)*(1 + c*ph)/9 - c*W1x(ph)/12) + t*t*((1 + c*ph)**2/36 - c*c/48)
        r72 = t*(Fr(1, 9) + c/12 + (c - 1)/9*ph - 7*c/36*ph**2 + ph**3/27 + c*ph**4/27) \
            + t*t*(Fr(1, 36) - c*c/48 + c*ph/18 + c*c*ph*ph/36)
        if lhs != r71: fails["E3_static_reduction_M462_7p1"] += 1
        if r71 != r72: fails["E3b_M462_7p2_expanded"] += 1
        # E4: y-chart
        m = runit(); p = (1 + m*m)/(1 - m*m); v = 2*m/(1 - m*m)          # p^2 - v^2 = 1
        rho, rp, th = rpos(), rq(), runit()
        al, be = rp/rho, th/rho
        X = al*p + be*v; H = al*v + be*p
        inv_a2 = (1 - th*th)/(rho*rho)                                    # a = rho cosh eta
        if H*H + inv_a2 != X*X + (1 - rp*rp)/(rho*rho): fails["E4_ychart_friedmann"] += 1
        # E5: chart equivalence  rho = Omega b, b = sqrt(sigma), rho' = 1 + 2 sigma w
        b, Om, w = rpos(), rpos(), rq()
        sg = b*b; rho = Om*b; rp = 1 + 2*sg*w
        ok = ((1 - rp*rp)/(rho*rho) == -4*w*(1 + sg*w)/(Om*Om))
        k = runit(); ch = (1 + k*k)/(1 - k*k); sh = 2*k/(1 - k*k)
        T, R = b*sh, b*ch
        bd = v/Om; etd = p/rho                                           # bdot = b ydot/rho = ydot/Omega
        Td = bd*sh + b*ch*etd; Rd = bd*ch + b*sh*etd
        ok &= (Om*Om*(Td*Td - Rd*Rd) == 1)
        Xc = Td/R + 2*w*(R*Td - T*Rd); Hc = Rd/R + 2*w*(R*Rd - T*Td)
        al, be = rp/rho, (sh/ch)/rho
        ok &= (Xc == al*p + be*v) and (Hc == al*v + be*p)
        if not ok: fails["E5_chart_equivalence"] += 1
        # E6: Codazzi / constraint propagation (off-shell):  Cdot = -(phidot/3) * S - H * C
        rho, rp, s, th = rpos(), rq(), rq(), runit()
        rm, weos, J, lp = rpos(), rq(-1, 1), rq(), rq()
        al, be = rp/rho, th/rho
        X = al*p + be*v; H = al*v + be*p
        lam = rq(); Cc = X - (lam + rm)/6                                 # C NOT imposed: off-shell identity
        Q = (lam - 2*rm - 3*weos*rm)/6
        psid = Q - al*p                                                   # tau-tau junction
        al_y = -1/(rho*rho) - s*s/3                                       # bulk eq (rho'/rho)' = -1/rho^2 - phi'^2/3
        be_y = -al*be; be_eta = (1 - th*th)/rho
        rmd = -3*H*(1 + weos)*rm + J*v*s                                  # energy equation
        lamd = lp*v*s                                                     # lam = lam(phi)
        Cd = al_y*v*p + be_y*v*v + be_eta*(p/rho)*v + H*psid - (lamd + rmd)/6
        S = p*s + (lp + J)/2
        if Cd != -(v*s/3)*S - H*Cc: fails["E6_codazzi_constraint_propagation"] += 1
        # E7: BPS no-force:  Lambda_eff = 0 identically at t = 0
        if Wx(ph)**2/9 + Ux(ph)/6 - W1x(ph)**2/12 != 0: fails["E7_BPS_no_force"] += 1
        # E8: S^th_th - S/3 = (lam+rho)/3 ; S^tau_tau - S/3 = (lam - 2 rho - 3 p)/3  (with sign)
        lam, r_, p_ = rq(), rq(), rq()
        Stt, Sthth = -lam - r_, -lam + p_; Str = Stt + 3*Sthth
        # kept-side K (normal out of kept bulk) = +(1/2)(S^mu_nu - S delta/3)
        if (Sthth - Str/3)/2 != (lam + r_)/6 or (Stt - Str/3)/2 != (lam - 2*r_ - 3*p_)/6:
            fails["E8_junction_trace_algebra"] += 1
        # E9: Maeda-Wands (hep-th/0008188 eq 22) Lambda_4 = (1/2)[V + lam^2/6 - lam'^2/8] = 3[lam^2/36 + V/6 - lam'^2/48]
        V, lpp = rq(), rq()
        if 3*(lam*lam/36 + V/6 - lpp*lpp/48) != Fr(1, 2)*(V + lam*lam/6 - lpp*lpp/8): fails["E9_maeda_wands_Lambda4"] += 1
        # E11: sum rule  (rho^4 Z)'/rho^4 = 4 alpha Z + Z' = 1/rho^2 + 2 Z^2 - Fd^2/6
        #      Z = alpha - W/3, Fd = phi' + W_phi, orientation phi' = -W_phi, alpha = +W/3 on the BPS flow
        Zd, Fd = rq(), rq(); Wv, W1v = Wx(ph), W1x(ph)
        al = Zd + Wv/3; s_ = Fd - W1v
        inv_r2 = al*al - s_*s_/12 + Ux(ph)/6                  # Hamiltonian constraint solved for 1/rho^2
        Zp = (-inv_r2 - s_*s_/3) - W1v*s_/3                   # Z' = alpha' - W_phi phi'/3
        if 4*al*Zd + Zp != inv_r2 + 2*Zd*Zd - Fd*Fd/6: fails["E11_sum_rule_integrand"] += 1
        # E12: rapidity law  psidot = beta sinh(psi) - (rho_m + p_m)/2   on C = 0
        rho, rp, th = rpos(), rq(), runit(); rm, weos = rpos(), rq(-1, 1)
        al, be = rp/rho, th/rho; lam = 6*(al*p + be*v) - rm
        if (lam - 2*rm - 3*weos*rm)/6 - al*p != be*v - (1 + weos)*rm/2: fails["E12_rapidity_law"] += 1
    return dict(samples=N, samples_E1b=n1b, failures=fails, all_pass=all(v == 0 for v in fails.values()))

# ----------------------------------------------------------------------------------------------
# PART B : floating checks on the registered background
# ----------------------------------------------------------------------------------------------
class Background:
    def __init__(self, phi_h, y_end, dv=2e-4):
        _, out = bg.integrate(phi_h, y_end, dv=dv, record=True)
        self.y, self.rho, self.phi, self.s = out[:, 0], out[:, 1], out[:, 2], out[:, 3]
        self.h = self.y[1] - self.y[0]
        U = bg.U(self.phi)
        self.rp = np.sqrt(1 + self.rho**2*(self.s**2/12 - U/6))
        self.sp = bg.U1(self.phi) - 4*(self.rp/self.rho)*self.s
        # L(y) = ln sigma - 2 ln y,   L' = 2 (1/rho - 1/y),  L(0) = 0  (Omega(0) = 1)
        f = 1/self.rho - 1/self.y
        fp = -self.rp/self.rho**2 + 1/self.y**2
        u0 = bg.U(phi_h)
        cum = np.concatenate([[0.0], np.cumsum(self.h/2*(f[:-1] + f[1:]) + self.h**2/12*(fp[:-1] - fp[1:]))])
        self.L = 2*(u0*self.y[0]**2/72 + cum)
        self.Lp = 2*f

    def _herm(self, F, Fp, yq):
        x = (yq - self.y[0])/self.h
        i = int(min(max(math.floor(x), 0), len(self.y) - 2)); t = x - i; h = self.h
        h00 = (1 + 2*t)*(1 - t)**2; h10 = t*(1 - t)**2; h01 = t*t*(3 - 2*t); h11 = t*t*(t - 1)
        return h00*F[i] + h10*h*Fp[i] + h01*F[i+1] + h11*h*Fp[i+1]

    def at(self, yq):
        rho = self._herm(self.rho, self.rp, yq); phi = self._herm(self.phi, self.s, yq)
        s = self._herm(self.s, self.sp, yq)
        rp = math.sqrt(1 + rho*rho*(s*s/12 - bg.U(phi)/6))
        return rho, rp, phi, s

    def sigma(self, yq):
        return yq*yq*math.exp(self._herm(self.L, self.Lp, yq))


def shell_rhs(B, st, weos, Jfun, lam_fun=None):
    """state st = (y, eta, psi, rho_m, lam).  lam_fun=None -> reconstructed tension (lamdot from scalar junction)."""
    y, eta, psi, rm, lam = st
    rho, rp, phi, s = B.at(y)
    ch, sh = math.cosh(psi), math.sinh(psi)
    al, be = rp/rho, math.tanh(eta)/rho
    H = al*sh + be*ch
    if Jfun == "junction":                  # reconstructed coupling: J := -lam'(phi) - 2 cosh(psi) phi'
        J = -lam_fun(phi)[1] - 2*ch*s
    else:
        J = Jfun(rm, phi)
    if lam_fun is not None:
        lam = lam_fun(phi)[0]; lamd = 0.0
    else:
        lamd = -(2*ch*s + J)*sh*s          # lam' = -2 (n.grad phi) - J ,  phidot = sinh(psi) phi'
    psid = (lam - 2*rm - 3*weos*rm)/6 - al*ch
    rmd = -3*H*(1 + weos)*rm + J*sh*s
    return np.array([sh, ch/rho, psid, rmd, lamd])


def run_traj(B, st0, weos, Jfun, tau_end, dt, lam_fun=None, y_lo=0.5, y_hi=None):
    st = np.array(st0, float); rows = []; n = int(round(tau_end/dt)); tau = 0.0
    y_hi = y_hi if y_hi is not None else B.y[-1] - 0.05
    for k in range(n + 1):
        if not (y_lo < st[0] < y_hi) or abs(st[2]) > 2.5: break      # stop before leaving the tabulated bulk / ultra-relativistic plunge
        rows.append((tau, *st))
        f = lambda z: shell_rhs(B, z, weos, Jfun, lam_fun)
        k1 = f(st); k2 = f(st + dt/2*k1); k3 = f(st + dt/2*k2); k4 = f(st + dt*k3)
        st = st + dt/6*(k1 + 2*k2 + 2*k3 + k4); tau += dt
    return np.array(rows)


def analyse(B, rows, weos, Jfun, lam_fun=None, conformal=True):
    tau = rows[:, 0]; dt = tau[1] - tau[0]
    n = len(tau); lna = np.zeros(n); Hf = np.zeros(n); C = np.zeros(n); S = np.zeros(n)
    rhsF = np.zeros(n); inva2 = np.zeros(n); phid = np.zeros(n); conf = np.zeros((n, 3)); lamv = np.zeros(n)
    Jrec = np.zeros(n); hyp = np.zeros((n, 3)); phib = np.zeros(n); rap = np.zeros(n)
    for i, (_, y, eta, psi, rm, lam) in enumerate(rows):
        rho, rp, phi, s = B.at(y); ch, sh = math.cosh(psi), math.sinh(psi)
        if lam_fun is not None: lam, lamp = lam_fun(phi)
        else: lamp = -2*ch*s - Jfun(rm, phi)
        Jv = (-lamp - 2*ch*s) if Jfun == "junction" else Jfun(rm, phi)
        Jrec[i] = Jv
        lamv[i] = lam
        al, be = rp/rho, math.tanh(eta)/rho
        X = al*ch + be*sh; Hf[i] = al*sh + be*ch
        a = rho*math.cosh(eta); lna[i] = math.log(a); inva2[i] = 1/a**2
        C[i] = X - (lam + rm)/6
        S[i] = ch*s + (lamp + Jv)/2
        rhsF[i] = (lam + rm)**2/36 + bg.U(phi)/6 - s*s/12
        phid[i] = sh*s; phib[i] = phi
        rap[i] = be*sh - (1 + weos)*rm/2 - shell_rhs(B, rows[i, 1:], weos, Jfun, lam_fun)[2]   # rapidity law (valid on C=0)
        if conformal:
            sg = B.sigma(y); b = math.sqrt(sg); Om = rho/b
            w = (rp - 1)/(2*sg)
            T, R = b*math.sinh(eta), b*math.cosh(eta)
            bd = sh/Om; etd = ch/rho
            Td = bd*math.sinh(eta) + b*math.cosh(eta)*etd; Rd = bd*math.cosh(eta) + b*math.sinh(eta)*etd
            Xc = Td/R + 2*w*(R*Td - T*Rd); Hc = Rd/R + 2*w*(R*Rd - T*Td)
            conf[i] = (Om*Om*(Td*Td - Rd*Rd) - 1, Xc - X, Hc - Hf[i])
            T0 = T - R*Rd/Td; b2 = R*R - (T - T0)**2
            hyp[i] = (T0, b2, lam/6 - (1 + w*(sg + b2 + T0*T0))/(Om*math.sqrt(abs(b2))))
    # independent Hubble rate: 4th-order centred finite difference of ln a(tau)
    Hfd = np.full(n, np.nan)
    Hfd[2:-2] = (-lna[4:] + 8*lna[3:-1] - 8*lna[1:-3] + lna[:-4])/(12*dt)
    Cfd = np.full(n, np.nan)
    Cfd[2:-2] = (-C[4:] + 8*C[3:-1] - 8*C[1:-3] + C[:-4])/(12*dt)
    sl = slice(2, n - 2)
    return dict(
        n_steps=int(n), tau_end=float(tau[-1]),
        y_range=[float(rows[:, 1].min()), float(rows[:, 1].max())],
        eta_end=float(rows[-1, 2]), psi_range=[float(rows[:, 3].min()), float(rows[:, 3].max())],
        a_ratio=float(math.exp(lna[-1] - lna[0])),
        lam_range=[float(lamv.min()), float(lamv.max())],
        H_range=[float(Hf.min()), float(Hf.max())],
        max_abs_theta_theta_constraint_C=float(np.abs(C).max()),
        max_abs_scalar_junction_S=float(np.abs(S).max()),
        max_abs_H_formula_minus_H_finite_difference=float(np.abs(Hf[sl] - Hfd[sl]).max()),
        max_abs_friedmann_residual_using_Hfd=float(np.abs(Hfd[sl]**2 + inva2[sl] - rhsF[sl]).max()),
        max_abs_friedmann_residual_using_Hformula=float(np.abs(Hf**2 + inva2 - rhsF).max()),
        friedmann_scale=float(np.abs(rhsF).max()),
        max_abs_Cdot_plus_phidot_S_over_3_plus_HC=float(np.abs(Cfd[sl] + phid[sl]*S[sl]/3 + Hf[sl]*C[sl]).max()),
        max_abs_Cdot_finite_difference=float(np.abs(Cfd[sl]).max()),
        J_range=[float(Jrec.min()), float(Jrec.max())],
        max_abs_rapidity_law_residual=float(np.abs(rap).max()),
        phi_b_monotone=bool(np.all(np.diff(phib) > 0) or np.all(np.diff(phib) < 0)),
        vacuum_shell_hyperboloid=(dict(T0_first=float(hyp[0, 0]), T0_spread=float(np.ptp(hyp[:, 0])),
                                       b2_first=float(hyp[0, 1]), b2_spread=float(np.ptp(hyp[:, 1])),
                                       max_abs_closed_form_tension_residual=float(np.abs(hyp[:, 2]).max()))
                                  if (conformal and np.all(rows[:, 4] == 0.0)) else None),
        conformal_chart_max_abs=dict(normalisation=float(np.abs(conf[:, 0]).max()),
                                     K_theta_theta=float(np.abs(conf[:, 1]).max()),
                                     Hubble=float(np.abs(conf[:, 2]).max())) if conformal else None)


def float_checks():
    out = {}
    reg = json.load(open(os.path.join(HERE, "..", "background", "REGISTERED_SHELL_FLOAT_SOLUTION.json")))
    t, c, phi_h, y_b = reg["t"], reg["c"], reg["phi_h"], reg["v_b"]
    Y_END = 9.6
    B = Background(phi_h, Y_END)
    sig_t = lambda ph: (2*bg.W(ph) + t*(1 + c*ph), 2*bg.W1(ph) + t*c)

    # --- B1: static registered shell
    rho, rp, phi, s = B.at(y_b); lam, lamp = sig_t(phi)
    h = 1/rho**2
    h71 = t*(bg.W(phi)*(1 + c*phi)/9 - c*bg.W1(phi)/12) + t*t*((1 + c*phi)**2/36 - c*c/48)
    out["B1_static_shell"] = dict(y_b=y_b, rho_b=rho, phi_b=phi, s_b=s, h=h,
        israel_residual=rp/rho - lam/6, scalar_residual=s + lamp/2,
        theorem_rhs=lam**2/36 + bg.U(phi)/6 - s*s/12, theorem_rhs_minus_h=lam**2/36 + bg.U(phi)/6 - s*s/12 - h,
        M462_7p1=h71, M462_7p1_minus_h=h71 - h, sigma_conformal_at_shell=B.sigma(y_b),
        Omega_at_shell=rho/math.sqrt(B.sigma(y_b)))

    # --- B2: conformal-chart background identities on a grid (finite differences in sigma)
    idx = np.arange(5000, len(B.y) - 5000, 2500)
    sg = B.y**2*np.exp(B.L); w = (B.rp - 1)/(2*sg); Hs = B.s*B.rho/(2*sg); Om2 = B.rho**2/sg
    dwdy = np.gradient(w, B.h, edge_order=2); wsig = dwdy*B.rho/(2*sg)          # dy/dsigma = rho/(2 sigma)
    r1 = wsig - (w*w - Hs*Hs/3)
    r2 = 4*w + sg*wsig + 3*sg*w*w + Om2*bg.U(B.phi)/6
    r3 = -4*w*(1 + sg*w)/Om2 - (bg.U(B.phi)/6 - B.s**2/12)                      # algebraic, no FD
    dOm = np.gradient(0.5*np.log(Om2), B.h, edge_order=2)*B.rho/(2*sg) - w       # w = Omega_sigma/Omega
    out["B2_conformal_background_identities"] = dict(
        grid_points=int(len(idx)), y_min=float(B.y[idx[0]]), y_max=float(B.y[idx[-1]]),
        max_abs_w_sigma_minus_w2_plus_Hs2_over_3=float(np.abs(r1[idx]).max()), scale_w2=float(np.abs(w[idx]**2).max()),
        max_abs_registered_constraint=float(np.abs(r2[idx]).max()), scale_4w=float(np.abs(4*w[idx]).max()),
        max_abs_F_conformal_minus_F_ychart=float(np.abs(r3[idx]).max()),
        max_abs_w_minus_dlnOmega_dsigma=float(np.abs(dOm[idx]).max()))

    # --- B3..: moving trajectories, reconstructed tension (exact solutions of all three junctions)
    noJ = lambda rm, ph: 0.0
    rho, rp, phi, s = B.at(y_b)
    def start(y0, eta0, psi0, rm0):
        r_, rp_, _, _ = B.at(y0)
        X = (rp_/r_)*math.cosh(psi0) + math.tanh(eta0)/r_*math.sinh(psi0)
        return (y0, eta0, psi0, rm0, 6*X - rm0)
    cases = {
        "B3_vacuum_shell_outward_kick":      dict(st=start(y_b, 0.0, 0.30, 0.0), w=0.0, J=noJ, tau=3.0),
        "B4_vacuum_shell_inward_kick":       dict(st=start(y_b, 0.2, -0.40, 0.0), w=0.0, J=noJ, tau=6.0),
        "B5_radiation_shell":                dict(st=start(6.0, 0.1, 0.25, 0.8), w=1/3, J=noJ, tau=4.0),
        "B6_dust_shell_with_coupling_J":     dict(st=start(7.0, -0.1, 0.20, 0.5), w=0.0, J=(lambda rm, ph: 0.3*rm), tau=4.0),
        "B7_stiff_matter_deep_bulk":         dict(st=start(3.0, 0.3, 0.50, 1.5), w=1.0, J=noJ, tau=3.0),
    }
    for name, cs in cases.items():
        rows = run_traj(B, cs["st"], cs["w"], cs["J"], cs["tau"], 1e-3)
        res = analyse(B, rows, cs["w"], cs["J"])
        # step-halving control on the constraint
        rows2 = run_traj(B, cs["st"], cs["w"], cs["J"], cs["tau"], 2e-3)
        res["max_abs_C_at_double_step"] = float(np.abs(analyse(B, rows2, cs["w"], cs["J"], conformal=False)
                                                       ["max_abs_theta_theta_constraint_C"]))
        res["initial_state_y_eta_psi_rhom_lam"] = [float(x) for x in cs["st"]]
        res["w_eos"] = cs["w"]
        # reconstructed tension sample lam(phi_b)
        k = np.linspace(0, len(rows) - 1, 7).astype(int)
        res["reconstructed_tension_samples_phi_lam"] = [[float(B.at(rows[i, 1])[2]), float(rows[i, 5])] for i in k]
        out[name] = res

    # --- B8: over-determination demo: registered sigma_t(phi) kept FIXED, shell kicked, frozen bulk.
    # Evolve with the tau-tau junction; the theta-theta constraint then drifts at rate -(phidot/3) S.
    rows = run_traj(B, (y_b, 0.0, 0.05, 0.0, 0.0), 0.0, noJ, 3.0, 1e-3, lam_fun=sig_t)
    res = analyse(B, rows, 0.0, noJ, lam_fun=sig_t, conformal=False)
    out["B8_fixed_registered_tension_kicked_shell_frozen_bulk"] = dict(
        note="C != 0 and S != 0 here: a generic fixed sigma_t(phi) is NOT compatible with motion in the frozen bulk; "
             "the drift law Cdot = -(phidot/3) S - H C is verified by finite differences.",
        psi0=0.05, max_abs_C=res["max_abs_theta_theta_constraint_C"], max_abs_S=res["max_abs_scalar_junction_S"],
        max_abs_Cdot_plus_phidot_S_over_3_plus_HC=res["max_abs_Cdot_plus_phidot_S_over_3_plus_HC"],
        max_abs_Cdot_finite_difference=res["max_abs_Cdot_finite_difference"], y_range=res["y_range"])

    # --- B11: reconstructed COUPLING family: registered sigma_t(phi) kept fixed, dust, J fixed by the scalar junction
    r_, rp_, ph_, s_ = B.at(y_b); psi0 = 0.30
    rm0 = 6*(rp_/r_)*math.cosh(psi0) - sig_t(ph_)[0]
    rows = run_traj(B, (y_b, 0.0, psi0, rm0, 0.0), 0.0, "junction", 3.0, 1e-3, lam_fun=sig_t)
    res = analyse(B, rows, 0.0, "junction", lam_fun=sig_t, conformal=True)
    res["rho_m_initial"] = rm0; res["rho_m_range"] = [float(rows[:, 4].min()), float(rows[:, 4].max())]
    out["B11_registered_tension_dust_shell_reconstructed_J"] = res

    # --- B10: exact sum rule  (rho^4 Z)' = rho^2 + rho^4 (2 Z^2 - Fd^2/6),  Z = rho'/rho - W/3, Fd = phi' + W_phi
    _, o = bg.integrate(phi_h, y_b, dv=2e-4, record=True)
    yy, r_, p_, s_ = o[:, 0], o[:, 1], o[:, 2], o[:, 3]
    rp_ = np.sqrt(1 + r_**2*(s_**2/12 - bg.U(p_)/6))
    Z = rp_/r_ - bg.W(p_)/3; Fd = s_ + bg.W1(p_)
    def simpson(f):
        m = len(f) - 1
        if m % 2: return simpson(f[:-1]) + (yy[-1] - yy[-2])*(f[-1] + f[-2])/2
        hh = yy[1] - yy[0]; return hh/3*(f[0] + f[-1] + 4*f[1:-1:2].sum() + 2*f[2:-1:2].sum())
    rb = r_[-1]
    I_kept = simpson((r_/rb)**2) + yy[0]**3/(3*rb**2)
    Q2 = simpson((r_/rb)**4*(2*Z**2 - Fd**2/6)) + 2*yy[0]**3/(3*rb**4)
    dlam = t*(1 + c*p_[-1])
    out["B10_planck_mass_sum_rule"] = dict(
        statement="delta_lambda(phi_b)/6 = h*I_kept + int (rho/rho_b)^4 (2 Z^2 - Fd^2/6) dy,  I_kept = int_0^{y_b} (rho/rho_b)^2 dy",
        Z_b=float(Z[-1]), delta_lambda_over_6=dlam/6, h_times_I_kept=float(I_kept/rb**2), quadratic_defect_integral=float(Q2),
        residual=float(dlam/6 - I_kept/rb**2 - Q2), I_kept=float(I_kept), I_plus_flat_BPS=reg["I_plus"],
        first_order_prediction_h=float(dlam/(6*I_kept)), h=float(1/rb**2),
        relative_error_of_first_order_prediction=float(dlam/(6*I_kept)*rb**2 - 1))

    # --- B9: Lambda_eff table along the registered background
    table = []
    for yq in [0.02, 0.5, 1, 2, 3, 4, 5, 6, 7, 7.5, 8, y_b, 8.5, 9, 9.5]:
        rho, rp, phi, s = B.at(yq); lam, lamp = sig_t(phi)
        F = bg.U(phi)/6 - s*s/12
        table.append(dict(y=yq, phi=phi, rho=rho, sigma_conf=B.sigma(yq), sigma_t=lam, tension_term=lam*lam/36,
                          F_bulk=F, F_bulk_check=(1 - rp*rp)/rho**2, Lambda_eff=3*(lam*lam/36 + F),
                          Lambda_eff_over_3h=(lam*lam/36 + F)*reg["rho_b"]**2,
                          israel_mismatch_sigma_t_over6_minus_alpha=lam/6 - rp/rho,
                          scalar_mismatch_phiprime_plus_half_dsigma=s + lamp/2,
                          defect_s_minus_Wphi=s - bg.W1(phi),
                          eightpiG_friedmann=lam/6))
    out["B9_Lambda_eff_table"] = table
    out["B9_reference"] = dict(three_h=3*reg["h"], eightpiG_zero_mode_1_over_2Iplus=1/(2*reg["I_plus"]),
                               eightpiG_friedmann_at_shell=sig_t(reg["phi_b"])[0]/6,
                               h_lin_coefficient_1_over_6Iplus=1/(6*reg["I_plus"]),
                               naive_coefficient_sigma_over_18=sig_t(reg["phi_b"])[0]/18)
    return out


class _Tee:
    def __init__(self, path):
        self.o = sys.stdout
        try: self.f = open(path, "w")
        except OSError: self.f = None          # log copy is optional
    def write(self, x):
        self.o.write(x)
        if self.f: self.f.write(x)
    def flush(self):
        self.o.flush()
        if self.f: self.f.flush()

if __name__ == "__main__":
    sys.stdout = _Tee(os.path.join(HERE, "MOVING_SHELL_VERIFY_LOG.txt"))
    t0 = time.time()
    RESULTS["exact"] = exact_checks()
    print("PART A exact rational controls:", json.dumps(RESULTS["exact"], indent=1))
    RESULTS["float"] = float_checks()
    for k, v in RESULTS["float"].items():
        if k == "B9_Lambda_eff_table":
            print("\nB9 Lambda_eff table")
            print("%8s %12s %10s %10s %13s %13s %13s %12s %12s" % ("y", "phi", "rho", "sigma_t", "sig_t^2/36", "F_bulk",
                  "Lambda_eff", "Isr.mism", "scal.mism"))
            for r in v:
                print("%8.4f %12.4e %10.4f %10.6f %13.6e %13.6e %13.6e %12.3e %12.3e" % (r["y"], r["phi"], r["rho"],
                      r["sigma_t"], r["tension_term"], r["F_bulk"], r["Lambda_eff"],
                      r["israel_mismatch_sigma_t_over6_minus_alpha"], r["scalar_mismatch_phiprime_plus_half_dsigma"]))
        else:
            print("\n" + k, json.dumps(v, indent=1))
    RESULTS["runtime_s"] = time.time() - t0
    json.dump(RESULTS, open(os.path.join(HERE, "MOVING_SHELL_CHECKS.json"), "w"), indent=1)
    print("\nruntime %.1f s ; wrote MOVING_SHELL_CHECKS.json" % RESULTS["runtime_s"])
