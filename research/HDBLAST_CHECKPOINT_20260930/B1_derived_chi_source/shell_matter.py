#!/usr/bin/env python3
"""B1 shell-matter sectors for the 5D roll-off solver (evolve_b1.py) and the reduced fixed-background estimate.

All shell-matter quantities enter the junctions only through the totals (model units, kappa5 = 1 convention of the solver)
    rho = kappa5^2 rho_tot,  p = kappa5^2 p_tot,  J0 = kappa5^2 j (velocity-independent part),  plus Y v (friction closure),
with (22 Sept matter extension, MATTER_EXTENSION_AND_RESIDUAL_VACUUM.md sections 2-3)
    n.A = (sigma + rho)/6,  n.B = (sigma - (2 rho + 3 p))/6,  n.phi = -(sigma' + J0 + Y v)/2,
    d rho_m/dtau + 3 H (rho_m + p_m) = j v  (Ward identity, summed over all shell matter).

Sectors
  Friction  (A1 closure, unchanged): state [R];  dR/dtau + 4 H R = Y v^2.
  Toy       (friction-equivalent, calibration only): state [R, Xt]; radiation-like fluid Xt fed by j = Y v and decaying into R
            at constant rate Gt (H0 units):  dXt/dtau + 4 H Xt = Y v^2 - Gt Xt,  dR/dtau + 4 H R = Gt Xt.
            The total R + Xt obeys the A1 ledger exactly; the split exercises the multi-component plumbing.
  Chi       (DERIVED source, MODEL CHANGE B1): a gas of chi quanta, m_chi^2 = mu0^2 + G^2 (phi_b - phi*)^2 (H0 units),
            produced at each crossing of phi* by instant preheating (optionally with the O(1/q) curved/expanding correction
            factor 1 + D/q of HDBLAST_CHECKPOINT_20260927/preheating CORRECTION_FIT.json), free thereafter, decaying with the
            rest-frame rate Gamma = y^2 m/(8 pi) (time-dilated per mode, Gamma m/omega) into massless radiation R.
            Per momentum node i of each cohort (physical momentum kappa_i and number density N_i at the crossing, H0 units):
               n_i = N_i a_rel^-3,  p_i = kappa_i / a_rel,  a_rel = exp(A_b - A_*),  omega_i = sqrt(p_i^2 + m^2)
               rho = sum n_i omega_i,  p = sum n_i p_i^2/(3 omega_i),  j = G^2 (phi - phi*) sum n_i/omega_i = sum n_i d omega_i/d phi
               dN_i/dtau = -Gamma (m/omega_i) N_i,   dR/dtau + 4 H R = kappa5^2 Gamma m n.
            Conversion to model units: kappa5^2 rho / H0 = b rho_hat, b = kappa5^2 H0^3 = (H0/M5)^3.
Units inside the chi sector: H0 units (hats); rb = 1/H0 in model units.
"""
import math
import numpy as np
from scipy.special import roots_genlaguerre

# O(1/q) correction coefficients (closed forms quoted in preheating/README.md section 5 and CORRECTION_FIT.json):
#   D = c_hh h0^2 + c_h1 h1 + c_AA (A/v)^2 + c_B (B/v) + c_hA h0 (A/v),  A = phi'', B = phi''' (s = H0 tau), v = phi'
C_HH = 9*math.pi/4 + 15/(2*math.pi)
C_H1 = 3*math.pi/2 - 15/(8*math.pi)
C_AA = 45/(32*math.pi) - math.pi/8
C_B = math.pi/8 - 45/(96*math.pi)
C_HA = 45/(8*math.pi)

def D_correction(v, A, B, h0, h1):
    a = A/v; bb = B/v
    return C_HH*h0*h0 + C_H1*h1 + C_AA*a*a + C_B*bb + C_HA*h0*a

def local_screen(v, A, h0, h1, eps=0.1):
    """minimum G of the 22 Sept local approximation screen (matter report section 7)"""
    v = abs(v)
    return max(h0*h0/(eps*eps*v), A*A/(eps*eps*v**3), abs(h1)/(eps*v))

def instant_nodes(q, nk=24, mu0=0.0, factor=1.0):
    """momentum nodes for n_k = factor exp(-pi (k^2 + mu0^2)/q):  number density per node (H0^3), sum = q^1.5/(8 pi^3) factor e^{-pi mu0^2/q}.
    Generalised Gauss-Laguerre (alpha = 1/2) in u = pi k^2/q:  int d^3k/(2 pi)^3 n_k f = (q/pi)^1.5/(4 pi^2) int u^1/2 e^-u f du."""
    u, w = roots_genlaguerre(nk, 0.5)
    kap = np.sqrt(q*u/math.pi)
    N = w*(q/math.pi)**1.5/(4*math.pi**2)*factor*math.exp(-math.pi*mu0*mu0/q)
    return kap, N


class Friction:
    kind = 'friction'
    pure_radiation = True
    def __init__(self, Y): self.Y = Y
    def init_state(self): return np.zeros(1)
    def totals(self, m, Ab, phib, T=None, phiT=None, eB=None):
        return m[0], m[0]/3, 0.0
    def rhs(self, m, Ab, phib, AT, phiT, eB, T=None):
        return np.array([-4*AT*m[0] + self.Y*phiT**2/eB])
    def totals_T(self, m, dm, Ab, phib, AT, phiT, eB, T=None):
        return dm[0], dm[0]/3, 0.0
    def extras(self, m, Ab, phib, AT, phiT, eB, T=None): return dict(X=0.0, Pchi=0.0, J0=0.0)


class Toy:
    """friction-equivalent toy through a separate radiation-like component that decays into R (calibration K2)"""
    kind = 'toy'
    pure_radiation = True
    def __init__(self, Y, Gt_hat, rb):
        self.Y = Y; self.Gt = Gt_hat/rb            # model-unit rate
    def init_state(self): return np.zeros(2)       # [R, Xt]
    def totals(self, m, Ab, phib, T=None, phiT=None, eB=None):
        r = m[0] + m[1]
        return r, r/3, 0.0
    def rhs(self, m, Ab, phib, AT, phiT, eB, T=None):
        dec = eB*self.Gt*m[1]
        return np.array([-4*AT*m[0] + dec, -4*AT*m[1] + self.Y*phiT**2/eB - dec])
    def totals_T(self, m, dm, Ab, phib, AT, phiT, eB, T=None):
        r = dm[0] + dm[1]
        return r, r/3, 0.0
    def extras(self, m, Ab, phib, AT, phiT, eB, T=None): return dict(X=float(m[1]), Pchi=float(m[1]/3), J0=0.0)


class ChiGas:
    """derived chi source (see module docstring).  State m = [R, (N_1..N_nk, u) per cohort].
    Production (method choice of B1, energy-consistent, field-space ramp): at a crossing of phi* a cohort with final comoving
    numbers N_fin (instant preheating x (1 + D/q)) is opened.  Its progress u in [0, 1] advances only while the scalar moves
    away from phi* (du/dtau = dir v / dphi_r, dir = sign of the crossing velocity, dphi_r = ramp_c |v*|/sqrt(q), i.e. the
    field excursion covered in ramp_c/sqrt(q) at the crossing velocity; the non-adiabatic region is |phi - phi*| ~ sqrt(|v|/G)).
    N = N_fin f(u), f(u) = 10u^3 - 15u^4 + 6u^5.  The production power P = sum dN/dtau|prod omega a^-3 = j_prod v is paid by the
    scalar through the bounded force j_prod = dir sum N_fin f'(u) omega a^-3 / dphi_r (zero when the scalar does not move
    outward), so rho_dot + 3H(rho+p) = (j + j_prod) v - Q holds at all times (B1_SYMBOLIC_CHECKS C4 with f = f(u(tau))).
    mode='instant' inserts N_fin at once (unpaid energy jump; control only)."""
    kind = 'chi'
    pure_radiation = False
    def __init__(self, G, phistar, y, b, rb, mu0=0.0, nk=24, apply_D=True, D_max=0.3, sign_j=1.0, ramp_c=3.0, mode='ramp', vs_frac=0.05):
        self.G, self.ps, self.y, self.b, self.rb, self.mu0, self.nk = G, phistar, y, b, rb, mu0, nk
        self.apply_D, self.D_max = apply_D, D_max
        self.sign_j = sign_j        # +1 physical; -1 = wrong-sign control (only for the ledger control, never for results)
        self.ramp_c, self.mode = ramp_c, mode
        self.vs_frac = vs_frac      # stall regularisation: v_s = vs_frac |v*|
        self.cohorts = []           # list of dict(Astar, kap, i0, iu, Nfin, dphi, dir, vs)
        self.events = []
        self.Y = 0.0
    def init_state(self): return np.zeros(1)
    @staticmethod
    def _f(u):
        if u <= 0: return 0.0, 0.0
        if u >= 1: return 1.0, 0.0
        return u**3*(10 - 15*u + 6*u*u), 30*u*u*(1 - u)**2
    # ---- per-node kinematics
    def _nodes(self, m, Ab, phib):
        if not self.cohorts: return None
        ns, ps = [], []
        for c in self.cohorts:
            ar = math.exp(Ab - c['Astar'])
            N = m[c['i0']:c['i0'] + self.nk]
            ns.append(N/ar**3); ps.append(c['kap']/ar)
        n = np.concatenate(ns); p = np.concatenate(ps)
        d = phib - self.ps
        m2 = self.mu0**2 + self.G**2*d*d
        om = np.sqrt(p*p + m2)
        return n, p, om, d, math.sqrt(m2)
    @staticmethod
    def _s(x):
        """C^1 switch: 0 for x <= 0, 3x^2 - 2x^3 on (0, 1), 1 for x >= 1"""
        if x <= 0: return 0.0
        if x >= 1: return 1.0
        return x*x*(3 - 2*x)
    def _active(self, m, vh):
        """cohorts in production: list of (cohort, f'(u), s) with u < 1; s = switch(dir vh / v_s) regularises the stop of production
        when the scalar stalls (du/dtau = dir v s/dphi_r, so j_prod = P/v = dir sum N_fin f' omega a^-3 s/dphi_r stays bounded and C^1)"""
        out = []
        if vh is None: return out
        for c in self.cohorts:
            if c['dphi'] <= 0: continue
            u = m[c['iu']]
            if u >= 1: continue
            sf = self._s(c['dir']*vh/c['vs'])
            if sf <= 0: continue
            f, fp = self._f(u)
            if fp > 0: out.append((c, fp, sf))
        return out
    def jprod_hat(self, m, Ab, phib, vh):
        """bounded production force (H0^4 units): dir sum N_fin f'(u) s omega a^-3 / dphi_r over active cohorts"""
        d = phib - self.ps; m2 = self.mu0**2 + self.G**2*d*d; J = 0.0
        for c, fp, sf in self._active(m, vh):
            ar = math.exp(Ab - c['Astar']); p = c['kap']/ar
            J += c['dir']*fp*sf*float(np.sum(c['Nfin']*np.sqrt(p*p + m2)))/ar**3/c['dphi']
        return J
    def hat_totals(self, m, Ab, phib, T=None, phiT=None, eB=None):
        """(rho, p, j_adiabatic, n, m_chi, Q, j_prod) in H0 units"""
        mc0 = math.sqrt(self.mu0**2 + self.G**2*(phib - self.ps)**2)
        nd = self._nodes(m, Ab, phib)
        if nd is None: return 0.0, 0.0, 0.0, 0.0, mc0, 0.0, 0.0
        n, p, om, d, mc = nd
        rho = float(np.sum(n*om)); pr = float(np.sum(n*p*p/om)/3); j = self.sign_j*self.G**2*d*float(np.sum(n/om))
        ntot = float(np.sum(n)); Q = self.y**2*mc*mc/(8*math.pi)*ntot
        vh = (self.rb*phiT/eB) if (phiT is not None and eB is not None) else None
        jp = self.sign_j*self.jprod_hat(m, Ab, phib, vh)
        return rho, pr, j, ntot, mc, Q, jp
    def totals(self, m, Ab, phib, T=None, phiT=None, eB=None):
        rho, pr, j, _, _, _, jp = self.hat_totals(m, Ab, phib, T, phiT, eB)
        f = self.b/self.rb
        return m[0] + f*rho, m[0]/3 + f*pr, f*(j + jp)
    def rhs(self, m, Ab, phib, AT, phiT, eB, T=None):
        out = np.zeros_like(m)
        nd = self._nodes(m, Ab, phib)
        Q = 0.0
        if nd is not None:
            n, p, om, d, mc = nd
            Gam = self.y**2*mc/(8*math.pi)                  # H0 units
            lapse = eB/self.rb                              # d tau_hat / dT
            for c in self.cohorts:
                sl = slice(c['i0'], c['i0'] + self.nk)
                ar = math.exp(Ab - c['Astar']); pc = c['kap']/ar
                omc = np.sqrt(pc*pc + mc*mc)
                out[sl] = -lapse*Gam*(mc/omc)*m[sl]
            for c, fp, sf in self._active(m, self.rb*phiT/eB):
                dudT = c['dir']*phiT*sf/c['dphi']           # s d|phi - phi*|/dT / dphi_r  (>= 0)
                out[c['iu']] = dudT
                out[c['i0']:c['i0'] + self.nk] += c['Nfin']*fp*dudT
            Q = Gam*mc*float(np.sum(n))
        out[0] = -4*AT*m[0] + eB*self.b*Q/self.rb**2
        return out
    def hat_rates(self, m, Ab, phib, hh, vh, T=None, lapse=None, phiT=None):
        """d/dtau_hat of (rho, p, j_adiabatic) (H0 units) along the flow (hh = H/H0, vh = dphi/dtau_hat), including decay and
        production (the derivative of j_prod itself is not included: it enters only the ghost-coupled shell diagnostics)"""
        nd = self._nodes(m, Ab, phib)
        if nd is None: return 0.0, 0.0, 0.0
        n, p, om, d, mc = nd
        G2 = self.G**2; g = G2*d
        Gam = self.y**2*mc/(8*math.pi)
        ndot = -3*hh*n - Gam*(mc/om)*n
        add = np.zeros_like(n)
        for c, fp, sf in self._active(m, vh):
            ci = self.cohorts.index(c)
            ar = math.exp(Ab - c['Astar'])
            add[ci*self.nk:(ci + 1)*self.nk] = c['Nfin']*fp*(c['dir']*vh*sf/c['dphi'])/ar**3
        ndot = ndot + add
        pdot = -hh*p
        omdot = (-hh*p*p + g*vh)/om
        rhodot = float(np.sum(ndot*om + n*omdot))
        prdot = np.sum(ndot*p*p/(3*om) + n*2*p*pdot/(3*om) - n*p*p*omdot/(3*om*om))
        jdot = np.sum(ndot*g/om + n*G2*vh/om - n*g*omdot/(om*om))
        return float(rhodot), float(prdot), self.sign_j*float(jdot)
    def totals_T(self, m, dm, Ab, phib, AT, phiT, eB, T=None):
        hh = self.rb*AT/eB; vh = self.rb*phiT/eB; lapse = eB/self.rb
        rd, pd, jd = self.hat_rates(m, Ab, phib, hh, vh, T, lapse)
        f = self.b/self.rb
        return dm[0] + f*lapse*rd, dm[0]/3 + f*lapse*pd, f*lapse*jd
    def extras(self, m, Ab, phib, AT, phiT, eB, T=None):
        rho, pr, j, ntot, mc, Q, jp = self.hat_totals(m, Ab, phib, T, phiT, eB)
        f = self.b/self.rb
        vh = self.rb*phiT/eB
        return dict(X=f*rho, Pchi=f*pr, J0=f*(j + jp), n_hat=ntot, m_hat=mc, rho_hat=rho, p_hat=pr, j_hat=j, jprod_hat=jp,
                    Q_hat=Q, P_hat=jp*vh, ncoh=len(self.cohorts), u_last=(float(m[self.cohorts[-1]['iu']]) if self.cohorts else 0.0))
    # ---- production
    def add_cohort(self, m, Astar, q, factor, v_star, phib_now):
        kap, N = instant_nodes(q, self.nk, self.mu0, factor)
        i0 = len(m); iu = i0 + self.nk
        dphi = self.ramp_c*abs(v_star)/math.sqrt(q) if self.mode == 'ramp' else 0.0
        c = dict(Astar=Astar, kap=kap, i0=i0, iu=iu, Nfin=N, dphi=dphi, dir=(1.0 if v_star > 0 else -1.0), vs=self.vs_frac*abs(v_star))
        self.cohorts.append(c)
        u0 = min(abs(phib_now - self.ps)/dphi, 1.0) if dphi > 0 else 1.0
        f0, _ = self._f(u0) if dphi > 0 else (1.0, 0.0)
        return np.concatenate([m, N*f0, [u0]])
