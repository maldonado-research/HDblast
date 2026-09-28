#!/usr/bin/env python3
"""VERIFICATION (math agent) - brute-force pointwise check of the linearised 5D Einstein-scalar equations.

Independent of any hand-derived perturbation formula: we build the full perturbed metric
    ds^2 = (1 + 2 eps xi Y) dy^2 + rho(y)^2 (1 + 2 eps psi Y) gamma,   gamma = -dtau^2 + e^{2 tau} dx^2  (unit dS_4)
    phi  = phi0(y) + eps chi(y) Y,     Y = F(tau) cos(k x1 + th),  Box_gamma Y = mu2 Y
and compute E_AB = R_AB - d_A phi d_B phi - (2/3) U g_AB and E_phi = Box phi - U_phi numerically:
  * derivative with respect to eps: complex step (eps = 1e-25 i) -> exact linearisation,
  * coordinate derivatives: 4th-order central differences of the LOCAL 2-jets (Taylor polynomials) of the 1D functions.
    (R_AB at a point depends only on the 2-jet of the metric, so truncating each 1D function to its Taylor
     polynomial at the point is exact; the jets are generated from the claimed ODEs.)
Claimed ODE system under test (orchestrator R1, written in first-order form):
    xi = -2 psi,
    psi' = -2 (rho'/rho) psi - phi' chi/3                                   (Codazzi, C)
    chi' = 3 (mu2 + 4) psi/(rho^2 phi') - 2 phi' psi + (phi''/phi') chi      (Gauss H combined with C)
 which is algebraically equivalent to the master equation (M_y) [shown in VERIFICATION_REPORT].
Random point data (rho, phi, phi', psi, chi, F, F', tau, k, mu2) and random POTENTIAL parameters are used.
Sabotage runs show that the test is sensitive to each coefficient.
"""
import numpy as np, math, sys
rng = np.random.default_rng(20260916)

def make_case(sab=None):
    sab = sab or {}
    # random cubic superpotential (the identities must hold for any U); registered W is used in case 0
    return sab

def run(seed, registered=True, sab=None, verbose=False):
    sab = sab or {}
    r = np.random.default_rng(seed)
    if registered:
        W  = lambda p: 1 - p + p**3/3
        W1 = lambda p: p*p - 1
        W2 = lambda p: 2*p
        U  = lambda p: 0.5*W1(p)**2 - (2/3)*W(p)**2
        U1 = lambda p: W1(p)*W2(p) - (4/3)*W(p)*W1(p)
    else:   # generic potential, not from a superpotential
        a0, a1, a2, a3 = r.uniform(-1, 1, 4); a0 = -abs(a0) - 0.5
        U  = lambda p: a0 + a1*p + a2*p*p + a3*p**3
        U1 = lambda p: a1 + 2*a2*p + 3*a3*p*p
    mu2 = r.uniform(-12, 3); k = r.uniform(0.3, 2.0); th = r.uniform(0, 6.28)
    tau0 = r.uniform(-0.5, 0.5); y0 = 0.0; xx0 = r.uniform(-1, 1, 3)
    rho0 = r.uniform(0.7, 2.5); ph0 = r.uniform(-0.9, 0.9); s0 = r.uniform(0.3, 1.2)*r.choice([-1, 1])
    psi0 = r.uniform(-1, 1); chi0 = r.uniform(-1, 1); F0 = r.uniform(0.5, 1.5); F1 = r.uniform(-1, 1)
    cXI = sab.get('xi', -2.0); cMU = sab.get('mu_shift', 4.0); cC = sab.get('codazzi', 1/3); cG = sab.get('gauss2', 2.0)

    def flow(st):
        rho, ph, s, psi, chi = st
        rp = np.sqrt(1 + rho*rho*(s*s/12 - U(ph)/6)); Hh = rp/rho
        spp = U1(ph) - 4*Hh*s
        dpsi = -2*Hh*psi - cC*s*chi
        dchi = 3*(mu2 + cMU)*psi/(rho*rho*s) - cG*s*psi + (spp/s)*chi
        return np.array([rp, s, spp, dpsi, dchi])
    st0 = np.array([rho0, ph0, s0, psi0, chi0], complex)
    d1 = flow(st0)
    hh = 1e-30
    d2 = np.imag(flow(st0 + 1j*hh*d1.real))/hh          # second derivatives via complex-step directional derivative
    d1 = d1.real
    # F jets
    F2 = -3*F1 - k*k*math.exp(-2*tau0)*F0 - mu2*F0
    poly = lambda c0, c1, c2: (lambda d: c0 + c1*d + 0.5*c2*d*d)
    rho_f = poly(rho0, d1[0], d2[0]); ph_f = poly(ph0, d1[1], d2[1])
    psi_f = poly(psi0, d1[3], d2[3]); chi_f = poly(chi0, d1[4], d2[4]); F_f = poly(F0, F1, F2)
    eps = 1e-25j

    def fields(x):
        dy = x[0] - y0; dt = x[1] - tau0
        Y = F_f(dt)*math.cos(k*x[2] + th)
        rho = rho_f(dy); psi = psi_f(dy)
        g = np.zeros((5, 5), complex)
        g[0, 0] = 1 + 2*eps*cXI*psi*Y
        conf = rho*rho*(1 + 2*eps*psi*Y)
        g[1, 1] = -conf
        e2 = math.exp(2*x[1])
        for i in (2, 3, 4): g[i, i] = conf*e2
        return g, ph_f(dy) + eps*chi_f(dy)*Y

    x0 = np.array([y0, tau0, *xx0])
    h = 2e-3
    st = [(-2, 1/12), (-1, -8/12), (1, 8/12), (2, -1/12)]
    st2 = [(-2, -1/12), (-1, 16/12), (0, -30/12), (1, 16/12), (2, -1/12)]
    cache = {}
    def ev(off):
        key = tuple(off)
        if key not in cache: cache[key] = fields(x0 + h*np.array(off, float))
        return cache[key]
    g0, p0 = ev([0]*5)
    dg = np.zeros((5, 5, 5), complex); dp = np.zeros(5, complex)
    ddg = np.zeros((5, 5, 5, 5), complex); ddp = np.zeros((5, 5), complex)
    for a in range(5):
        for j, w in st:
            off = [0]*5; off[a] = j; G, P = ev(off); dg[a] += w*G/h; dp[a] += w*P/h
        for j, w in st2:
            off = [0]*5; off[a] = j; G, P = ev(off); ddg[a, a] += w*G/h**2; ddp[a, a] += w*P/h**2
        for b in range(a + 1, 5):
            for j, w in st:
                for j2, w2 in st:
                    off = [0]*5; off[a] = j; off[b] = j2; G, P = ev(off)
                    ddg[a, b] += w*w2*G/h**2; ddp[a, b] += w*w2*P/h**2
            ddg[b, a] = ddg[a, b]; ddp[b, a] = ddp[a, b]
    gi = np.linalg.inv(g0)
    dgi = np.array([-gi @ dg[e] @ gi for e in range(5)])
    # S[e,b,c,d] = d_e ( d_b g_dc + d_c g_db - d_d g_bc ),  T[b,c,d] = d_b g_dc + d_c g_db - d_d g_bc
    T = np.einsum('bdc->bcd', dg) + np.einsum('cdb->bcd', dg) - np.einsum('dbc->bcd', dg)
    S = np.einsum('ebdc->ebcd', ddg) + np.einsum('ecdb->ebcd', ddg) - np.einsum('edbc->ebcd', ddg)
    Gam = 0.5*np.einsum('ad,bcd->abc', gi, T)
    dGam = 0.5*np.einsum('ead,bcd->eabc', dgi, T) + 0.5*np.einsum('ad,ebcd->eabc', gi, S)   # dGam[e,a,b,c] = d_e Gamma^a_bc
    Ric = np.einsum('aabc->bc', dGam) - np.einsum('caba->bc', dGam) \
        + np.einsum('aad,dbc->bc', Gam, Gam) - np.einsum('acd,dba->bc', Gam, Gam)
    Uc = U(p0); U1c = U1(p0)
    E = Ric - np.outer(dp, dp) - (2/3)*Uc*g0
    boxphi = np.einsum('ab,ab->', gi, ddp) - np.einsum('ab,cab,c->', gi, Gam, dp)
    Ephi = boxphi - U1c
    # normalise components with the orthonormal-frame scale
    sc = np.sqrt(np.abs(np.diag(g0.real)))
    Ebg = (E.real/np.outer(sc, sc)); Elin = (E.imag/abs(eps))/np.outer(sc, sc)
    out = dict(bg=np.abs(Ebg).max(), lin=np.abs(Elin).max(), bg_phi=abs(Ephi.real), lin_phi=abs(Ephi.imag/abs(eps)),
               scale=np.abs((Ric.imag/abs(eps))/np.outer(sc, sc)).max(), mu2=mu2)
    if verbose:
        np.set_printoptions(precision=2, linewidth=150); print(Elin)
    return out

if __name__ == "__main__":
    print("== correct equations, registered W: max |background residual|, max |linear residual| (frame comps), scalar eq, scale of dR")
    worst = 0
    for seed in range(12):
        o = run(seed, True); worst = max(worst, o['lin'], o['lin_phi'], o['bg'], o['bg_phi'])
        print("seed %2d mu2=%+7.3f  bg %.1e  lin %.1e  bgphi %.1e  linphi %.1e   (|dRic| ~ %.1e)" % (seed, o['mu2'], o['bg'], o['lin'], o['bg_phi'], o['lin_phi'], o['scale']))
    print("== correct equations, GENERIC potential U (no superpotential)")
    for seed in range(100, 108):
        o = run(seed, False); worst = max(worst, o['lin'], o['lin_phi'], o['bg'], o['bg_phi'])
        print("seed %2d mu2=%+7.3f  bg %.1e  lin %.1e  bgphi %.1e  linphi %.1e   (|dRic| ~ %.1e)" % (seed, o['mu2'], o['bg'], o['lin'], o['bg_phi'], o['lin_phi'], o['scale']))
    print("WORST residual with the claimed equations: %.2e" % worst)
    print("== sabotage runs (seed 3): each must give O(1e-2..1) residuals")
    for name, sab in [("xi = -2.05 psi", dict(xi=-2.05)), ("(mu2+4) -> (mu2+4.1)  [i.e. (2+mu2)->(2.1+mu2) in M_y]", dict(mu_shift=4.1)),
                      ("Codazzi 1/3 -> 0.35", dict(codazzi=0.35)), ("-2 phi' psi -> -2.1 phi' psi", dict(gauss2=2.1))]:
        o = run(3, True, sab)
        print("  %-60s lin %.2e  linphi %.2e" % (name, o['lin'], o['lin_phi']))
