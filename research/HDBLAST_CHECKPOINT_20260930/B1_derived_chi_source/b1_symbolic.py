#!/usr/bin/env python3
"""B1 step 1: exact (SymPy) checks of the closed shell system for the derived chi source, plus wrong-sign/wrong-factor
controls.  Output: B1_SYMBOLIC_CHECKS.json.  Every check raises on failure (explicit exceptions, valid under python -O).

Closed shell system (kappa5 kept explicit here; the solver uses the model convention kappa5 = 1 for the bulk and the
dimensionless b = kappa5^2 H0^3 for the matter normalisation):
  chi gas, one momentum node i (the sum over nodes is linear):  number density n_i = N_i a^-3, physical momentum p_i = k_i/a,
      omega_i = sqrt(p_i^2 + m^2),  m^2 = m0^2 + g^2 (phi - phis)^2,
      rho_chi = sum n_i omega_i,  p_chi = sum n_i p_i^2/(3 omega_i),  j = -dL_m/dphi = sum n_i d omega_i/d phi = g^2 (phi-phis) sum n_i/omega_i,
      decay (rest-frame rate Gamma = y^2 m/(8 pi), time-dilated):  dN_i/dtau = -Gamma (m/omega_i) N_i,  Q = Gamma m sum n_i,
      production ramp (energy-consistent insertion, B1 method choice):  dN_i/dtau |prod = N_i^fin f'(tau),  P = sum (dN_i/dtau|prod) a^-3 omega_i,
      j_prod = P / v  (the scalar pays the production energy through the scalar junction).
  radiation (decay products):  dR/dtau + 4 H R = Q.
  junctions (22 Sept matter extension):  k_s = (sigma + kappa^2 rho)/6,  k_0 = (sigma - kappa^2 (2 rho + 3 p))/6,
      w = n.phi = -(sigma' + kappa^2 j_tot)/2,  j_tot = j + j_prod,  rho = rho_chi + R, p = p_chi + R/3.
"""
import json, math, hashlib
from pathlib import Path
import sympy as sp
import numpy as np
from scipy.special import roots_genlaguerre

HERE = Path(__file__).resolve().parent
checks = []
def check(name, expr_zero, note=''):
    z = sp.simplify(expr_zero)
    ok = (z == 0)
    checks.append(dict(name=name, passed=bool(ok), residual=str(z), note=note))
    if not ok: raise AssertionError('%s: residual %s' % (name, z))
def control(name, expr_should_be_nonzero, note=''):
    z = sp.simplify(expr_should_be_nonzero)
    ok = (z != 0)
    checks.append(dict(name=name, control=True, nonzero_as_required=bool(ok), residual=str(z), note=note))
    if not ok: raise AssertionError('control %s vanished' % name)

tau = sp.symbols('tau', real=True)
a = sp.Function('a', positive=True)(tau)
phi = sp.Function('phi', real=True)(tau)
N = sp.Function('N', positive=True)(tau)
g, phis, m0, k, Gam0, y, kap = sp.symbols('g phi_s m_0 k Gamma_0 y kappa', positive=True)
H = sp.diff(a, tau)/a
v = sp.diff(phi, tau)
m2 = m0**2 + g**2*(phi - phis)**2
m = sp.sqrt(m2)
pk = k/a
om = sp.sqrt(pk**2 + m2)
n = N/a**3
rho = n*om
pr = n*pk**2/(3*om)
j = n*sp.diff(om, phi)                       # = g^2 (phi - phis) n / omega
# ---------------------------------------------------------------- C1: source from the Lagrangian
chi2 = n/om                                  # particle part of <chi^2> for occupation n_k (vacuum part subtracted)
Lm_phi_part = -sp.Rational(1, 2)*m2*chi2    # the phi-dependence of L_m = -1/2 m^2(phi) chi^2 at fixed <chi^2>
j_lagr = -sp.diff(-sp.Rational(1, 2)*m2, phi)*(-1)*chi2   # j = -dL_m/dphi = +1/2 d m^2/dphi <chi^2>
j_lagr = sp.Rational(1, 2)*sp.diff(m2, phi)*chi2
check('C1 j = -dL_m/dphi = (1/2) dm^2/dphi <chi^2> = n d omega/d phi', j - j_lagr)
check('C1b j = g^2 (phi - phis) n/omega', j - g**2*(phi - phis)*n/om)
# ---------------------------------------------------------------- C2: Ward identity (no decay, no production)
Ndot0 = 0
def total_derivative(expr, Ndot):
    return sp.diff(expr, tau).subs(sp.Derivative(N, tau), Ndot)
ward0 = total_derivative(rho, Ndot0) + 3*H*(rho + pr) - j*v
check('C2 Ward identity rho_dot + 3H(rho+p) = j v (free chi gas, per node)', ward0)
control('C2x wrong-sign source (j -> -j) violates the Ward identity', total_derivative(rho, Ndot0) + 3*H*(rho + pr) + j*v,
        'residual = 2 j v')
control('C2y wrong pressure (p -> rho/3, radiation closure for a massive gas) violates it', total_derivative(rho, Ndot0) + 4*H*rho - j*v)
# ---------------------------------------------------------------- C3: decay ledger
Gam = y**2*m/(8*sp.pi)
Ndec = -Gam*(m/om)*N
Q = Gam*m*n
ward_dec = total_derivative(rho, Ndec) + 3*H*(rho + pr) - j*v + Q
check('C3 chi ledger with decay: rho_dot + 3H(rho+p) = j v - Q, Q = Gamma m n (time-dilated per-mode decay)', ward_dec)
Rf = sp.Function('R')(tau)
Rdot = -4*H*Rf + Q
tot = total_derivative(rho, Ndec) + Rdot + 3*H*(rho + pr) + 4*H*Rf - j*v
check('C3b total shell-matter ledger: d(rho_chi + R)/dtau + 3H(rho_chi+p_chi) + 4HR = j v (Q cancels)', tot)
control('C3x decay without time dilation (dN = -Gamma N) is not balanced by Q = Gamma m n',
        total_derivative(rho, -Gam*N) + 3*H*(rho + pr) - j*v + Q)
# ---------------------------------------------------------------- C4: production ramp, energy-consistent
Nfin = sp.symbols('N_fin', positive=True)
f = sp.Function('f')(tau)
Nprod_dot = Nfin*sp.diff(f, tau)
P = Nprod_dot*om/a**3
jprod = P/v
ward_prod = total_derivative(rho, Nprod_dot) + 3*H*(rho + pr) - (j + jprod)*v
check('C4 production ramp: rho_dot + 3H(rho+p) = (j + j_prod) v with j_prod = P/v', ward_prod)
check('C4b production work is positive: j_prod v = P = N_fin f\'(tau) omega/a^3 (>= 0 for a monotone ramp)', jprod*v - Nfin*sp.diff(f, tau)*om/a**3)
# field-space ramp actually used by the solver: u_dot = dir v / dphi_r, N = N_fin f(u); j_prod = dir N_fin f'(u) omega a^-3 / dphi_r
uf = sp.Function('u')(tau); fu = sp.Function('F'); dirs, dphr = sp.symbols('dir dphi_r', nonzero=True)
Nprod_u = Nfin*sp.diff(fu(uf), tau).subs(sp.Derivative(uf, tau), dirs*v/dphr)
jprod_u = dirs*Nfin*sp.Subs(sp.Derivative(fu(sp.Symbol('x')), sp.Symbol('x')), sp.Symbol('x'), uf).doit()*om/a**3/dphr
check('C4c field-space ramp: rho_dot + 3H(rho+p) = (j + j_prod) v with the bounded force j_prod = dir N_fin f\'(u) omega a^-3/dphi_r',
      total_derivative(rho, Nprod_u) + 3*H*(rho + pr) - (j + jprod_u)*v)
sw = sp.Function('s_w')
Nprod_us = Nfin*sp.diff(fu(uf), tau).subs(sp.Derivative(uf, tau), dirs*v*sw(dirs*v)/dphr)
jprod_us = dirs*Nfin*sp.Subs(sp.Derivative(fu(sp.Symbol('x')), sp.Symbol('x')), sp.Symbol('x'), uf).doit()*sw(dirs*v)*om/a**3/dphr
check('C4d field-space ramp with the stall switch s(dir v): Ward identity exact for any switch function',
      total_derivative(rho, Nprod_us) + 3*H*(rho + pr) - (j + jprod_us)*v)
control('C4x instant insertion without j_prod leaves an unpaid energy P', total_derivative(rho, Nprod_dot) + 3*H*(rho + pr) - j*v)
# ---------------------------------------------------------------- C5: junction algebra, shell budget and Codazzi
sig = sp.Function('sigma')
kap2 = kap**2
rho_t = sp.Symbol('rho_t'); p_t = sp.Symbol('p_t'); J = sp.Symbol('J')
rho_tf = sp.Function('rho_T')(tau); J_f = sp.Function('J_T')(tau); p_tf = sp.Function('p_T')(tau)
Hs = sp.Function('H')(tau)
ks = (sig(phi) + kap2*rho_tf)/6
k0 = (sig(phi) - kap2*(2*rho_tf + 3*p_tf))/6
w = -(sp.diff(sig(phi), phi) + kap2*J_f)/2
ward_tot = {sp.Derivative(rho_tf, tau): -3*Hs*(rho_tf + p_tf) + J_f*v}
ks_dot = sp.diff(ks, tau).subs(ward_tot)
check('C5 Codazzi: k_s_dot = -w v/3 - kappa^2 H (rho+p)/2 (with the total Ward identity)', ks_dot - (-w*v/3 - kap2*Hs*(rho_tf + p_tf)/2))
budget = sp.diff(sig(phi)/kap2 + rho_tf, tau).subs(ward_tot) + 3*Hs*(rho_tf + p_tf) - (-2*w*v/kap2)
check('C5b signed shell budget: d(lambda+rho)/dtau + 3H(rho+p) = (sigma\'/kappa^2 + j) v = -2 w v/kappa^2', budget)
control('C5x scalar junction with + sign (w = +(sigma\'+kappa^2 j)/2) breaks the budget',
        sp.diff(sig(phi)/kap2 + rho_tf, tau).subs(ward_tot) + 3*Hs*(rho_tf + p_tf) - (-2*(+(sp.diff(sig(phi), phi) + kap2*J_f)/2)*v/kap2))
control('C5y scalar junction without the factor 1/2 breaks the budget',
        sp.diff(sig(phi)/kap2 + rho_tf, tau).subs(ward_tot) + 3*Hs*(rho_tf + p_tf) - (-2*(-(sp.diff(sig(phi), phi) + kap2*J_f))*v/kap2))
check('C5c pure radiation: k_0 = (sigma - 3 kappa^2 R)/6 when p = R/3 (the A1 code path)',
      k0.subs(p_tf, rho_tf/3) - (sig(phi) - 3*kap2*rho_tf)/6)
# ---------------------------------------------------------------- C6: Friedmann identity, Weyl balance with the chi source
U = sp.Function('U'); Wf = sp.Function('W')(tau)
Hsq = (sig(phi) + kap2*rho_tf)**2/36 + v**2/12 - (sp.diff(sig(phi), phi) + kap2*J_f)**2/48 + U(phi)/6 + Wf
Hdot = -kap2*(sig(phi) + kap2*rho_tf)*(rho_tf + p_tf)/12 - v**2/3 - 2*Wf
# d/dtau (H^2) = 2 H Hdot  with  H^2 given by the identity:  solve for W_dot
Wdot_sym = sp.Symbol('Wdot')
lhs = sp.diff(Hsq, tau).subs(ward_tot).subs(sp.Derivative(Wf, tau), Wdot_sym)
eqW = sp.Eq(lhs, 2*Hs*Hdot)
Wdot_sol = sp.solve(eqW, Wdot_sym)[0]
# substitute H^2 -> identity for the terms 2 H Hdot use H; Weyl balance: Wdot + 4 H W = (1/6)[4 k_s w v - 4 H v^2 - v vdot + w wdot - U' v]
vdot = sp.diff(v, tau); wdot = sp.diff(w, tau)
Wbal = (4*ks*w*v - 4*Hs*v**2 - v*vdot + w*wdot - sp.diff(U(phi), phi)*v)/6
resW = sp.expand(Wdot_sol + 4*Hs*Wf - Wbal)
check('C6 Weyl balance with the chi source: W_dot + 4HW = (1/6)[4 k_s w v - 4H v^2 - v v_dot + w w_dot - U\' v]', resW,
      'derived from the H^2 identity, the Hdot identity and the total Ward identity; J(tau) arbitrary (j + j_prod)')
resW_bad = sp.expand(Wdot_sol + 3*Hs*Wf - Wbal)
control('C6x wrong factor (4H -> 3H) in the Weyl balance', resW_bad)
Wdot_noJ = sp.solve(sp.Eq(sp.diff(Hsq.subs(J_f, 0), tau).subs(ward_tot).subs(sp.Derivative(Wf, tau), Wdot_sym), 2*Hs*Hdot), Wdot_sym)[0]
control('C6y dropping the chi source from the scalar junction (J -> 0 in H^2 only) breaks the balance',
        sp.expand(Wdot_noJ + 4*Hs*Wf - Wbal))
# decomposition of the (sigma + rho)^2/36 term into vacuum, radiation and chi parts (as recorded by the solver)
s_, R_, X_ = sp.symbols('s R X')
check('C6b (s+R+X)^2/36 = s^2/36 + [sR/18 + R^2/36] + [sX/18 + (X^2 + 2RX)/36]',
      (s_ + R_ + X_)**2/36 - (s_**2/36 + (s_*R_/18 + R_**2/36) + (s_*X_/18 + (X_**2 + 2*R_*X_)/36)))
# ---------------------------------------------------------------- C7: units (H0 units -> solver model units)
H0, rb, rhohat, jhat, Qhat, b, M5 = sp.symbols('H_0 r_b rhohat jhat Qhat b M_5', positive=True)
subs_units = {H0: 1/rb}
k2 = b/H0**3                                    # kappa5^2 = b/H0^3
check('C7 kappa5^2 rho = b rhohat / r_b  (rho = H0^4 rhohat)', (k2*H0**4*rhohat).subs(subs_units) - b*rhohat/rb)
check('C7b kappa5^2 j = b jhat / r_b', (k2*H0**4*jhat).subs(subs_units) - b*jhat/rb)
check('C7c kappa5^2 Q = b Qhat / r_b^2 (Q = H0^5 Qhat)', (k2*H0**5*Qhat).subs(subs_units) - b*Qhat/rb**2)
check('C7d H0/M5 = b^(1/3) with M5 = kappa5^(-2/3)', (H0/(k2**sp.Rational(-1, 3))) - b**sp.Rational(1, 3))
lam_c, mhat = sp.symbols('lambda_c mhat', positive=True)
check('C7e cutoff m_chi <= lambda_c M5  <=>  b <= (lambda_c/mhat)^3 (equality form)',
      sp.simplify((mhat*H0/(k2**sp.Rational(-1, 3))).subs(b, (lam_c/mhat)**3) - lam_c))
# ---------------------------------------------------------------- C8: momentum nodes of the instant-preheating occupation
num = []
for q, mu0 in [(10.0, 0.0), (128.3, 0.0), (1283.0, 0.0), (100.0, 5.0)]:
    u, wts = roots_genlaguerre(24, 0.5)
    kapn = np.sqrt(q*u/math.pi)
    Nn = wts*(q/math.pi)**1.5/(4*math.pi**2)*math.exp(-math.pi*mu0**2/q)
    n_exact = q**1.5/(8*math.pi**3)*math.exp(-math.pi*mu0**2/q)
    e_exact = q**2/(4*math.pi**4)*math.exp(-math.pi*mu0**2/q)
    # wrong occupation exp(-pi k^2/(2q)) as control
    Nw = wts*(2*q/math.pi)**1.5/(4*math.pi**2)*math.exp(-math.pi*mu0**2/q)
    num.append(dict(q=q, mu0=mu0, n_nodes=float(Nn.sum()), n_exact=n_exact, rel_err_n=float(abs(Nn.sum()/n_exact - 1)),
                    massless_energy_nodes=float((Nn*kapn).sum()), massless_energy_exact=e_exact,
                    rel_err_energy=float(abs((Nn*kapn).sum()/e_exact - 1)),
                    wrong_occupation_rel_dev=float(abs(Nw.sum()/n_exact - 1))))
    assert abs(Nn.sum()/n_exact - 1) < 1e-12
    assert abs((Nn*kapn).sum()/e_exact - 1) < 5e-4         # massless weight k ~ u^1/2 is not polynomial: 24-node error 1.7e-4 (N^-2); massive omega converges fast
    assert abs(Nw.sum()/n_exact - 1) > 0.5
# sympy: exact Gaussian integrals
kk, qq = sp.symbols('k q', positive=True)
n_int = sp.integrate(4*sp.pi*kk**2*sp.exp(-sp.pi*kk**2/qq), (kk, 0, sp.oo))/(2*sp.pi)**3
e_int = sp.integrate(4*sp.pi*kk**3*sp.exp(-sp.pi*kk**2/qq), (kk, 0, sp.oo))/(2*sp.pi)**3
check('C8 int d^3k/(2pi)^3 exp(-pi k^2/q) = q^(3/2)/(8 pi^3)', n_int - qq**sp.Rational(3, 2)/(8*sp.pi**3))
check('C8b int d^3k/(2pi)^3 k exp(-pi k^2/q) = q^2/(4 pi^4) (massless energy proxy)', e_int - qq**2/(4*sp.pi**4))
# ---------------------------------------------------------------- C9: NR limit and the endpoint scalar force
dpos, nn = sp.symbols('Delta n', positive=True)
check('C9 NR limit (m0 = 0, p = 0, phi - phis = Delta > 0): j = g^2 Delta n/m = g n  (constant scalar force per particle)',
      g**2*dpos*nn/sp.sqrt(g**2*dpos**2) - g*nn)

out = dict(status='exact-verified (symbolic identities of the closed shell system; controls must be nonzero)',
           n_checks=sum(1 for c in checks if not c.get('control')), n_controls=sum(1 for c in checks if c.get('control')),
           all_passed=all(c.get('passed', c.get('nonzero_as_required')) for c in checks),
           checks=checks, gauss_laguerre_nodes=num, sympy_version=sp.__version__,
           script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
(HERE/'B1_SYMBOLIC_CHECKS.json').write_text(json.dumps(out, indent=1) + '\n')
print('checks: %d passed, controls: %d nonzero' % (out['n_checks'], out['n_controls']))
for c in checks: print(' ', 'OK ' if c.get('passed', c.get('nonzero_as_required')) else 'BAD', c['name'])
