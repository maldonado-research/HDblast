#!/usr/bin/env python3
"""HDBLAST Chat 9 - analytic O(t) coefficient of the shell tachyon mass mu^2(t) = mu0 + s1 t + ...
Curvature expansion (X = 1/rho^2) of the 5D boundary-value problem with the bulk scalar as radial coordinate.
See SLOPE_DERIVATION.md.  numpy + stdlib only.
Orientation: y increases from the AdS_- throat (phi=-1) to the shell; BPS wall phi = tanh(y), A' = W/3, phi' = -W_phi.
"""
import math, json, sys
import numpy as np
W   = lambda p: 1 - p + p**3/3
W1  = lambda p: p*p - 1
W2  = lambda p: 2*p
W3  = lambda p: 2.0
LM  = 9.0/5.0                  # AdS_- radius: 1/l = W(-1)/3 = 5/9
Z2C = -LM**3/8                 # z2 = -l^3/8 = -729/1000  (no bulk black-hole mass / regular cone)

def _rhs(y, S):
    p = math.tanh(y); w, w1, w2 = W(p), -1.0/math.cosh(y)**2, W2(p)     # W_phi = -sech^2 (stable in the throat)
    Q = S[0]
    I = 1.5*(1 - w*Q)/w; Ip = -w*Q/w1
    return np.array([w1*w1/(w*w) - (2*w/3)*Q, Ip*Ip, (w - w2)*I*Ip*Ip])

def wall(phis, L=20.0, dy=1e-3):
    """Returns dict phi -> (I, K1, K2):  I = e^{-2A} int_{-inf}^y e^{2A},  K1 = int I_phi^2 dy,  K2 = int (W - W_phiphi) I I_phi^2 dy."""
    out = {}; y = -L; S = np.zeros(3)
    for p in sorted(phis):
        yt = math.atanh(p); n = max(2, int(math.ceil((yt - y)/dy))); h = (yt - y)/n
        for _ in range(n):
            k1 = _rhs(y, S); k2 = _rhs(y + h/2, S + h/2*k1); k3 = _rhs(y + h/2, S + h/2*k2); k4 = _rhs(y + h, S + h*k3)
            S = S + h/6*(k1 + 2*k2 + 2*k3 + k4); y += h
        w = W(p); out[p] = (1.5*(1 - w*S[0])/w, S[1], S[2])
    return out

def local(p, I, K1, K2, Vf):
    """All curvature-expansion coefficient functions at scalar value p.  Vf(p) -> (V, V', V'') detuning profile."""
    w, w1, w2, w3 = W(p), W1(p), W2(p), W3(p)
    V, Vp, Vpp = Vf(p)
    I1 = (2*w*I/3 - 1)/w1
    I2 = (-w2*I1 + (2/3)*(w1*I + w*I1))/w1
    z2 = Z2C
    p1 = -6*I1; p1d = -6*I2
    p2 = (6/w1)*(-2*w*z2/3 - I*I + p1*p1/12)
    p2d = -(w2/w1)*p2 + (6/w1)*(-2*w1*z2/3 - 2*I*I1 + p1*p1d/6)
    k0 = 2*w/(3*w1); k1 = k0*(3*I/w + p1/w1)
    G0 = -w1*(k0*p1 + p1d)
    G1 = p1*(k0*p1 + p1d) - w1*(k1*p1 + 2*k0*p2 + p2d)
    T0 = G0 - 3*w1*Vpp*I/V
    T1 = G1 + 3*Vpp*(p1*I - w1*z2)/V
    s0 = -I1/(2*I)
    s1a = ((I**3 + z2)/3 - 6*K2)/(I*I*w1)
    s1b = K1/(4*I*I*w1)
    mu0 = -4 - T0/(3*s0)
    return dict(I=I, I1=I1, I2=I2, p1=p1, p2=p2, z2=z2, G0=G0, G1=G1, T0=T0, T1=T1, s0=s0, s1a=s1a, s1b=s1b,
                mu0=mu0, ratio=p1/I + 3*Vp/V, x1=V/(6*I), K1=K1, K2=K2)

def slope(phi0, Vf_maker, e=2e-3):
    """phi0 = leading-order equilibrium; Vf_maker(phi0, I, I1) returns Vf tuned so that phi0 is the equilibrium."""
    pts = [phi0 + k*e for k in (-2, -1, 0, 1, 2)]
    tab = wall(pts)
    I0 = tab[phi0][0]; I10 = (2*W(phi0)*I0/3 - 1)/W1(phi0)
    Vf = Vf_maker(phi0, I0, I10)
    loc = [local(p, *tab[p], Vf) for p in pts]
    d = lambda key: (loc[0][key] - 8*loc[1][key] + 8*loc[3][key] - loc[4][key])/(12*e)
    L0 = loc[2]
    beta = -L0['x1']*(L0['p2']/L0['I'] - L0['p1']*L0['z2']/L0['I']**2)/d('ratio')
    s1 = L0['s1a'] + (3*L0['mu0'] + 12)*L0['s1b']
    sl = beta*d('mu0') - L0['x1']*(L0['T1']/(3*L0['s0']) - L0['T0']*s1/(3*L0['s0']**2))
    x2 = None
    L0.update(beta=beta, s1_total=s1, slope=sl, equilibrium_residual=L0['ratio'], Vf=Vf)
    return L0

def linear_detuning(phi0, I0, I10, d=0.0):
    g = 2*I10/I0                                  # V'/V = 2 I'/I at equilibrium
    c = (g*(1 + d*phi0*phi0/2) - d*phi0)/(1 - g*phi0)
    Vf = lambda p: (1 + c*p + d*p*p/2, c + d*p, d)
    Vf.c = c; Vf.d = d
    return Vf

if __name__ == "__main__":
    r = slope(0.0, linear_detuning)
    Ip = r['I']; c = r['Vf'].c
    print("I_plus = %.13f  c_star = %.13f  (2/I-4/3 = %.13f)" % (Ip, c, 2/Ip - 4/3))
    print("K1 = %.12f   K2 = %.12f   z2 = %.6f" % (r['K1'], r['K2'], r['z2']))
    print("mu0 = %.10f  closed form %.10f" % (r['mu0'], -4*(3*c*c - 4*c + 8)/(c*(3*c + 4))))
    print("beta = phi_b/t = %.10f   (5D float: 2.279e-5/1e-3 = 0.02279)" % r['beta'])
    print("s0 = %.10f  s1a = %.10f  s1b = %.10f" % (r['s0'], r['s1a'], r['s1b']))
    print("SLOPE s1 = d mu^2/dt = %.10f   (5D fit: 1.9243)" % r['slope'])
