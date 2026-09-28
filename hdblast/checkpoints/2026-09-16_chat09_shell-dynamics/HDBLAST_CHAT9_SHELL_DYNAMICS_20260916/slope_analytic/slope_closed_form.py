#!/usr/bin/env python3
"""Closed form of s1 = d mu^2/dt|_{t=0} for the registered shell (phi_b -> 0, linear detuning, c = c_star) in terms of
c, zeta = l_-^3/8 = 729/1000 and the two wall quadratures K1, K2;  cross-check against the general routine; table over c."""
import json, math
import slope_formula as sf

def closed_form(c, K1, K2, zeta=729/1000):
    a = 6/(3*c + 4)                                  # a = I_plus
    Q = 3*c*c - 4*c + 8
    mu0 = -4*Q/(c*(3*c + 4))
    beta = (a**3*(4 - 3*c*c) - 4*zeta)/(2*a**3*Q)    # phi_b = beta t + O(t^2)
    dmu0 = 8*a*c/3 - 56*a/3 + (32*a/9)*(2 - c)/c**2  # d mu0(phi)/d phi at phi = 0 along the equilibrium family (V = 1 + c phi fixed)
    G0 = 4*a*(c - 1)
    G1 = a*a*(6*c*c + 12*c - 8) + 28*zeta/3
    S = -((a**3 - zeta)/3 - 6*K2)/a**2 - (3*mu0 + 12)*K1/(4*a*a)      # s1(0) = first curvature correction of psi/chi
    s0 = -c/4
    parts = dict(position_shift=beta*dmu0, junction_curvature=-(1/(6*a))*G1/(3*s0), mode_profile=(1/(6*a))*G0*S/(3*s0*s0))
    return sum(parts.values()), dict(mu0=mu0, beta=beta, dmu0_dphi=dmu0, S=S, **parts)

if __name__ == "__main__":
    res = {}
    for dy in (1e-3, 5e-4):
        sf_wall = sf.wall
        sf.wall = lambda phis, L=20.0, dy=dy: sf_wall(phis, L, dy)
        r = sf.slope(0.0, sf.linear_detuning); sf.wall = sf_wall
        c = r['Vf'].c; s, parts = closed_form(c, r['K1'], r['K2'])
        print("dy=%g: I+=%.14f c*=%.14f K1=%.14f K2=%.14f" % (dy, r['I'], c, r['K1'], r['K2']))
        print("   general routine slope = %.12f    closed form = %.12f" % (r['slope'], s)); print("  ", parts)
        res[str(dy)] = dict(I_plus=r['I'], c_star=c, K1=r['K1'], K2=r['K2'], slope_general=r['slope'], slope_closed=s, **parts)
    print("\nTable: linear detuning, equilibrium at phi0 (c fixed by phi0):   phi0   c   mu0   beta   s1")
    tab = []
    for p0 in (-0.9, -0.8, -0.6, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8):
        r = sf.slope(p0, sf.linear_detuning)
        print("  %+.2f  %.10f  %+.10f  %+.10f  %+.10f" % (p0, r['Vf'].c, r['mu0'], r['beta'], r['slope']))
        tab.append(dict(phi0=p0, c=r['Vf'].c, mu0=r['mu0'], beta=r['beta'], slope=r['slope'], K1=r['K1'], K2=r['K2']))
    res['table'] = tab
    json.dump(res, open("SLOPE_CLOSED_FORM.json", "w"), indent=1)
