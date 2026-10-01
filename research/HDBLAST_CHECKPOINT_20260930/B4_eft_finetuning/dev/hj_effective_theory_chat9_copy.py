#!/usr/bin/env python3
"""HDBLAST Chat 9 - Hamilton-Jacobi (holographic-RG) 4D effective theory of the Z2 shell in the BPS bulk.
  S_4 = int sqrt(-g) [ f(phi) R/2 - (1/2) Z(phi) (d phi)^2 - V(phi) ],   phi = bulk scalar value at the shell,
  f = 2 I,   W_phi I' - (2/3) W I = -1  (IR-regular solution: I = int_{y_b}^inf exp(2(A-A_b)) dy along phi=-tanh y),
  Z = - W f'/W_phi = 2 W (1 - 2 W I/3)/W_phi^2,     V = sigma_t - 2W  (detuning).
Einstein frame: V_E = V/f^2,  Z_E = Z/f + (3/2)(f'/f)^2;   at a de Sitter stationary point  m^2/H^2 = 3 (ln V_E)''/Z_E.
"""
import numpy as np, json, math
W  = lambda p: 1 - p + p**3/3
W1 = lambda p: p*p - 1
W2 = lambda p: 2*p
def A_of_y(y): return -(1/3)*(y + (2/3)*np.log(np.cosh(y)) - 1/(6*np.cosh(y)**2))
_x, _w = np.polynomial.legendre.leggauss(80)
def I_of_phi(pb, L=90.0, panels=900):
    yb = -np.arctanh(pb); ed = np.linspace(yb, yb + L, panels + 1); tot = 0.0
    for i in range(panels):
        m = 0.5*(ed[i] + ed[i+1]); hh = 0.5*(ed[i+1] - ed[i])
        tot += hh*np.sum(_w*np.exp(2*(A_of_y(m + hh*_x) - A_of_y(yb))))
    return tot
def eft(pb):
    I = I_of_phi(pb); Wv, Wp, Wpp = W(pb), W1(pb), W2(pb)
    I1 = (-1 + 2*Wv*I/3)/Wp
    I2 = ((2*Wp*I/3 + 2*Wv*I1/3)*Wp - (-1 + 2*Wv*I/3)*Wpp)/Wp**2
    f, f1, f2 = 2*I, 2*I1, 2*I2
    Z = -Wv*f1/Wp
    ZE = Z/f + 1.5*(f1/f)**2
    return dict(phi_b=pb, I=I, f=f, f1=f1, f2=f2, Z=Z, Z_E=ZE)
def mu2_linear_detuning(pb, d=0.0):
    """sigma_t = 2W + t(1 + c phi + d phi^2/2); stationarity fixes c given phi_b; returns c, mu2 (t->0)."""
    e = eft(pb); g = 2*e['f1']/e['f']                     # need V'/V = g
    # V = 1 + c p + d p^2/2 ; V' = c + d p  => (c + d p) = g (1 + c p + d p^2/2)
    c = (g*(1 + d*pb*pb/2) - d*pb)/(1 - g*pb)
    V = 1 + c*pb + d*pb*pb/2; Vp = c + d*pb
    lnVE2 = d/V - (Vp/V)**2 - 2*(e['f2']/e['f'] - (e['f1']/e['f'])**2)
    e.update(c=c, d=d, lnVE_pp=lnVE2, mu2=3*lnVE2/e['Z_E'], VE_over_t=V/e['f']**2, h_over_t=V/(3*e['f']))
    return e
if __name__ == "__main__":
    r0 = mu2_linear_detuning(0.0)
    print("registered shell (phi_b -> 0, d = 0):"); print(json.dumps(r0, indent=1))
    c = r0['c']; print("closed form  -4(3c^2-4c+8)/(c(3c+4)) =", -4*(3*c*c - 4*c + 8)/(c*(3*c + 4)))
    print("\nscan of equilibrium position (linear detuning, d = 0):")
    print("  phi_b        I          Z          Z_E        c          (lnV_E)''    mu2=m^2/H^2   (I^2)''>0?")
    rows = []
    for pb in [-0.95, -0.9, -0.8, -0.6, -0.4, -0.2, -0.1, 0.0, 0.1, 0.2, 0.4, 0.6, 0.8, 0.9, 0.95]:
        e = mu2_linear_detuning(pb); rows.append(e)
        conc = 2*(e['f1']/2)**2 + 2*e['I']*e['f2']/2
        print("  %+.2f   %.6f   %+.6f  %+.6f  %+.6f  %+.6f   %+.6f    %+.4f" % (pb, e['I'], e['Z'], e['Z_E'], e['c'], e['lnVE_pp'], e['mu2'], conc))
    # quadratic detuning needed for slow-roll hilltop at phi_b = 0
    e0 = eft(0.0)
    for target in (-0.0525, -0.06, -0.03):
        # mu2 = 3 (d - c^2 - 2(ln f)'')/Z_E  with c fixed by stationarity (independent of d at phi_b = 0)
        lnf2 = e0['f2']/e0['f'] - (e0['f1']/e0['f'])**2
        dneed = target*e0['Z_E']/3 + c*c + 2*lnf2
        print("target mu2 = %+.4f (eta_V = %+.5f, n_s ~ %.4f): d = %.8f" % (target, target/3, 1 + 2*target/3, dneed))
    json.dump(dict(registered=r0, scan=rows), open("HJ_EFT_TABLE.json", "w"), indent=1)
