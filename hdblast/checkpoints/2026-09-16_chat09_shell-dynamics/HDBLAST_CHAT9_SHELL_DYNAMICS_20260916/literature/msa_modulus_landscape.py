#!/usr/bin/env python3
"""MSA (moduli-space approximation, Brax et al. hep-th/0209158 method) landscape of the single-brane modulus for the
registered HDBLAST wall.  Leading order in the detuning t.  HEURISTIC (4D effective theory), not a 5D result.
   z = BPS proper coordinate, phi = -tanh z, a = e^A, A(0) = 0, bulk kept on z > z_b (the phi -> -1 side).
   f(z) = int_z^inf a^2,  V_E/t = (1 + c phi) a^4/(4 f^2),  G_zz = (3/2) a^4/f^2 - a^2 W/f,  chi = int sqrt(G_zz) dz.
   eta_V = V_E''/(V_E) in canonical field units (M_pl = 1) = [ (lnV)_zz + (lnV)_z^2 - (lnV)_z G_z/(2G) ] / G  ; eps_V = (lnV)_z^2/(2G).
run: python3 msa_modulus_landscape.py"""
import numpy as np, json, os
A = lambda z: -(1/3)*(z + (2/3)*np.log(np.cosh(z)) - 1/(6*np.cosh(z)**2) + 1/6)
Wf = lambda p: 1 - p + p**3/3
z = np.linspace(-14, 40, 540001); dz = z[1] - z[0]
a2 = np.exp(2*A(z))
# f(z) = int_z^inf a^2 : cumulative trapezoid from the right + analytic tail (a^2 ~ e^{-(10/9) z})
tail = a2[-1]/(10/9)
trap = 0.5*(a2[1:] + a2[:-1])*dz
f = np.concatenate((np.cumsum(trap[::-1])[::-1], [0.0])) + tail      # summed from the right: no cancellation on the plateau
I = np.interp(0.0, z, f); c = 2/I - 4/3
phi = -np.tanh(z)
VE = (1 + c*phi)*a2**2/(4*f**2)
G_naive = 1.5*a2**2/f**2 - a2*Wf(phi)/f
# cancellation-free identity (integration by parts, uses a^2 = int_z^inf (2W/3) a^2 and W_z = W_phi^2 >= 0):
#     G_zz = (a^2/f^2) int_z^inf W_phi(z')^2 f(z') dz'   >= 0   (vanishes only for a constant scalar = Randall-Sundrum)
w1f = (phi**2 - 1)**2*f
trapG = 0.5*(w1f[1:] + w1f[:-1])*dz
G = a2/f**2*np.concatenate((np.cumsum(trapG[::-1])[::-1], [0.0]))
print("max |G_identity - G_naive| for |z|<6 : %.2e" % np.max(np.abs(G - G_naive)[np.abs(z) < 6]))
lnV = np.log(VE); d1 = np.gradient(lnV, dz); d2 = np.gradient(d1, dz)
i0 = np.argmin(abs(z))
print("I_plus = %.10f, c_star = %.10f" % (I, c))
print("at z=0: (lnV)_z = %.3e, (lnV)_zz = %.6f, G_zz = %.6f, mu^2 = %.5f" % (d1[i0], d2[i0], G[i0], 3*d2[i0]/G[i0]))
dd = d1[2000:-2000]; sgn = np.sign(dd)
ch = [k for k in np.where(sgn[1:]*sgn[:-1] < 0)[0] if max(abs(dd[k]), abs(dd[k+1])) > 1e-7]   # ignore rounding noise on the flat plateaus
print("stationary points of V_E (|slope| above rounding noise) at z =", [round(float(z[2000 + k]), 4) for k in ch])
chi = np.concatenate(([0.0], np.cumsum(0.5*(np.sqrt(np.abs(G[1:])) + np.sqrt(np.abs(G[:-1])))*dz))); chi -= chi[i0]
print("min G_zz on grid (z<39) = %.3e (positivity of the Einstein-frame kinetic term)" % G[z < 39].min())
print("canonical distance hilltop -> phi_b=-1 end: %.5f ;  hilltop -> phi_b=+1 end: %.5f (grid to z=-14) + %.5f (tail 9 sqrt(G))" % (
      chi[-1], -chi[0], 9*np.sqrt(G[0])))
rows = []
print("   z      phi_b     V_E/t      G_zz       chi(canonical, M_pl=1)   H^2/H^2(hilltop)")
for zz in [-12, -8, -4, -2, -1, -0.5, 0, 0.5, 1, 2, 4, 8, 12, 20]:
    k = np.argmin(abs(z - zz))
    print("%6.1f  %8.5f  %9.6f  %9.3e  %9.5f   %8.5f" % (z[k], phi[k], VE[k], G[k], chi[k], VE[k]/VE[i0]))
    rows.append(dict(z=float(z[k]), phi_b=float(phi[k]), VE_over_t=float(VE[k]), G_zz=float(G[k]), chi=float(chi[k])))
k5 = (5/9); k1 = 1/9
print("analytic plateaus: V_E/t -> (1-c) k^2 = %.6f (phi_b -> -1, k = 5/9);   (1+c) k'^2 = %.6f (phi_b -> +1, k' = 1/9)" % ((1 - c)*k5**2, (1 + c)*k1**2))
json.dump(dict(I_plus=float(I), c_star=float(c), rows=rows), open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "msa_modulus_landscape_output.json"), "w"), indent=1)
