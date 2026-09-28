#!/usr/bin/env python3
"""Check (floating point) that the HDBLAST radial-HJ transport equations are the single-field specialisation of
de Boer-Verlinde-Verlinde hep-th/9912012 eqs (16),(21),(23),(24), under the dictionary
    phi_d = sqrt(2) phi,  U_d = -2 W,  V_d = -2 U_hdblast,  Phi_d = 2 Phi = I,  M_d = M = W I'/W_phi,
    beta^phi d/dphi = (6/U_d) dU_d/dphi_d d/dphi_d = (3 W_phi/W) d/dphi.
(23):  beta Phi_d' = 2 Phi_d + 6/U_d            <=>  W_phi I' - (2/3) W I = -1
(24) is NOT verified here (see the message printed at the end).
run: python3 dbvv_transport_check.py"""
import numpy as np, warnings
warnings.simplefilter('ignore')
z = np.linspace(-4, 40, 880001); dz = z[1]-z[0]
A = -(1/3)*(z + (2/3)*np.log(np.cosh(z)) - 1/(6*np.cosh(z)**2) + 1/6)
a2 = np.exp(2*A); phi = -np.tanh(z)
W = 1 - phi + phi**3/3; W1 = phi**2 - 1; W2 = 2*phi
tail = a2[-1]/(10/9)
trap = 0.5*(a2[1:]+a2[:-1])*dz
F = np.concatenate((np.cumsum(trap[::-1])[::-1], [0.0])) + tail
I = F/a2                                   # I(phi_b) = int_{z_b}^inf e^{2(A-A_b)}
# d/dphi = (1/phi_z) d/dz, phi_z = W_phi on the BPS flow
ddphi = lambda g: np.gradient(g, dz)/W1
Ip = (-1 + (2/3)*W*I)/W1          # transport equation used analytically for I'; (23) is tested with the numerical derivative below
Ip_num = ddphi(I)
sel = (np.abs(z) < 1.5)
print("(23)  max |W_phi I' - (2/3) W I + 1| on |z|<1.5 : %.2e" % np.max(np.abs(W1*Ip_num - (2/3)*W*I + 1)[sel]))
Ud = -2*W
print("(16)  max |V_d - (U_d^2/3 - (1/2)(dU_d/dphi_d)^2)| with V_d=-2U : %.2e" % np.max(np.abs(-2*(0.5*W1**2-(2/3)*W**2) - (Ud**2/3 - 0.5*(2*W1**2)))[sel]))
print("(24)  NOT tested: the index structure of dBVV (24) could not be read reliably from the HTML rendering; three trial readings\n      did not match M = W I'/W_phi, so no claim is made about (24). HDBLAST uses its own (box phi)-coefficient equation M = 2 W Phi'/W_phi.")
