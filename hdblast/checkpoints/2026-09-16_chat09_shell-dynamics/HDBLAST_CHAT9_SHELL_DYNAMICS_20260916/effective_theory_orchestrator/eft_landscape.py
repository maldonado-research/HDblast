#!/usr/bin/env python3
"""Shape of the shell-modulus potential in the canonical Einstein-frame field Theta (registered linear detuning).
Coordinate: y_b (bulk proper position of the shell along the BPS wall, phi_b = -tanh y_b)."""
import numpy as np, json
W  = lambda p: 1 - p + p**3/3
W1 = lambda p: p*p - 1
Ip = 1.0357712571566782; c = 2/Ip - 4/3
# integrate dI/dy = -1 + (2W/3) I backwards from y=+30 (I -> 3/(2W(-1)) = 0.9)
N = 600001; ys = np.linspace(30, -30, N); dy = ys[1] - ys[0]
I = np.empty(N); I[0] = 0.9
def dI(y, I_): return -1 + (2*W(-np.tanh(y))/3)*I_
for i in range(N-1):
    y = ys[i]; k1 = dI(y, I[i]); k2 = dI(y+dy/2, I[i]+dy/2*k1); k3 = dI(y+dy/2, I[i]+dy/2*k2); k4 = dI(y+dy, I[i]+dy*k3)
    I[i+1] = I[i] + dy/6*(k1+2*k2+2*k3+k4)
phi = -np.tanh(ys); f = 2*I; fy = 2*dI(ys, I)                 # df/dy_b
Zy = -W(phi)*fy*1.0*W1(phi)/W1(phi)                           # Z_y = Z W_phi^2 = -W f' W_phi = -W df/dy_b
Zy = -W(phi)*fy
ZEy = Zy/f + 1.5*(fy/f)**2                                    # Einstein-frame kinetic coefficient for y_b
V = (1 + c*phi)/f**2                                          # V_E / t
i0 = np.argmin(np.abs(ys)); print("check I(0) =", I[i0], " (I_plus = %.10f)" % Ip, " min Z_E,y =", ZEy.min())
# canonical field Theta(y_b) = int sqrt(ZEy) dy from y=0
s = np.sqrt(np.maximum(ZEy, 0)); Theta = np.concatenate([[0], np.cumsum(0.5*(s[1:]+s[:-1])*np.abs(dy))])
Theta -= Theta[i0]; Theta = -Theta          # sign convention: Theta > 0 toward the AdS_- throat (phi_b -> -1), Theta < 0 toward phi_b -> +1
# NOTE (post-review): phi_b = -1 is NOT the end of field space (Z(-1) = 9/23 finite; EFT continues on the phi = -coth y branch).
print("field range: Theta(y=+30; phi=-1 throat) = %.6f ,  Theta(y=-30; phi=+1) = %.6f" % (Theta[0], Theta[-1]))
print("V_E/t: hilltop %.6f ; throat end %.6f ; phi=+1 end %.6f" % (V[i0], V[0], V[-1]))
dVdTh = np.gradient(V, Theta); d2 = np.gradient(dVdTh, Theta)
print("hilltop check: V''/V*3 = mu2 ->", 3*d2[i0]/V[i0])
print("\n  y_b      phi_b      Theta     V_E/t      dlnV/dTheta   eps_V=0.5(V'/V)^2   eta_V=V''/V")
for yy in [12, 8, 6, 4, 3, 2, 1.5, 1, 0.5, 0.2, 0, -0.2, -0.5, -1, -1.5, -2, -3, -4, -6, -8, -12, -20]:
    i = np.argmin(np.abs(ys - yy))
    print(" %+6.2f  %+.6f  %+.5f  %.6f   %+.5f      %.5f        %+.5f" % (ys[i], phi[i], Theta[i], V[i], dVdTh[i]/V[i], 0.5*(dVdTh[i]/V[i])**2, d2[i]/V[i]))
np.savez("eft_landscape.npz", y=ys, phi=phi, I=I, f=f, ZEy=ZEy, Theta=Theta, V=V)
