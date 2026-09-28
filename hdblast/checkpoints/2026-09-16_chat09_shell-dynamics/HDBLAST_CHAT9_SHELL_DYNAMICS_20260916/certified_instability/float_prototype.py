#!/usr/bin/env python3
"""NON-RIGOROUS floating-point prototype of the polynomial/Riccati formulation used by the certificate.
Checks that  b(mu2) = (mu2+4) w + (G - 4H + 2 phi) R  has its zero at the orchestrator's root mu2 = -7.7178716."""
import math, json, sys
import numpy as np
def U1(p): return (p*p-1)*(-4/3+(10/3)*p-(4/9)*p**3)
def U(p):
    W=1-p+p**3/3; return 0.5*(p*p-1)**2-(2/3)*W*W
def cone(p, y0, N=40):
    # series beta, s, phi in y
    import numpy.polynomial.polynomial as P
    be=np.zeros(N+1); s=np.zeros(N+1); ph=np.zeros(N+1); ph[0]=p
    Uc=[-1/6,4/3,-5/3,-4/9,17/18,0,-2/27]
    U1c=[(k+1)*Uc[k+1] for k in range(6)]
    def comp(c, ser, n):
        out=np.zeros(n+1); pw=np.zeros(n+1); pw[0]=1
        for k,ck in enumerate(c):
            out+=ck*pw; pw=np.convolve(pw,ser[:n+1])[:n+1]
        return out
    for n in range(1,N+1):
        Us=comp(Uc,ph,n-1); U1s=comp(U1c,ph,n-1)
        f=-np.convolve(be[:n],be[:n])[n-1]-np.convolve(s[:n],s[:n])[n-1]/4-Us[n-1]/6
        g=U1s[n-1]-4*np.convolve(be[:n],s[:n])[n-1]
        be[n]=f/(n+2); s[n]=g/(n+4); ph[n]=s[n-1]/n
    ev=lambda c: sum(c[k]*y0**k for k in range(N+1))
    return ev(be), ev(s), ev(ph)
def rhs(Y, mu2):
    phi,s,H,w,R=Y
    G=U1(phi)/s
    return np.array([s, U1(phi)-4*H*s, -w-s*s/3, -2*H*w, R*R+(2*G-6*H)*R+(mu2+4)*w-(2/3)*s*s])
def run(mu2, p, yb, y0=0.05, n=40000):
    be,s,ph=cone(p,y0); H=1/y0+be; w=H*H-s*s/12+U(ph)/6
    rp=-2.5-math.sqrt(2.25-mu2)
    Y=np.array([ph,s,H,w,rp/y0]); h=(yb-y0)/n
    for i in range(n):
        k1=rhs(Y,mu2);k2=rhs(Y+h/2*k1,mu2);k3=rhs(Y+h/2*k2,mu2);k4=rhs(Y+h*k3,mu2)
        Y=Y+h/6*(k1+2*k2+2*k3+k4)
    phi,s,H,w,R=Y; G=U1(phi)/s
    return (mu2+4)*w+(G-4*H+2*phi)*R, Y
sol=json.load(open("../background/REGISTERED_SHELL_FLOAT_SOLUTION.json"))
for m in (-8,-7.7178716,-7.5):
    b,Y=run(m,sol["phi_h"],sol["v_b"]); print(m,b,Y, "h=",sol["h"])
