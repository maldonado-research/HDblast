"""Independent fixed-step Gaussian loss channel and a static damped-mode benchmark."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

spec = importlib.util.spec_from_file_location("inputs", Path(__file__).with_name("reservoir.py"))
inputs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inputs)


def step(u,c,v,a,mu,k,y,dt):
    e = np.hypot(k/a, mu)
    b = a**3*e
    gamma = y*y*mu*mu/(8*math.pi*e)
    loss = -np.expm1(-gamma*dt)
    co,sn = np.cos(e*dt), np.sin(e*dt)
    # Symplectic rotation in the frozen Hamiltonian, followed by exact vacuum loss.
    u1 = co*co*u+2*co*sn*c/b+sn*sn*v/b**2
    c1 = -co*sn*u*b+(co*co-sn*sn)*c+co*sn*v/b
    v1 = sn*sn*u*b*b-2*co*sn*b*c+co*co*v
    n = (b*u+v/b-1)/2
    return ((1-loss)*u1 + loss/(2*b), (1-loss)*c1,
            (1-loss)*v1 + loss*b/2, a*e*n*loss)


def static_check():
    a,mu,y,dt = 1.2,-3.7,1.1,.003
    k=np.array([1.,3.,5.])
    e=np.hypot(k/a,mu)
    b=a**3*e
    n0=np.array([.4,.1,.3])
    u=(n0+.5)/b
    c=np.zeros_like(k)
    v=(n0+.5)*b
    z=np.zeros_like(k)
    for _ in range(700):
        u,c,v,dz=step(u,c,v,a,mu,k,y,dt)
        z+=dz
    n=(b*u+v/b-1)/2
    exact=n0*np.exp(-y*y*mu*mu/(8*math.pi*e)*700*dt)
    defect=max(np.max(abs(n-exact)),np.max(abs(z-a*e*(n0-exact))))
    if defect>5e-11:raise ValueError("Static damping benchmark")
    return float(defect)


def advance(geom, y, steps):
    k,w=inputs.rule(48)
    af,_,mf,_=geom(2.1)
    a0,_,m0,_=geom(.8)
    ef=np.hypot(k/af,mf)
    b0=a0**3*np.hypot(k/a0,m0)
    u,c,v=.5/b0,np.zeros_like(k),.5*b0
    dt=(2.1-.8)/steps
    a,_,mu,_=geom(.8+(np.arange(steps)+.5)*dt)
    z=0.
    for ai,mi in zip(a,mu,strict=True):
        u,c,v,dz=step(u,c,v,ai,mi,k,y,dt)
        z+=w@dz
    n=(af**3*ef*u+v/(af**3*ef)-1)/2
    if min(u*v-c*c)<.25-1e-8 or min(n)<-1e-8:raise ValueError("Nonphysical final field moments")
    return dict(band_comoving_number=w@n,chi_comoving_energy=w@(ef*n),daughter_physical_density=z/af**4)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("previous_capsule",type=Path)
    b=p.parse_args().previous_capsule
    out={"static_damping_defect":static_check()}
    for steps in (1800,3600):
        out[str(steps)]={arm: advance(inputs.background(b,"fine",arm),1.,steps) for arm in inputs.BRANCHES}
    print(json.dumps(out, indent=2,allow_nan=False))
