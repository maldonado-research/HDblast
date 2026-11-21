"""Fixed-energy Born comparison of the crossed channels of one local diagonal Dirac vector."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad

M=8300.
WIDTH=.30
G2=12*np.pi*WIDTH
ID=np.eye(2);ZERO=np.zeros((2,2))
PAULI=np.array([[[0,1],[1,0]],[[0,-1j],[1j,0]],[[1,0],[0,-1]]])
GAMMA=np.array([np.block([[ID,ZERO],[ZERO,-ID]])]+[np.block([[ZERO,sig],[-sig,ZERO]]) for sig in PAULI])
ETA=np.array([1.,-1.,-1.,-1.])


def blocks(x,y,z,m=1.):
    """Two spin sums and the interference invariant; x,y are exchange variables."""
    m2=m*m
    first=2*((y-2*m2)**2+(z-2*m2)**2+4*m2*x)
    second=2*((x-2*m2)**2+(z-2*m2)**2+4*m2*y)
    cross=4*(z*z-8*m2*z+12*m2*m2)
    return first,second,cross


def propagators(s,t,u,m=1.):
    ds=1./(s-M*M+1j*WIDTH*(s+2*m*m)*np.sqrt(1.-4.*m*m/s))
    return ds,1./(t-M*M),1./(u-M*M)


def terms(s,z,m=1.):
    t=-.5*(s-4*m*m)*(1.-z)
    u=-.5*(s-4*m*m)*(1.+z)
    ds,dt,du=propagators(s,t,u,m)
    a,b,c=blocks(s,t,u,m); mt,mu,mi=blocks(t,u,s,m)
    return np.array([a*abs(ds)**2,b*dt**2,c*ds.real*dt,mt*dt**2,mu*du**2,mi*dt*du])


def spin(p,m,anti):
    n=np.sqrt(p[0]+m);large=n*ID;small=np.einsum("i,ijk->jk",p[1:],PAULI)/n
    return np.concatenate((small,large) if anti else (large,small))


def amp(p,q):
    return (p.conj().T@GAMMA[0]) @ GAMMA @ q


def spin_contract(s,z,m):
    E=np.sqrt(s)/2.;k=np.sqrt(E*E-m*m);v=k*np.sqrt((1.-z)*(1.+z))
    p1,p2,p3,p4=np.array([[E,0.,0.,k],[E,0.,0.,-k],[E,v,0.,k*z],[E,-v,0.,-k*z]])
    U1,U2,U3,U4=[spin(p,m,False) for p in (p1,p2,p3,p4)]
    V2,V4=spin(p2,m,True),spin(p4,m,True)
    S=np.einsum("mba,mcd,m->abcd",amp(V2,U1),amp(U3,V4),ETA)
    T=np.einsum("mca,mbd,m->abcd",amp(U3,U1),amp(V2,V4),ETA)
    TT=np.einsum("mca,mdb,m->abcd",amp(U3,U1),amp(U4,U2),ETA)
    UU=np.einsum("mda,mcb,m->abcd",amp(U4,U1),amp(U3,U2),ETA)
    x=-.5*(s-4*m*m)*(1.-z);y=-.5*(s-4*m*m)*(1.+z)
    ds,dt,du=propagators(s,x,y,m)
    a,b,c,d=ds*S,dt*T,dt*TT,du*UU
    return np.array([np.sum(abs(a)**2),np.sum(abs(b)**2),-2*np.vdot(b,a).real,
                     np.sum(abs(c)**2),np.sum(abs(d)**2),-2*np.vdot(d,c).real])/4.


def summary(s):
    def integrate(i,weight):
        def elem(z):
            a,b,c,d,e,f=terms(s,z)
            value=(a,a+b+c,(d+e+f)/2.)[i]
            return (G2**2/(32*np.pi*s))*weight(z)*value
        return quad(elem,-1.,1.,epsabs=0.,epsrel=1e-10,limit=160)[0]
    out={}
    for key,w in (("total",lambda z:1.),("transfer",lambda z:1.-z),
                  ("quadrupole",lambda z:1.5*(1.-z*z)),("identical_transfer",lambda z:1.-abs(z))):
        pole,unlike,like=[integrate(i,w) for i in range(3)]
        out[key]=dict(pole=pole,unlike=unlike,like=like,unlike_over_pole=unlike/pole,like_over_pole=like/pole)
    return out


def main(output):
    checks=[]
    for m,s,z in ((1.,4.01,-.15),(1.,18.,.4),(1.,4.*M*M,-.79),(1.,M*M,.12),(0.,1.4*M*M,-.6)):
        contracted=spin_contract(s,z,m);analytic=terms(s,z,m)
        scale=np.maximum(abs(analytic),1e-8*np.max(abs(analytic)))
        delta=float(np.max(abs(contracted-analytic)/scale))
        if delta>1.e-9: raise ValueError(("Spin contractions disagree",m,s,z,delta))
        checks.append(dict(m=m,s=s,z=z,spin_contraction_error=delta))
    grid=[4.001,0.25*M*M,M*M,8.63903342e7,4.*M*M,4.81731786e8,3.61011614e9,100.*M*M]
    rows=[dict(s=s,root_s_over_M=np.sqrt(s)/M,rates=summary(s)) for s in grid]
    for r in rows:
        s=r["s"];v=WIDTH*(s+2.);beta=np.sqrt(1.-4./s)
        frozen=12.*np.pi*v*v/(s*((s-M*M)**2+(beta*v)**2))
        diff=abs(r["rates"]["total"]["pole"]/frozen-1.)
        r["parent_total_relative_error"]=diff
        if diff>1e-12:raise ValueError("Unmatched isolated pole normalization")
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(dict(controls=checks,rows=rows),indent=2)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output",type=Path,required=True)
    main(parser.parse_args().output)
