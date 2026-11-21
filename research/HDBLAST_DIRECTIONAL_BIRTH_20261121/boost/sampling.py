"""Independent two-bath collision integral using Lorentz-boosted Maxwell–Jüttner samples."""
import argparse
import json
from pathlib import Path
import math
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.stats import qmc

HERE=Path(__file__).resolve().parent
POLE_MASS=8300.
WIDTH=.30
PROBES={'cool':(.03,.2),'resonant_cold_beams':(.0001,9.),'moderate':(10.,4.),
        'hot':(3000.,.7),'hot_slow':(3000.,.05),'boosted_hot':(100.,5.)}

def momenta(theta,uniform):
    # In kinetic-energy variable x=(E/m-1)/Theta, p^2 dp is proportional to
    # (1+Theta*x)*sqrt(x*(2+Theta*x)) dx. Use sqrt(x) for a regular endpoint.
    t=np.linspace(0,math.sqrt(55.),12001)
    x=t*t
    pdf=2*x*(1+theta*x)*np.sqrt(2+theta*x)*np.exp(-x)
    cdf=cumulative_trapezoid(pdf,t,initial=0)
    xi=np.interp(uniform,cdf/cdf[-1],t)**2
    energy=1+theta*xi
    momentum=np.sqrt(theta*xi*(2+theta*xi))
    return energy,momentum

def rate(theta,rapidity,sequence):
    E1,p1=momenta(theta,sequence[:,0])
    E2,p2=momenta(theta,sequence[:,1])
    u1,u2=2*sequence[:,2]-1,2*sequence[:,3]-1
    transverse=p1*p2*np.sqrt((1-u1*u1)*(1-u2*u2))*np.cos(2*np.pi*sequence[:,4])
    v,gamma=math.tanh(rapidity),math.cosh(rapidity)
    enthalpy=2.5*theta+math.hypot(1.,1.5*theta)
    # Parity partners reduce cancellation of the signed momentum change.
    combined=0.
    for sign in (1.,-1.):
        z1,z2=sign*u1*p1,sign*u2*p2
        lab1,lab2=gamma*(E1+v*z1),gamma*(E2-v*z2)
        pz1,pz2=gamma*(z1+v*E1),gamma*(z2-v*E2)
        invariant=np.maximum(4.,2+2*(lab1*lab2-pz1*pz2-transverse))
        beta=np.sqrt(np.maximum(0.,(invariant-4.)/invariant))
        running=WIDTH*(invariant+2)*beta
        # m² times Møller invariant relative flux times transfer cross section.
        # Writing BW/beta as beta*[WIDTH*(s+2)]^2/den avoids 0/0 at threshold.
        flux_xs=6*np.pi*beta*np.square(WIDTH*(invariant+2))/(np.square(invariant-POLE_MASS**2)+np.square(running))
        combined+=np.mean((pz1-pz2)*flux_xs/(E1*E2))/2
    return float(combined/(2*gamma**2*enthalpy*v))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    sequences=[qmc.Sobol(d=5,scramble=True,seed=s).random_base2(16) for s in (113,194,271,349)]
    result={}
    for name,(theta,r) in PROBES.items():
        rates=[rate(theta,r,seq) for seq in sequences]
        result[name]=dict(theta=theta,rapidity=r,mean=float(np.mean(rates)),draws_per_seed=65536,
                    scrambles=rates,relative_standard_error=float(np.std(rates,ddof=1)/(2*np.mean(rates))))
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+'\n')

if __name__=='__main__':
    main()
