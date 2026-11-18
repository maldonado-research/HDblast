"""Independent metric shooting check from the original Hamiltonian and momentum equations."""
import argparse
import hashlib
import json
import math
from pathlib import Path

from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def main(source):
    raw=source.read_bytes()
    if hashlib.sha256(raw).hexdigest()!="bf162614a84fd207f2f2736b8b77cd1cb173c2db012469cf1d4fe1ccc50c8d0f":
        raise ValueError("Unexpected B1 reference")
    cfg=json.loads(raw)
    phi_b=cfg["phi_b_static"]
    pars=cfg["params"]
    delta,c,d=(pars[key] for key in ("delta","c","d"))

    def u(phi):
        return -2*phi**6/27+17*phi**4/18-4*phi**3/9-5*phi**2/3+4*phi/3-1/6

    def up(phi):
        return -4*phi**5/9+34*phi**3/9-4*phi**2/3-10*phi/3+4/3

    sig=2*(1-phi_b+phi_b**3/3)+delta*(1+c*phi_b+d*phi_b**2/2)
    dsig=2*(phi_b**2-1)+delta*(c+d*phi_b)
    rho=(sig**2/36-dsig**2/48+u(phi_b)/6)**(-.5)

    def quiet(z,state):
        b,q,f,r=state
        return (q,1-q*q-r*r/6-math.exp(2*b)*u(f)/3,r,math.exp(2*b)*up(f)-3*q*r)

    opts=dict(method="RK45",rtol=2e-12,atol=2e-15,max_step=.00025,dense_output=True)
    bg=solve_ivp(quiet,(0,-.2),[math.log(rho),rho*sig/6,phi_b,-rho*dsig/2],**opts)
    if not bg.success:
        raise RuntimeError(bg.message)

    def end(amplitude):
        def rhs(z,state):
            b0,q,f,r0=bg.sol(z)
            x=(z+.12)/.02
            if abs(x)<1:
                bump=math.exp(1-1/(1-x*x))
                df=-100*x*bump/(1-x*x)**2
            else:
                bump,df=0.,0.
            f1=f+amplitude*bump
            p=-amplitude*df
            r=r0+amplitude*df
            # Eliminate B_T from the two original Einstein constraints.
            new_pot=math.exp(2*(b0+state[0]))*u(f1)
            old_pot=math.exp(2*b0)*u(f)
            return [(2*(new_pot-old_pot)+p*p+r*r-r0*r0-2*p*r/q)/(6*(q-1/q))]
        sol=solve_ivp(rhs,(-.10,-.14),[0.],**opts)
        if not sol.success:
            raise RuntimeError(sol.message)
        return float(sol.y[0,-1])

    root=brentq(end,.001,.003,xtol=1e-15)
    result=dict(method="RK45 Hamiltonian-form static background and original constraint equations for lapse factor B",
                endpoint_delta_B_at_001=end(.001),endpoint_delta_B_at_003=end(.003),
                amplitude=root,endpoint_delta_B_at_root=end(root))
    print(json.dumps(result,indent=2,allow_nan=False))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_summary",type=Path)
    main(parser.parse_args().source_summary)
