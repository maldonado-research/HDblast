#!/usr/bin/env python3
"""Pure rational analytic ceilings; no pulse, response or floating evaluation."""
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path


def main():
    # Uniform derivative total variations have independent exact rational
    # certificates in pulse-norms/GLOBAL_RATIONAL_ENVELOPES.json.
    L=F(2,3)
    a2=L*L
    mass2=2*a2
    K=F(256)
    N=(F(97),F(2926),F(161866))
    jets=(F(1),F(4),F(49))
    # |B'|<=8/e<3; |(uB)'|<4; endpoint-zero f'' has sup<=TV/2<49.
    # pi>3, 1-v<=M^2/(2K^2), 1-v^3<=3M^2/(2K^2),
    # 1+2v^3-3v^5<=9(1-v), and 0<=v<=1.
    u=mass2/(2*K*K)
    d3=3*u
    g2=9*u
    basic=[(3*mass2*jets[j]/2+N[j]/4)/(144*K*K) for j in range(3)]
    T0=basic[0]
    T1=basic[1]+L*jets[0]*d3/72
    T2=basic[2]+(2*L*jets[1]*d3+L*L*jets[0]*g2)/72
    ref=60*u*u/(96*9)
    density_local=L*(3*L*jets[0]+jets[1])/(96*9)
    density_tail=density_local*d3
    moment5=d3/3
    moment7=u*u
    moment9=F(4,3)*u**3
    anomaly=((jets[2]+12*L*L*jets[0])*moment5/16
             +(30*L*L*jets[0]+10*L*jets[1])*moment7/32
             +70*L*L*jets[0]*moment9/64)/18
    rho=(3*L*L*T0+L*T1)/2+a2*jets[0]*ref/2+density_tail
    pressure=(T2+3*L*T1+3*L*L*T0)/6+a2*jets[0]*ref/6+(anomaly+density_tail)/3
    current=T0+jets[0]*ref/4
    ceilings={"rho":F(3,100000),"p":F(1,1000),"current":F(3,1000000)}
    computed={"rho":rho,"p":pressure,"current":current}
    for key,value in computed.items():
        if not value<ceilings[key]:
            raise RuntimeError("Global rational stress ceiling failed: "+key)
    directory=Path(__file__).parent
    proof=directory/"pulse-norms"/"GLOBAL_RATIONAL_ENVELOPES.json"
    data={"status":"PASS","scope":"Global analytic upper ceilings from rational inequalities; no source values or response numerics evaluated",
          "K":"256","domain":"inherited observation times; 0<L<=2/3; H=1,r=2",
          "normalized_derivative_budget_caps":[str(v) for v in N],
          "local_absolute_jet_caps":[str(v) for v in jets],
          "strict_ceilings":{k:str(v) for k,v in ceilings.items()},
          "rational_upper_enclosures":{k:str(v) for k,v in computed.items()},
          "derivative_certificate_sha256":hashlib.sha256(proof.read_bytes()).hexdigest(),
          "note":"These are analytic preparation results, not newly selected numerical acceptance tolerances. The registered run uses the tighter actual deferred interval bounds."}
    path=directory/"GLOBAL_STRESS_CEILINGS.json"
    with path.open("x") as handle:
        json.dump(data,handle,indent=2)
        handle.write("\n")
    print(json.dumps(data,indent=2))


if __name__=="__main__":
    main()
