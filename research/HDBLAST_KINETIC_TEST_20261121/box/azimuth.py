"""Angular fourth moment of the first products from a common three-flavor pole."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np
from scipy.special import eval_legendre, roots_legendre

HERE = Path(__file__).resolve().parent
kernel = runpy.run_path(str(HERE / "run.py"))["kernel"]


def azimuthal_moment(energy, split, n):
    energies = [1+.15*(energy-1), 1+1.85*(energy-1)] if split else [energy]
    c, w = roots_legendre(n)
    rates, moments = [], []
    for ea in energies:
        for eb in energies:
            pa, pb = math.sqrt(ea**2-1), math.sqrt(eb**2-1)
            s = 2+2*(ea*eb+pa*pb)
            beta2 = 1-4/s
            a, b = 3-2*beta2, beta2**2
            angular = (a+b*c**2)/(2*(a+b/3))
            velocity = (pa-pb)/(ea+eb)
            gamma = (ea+eb)/math.sqrt(s)
            pcm = math.sqrt(s/4-1)
            lab_z = gamma*(pcm*c+velocity*math.sqrt(s)/2)
            lab_perp = pcm*np.sqrt(1-c*c)
            lab_cosine = lab_z/np.hypot(lab_z,lab_perp)
            rates.append(kernel(ea, eb))
            moments.append(w @ (angular*eval_legendre(4, lab_cosine)))
    mean_p4 = float(np.average(moments, weights=rates))
    # For H(n)=Re[(n_x+i n_y)^4], zonal averaging about axes x,y,z
    # gives H(axis) <P4(n.axis)> with H(x)=H(y)=1 and H(z)=0.
    return dict(input_axis_mean_p4=mean_p4, destination_z_transverse_fourth=2*mean_p4/3)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    e=json.loads(args.input.read_text())["initial_energies"]
    result={}
    for label,energy in e.items():
        for split in (False,True):
            main=azimuthal_moment(energy,split,96)
            check=azimuthal_moment(energy,split,48)
            main["max_quadrature_change"]=max(abs(main[k]-check[k]) for k in main)
            result[label+("_split" if split else "_mono")]=main
    args.output.write_text(json.dumps(result,sort_keys=True,indent=2,allow_nan=False)+"\n")
