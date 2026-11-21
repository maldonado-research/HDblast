"""Independent checks of the 3-D coordinate transforms and magnetic angle integral."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.integrate import quad
from scipy.special import roots_legendre

HERE=Path(__file__).resolve().parent
PROBE=r"""
#include <cmath>
#include <iostream>
#include <iomanip>
struct Vector {double x,y,z;};
struct Particle {double qperp,qz;int hits;double qx,qy;};
#include "oriented.hpp"
int main(){
    std::cout<<std::setprecision(16);
    int type;double a,b,c,d;
    while(std::cin>>type>>a>>b>>c>>d){
        if(type==0){std::cout<<average_inverse_energy(a,b,c,d)<<'\n';continue;}
        if(type==1){
            const Vector n=final_basis({a,b,c},int(d)), m=charge_partner(n);
            std::cout<<n.x<<' '<<n.y<<' '<<n.z<<' '<<m.x<<' '<<m.y<<' '<<m.z<<'\n';continue;
        }
        if(type==2){
            Particle p=canonical_particle({a,b,c},1.,0.,0);
            rotate_gyro(p,d);
            std::cout<<p.qx<<' '<<p.qy<<' '<<p.qz<<'\n';continue;
        }
        Vector n=quarter_turn({a,b,c},int(d));
        std::cout<<n.x<<' '<<n.y<<' '<<n.z<<'\n';
    }
}
"""


def run():
    average_probes=[(0.,0.,0.,1.), (1.,3.,3.,2.), (.001,-450.,.003,2.),
                    (.3,-1e4,2e5,3.),(8e3,2e5,-2e5,2.),(.13,65.,65.00001,5.),
                    (0.,23000.,23100.,.9)]
    rows=[(0,*v) for v in average_probes]
    phi=.916
    vector=np.array([.4,-.31,.78])
    reflect=np.array([1.,-1.,-1.])
    rows.extend([(2,*vector,phi),(2,*(reflect*vector),-phi)])
    rows.extend([(3,*vector,j) for j in range(4)])
    rows.extend([(1,*vector,j) for j in range(3)])
    with tempfile.TemporaryDirectory(prefix="hdblast-oriented-check-") as tmp:
        dest=Path(tmp); src=dest/"check.cpp";exe=dest/"check"
        src.write_text(PROBE)
        subprocess.run(["g++","-O3","-std=c++17","-I",str(HERE),str(src),"-o",str(exe)],check=True)
        lines="".join(" ".join(map(str,row))+"\n" for row in rows)
        run=subprocess.run([str(exe)],input=lines,capture_output=True,check=True,text=True)
    out=[np.fromstring(x,sep=" ") for x in run.stdout.splitlines()]
    ratios=[]
    for index,(q,z0,z1,r) in enumerate(average_probes):
        reference=quad(lambda t:1/math.sqrt(1+(q/r)**2+((z0+(z1-z0)*t)/r)**2),
                       0,1,points=[min(1.,max(0.,-z0/(z1-z0)))] if z1!=z0 else [],epsabs=1e-14)[0]
        ratios.append(float(out[index][0]/reference-1))
    i=len(average_probes)
    expected=np.array([math.cos(phi)*vector[0]+math.sin(phi)*vector[1],
                       -math.sin(phi)*vector[0]+math.cos(phi)*vector[1],vector[2]])
    signed_error=float(np.max(abs(out[i]-expected)))
    mirror_error=float(np.max(abs(reflect*out[i]-out[i+1])))
    rots=out[i+2:i+6]
    mean_transverse=np.mean(rots,axis=0)[:2]
    harmonics=[(v[0]+1j*v[1])**4 for v in rots]
    fourth_error=float(max(abs(k-harmonics[0]) for k in harmonics))
    mapped=out[i+6:i+9]
    basis_error=float(max(np.max(abs(v[:3]-vector[idx])) for v,idx in zip(mapped,([1,2,0],[2,0,1],[0,1,2]))))
    if max(max(abs(x) for x in ratios),signed_error,mirror_error,fourth_error,basis_error)>1e-10:
        raise ValueError("3-D transform or magnetic-advance check failed")
    return dict(gyro_integral_probes=average_probes,relative_errors_against_direct_quadrature=ratios,
                positive_rotation_absolute_error=signed_error,charge_conjugate_rotation_absolute_error=mirror_error,
                quarter_turn_current_residual=float(np.max(abs(mean_transverse))),
                quarter_turn_fourth_harmonic_error=fourth_error,cyclic_destination_basis_error=basis_error)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(run(),allow_nan=False,sort_keys=True,indent=2)+"\n")
