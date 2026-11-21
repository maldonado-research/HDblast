"""Test monoenergetic-shell drag and Synge primitive recovery by independent quadrature."""
import argparse
import json
import math
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np
from scipy.special import kve,roots_legendre
HERE=Path(__file__).resolve().parent


def direct(t,u,order):
    # Two independent directions in the two rest frames; integrate the relative
    # azimuth explicitly. This calculation has no primitive table or rapidity transform.
    e=1.5*t+np.hypot(1.,1.5*t); k=np.sqrt(3*t*e); h=e+t; gamma=np.hypot(1.,u)
    x,w=roots_legendre(order); z1=x[:,None,None]; z2=x[None,:,None]
    phi=(np.arange(2*order)+.5)*np.pi/order
    et1=gamma*e+u*k*z1;et2=gamma*e-u*k*z2
    pz1=gamma*k*z1+u*e;pz2=gamma*k*z2-u*e
    side=k*k*np.sqrt((1-z1*z1)*(1-z2*z2))*np.cos(phi)
    s=np.maximum(4.,2*(1+et1*et2-pz1*pz2-side)); beta=np.sqrt(1-4/s)
    g=.30*(s+2); flux=6*np.pi*beta*g*g/((s-8300.**2)**2+beta*beta*g*g)
    val=((pz1-pz2)*flux).mean(axis=2)
    return float(w @ val @ w/(8*gamma*u*h*e**2))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    with tempfile.TemporaryDirectory(prefix='hdblast-closure-tests-') as tmp:
        d=Path(tmp)
        runpy.run_path(str(HERE/'flow.py'))['eos_table'](d/'synge_table.inc')
        cpp=d/'check.cpp'
        cpp.write_text('''#include "synge.hpp"
#include "shell.hpp"
#include <iomanip>
#include <iostream>
int main(){double t,u,extra; while(std::cin>>t>>u>>extra) {
 const auto e=synge::thermal(t);
 std::cout<<std::setprecision(17)<<shellbath::rate(t,u,24)<<" "<<e.specific<<" "<<
 synge::invert((1.+e.specific+t)*u,extra,1.01*t)<<"\\n";}}''')
        subprocess.run(['g++','-std=c++17','-O3','-I',str(HERE),'-I',str(d),str(cpp),'-o',str(d/'check')],check=True)
        pts=[(.01,.2),(16.90075,30.112),(4100.,.4),(4765.,.001),(720.,1.1),(1e-8,200.),(.5,0.),(1.,4.),(1e5,10.),(1e-5,.0002)]
        data=[]; eos=[]
        for t,u in pts:
            gamma=np.hypot(1.,u)
            ej=1.5*t+1.875*t*t-1.875*t**3 if t<1e-5 else 3*t+kve(1,1/t)/kve(2,1/t)-1
            h=1+ej+t; mom=h*u
            # Rationalized lab excess avoids subtracting large equal ultra-cold energies.
            delta=(ej+t)*(2.+ej+t)/(h*gamma+np.hypot(1,mom))-t/gamma
            data.append(f'{t:.17g} {u:.17g} {delta:.17g}\n');eos.append(ej)
        output=subprocess.check_output([str(d/'check')],input=''.join(data),text=True)
    result=[]
    for i, line in enumerate(output.splitlines()):
        rate,eth,rec=map(float,line.split());t,u=pts[i]
        row=dict(theta=t,proper_speed=u,kernel=rate,specific_energy_relative_difference=eth/eos[i]-1.,
                 temperature_recovery_relative_difference=rec/t-1)
        if 1e-8<t<1e5 and u>0:
            r48=direct(t,u,48);r72=direct(t,u,72)
            row['angular_reference']=r72;row['angular_order_change']=r72/r48-1
            row['rate_relative_difference']=rate/r72-1.
            if abs(row['rate_relative_difference'])>.0005:raise ValueError(row)
        if max(abs(row[k]) for k in ('temperature_recovery_relative_difference','specific_energy_relative_difference'))>1e-7:
            raise ValueError(row)
        result.append(row)
    a.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')

if __name__=='__main__':main()
