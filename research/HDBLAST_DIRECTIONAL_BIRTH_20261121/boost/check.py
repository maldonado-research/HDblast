"""Check the finite-mass collision integral and the same C++ table reader used in the live evolution."""
import argparse
import json
import math
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np
from scipy.stats import qmc

HERE=Path(__file__).resolve().parent
T=runpy.run_path(str(HERE/'thermal.py'))
rate=T['rate']


def compile_reader(points,table):
    cpp='''#include "lookup.hpp"
#include <iostream>
#include <iomanip>
int main(int argc,char** argv) {
    load_thermal_table(argv[1]);
    double theta, r;
    while(std::cin >> theta >> r)
        std::cout << std::setprecision(17) << boosted_rate(theta,std::sinh(r)) << "\\n";
}
'''
    with tempfile.TemporaryDirectory() as scratch:
        source=Path(scratch)/'probe.cpp'; exe=source.with_suffix('')
        source.write_text(cpp)
        subprocess.run(['g++','-std=c++17','-O2','-I',str(HERE),str(source),'-o',str(exe)],check=True)
        text=''.join(f'{t:.17g} {r:.17g}\n' for t,r in points)
        values=subprocess.check_output([str(exe),str(table)],input=text,text=True)
    return list(map(float,values.splitlines()))


def check(table, sampling_path):
    points=list(map(tuple,np.column_stack((10**(-5+11*(a:=qmc.Sobol(2,scramble=True,seed=901).random_base2(8))[:,0]),12*a[:,1]))))
    pole=math.acosh(T['POLE']/2)
    points += [(theta,pole+off) for theta in (1e-5,3.9e-5,.0009,.01,.055,.45)
               for off in np.linspace(-.7,.7,57)]
    points += [(theta,r) for theta in (1e-10,1e-7,9.99e-6) for r in (0.,1e-6,.005,.05,.5,3.,8.8,9.,pole,9.3,11.)]
    points += [(1.,0.),(10.,1e-10),(100.,1e-8),(1e5,.001),(3000.,.01)]
    values=compile_reader(points,table)
    comparisons=[]
    for (theta,r),interpolated in zip(points,values):
        exact=rate(theta,r,order=32)
        default=rate(theta,r,order=16)
        comparisons.append(dict(theta=theta,rapidity=r,interpolation_rel=abs(interpolated/exact-1),
                                quadrature_rel=abs(default/exact-1)))
    reference=json.loads(sampling_path.read_text())
    sampling={}
    for name,row in reference.items():
        exact=rate(row['theta'],row['rapidity'],order=32)
        sampling[name]=dict(exact=exact,sobol_mean=row['mean'],sobol_relative_standard_error=row['relative_standard_error'],
                            relative_difference=exact/row['mean']-1)
    sigma0=float(T['m2sigma'](0.))
    controls={}
    for t in (1e-8,1e-6):
        controls[f'nonrelativistic_zero_drift_{t:.0e}']=rate(t,0)/(16*sigma0*math.sqrt(t)/(3*math.sqrt(math.pi)))
    for r in (.3,3.,pole,11.):
        controls[f'cold_beam_at_r_{r:.6g}']=rate(1e-9,r)/rate(0.,r)
    controls['ultrarelativistic_constant_cross_section_over_5_4']=rate(5000.,0.,cross_section=np.ones_like)/1.25
    controls['ultrarelativistic_s_wave_unitarity_over_pi_over_theta_squared']=rate(5000.,0.,cross_section=lambda e:4*np.pi/np.sinh(e)**2)*5000.**2/np.pi
    controls['low_drift_continuity_at_theta_3000']=rate(3000.,1e-9)/rate(3000.,0.)
    # Independent massless pole average using the published low-drift kernel in z=sqrt(s)/T.
    from scipy.integrate import quad
    from scipy.special import kv
    for t in (500.,2000.,10000.):
        ratio=T['POLE']/t
        def kernel(z):
            return (z**3*kv(2,z)+z**4*kv(1,z))/(24*(1+((ratio**2/z**2-1)/T['FRACTION'])**2))
        fraction=sum(quad(kernel,lo,hi,epsabs=1e-13)[0] for lo,hi in ((0,ratio),(ratio,max(ratio+65,100))))
        controls[f'massless_vector_zero_drift_{t:g}']=rate(t,0.)/(.75*np.pi*fraction/t**2)
    stats=dict(count=len(points),
        max_quadrature=max(comparisons,key=lambda q:q['quadrature_rel']),
        max_interpolation=max((q for q in comparisons if q['theta']>=1e-5),key=lambda q:q['interpolation_rel']),
        max_cold_continuation=max((q for q in comparisons if q['theta']<1e-5),key=lambda q:q['interpolation_rel']),
        p95_interpolation=float(np.quantile([x['interpolation_rel'] for x in comparisons],.95)))
    if stats['max_quadrature']['quadrature_rel']>1e-5 or stats['max_interpolation']['interpolation_rel']>.01 or stats['max_cold_continuation']['interpolation_rel']>.002:
        raise ArithmeticError(stats)
    if any(abs(x-1)>.0001 for x in controls.values()) or max(abs(x['relative_difference']) for x in sampling.values())>.015:
        raise ArithmeticError((controls,sampling))
    return dict(table=stats,limit_control_ratio=controls,sobol=sampling)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--table',type=Path,required=True)
    p.add_argument('--sampling',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(check(args.table,args.sampling),sort_keys=True,indent=2,allow_nan=False)+'\n')
