"""Independent null-characteristic propagation cross-check for the primary matched pulse."""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np


def forces(fields, plus, minus):
    # plus = u_T + u_Z propagates left. minus = u_T - u_Z propagates toward the wall.
    b,phi=fields[1:]
    w=1-phi+phi**3/3
    dw=phi**2-1
    ew=np.exp(2*b)
    v=(dw*dw/2-2*w*w/3)*ew
    vp=(2*phi*dw-4*w*dw/3)*ew
    # A_T^2-A_Z^2 = A_+ A_-; same polarized identity for A and phi.
    return np.array([-3*plus[0]*minus[0]+2*v/3,
                     3*plus[0]*minus[0]-plus[2]*minus[2]/2-v/3,
                     -1.5*(plus[0]*minus[2]+minus[0]*plus[2])-vp])


def wall_slope(fields,pars):
    b,phi=fields[1:,-1]
    delta,c,d=(pars[key] for key in ('delta','c','d'))
    sigma=2*(1-phi+phi**3/3)+delta*(1+c*phi+d*phi*phi/2)
    sigma1=2*(phi*phi-1)+delta*(c+d*phi)
    return np.exp(b)*np.array([sigma/6,sigma/6,-sigma1/2])


def evolve(fields, velocity, spacing, pars, end):
    # Second-order initial slopes may draw on one extra grid point.
    slope=np.gradient(fields,spacing,axis=1,edge_order=2)
    slope[:,-1]=wall_slope(fields,pars)
    plus=velocity+slope
    minus=velocity-slope
    results=[fields[2,-1]]
    for n in range(round(end/spacing)):
        a=forces(fields,plus,minus)
        from_left=minus[:,:-1]+spacing*a[:,:-1]/2
        from_right=plus[:,2:]+spacing*a[:,2:]/2
        old=fields[:,1:]
        trapezoid=old+spacing*(plus[:,1:]+minus[:,1:])/4
        fn=old.copy()
        an=a[:,1:]
        for _ in range(5):
            mn=from_left+spacing*an/2
            pn=np.empty_like(mn)
            pn[:,:-1]=from_right+spacing*an[:,:-1]/2
            pn[:,-1]=mn[:,-1]+2*wall_slope(fn,pars)
            fn=trapezoid+spacing*(pn+mn)/4
            an=forces(fn,pn,mn)
        if not np.isfinite(fn).all() or fn[-1,-1]>=.5:
            raise ValueError('Lost preproduction branch or finite solution')
        fields,plus,minus=fn,pn,mn
        results.append(fn[2,-1])
    return np.array(results)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary',type=Path)
    parser.add_argument('profiles',type=Path)
    parser.add_argument('runs',type=Path)
    plist=parser.parse_args()
    roots=Path(__file__).resolve().parent
    spec=importlib.util.spec_from_file_location('initial',roots/'evolve_patch.py')
    init=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(init)
    raw=plist.summary.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=init.SOURCE_HASH:
        raise ValueError('Background summary hash mismatch')
    summary=json.loads(raw)
    profiles=json.loads(plist.profiles.read_text())
    report={}
    for grid_name in ('g2','g4','g8'):
        meta=json.loads((plist.runs/grid_name/'metadata.json').read_text())
        h=meta['spacing']
        z=np.arange(-round(-init.LEFT/h),1)*h
        chosen={n:profiles[n] for n in ('quiet','delayed','uncompensated','delayed_negative')}
        initial,_=init.prepare(summary,chosen,z)
        fields={}
        for name,(u,velocity,_) in initial.items():
            fields[name]=evolve(u,velocity,h,summary['params'],init.END)
        T=np.arange(len(fields['quiet']))*h
        cross={}
        for name in ('delayed','uncompensated','delayed_negative'):
            response=fields[name]-fields['quiet']
            def series(n):
                with (plist.runs/grid_name/(n+'.csv')).open() as f:
                    return np.array([float(row['phi_wall']) for row in csv.DictReader(f)])
            spatial_response=(series(name)-series('quiet'))[::2]
            if len(response)!=len(spatial_response):
                raise ValueError('Unmatched output cadence')
            j=np.argmax(abs(response))
            cross[name]=dict(char_peak_delta_phi=float(response[j]),T_peak=float(T[j]),
                             max_residual_vs_position_scheme=float(max(abs(response-spatial_response))),
                             char_prearrival_max=float(max(abs(response[T<=.095]))))
        report[grid_name]=dict(spacing=h,cases=cross,
                               causal_comparator_max_through_T_015=float(max(abs((fields['delayed']-fields['uncompensated'])[T<=.15]))))
    print(json.dumps(report,indent=2,allow_nan=False))


if __name__=='__main__':
    main()
