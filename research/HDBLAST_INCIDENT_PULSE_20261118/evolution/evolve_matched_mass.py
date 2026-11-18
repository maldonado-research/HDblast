"""Propagate the near-root mass-matched scalar and its canonical-radius quiet control."""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


def module(path):
    spec=importlib.util.spec_from_file_location(path.stem,path)
    result=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('summary',type=Path)
    parser.add_argument('mass_module',type=Path)
    parser.add_argument('output',type=Path)
    parser.add_argument('--spacing',type=float,required=True)
    args=parser.parse_args()
    eq=module(Path(__file__).resolve().with_name('evolve_patch.py'))
    mass=module(args.mass_module.resolve())
    h=args.spacing
    z=np.arange(-round(-eq.LEFT/h),1)*h
    patch=mass.MatchedPulse(args.summary)
    amplitude=patch.root()
    mass_solution=patch.mass_solution(amplitude)

    def background(position, fields):
        b,q,phi,r=fields
        U,Up=eq.potential(phi)
        return q,3-3*q*q-2*np.exp(2*b)*U/3,r,np.exp(2*b)*Up-3*q*r

    inside=z>=mass.PATCH_LEFT
    quiet=np.empty((4,len(z)))
    quiet[:,inside]=patch.static(z[inside])
    extended=solve_ivp(background,(mass.PATCH_LEFT,eq.LEFT),patch.static(mass.PATCH_LEFT),
                       method='DOP853',max_step=.0005,rtol=2e-13,atol=2e-15,dense_output=True)
    if not extended.success:
        raise ArithmeticError(extended.message)
    quiet[:,~inside]=extended.sol(z[~inside])
    B,Q,phi,r=quiet
    fields=np.array([B,B,phi])
    speeds=np.zeros_like(fields)
    speeds[0]=1
    corrected_fields=fields.copy()
    corrected_speeds=speeds.copy()
    for i in np.flatnonzero((z>mass.PULSE_LEFT)&(z<mass.PULSE_RIGHT)):
        d=patch.fields(z[i],amplitude,mass_solution)
        corrected_fields[1,i]=d['B']
        corrected_fields[2,i]=d['phi']
        corrected_speeds[1,i]=d['BT']
        corrected_speeds[2,i]=d['p']
    summary=json.loads(args.summary.read_text())
    args.output.mkdir(parents=True,exist_ok=True)
    meta=dict(archived_radius=summary['rho_b'],canonical_radius=patch.radius,amplitude=amplitude,
              residual_mass_at_left=float(mass_solution(mass.PULSE_LEFT)[0]),
              summary_sha256=eq.SOURCE_HASH,spacing=h,cases={})
    for name,f,v in [('quiet',fields,speeds),('mass_matched',corrected_fields,corrected_speeds)]:
        rows=eq.evolve(f,v,z,summary['params'],h)
        path=args.output/(name+'.csv')
        with path.open('w',newline='') as handle:
            w=csv.writer(handle,lineterminator='\n')
            w.writerow(eq.HEADER)
            w.writerows([[format(val,'.17g') for val in row] for row in rows])
        meta['cases'][name]={'rows':len(rows),'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
    (args.output/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')


if __name__=='__main__':
    main()
