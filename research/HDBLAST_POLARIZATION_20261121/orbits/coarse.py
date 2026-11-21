"""Measure broad energy-direction cells of the collisionless daughter streams."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import runpy
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.integrate import quad
from scipy.special import kve

HERE=Path(__file__).resolve().parent
parent=runpy.run_path(str(HERE/'pairs.py'))
field_table=parent['field_table'];seeds=parent['seeds'];SCRAMBLE=parent['SCRAMBLE']


ENERGY=np.array([0.,1.,2.,3.,5.,np.inf])
COSINE=np.linspace(-1.,1.,5)


def thermal_cells(theta,drift):
    """Lab-count probabilities for bins in the reference stream's rest coordinates."""
    beta=abs(drift)/np.hypot(1.,drift)
    normalization=kve(2,1./theta)
    rows=[]
    for lo,hi in zip(ENERGY[:-1],ENERGY[1:]):
        def density(x):
            return np.exp(-x)*(1.+theta*x)*np.sqrt(theta*x*(2.+theta*x))/normalization
        def current(x):
            return np.exp(-x)*theta*x*(2.+theta*x)/normalization
        density_mass=quad(density,lo,hi,epsabs=1.e-12)[0]
        current_mass=quad(current,lo,hi,epsabs=1.e-12)[0]
        rows.append(.5*density_mass*np.diff(COSINE)+.25*beta*current_mass*np.diff(COSINE**2))
    result=np.array(rows)
    if abs(result.sum()-1.)>2.e-10 or np.min(result)<=0.:raise ValueError('Thermal normalization')
    return result


def compare(h,marker,raw):
    drift=np.interp(marker['time'],h['t'],h['u'])
    boost=np.hypot(1.,drift);theta=marker['theta']
    p=raw[:,5:];gamma=np.hypot(1.,np.linalg.norm(p,axis=1))
    along=p[:,2]*raw[:,0]*np.sign(drift)
    rest_gamma=boost*gamma-abs(drift)*along
    rest_along=boost*along-abs(drift)*gamma
    rest_magnitude=np.hypot(rest_along,np.linalg.norm(p[:,:2],axis=1))
    mu=np.divide(rest_along,rest_magnitude,out=np.zeros_like(along),where=rest_magnitude>0.)
    observed=np.histogram2d(np.maximum(0.,(rest_gamma-1.)/theta),mu,bins=(ENERGY,COSINE))[0]/len(p)
    if abs(observed.sum()-1.)>1.e-14:raise ValueError('Unbinned orbit')
    target=thermal_cells(theta,drift);delta=observed-target
    cut=tuple(map(int,np.unravel_index(np.argmax(delta),delta.shape)))
    return dict(coarse_total_variation=float(.5*np.abs(delta).sum()),
                energy_only_total_variation=float(.5*np.abs(delta.sum(axis=1)).sum()),
                direction_only_total_variation=float(.5*np.abs(delta.sum(axis=0)).sum()),
                observed=observed.tolist(),thermal=target.tolist(),
                largest_excess_cell=list(cut),largest_cell_excess=float(delta[cut]))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('history','fields','scan','pairs','output'):
        parser.add_argument('--'+name,type=Path,required=True)
    parser.add_argument('--jobs',type=int,default=1)
    a=parser.parse_args()
    markers=json.loads((a.history/'markers.json').read_text())['runs']
    runs=json.loads(a.pairs.read_text())['runs']
    inputs={}; cases=[]
    for r in runs:
        if r['case'] in ('uniform_local_pressure_match','pooled_pressure_match','low_cdf_discrepancy_pilot_selection') or (
                r['case']=='distributed_storage' and r['k2']==515 and r['history'].endswith('_half')):
            name=r['history']; label=name.removeprefix('scan_')
            if name not in inputs:
                with np.load(a.history/(label+'.npz')) as z:h=dict(z)
                folder=a.scan if name.startswith('scan_') else a.fields
                with np.load(folder/(label+'.npz')) as z:f=dict(z)
                inputs[name]=(h,f,markers[label]['markers']['0.01'])
            for power,step,seed in [(13,256,SCRAMBLE+19),(12,192,SCRAMBLE+37)]:
                cases.append((r,power,step,seed))
    with tempfile.TemporaryDirectory(prefix='polar-counts-') as wd:
        tmp=Path(wd);exe=tmp/'sphere'
        subprocess.run(['g++','-std=c++17','-O3',str(HERE/'sphere.cpp'),'-o',str(exe)],check=True)
        def run(case):
            r,power,steps,seed=case
            h,f,m=inputs[r['history']]
            t=field_table(h,f,m,r['k1'],r['k2'],r['profile1'],r['profile2'],steps)
            with tempfile.TemporaryDirectory(dir=tmp) as td:
                d=Path(td)
                t.tofile(d/'fields');seeds(h,m,power,seed).tofile(d/'births')
                args=[r[k] for k in ('k1','k2','amplitude1','amplitude2','phase2','second_shape')]
                subprocess.run([str(exe),str(d/'fields'),str(d/'births'),*map(str,args),str(d/'output')],check=True)
                raw=np.fromfile(d/'output',dtype='=f8').reshape(-1,8)
            if len(raw)!=2**(power+1) or not np.isfinite(raw).all():raise ValueError('Missing orbit')
            return dict(case=r['case'],history=r['history'],pairs=2**power,steps=steps,seed=seed,**compare(h,m,raw))
        with ThreadPoolExecutor(a.jobs) as pool:results=list(pool.map(run,cases))
    # A zero drift gives a uniformly distributed polar cosine in every energy bin.
    at_rest=thermal_cells(1.,0.)
    if np.max(np.ptp(at_rest,axis=1))>1.e-14:raise ValueError('Zero-drift target')
    out=dict(kinetic_energy_over_temperature_edges=[0.,1.,2.,3.,5.,'infinity'],
             reference_frame_cosine_edges=COSINE.tolist(),runs=results)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(out,allow_nan=False,sort_keys=True,indent=2)+'\n')


if __name__=='__main__':main()
