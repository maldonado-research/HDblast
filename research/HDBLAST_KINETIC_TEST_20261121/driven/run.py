"""Drive a gyrotropic relativistic particle ensemble from the accepted wall birth history."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
import os
from pathlib import Path
import struct
import subprocess
import tempfile

import numpy as np
from scipy.integrate import cumulative_trapezoid

HERE=Path(__file__).resolve().parent
FIELDS='t n charge_velocity landau_velocity mean_lab_energy eckart_energy landau_energy eckart_density landau_density replaced_total_density field_change gas_comoving midpoint_work midpoint_redshift injection pair_error hits_per_particle never_hit longitudinal_share proposals electric_work redshift_loss rest_fit_theta energy_cdf_distance rest_longitudinal_share tube_fraction thermal_tube_fraction averaged_landau_fraction averaged_longitudinal_share'.split()
MARKERS=('half_work','nine_tenths_work','reference_first_hot','radius_doubling','turnaround')
SEEDS=(6247,8693,18127,27361)

def prepare(native,keys,channel,count,mesh,path,dt):
    code=channel+'_'+mesh
    meta=keys['runs'][code]
    a=np.memmap(native/(code+'.bin'),dtype=np.float64,mode='r').reshape(-1,len(keys['columns']))
    t,rad,H,X,D,theta,u,R,tried,e,intrinsic,c,d,target=a.T
    opening=meta['opening_time']
    ih=np.searchsorted(t,meta['first_hot_time'])+1
    abs_work=cumulative_trapezoid(np.abs(keys['charge_to_mass']*D[:ih]*X[:ih]*u[:ih]/np.hypot(1.,u[:ih])),t[:ih],initial=0.)
    # Markers use absolute work before first hot, allowing repeated reversals of flow.
    ends=list(np.interp([.5*abs_work[-1],.9*abs_work[-1]],abs_work,t[:ih]))
    ends.extend(meta[k] for k in ('first_hot_time','doubling_time','maximum_tau'))
    t0=opening-4e-5
    steps=[t0];markers=[]
    for end in ends:
        # The short startup has fast electric motion. Later macro steps are bounded
        # separately; all endpoints and both half-step midpoints share native interpolation.
        step=dt if end<=ends[2] else 6*dt
        n=max(1,math.ceil((end-steps[-1])/step))
        steps.extend(np.linspace(steps[-1],end,4*n+1)[1:])
        markers.append(len(steps)-1)
    steps=np.array(steps)
    def interp(x,y):return np.interp(x,t,y)
    aa=interp(steps,rad)
    flux=cumulative_trapezoid(rad*X,t,initial=0.)
    Ii=keys['charge_to_mass']*(interp(steps,flux)-interp(t0,flux))
    Ds=interp(steps,D); us=interp(steps,u); ts=interp(steps,theta); es=interp(steps,e)
    gs=Ds/aa**3*((2.5*ts+np.hypot(1,1.5*ts))*np.hypot(1,us)-ts/np.hypot(1,us))
    rows=np.stack((steps,aa,Ii,Ds,us/np.hypot(1.,us),es,gs,2*interp(steps,flux)/aa**2,interp(steps,X)),axis=1)
    w=meta['end_rest']/count
    mass=(np.arange(count)+.5)*w
    bt=np.interp(mass,D,t)
    RR=interp(bt,R)
    births=np.stack((bt,interp(bt,rad),np.sqrt(np.maximum(0.,(RR-1)*(RR+1))),interp(bt,c),interp(bt,d),interp(bt,target)),axis=1)
    with path.open('wb') as out:
        out.write(struct.pack('=iiidd',len(rows),count,len(markers),w,keys['action']/(6*keys['daughter_mass_L']**3)))
        rows.tofile(out);births.tofile(out);np.array(markers,dtype=np.int32).tofile(out)
    return {name:float(z) for name,z in zip(MARKERS,ends)}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--native',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--jobs',type=int,default=2)
    p.add_argument('--pilot',action='store_true')
    p.add_argument('--large',action='store_true',help='Higher particle count and smaller time steps on both channel choices')
    args=p.parse_args()
    keys=json.loads((args.native/'metadata.json').read_bytes())
    configs=[('baseline',2048,'native',.001),('finer_time',2048,'native',.0005),('more_particles',4096,'half',.0005)]
    if args.large: configs=[('large',8192,'half',.00025),('largest',16384,'half',.00025)]
    if args.pilot: configs=[('pilot',1024,'native',.0015)]
    with tempfile.TemporaryDirectory(prefix='hdblast-driven-') as temp:
        dest=Path(temp); exe=dest/'ensemble'
        subprocess.run(['g++','-O3','-std=c++17',str(HERE/'ensemble.cpp'),'-o',str(exe)],check=True)
        jobs=[]
        anchors={}
        for name,n,mesh,dt in configs:
            for channel in ('dedicated','shared'):
                trace=dest/(name+'_'+channel)
                anchors[channel+'_'+mesh]=prepare(args.native,keys,channel,n,mesh,trace,dt)
                modes=['kinetic']
                if name in ('baseline','pilot'):
                    modes.extend(['fixed_field','no_scattering'])
                for mode in modes:
                    for seed in (SEEDS[:1] if args.pilot else SEEDS[:2] if mode!='kinetic' or args.large else SEEDS):
                        jobs.append((name,channel,mode,seed,trace))
        def job(task):
            name,channel,mode,seed,trace=task
            cmd=[str(exe),str(trace),str(seed),str(int(mode!='fixed_field')),str(int(mode!='no_scattering')),str(1. if channel=='dedicated' else 1./3)]
            table=np.array([[float(x) for x in line.split()] for line in subprocess.check_output(cmd,text=True).splitlines()])
            if table.shape!=(5,len(FIELDS)+1) or not np.isfinite(table).all():
                raise ValueError('Invalid driven ensemble')
            out={}
            for tag,row in zip(MARKERS,table):
                rec=dict(zip(FIELDS,map(float,row[1:])))
                rec['energy_residual']=abs(rec['gas_comoving']-rec['injection']-rec['electric_work']+rec['redshift_loss'])/rec['gas_comoving']
                rec['landau_fraction']=rec['landau_density']/rec['replaced_total_density']
                rec['eckart_fraction']=rec['eckart_density']/rec['replaced_total_density']
                if rec['energy_residual']>2e-6 or rec['pair_error']>1e-7:raise ValueError((tag,rec))
                out[tag]=rec
            print(name,channel,mode,seed,'first_hot',out['reference_first_hot']['landau_fraction'],'doubling',out['radius_doubling']['landau_fraction'],flush=True)
            return (name,channel,mode),out
        data={}
        with ThreadPoolExecutor(args.jobs) as pool:
            for group,out in pool.map(job,jobs):
                data.setdefault('_'.join(group),[]).append(out)
    rows={}
    for group,runs in data.items():
        rows[group]={}
        for marker in MARKERS:
            vals=[r[marker] for r in runs]
            rows[group][marker]={'mean':{k:float(np.mean([v[k] for v in vals])) for k in vals[0]},
                 'standard_error':{k:(float(np.std([v[k] for v in vals],ddof=1)/math.sqrt(len(vals))) if len(vals)>1 else 0.) for k in vals[0]},
                 'seed_range':{k:[min(v[k] for v in vals),max(v[k] for v in vals)] for k in ('eckart_fraction','landau_fraction','never_hit','electric_work','energy_residual','energy_cdf_distance','averaged_landau_fraction','tube_fraction')}}
    args.output.write_text(json.dumps(dict(seeds=SEEDS[:1] if args.pilot else SEEDS,markers=anchors,
            groups=rows,configurations=configs),indent=2,sort_keys=True,allow_nan=False)+'\n')

if __name__=='__main__':
    main()
