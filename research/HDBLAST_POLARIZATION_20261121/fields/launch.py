"""Time-patterned coexact Maxwell waves launched at the outer tube on the #197 geometry."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import shlex
import subprocess
from tempfile import TemporaryDirectory

import numpy as np

HERE = Path(__file__).resolve().parent
KAPPAS = (2,5,7,40,258)
PROFILES = ('early','prepare','release','birth','late_birth','ramp','bipolar')
CHANNELS = tuple((k,f) for k in KAPPAS for f in PROFILES)


def build(parent, temp, channels=CHANNELS, delays=None):
    prior = parent/'parent/parent/prior'
    ancestor = runpy.run_path(str(parent/'geometry/probe.py'))
    ancestor['build'](prior,temp,tuple(k for k,p in channels))
    src=temp/'probe.cpp'; text=src.read_text()

    def change(old,new):
        nonlocal text
        if text.count(old)!=1: raise ValueError('Unexpected parent at '+old[:100])
        text=text.replace(old,new)

    names=list(dict.fromkeys(f for _,f in channels)) if delays is not None else PROFILES
    assignments=','.join(str(names.index(p)) for k,p in channels)
    waveform=(HERE/'waveforms.hpp').read_text()
    if delays is not None:
        offsets=','.join(str(delays[n]) for n in names)
        waveform='double shaped(double tau,int profile,double early) const { constexpr double shift[]={'+offsets+'}; return incoming(tau-shift[profile],dose,norm,stop); }'
    change('constexpr int channels=','constexpr int profile_id[probe_count]={'+assignments+',0};\nconstexpr int channels=')
    change('COUNT=PROBE_DRIVER+probe_count','PROBE_INPUT=PROBE_DRIVER+probe_count, PROBE_SQUARE=PROBE_INPUT+probe_count, PROBE_GRAM=PROBE_SQUARE+probe_count, COUNT=PROBE_GRAM+'+str(len(names)**2))
    change('out["tube_proper_time"] = y[SU];', 'out["tube_proper_time"] = y[SU];\n    for (int j=0;j<'+str(len(names)**2)+';++j) out["wavegram"+std::to_string(j)]=y[PROBE_GRAM+j];')
    change('void probe_cell(', waveform+'\nvoid probe_cell(')
    change('(m==unforced ? 0. : emitted)', '(m==unforced ? 0. : shaped(y[SU],profile_id[m],emitted))')
    change('dy[SU] = g.n;', '''dy[SU] = g.n;
        for(int m=0; m<probe_count; ++m) {
            const double wave = m==unforced ? 0. : shaped(y[SU],profile_id[m],emitted);
            dy[PROBE_SQUARE+m]=g.n*wave*wave;
            dy[PROBE_INPUT+m]=2*r*r*r*dy[PROBE_SQUARE+m];
        }''')
    change('dy[SU] = g.n;', '''dy[SU] = g.n;
        constexpr int nf=PROFILE_NUMBER;
        std::array<double,nf> input_shape;
        for(int f=0;f<nf;++f)input_shape[f]=shaped(y[SU],f,emitted);
        for(int f=0;f<nf;++f)for(int l=0;l<nf;++l)
            dy[PROBE_GRAM+f*nf+l]=2.*g.n*r*r*r*input_shape[f]*input_shape[l];'''.replace('PROFILE_NUMBER',str(len(names))))
    change('gas.comoving_rest,y[PHI]};','gas.comoving_rest,y[PHI],y[SU],y[EMITTED]};')
    change('y[PROBE_NET+m],y[PROBE_DRIVER+m]};','y[PROBE_NET+m],y[PROBE_DRIVER+m],y[PROBE_INPUT+m],y[PROBE_SQUARE+m]};')
    # Retain source arrival sampling through the sharp late feed and the K transition.
    change('t>6.24','t>6.18')
    src.write_text(text); exe=temp/'profiled'
    subprocess.run(['g++','-std=c++17','-O3','-I',str(prior/'boost'),str(src),'-o',str(exe)],check=True)
    return prior,exe


def solve_case(name, prior, exe, history, output, temp, channels=CHANNELS):
    spatial=name=='dedicated_spatial'
    key='dedicated_native' if spatial else name
    meta=json.loads((history/'markers.json').read_bytes())
    spec=meta['runs'][key]
    with np.load(history/(key+'.npz')) as dat:
        ht=dat['t']; baseE=dat['X']; baseB=dat['B']
    cfg=json.loads((prior/'prior/envelope/ENVELOPE.json').read_bytes())['parameters']['mass_4e12']
    trace=temp/(name+'.bin')
    args=['2',str(meta['action']),str(prior/'boost/BOOST_TABLE.dat'),'1',str(spec['channel_fraction']),str(trace)]
    command=temp/name
    command.write_text('#!/bin/sh\nexec '+shlex.quote(str(exe))+' "$@" '+shlex.join(args)+'\n')
    command.chmod(0o755)
    end=spec['markers']['0.05']['time']+.0004
    settings=dict(cfg,end=end,cfl=cfg['cfl']*(.5 if name.endswith('half') else 1.),sample=0,
                  early_cells=cfg['early_cells']*(2 if spatial else 1),late_cells=cfg['late_cells']*(2 if spatial else 1))
    diag=runpy.run_path(str(prior/'prior/coupled/joined.py'))['run'](command,settings)
    n=len(channels)+1
    rows=np.memmap(trace,dtype='=f8',mode='r').reshape(-1,10+7*n)
    t,a,K,X,XB,I,D,phi,tube,emitted=rows[:,:10].T
    E,B,I_mode,net,driver,energy,square=rows[:,10:].reshape(-1,n,7).transpose(2,1,0)
    if np.max(abs(rows[:,-7:]))>0: raise ValueError('Unforced mode changed')
    if not np.allclose(energy[0,-1],emitted[-1],rtol=1e-11): raise ValueError('Emitted energy normalization')
    scale=np.sqrt(emitted[-1]/energy[:-1,-1])
    t0=ht[0]; windows={k:mark['time'] for k,mark in spec['markers'].items()}
    present=(t<=min(t0,t[np.flatnonzero(D>0)[0]-1]))
    original_error=np.max(abs(E[0,present]-X[present]))
    if original_error>2e-9:raise ValueError('Lowest angular mode differs')
    if not spatial and np.max(abs(np.interp(ht,t,X)-baseE))>1e-5:raise ValueError('Background differs')
    traces={}
    for name_f,table,power in [('E',E,1),('B',B,1),('I',I_mode,1),('net',net,2),('driver',driver,2)]:
        traces[name_f]=np.array([s**power*np.interp(ht,t,row) for s,row in zip(scale,table)])
    np.savez_compressed(output/'traces'/(name+'.npz'),t=ht,kappa=[k for k,f in channels],profile=[f for k,f in channels],**traces)
    born_peak=[]
    reference=float(emitted[-1])
    res={}
    for m,(k,f) in enumerate(channels):
        factor=scale[m]; e=E[m]; b=B[m]; incoming=I_mode[m]; work=net[m]+driver[m]
        physical=.5*K*a**4*(e*e+b*b)
        denom=np.max(abs(net[m])+abs(driver[m])+physical)
        err=float(np.max(abs(physical-work))/denom)
        if err>7.e-5:raise ValueError(('Gauge and boundary energy balance',name,k,f,err))
        peak=int(np.argmax(abs(incoming)))
        sample=lambda when,val:float(np.interp(when,t,val))
        markers={w:dict(E=sample(v,e)*factor, B=sample(v,b)*factor,
                     total_stored=sample(v,physical)*factor**2,net_boundary_work=sample(v,net[m])*factor**2,
                     gauge_work=sample(v,driver[m])*factor**2) for w,v in dict(birth=t0,**windows).items()}
        res[f'{k}_{f}']=dict(kappa=k,profile=f,wall_incoming_peak_time=float(t[peak]),
             wall_incoming_peak_magnitude=float(abs(incoming[peak])*factor),
             flux_normalizing_multiplier=float(factor), tube_characteristic_squared_integral=float(square[m,-1]*factor**2),
             comoving_wall_ledger_error_over_peak=err,markers=markers)
    labels=list(dict.fromkeys(f for _,f in channels));nf=len(labels)
    gram=np.array([diag['wavegram'+str(j)] for j in range(nf**2)]).reshape(nf,nf)
    normalizing=np.array([scale[next(j for j,(_,f) in enumerate(channels) if f==label)] for label in labels])
    gram=gram*np.outer(normalizing,normalizing)/reference
    if np.max(abs(gram.T-gram))>1.e-9 or np.max(abs(gram.diagonal()-1.))>2.e-6 or np.linalg.eigvalsh(gram)[0]<1.e-7:
        raise ValueError(('Singular or misnormalized launch Gram',gram))
    return dict(profile_order=labels, emitted_energy_gram=gram.tolist(),first_birth_time=float(t0),tube_proper_time_at_birth=sample(t0,tube),unit_emitted_energy=reference,
                original_zero_mode_prebirth_error=float(original_error),channels=res)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--parent',type=Path,required=True);p.add_argument('--history',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True);p.add_argument('--jobs',type=int,default=2)
    p.add_argument('--scan',action='store_true',help='Scan narrow pulses shifted toward preparation, at higher angular wavenumbers')
    p.add_argument('--cases',nargs='+')
    args=p.parse_args();out=args.output_dir.resolve(); (out/'traces').mkdir(parents=True,exist_ok=True)
    with TemporaryDirectory(prefix='hd-two-pulse-') as raw:
        temp=Path(raw)
        delays=None;channels=CHANNELS
        if args.scan:
            delays={'early':0.,'before':4.989,'t5.085':5.085,'t5.155':5.155,'t5.205':5.205,'t5.225':5.225,
                    't5.235':5.235,'t5.245':5.245,'t5.255':5.255}
            channels=((2,'early'),)+tuple((k,p) for k in (7,20,40,55,72,90,120,160,200,258,360,515)
                                          for p in delays if p!='early')
        prior,exe=build(args.parent.resolve(),temp,channels,delays)
        def one(case):
            val=solve_case(case,prior,exe,args.history,out,temp,channels)
            print('Completed',case,flush=True)
            return case,val
        cases=args.cases or (['dedicated_native','dedicated_half','dedicated_spatial','shared_half'] if args.scan else
                             ['dedicated_native','shared_native','dedicated_half','shared_half','dedicated_spatial'])
        with ThreadPoolExecutor(max_workers=args.jobs) as pool: result=dict(pool.map(one,cases))
    product='STORAGE_SCAN.json' if args.scan else 'FIELDS.json'
    (out/product).write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')

if __name__=='__main__':main()
