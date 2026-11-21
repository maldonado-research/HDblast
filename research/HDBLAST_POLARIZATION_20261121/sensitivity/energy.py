"""Recheck the timing comparison at nearby source energies and a later wall time."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy
import shlex
import subprocess
import tempfile

import numpy as np

HERE=Path(__file__).resolve().parent
SOURCE=runpy.run_path(str(HERE.parent/'source/backreaction.py'))


def run(parent, calibration, jobs):
    original=parent/'parent/parent/prior'
    cfg=json.loads((original/'prior/envelope/ENVELOPE.json').read_bytes())['parameters']['mass_4e12']
    meta=json.loads((parent/'parent/history/markers.json').read_bytes())
    evolve=runpy.run_path(str(original/'prior/coupled/joined.py'))['run']
    columns=SOURCE['BASE']
    width=len(columns)+3*len(SOURCE['MODE'])
    cases=[(mesh,dose) for mesh in ('native','half','radial') for dose in (.999,1.,1.001)]
    with tempfile.TemporaryDirectory(prefix='pulse-dose-') as wd:
        tmp=Path(wd)
        exe=SOURCE['build'](parent,tmp)
        def case(item):
            mesh,dose=item
            settings=dict(cfg,dose=cfg['dose']*dose,end=11.,sample=0,
                cfl=cfg['cfl']*(.5 if mesh=='half' else 1.),
                early_cells=cfg['early_cells']*(2 if mesh=='radial' else 1),
                late_cells=cfg['late_cells']*(2 if mesh=='radial' else 1))
            with tempfile.TemporaryDirectory(dir=tmp) as td:
                folder=Path(td); data=folder/'trace'; cmd=folder/'wrapped'
                def trial(late,scale):
                    args=[2,meta['action'],original/'boost/BOOST_TABLE.dat',1,1,data,1,.025,.025,7,7,2 if late else 0,scale]
                    cmd.write_text('#!/bin/sh\nexec '+shlex.quote(str(exe))+' "$@" '+shlex.join(list(map(str,args)))+'\n')
                    cmd.chmod(0o755)
                    try: report=evolve(cmd,settings)
                    except subprocess.CalledProcessError as exc: raise RuntimeError(exc.stderr) from exc
                    raw=np.fromfile(data,dtype='=f8').reshape(-1,width)
                    get=lambda k:raw[:,columns.index(k)]
                    t=get('t'); D=get('D'); phi=get('phi')
                    if not np.isfinite(raw).all() or abs(t[-1]-11.)>.00031 or abs(report['end_tau']-11.)>1.e-10:
                        raise ValueError('Run did not reach common time')
                    positive=np.flatnonzero(D>0.)
                    e=np.array([get(k)[-1] for k in ('emitted0','emitted1','emitted2')]); j=int(np.argmax(phi))
                    return dict(emitted=float(e.sum()),second_energy=float(e[-1]),end_radius=float(get('a')[-1]),
                                first_birth=float(t[positive[0]-1]) if positive.size else None,first_birth_by_8=bool(np.any(D[t<=8.])),
                                max_phi=float(phi[j]),time_of_max_phi=float(t[j]),final_phi=float(phi[-1]),
                                final_D=float(D[-1]),source_end=float(t[-1]),radius_cap=cfg['radius_limit'])
                a=trial(False,1.)
                scale=calibration['delayed_matched_'+mesh]['second_scale']*np.sqrt(dose)
                for _ in range(4):
                    b=trial(True,scale)
                    err=(b['emitted']-a['emitted'])/a['emitted']
                    if abs(err)<5.e-11: break
                    scale*=np.sqrt(1+(a['emitted']-b['emitted'])/b['second_energy'])
                else: raise ValueError('Delayed energy match failed')
                print(mesh,dose,a['first_birth'],b['first_birth'],b['max_phi'],flush=True)
                return dict(mesh=mesh,input_dose_scale=dose,both_early=a,one_late=b,matched_energy_relative_error=err)
        with ThreadPoolExecutor(jobs) as pool:return list(pool.map(case,cases))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--parent',type=Path,required=True); parser.add_argument('--source',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True); parser.add_argument('--jobs',type=int,default=2)
    opt=parser.parse_args(); results=run(opt.parent.resolve(),json.loads(opt.source.read_bytes()),opt.jobs)
    opt.output.write_text(json.dumps(dict(runs=results),indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__': main()
