"""Compare the kinetic solver to the independent exact plane-wave particle solution."""
import argparse
import json
import math
from pathlib import Path
import runpy
import subprocess
import tempfile

import numpy as np
from scipy.optimize import brentq

HERE=Path(__file__).resolve().parent
run_case=runpy.run_path(str(HERE/'run.py'))['run_case']


def plane_wave(k, amplitude, time, nx):
    strength=1.e7*amplitude/k
    initial=2*math.pi*(np.arange(nx)+.61803398874989485)/nx
    pz=[]
    for start in initial:
        sn=math.sin(start)
        def primitive(x):
            return x+.5*strength**2*((.5+sn**2)*x-.25*math.sin(2*x)+2*sn*math.cos(x))
        phase=brentq(lambda p:primitive(start)-primitive(p)-k*time,start-k*time-1,start)
        pz.append(strength*(sn-math.sin(phase)))
    pz=np.array(pz); px=.5*pz**2;energy=1+px
    return {'gas_conformal':float(np.mean(energy)), 'transverse_rms':float(np.sqrt(np.mean(px**2))),
            'wall_axial_stress_share':float(np.sum(pz**2/energy)/np.sum((px**2+pz**2)/energy))}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args(); k=100.;amp=1.e-5;duration=.038;nx=24
    exact=plane_wave(k,amp,duration,nx)
    with tempfile.TemporaryDirectory() as folder:
        tmp=Path(folder); exe=tmp/'patch'; rows={}
        subprocess.run(['g++','-std=c++17','-O3',str(HERE/'patch.cpp'),'-o',str(exe)],check=True)
        for steps in (96,192,384):
            path=tmp/'constant.txt'
            with path.open('w') as f:
                print(steps,1,1,1.,1.,1.,file=f)
                print(0,1,0,0,0,file=f)
                for _ in range(steps): print(1,1,1,duration/steps,0,0,1,file=f)
            row=run_case(exe,path,(1,1,1,nx,1,steps),k,amp,False,3)
            rows[steps]={'relative_error':{key:abs(row[key]/value-1.) for key,value in exact.items()},
                         'particle_ledger_residual':row['particle_energy_residual_relative']}
    for key in exact:
        if rows[384]['relative_error'][key]>4e-5 or rows[192]['relative_error'][key]>0.30*rows[96]['relative_error'][key]:
            raise ValueError('Plane-wave trajectory convergence failed: '+key)
    args.output.write_text(json.dumps({'analytic':exact,'steps':rows,'duration':duration,'wave_number':k,
                                      'electric_peak':amp},indent=2,sort_keys=True)+'\n')

if __name__=='__main__':
    main()
