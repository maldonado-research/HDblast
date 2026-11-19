"""Validate the native-shell endpoints, integrated energy, controls and sensitivity."""
import argparse
import csv
import io
import json
import math
from pathlib import Path
import subprocess
import sys
import zipfile

HERE=Path(__file__).resolve().parent


def require(ok, description):
    if not ok:
        raise ValueError(description)


def near(x,y,description,rel=2e-9,absolute=3e-11):
    require(math.isclose(x,y,rel_tol=rel,abs_tol=absolute),description)


def verify():
    report=json.loads((HERE/'COUPLED.json').read_text())
    sensitivities=json.loads((HERE/'SENSITIVITY.json').read_text())
    saved=report['cases']
    require(len(saved)==32 and len(sensitivities)==5, 'Incomplete experiment')
    with zipfile.ZipFile(HERE/'WALLS.zip') as z:
        receipts=json.loads(z.read('initial_data.json'))
        for label,ref in saved.items():
            grid,config,arm=label.split('/')
            cut=report['cut_s']
            pre=receipts[grid+'/'+config+'/prepare']
            require(max(abs(pre[k]) for k in ('energy_density','pressure','force'))<1e-14,'Nonempty initial junction')
            require(abs(pre['derivatives'][0])<1e-14,'Initial density derivative')
            for field in ('energy_density','pressure','force','derivatives','s','T'):
                require(report['preparations'][grid+'/'+config+'/prepare'][field]==pre[field], field)
            with z.open(label+'.csv') as f:
                rows=list(csv.DictReader(io.TextIOWrapper(f,encoding='ascii')))
            lo,hi=next((a,b) for a,b in zip(rows,rows[1:]) if float(a['H0tau'])<=cut<=float(b['H0tau']))
            time1,time2=float(lo['H0tau']),float(hi['H0tau'])
            fraction=(cut-time1)/(time2-time1)
            def sample(row,key):
                if key=='rad_fraction':
                    return float(row['rad'])/float(row['H_over_H0'])**2
                return float(row[key])
            for key in ('phi_b','H_over_H0','rad_hat','rho_hat','comoving_number','scalar_work','pressure_work','comoving_energy','rad_fraction'):
                v0,v1=sample(lo,key),sample(hi,key)
                near(v0+(v1-v0)*fraction,ref[key],label+key)
            require(ref['min_uncertainty_through_cut']>=.2499,'Uncertainty inequality')
            require(ref['max_scaled_ward_error']<1e-7,'Integrated Ward error')
            require(ref['max_physical_frequency_through_cut']<report['M5_over_H0'],'Band exceeds admitted frequency')
            if config.startswith('y0'):
                require(ref['rad_hat']==0., 'Stable sector generated radiation')
            if arm!='quiet':
                matched=receipts[label]
                require(abs(matched['matched_delta_y_left'])<2e-13,'Pulse boundary mismatch')
        for g in ('coarse','fine'):
            for branch,sign in [('plus',-1),('minus',1),('far',1)]:
                c=report['contrasts'][g+'/y1/'+branch]['rad_hat']
                for feed,k in [(0,'passive_delta'),(1,'coupled_delta')]:
                    name=g+f'/y1_b{feed}/'
                    d=saved[name+branch]['rad_hat']-saved[name+'quiet']['rad_hat']
                    near(d,c[k],branch+k)
                    require(sign*d>0.,'Pulse sign changed')
    for name,r in sensitivities.items():
        v=r['cases']
        for branch,sgn in [('plus',-1),('minus',1),('far',1)]:
            contrast=100*(v[branch]['rad_hat']/v['quiet']['rad_hat']-1)
            near(contrast,v[branch]['radiation_contrast_percent'],name+branch)
            require(sgn*contrast>0.,'Sensitivity sign changed')
    print('Verified 32 native baseline histories, common preparations, signed sources and 12 sensitivity endpoints.')


if __name__=='__main__':
    par=argparse.ArgumentParser(description=__doc__)
    par.add_argument('--full',action='store_true',help='Regenerate the full NumPy history analysis')
    args=par.parse_args()
    verify()
    if args.full:
        flags=['-O'] if sys.flags.optimize else []
        rebuilt=json.loads(subprocess.check_output([sys.executable,'-I','-B',*flags,HERE/'summarize.py']))
        stored=json.loads((HERE/'COUPLED.json').read_bytes())
        require(rebuilt.keys()==stored.keys() and rebuilt['cases'].keys()==stored['cases'].keys(),'Report structure')
        def compare(a,b):
            if isinstance(a,dict):
                require(a.keys()==b.keys(),'Missing keys')
                for name in a:
                    compare(a[name],b[name])
            elif isinstance(a,list):
                for x,y in zip(a,b,strict=True):
                    compare(x,y)
            elif isinstance(a,float):
                near(a,b,'Numerical report',rel=3e-7,absolute=5e-11)
            else:
                require(a==b,'Changed report')
        compare(rebuilt,stored)
        print('Regenerated the coupled-shell report from the archived walls.')
