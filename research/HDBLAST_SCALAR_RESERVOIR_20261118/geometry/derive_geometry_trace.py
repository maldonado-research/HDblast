#!/usr/bin/env python3
"""Private extraction of two limited B1 geometry/matter comparison tables."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import subprocess

import numpy as np

HERE = Path(__file__).resolve().parent
PUBLIC_COMMIT = 'e1d1285c65675899b9a1e3958171c11a8fbab961'
ARCHIVE_COMMIT = '140a0e9c672b03e930400d5109213f662842bdf4'
BASE = 'research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source'
HEADER = ('x', 'tension_energy', 'tension_expansion', 'tension_work',
          'matter_energy', 'matter_expansion', 'matter_work')


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def git(repo, *args):
    return subprocess.check_output(['git', '-C', str(repo), *args], stderr=subprocess.PIPE)


def raw_pin(path):
    data = (ARCHIVE/path).read_bytes()
    blob = git(ARCHIVE, 'rev-parse', ARCHIVE_COMMIT+':'+path).decode().strip()
    require(data == git(ARCHIVE, 'cat-file', 'blob', blob), 'archive byte mismatch: '+path)
    return dict(path=path, sha256=sha(data), git_blob=blob, commit=ARCHIVE_COMMIT)


def main():
    (HERE/'data').mkdir(exist_ok=True)
    sources=[]
    source_paths = [BASE+'/'+s for s in ('evolve_b1.py','shell_matter.py','static_w.py','README.md','AUDIT.md')]
    source_paths += [
        'hdblast/checkpoints/2026-09-22_scalar-profile-branch-and-matter_zenodo-22922928/HDBLAST_SCALAR_PROFILE_BRANCH_AND_MATTER_20260922/matter/MATTER_EXTENSION_AND_RESIDUAL_VACUUM.md',
        'research/HDBLAST_CHECKPOINT_20261002_JUNCTIONS/ACTION_AND_JUNCTION_SPECIFICATION.md']
    for rel in source_paths:
        data=(ARCHIVE/rel).read_bytes()
        require(data == git(PUBLIC,'show',PUBLIC_COMMIT+':'+rel), 'public source mismatch: '+rel)
        sources.append(dict(path=rel,sha256=sha(data),bytes=len(data),
                            url='https://github.com/maldonado-research/HDblast/blob/'+PUBLIC_COMMIT+'/'+rel))
    (HERE/'SOURCES.json').write_text(json.dumps(dict(public_commit=PUBLIC_COMMIT,sources=sources),indent=2,sort_keys=True)+'\n')
    receipt=dict(scope='Private provenance for limited derived tables; no bulk reservoir balance is measured.',runs=[])
    for spacing,label in (('1e-3','coarse'),('5e-4','fine')):
        tag='main_ps0.5_G100_y1_bcut_dzf'+spacing
        prefix=BASE+'/runs/main/'+tag
        summary=json.loads((ARCHIVE/(prefix+'_summary.json')).read_text())
        for field,file in (('code_sha256','evolve_b1.py'),('matter_sha256','shell_matter.py'),('static_w_sha256','static_w.py')):
            require(summary[field] == sha((ARCHIVE/BASE/file).read_bytes()), 'recorded producer mismatch')
        with np.load(ARCHIVE/(prefix+'_timeseries.npz'),allow_pickle=False) as z:
            d={str(c):z['rec'][:,i].copy() for i,c in enumerate(z['cols'])}
            inventory={k:dict(shape=list(z[k].shape),dtype=str(z[k].dtype)) for k in z.files}
        x=d['H0tau']; require(np.all(np.diff(x)>0), 'nonincreasing clock')
        # Same predeclared [1,4] interval as the contemporaneous matter-ledger audit.
        lo,hi=1.,4.; require(x[0]<=lo<hi<=x[-1], 'unavailable window')
        bad=np.flatnonzero((np.nan_to_num(d['Hmax_near'],nan=0)>.05)|(np.nan_to_num(d['Mmax_near'],nan=0)>.05))
        end=int(bad[0])-1 if len(bad) else len(x)-1
        require(hi<=x[end], 'window beyond original reliable prefix')
        for k in ('ln_a','sigma','sigma1','H_over_H0','v_over_H0','J','R','X','P'):
            require(np.all(np.isfinite(d[k][(x>=lo)&(x<=hi)])), 'nonfinite selected field')
        a3=np.exp(3*d['ln_a']);h=d['H_over_H0'];v=d['v_over_H0']
        arrays=(a3*d['sigma'],3*h*a3*d['sigma'],a3*d['sigma1']*v,
                a3*(d['R']+d['X']),3*h*a3*d['P'],a3*d['J']*v)
        xx=np.r_[lo,x[(x>lo)&(x<hi)],hi]
        vals=[np.interp(xx,x,a) for a in arrays]
        path=HERE/'data'/('b1_'+label+'.csv')
        with path.open('w',newline='') as f:
            w=csv.writer(f,lineterminator='\n');w.writerow(HEADER)
            w.writerows(tuple(format(float(q),'.17g') for q in row) for row in zip(xx,*vals))
        # Independent evaluation of the explicit archived tension law, not a producer import.
        phi=d['phi_b'];p=summary['params'];W=1-phi+phi**3/3
        sig=2*W+p['delta']*(1+p['c']*phi+p['d']*phi**2/2)
        sig1=2*(-1+phi**2)+p['delta']*(p['c']+p['d']*phi)
        select=(x>=lo)&(x<=hi)
        require(np.max(np.abs(sig[select]-d['sigma'][select]))<2e-14,'tension value mismatch')
        require(np.max(np.abs(sig1[select]-d['sigma1'][select]))<2e-14,'tension derivative mismatch')
        event=summary['matter']['events'][0]
        receipt['runs'].append(dict(tag=tag,npz=raw_pin(prefix+'_timeseries.npz'),
            summary=raw_pin(prefix+'_summary.json'),rho_b=summary['rho_b'],b=summary['matter']['b'],
            csv=dict(path='data/'+path.name,sha256=sha(path.read_bytes()),rows=len(xx)),
            interval=[lo,hi],original_reliable_end=float(x[end]),inventory=inventory,
            columns=list(d),explicit_tension_value_max_difference=float(np.max(np.abs(sig[select]-d['sigma'][select]))),
            explicit_tension_derivative_max_difference=float(np.max(np.abs(sig1[select]-d['sigma1'][select]))),
            event=event,unpaid_jump_over_reference_crossing_energy=event['kappa5sq_rho_jump_over_H0']/event['massless_energy_proxy_over_H0']))
    (HERE/'PRIVATE_EXTRACTION.json').write_text(json.dumps(receipt,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path,help='Versioned archive checkout holding the original NPZ files')
    parser.add_argument('public',type=Path,help='Public HDblast checkout')
    args=parser.parse_args()
    ARCHIVE, PUBLIC=args.archive,args.public
    main()
