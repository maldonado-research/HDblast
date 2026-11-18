#!/usr/bin/env python3
"""Independent finite-collar diagnostic from saved six-field rerun snapshots.

No producer imports or imposed Neumann derivative values are used.
"""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path

import numpy as np

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent/'interior'
FIELDS=('A','B','phi','A_T','B_T','phi_T')
HEADER=('T','s','collar_energy','inner_flux','wall_flux','metric_work','junction_wall_flux')


def require(ok,message):
    if not ok:raise RuntimeError(message)


def weights(z,point,degree,derivative=0,strict_interior=False):
    allowed=np.flatnonzero(z<0) if strict_interior else np.arange(len(z))
    inds=allowed[np.argsort(np.abs(z[allowed]-point),kind='stable')[:degree+1]]
    h=float(np.max(np.abs(z[inds]-point)))
    require(h>0,'zero stencil diameter')
    q=(z[inds]-point)/h
    mat=np.array([q**j for j in range(degree+1)])
    rhs=np.zeros(degree+1);rhs[derivative]=math.factorial(derivative)
    result=np.linalg.solve(mat,rhs)/h**derivative
    require(np.max(np.abs(mat@result-rhs/h**derivative))<1e-8,'stencil polynomial check failed')
    return inds,result


def trap(x,y):
    return math.fsum(float((x[i+1]-x[i])*(y[i]+y[i+1])/2) for i in range(len(x)-1))


def extract(grid,degree,snapshot_dir=None,output_dir=None):
    snapshot_dir=ROOT/f'run_{grid}'/'run/collar' if snapshot_dir is None else snapshot_dir
    output_dir=HERE if output_dir is None else output_dir
    paths=sorted(snapshot_dir.glob('sample_*.npz'))
    require(paths,'no snapshots')
    records=[];pins=[];z0=None;A0=None;mat0=None;mat1=None;zz=None
    for path in paths:
        raw=path.read_bytes()
        pins.append(dict(file=path.name,sha256=hashlib.sha256(raw).hexdigest()))
        with np.load(io.BytesIO(raw),allow_pickle=False) as f:d={k:f[k].copy() for k in f.files}
        z=d['z'];require(np.all(np.diff(z)>0) and z[-1]==0,'invalid spatial grid')
        require(all(np.all(np.isfinite(d[k])) for k in FIELDS),'nonfinite state')
        if z0 is None:
            z0=z.copy();A0=float(d['A'][-1]);require(float(d['T'])==0,'missing initial normalization')
            zz=np.r_[-.04,z[(z>-.04+1e-12)&(z<-1e-12)],0.]
            require(z[0]<-.04 and np.min(np.diff(zz))>1e-8,'collar geometry invalid')
            mat0=np.zeros((len(zz),len(z)));mat1=mat0.copy()
            for j,point in enumerate(zz):
                exact=np.flatnonzero(np.abs(z-point)<1e-12)
                if len(exact):mat0[j,exact[0]]=1
                else:
                    inds,w=weights(z,point,degree);mat0[j,inds]=w
                inds,w=weights(z,point,degree,1,strict_interior=(j==len(zz)-1))
                mat1[j,inds]=w
            require(mat1[-1,-1]==0,'wall derivative used boundary node')
            require(np.max(np.abs(mat1@z-1))<1e-8,'linear derivative witness failed')
        else:require(np.array_equal(z,z0),'spatial grid changed')
        a,b,phi,at,bt,pt=[mat0@d[k] for k in FIELDS]
        pz=mat1@d['phi'];bz=mat1@d['B']
        W=1-phi+phi**3/3;U=.5*(-1+phi**2)**2-(2/3)*W**2
        weight=np.exp(3*(a-A0)-b);potential=np.exp(2*b)*U
        energy=weight*(.5*(pt**2+pz**2)+potential)
        flux=-weight*pt*pz
        work=weight*(.5*(3*at-bt)*pz**2-.5*(3*at+bt)*pt**2+(3*at+bt)*potential+bz*pt*pz)
        inferred=.5*math.exp(3*(float(d['A'][-1])-A0))*(float(d['sigma1'])+float(d['J']))*float(d['phi_T'][-1])
        records.append((float(d['T']),float(d['s']),trap(zz,energy),float(flux[0]),float(flux[-1]),trap(zz,work),inferred))
        # Only the first bracket beyond the fixed window is needed, even for full runs.
        if float(d['s'])>=4:break
    rows=np.asarray(records);T=rows[:,0];clock=rows[:,1]
    require(np.all(np.diff(T)>0) and np.all(np.diff(clock)>0),'invalid time grid')
    require(clock[0]<=1<4<=clock[-1],'rerun has not completed selected window')
    endpoints=np.interp([1.,4.],clock,T)
    tt=np.r_[endpoints[0],T[(T>endpoints[0])&(T<endpoints[1])],endpoints[1]]
    out=np.column_stack([np.interp(tt,T,rows[:,k]) for k in range(len(HEADER))])
    out[:,0]=tt;out[0,1]=1.;out[-1,1]=4.
    folder=output_dir/'collar_data';folder.mkdir(parents=True,exist_ok=True)
    file=folder/f'{grid}_d{degree}.csv'
    with file.open('w',newline='') as f:
        writer=csv.writer(f,lineterminator='\n');writer.writerow(HEADER)
        writer.writerows(tuple(format(float(v),'.17g') for v in row) for row in out)
    receipt=dict(grid=grid,derivative_degree=degree,collar=[-.04,0.],collar_nodes=len(zz),
                 source_nodes=len(z0),source_snapshots=len(pins),A0=A0,
                 selected_T=endpoints.tolist(),selected_s=[1,4],selected_rows=len(tt),
                 public_csv=file.name,csv_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
                 wall_derivative_excludes_boundary_node=True,
                 source_pins=pins)
    (output_dir/f'COLLAR_EXTRACTION_{grid}_d{degree}.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n')
    print(file,len(tt),receipt['csv_sha256'])


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--grid',choices=['coarse','fine'],required=True)
    p.add_argument('--snapshot-dir',type=Path,
                   help='directory containing sample_*.npz (default: ../interior/run_GRID/run/collar)')
    p.add_argument('--output-dir',type=Path,
                   help='root for collar_data/ and extraction receipts (default: this script directory)')
    a=p.parse_args()
    for degree in (2,4):extract(a.grid,degree,a.snapshot_dir,a.output_dir)


if __name__=='__main__':main()
