#!/usr/bin/env python3
"""Independent stereographic coordinate check of all producer frame slots."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import sys, argparse, json, hashlib
sys.path.insert(0,str(Path(__file__).resolve().parent))
import independent_coordinate_oracle as ref

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'geometry'))
import geometric_sources as prod

def transition(point,b):
    den=1+sum(x*x for x in point)
    n=[2*x/den for x in point]+[(2-den)/den]
    matrix=[[F() for _ in range(6)] for _ in ref.I]
    matrix[0][0]=F(1);matrix[1][1]=b['rho']
    for i in range(3):
        for j in range(3):matrix[i+2][j+2]=b['rho']*b['c']*(2*(i==j)/den-4*point[i]*point[j]/den**2)
        matrix[i+2][5]=b['rho']*b['c']*(-4*point[i]/den**2)
    return n,matrix

def evaluate(poly,n):
    return sum((coef*product_scalar(n[i]**power for i,power in enumerate(exponent)) for exponent,coef in poly.items()),F())
def product_scalar(values):
    out=F(1)
    for value in values:out*=value
    return out

def pullback(tensor,n,matrix):
    rank=len(next(iter(tensor)))
    work={idx:evaluate(poly,n) for idx,poly in tensor.items()}
    dims=[range(6)]*rank
    for slot in range(rank):
        dims[slot]=range(5);new={}
        for idx in product(*dims):
            val=F()
            for old_index in range(6):
                weight=matrix[idx[slot]][old_index]
                if weight:
                    old=list(idx);old[slot]=old_index
                    val+=weight*work[tuple(old)]
            new[idx]=val
        work=new
    return work

def compare_frame(case):
    g,q,phi,eta,m,s,b,point=ref.manufactured(case)
    direct,Gamma,gi=ref.direct(g,q,phi,eta,b)
    own,aux=ref.indexed(g,q,phi,eta,b,Gamma,gi)
    sources=prod.evaluate_sources(m,s,b,harmonic=False,include_jets=True)
    n,matrix=transition(point,b)
    fields={key:pullback(sources[key],n,matrix) for key in ('QRic','QE','QEred_general','QC_actual','QC_background','q','dq','ddq','deta','ddeta','dQC')}
    expected_jets={
        'q':{(a,c):q[a][c].v for a,c in ref.PAIRS},
        'dq':{idx:v.v for idx,v in aux['dq'].items()},'ddq':aux['ddq'],
        'deta':{(a,):eta.d[a] for a in ref.I},
        'ddeta':{(a,c):eta.h[a][c]-sum((Gamma[d,a,c][0].v*eta.d[d] for d in ref.I),F()) for a,c in ref.PAIRS},
        'dQC':aux['DQC_lower_actual']}
    ref.compare(expected_jets,fields,'FRAME_JETS_TO_COORDINATE_JETS')
    got={'QRic':fields['QRic'],'QE':fields['QE'],'QEred_general':fields['QEred_general'],
        'QC_lower_actual':{a:fields['QC_actual'][a,] for a in ref.I},
        'QC_upper':{a:gi[a][a].v*fields['QC_background'][a,] for a in ref.I},
        'Qphi_general':evaluate(sources['Qphi_general'],n)}
    ref.compare(direct,got,'FULL_PRODUCER_FRAME_VS_DIRECT_COORDINATE')
    # Both independent coordinate and frame constructions expose nonzero
    # gauge constraints. Their agreement therefore exercises general 2PC1.
    if all(not p for p in sources['C1'].values()):raise ValueError('GENERAL_FRAME_C1_MISSING')
    return {'case':case,'status':'PASS_FULL_FRAME_AND_SOURCE_AGAINST_COORDINATE_FAMILY',
        'verified':'q,dq,ddq,deta,ddeta,dQC,QRic,QE,QEred_general,QC_actual,QC_background,Qphi_general'}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    if args.output.exists():raise ValueError('OUTPUT_ALREADY_EXISTS')
    cases=[compare_frame(k) for k in range(2)]
    paths=[Path(__file__),ROOT/'review/independent_coordinate_oracle.py',ROOT/'geometry/indexed_geometry.py',ROOT/'geometry/geometric_sources.py']
    out={'status':'PASS_INDEPENDENT_COMPLETE_FRAME_SOURCE_REVIEW','cases':cases,
        'sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        'scope':'Finite exact genuine warped metric families, never a claim of proof from angular or spacetime samples.'}
    args.output.write_text(json.dumps(out,indent=2)+'\n');print(out['status'])

if __name__=='__main__':main()
