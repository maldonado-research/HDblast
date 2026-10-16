#!/usr/bin/env python3
"""Read-only comparison of producer core with independent coordinate oracle."""
from fractions import Fraction as F
from pathlib import Path
from itertools import product
import sys, argparse, json, hashlib
sys.path.insert(0,str(Path(__file__).resolve().parent))
import independent_coordinate_oracle as ref

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'geometry'))
import indexed_geometry as prod

def compare_core(case):
    g,q,phi,eta,m,s,b,point=ref.manufactured(case)
    direct,Gamma,gi=ref.direct(g,q,phi,eta,b)
    own,aux=ref.indexed(g,q,phi,eta,b,Gamma,gi)
    q0={(a,c):q[a][c].v for a,c in ref.PAIRS}
    dq={k:v.v for k,v in aux['dq'].items()}
    metric={(a,c):g[a][c].v for a,c in ref.PAIRS}
    sig=tuple(gi[a][a].v for a in ref.I)
    core=prod.indexed_geometry(q0,dq,aux['ddq'],sig)
    Hphi={(a,c):phi.h[a][c]-sum((Gamma[d,a,c][0].v*phi.d[d] for d in ref.I),F()) for a,c in ref.PAIRS}
    Heta={(a,c):eta.h[a][c]-sum((Gamma[d,a,c][0].v*eta.d[d] for d in ref.I),F()) for a,c in ref.PAIRS}
    source=prod.natural_sources(q0,eta.v,{(a,):eta.d[a] for a in ref.I},Heta,Hphi,
        {(a,):phi.d[a] for a in ref.I},metric,sig,core,b['Vp'],b['Vpp'],b['Vppp'],harmonic=False)
    got={'QRic':core['QRic'],'QE':source['QE'],'QEred_general':source['QEred'],
        'QC_lower_actual':{a:core['QC_actual'][a,] for a in ref.I},
        'QC_upper':{a:sig[a]*core['QC_background'][a,] for a in ref.I},'Qphi_general':source['Qphi']}
    ref.compare(direct,got,'PRODUCER_CORE_VS_DIRECT_COORDINATE')
    controls=[]
    bad=prod.indexed_geometry(q0,dq,aux['ddq'],sig,mutation='ricci_square_half')
    if bad['QRic']==direct['QRic']:raise ValueError('RICCI_FACTOR_CONTROL_NOT_DETECTED')
    controls.append('ricci_square_half')
    badmatter=prod.natural_sources(q0,eta.v,{(a,):eta.d[a] for a in ref.I},Heta,Hphi,
        {(a,):phi.d[a] for a in ref.I},metric,sig,core,b['Vp'],b['Vpp'],b['Vppp'],harmonic=False,mutation='matter_half')
    if badmatter['QE']==direct['QE']:raise ValueError('MATTER_FACTOR_CONTROL_NOT_DETECTED')
    controls.append('matter_half')
    return {'case':case,'status':'PASS_PRODUCER_CORE_VS_INDEPENDENT_COORDINATE_FAMILY','controls':controls}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    if args.output.exists():raise ValueError('OUTPUT_ALREADY_EXISTS')
    results=[compare_core(k) for k in range(2)]
    pins={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
          [Path(__file__),ROOT/'review/independent_coordinate_oracle.py',ROOT/'geometry/indexed_geometry.py']}
    report={'status':'PASS_PRODUCER_CORE_INDEPENDENT_COORDINATE_REVIEW','input_sha256':pins,'cases':results,
        'scope':'Finite exact coordinate families verify implementation; universal identities are independently derived in INDEPENDENT_DERIVATION.md.'}
    args.output.write_text(json.dumps(report,indent=2)+'\n');print(report['status'])

if __name__=='__main__':main()
