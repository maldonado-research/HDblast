"""Check that simple published constants round the exact report outward."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--input',type=Path,required=True)
    ap.add_argument('--input-sha256',required=True)
    ap.add_argument('--out',type=Path,required=True)
    args=ap.parse_args()
    raw=args.input.read_bytes()
    sha=hashlib.sha256(raw).hexdigest()
    if sha!=args.input_sha256: raise ValueError('Exact report hash mismatch')
    data=json.loads(raw)
    rows=[]
    proposals=[('registered_exact_decimal_model',Q(1,10**7),Q(4318),Q(335),Q(89,10**5)),
               ('simple_polynomial_model',Q(1,1000),Q(1021,5000),Q(59,2500),Q(69,5000))]
    for m,(name,du,ku,hu,sl) in zip(data['models'],proposals,strict=True):
        if m['name']!=name: raise ValueError('Model identity mismatch')
        k=Q(m['exact_constants']['K_eta']); h=Q(m['strict_curvature_suppression']['best_K_H'])
        d=Q(m['delta_star']); q=Q(m['strict_curvature_suppression']['q'])
        checks=dict(delta_shrunk=0<du<=d,scalar_constant_rounded_up=ku>=k,
                    curvature_constant_rounded_up=hu>=h,suppression_lower_bound=q-hu*du>sl>0)
        if not all(checks.values()): raise ArithmeticError('Inward-rounded presentation constant')
        rows.append(dict(name=name,delta_readable=str(du),K_eta_upper=str(ku),K_H_upper=str(hu),
                         suppression_coefficient_lower=str(sl),checks=checks))
    out=dict(status='PASS_EXACT_OUTWARD_PRESENTATION',input_sha256=sha,
             source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),rows=rows)
    args.out.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'])


if __name__=='__main__': main()
