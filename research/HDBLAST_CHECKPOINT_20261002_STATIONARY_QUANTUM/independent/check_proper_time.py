#!/usr/bin/env python3
"""Independent direct metric integrals at four preregistered stationary roots."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import proper_time_source as pt


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--roots',type=Path,required=True)
    parser.add_argument('--registration',type=Path,required=True)
    parser.add_argument('--registration-sha256',required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    require(hashlib.sha256(args.registration.read_bytes()).hexdigest()==args.registration_sha256,'Registration hash mismatch')
    require(not args.output.exists(),'Refusing to overwrite prior independent source output')
    roots=json.loads(args.roots.read_text())
    require(roots['status']=='PASS','Independent roots did not pass')
    selected=[r for r in roots['runs'] if r['setting']['rtol']==3e-14 and r['b'] in (-1.,1.) and r['gamma'] in (100.,1e6)]
    require(len(selected)==4,'Exactly four registered independent endpoints required')
    settings=[dict(dps=40,join='.01',upper='60',order=10),dict(dps=55,join='.005',upper='80',order=12)]
    out=dict(started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
             registration_sha256=args.registration_sha256,
             root_input_sha256=hashlib.sha256(args.roots.read_bytes()).hexdigest(),
             source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(pt.__file__))},
             independence='Direct proper-time density from separate metric variation and scalar source; no producer imports.',
             results=[],gates=[],total_error_enclosure=False)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    for root in selected:
        # Scale the full action by the FIXED reference: xbar=x/r,zbar=z/r,rbar=1.
        # Then rho=r²*rhobar and Q=r*Qbar. This is a units change, not a variation.
        with mp.workdps(65):
            ref=mp.mpf(str(root['r']))
            x=mp.mpf(str(root['x']))/ref
            z=mp.mpf(str(root['H2']))/ref
        point=[]
        for setting in settings:
            start=time.monotonic()
            result=pt.sources(x,z,1,setting)
            result.update(b=root['b'],gamma=root['gamma'],seconds=time.monotonic()-start,
                          units='xbar=x/r, zbar=z/r, rhobar=rho/r², Qbar=Q/r')
            point.append(result)
            out['results'].append(result)
            with mp.workdps(60):
                # Trace and tail bounds are independent acceptance gates, not
                # merely fields copied from the inner source calculation.
                rho_value=mp.mpf(result['rho'])
                Q_value=mp.mpf(result['Q'])
                trace=abs(mp.mpf(result['trace_residual']))/z**2
                trace_tolerance=mp.mpf('1e-9')+mp.mpf('1e-7')*max(abs(4*rho_value/z**2),abs(x*Q_value/z**2))
                out['gates'].append(dict(b=root['b'],gamma=root['gamma'],setting=setting,
                    diagnostic='matched_trace',normalized_defect=str(trace),
                    normalized_tolerance=str(trace_tolerance),passed=bool(trace<=trace_tolerance)))
                for index,name,power in [(0,'W',2),(1,'rho',2),(2,'Q',1)]:
                    tail=(mp.mpf(result['spectral_tail_bounds'][index])+mp.mpf(result['IR_tail_bounds'][index]))/z**power
                    tail_tolerance=mp.mpf('1e-9')+mp.mpf('1e-7')*abs(mp.mpf(result[name])/z**power)
                    out['gates'].append(dict(b=root['b'],gamma=root['gamma'],setting=setting,
                        diagnostic='bounded_spectral_plus_IR_tail',source=name,normalized_bound=str(tail),
                        normalized_tolerance=str(tail_tolerance),passed=bool(tail<=tail_tolerance)))
            print(json.dumps({k:result[k] for k in ('b','gamma','settings','rho','Q','seconds')}),flush=True)
            args.output.write_text(json.dumps(out,indent=2)+'\n')
        with mp.workdps(55):
            for name,power in [('rho',2),('Q',1)]:
                expected=mp.mpf(str(root[name]))/ref**power
                actual=mp.mpf(point[1][name])
                normalized_expected=expected/z**power
                normalized_actual=actual/z**power
                tolerance=mp.mpf('1e-9')+mp.mpf('1e-7')*abs(normalized_actual)
                discrepancy=abs(normalized_actual-normalized_expected)
                spread=abs(mp.mpf(point[0][name])-actual)/z**power
                out['gates'].append(dict(b=root['b'],gamma=root['gamma'],source=name,
                    normalized_difference=str(discrepancy), normalized_refinement_spread=str(spread),
                    normalized_tolerance=str(tolerance),passed=bool(discrepancy<=tolerance and spread<=tolerance)))
            rho=mp.mpf(point[1]['rho']); W=mp.mpf(point[1]['W'])
            quantum_tolerance=mp.mpf('1e-9')+mp.mpf('1e-7')*abs(rho/z**2)
            for name,wrongrho in [('density_as_action',W),('metric_curvature_sign',2*W-rho)]:
                difference=abs(wrongrho-rho)/z**2
                out['gates'].append(dict(b=root['b'],gamma=root['gamma'],negative_control=name,
                                       normalized_defect=str(difference),passed=bool(difference>quantum_tolerance)))
        args.output.write_text(json.dumps(out,indent=2)+'\n')
    out['status']='PASS' if all(g['passed'] for g in out['gates']) else 'FAIL'
    out['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    args.output.write_text(json.dumps(out,indent=2)+'\n')
    require(out['status']=='PASS','Independent proper-time gates failed')


if __name__=='__main__':
    main()
