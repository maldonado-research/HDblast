"""Weak finite-band fermion response near the short high-mass threshold crossing."""
import argparse
import json
from pathlib import Path
import runpy

HERE=Path(__file__).resolve().parent
solver=runpy.run_path(str(HERE/'fermions.py'))


def calculate(doc,case):
    bg=solver['Parent'](doc,case);t=bg.opening
    runs={}
    for label,lo,hi,ramp in [('crossing',t-.035,t+.035,.0125),('longer',t-.12,t+.50,.045),
                              ('gentler',t-.12,t+.50,.085),('closed',t-.12,t-.05,.0125)]:
        if case=='prior_400' and label not in ('crossing','longer'):
            continue
        n=12288 if label in ('crossing','closed') else 6144
        dt=.000125
        check=label in ('crossing','longer') and case=='matched_400'
        value=solver['shells'](bg,lo,hi,ramp,n,dt,lam=1.,nonlinear=check)
        if check:
            fine=solver['shells'](bg,lo,hi,ramp,n,dt/2)
            wide=solver['shells'](bg,lo,hi,ramp,2*n,dt)
            value.update(halved_step_relative_change=fine['pair_response_per_lambda2']/value['pair_response_per_lambda2']-1,
                         doubled_band_relative_change=wide['pair_response_per_lambda2']/value['pair_response_per_lambda2']-1)
        if check and (abs(value['halved_step_relative_change'])>.005 or abs(value['doubled_band_relative_change'])>.05):
            raise ValueError(f'Band or step sensitivity in {label}: {value}')
        runs[label]=dict(start=lo,end=hi,ramp=ramp,shell_count=n,**value)
    return dict(mu=bg.mu,opening=t,K_late=doc['cases'][case]['K_late'],relative_history_error=bg.max_history_error,windows=runs)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--data',type=Path,default=HERE/'NARROW_HISTORY.json')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    doc=json.loads(args.data.read_text())
    result={case:calculate(doc,case) for case in doc['cases']}
    args.output.write_text(json.dumps(result,
                                      sort_keys=True,indent=2,allow_nan=False)+'\n')
