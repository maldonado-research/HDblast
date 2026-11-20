"""Continue the pinned #191 high-mass case and sample a short interval round its opening."""
import argparse
import hashlib
import json
from pathlib import Path
import runpy
import subprocess
import tempfile


def run(directory):
    case_doc=json.loads((directory/'narrow_cases.json').read_text())
    cfg=case_doc['defaults']|case_doc['cases']['fine']|dict(continue_stable=1,end=6.755,sample=1)
    src=(directory/'production_forward.cpp').read_bytes()
    digest=hashlib.sha256(src).hexdigest()
    if digest!='dd305f17010ce54a59799db07e51060e0d299f880a6226ed4ad1f51e0f4fe600':
        raise ValueError('Unrecognized #191 parent evolution.')
    old='int samples_per_time=(sampling==2?1000:50);'
    new='int samples_per_time=5000;'
    if src.decode().count(old)!=1:raise ValueError('Expected one dense-output site.')
    # Restrict dense output to this calculation; integration is unaffected.
    tmp=src.decode().replace(old,new).replace('for(int i=4*samples_per_time;i<=13*samples_per_time;i++)',
                          'for(int i=30400;i<=33775;i++)')
    runner=runpy.run_path(str(directory/'production.py'))
    with tempfile.TemporaryDirectory() as folder:
        altered=Path(folder)/'sample.cpp'; altered.write_text(tmp)
        exe=Path(folder)/'forward'
        subprocess.run(['g++','-O3','-std=c++17',str(altered),'-o',str(exe)],check=True)
        cases={}
        columns='t a H psi pi chi chi_dot'.split()
        ids=[runner['COLUMNS'].index(k) for k in columns]
        for key,K in [('prior_400',26.),('matched_400',2.5)]:
            data=runner['run'](exe,cfg | {'final_gauge_increment':K})
            o=data['opening_time']
            rows=[[r[k] for k in ids] for r in data['trace'] if o-.14<=r[0]<=o+.51]
            rows=list({r[0]:r for r in rows}.values())
            if max(b[0]-a[0] for a,b in zip(rows,rows[1:]))>.000201:
                raise ValueError('Missing dense wave samples.')
            cases[key]=dict(rows=rows,opening=o,turn=None,mu=400.,K_late=.04+K,
                initial_ratio=.95,final_ratio=1.05,gate=.53,width=.003)
    return dict(columns=columns,source_solver_sha256=digest,cases=cases)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('previous_coupled_directory',type=Path);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(run(args.previous_coupled_directory),sort_keys=True,
                                      separators=(',',':'),allow_nan=False)+'\n')
