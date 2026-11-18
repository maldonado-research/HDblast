"""Run the published B1 RK4 solver with a matched pulse initial condition."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
CAP = HERE.parent
B1 = 'research/HDBLAST_CHECKPOINT_20260930/B1_derived_chi_source/'
PINS = {
    B1+'evolve_b1.py':'a0ae79e12ff90f549a0d1ca674fa3285b8ad456b780880df405b9597d3434f0f',
    B1+'static_w.py':'ac4f5df5d8a1fcedfbd9503caf4717c1473fde6699388e095618f269b6577b1f',
    B1+'shell_matter.py':'58c4c4cc18132f86aaa45ac9f0930490ef3debe160fe8f6ed59274451aadc7d1',
    'research/HDBLAST_CHECKPOINT_20260927/mechanisms/pilot_5d/rolloff5d_v1.py':
        '9ea4d8d6ef36ce3d242fbb0e09f27a944b571fdfb8646e3d71195b60b089080b',
}
# Only the initial deviation vector is supplied to the producer from outside.
BEFORE = '        proj_log = None'
AFTER = '        Yv = self.initial_deviations.copy()\n' + BEFORE
DRIVER = """
import json
import sys
import numpy as np
from pathlib import Path
from evolve_b1 import Run
from shell_matter import ChiGas
from matched_pulse import initial_data

summary,output,h = Path(sys.argv[1]),Path(sys.argv[2]),float(sys.argv[3])
model = json.loads(summary.read_text())
matter = lambda rb: ChiGas(100., .5, 1., 6.49e-6, rb)
setup = dict(delta=.1,c=model['params']['c'],d=model['params']['d'],dc=0.,matter=matter,
             dzf=h,dzc=h,L=.5,xc=np.inf,project=None,kappa=0.,order=4,cfl=.5,log=lambda *x:None)
shared=None
for case in ('quiet','delayed'):
    run = Run(**setup,shells=shared)
    if shared is None:
        shared=(run.bg,run.seed)
        prepared=initial_data(run.z,summary)
    amp=0 if case=='quiet' else 1
    run.initial_deviations=amp*np.array([prepared['delta'][k] for k in ('A','B','phi','A_T','B_T','phi_T')])
    run.evolve(T_final=.18,rec_dT=.00025,tag=case,out_dir=str(output),wall_limit=240)
"""

def compare(repo, output, spacing, constructor):
    import numpy as np
    output.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        for path,digest in PINS.items():
            data=(repo/path).read_bytes()
            if hashlib.sha256(data).hexdigest()!=digest:
                raise ValueError('Historical source mismatch: '+path)
            if path.endswith('/evolve_b1.py'):
                if data.count(BEFORE.encode())!=1:
                    raise ValueError('Unexpected initialization in producer')
                data=data.replace(BEFORE.encode(),AFTER.encode())
            (work/Path(path).name).write_bytes(data)
        shutil.copyfile(constructor,work/'matched_pulse.py')
        (work/'driver.py').write_text(DRIVER)
        rawout = work/'output'
        rawout.mkdir()
        env = dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
        subprocess.run([sys.executable,'-B','driver.py',str(CAP/'B1_summary.json'),str(rawout),str(spacing)],
                       cwd=work,check=True,env=env)
        for case in ('quiet','delayed'):
            meta=json.loads((rawout/(case+'_summary.json')).read_text())
            if meta['stop_reason']!='T_final' or meta['T_end'] < .18-1e-11 or meta['matter']['events']:
                raise ValueError('Incomplete integration or particle event')
            with np.load(rawout/(case+'_timeseries.npz'), allow_pickle=False) as record:
                names=record['cols'].tolist()
                keys=('T','H0tau','phi_b','A_b','B_b','v_over_H0','R','X')
                rows=record['rec'][:,[names.index(k) for k in keys]]
            with (output/(case+'.csv')).open('w',newline='') as f:
                writer=csv.writer(f,lineterminator='\n')
                writer.writerow(keys)
                writer.writerows([format(x,'.17g') for x in row] for row in rows)
            print(case,meta['n_points'],len(rows),meta['rho_b'],meta['phi_b_static'])


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('repo',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--spacing',type=float,required=True)
    p.add_argument('--constructor',type=Path,default=CAP/'controlled/matched_pulse.py')
    options=p.parse_args()
    compare(options.repo.resolve(), options.output.resolve(), options.spacing, options.constructor)
