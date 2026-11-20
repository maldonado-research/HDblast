"""Build and run a finite-source charged-gas wall study using only Python and C++17."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = Path(__file__).resolve().parent
COLUMNS = "t a H psi pi rad sheet e V source_time emitted returned source_radius K B X mass rest u theta gas pgas random intrinsic_random charge_work drag_heat".split()


def waveform_constants():
    stop = 3 + math.log(39) / 4 + .0625
    def logistic(z):
        val = math.exp(-abs(z))
        return 1/(1+val) if z >= 0 else val/(1+val)
    def step(v):
        if v <= 0: return 0.
        if v >= 1: return 1.
        return logistic(1/(1-v)-1/v)
    def square(t):
        z=4*(t-3)
        return logistic(z)*logistic(-z)**2*step(t/.35)**2*step((stop-t)/.125)**2
    n=32768
    total=math.fsum((2 if j%2==0 else 4)*square(j*stop/n) for j in range(1,n))
    return math.sqrt(2*stop*total/(3*n)),stop


def run(exe, params):
    norm,stop=waveform_constants()
    cells=params['cells']
    late_cells=params.get('late_cells',6)
    args=(cells,.002,params.get('dose',.123),params.get('end',16),norm,stop,
          params.get('cfl',.08),late_cells,.08,13200,int(params.get('trace',False)),
          params['D'],params['Q'],params['nu'],params.get('activation',0.))
    output=subprocess.check_output([str(exe),*map(str,args)],text=True)
    scalars={}
    trace=[]
    for row in output.splitlines():
        name,*values=row.split()
        if name=='@':
            if len(values)!=len(COLUMNS): raise ValueError('Invalid trace length')
            trace.append([float(x) for x in values])
        else:
            scalars[name]=float(values[0])
    if trace: scalars['trace']=trace
    for k,v in scalars.items():
        if k!='trace' and not math.isfinite(v):raise ValueError((k,v))
    return scalars


def compile_solver(exe):
    subprocess.run(['g++','-O3','-std=c++17',str(HERE/'charged_forward.cpp'),'-o',str(exe)],check=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--cases',type=Path,default=HERE/'cases.json')
    parser.add_argument('--jobs',type=int,default=2)
    args=parser.parse_args()
    cases=json.loads(args.cases.read_text())
    with tempfile.TemporaryDirectory(prefix='charged-wall-') as temp:
        exe=Path(temp)/'charged'
        compile_solver(exe)
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            rows=dict(zip(cases,pool.map(lambda name:run(exe,cases[name]),cases)))
    data={'parameters':cases,'trace_columns':COLUMNS,'runs':rows}
    args.output.write_text(json.dumps(data,sort_keys=True,indent=2,allow_nan=False)+'\n')
