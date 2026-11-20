"""Classical closed-wall spectator propagation through the two parity-even scalar steps."""
import argparse
from bisect import bisect_right
import hashlib
import json
import math
from pathlib import Path
import runpy
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
mode = runpy.run_path(str(HERE/'helicity.py'))
smooth = mode['smooth']
HARMONICS = mode['HARMONICS']
PREP = (4.0, 4.4)
RELEASE = (6.56, 6.67)


def load_source():
    case_dir = HERE.parent/'exit'
    source = case_dir/'growth_forward.cpp'
    runner = runpy.run_path(str(case_dir/'growth.py'))
    cases = json.loads((case_dir/'cases.json').read_text())
    params = cases['defaults'] | cases['cases']['primary']
    # The existing dense-output mode samples the original solver every .001 proper-time units.
    with tempfile.TemporaryDirectory() as folder:
        exe = Path(folder)/'forward'
        subprocess.run(['g++','-std=c++17','-O3',str(source),'-o',str(exe)], check=True)
        run = runner['run'](exe,params | {'end':6.7,'sample':2})
    return [dict(zip(runner['COLUMNS'],row,strict=True)) for row in run['trace']
            if PREP[0]-1e-5 <= row[0] <= RELEASE[1]+.001],params,hashlib.sha256(source.read_bytes()).hexdigest()


def field_functions(rows,params):
    times = [r['t'] for r in rows]
    def K(psi):
        up = smooth((psi-params['preparation_field'])/params['preparation_width'])
        down = smooth((psi-params['exit_field'])/params['exit_width'])
        k = rows[0]['K']+params['preparation_gauge_increment']*up
        return 1/((1-down)/k + down/(rows[0]['K']+params['final_gauge_increment']))

    def state(t):
        j = min(len(rows)-2,max(0,bisect_right(times,t)-1))
        l,r = rows[j],rows[j+1]
        dt = r['t']-l['t']
        u=(t-l['t'])/dt
        def cubic(y0,y1,d0,d1):
            return y0+u*(dt*d0+u*(3*(y1-y0)-dt*(2*d0+d1)+u*(dt*(d0+d1)-2*(y1-y0))))
        a=cubic(l['a'],r['a'],l['H']*l['a'],r['H']*r['a'])
        psi=cubic(l['psi'],r['psi'],l['pi'],r['pi'])
        return a,K(psi)
    error = max(abs(K(row['psi'])/row['K']-1) for row in rows)
    if error>1e-8 or abs(state(PREP[1])[1]/state(RELEASE[0])[1]-1)>1e-12:
        raise ValueError('The fixed windows or constitutive law disagree with the trajectory.')
    return state,error


def segment(state,window,n):
    start,stop = window
    step = (stop-start)/n
    grid = []
    for i in range(n):
        a,K = state(start+step*(i+.5))
        grid.append((step/a,K))
    return grid


def multiply(l,r):
    a,b,c,d=l
    e,f,g,h=r
    return (a*e+b*g,a*f+b*h,c*e+d*g,c*f+d*h)


def waves(grid,k,kref):
    w=(1.,0.,0.,1.)
    for deta,K in grid:
        co,si=math.cos(k*deta),math.sin(k*deta)
        w=multiply((co,si*kref/K,-si*K/kref,co),w)
    return w


def gain(matrix,kin,kout):
    a,b,c,d=matrix
    ratio=math.sqrt(kout/kin)
    a,b,c,d=a*ratio,b*ratio,c/ratio,d/ratio
    # Singular value formula avoids subtracting nearly equal energies at unit gain.
    largest=(math.hypot(a+d,c-b)+math.hypot(a-d,b+c))/2
    return dict(max_energy_gain=largest**2,phase_averaged_energy_gain=(a*a+b*b+c*c+d*d)/2,
                determinant_minus_one=a*d-b*c-1)


def evaluate(rows,params,n):
    state,law_error = field_functions(rows,params)
    left = segment(state,PREP,n)
    right = segment(state,RELEASE,n)
    flat = segment(state,(PREP[1],RELEASE[0]),n)
    bridge = math.fsum(deta for deta,K in flat)
    k0,km,kf=(state(t)[1] for t in (PREP[0],PREP[1],RELEASE[1]))
    out={}
    for k in HARMONICS:
        first=waves(left,k,k0)
        flat_step=waves([(bridge,km)],k,k0)
        second=waves(right,k,k0)
        total=multiply(second,multiply(flat_step,first))
        # Change reference impedance for the isolated second step to its own initial value.
        a,b,c,d=second
        out[str(k)] = dict(combined=gain(total,k0,kf),preparation=gain(first,k0,km),
                           release=gain((a,b*km/k0,c*k0/km,d),km,kf))
    return dict(modes=out,conformal_preparation=math.fsum(x[0]for x in left),
                conformal_free_interval=bridge,conformal_release=math.fsum(x[0]for x in right),
                K_values=[k0,km,kf],max_K_law_relative_error=law_error)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    rows,params,source_hash=load_source()
    n=16384
    primary=evaluate(rows,params,n)
    half=evaluate(rows,params,n//2)
    sparser=rows[::2]+([rows[-1]] if len(rows)%2==0 else [])
    coarse_source=evaluate(sparser,params,n)
    def difference(other):
        return max(abs(r['max_energy_gain']/other['modes'][k][stage]['max_energy_gain']-1)
                   for k,detail in primary['modes'].items() for stage,r in detail.items())
    worst=max(abs(r['determinant_minus_one']) for detail in primary['modes'].values() for r in detail.values())
    end = next(row for row in rows if abs(row['t']-RELEASE[1]) < 1e-11)
    endpoint = {key:end[key] for key in ('t','a','K','B','X')}
    density = .5*end['K']*(end['B']**2+end['X']**2)
    if abs(density/end['sheet']-1)>1e-12:
        raise ValueError('Endpoint fields and archived sheet density do not agree.')
    endpoint.update(sheet_density=density, comoving_sheet_energy=end['a']**4*density)
    primary.update(endpoint_reference=endpoint,numerical_checks=dict(steps_per_transition=n,source_samples=len(rows),
                                         half_steps_relative_gain_change=difference(half),
                                         half_source_relative_gain_change=difference(coarse_source),
                                         max_abs_wronskian_error=worst),
                   source_solver_sha256=source_hash,preparation_proper_time=PREP,release_proper_time=RELEASE)
    if worst>1e-6 or difference(half)>.01 or difference(coarse_source)>.01:
        raise ValueError(primary['numerical_checks'])
    args.output.write_text(json.dumps(primary,sort_keys=True,indent=2,allow_nan=False)+'\n')


if __name__ == '__main__':
    main()
