"""Count birth-accessible Dirac states and test a single thermal elastic s-wave channel."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
FRACTIONS=(.01,.1,.25,.5,.75,.9,.99)


def triangular(n):
    return n*(n+1)//2


def shell_subspace(ell, b):
    # Direct diagonalization of each uncoupled edge and spin-orbit 2x2 block.
    eigen=[ell+.5-b,ell+.5+b]
    spin=(ell-1)/2
    for twice_m in range(1-ell,ell-1,2):
        m=twice_m/2
        a=2*m+1.5-b
        d=-2*(m+1)+1.5+b
        off=2*math.sqrt((spin-m)*(spin+m+1))
        root=math.hypot((a-d)/2,off)
        eigen.extend(((a+d)/2-root,(a+d)/2+root))
    return eigen


def complete_count(b, cutoff):
    """One charge and flavor. Each admissible intrinsic eigenspinor gives one positive 4-spinor."""
    shift=abs(b)
    def count_ell(lo,hi):
        first=max(1,math.ceil(lo));last=math.floor(hi)
        return triangular(last)-triangular(first-1) if last>=first else 0
    edge=sum(count_ell(-cutoff-.5+s,cutoff-.5+s) for s in (shift,-shift))
    if shift==0.:
        n=math.floor(cutoff-1.5)
        return 2*(n+1)*(n+2)*(n+3)//3 if n>=0 else 0
    total=edge
    for lim in (cutoff+.5, cutoff-.5):
        if lim<0:continue
        top=math.floor(min(lim**2/(4*shift),shift+lim-1))
        for j in range(1,top+1):
            half=math.sqrt(max(0.,lim**2-4*j*shift))
            total+=count_ell(max(j+1,shift-half),shift+half)
    return total


def checks():
    tests=0
    for b in (0.,.13,1.,8.7,12.15):
        for cut in (.23,1.51,2.46,7.23,17.5):
            direct=sum(ell for ell in range(1,45) for val in shell_subspace(ell,b) if abs(val)<=cut)
            assert direct==complete_count(b,cut),(b,cut,direct,complete_count(b,cut))
            tests+=1
    return tests


def quantify(readout, p, final_rest):
    mu=p['daughter_mass'];nu=p['drag'];point=readout['spinor_drag_floor']
    duration=readout['spinor_hot_duration']
    integrated=nu*duration/readout['spinor_hot_ceiling_integral']
    summary=[]
    for j,(f0,f1) in enumerate(zip(FRACTIONS,FRACTIONS[1:])):
        row={k.removeprefix(f'spinor_{j}_'):v for k,v in readout.items() if k.startswith(f'spinor_{j}_')}
        lo,hi=row['bmin'],row['bmax']
        if lo<=0<=hi: raise ValueError('Magnetic sign changes in selected production bin')
        low,high=sorted((abs(lo),abs(hi)))
        # Enlarge the range by two native increments and 1% of its measured span,
        # and allow the radius of the birth shell to exceed the recorded maximum by 1%.
        guard=2*row['max_step_change']+.01*(high-low)
        low=max(0.,low-guard);high+=guard
        k=1.01*row['kmax']
        gap=math.sqrt(4*low)-.5
        if gap<=k or low<k+1:raise ValueError('More magnetic branches can be occupied')
        first=max(1,math.ceil(low-.5-k));last=math.floor(high-.5+k)
        cap=6*(triangular(last)-triangular(first-1)) # both signs, three Dirac flavors
        requested_per_action=2*math.pi**2*(f1-f0)*final_rest/mu
        ceiling=cap/requested_per_action
        elapsed=row['end']-row['start']
        summary.append(dict(fraction_bin=[f0,f1],**row,magnetic_guard=guard,conservative_gap_over_birth_cut=gap/k,
            first_ell=first,last_ell=last,count_all_six=cap,requested_daughters_per_action=requested_per_action,
            no_recycling_action_ceiling=ceiling,
            required_over_available_at_pointwise_drag_floor=point/ceiling,
            minimum_clearance_rate_at_pointwise_drag_floor=max(0.,point/ceiling-1)/elapsed,
            clearance_limited_action_ceiling_at_fluid_rate=ceiling*(1.+nu*elapsed)))
    return dict(rest_target=final_rest,drag_target=nu,hot_low_drift_duration=duration,
        pointwise_swave_action_floor=point,hot_integrated_swave_action_floor=integrated,
        pointwise_floor_witness={k.removeprefix('spinor_drag_floor_'):v for k,v in readout.items() if k.startswith('spinor_drag_floor_')},
        delivery_windows=summary,largest_window_occupancy_deficit=max(row['required_over_available_at_pointwise_drag_floor'] for row in summary))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study',type=Path,default=HERE.parent/'coupled/JOINED.json')
    parser.add_argument('--cases',nargs='+',default=['primary','fine'])
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    joined=HERE.parent/'coupled'
    runner=runpy.run_path(str(joined/'joined.py'))
    source=(joined/'joined_forward.cpp').read_text()
    anchor='        const double margin=gas.intrinsic_random-.9*wnew.e;'
    if source.count(anchor)!=1:raise ValueError('Unexpected diagnostic insertion point')
    study=json.loads(args.study.read_text())
    if hashlib.sha256(source.encode()).hexdigest()!=study['solver_sha256']:
        raise ValueError('Coupled data and solver differ')
    runs={}
    with tempfile.TemporaryDirectory(prefix='hdblast-landau-') as scratch:
        cpp=Path(scratch)/'measure.cpp';exe=Path(scratch)/'measure'
        for name in args.cases:
            p=study['parameters'][name]; rest=study['runs'][name]['end_rest']
            extra=(HERE/'readout.inc').read_text().replace('KINETIC_FINAL_COMOVING_REST',repr(rest))
            cpp.write_text(source.replace(anchor,extra+anchor))
            subprocess.run(['g++','-std=c++17','-O3',str(cpp),'-o',str(exe)],check=True)
            obs=runner['run'](exe,dict(p,sample=0))
            if not math.isclose(obs['end_rest'],rest,rel_tol=3e-8):raise ValueError('Changed reference number')
            if not math.isclose(obs['hot_efolds'],study['runs'][name]['hot_efolds'],rel_tol=3e-7):
                raise ValueError('Changed hot interval')
            runs[name]=quantify(obs,p,rest)
    args.output.write_text(json.dumps({'magnetic_spectrum_count_checks':checks(),'solver_sha256':study['solver_sha256'],
                                     'runs':runs},allow_nan=False,sort_keys=True,indent=2)+'\n')

if __name__=='__main__':
    main()
