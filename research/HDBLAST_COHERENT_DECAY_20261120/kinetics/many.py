"""Finite many-branch, instantaneous magnetic Dirac-state diagnostic for averaged wall runs."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import runpy
import subprocess
import tempfile

HERE=Path(__file__).resolve().parent
KIN=runpy.run_path(str(HERE/'kinetics.py'))
triangular=KIN['triangular']
FRACTIONS=KIN['FRACTIONS']


def count_range(low,high):
    first=max(1,low)
    return triangular(high)-triangular(first-1) if high>=first else 0


def upper_count(low, high, kmax, groups):
    """Outward-rounded count of an enclosing set of labels, per sign and flavor."""
    b0=max(0,math.floor(low))
    b1=math.ceil(high)
    k=math.ceil(kmax+.5)
    # For an edge one only needs |ell + 1/2 +/- b| <= kmax.
    edges=count_range(b0-k-1,b1+k+1)+count_range(-b1-k-1,-b0+k+1)
    jmax=b1+k-1
    if b0:
        jmax=min(jmax,k*k//(4*b0))
    def labels(j):
        square=k*k-4*j*b0
        if square<0:return 0
        z=math.isqrt(square)
        if z*z<square:z+=1
        # Count BOTH positive-energy states of the interior block whenever either
        # is accessible, thus allowing unitary intrablock mixing.
        return 2*count_range(max(j+1,b0-z),b1+z)
    total=edges
    lower=edges
    if jmax>0:
        step=max(1,(jmax+groups-1)//groups)
        for j in range(1,jmax+1,step):
            n=min(step,jmax-j+1)
            # The enclosing interval can only shrink as j increases.
            total+=n*labels(j)
            lower+=n*labels(j+n-1)
    return total,lower,jmax


def measure(src, runner, parameters, reference, screen):
    anchor='        const double margin=gas.intrinsic_random-.9*wnew.e;'
    if src.count(anchor)!=1:raise ValueError('Missing or ambiguous observer insertion point')
    rest=reference['end_rest']
    extra=(HERE/'readout.inc').read_text().replace('KINETIC_FINAL_COMOVING_REST',repr(rest))
    with tempfile.TemporaryDirectory(prefix='hdblast-many-') as tmp:
        cpp=Path(tmp)/'measure.cpp';exe=Path(tmp)/'measure'
        cpp.write_text(src.replace(anchor,extra+anchor))
        subprocess.run(['g++','-std=c++17','-O3',str(cpp),'-o',str(exe)],check=True)
        d=runner['run'](exe,dict(parameters,sample=0))
    for key in ('end_rest','hot_efolds'):
        if not math.isclose(d[key],reference[key],rel_tol=2e-8):
            raise ValueError(f'Observer changed {key}')
    floor=d['spinor_drag_floor']
    if not math.isclose(floor,screen['elastic_drag_action_floor'],rel_tol=1e-10):
        raise ValueError('Thermal action floor mismatch')
    mu=parameters['daughter_mass'];nu=parameters['drag']
    results=[]
    for j,(f0,f1) in enumerate(zip(FRACTIONS,FRACTIONS[1:])):
        row={k.removeprefix(f'spinor_{j}_'):v for k,v in d.items() if k.startswith(f'spinor_{j}_')}
        neg,pos=row['bmin'],row['bmax']
        low=0. if neg<=0<=pos else min(abs(neg),abs(pos))
        high=max(abs(neg),abs(pos))
        guard=.01*(high-low)+2*row['max_step_change']+8*math.ulp(high)
        low=max(0.,low-guard);high+=guard
        k=1.01*row['kmax']
        up,down,jmax=upper_count(low,high,k,4096)
        up2,down2,_=upper_count(low,high,k,512)
        if not (down2<=down<=up<=up2 or (up-down)/up<.001):
            raise ValueError('Unexpected enclosure discretization')
        cap=6*up
        # Weyl phase volume of a ball of radius k swept along a segment of length high-low.
        # Six charge/flavor groups with spin two on a unit S3 of volume 2*pi^2.
        continuum=4*k**3+3*k*k*(high-low)
        req=2*math.pi**2*(f1-f0)*rest/mu
        ceiling=cap/req
        item=dict(fraction_bin=[f0,f1],start=row['start'],end=row['end'],max_interior_branch_index=jmax,
                  label_upper=cap,weyl_swept_ball_count=continuum,upper_to_weyl_swept_ball=cap/continuum,
                  enclosure_sum_relative_bracket=(up-down)/up,
                  coarser_enclosure_relative_change=(up2-up)/up,
                  no_recycling_S0_ceiling=ceiling,requested_over_capacity_at_drag_floor=floor/ceiling,
                  turn_over_at_drag_rate_S0_ceiling=ceiling*(1+nu*(row['end']-row['start'])))
        if 'witness_action' in screen:
            item['requested_over_capacity_at_screen_witness']=screen['witness_action']/ceiling
        results.append(item)
    return dict(mu=mu,charge_to_mass=parameters['charge_to_mass'],hot_radius_factor=math.exp(d['hot_efolds']),
                pointwise_drag_floor=floor,integrated_drag_floor=nu*d['spinor_hot_duration']/d['spinor_hot_ceiling_integral'],
                original_screen=screen,windows=results,
                lowest_collisionless_S0_ceiling=min(x['no_recycling_S0_ceiling'] for x in results),
                lowest_optimistic_S0_ceiling_if_states_clear_at_drag_rate=min(x['turn_over_at_drag_rate_S0_ceiling'] for x in results))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--envelope',type=Path,default=HERE.parent/'envelope/ENVELOPE.json')
    parser.add_argument('--cases',nargs='+',default=['mass_4e8','mass_4e10'])
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    env=json.loads(args.envelope.read_text())
    coupled=HERE.parent/'coupled'
    runner=runpy.run_path(str(coupled/'joined.py'))
    transform=runpy.run_path(str(HERE.parent/'envelope/build.py'))['transform']
    src=transform((coupled/'joined_forward.cpp').read_text())
    digest=hashlib.sha256(src.encode()).hexdigest()
    if digest!=env['transformed_solver_sha256']:raise ValueError('Averaged solver source mismatch')
    for b in (0.,.13,12.7,120.,1e3):
        for k in (1.51,10.4,37.8):
            assert upper_count(b,b,k,1024)[0]>=KIN['complete_count'](b,k)
    runs={n:measure(src,runner,env['parameters'][n],env['runs'][n],env['screen'][n]) for n in args.cases}
    args.output.write_text(json.dumps({'transformed_solver_sha256':digest,'runs':runs},
                                     sort_keys=True,allow_nan=False,indent=2)+'\n')
