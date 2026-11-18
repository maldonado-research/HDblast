"""Verify the signed remote-control receipts and optionally evolve all five spectra anew."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

HERE=Path(__file__).resolve().parent


def read(p):
    return json.loads(p.read_bytes())


def require(ok, message):
    if not ok:
        raise ValueError(message)


def close(x,y,label):
    require(math.isfinite(x) and math.isclose(x,y,rel_tol=8e-7,abs_tol=3e-9),label)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--full',action='store_true')
    p.add_argument('--classical',type=Path,default=HERE.parent/'classical')
    args=p.parse_args()
    r=read(HERE/'FAR_RESULTS.json')
    names={'positive':('0.001','plus'),'half':('0.0005','half'),'shift':('shifted','shift')}
    for label,case in r['cases'].items():
        grid,kind=label.split('_')
        branch,stem=names[kind]
        a=read(HERE/'runs'/grid/branch/'initial.json')
        e=read(HERE/'runs'/grid/branch/'rolling_summary.json')
        m=read(HERE/'modes'/f'{stem}-{grid}.json')
        q=read(args.classical/'modes'/f'base-{grid}-quiet.json')
        require(e['stop_reason']=='T_final' and len(e['matter']['events'])==1,label+' native event')
        contact=a['T']-a['control_center']-a['control_halfwidth']
        close(contact,case['earliest_T_compensator'],label+' remote front')
        require(abs(a['matched_delta_y_left'])<1e-10 and a['min_abs_Q']>.9 and a['min_abs_D']>.4,
                label+' compact constraint solve')
        close(100*(m['event_comoving_N_8']/q['event_comoving_N_8']-1),case['event_N8_change_percent'],label+' event change')
        for cut,measurement in case['cuts'].items():
            z=m['ends'][cut]
            close(100*(z['comoving_N_8']/q['ends'][cut]['comoving_N_8']-1),measurement['change_percent'],label+cut)
            require(contact-measurement['wall_T']>.09 and z['band_E_upper_H0']<m['M5_over_H0'] and
                    0<measurement['change_percent']<.055,label+cut+' causal/band check')
    shift=r['shifted_wall_max_abs_difference']
    require(shift['phi_b']<1e-10 and shift['H_over_H0']<1e-8,'relocation equivalence')
    if args.full:
        command=[sys.executable,'-I','-B']
        if sys.flags.optimize:
            command.append('-O')
        actual=json.loads(subprocess.check_output(command+[str(HERE/'report_far.py'),'--classical',
                                   str(args.classical),'--replay-modes']))
        require(actual==r,'history recomputation disagreeing with receipts')
    print('Passed five remote-pulse comparisons and four outgoing cuts (independently reran modes).'
          if args.full else 'Passed remote-pulse receipts and causal margins.')

if __name__=='__main__':
    main()
