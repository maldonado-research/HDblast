"""Check the independent time refinements, bath controls and optional replay."""
import argparse
import json
import math
from pathlib import Path

REF=Path(__file__).resolve().with_name('TRANSFER.json')


def check(rows, expected):
    if set(rows)!=set(expected):
        raise ValueError('Three resolution levels are required')
    bounds={'s_min_R_H':0.001,'s_min_a4_R_H':0.001}
    for level, cases in rows.items():
        if set(cases)!=set(expected[level]):
            raise ValueError('The two histories and three preparation-target pairs are required')
        for name, result in cases.items():
            ref=expected[level][name]
            if set(result)!=set(ref):
                raise ValueError(name)
            for key, target in ref.items():
                val=result[key]
                if type(target) is not float:
                    if target != val:
                        raise ValueError((name,key))
                elif not math.isclose(val,target,abs_tol=bounds.get(key,3e-8),rel_tol=4e-7):
                    raise ValueError((name,key,val,target))
            if result['R_count']<=0 or result['cut_s']!=2.1:
                raise ValueError('Incomplete deposition')
            if name.endswith('_H'):
                if abs(result['R_count']-result['R_H'])>1e-11 or result['min_a4_R_H'] < -1e-10:
                    raise ValueError('Physical bath failed to dissipate passively')
            elif result['min_R_H']>=-.06 or abs(result['R_H']-result['R_count'])<=1e-3:
                raise ValueError('Missing signed loss or distinct final radiation budgets')
    fine=rows['joint_refined']
    medium=rows['time_refined']
    for name,row in fine.items():
        for key in ('R_H','R_count'):
            if abs(row[key]-medium[name][key])>6.1e-5:
                raise ValueError(f'Resolution changed the quoted endpoint {name}/{key}')
    for pair in ('H','B','B_prepared'):
        for key in ('R_H','R_count'):
            contrast=100*(fine['plus_'+pair][key]/fine['quiet_'+pair][key]-1)
            if not -.176<contrast<-.170:
                raise ValueError('Pulse contrast changed')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('report',type=Path,nargs='?',default=REF)
    args=parser.parse_args()
    check(json.loads(args.report.read_bytes()),json.loads(REF.read_bytes()))
    print('Verified independent channel transfers, both preparations and signed radiation minima.')
