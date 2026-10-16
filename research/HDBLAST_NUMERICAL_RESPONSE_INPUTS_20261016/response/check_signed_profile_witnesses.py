#!/usr/bin/env python3
"""Exact consequences of authenticated whole-cell initial-profile enclosures."""
import argparse
from fractions import Fraction as Q
import hashlib
import json
from pathlib import Path


def need(ok, reason):
    if not ok:
        raise ValueError(reason)


def run():
    here = Path(__file__).resolve().parent
    pin = json.loads((here/'PINS.json').read_text())['initial']
    raw = (here/pin['path']).read_bytes()
    need(hashlib.sha256(raw).hexdigest() == pin['sha256'], 'initial_hash')
    data = json.loads(raw)
    need(data['status'] == 'PASS_CONTINUOUS_ACTUAL_SECOND_INITIAL_PROFILE_ENCLOSURES', 'initial_status')
    rows = data['profile_cells']
    need(Q(rows[0]['r'][0]) == 1 and Q(rows[-1]['r'][1]) == 2, 'complete_endpoints')
    need(all(Q(a['r'][1]) == Q(b['r'][0]) for a,b in zip(rows,rows[1:])), 'complete_cell_joins')
    positive = [v for v in rows if Q(v['scalar_mean'][0]) > 0]
    negative = [v for v in rows if Q(v['scalar_mean'][1]) < 0]
    alternating = []
    for row in rows:
        lo,hi = map(Q,row['scalar_mean'])
        sign = 1 if lo > 0 else -1 if hi < 0 else 0
        if sign and (not alternating or alternating[-1]['sign'] != sign):
            alternating.append({'sign':sign,'r':row['r'],
                                'scalar_mean_interval':row['scalar_mean']})
    need(len(alternating) >= 6, 'six_alternating_signed_cells')
    gaps=[]
    for left,right in zip(alternating,alternating[1:]):
        need(Q(left['r'][1]) < Q(right['r'][0]), 'disjoint_zero_bracket')
        gaps.append([left['r'][1],right['r'][0]])
    witnesses = []
    for name,left,right,sign in [('positive','73/64','585/512',1),
                                  ('negative','951/512','119/64',-1)]:
        found = [v for v in rows if list(map(Q,v['r'])) == [Q(left),Q(right)]]
        need(len(found) == 1, 'unique_witness_' + name)
        row = found[0]
        lo,hi = map(Q,row['scalar_mean'])
        if sign > 0:
            need(lo > 7*10**30 and hi < Q('1.4e31'), 'positive_signed_witness')
        else:
            need(hi < -7*10**30 and lo > -Q('1.4e31'), 'negative_signed_witness')
        witnesses.append({'sign':name,'r':[left,right],
                          'certified_scalar_mean_interval':row['scalar_mean'],
                          'magnitude_exceeds':'7000000000000000000000000000000'})
    return {'status':'PASS_ACTUAL_SECOND_INITIAL_SCALAR_MEAN_SIGNED_WITNESSES',
            'initial_sha256':pin['sha256'],
            'whole_cells':len(rows),'positive_cells':len(positive),'negative_cells':len(negative),
            'witnesses':witnesses,
            'ordered_alternating_signed_cells':alternating,
            'disjoint_open_zero_brackets':gaps,
            'distinct_interior_zeros_at_least':len(gaps),
            'consequence':'Continuity and the ordered alternating signed whole cells give a distinct zero in each disjoint open bracket.',
            'scope':'Actual selected initial scalar second parameter derivative, in the registered initial prescription and gauge. This is not an evolved response, force or invariant observation.'}


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    need(not a.output.exists(),'refuse_existing_output')
    report=run()
    a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
    print(report['status'])


if __name__=='__main__':
    main()
