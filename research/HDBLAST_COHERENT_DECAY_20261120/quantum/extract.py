"""Extract the stable parent trajectories from the pinned #191 release."""
import argparse
import hashlib
import gzip
import json
from pathlib import Path

COLUMNS = 't a H psi pi chi chi_dot'.split()
SOURCE = '17673e63c8e6a7393326f0f623738f4aa532263df0765839e7f177a196670acf'


def extract(path):
    raw=gzip.decompress(path.read_bytes()) if path.suffix=='.gz' else path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=SOURCE:
        raise ValueError('The archive differs from the studied #191 source.')
    d=json.loads(raw)
    ids=[d['trace_columns'].index(key) for key in COLUMNS]
    cases={}
    for name in ('stable','heavy_stable'):
        model=d['runs'][name]; params=d['parameters'][name]
        rows=[row for row in model['trace'] if model['opening_time']-1.5<row[0]<model['maximum_tau']+.15]
        cases[name]=dict(rows=[[r[j] for j in ids] for r in rows],
            opening=model['opening_time'],turn=model['maximum_tau'],mu=params['daughter_mass'],
            initial_ratio=params['initial_mass_ratio'],final_ratio=params['final_mass_ratio'],
            gate=params['decay_threshold'],width=params['decay_width'])
    return dict(source_sha256=SOURCE,source='research/HDBLAST_DRIVEN_PARENT_SUPPLY_20261120/coupled/PARENT.json',
                commit='dfdbad8b4859d331fab1dcfa4a09051d2113d906',
                columns=COLUMNS,cases=cases)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('archive',type=Path);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(extract(args.archive),sort_keys=True,separators=(',',':'),allow_nan=False)+'\n')
