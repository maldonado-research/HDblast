"""Write compact per-sample frequency and energy diagnostics for the wider-band controls."""
import argparse
import csv
import io
import json
from pathlib import Path
import zipfile

import numpy as np


def write(output,root,previous):
    with zipfile.ZipFile(output,'w',compression=zipfile.ZIP_DEFLATED) as z:
        for edge in (12,14):
            for grid in ('coarse','fine'):
                p=root/f'K{edge}'/grid/'k32_y1_b1_t0.700'
                Astar=json.loads((previous/'classical/runs'/grid/'quiet/rolling_summary.json').read_bytes())['matter']['events'][0]['A_star']
                sources={'pre':Path(str(p/'fork.npz')+'.history.npz')}
                sources.update({arm:p/arm/'wall_timeseries.npz' for arm in ('quiet','plus','minus','far')})
                for arm,path in sources.items():
                    with np.load(path,allow_pickle=False) as f:
                        raw=dict(zip(f['cols'],f['rec'].T))
                    def field(key):return np.asarray(raw[key],dtype=float)
                    samples=(field('H0tau'),field('A_b')-Astar,field('phi_b'),field('rad_hat'),field('rho_hat'),
                             field('rad')/field('H_over_H0')**2)
                    data=io.StringIO()
                    out=csv.writer(data,lineterminator='\n')
                    out.writerow(('s','log_a','phi','rad_hat','massive_hat','rad_fraction'))
                    for row in zip(*samples,strict=True):out.writerow([format(v,'.16g') for v in row])
                    z.writestr(f'{grid}/K{edge}/{arm}.csv',data.getvalue())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('native_work',type=Path)
    p.add_argument('previous_capsule',type=Path)
    p.add_argument('archive',type=Path)
    args=p.parse_args()
    write(args.archive,args.native_work,args.previous_capsule)
