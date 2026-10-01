#!/usr/bin/env python3
"""Print the restart job line of REGISTRATION section 4: parent run's saved state at the last saved T with shell lapse <= 3.
usage: make_restart_job.py <parent_tag> <out_subdir> <new_tag> <dzf> <dzc> [extra args]"""
import sys, glob, re, json
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
parent, outdir, newtag, dzf, dzc = sys.argv[1:6]; extra = ' '.join(sys.argv[6:])
sub = 'runs/main'
s = json.load(open(HERE/sub/(parent + '_summary.json')))
d = np.load(HERE/sub/(parent + '_timeseries.npz')); cols = list(d['cols']); rec = d['rec']
T = rec[:, cols.index('T')]; lapse = rec[:, cols.index('lapse')]
best = None
for f in sorted(glob.glob(str(HERE/sub/'states'/(parent + '_state_T*.npz')))):
    Ts = float(re.search(r'_state_T([0-9.]+)\.npz', f).group(1))
    lp = float(np.interp(Ts, T, lapse))
    if lp <= 3.0 and (best is None or Ts > best[0]): best = (Ts, f, lp)
p = s['params']
args = '--dstar --Y %g --dc %g --dzf %s --dzc %s --L %g --Tf 20 --kappa %g --xc %g --stop_recollapse --cfl %g --restart %s %s' % (
    p['Y'], p['dc'], dzf, dzc, p['L'], p['kappa'], p['xc'], p['cfl'], str(Path(best[1]).relative_to(HERE)), extra)
print(('%s %s %s' % (outdir, newtag, args)).rstrip())
sys.stderr.write('restart at T=%.2f (parent lapse %.3f)\n' % (best[0], best[2]))
