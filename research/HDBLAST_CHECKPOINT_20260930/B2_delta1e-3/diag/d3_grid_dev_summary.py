"""B2 dev D3: grid comparison for the tuned model (d*), delta = 1e-3, Y = 0, dc = 1e-2, old chart, T <= 2.5 (runs/dev/g*).
Output diag/D3_GRID_DEV.json: growth rates of |phi_b - phi_b(static)| in chart time and phi_b, H at T = 2.  Numerical."""
import numpy as np, glob, json
out = {}
for f in sorted(glob.glob('runs/dev/g*_summary.json')):
    s = json.load(open(f)); d = np.load(f.replace('_summary.json', '_timeseries.npz')); cols = list(d['cols']); rec = d['rec']
    ts = {c: rec[:, i] for i, c in enumerate(cols)}
    dev = np.abs(ts['phi_b'] - s['phi_b_static']); T = ts['T']
    def sl(a, b):
        m = (T >= a) & (T <= b); return float(np.polyfit(T[m], np.log(dev[m]), 1)[0]) if m.sum() > 10 else None
    tag = s['tag']
    out[tag] = dict(n_points=s['n_points'], grid=s['params'].get('grid_kind'), dzf=s['params']['dzf'], dzc=s['params']['dzc'], L=s['params']['L'],
                    zfine=s['params']['zfine'], points_per_wall=s.get('points_per_wall'), rate_T0p2_1p0=sl(0.2, 1.0), rate_T0p5_1p5=sl(0.5, 1.5),
                    phi_b_T2=float(np.interp(2.0, T, ts['phi_b'])), H_T2=float(np.interp(2.0, T, ts['H_over_H0'])),
                    lapse_T2=float(np.interp(2.0, T, ts['lapse'])), max_Hnear=float(np.nanmax(ts['Hmax_near'])), T_end=float(T[-1]), runtime_s=s['runtime_s'])
    print(tag, json.dumps(out[tag]))
json.dump(dict(status='numerical', purpose=__doc__, runs=out), open('diag/D3_GRID_DEV.json', 'w'), indent=1)
