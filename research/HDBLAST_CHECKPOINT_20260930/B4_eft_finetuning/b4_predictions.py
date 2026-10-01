#!/usr/bin/env python3
"""B4 step 3 (REGISTRATION.md section 3.1-3.2): reduced-model Y scan and tuning window in d.
The registered reduced model was NOT VALIDATED (B4_VALIDATION.json), so every number here is CONDITIONAL on the reduced
model; r(Y) from the registered matching is reported but is known to be ~12x too large at Y = 1 and is not a prediction.
The radiation-only variant (post-registration) is used for the window because it passed the timing targets T5-T8 and T3.
Output: B4_PREDICTIONS.json."""
import json, math, hashlib, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import eft_reduced as M

def rd_efolds(res, an_tau_m):
    """longest phase-II interval (after the match) with Omega_r >= 0.9, |Omega_vac| <= 0.1 and H > 0 (radiation dominated)."""
    D = res['D']; m0 = D['H0tau'] >= an_tau_m
    H = D['H_over_H0'][m0]; lna = D['ln_a'][m0]
    with np.errstate(divide='ignore', invalid='ignore'):
        om_r = D['rad'][m0]/H**2; om_v = D['vac'][m0]/H**2
    ok = (om_r >= 0.9) & (np.abs(om_v) <= 0.1) & (H > 0)
    best = 0.0; st = None
    for i in range(len(ok)):
        if ok[i] and st is None: st = i
        if (not ok[i] or i == len(ok) - 1) and st is not None:
            en = i if ok[i] else i - 1; best = max(best, float(lna[en] - lna[st])); st = None
    return best

def one(e, Y, q, variant, dc=1e-2):
    res = M.run(e, Y, q*M.DSTAR, dc, 0.95, variant=variant, tmax_E=3000.0)
    if not res['ok']: return dict(ok=False, msg=res.get('msg'))
    an = M.analyse(res); mt = res['match']
    return dict(ok=True, H0tau_match=mt['H0tau'], h_match=mt['h'], R_over_H0sq_match=mt['R']/mt['H0_EFT']**2, rad_match=mt['rad'],
                r_match=mt['r'], lna_match=mt['lna'], classification=an['classification'], plateau=an['plateau'],
                H_zero=an['H_zero_H0tau'], recollapse=an['recollapse_H0tau'], rad_era=an['rad_era'],
                efolds_RD=rd_efolds(res, mt['H0tau']))

def main():
    e = M.EFT(); out = dict(status='CONDITIONAL (reduced model not validated; see B4_VALIDATION.json)')
    Ys = [0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
    out['Y_scan'] = {}
    for Y in Ys:
        out['Y_scan']['%g' % Y] = {v: one(e, Y, 1.0, v) for v in ('registered', 'radonly')}
        a = out['Y_scan']['%g' % Y]
        print('Y=%g tau_m=%.2f R_m/H0^2=%.4f rad_m=%.4f r_reg=%.3f cls_reg=%s | radonly: Om_vac_pl=%s efolds_RD=%.2f' % (
            Y, a['registered']['H0tau_match'], a['registered']['R_over_H0sq_match'], a['registered']['rad_match'], a['registered']['r_match'],
            a['registered']['classification'], a['radonly']['plateau'] and round(a['radonly']['plateau']['Omega_vac'], 4), a['radonly']['efolds_RD']), flush=True)
    # tuning window
    qs = np.round(np.arange(0.97, 1.0301, 0.0025), 6)
    out['window'] = {}
    for Y in (0.3, 1.0, 3.0):
        for v in ('radonly', 'registered'):
            rows = [dict(q=float(q), **{k: x for k, x in one(e, Y, q, v).items() if k in ('rad_era', 'efolds_RD', 'classification', 'H_zero', 'recollapse')}) for q in qs]
            era = [r['q'] for r in rows if r.get('rad_era')]
            rd1 = [r['q'] for r in rows if r.get('efolds_RD', 0) >= 1.0]
            out['window']['Y=%g %s' % (Y, v)] = dict(rows=rows, q_range_rad_era=[min(era), max(era)] if era else None,
                                                    q_range_1efold_RD=[min(rd1), max(rd1)] if rd1 else None,
                                                    width_rad_era=(max(era) - min(era) + 0.0025) if era else 0.0,
                                                    width_1efold_RD=(max(rd1) - min(rd1) + 0.0025) if rd1 else 0.0, grid_step=0.0025)
            w = out['window']['Y=%g %s' % (Y, v)]
            print('window Y=%g %s rad_era %s 1-efold-RD %s' % (Y, v, w['q_range_rad_era'], w['q_range_1efold_RD']), flush=True)
    out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    def clean(v):
        if isinstance(v, dict): return {str(k): clean(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)): return [clean(x) for x in v]
        if isinstance(v, (np.floating, float)): return None if not math.isfinite(float(v)) else float(v)
        if isinstance(v, (np.integer,)): return int(v)
        if isinstance(v, np.bool_): return bool(v)
        return v
    (HERE/'B4_PREDICTIONS.json').write_text(json.dumps(clean(out), indent=1) + '\n')

if __name__ == '__main__':
    main()
