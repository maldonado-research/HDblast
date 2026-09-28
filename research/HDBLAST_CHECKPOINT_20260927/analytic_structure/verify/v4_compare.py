"""Compare the auditor's independent BVP solutions (runs/V2_BVP_*.json) with (a) the audited BVP values,
(b) the audited exact series through delta^8 and (c) the auditor's exact local delta^9 coefficients
(V1_SERIES_INDEPENDENT.json).  Also re-checks selected numbers quoted in the audited README against its JSON.
Writes verify/V4_COMPARE.json
"""
import json
from pathlib import Path
import mpmath as mp
import sympy as sp

mp.mp.dps = 60
HERE = Path(__file__).resolve().parent; AUD = HERE.parent
ser = json.loads((AUD / 'SERIES_COEFFICIENTS.json').read_text())
v1 = json.loads((HERE / 'V1_SERIES_INDEPENDENT.json').read_text())
cs = sp.Symbol('c')
def coeffs(exprs, cval):
    return [mp.mpf(str(sp.N(sp.sympify(e, locals={'c': cs}).subs(cs, sp.Float(str(cval), 70)), 60))) for e in exprs]
out = {}
aud_rows = {}
for f in (AUD / 'runs').glob('BVP_SCAN_*.json'):
    s = json.loads(f.read_text())
    for r in s['rows']:
        aud_rows[(s['c'][:20], r['delta'])] = r
for f in sorted((HERE / 'runs').glob('V2_BVP_*.json')):
    s = json.loads(f.read_text()); cval = mp.mpf(s['c'])
    q8 = {n: coeffs(ser[n], cval) for n in ('eta_b', 'H2')}
    q9 = {n: coeffs(v1[n], cval) for n in ('eta_b', 'H2')}   # auditor's series incl. local delta^9 term
    rows = []
    for r in s['rows']:
        d = mp.mpf(r['delta'])
        rec = dict(delta=r['delta'], J2=r['J2'])
        for n in ('eta_b', 'H2'):
            Q = mp.mpf(r[n])
            S8 = mp.fsum(q8[n][i] * d**i for i in range(9)); S9 = S8 + q9[n][9] * d**9
            rec[n] = dict(value=r[n], minus_series8_over_d9=mp.nstr((Q - S8) / d**9, 12), exact_local_q9=mp.nstr(q9[n][9], 12),
                          abs_diff_series8=mp.nstr(abs(Q - S8), 5), abs_diff_series9=mp.nstr(abs(Q - S9), 5),
                          ratio_diff9_to_d10=mp.nstr(abs(Q - S9) / d**10, 5))
        key = (mp.nstr(cval, 20)[:20], r['delta'])
        # match to audited BVP rows with same c and delta
        match = [v for k, v in aud_rows.items() if k[1] == r['delta'] and abs(mp.mpf(k[0]) - cval) < mp.mpf('1e-15')]
        if match:
            a = match[0]
            rec['vs_audited_BVP_rel'] = {n: (mp.nstr(abs(mp.mpf(r[n]) - mp.mpf(a[n])) / abs(mp.mpf(a[n])), 4) if mp.mpf(a[n]) != 0 else 'abs ' + mp.nstr(abs(mp.mpf(r[n]) - mp.mpf(a[n])), 4)) for n in ('eta_b', 'H2', 'rho_b', 'eta_h')}
        rows.append(rec)
    out[f.stem] = dict(c=s['c'], dps=s['dps'], u0=s['u0'], rows=rows)
# cross-setting convergence of the auditor's own solver (reg40 vs reg50)
if 'V2_BVP_reg40' in out and 'V2_BVP_reg50' in out:
    A = {r['delta']: r for r in json.loads((HERE / 'runs' / 'V2_BVP_reg40.json').read_text())['rows']}
    B = {r['delta']: r for r in json.loads((HERE / 'runs' / 'V2_BVP_reg50.json').read_text())['rows']}
    out['auditor_solver_convergence_dps40_u0_1e-3_vs_dps50_u0_5e-4'] = [
        dict(delta=k, **{n: mp.nstr(abs(mp.mpf(A[k][n]) - mp.mpf(B[k][n])) / abs(mp.mpf(B[k][n])), 4) for n in ('eta_b', 'H2', 'eta_h')})
        for k in sorted(set(A) & set(B))]
# eta_h of the auditor's BVP vs audited eta_h series (through delta^7 relative), reg only
for lab in ('V2_BVP_reg40',):
    if lab in out:
        s = json.loads((HERE / 'runs' / f'{lab}.json').read_text()); cval = mp.mpf(s['c'])
        qe = coeffs(ser['eta_h_over_delta8'], cval)
        out[lab]['eta_h_check'] = [dict(delta=r['delta'], rel_dev_from_series=mp.nstr(abs(mp.mpf(r['eta_h']) / mp.mpf(r['delta'])**8 - mp.fsum(qe[i] * mp.mpf(r['delta'])**i for i in range(8))) / abs(qe[0]), 5))
                                   for r in s['rows']]
# audited claim: scaled remainders deviate proportionally to delta (ratio ~2 per doubling) -- read from audited JSON
sv = json.loads((AUD / 'SERIES_VS_BVP.json').read_text())
ratios = {}
for scan in ('BVP_SCAN_reg', 'BVP_SCAN_cp13', 'BVP_SCAN_cm04'):
    per = {r['delta']: r for r in sv['scans'][scan]['per_delta']}
    rr = {}
    for nme in ('eta_b', 'H2', 'rho_b_factor', 'y_b_correction', 'eta_h_over_delta8'):
        lst = []
        for n in range(len(per['0.0008'][nme]['remainders'])):
            a = mp.mpf(per['0.0008'][nme]['remainders'][n]['rel_dev']); b_ = mp.mpf(per['0.0016'][nme]['remainders'][n]['rel_dev'])
            lst.append(mp.nstr(b_ / a, 4) if a != 0 else 'nan')
        rr[nme] = lst
    ratios[scan] = rr
out['audited_scaled_remainder_ratio_dev(0.0016)/dev(0.0008)_expect_2'] = ratios
(HERE / 'V4_COMPARE.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out, indent=1))
