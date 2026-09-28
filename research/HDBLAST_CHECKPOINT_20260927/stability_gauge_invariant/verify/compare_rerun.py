"""Auditor: compare re-run outputs (verify/rerun/*.json) with the audited JSON files, and the independent
Pruefer-angle solver (indep_spectrum.json) with the audited spectrum.  Output: compare_rerun.json"""
import json, math
from pathlib import Path
H = Path(__file__).resolve().parent; A = H.parent; R = H/'rerun'
L = lambda p: json.loads(p.read_text())
out = {}
for f in ('BACKGROUNDS.json', 'DERIVATION_RESULTS.json', 'SPECTRUM_RESULTS.json', 'CONTROLS_RESULTS.json'):
    a, b = L(A/f), L(R/f)
    keys = [k for k in a if k not in ('runtime_s',)]
    diff = [k for k in keys if a[k] != b.get(k)]
    if f == 'DERIVATION_RESULTS.json':
        diff = [k for k in diff if k != 'runtime_s']
    out[f] = dict(identical_keys=len(keys) - len(diff), differing_keys=diff)
bg = L(A/'BACKGROUNDS.json')
out['max_rel_diff_rho_b_plus_vs_package'] = max(r['rel_diff_rho_b'] for r in bg['plus'])
out['README_states_max_rel_diff'] = 6e-15
sp_ = L(A/'SPECTRUM_RESULTS.json'); grid = sp_['grid']
out['grid_point_used_for_Mhat_0'] = min(grid, key=abs)
ind = L(H/'indep_spectrum.json')
rows = []
for r in ind['branches']:
    if r['branch'] != 'plus':
        continue
    a = [s for s in sp_['plus_summary'] if s['delta'] == r['delta']][0]
    rows.append(dict(delta=r['delta'], indep_min_Mhat=r['scalar']['min_Mhat'], audited_min_Mhat=a['min_Mhat'],
                     abs_diff=abs(r['scalar']['min_Mhat'] - a['min_Mhat']), indep_roots=r['scalar']['roots'],
                     audited_roots=a['scalar_roots'], B_indep=r['shell']['B'], B_audited=a['B'],
                     Bstar_range_mu2_neg=[r['Bstar_control']['Bstar_min_mu2_neg'], r['Bstar_control']['Bstar_max_mu2_neg']],
                     stability_margin_B_minus_maxBstar=r['Bstar_control']['margin'], tensor_roots_indep=r['tensor']['roots']))
out['plus_indep_vs_audited'] = rows
o = {r['delta']: r['scalar']['roots'][0] for r in ind['branches'] if r['branch'] == 'original'}
c = 0.5975949350280132
cf = -4*(3*c*c - 4*c + 8)/(c*(3*c + 4)); sl = (o[0.003] - o[0.001])/0.002
out['calibration_indep'] = dict(roots=o, audited_finest=sp_['calibration']['mu2'], diff_indep_vs_audited=o[0.001] - sp_['calibration']['mu2'],
                                diff_indep_vs_chat9_GN=o[0.001] + 7.717871625176294,
                                C4_intercept_indep=o[0.001] - 0.001*sl, closed_form_4D=cf, C4_diff=o[0.001] - 0.001*sl - cf)
(H/'compare_rerun.json').write_text(json.dumps(out, indent=2) + '\n')
print(json.dumps(out, indent=1))
