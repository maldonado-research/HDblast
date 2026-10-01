#!/usr/bin/env python3
"""B4 step 2: registered validation of reduced model E against the A1 5D numbers (REGISTRATION.md section 2),
controls C1, C3, systematics (phi_m, rtol), post-registration variants, and diagnostics against the A1 time series
(read-only).  Output: B4_VALIDATION.json.  Numerical."""
import json, math, hashlib, sys
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import eft_reduced as M
A1 = Path('/home/user/unified-theory-maldonado/research/HDBLAST_CHECKPOINT_20260928/A1_tuned_vacuum')

TARGETS = dict(T1=('r', 0.3, 1.0, 0.463, 'rel', 0.30), T2=('r', 1.0, 1.0, 0.0734, 'rel', 0.30),
               T3=('Omega_vac', 0.3, 1.0, -0.148, 'abs', 0.05), T4=('Omega_vac', 1.0, 1.0, -0.207, 'abs', 0.05),
               T5=('H_zero', 1.0, 1.01, 13.54, 'rel', 0.15), T6=('recollapse', 1.0, 1.01, 18.85, 'rel', 0.15),
               T7=('H_zero', 1.0, 1.05, 6.93, 'rel', 0.15), T8=('recollapse', 1.0, 1.05, 8.20, 'rel', 0.15))

def get(an, what):
    if what == 'r': return an['plateau']['r'] if an['plateau'] else None
    if what == 'Omega_vac': return an['plateau']['Omega_vac'] if an['plateau'] else None
    if what == 'H_zero': return an['H_zero_H0tau']
    if what == 'recollapse': return an['recollapse_H0tau']

def evaluate(e, variant='registered', phi_m=0.95, rtol=1e-10, dcs=(1e-2, 1e-4), wrong=None):
    cache = {}; rows = {}; allpass = True
    for k, (what, Y, q, tgt, kind, tol) in TARGETS.items():
        vals = []
        for dc in (dcs if q == 1.0 else dcs[:1]):     # A1 detuned controls (1.01 d*, 1.05 d*) exist only for dc = 1e-2
            key = (Y, q, dc)
            if key not in cache:
                res = M.run(e, Y, q*M.DSTAR, dc, phi_m, rtol=rtol, variant=variant, wrong=wrong)
                cache[key] = (M.analyse(res), res['match']) if res['ok'] else (None, res)
            an = cache[key][0]
            vals.append(get(an, what) if an else None)
        ok = []
        for v in vals:
            if v is None: ok.append(False); continue
            dev = abs(v - tgt)/abs(tgt) if kind == 'rel' else abs(v - tgt)
            side = True
            if k == 'T1': side = v > 0.1
            if k == 'T2': side = v <= 0.1
            ok.append(bool(dev <= tol and side))
        rows[k] = dict(quantity=what, Y=Y, d_over_dstar=q, target=tgt, tolerance=('%s %.2f' % (kind, tol)), values_dc=vals, PASS=all(ok))
        allpass &= all(ok)
    late = all(rows[k]['PASS'] for k in ('T3', 'T4', 'T5', 'T6', 'T7', 'T8'))
    verdict = 'VALIDATED' if allpass else ('PARTIAL (late sector validated)' if late else 'NOT VALIDATED')
    matches = {('Y=%g q=%g dc=%g' % k): v[1] for k, v in cache.items()}
    classes = {('Y=%g q=%g dc=%g' % k): (v[0]['classification'] if v[0] else None) for k, v in cache.items()}
    return dict(variant=variant, phi_m=phi_m, rtol=rtol, wrong=wrong, targets=rows, verdict=verdict, all_pass=allpass,
                matches=matches, classifications=classes)

def load5d(tag):
    sub = 'ctl' if tag.startswith('ctl') else 'main'
    d = np.load(A1/'runs'/sub/(tag + '_timeseries.npz')); cols = list(d['cols']); R = d['rec']
    return {c: R[:, i] for i, c in enumerate(cols)}

def diag_recipe_on_5d(phi_m=0.95):
    """Post-registration diagnostic: apply the registered matching recipe to the A1 5D trajectory itself."""
    out = {}
    lam = float(np.polyval(M.FIT['coef_Lambda'], M.DSTAR)); sig_s = float(np.polyval(M.FIT['coef_sigma'], M.DSTAR))
    rb = 7.838285538074403
    for tag, rpl in [('main_dstar_Y1_dc1e-2_dzf5e-4', 0.07337), ('main_dstar_Y0.3_dc1e-2_dzf5e-4', 0.4629)]:
        D = load5d(tag); i = int(np.flatnonzero(D['phi_b'] >= phi_m)[0])
        h = D['H_over_H0'][i]; R = D['R'][i]; rad = rb**2*(sig_s*R/18 + R*R/36)
        w = h*h - lam - rad
        out[tag] = dict(H0tau=float(D['H0tau'][i]), phi_b=float(D['phi_b'][i]), recipe_r=float(w/rad), plateau_r_5D=rpl,
                        Wa4_ratio_plateau_over_match=None)
    return out

def diag_trajectory(e, phi_pts=(0.1, 0.3, 0.5, 0.7, 0.9, 0.95)):
    """Post-registration diagnostic: EFT phase-I trajectory against the A1 5D trajectory at equal phi_b (dc = 1e-2)."""
    out = {}
    for Y, tag in [(1.0, 'main_dstar_Y1_dc1e-2_dzf5e-4'), (0.3, 'main_dstar_Y0.3_dc1e-2_dzf5e-4'), (0.0, 'main_dstar_Y0_dc1e-2_dzf5e-4')]:
        res = M.run(e, Y, M.DSTAR, 1e-2, 0.95, variant='registered')
        D = res['D']; F = load5d(tag)
        rows = []
        for p in phi_pts:
            iE = int(np.flatnonzero(D['phi_b'] >= p)[0]); i5 = int(np.flatnonzero(F['phi_b'] >= p)[0])
            rbsq = 7.838285538074403**2
            rows.append(dict(phi=p, H0tau_EFT=float(D['H0tau'][iE]), H0tau_5D=float(F['H0tau'][i5]),
                             H_EFT=float(D['H_over_H0'][iE]), H_5D=float(F['H_over_H0'][i5]),
                             lna_EFT=float(D['ln_a'][iE]), lna_5D=float(F['ln_a'][i5]),
                             R_over_H0sq_EFT=float(D['R'][iE]/res['match']['H0_EFT']**2), R_over_H0sq_5D=float(F['R'][i5]*rbsq),
                             v_EFT=float(D['v_over_H0'][iE]), v_5D=float(F['v_over_H0'][i5])))
        out['Y=%g' % Y] = rows
    return out

def main():
    e = M.EFT()
    out = dict(status='numerical; reduced model (conditional unless validated)', registration='REGISTRATION.md')
    g = M.growth_rate(e)
    out['C1_growth'] = dict(g, target_5D_delta_1em3=1.65719, closed_form_mu2=-7.719795918,
                            PASS=bool(abs(g['rate_over_H'] - 1.65719) <= 2e-3 and abs(g['mu2'] + 7.719795918) <= 2e-3))
    gd = M.growth_rate(e, d=M.DSTAR); out['growth_at_dstar'] = gd
    print('C1', out['C1_growth'], 'd*', gd, flush=True)
    out['primary'] = evaluate(e)
    print('primary verdict', out['primary']['verdict'], {k: (v['values_dc'], v['PASS']) for k, v in out['primary']['targets'].items()}, flush=True)
    out['rtol_check'] = evaluate(e, rtol=1e-12, dcs=(1e-2,))
    a = out['primary']['targets']; b = out['rtol_check']['targets']
    out['rtol_max_rel_diff'] = max(abs(a[k]['values_dc'][0] - b[k]['values_dc'][0])/abs(b[k]['values_dc'][0])
                                   for k in a if a[k]['values_dc'][0] is not None and b[k]['values_dc'][0] is not None)
    out['systematics_phi_m'] = {('%.2f' % pm): evaluate(e, phi_m=pm, dcs=(1e-2,)) for pm in (0.90, 0.98)}
    out['C3_wrong_gamma'] = evaluate(e, wrong='gamma', dcs=(1e-2,))
    out['C3_wrong_nofdot'] = evaluate(e, wrong='nofdot', dcs=(1e-2,))
    p0 = out['primary']['targets']['T2']['values_dc'][0]
    out['C3_sensitivity'] = {k: dict(T2=out[k]['targets']['T2']['values_dc'][0],
                                     changes_T2_by_more_than_tol=bool(abs(out[k]['targets']['T2']['values_dc'][0] - p0) > 0.3*0.0734))
                             for k in ('C3_wrong_gamma', 'C3_wrong_nofdot')}
    out['post_registration_variants'] = {v: evaluate(e, variant=v) for v in ('einstein', 'radonly')}
    out['post_registration_diag_recipe_on_5D'] = diag_recipe_on_5d()
    out['post_registration_diag_trajectory'] = diag_trajectory(e)
    out['script_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out['model_sha256'] = hashlib.sha256((HERE/'eft_reduced.py').read_bytes()).hexdigest()
    def clean(v):
        if isinstance(v, dict): return {str(k): clean(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)): return [clean(x) for x in v]
        if isinstance(v, (np.floating, float)): return None if not math.isfinite(float(v)) else float(v)
        if isinstance(v, (np.integer,)): return int(v)
        if isinstance(v, np.bool_): return bool(v)
        return v
    (HERE/'B4_VALIDATION.json').write_text(json.dumps(clean(out), indent=1) + '\n')
    for k in ('einstein', 'radonly'):
        print(k, out['post_registration_variants'][k]['verdict'], {t: (v['values_dc'], v['PASS']) for t, v in out['post_registration_variants'][k]['targets'].items()})
    for pm, v in out['systematics_phi_m'].items(): print('phi_m', pm, v['verdict'], {t: x['values_dc'] for t, x in v['targets'].items()})
    print('rtol', out['rtol_max_rel_diff']); print('C3', out['C3_sensitivity'])
    print(json.dumps(clean(out['post_registration_diag_recipe_on_5D']), indent=1))
    print(json.dumps(clean(out['post_registration_diag_trajectory']), indent=0)[:4000])

if __name__ == '__main__':
    main()
