"""Audit V4: independent arithmetic re-computation of headline numbers from first principles / saved inputs.
Output: V4_ARITHMETIC_CHECKS.json"""
import json, math
from pathlib import Path
HERE = Path(__file__).resolve().parent; M = HERE.parent
m2 = json.load(open(M/'M2_ENERGY_BUDGET_AND_SCALES.json')); m1 = json.load(open(M/'M1_WEYL_ALONG_CHAT14.json'))
c = 2/1.0357712571566784 - 4/3
out = {}
# f_H leading-order limit: H_vac^2/delta -> (1+c)/27, H0^2/delta -> 0.1609 (M2 smallest delta)
row = m2['rows'][0]
out['fH_leading_limit'] = ((1+c)/27)/row['budget']['H02_over_delta']
out['fH_by_delta'] = {str(r['delta']): r['budget']['H_vac2_over_H02'] for r in m2['rows']}
# f_E and N_max from f_H and M_i^2/M_f^2 (V_E ∝ H^2/M^2)
reg = [r for r in m2['rows'] if r['delta'] == 1e-3][0]['budget']
fE = reg['H_vac2_over_H02']*reg['M_i2_over_M_f2']; out['fE_recomputed'] = fE
out['N_max_recomputed'] = 0.25*math.log((1-fE)/fE)
out['wrong_direction_control_fE_with_inverted_mass_ratio'] = reg['H_vac2_over_H02']/reg['M_i2_over_M_f2']
# target: rho_r/rho_Lambda at T = 5 MeV, g* = 10.75, h = 0.674, Omega_L = 0.685, reduced Planck mass 2.435e27 eV
H0 = 2.1332e-33*0.674; rho_L = 3*(2.435e27)**2*0.685*H0**2; rho_r = math.pi**2/30*10.75*(5e6)**4
out['target_ratio'] = rho_r/rho_L; out['target_N'] = 0.25*math.log(rho_r/rho_L)
# Delta N_eff -> rho_DR/rho_SM
for dn in (0.3, 0.09):
    r_bbn = dn*(7/8)*2/10.75; out['dNeff_%g' % dn] = dict(bbn=r_bbn, production=r_bbn*(106.75/10.75)**(1/3))
# delta-Lambda lock
HL = math.sqrt(0.685)*H0; hbarc = 1.973269804e-7
for ell in (38.6, 0.1):
    invL = hbarc/(ell*1e-6/9)
    out['delta_required_ell_%g_um' % ell] = (HL/invL)**2/reg['Hvac2_over_delta']
gamma = 1.657193663767*HL/math.sqrt(reg['H_vac2_over_H02'])
out['efold_time_Gyr'] = 6.582119569e-16/gamma/3.15576e16
out['bounce_H_over_Hvac'] = math.sqrt(out['target_ratio'])
# Weyl peak share of H^2 in the Chat14 +1 runs
pk = {}
for k, v in m1['runs'].items():
    a = v['analysis']
    if 'plus' in k:
        import numpy as np
        s = np.array(a['s_sample']); W = np.array(a['Wy_sample'])
        pk[k] = dict(peak=a['Wy_over_H0sq_max'], at=a['s_at_Wy_max'], p95_identity=a['weyl_identity_p95_rel_residual'], median_identity=a['weyl_identity_median_rel_residual'], end=a['Wy_over_H0sq_end'])
out['chat14_plus_runs'] = pk
(HERE/'V4_ARITHMETIC_CHECKS.json').write_text(json.dumps(out, indent=1) + '\n')
print(json.dumps(out, indent=1))
