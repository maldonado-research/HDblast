"""Audit V2: independent consistency checks on the workstream's coupled-5D pilot time series (read-only).

For each pilot run in ../pilot_5d/runs (and the audit's own re-runs in ./pilot):
 (1) Weyl transport identity WITH shell radiation (re-derived in V1, check D2):
        dWy/dtau + 4 H Wy = (1/6)[4 k w v - 4 H v^2 - v vdot + w wdot - U' v],  k = (sigma+R)/6,  w = -(sigma' + Y v)/2
     evaluated on the saved samples (hatted units, H0 = 1/rho_b).  This is NOT enforced by the solver: it holds only
     if the bulk evolution, the three junctions and the radiation ODE are mutually consistent (4D Bianchi identity).
     Controls: (a) k = sigma/6 (radiation dropped from nA), (b) w = -sigma'/2 (friction dropped from the scalar junction),
     (c) 4H -> 3H.
 (2) H^2 budget decomposition at the end of the reliable window (same window rule as m6: first t>1 with d(H0 tau)/dt<0.2),
     in H0^2 units:  vacuum part sigma^2/36 - sigma'^2/48 + U/6,  radiation terms sigma R/18 + R^2/36,  scalar kinetic v^2/12,
     friction cross terms -(2 sigma' Y v + Y^2 v^2)/48,  Weyl Wy.  Reports the fraction of H^2 in each.
 (3) Weyl a^4 and R a^4 over the last 0.2 H0 tau of the window (the workstream's 'still falling' statement).
Output: V2_PILOT_CONSISTENCY.json
"""
import json, math, glob
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
RUNS = sorted(glob.glob(str(HERE.parent/'pilot_5d/runs/*_summary.json'))) + sorted(glob.glob(str(HERE/'pilot/*_summary.json')))

def W(p): return 1 - p + p**3/3
def Wp(p): return p*p - 1
def U(p): return 0.5*(p*p-1)**2 - (2/3)*W(p)**2
def Up(p): return (p*p-1)*2*p - (4/3)*W(p)*(p*p-1)

def analyse(fsum):
    s = json.load(open(fsum)); rec = np.load(fsum.replace('_summary.json', '_timeseries.npz'))['rec']
    td, c, dq, Y, rb = s['t_det'], s.get('c', 2/1.0357712571566784 - 4/3), s.get('d_quad', 0.0), s['Y'], s['rho_b']
    sig = lambda p: 2*W(p) + td*(1 + c*p + dq*p*p/2)
    sig1 = lambda p: 2*Wp(p) + td*(c + dq*p)
    sig2 = lambda p: 4*p + td*dq
    t = rec[:, 0]; ss = rec[:, 8]
    lapse = np.gradient(ss, t); bad = np.where((lapse < 0.2) & (t > 1.0))[0]
    ie = int(bad[0]) - 1 if len(bad) else len(t) - 1
    keep = np.concatenate(([True], np.diff(ss) > 1e-12)); keep[ie+1:] = False
    phi, h, vh, R, Wy, s_ = rec[keep, 1], rec[keep, 2], rec[keep, 5], rec[keep, 6], rec[keep, 7], ss[keep]
    v = vh/rb
    # hatted quantities (units of H0 = 1/rb)
    k = rb*(sig(phi) + R)/6
    wv = -(rb*sig1(phi) + Y*vh)/2
    wdot = -(rb*sig2(phi)*vh + Y*np.gradient(vh, s_, edge_order=2))/2
    vdot = np.gradient(vh, s_, edge_order=2)
    dWy = np.gradient(Wy, s_, edge_order=2)
    def resid(kk, ww, wwd, fac=4):
        rhs = (4*kk*ww*vh - 4*h*vh*vh - vh*vdot + ww*wwd - rb**2*Up(phi)*vh)/6
        lhs = dWy + fac*h*Wy
        scale = (np.abs(4*kk*ww*vh) + np.abs(4*h*vh*vh) + np.abs(vh*vdot) + np.abs(ww*wwd) + np.abs(rb**2*Up(phi)*vh))/6 + np.abs(dWy) + np.abs(fac*h*Wy) + 1e-300
        r = (np.abs(lhs - rhs)/scale)[3:-3]
        return dict(median=float(np.median(r)), p95=float(np.percentile(r, 95)), max=float(np.max(r)))
    out = dict(tag=Path(fsum).name.replace('_summary.json', ''), Y=Y, c=c, d_quad=dq, dc=s['dc'], dz_fine=s['dz_fine'],
               window_end_t=float(t[ie]), window_end_H0tau=float(ss[ie]), n_used=int(keep.sum()))
    out['weyl_transport_with_matter'] = resid(k, wv, wdot)
    out['control_nA_without_R'] = resid(rb*sig(phi)/6, wv, wdot)
    wv0 = -(rb*sig1(phi))/2; wd0 = -(rb*sig2(phi)*vh)/2
    out['control_scalar_junction_without_friction'] = resid(k, wv0, wd0) if Y > 0 else None
    out['control_4H_to_3H'] = resid(k, wv, wdot, fac=3)
    # budget at window end (index -1 of kept arrays)
    j = -1; p = phi[j]; Rv = R[j]; vv = v[j]
    terms = dict(vacuum_sigma2_minus_sigma1sq_plus_U=rb*rb*(sig(p)**2/36 - sig1(p)**2/48 + U(p)/6),
                 radiation_terms=rb*rb*(sig(p)*Rv/18 + Rv*Rv/36),
                 scalar_kinetic=rb*rb*vv*vv/12,
                 friction_cross_terms=-rb*rb*(2*sig1(p)*Y*vv + Y*Y*vv*vv)/48,
                 weyl=float(Wy[j]))
    H2 = float(h[j]**2); tot = sum(terms.values())
    out['budget_at_window_end_H0sq_units'] = dict(H2=H2, sum_of_terms=float(tot), closure_rel=float(abs(tot - H2)/H2),
                                                   fractions_of_H2={kk: float(vv_/H2) for kk, vv_ in terms.items()},
                                                   phi_b=float(p), H_over_H0=float(h[j]), R=float(Rv))
    # a^4-scaled quantities over the last 0.2 H0 tau
    lna = np.concatenate(([0.], np.cumsum(0.5*(h[1:] + h[:-1])*np.diff(s_))))
    i0 = int(np.searchsorted(s_, s_[-1] - 0.2))
    wa4 = Wy*np.exp(4*lna); ra4 = R*np.exp(4*lna)
    out['last_0.2_H0tau'] = dict(Wy_a4_rel_change=float((wa4[-1] - wa4[i0])/wa4[i0]) if wa4[i0] != 0 else None,
                                 R_a4_rel_change=float((ra4[-1] - ra4[i0])/ra4[i0]) if ra4[i0] > 0 else None,
                                 ln_a_advance=float(lna[-1] - lna[i0]),
                                 vacuum_fraction_start_end=None)
    return out

def main():
    res = [analyse(f) for f in RUNS]
    for r in res:
        b = r['budget_at_window_end_H0sq_units']
        print('%-45s Y=%-4g d=%-7.3f win=%.3f  weyl-id med=%.1e p95=%.1e | noR=%.1e nofric=%s 3H=%.1e | fr: vac=%.3f rad=%.3f kin=%.3f fric=%.3f Wy=%.3f clos=%.1e' % (
            r['tag'], r['Y'], r['d_quad'], r['window_end_H0tau'], r['weyl_transport_with_matter']['median'], r['weyl_transport_with_matter']['p95'],
            r['control_nA_without_R']['median'], ('%.1e' % r['control_scalar_junction_without_friction']['median']) if r['control_scalar_junction_without_friction'] else '-',
            r['control_4H_to_3H']['median'], *[b['fractions_of_H2'][k] for k in ['vacuum_sigma2_minus_sigma1sq_plus_U', 'radiation_terms', 'scalar_kinetic', 'friction_cross_terms', 'weyl']], b['closure_rel']))
    (HERE/'V2_PILOT_CONSISTENCY.json').write_text(json.dumps(dict(
        status='numerical (derived from saved floating-point pilot time series; centred finite differences on ~0.01-t samples)',
        runs=res), indent=1) + '\n')

if __name__ == '__main__':
    main()
