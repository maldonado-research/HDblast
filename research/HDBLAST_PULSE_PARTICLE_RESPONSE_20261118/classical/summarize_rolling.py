"""Event, timing, proper-time and mesh checks for the matched rolling runs."""
import argparse
import json
from pathlib import Path
import numpy as np

def analyze(root):
    output={}
    for grid in ('coarse','fine'):
        quiet=None; output[grid]={}
        for branch in ('quiet','0.001','-0.001'):
            folder=root/grid/branch
            raw=json.loads((folder/'rolling_summary.json').read_text())
            events=raw['matter']['events']
            assert raw['stop_reason']=='T_final' and len(events)==1
            e=events[0]
            with np.load(folder/'rolling_timeseries.npz') as data:
                col={k:data['rec'][:,i] for i,k in enumerate(data['cols'])}
            if quiet is None: quiet=col
            t=col['T'];s=col['H0tau'];v=col['v_over_H0']
            pre=t<=1.20
            sfront=np.interp([1.21,1.25],t,s)
            sstar=e['s_star']; scale=e['q']**-.5
            in_scale=(s>=sstar-scale)&(s<=sstar+scale)
            read=lambda c:float(np.interp(2.1,s,col[c]))
            source=read('Q_hat')
            summary=dict(T_star=e['T_star'],s_star=sstar,v_star=e['v_star'],q=e['q'],v_fit=e['v_fit'],
                n_target=e['n_hat_at_crossing'],n_leading=e['q']**1.5/(8*np.pi**3),
                derivative_factor=e['number_factor'],D_over_q=e['D_over_q'],
                G_min_screen=e['G_min_screen'],passes_screen=e['screen_pass'],
                window_proper_width=float(np.diff(sfront)[0]),intrinsic_proper_width=scale,
                max_rel_velocity_change_in_plusminus_q_mhalf=float(np.max(abs(v[in_scale]/e['v_star']-1))),
                pre_arrival_max_abs_phi_difference=float(max(abs(col['phi_b'][pre]-quiet['phi_b'][pre]))),
                pre_arrival_max_abs_v_difference=float(max(abs(v[pre]-quiet['v_over_H0'][pre]))),
                R_at_s2_1=read('R'),Qhat_at_s2_1=source,proper_nhat_at_s2_1=read('n_hat'),
                source_event_jump_kappa5sq_rho_over_H0=e['kappa5sq_rho_jump_over_H0'])
            if branch!='quiet':
                summary['construction']=json.loads((folder/'initial.json').read_text())
            output[grid][branch]=summary
    return output

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('runs',type=Path);p.add_argument('output',type=Path)
    a=p.parse_args();a.output.write_text(json.dumps(analyze(a.runs),indent=2)+'\n')
