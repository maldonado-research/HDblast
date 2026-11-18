"""Inspect the remote-control histories and optionally recompute the band mode spectra."""
import argparse
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import numpy as np

CASES = [('coarse','positive','0.001','plus'),('fine','positive','0.001','plus'),
         ('coarse','half','0.0005','half'),('fine','half','0.0005','half'),
         ('fine','shift','shifted','shift')]


def read(path):
    return json.loads(path.read_text())


def rows(folder):
    with np.load(folder/'rolling_timeseries.npz',allow_pickle=False) as data:
        return dict(zip(data['cols'],data['rec'].T,strict=True))


def require(ok, message):
    if not ok:
        raise ValueError(message)


def close(x, ref):
    require(math.isclose(x,ref,rel_tol=7e-8,abs_tol=3e-10),f'Replay mismatch {x} and {ref}')


def compute(data, classical, replay):
    out={}
    for grid,name,branch,mode in CASES:
        path=data/'runs'/grid/branch
        initial=read(path/'initial.json')
        native=read(path/'rolling_summary.json')
        event=native['matter']['events'][0]
        qpath=classical/'runs'/grid/'quiet'
        eventq=read(qpath/'rolling_summary.json')['matter']['events'][0]
        quiet=read(classical/'modes'/f'base-{grid}-quiet.json')
        meas=read(data/'modes'/f'{mode}-{grid}.json')
        history=rows(path)
        if replay:
            with tempfile.TemporaryDirectory() as scratch:
                tmp=Path(scratch)
                rg=tmp/grid
                rg.mkdir()
                (rg/branch).symlink_to(path.resolve(),target_is_directory=True)
                (rg/'quiet').symlink_to(qpath.resolve(),target_is_directory=True)
                (rg/'checkpoint.npz').symlink_to((classical/'runs'/grid/'checkpoint.npz').resolve())
                command=[sys.executable,'-I','-B']
                if sys.flags.optimize:
                    command.append('-O')
                subprocess.run(command+[str(classical/'rolling_modes.py'),str(tmp),str(tmp/'new.json'),
                               '--grid',grid,'--branch',branch],check=True,
                               env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
                new=read(tmp/'new.json')
                for cut in meas['ends']:
                    for k in ('comoving_N_8','band_E_upper_H0'):
                        close(new['ends'][cut][k],meas['ends'][cut][k])
                close(new['event_comoving_N_8'],meas['event_comoving_N_8'])
        contact=initial['T']-(initial['control_center']+initial['control_halfwidth'])
        cuts={}
        for cut,m in meas['ends'].items():
            quiet_n=quiet['ends'][cut]['comoving_N_8']
            wall_t=float(np.interp(float(cut),history['H0tau'],history['T']))
            require(native['stop_reason']=='T_final' and history['T'][0]<wall_t<contact and
                    history['H0tau'][0]<float(cut)<history['H0tau'][-1] and
                    m['band_E_upper_H0']<meas['M5_over_H0'],f'Invalid readout {grid}/{branch}/{cut}')
            cuts[cut]={'wall_T':wall_t, 'time_margin_to_compensator':contact-wall_t,
                       'comoving_N_8':m['comoving_N_8'], 'quiet_N_8':quiet_n,
                       'change_percent':100*(m['comoving_N_8']/quiet_n-1),
                       'frequency_upper_H0':m['band_E_upper_H0']}
        require(abs(initial['matched_delta_y_left'])<1e-10 and initial['min_abs_D']>.4,
                'Constraint shooting failure')
        out[f'{grid}_{name}']=dict(remote_center=initial['control_center'],
                               remote_width=initial['control_halfwidth'],remote_coefficient=initial['compensator_coefficient'],
                               earliest_T_compensator=contact,core_incident_front_T=initial['T']+.10,
                               crossing_T=event['T_star'],quiet_crossing_T=eventq['T_star'],
                               event_N8_change_percent=100*(meas['event_comoving_N_8']/quiet['event_comoving_N_8']-1),
                               solved_left_residual=initial['matched_delta_y_left'],cuts=cuts)
    a=rows(data/'runs'/'fine'/'0.001');b=rows(data/'runs'/'fine'/'shifted')
    shift_max={k:float(np.max(abs(a[k]-b[k]))) for k in ('T','phi_b','H_over_H0','A_b','H0tau','B_b','n_hat','rho_hat')}
    require(shift_max['phi_b']<1e-10 and shift_max['H_over_H0']<1e-8,'Remote relocation failed')
    failed={}
    for kind in ('incident','velocity','metric','metric_field'):
        v=read(data/f'profile-{kind}.json')[0]
        require(v['minimum']>.17 and v['admissible_samples']>50,'Changed negative profile')
        failed[kind]=v
    return dict(cases=out,shifted_wall_max_abs_difference=shift_max,
                negative_compact_endpoint_profiles=failed)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--classical',type=Path)
    parser.add_argument('--replay-modes',action='store_true')
    args=parser.parse_args()
    here=Path(__file__).resolve().parent
    print(json.dumps(compute(here,args.classical or here.parent/'classical',args.replay_modes),
                     indent=2,allow_nan=False))
