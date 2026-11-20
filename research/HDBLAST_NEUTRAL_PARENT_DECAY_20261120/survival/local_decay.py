"""Two-body phase space, survival and daughter injection on the neutral wall trajectory."""
import argparse
import importlib.util
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicHermiteSpline
from scipy.optimize import brentq
from scipy.special import expit

HERE=Path(__file__).resolve().parent
REFERENCE='5dbed28e158584403a4cd238d392aab15c00480522cf489d1ee0a36f1d550f90'
PRIMARY_STAR=.77
WIDTH=.03
LANDMARKS=(9.5,10.02608,10.5017,10.82926,11.0)
Q=1e8


def history(path,stride=1):
    raw=path.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=REFERENCE:
        raise ValueError("The pinned charged-gas history changed.")
    data=json.loads(raw)
    row=data['runs']['neutral_precursor']
    grid=np.asarray(row['trace'])
    idx=np.unique(np.r_[np.arange(0,len(grid),stride),len(grid)-1])
    col={key:grid[idx,i] for i,key in enumerate(data['trace_columns'])}
    t=col['t']
    psi=CubicHermiteSpline(t,col['psi'],col['pi'])
    a=CubicHermiteSpline(t,col['a'],col['H']*col['a'])
    return psi,a,row['maximum_tau']


def gate(psi,star=PRIMARY_STAR,width=WIDTH):
    z=(np.asarray(psi)-star)/width+.5
    w=np.ones_like(z)
    inner=(z>0)&(z<1)
    w[z<=0]=0.
    w[inner]=expit(1/(1-z[inner])-1/z[inner])
    return .95+.1*w


def interval(psi,a,turn,star=PRIMARY_STAR,width=WIDTH,dt=.001):
    t=np.unique(np.r_[np.linspace(0,turn,math.ceil(turn/dt)+1),LANDMARKS,turn])
    t=t[(t<=turn)&(t>=0)]
    R=gate(psi(t),star,width)
    factor=R*np.maximum(0.,1.-1./(R*R))**1.5
    cum=cumulative_trapezoid(factor,t,initial=0.)
    start=brentq(lambda x:psi(x)-star,0.,turn,xtol=1e-13)
    fully_open=brentq(lambda x:psi(x)-star-.5*width,start,turn,xtol=1e-13)
    return dict(t=t,a=a(t),R=R,shape=factor,J=cum,start=start,full=fully_open,turn=turn)


def case(win,kappa,full=False):
    t,R,a,J=map(win.__getitem__,('t','R','a','J'))
    S=np.exp(-kappa*J)
    F=-np.expm1(-kappa*J)
    births=S[:-1]*(-np.expm1(-kappa*np.diff(J)))
    mass_mid=.5*(R[:-1]+R[1:])
    energy_in=.5*(S[1:]+S[:-1])@np.diff(R)
    r_initial=float(R[0])
    kin_birth=float(np.dot(mass_mid-1.,births))
    mom2=np.maximum(0,mass_mid**2-1.)*((a[:-1]+a[1:])/(2*a[-1]))**2
    kin_end=float(np.dot(mom2/(1+np.sqrt(1+mom2)),births))
    result=dict(kappa=kappa,final_fraction=float(F[-1]),field_work_in_units_2m=float(energy_in),
                daughter_birth_kinetic_in_units_2m=kin_birth,daughter_turn_kinetic_in_units_2m=kin_end,
                expansion_work_in_units_2m=kin_birth-kin_end,
                comoving_energy_balance_error=abs(float(S[-1]*R[-1]+np.dot(mass_mid,births)-r_initial-energy_in)))
    if result['comoving_energy_balance_error']>5e-12:
        raise ValueError('Parent/daughter mass-work identity failed.')
    if full:
        result['landmarks']={str(l):float(np.interp(l,t,F)) for l in LANDMARKS}
        result['first_arrivals']={}
        for f in (.01,.1,.5,.9,.99):
            result['first_arrivals'][str(f)]=float(np.interp(-math.log1p(-f)/kappa,J,t)) if kappa and F[-1]>=f else None
    return result


def calculate(source):
    psi,a,turn=history(source)
    main=interval(psi,a,turn)
    t,J,R=main['t'],main['J'],main['R']
    base={str(k):case(main,k,True) for k in (0,1,10,50,100,1000)}
    trial={}
    for star in (.42,.64,.74,.77,.79):
        win=interval(psi,a,turn,star)
        trial[str(star)]=dict(threshold_time=win['start'],terminal_mass_time=win['full'],
            available_exposure=float(win['J'][-1]),
            kappa_for_half_before_turn=math.log(2)/win['J'][-1],
            kappa_for_90pct_before_turn=math.log(10)/win['J'][-1])
    width_test={}
    for width in (.015,.03,.06):
        win=interval(psi,a,turn,PRIMARY_STAR,width)
        width_test[str(width)]=dict(threshold_time=win['start'],terminal_mass_time=win['full'],
            available_exposure=float(win['J'][-1]),kappa_for_90pct_before_turn=math.log(10)/win['J'][-1])
    requirements={}
    for l in (*LANDMARKS,turn):
        exposure=float(np.interp(l,t,J))
        requirements[str(l)]=dict(exposure=exposure,kappa_for_90pct=math.log(10)/exposure if exposure>0 else None)
    plateau=1.05*(1-1/1.05**2)**1.5
    # Q is the force parameter of the charge-projected triad. g^2/(4pi) is
    # the neutral Yukawa coupling and does not equal the electromagnetic one.
    matched=[]
    for q in (.1,1.,10.):
        mu=q/(math.sqrt(3)*Q)
        matched.append(dict(charge_norm=q,daughter_compton_mass=mu,
            weak_yukawa_fraction_upper=float(-np.expm1(-mu*J[-1])),
            yukawa_g2_over_4pi_for_half=math.log(2)/(mu*J[-1]),
            yukawa_g2_over_4pi_for_90pct=math.log(10)/(mu*J[-1])))
    fine=interval(psi,a,turn,dt=.0005)
    p2,a2,tr=history(source,2)
    coarse=interval(p2,a2,tr)
    # Accumulated detuning phase of the zero-momentum pair after opening.
    phase_per_mu=float(np.trapezoid(2*np.maximum(0,R-1),t))
    spec=importlib.util.spec_from_file_location('weak_yukawa_window',HERE/'finite_time.py')
    finite=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(finite)
    finite_time=finite.calculate(psi,a,turn)
    return dict(history_sha256=REFERENCE,finite_time=finite_time,parameters=dict(Q=Q,star=PRIMARY_STAR,width=WIDTH,
            initial_mass_over_2m=.95,final_mass_over_2m=1.05,initial_radius=2.),
        threshold_time=main['start'],fully_open_time=main['full'],turn=turn,
        peak_rate_over_kappa=plateau,integrated_exposure=float(J[-1]),
        available_superthreshold_phase_per_daughter_mass=phase_per_mu,
        width_over_parent_mass_per_yukawa_alpha=plateau/(2*1.05),
        primary=base,scalar_threshold_scan=trial,gate_width_scan=width_test,
        required_rates= requirements,ordinary_charge_normalization_controls=matched,
        errors=dict(halved_timestep_exposure=float(fine['J'][-1]-J[-1]),
                    decimated_history_exposure=float(coarse['J'][-1]-J[-1]),
                    decimated_history_threshold=coarse['start']-main['start']))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source',type=Path,default=HERE.parent/'prior/CHARGED.json')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(calculate(args.source),indent=2,allow_nan=False)+'\n')
