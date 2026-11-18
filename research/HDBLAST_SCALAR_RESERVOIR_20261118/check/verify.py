"""Audit the archived paid cohort and a frozen-time full-insertion control, without running the producer."""
import argparse
import hashlib
import io
import json
import math
from pathlib import Path

import numpy as np
from scipy.special import roots_genlaguerre

PINS = {
    'dzf1e-3': ('2f2b61e745c438ef057533d75329f998621983a3e565226c9d725e4925e34716','3c9e732322715bc33d09eff6212066b15448da70c7ea10a177edc39e9f6ac9b2'),
    'dzf5e-4': ('bf162614a84fd207f2f2736b8b77cd1cb173c2db012469cf1d4fe1ccc50c8d0f','9eae8828655af2ec8cc91bb658ecf0e768ca9f54215b8927313937064f1b1646'),
}

def file(path, pin):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != pin:
        raise ValueError('source checksum: ' + path.name)
    return raw

def trapezoid(x, y):
    return math.fsum((u-t)*(q+p)/2 for t,u,p,q in zip(x,x[1:],y,y[1:]))

def signed_parts(x,y):
    pos, neg = [], []
    for left,right,p,q in zip(x,x[1:],y,y[1:]):
        if p*q<0:
            w=abs(p)/(abs(p)+abs(q))
            areas=[(right-left)*w*p/2,(right-left)*(1-w)*q/2]
        else:
            areas=[(right-left)*(p+q)/2]
        pos.extend(v for v in areas if v>0)
        neg.extend(v for v in areas if v<0)
    return math.fsum(pos),math.fsum(neg)

def one(root, spacing):
    tag = "main_ps0.5_G100_y1_bcut_" + spacing
    prefix = root/"runs/main"
    sm = json.loads(file(prefix/(tag+"_summary.json"),PINS[spacing][0]))
    with np.load(io.BytesIO(file(prefix/(tag+"_timeseries.npz"),PINS[spacing][1])),allow_pickle=False) as z:
        v = dict(zip(z["cols"].tolist(), z["rec"].T, strict=True))
    f=sm['params']['matter']; b=f['b']; rb=sm['rho_b']; e,=sm['matter']['events']
    if (f['mode'],f['sign_j'],f['G'],f['y'],f['phistar'],f['nk'],f['mu0'])!=('ramp',1,100,1,.5,24,0):
        raise ValueError(tag)
    s=v['H0tau']; a3=np.exp(3*(v['A_b']-v['A_b'][0]));a3star=math.exp(3*(e['A_star']-v['A_b'][0]))
    if not np.all(np.diff(s)>0):
        raise ValueError('invalid clock')
    nodes=np.r_[1.,s[(s>1)&(s<4)],4.]
    clip=lambda y: np.interp(nodes,s,y)
    def integral(y):return trapezoid(nodes,clip(y))
    def increment(y):return clip(y)[-1]-clip(y)[0]
    phi=v['phi_b']; vel=v['v_over_H0']; h=v['H_over_H0']; er=rb*a3*v['R']; ex=b*a3*v['rho_hat']
    w_mass=a3*b*v['j_hat']*vel; w_prod=a3*b*v['P_hat']
    source=a3*b*v['Q_hat']; c_r=h*er;c_x=3*h*a3*b*v['p_hat']
    im=integral(w_mass); ip=integral(w_prod); iq=integral(source); cr=integral(c_r); cx=integral(c_x)
    mpos,mneg=signed_parts(nodes,clip(w_mass)); ppos,pneg=signed_parts(nodes,clip(w_prod))
    # No use of dphi/ds: integrate J dphi as a Riemann-Stieltjes integral in time order.
    im_s=trapezoid(clip(phi),clip(a3*b*v['j_hat']))
    ip_s=trapezoid(clip(phi),clip(a3*b*v['jprod_hat']))
    jtot=integral(a3*rb*v['J0']*vel)
    u0=abs(e['phi_after']-f['phistar'])/e['ramp_dphi']; initial=u0**3*(10-15*u0+6*u0*u0)
    # Exact for the algorithm at the same post-crossing solver step: totals are linear
    # in the cohort occupation vector. Divide out f(u0); no alternative PDE run.
    recorded_jump=e['kappa5sq_rho_jump_over_H0']; instant_jump=recorded_jump/initial
    # At the precise crossing compute all quadrature mode momenta and full occupation energies.
    roots,weights=roots_genlaguerre(f['nk'],.5); q=e['q']
    momenta=np.sqrt(q*roots/math.pi)
    Ns=weights*(q/math.pi)**1.5/(4*math.pi**2)*e['number_factor']*math.exp(-math.pi*f['mu0']**2/q)
    ideal_24=b*float(np.dot(Ns,momenta)); ideal_continuum=b*e['number_factor']*q*q/(4*math.pi**4)
    ncal=a3star*e['n_hat_at_crossing']; nstart=increment(a3*v['n_hat'])
    nloss=integral(a3*f['y']**2*(v['rho_hat']-3*v['p_hat'])/(8*math.pi))
    gain=increment(er); chi=increment(ex);w=im+ip
    result={
     'window':[1,4], 'points':len(nodes), 'crossing':e['s_star'],
     'a_star_over_a0_cubed':a3star,'total_comoving_work':w,
     'mass_work':im,'mass_work_positive':mpos,'mass_work_negative':mneg,
     'production_work':ip,'production_work_positive':ppos,'production_work_negative':pneg,
     'mass_to_production':im/ip,'mass_fraction_of_work':im/w,
     'from_total_shell_J':jtot, 'from_Jdphi_mass_work':im_s,'from_Jdphi_production_work':ip_s,
     'from_Jdphi_total_relative_to_work':(im_s+ip_s-w)/w,
     'decay_transfer':iq, 'gross_decay_transfer_over_work':iq/w,'radiation_retained_gain_over_work':gain/w,
     'remaining_particle_gain_over_work':chi/w,'particle_pressure_work_over_work':cx/w,
     'radiation_pressure_work_over_work':cr/w,'integrated_all_budget_residual_over_work':(gain+chi+cx+cr-w)/w,
     'actual_ramp_u_first_post_crossing':u0, 'actual_ramp_occupation_fraction_at_first_step':initial,
     'recorded_unpaid_local_jump':recorded_jump,'same_step_unpaid_full_insertion_local_jump':instant_jump,
     'ideal_crossing_unpaid_full_insertion_local_24_nodes':ideal_24,'ideal_crossing_unpaid_continuum_massless_proxy':ideal_continuum,
     'paid_production_comoving_work_over_ideal_instant_comoving_jump':ip/(a3star*ideal_24),
     'after_completion_u_last':v['u_last'][-1],
     'last_sample_before_u1':float(s[np.flatnonzero(v['u_last']>=1)[0]-1]),'first_sample_after_u1':float(s[np.flatnonzero(v['u_last']>=1)[0]]),
     'comoving_theory_full_births':ncal,'net_comoving_number':nstart,'inferred_comoving_decays':nloss,
     'recovered_birth_count_over_full_minus_one':(nstart+nloss)/ncal-1,
     'radiation_gain':gain,'particle_gain':chi,'radiation_pressure_work':cr,'particle_pressure_work':cx
    }
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('archive',type=Path,help='Directory of the original B1 simulation archive')
    args=parser.parse_args()
    print(json.dumps({s:one(args.archive,s) for s in PINS},indent=2,allow_nan=False))
