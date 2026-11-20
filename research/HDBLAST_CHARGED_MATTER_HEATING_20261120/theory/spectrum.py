"""Finite-temperature, charged spectator distribution on a fixed HDblast history."""
import argparse
import csv
import gzip
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import entr, expit

DEFAULT_TRACE = Path(__file__).resolve().parents[1] / 'background' / 'trace.csv.gz'
# F_j(eta) = integral y^j / (exp(y-eta)+1) dy.
def moment(j, eta=0.0):
    return quad(lambda y: y**j * expit(eta-y), 0., np.inf, epsabs=2e-12, epsrel=5e-13)[0]


F2 = moment(2)
F3 = moment(3)


def integrate_radial(fn, s):
    integrand = lambda y: y*y*expit(-y)*fn(y)
    return (quad(integrand, 0., s, epsabs=2e-12, epsrel=1e-12)[0]
            + quad(integrand, s, np.inf, epsabs=2e-12, epsrel=1e-12)[0])


def drift(s):
    """Energy / initial energy and signed current / maximal charge current."""
    if s == 0.:
        return 1., 0.
    u = abs(s)
    def energy(y):
        return u + y*y/(3*u) if y <= u else y + u*u/(3*y)
    def current(y):
        return 1-y*y/(3*u*u) if y <= u else 2*u/(3*y)
    return integrate_radial(energy, u)/F3, math.copysign(integrate_radial(current, u)/F2, s)


def target(energy_ratio):
    # Charge symmetry and separate number conservation give a common kinetic fugacity.
    def T(eta):
        return (F2/moment(2, eta))**(1./3.)
    def excess(eta):
        return T(eta)**4*moment(3, eta)/F3-energy_ratio
    if abs(energy_ratio-1.) < 2e-14:
        return 0., 1.
    eta = brentq(excess, -45., 0., xtol=1e-13)
    return eta, T(eta)


def angular_witness(s, energy_ratio, order):
    """Comoving radial-particle-number total variation. The thermal fit matches E and n."""
    chemical, temp = target(energy_ratio)
    mu, weight = np.polynomial.legendre.leggauss(order)
    def occupations(y):
        return expit(-np.sqrt((y-s)**2+2*y*s*(1-mu)))
    def mean(y):
        return .5*np.dot(weight,occupations(y))
    def deviation(y):
        return y*y*abs(mean(y)-expit(chemical-y/temp))/(2*F2)
    split = sorted(set([0.,s,*[1.*i for i in range(1,31)],60.,np.inf]))
    tv = sum(quad(deviation, l, h, epsrel=5e-10,epsabs=2e-11)[0] for l,h in zip(split,split[1:]))
    # Independent energy/current integrals in final (rather than initial) momenta.
    def total(fn):
        return quad(fn, 0., s,epsabs=2e-11,epsrel=1e-10)[0]+quad(fn,s,np.inf,epsabs=2e-11,epsrel=1e-10)[0]
    e = total(lambda y:y**3*mean(y))/F3
    n = total(lambda y:y*y*mean(y))/F2
    j = total(lambda y:y*y*.5*np.dot(weight*mu,occupations(y)))/F2
    return dict(relative_number_TV=tv,kinetic_chemical_potential_over_T=chemical,
                kinetic_comoving_T_over_initial=temp,number_ratio_num=n,
                energy_ratio_num=e,current_ratio_num=j)


def entropy_per_particle(s, eta, temperature, angular_depth, kinetic_depth):
    """Binary Fermi entropy in units of k_B, summed over occupied phase-space cells."""
    mu,weights=np.polynomial.legendre.leggauss(96)
    w_a,w_k=math.exp(-angular_depth),math.exp(-kinetic_depth)
    def integrand(y):
        f=expit(-np.sqrt((y-s)**2+2*y*s*(1-mu)))
        average=.5*np.dot(weights,f)
        relaxed=w_k*(w_a*f+(1-w_a)*average)+(1-w_k)*expit(eta-y/temperature)
        return y*y*.5*np.dot(weights,entr(relaxed)+entr(1-relaxed))
    return (quad(integrand,0,s,epsrel=2e-9,epsabs=1e-12)[0]
            +quad(integrand,s,np.inf,epsrel=2e-9,epsabs=1e-12)[0])/F2


def calculate(path):
    with gzip.open(path,'rt',newline='') as saved:
        trace=[{k:float(r[k]) for k in ('t','a','H','Bmag','sheet','e')} for r in csv.DictReader(saved)]
    def A(row):
        return row['a']**2*row['Bmag']/2.
    def sample(left,right,x1,x2):
        u=-x1/(x2-x1)
        return {k:left[k]+u*(right[k]-left[k]) for k in ('t','a')}
    start=end=None
    crossings=[]
    for left,right in zip(trace,trace[1:]):
        if left['t']>=4.51 and A(left)*A(right)<0:
            crossings.append(dict(bracket=[left['t'],right['t']],linear_time=sample(left,right,A(left),A(right))['t']))
        if right['t']<4.51:  # select expanding interval after completed scalar exit
            continue
        d1=left['sheet']-.9*left['e']
        d2=right['sheet']-.9*right['e']
        if d1<0<=d2 and right['H']>0 and start is None:
            start=sample(left,right,d1,d2)
        if left['H']>0>=right['H']:
            end=sample(left,right,left['H'],right['H'])
    if start is None or end is None or not crossings:
        raise ValueError('Cannot find the late radiation/turnaround landmarks')
    N=math.log(end['a']/start['a'])
    selected={}
    for t in (.2,4.5,8.1,12.4):
        row=min(trace,key=lambda r:abs(r['t']-t))
        assert abs(row['t']-t)<5e-4
        av=A(row)-A(trace[0])
        controls=[]
        for charge_per_T in (1.,10.,100.,300.):
            shift=charge_per_T*av
            ratio,current=drift(shift)
            controls.append(dict(qeff_over_initial_comoving_T=charge_per_T,
                       signed_comoving_momentum_shift_over_initial_T=shift,energy_factor=ratio,
                       current_over_q_times_total_number=current,
                       number_relative_to_zero_mu_thermal_at_same_energy=ratio**(-.75)))
        selected[str(t)]=dict(t=row['t'],a=row['a'],A=av,charge_controls=controls)
    s=selected['8.1']['charge_controls'][2]['signed_comoving_momentum_shift_over_initial_T']
    L,J=drift(s)
    result = angular_witness(s,L,96)
    independent = angular_witness(s,L,48)
    assert max(abs(result['number_ratio_num']-1),abs(result['energy_ratio_num']-L),
               abs(result['current_ratio_num']-J),abs(result['relative_number_TV']-independent['relative_number_TV']))<2e-6
    dL=(drift(s+1e-4)[0]-drift(s-1e-4)[0])/(2e-4)
    assert abs(dL-F2/F3*J)<1e-9
    half=end['t']-selected['8.1']['t']
    cases={}
    eta=result['kinetic_chemical_potential_over_T']
    temperature=result['kinetic_comoving_T_over_initial']
    entropy_zero=4*F3/(3*F2)
    entropy_eq=4*F3*L/(3*F2*temperature)-eta
    for name,ang,kin in [('collisionless',0.,0.),('angular_only',math.log(100),0.),
                    ('kinetic_only',0.,math.log(100)),('same_current_two_rates',math.log(10),math.log(10))]:
        entropy=entropy_per_particle(s,eta,temperature,ang,kin)
        if not entropy_zero-1e-6<entropy<=entropy_eq+1e-6:
            raise ValueError('Fermi entropy outside the kinetic bounds')
        cases[name]=dict(fermi_entropy_per_particle=entropy,
            remaining_entropy_deficit_per_particle=entropy_eq-entropy,angular_optical_depth=ang,kinetic_optical_depth=kin,
            current_relative_to_initial=math.exp(-ang-kin),remaining_relative_radial_error=math.exp(-kin),
            final_current_over_q_times_total_number=J*math.exp(-ang-kin),
            final_radial_number_TV=result['relative_number_TV']*math.exp(-kin))
    return dict(collisionless_zero_shift_crossings=crossings,
        input_source='dense forward wall trace of the published curved Maxwell wall model',
        source_time_unit='bulk curvature length at c=1',
        model='separate fixed-background test-charge sectors, initial zero-chemical-potential massless Fermi gas',
        selected_collisionless=selected,preparation=dict(t=selected['8.1']['t'],qeff_over_initial_comoving_T=100.,
            full_expanding_em90_elapsed_efolds=N,em90_start_time_linear_density=start['t'],
            maximum_time=end['t'],elapsed_post_preparation_time=half,
            full_em90_residual_if_Gamma_equals_H=math.exp(-N),
            mean_Gamma_over_H_for_99pct_full_em90=math.log(100)/N,
            constant_rate_for_99pct_post_preparation=math.log(100)/half,
            charge_shift=s,number_relative_to_zero_mu_thermal_at_same_energy=L**(-.75),
            fractional_extra_pairs_for_zero_mu_at_fixed_energy=L**.75-1,
            initial_exact_fermi_entropy_per_particle=entropy_zero,
            matched_equilibrium_entropy_per_particle=entropy_eq,
            fitted_and_checked=result,quadrature_order_sensitivity=abs(result['relative_number_TV']-independent['relative_number_TV']),
            work_identity_error=abs(dL-F2/F3*J)),controlled_relaxation=cases)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',default=DEFAULT_TRACE,type=Path,help='Dense primary trace from the study background directory')
    parser.add_argument('--output',required=True,type=Path)
    cfg=parser.parse_args()
    cfg.output.write_text(json.dumps(calculate(cfg.input),indent=2,sort_keys=True)+'\n')
