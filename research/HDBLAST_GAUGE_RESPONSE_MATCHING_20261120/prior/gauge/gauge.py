"""Vacuum-to-photon transfer from a prescribed gauge coefficient on a closed FLRW roll."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq
from scipy.special import roots_legendre

HERE = Path(__file__).resolve().parent
DELTAS = np.array([0., .01, .1, .5, 1., -.1])
PRIMARY = 2
SELECTED = [2, 3, 4, 8, 16, 32, 48, 64, 80, 96, 128, 192, 256, 384, 512, 768, 1024]
CUTOFFS = [32, 64, 128, 256, 512, 768, 1024, 1536, 2048]
FRACTIONS = [.5, .9, .99, .999]


def load(prior):
    spec = importlib.util.spec_from_file_location('historical_exchange', prior/'homogeneous/exchange.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    archive = json.loads((prior/'homogeneous/results.json').read_text())['branches']
    return mod.branches(prior/'prior/gravity'), archive


def raw_bump(u):
    if u <= 0. or u >= 1.:
        return 0.
    return np.exp(-1/(u*(1-u)))


NORM = quad(raw_bump, 0., 1., epsabs=3e-16, epsrel=8e-14)[0]


class Roll:
    def __init__(self, branch, archive, mesh, limits=(.2,.8)):
        lower,upper=limits
        if not 0<lower<upper<1:
            raise ValueError("Switch limits must lie inside the directed roll")
        self.fit = branch.fitted
        self.Hpl = self.fit.Hpl
        self.q0 = self.fit.sign*self.fit.dq
        self.qmin = self.fit.tmin-branch.top

        def parent_potential(q):
            return self.fit.at(branch.top+q)[0]-branch.du_top*q

        def rhs(_, y):
            q, v, a, h, eta = y
            slope = self.fit.at(branch.top+q)[1]-branch.du_top
            return v, -3*h*v-3*slope, a*h, 1/a**2-v*v/2, 1/a

        initial = [self.q0, 0., 1/np.sqrt(parent_potential(self.q0)), 0., 0.]
        target = archive['runs']['inverse_f_seed_0']['first_minimum_crossing']
        cross = target['x']
        bg = solve_ivp(rhs, (0., cross+.001), initial, method='DOP853', max_step=.002,
                       rtol=2e-13, atol=1e-15, dense_output=True)
        if not bg.success:
            raise ArithmeticError(bg.message)
        self.parent_error = float(max(abs(bg.sol(cross)[:4]-[target[k] for k in ('q','v','a','h')])))
        if self.parent_error > 2e-8:
            raise ArithmeticError('The reconstructed roll misses the archived first minimum')

        def progress(x):
            return (bg.sol(x)[0]-self.q0)/(self.qmin-self.q0)

        def time_at(frac):
            return brentq(lambda x: progress(x)-frac, 0., cross+.001, xtol=6e-15)

        self.points = {}
        for name, fraction in [('launch',0.), ('switch_start',lower), ('switch_end',upper), ('first_minimum',1.)]:
            x = time_at(fraction) if fraction else 0.
            q,v,a,h,eta = bg.sol(x)
            parent = parent_potential(q)+v*v/6
            self.points[name] = dict(x=x, eta=eta, q=q, a=a, h=h, velocity=v, parent_over_V0=parent,
                                     photon_to_parent_per_D=self.Hpl**2/(3*a**4*parent))
        self.xstar_error = abs(self.points['first_minimum']['x']-cross)
        x = np.linspace(self.points['switch_start']['x'],self.points['switch_end']['x'],mesh)
        u = (progress(x)-lower)/(upper-lower)
        self.eta = CubicSpline(u,bg.sol(x)[4]-self.points['switch_start']['eta'])
        xm = (x[1:]+x[:-1])/2
        um = (progress(xm)-lower)/(upper-lower)
        self.phase_error = float(np.max(abs(self.eta(um)-(bg.sol(xm)[4]-self.points['switch_start']['eta']))))
        if self.phase_error > 1e-9:
            raise ArithmeticError('Inverse background interpolation failed')


def bogoliubov(roll, labels, deltas, rtol):
    """Integrate a dimensionless squeezing rate; alpha,beta phases use eta at switch-on."""
    shape = (2,len(deltas),len(labels))
    launch = np.zeros(shape, dtype=complex)
    launch[0] = 1.

    def rhs(u, y):
        alpha,beta = y.reshape(shape)
        osc = np.exp(2j*labels*roll.eta(u))
        pump = (deltas[:,None]*raw_bump(u)/NORM)*osc
        return np.stack((pump*beta.conj(),pump*alpha.conj())).ravel()

    phase_span = labels[-1]*(roll.eta(1.)-roll.eta(0.))
    sol = solve_ivp(rhs,(0.,1.),launch.ravel(),method='DOP853',rtol=rtol,atol=rtol*1e-4,
                    max_step=min(.025,1.7/phase_span),t_eval=[1.])
    if not sol.success:
        raise ArithmeticError(sol.message)
    return sol.y[:,-1].reshape(shape)


def direct_fields(roll, labels, delta):
    """Independent check with the vector potential Q and canonical momentum P=I^2 Q_eta."""
    labels = np.asarray(labels,dtype=float)
    size = len(labels)
    y0 = np.r_[np.array([0j]), 1/np.sqrt(2*labels), -1j*np.sqrt(labels/2)]
    def rhs(u,y):
        log_i = y[0].real
        q,p = y[1:1+size],y[1+size:]
        slope = roll.eta(u,1)
        i2 = np.exp(2*log_i)
        return np.r_[delta*raw_bump(u)/NORM,slope*p/i2,-slope*i2*labels**2*q]
    sol = solve_ivp(rhs,(0.,1.),y0,method='DOP853',rtol=2e-13,atol=1e-15,
                    max_step=.008,t_eval=[1.])
    if not sol.success:
        raise ArithmeticError(sol.message)
    y = sol.y[:,-1]
    if abs(y[0]-delta)>2e-12:
        raise ArithmeticError('Incorrect integrated log I')
    q,p = y[1:1+size],y[1+size:]
    return abs(np.sqrt(labels/2)*np.exp(delta)*q-1j*p/(np.exp(delta)*np.sqrt(2*labels)))**2


def born(roll,labels,order):
    z, w = roots_legendre(order)
    u, weights = (1+z)/2, w*np.exp(-1/(((1+z)/2)*((1-z)/2)))/(2*NORM)
    if abs(np.sum(weights)-1)>2e-12:
        raise ArithmeticError('Incorrect bump normalization')
    return weights@np.exp(2j*roll.eta(u)[:,None]*labels)


def spectrum(roll,nmax,rtol,block,deltas):
    labels=np.arange(2,nmax+1)
    mix=np.empty((len(deltas),len(labels)),dtype=complex)
    born_unit=np.empty(len(labels),dtype=complex)
    w_error=0.
    for start in range(0,len(labels),block):
        part=slice(start,start+block)
        alpha,beta=bogoliubov(roll,labels[part],deltas,rtol)
        mix[:,part]=beta
        born_unit[part]=born(roll,labels[part],max(160,int(nmax*.3)+80))
        w_error=max(w_error,float(np.max(abs(abs(alpha)**2-abs(beta)**2-1))))
    if w_error>1e-9:
        raise ArithmeticError('Bogoliubov normalization failed')
    return labels,mix,born_unit,w_error


def summarize(labels, occ, born_occ, roll):
    weights=labels*(labels**2-1)/np.pi**2
    d=weights*occ
    b=weights*born_occ
    cumulative=np.cumsum(d)
    total=float(cumulative[-1])
    quantiles={str(f):int(labels[np.searchsorted(cumulative,f*total)]) for f in FRACTIONS} if total else {}
    return dict(D=total,D_born=float(b.sum()),D_relative_Born_error=float((b.sum()-total)/total) if total else None,
                number_per_comoving_unit_volume=float(np.sum((labels**2-1)*occ)/np.pi**2),
                peak_energy_shell=int(labels[np.argmax(d)]) if total else None,
                cumulative_energy_shells=quantiles,
                D_by_cutoff={str(k):float(cumulative[k-2]) for k in sorted(set(CUTOFFS+[int(labels[-1])])) if k<=labels[-1]},
                photon_to_parent={k:total*roll.points[k]['photon_to_parent_per_D'] for k in ('switch_end','first_minimum')},
                selected_occupations={str(L):{'exact':float(occ[L-2]),'born':float(born_occ[L-2])}
                                      for L in SELECTED if L<=labels[-1]})


def run(prior,nmax,mesh,rtol,block):
    if nmax<128 or mesh<33 or rtol<=0 or block<1:
        raise ValueError('Use nmax >= 128, mesh >= 33 and positive solver controls')
    branches, archives = load(prior)
    result=dict(maximum_shell=nmax,background_mesh=mesh,rtol=rtol,block_size=block,
                bump_normalization=NORM,switch_field_fractions=[.2,.8],primary_delta=.1,branches={})
    for name,branch in branches.items():
        roll=Roll(branch,archives[name],mesh)
        labs,mix,born_unit,wr=spectrum(roll,nmax,rtol,block,DELTAS)
        occ=abs(mix)**2
        direct_labels=np.array([2,4,8,16,32,64,128])
        direct=direct_fields(roll,direct_labels,.1)
        error=float(max(abs(direct-occ[PRIMARY,direct_labels-2])))
        if error>1e-11 or np.any(mix[0]) or np.max(abs(mix[-1]+mix[PRIMARY]))>1e-12:
            raise ArithmeticError('Direct oscillator, zero-source or sign-reciprocity check failed')
        rows={str(d):summarize(labs,o,abs(d*born_unit)**2,roll) for d,o in zip(DELTAS,occ)}
        primary_spectrum=dict(first_shell=2,last_shell=int(labs[-1]),
                              occupation=list(occ[PRIMARY]),occupation_born=list(abs(.1*born_unit)**2),
                              cumulative_D=list(np.cumsum(labs*(labs**2-1)*occ[PRIMARY]/np.pi**2)))
        sensitivity={}
        # Supplementary cases change only the support width, around the same midpoint.
        for low,high in ((.15,.85),(.25,.75)):
            other=Roll(branch,archives[name],mesh,limits=(low,high))
            ls,betas,seed,_=spectrum(other,nmax,rtol,block,np.array([.1]))
            info=summarize(ls,abs(betas[0])**2,abs(.1*seed)**2,other)
            sensitivity[f'{low:.2f}-{high:.2f}']=dict(
                limits=[low,high],conformal_duration=other.points['switch_end']['eta']-other.points['switch_start']['eta'],
                D=info['D'],D_by_cutoff=info['D_by_cutoff'],peak_energy_shell=info['peak_energy_shell'],
                cumulative_energy_shells=info['cumulative_energy_shells'],photon_to_parent=info['photon_to_parent'])
        photon_unit=roll.points['first_minimum']['photon_to_parent_per_D']
        targets=dict(half_fixed_background_budget_D=.5/photon_unit,
                     equality_to_original_parent_density_D=1/photon_unit)
        result['branches'][name]=dict(points=roll.points,archived_crossing_error=roll.parent_error,
                       archived_crossing_time_error=roll.xstar_error,phase_interpolation_error=roll.phase_error,
                       H0_over_Mpl=roll.Hpl,wronskian_error=wr,sign_reciprocity_beta_error=float(np.max(abs(mix[-1]+mix[PRIMARY]))),
                       direct_mode_occupation_error=error,scan=rows,supplementary_support_test=sensitivity,first_minimum_target=targets,primary_spectrum=primary_spectrum)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior',type=Path,default=HERE.parent/'prior')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--nmax',type=int,default=1536)
    parser.add_argument('--mesh',type=int,default=2049)
    parser.add_argument('--rtol',type=float,default=2e-11)
    parser.add_argument('--block',type=int,default=96)
    args=parser.parse_args()
    data=run(args.prior,args.nmax,args.mesh,args.rtol,args.block)
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
