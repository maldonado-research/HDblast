"""Prepare time-smeared low-energy states on the closed rolling backgrounds."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre
from scipy.optimize import brentq
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parent
WIDTHS = (.10, .20, .30)
REF = .20


def load(prior):
    spec = importlib.util.spec_from_file_location('archived_homogeneous', prior/'homogeneous/exchange.py')
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.branches(prior/'prior/gravity'), json.loads((prior/'homogeneous/results.json').read_bytes())['branches']


class Background:
    def __init__(self, branch, archive):
        self.branch = branch
        self.arch = archive
        self.fit = branch.fitted
        self.ftop = self.fit.spline(branch.top)[1]
        offset = self.fit.sign*self.fit.dq
        self.f0=self.fit.spline(branch.top+offset)[1]
        self.a0 = 1/np.sqrt(self.potential(offset)[0])
        self.xstar = archive['runs']['inverse_f_seed_0']['first_minimum_crossing']['x']
        def dynamics(t, y):
            q, v, a, H, eta, tau = y
            u, du = self.potential(q)
            return v, -3*H*v-3*du, a*H, 1/a**2-v*v/2, 1/a, np.sqrt(self.f0/self.fit.spline(branch.top+q)[1])
        self.solve = solve_ivp(dynamics,(0.,self.xstar+.012),(offset,0.,self.a0,0.,0.,0.),
                               method='DOP853',max_step=.002,rtol=2e-13,atol=1e-15,dense_output=True)
        if not self.solve.success:
            raise ArithmeticError(self.solve.message)
        q,v,a,H,_,_=self.solve.sol(self.xstar)
        cross=archive['runs']['inverse_f_seed_0']['first_minimum_crossing']
        self.cross_error=float(max(abs(q-cross['q']),abs(v-cross['v']),abs(a-cross['a']),abs(H-cross['h'])))
        if self.cross_error > 2e-8:
            raise ArithmeticError('Background disagrees with archived crossing')
        # Supply the oscillator RHS from a separate dense table; tested against the
        # integrated background at interleaved times. Rows are eta,H,s=d(log Z)/dx.
        table_t=np.linspace(0.,self.xstar+.011,4097)
        _,H,eta,s=self.terms(table_t,True)
        self.fast=CubicSpline(table_t,np.column_stack((eta,H,s)),axis=0)
        mid=(table_t[:-1]+table_t[1:])/2
        _,HH,ee,ss=self.terms(mid,True)
        self.table_error=float(np.max(abs(self.fast(mid)-np.column_stack((ee,HH,ss)))))
        if self.table_error>5e-9:
            raise ArithmeticError('Oscillator background interpolation failed')

    def potential(self,q):
        u,du,_ = self.fit.at(self.branch.top+q)
        return u-self.branch.du_top*q,du-self.branch.du_top

    def terms(self,t,matter):
        q,v,a,H,eta,_=self.solve.sol(t)
        if matter:
            f=self.fit.spline(self.branch.top+q)[...,1]
            s=-v*self.fit.spline(self.branch.top+q,1)[...,1]/f
        else:
            s=0.*v
        return a,H,eta,s


def quadrature(lo,hi,count):
    z,w=roots_legendre(count)
    return (lo+hi)/2+(hi-lo)*z/2,(hi-lo)*w/2


def grids(bg,matter,order):
    xstar=bg.xstar
    samples={'bounce':0.,'probe':.12,'before':xstar-.005,'first_passage':xstar,'after':xstar+.005}
    smear={}
    for T in WIDTHS:
        x,w=quadrature(0.,T,order)
        weights=w*np.exp(-1/(1-(x/T)**2))
        smear[f'{T:.2f}']=(x,weights/np.sum(weights),np.ones_like(x))
    # Normalize the matter metric to agree with the Einstein metric at the bounce.
    # w d(tau_J) times the J-frame physical density is w sqrt(B) dx times rho_E/B².
    if matter:
        end=brentq(lambda x:bg.solve.sol(x)[5]-REF,0.,.4)
        x,w=quadrature(0.,end,order)
        q,_,_,_,_,tau=bg.solve.sol(x)
        ratio=bg.f0/bg.fit.spline(bg.branch.top+q)[...,1]
        weights=w*np.sqrt(ratio)*np.exp(-1/(1-(tau/REF)**2))
        smear['matter_clock_0.20']=(x,weights/np.sum(weights),1/ratio**2)
    ledger={}
    for label,(lo,hi) in {'bounce_to_passage':(0.,xstar),'passage_window':(xstar-.005,xstar+.005)}.items():
        ledger[label]=quadrature(lo,hi,order*2 if lo==0. else max(64,order//4))
    times=np.unique(np.concatenate([list(samples.values()),*[v[0] for v in smear.values()],
                                    *[v[0] for v in ledger.values()]]))
    return times,samples,smear,ledger


def propagate(bg,matter,labels,times,rtol):
    # Factor out the rapid phase to integrate Bogoliubov amplitudes with small error in B.
    k=np.sqrt(labels*(labels+2.))
    size=len(k)
    def rhs(t,coeffs):
        eta,H,s=bg.fast(t)
        factor=(H+(s/2 if matter else 0.))*np.exp(2j*k*eta)
        A,B=coeffs[:size],coeffs[size:]
        return np.r_[factor*B.conj(),factor*A.conj()]
    start=np.r_[np.ones(size,dtype=complex),np.zeros(size,dtype=complex)]
    sol=solve_ivp(rhs,(0.,times[-1]),start,method='DOP853',rtol=rtol,atol=rtol*1e-3,
                  t_eval=times,max_step=min(.01,2.8/k[-1]))
    if not sol.success:
        raise ArithmeticError(sol.message)
    A,B=sol.y[:size].T,sol.y[size:].T
    return A,B,float(np.max(abs(abs(A)**2-abs(B)**2-1)))


def shells(bg,matter,labels,times,smear,rtol):
    A,B,error=propagate(bg,matter,labels,times,rtol)
    k=np.sqrt(labels*(labels+2.))
    a,H,eta,s=bg.terms(times,matter)
    C=1+2*abs(B)**2
    D=2*(A*B.conj()).real
    V=2*np.pi**2*a**3
    # Both energy terms have common shell-dependent k/(2a), hence weight is w/a^4.
    u={}
    for T,(x,w,density_factor) in smear.items():
        i=np.searchsorted(times,x)
        measure=w*density_factor/a[i]**4
        u[T]=-np.arctanh((measure@D[i])/(measure@C[i]))
    # Per-mode oscillator energy and potential-minus-kinetic imbalance (expectation differences).
    phase=np.exp(-2j*eta[:,None]*k)
    X=(phase*A*B).real
    Y=.5*(phase*(A*A+B*B)).real
    omega=k/a[:,None]
    deg=(labels+1)**2
    results={}
    for T in smear:
        if T==f'{REF:.2f}':
            continue
        delta=(u[T]-u[f'{REF:.2f}'])/2
        mean=(u[T]+u[f'{REF:.2f}'])/2
        dh=2*np.sinh(delta)*np.sinh(mean)
        ds=2*np.sinh(delta)*np.cosh(mean)
        eng=.5*omega*(C*dh+D*ds)*deg
        imbalance=omega*(X*dh+Y*ds)*deg
        results[T]={'rho':eng/V[:,None], 'pressure':(eng-2*imbalance)/(3*V[:,None]),
                    'comoving_energy':eng, 'power':s[:,None]*imbalance,
                    'metric_power':-H[:,None]*(eng-2*imbalance)}
    return u,results,error


def run_case(bg,matter,nmax,order,rtol,cutoffs,block=128):
    times,samples,smear,ledger=grids(bg,matter,order)
    at={name:int(np.searchsorted(times,t)) for name,t in samples.items()}
    max_w=0.
    all_rows={key:{} for key in smear if key!=f'{REF:.2f}'}
    low={}
    logwidths={T:[] for T in smear}
    quantities={}
    for lo in range(1,nmax+1,block):
        labels=np.arange(lo,min(lo+block,nmax+1))
        u,data,err=shells(bg,matter,labels,times,smear,rtol)
        max_w=max(max_w,err)
        for T in smear: logwidths[T].extend(u[T])
        for T,values in data.items():
            if T not in quantities:
                quantities[T]={}
            series=quantities[T]
            for tag,i in at.items():
                for field in ('rho','pressure','comoving_energy'):
                    series.setdefault((tag,field),[]).extend(values[field][i])
            for tag,(x,w) in ledger.items():
                indices=np.searchsorted(times,x)
                for field in ('power','metric_power'):
                    series.setdefault((tag,field),[]).extend(w@values[field][indices])
            for tlabel,(x,w,fac) in smear.items():
                i=np.searchsorted(times,x)
                series.setdefault((tlabel,'smear'),[]).extend((w*fac)@values['rho'][i])
    for label in (1,2,4,8,16,32,64,128,256,512):
        if label<=nmax:
            low[str(label)]={T:float(logwidths[T][label-1]) for T in smear}
    for T,series in quantities.items():
        output=all_rows[T]
        for N in cutoffs:
            row={'samples':{},'work':{},'sampled_density':{}}
            def total(stage,key):
                return float(np.sum(series[stage,key][:N]))
            for tag,t in samples.items():
                q,v,*_=bg.solve.sol(t)
                frame_ratio=bg.f0/bg.fit.spline(bg.branch.top+q)[1] if matter else 1.
                rho=total(tag,'rho')
                parent=bg.potential(q)[0]+v*v/6
                row['samples'][tag]=dict(x=t,rho_difference_over_H0_4=rho,
                    pressure_difference_over_H0_4=total(tag,'pressure'),
                    matter_measured_density_difference_over_H0_4=rho/frame_ratio**2,
                    density_difference_over_parent=rho*bg.fit.Hpl**2/(3*parent))
            for tag,(first,last) in {'bounce_to_passage':('bounce','first_passage'),
                                     'passage_window':('before','after')}.items():
                energy=total(last,'comoving_energy')-total(first,'comoving_energy')
                Wz=total(tag,'power')
                Wg=total(tag,'metric_power')
                row['work'][tag]=dict(comoving_energy_change_over_H0=energy,
                    coupling_work_over_H0=Wz,geometry_work_over_H0=Wg,residual=energy-Wz-Wg)
            for T2 in smear:
                row['sampled_density'][T2]=total(T2,'smear')
            # Integrated energy of minimizing state minus same functional on reference is nonpositive.
            if row['sampled_density'][T]>2e-12 or row['sampled_density'][f'{REF:.2f}'] < -2e-12:
                raise ArithmeticError('State failed its energy-minimization check')
            output[str(N)]=row
    if max_w>1e-9:
        raise ArithmeticError('Canonical evolution failed')
    return dict(maximum_mode=nmax,quadrature_order=order,rtol=rtol,wronskian_maximum=max_w,
                selected_log_covariance_widths=low,relative_to_width_0p20=all_rows)



def check_de_sitter():
    """Compare the rolling solver's mode convention to exact unit-H closed de Sitter."""
    class ClosedDeSitter:
        @staticmethod
        def fast(t):
            return np.arctan(np.sinh(t)),np.tanh(t),0.

        @staticmethod
        def terms(t,matter):
            return np.cosh(t),np.tanh(t),np.arctan(np.sinh(t)),np.zeros_like(t)

    labels=np.array([1,2,3,4,7])
    k=np.sqrt(labels*(labels+2.))
    L=labels+1
    bg=ClosedDeSitter()
    windows={}
    for T in WIDTHS:
        x,w=quadrature(0.,T,96)
        prob=w*np.exp(-1/(1-(x/T)**2))
        windows[f'{T:.2f}']=x,prob/np.sum(prob),np.ones_like(x)
    t=np.unique(np.concatenate([x for x,_,_ in windows.values()]))
    u,_,wronskian=shells(bg,False,labels,t,windows,2e-12)
    rows={}
    error=0.
    for key,(x,w,_) in windows.items():
        a=np.cosh(x)[:,None]
        eta=np.arctan(np.sinh(x))[:,None]
        phase=np.exp(-1j*eta*L)/np.sqrt(2*L*k**2)
        field=(L*np.cos(eta)+1j*np.sin(eta))*phase
        derivative=-1j*k**2*np.cos(eta)*phase/a
        omega=k/a
        qquad=w@(derivative.real**2+(omega*field.real)**2)
        pquad=w@(derivative.imag**2+(omega*field.imag)**2)
        exact_u=np.log(L/k)+.5*np.log(pquad/qquad)
        error=max(error,float(np.max(abs(exact_u-u[key]))))
        rows[key]={'from_closed_modes':list(exact_u),'from_integrator':list(u[key])}
    if error>2e-12:
        raise ArithmeticError('Physical-density smearing failed closed de Sitter control')
    return dict(shells=list(map(int,labels)),h=1.,log_widths=rows,
                maximum_absolute_error=error,wronskian_maximum=wronskian)


def run(prior,nmax=2048,order=584,rtol=2e-11,block=128):
    if nmax<1 or order<8 or block<1 or rtol<=0:
        raise ValueError("Expected positive nmax, block and tolerance and at least eight quadrature points")
    branches,archive=load(prior)
    out={'sampling':'Even compact proper-time density weights proportional to exp(-1/(1-(x/T)^2))',
         'time_units':'x=H0 t','half_widths':list(WIDTHS),'reference_half_width':REF,
         'zero_mode':'identical independent state across preparations, cancels in the differences',
         'closed_de_sitter_check':check_de_sitter(),'branches':{}}
    for name,b in branches.items():
        bg=Background(b,archive[name])
        data={'H0_over_Mpl':bg.fit.Hpl,'archived_crossing_x':bg.xstar,
              'maximum_archived_crossing_error':bg.cross_error,'background_table_error':bg.table_error,'couplings':{}}
        for matter in (True,False):
            cutoffs=sorted({n for n in (32,64,128,256,512,768,1024,1536,2048,nmax) if n<=nmax})
            run=run_case(bg,matter,nmax,order,rtol,cutoffs,block)
            data['couplings']['inverse_f' if matter else 'constant']=run
        out['branches'][name]=data
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior',type=Path,default=HERE.parent/'prior')
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--nmax',type=int,default=2048)
    parser.add_argument('--order',type=int,default=584)
    parser.add_argument('--rtol',type=float,default=2e-11)
    parser.add_argument('--block',type=int,default=128)
    args=parser.parse_args()
    result=run(args.prior,args.nmax,args.order,args.rtol,args.block)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
