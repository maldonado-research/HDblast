"""Static Euclidean kernel enclosures on the tuned leaves for scalar localization weights."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq
from scipy.special import betainc, expit, hyp2f1

ROOT = Path(__file__).resolve().parents[1]
K = 5/9
EXPONENTS = (0.25, K, 1., 2., 4.)
MOMENTA = (0.0001, 0.001, 0.01, 0.1)
LOG2 = math.log(2.)


def integrate(fn, lo, hi):
    return float(quad(fn, lo, hi, epsabs=3e-13, epsrel=2e-11, limit=200)[0])


def log_weight(y, r):
    return 2*r*(LOG2-float(np.logaddexp(y,-y)))


def log_tail(y, r):
    # Incomplete beta written with the explicit power outside to avoid underflow.
    u = float(expit(-2*y))
    logu = -float(np.logaddexp(0., 2*y))
    return (2*r-1)*LOG2 + r*logu - math.log(r) + math.log(hyp2f1(r,1-r,r+1,u))


def warp(y):
    return -(y+(2/3)*(float(np.logaddexp(y,-y))-LOG2)+math.tanh(y)**2/6)/3


def selected_log_tail(y):
    # v=(1-tanh y)/2, then v=u t^(1/k) removes the integrable endpoint singularity.
    u=float(expit(-2*y))
    logu=-float(np.logaddexp(0.,2*y))
    def integrand(t):
        v=u*t**(1/K)
        return (1-v)**(-10/9)*math.exp(-(1-2*v)**2/9)
    return -K*LOG2+K*logu-math.log(K)+math.log(integrate(integrand,0.,1.))


def finite_bounds(b, r, q, n, j):
    # R(q) places the test layer a logarithmic proper depth into the throat.
    end = b+math.log1p(1/q)/K
    ab = warp(b)
    # r=None denotes the separate null-lift weight in the same normalized metric.
    weight = (lambda y:2*warp(y)) if r is None else (lambda y:log_weight(y,r))
    tail = selected_log_tail if r is None else (lambda y:log_tail(y,r))
    def log_stiffness(y):
        return weight(y)+2*(warp(y)-ab)
    def j_integrand(y):
        return math.exp(2*tail(y)-log_stiffness(y))
    def h_integrand(y):
        return math.exp(tail(y)-log_stiffness(y))
    jt = integrate(j_integrand,b,end)
    h = integrate(h_integrand,b,end)
    omitted = math.exp(tail(end))
    nr = n-omitted
    # The finite-Neumann bound only needs J_R <= jt. This form also avoids
    # cancellation between two nearly equal unit-normalized kernels at small q.
    upper = (omitted+nr*(q*q*jt/(nr+q*q*jt)))/n
    if j is not None:
        upper = min(upper,q*q*j/(n+q*q*j))
    lower = q*q*jt/(n*(1+q*q*h))
    # A slightly shallower test gives the exact leading coefficient in the marginal case.
    log_scale=math.log(1/q)
    if log_scale>1.:
        test_end=b+(log_scale-.5*math.log(log_scale))/K
        test_j=integrate(j_integrand,b,test_end)
        test_h=integrate(h_integrand,b,test_end)
        lower=max(lower,q*q*test_j/(n*(1+q*q*test_h)))
    if not 0 < lower <= upper < 1:
        raise ArithmeticError('Invalid analytic enclosure')
    return dict(q=q,end=end,relative_loss_lower=lower,relative_loss_upper=upper)


def independent_enclosure(b, r, q, n):
    """Dirichlet and Neumann bracketing at a distant wall via bounded reflection variable."""
    ab = warp(b)
    end = b+math.log1p(12*K/q)/K+2.
    def rhs(y,s):
        phi = -math.tanh(y)
        aprime = -(1-phi+phi**3/3)/3
        log_derivative = 3*aprime if r is None else -2*r*math.tanh(y)+aprime
        return 2*q*math.exp(ab-warp(y))*s-log_derivative*(1-s*s)/2
    loss = []
    for value in (1.,-1.):
        out = solve_ivp(rhs,(end,b),[value],method='DOP853',rtol=2e-12,atol=1e-14,
                        max_step=.07)
        if not out.success:
            raise ArithmeticError(out.message)
        s = out.y[0,-1]
        weight = math.exp(2*ab if r is None else log_weight(b,r))
        kernel = q*weight*(1+s)/(1-s)
        loss.append(1-kernel/(q*q*n))
    return dict(cutoff=end,relative_loss_dirichlet=loss[0],relative_loss_neumann=loss[1])


def leaf(phi,r):
    b = -math.atanh(phi)
    ab = warp(b)
    n = math.exp(log_tail(b,r))
    nq = integrate(lambda y:math.exp(log_weight(y,r)),b,np.inf)
    if abs(n-nq)>2e-11*n:
        raise ArithmeticError('Incomplete beta and direct gauge norm differ')
    j = None
    if r>K:
        j = integrate(lambda y:math.exp(2*log_tail(y,r)-log_weight(y,r)-2*(warp(y)-ab)),
                      b,np.inf)
    points = {}
    for q in MOMENTA:
        row = finite_bounds(b,r,q,n,j)
        if q>=.01:
            direct=independent_enclosure(b,r,q,n)
            lo,hi = direct['relative_loss_dirichlet'],direct['relative_loss_neumann']
            # Numerical integration tolerance here allows for reflection-ratio cancellation.
            if hi-lo > 2e-8 or lo < row['relative_loss_lower']-2e-8 or hi > row['relative_loss_upper']+2e-8:
                raise ArithmeticError('Independent Dirichlet/Neumann check outside enclosure')
            row['independent']=direct
        points[str(q)]=row
    marginal=None
    if r==K:
        warp_offset=(2/9)*LOG2-1/18
        coefficient=math.exp(2*r*LOG2+2*(ab-warp_offset))/(4*K**3*n)
        marginal=marginal_check(b,r,n,coefficient)
    return dict(phi=phi,single_copy_norm=n,marginal=marginal,check_norm_absolute_error=abs(n-nq),
                j_over_norm=j/n if j is not None else None,
                one_percent_kernel_sufficient_q=math.sqrt(n/(99*j)) if j is not None else None,
                points=points)


def contrast(r,phi_i,phi_f):
    return .5*(math.log(betainc(r,r,(1+phi_f)/2))-math.log(betainc(r,r,(1+phi_i)/2)))


def marginal_check(b,r,n,coefficient):
    q=1e-40
    scale=q*q*math.log(1/q)
    limits=finite_bounds(b,r,q,n,None)
    lo,hi=limits['relative_loss_lower']/scale,limits['relative_loss_upper']/scale
    if max(abs(lo/coefficient-1),abs(hi/coefficient-1))>.1:
        raise ArithmeticError('Deep infrared coefficient disagrees with variational bounds')
    return dict(coefficient=coefficient,q=q,scaled_loss_lower=lo,scaled_loss_upper=hi)


def null_lift(history):
    branches={}
    for name,branch in history.items():
        ends={}
        for event in ('bounce','first_minimum'):
            entry=branch[event]
            b=-math.atanh(entry['phi_b'])
            ab=warp(b)
            n=math.exp(selected_log_tail(b))
            direct=integrate(lambda y:math.exp(2*warp(y)),b,np.inf)
            f=2*math.exp(-2*ab)*n
            if abs(n-direct)>2e-11 or abs(f-entry['M4_squared_rolling_leaf'])>1e-9:
                raise ArithmeticError('Warp-selected kinetic norm differs from background')
            coefficient=math.exp(2*ab)/(4*K**3*n)
            if abs(coefficient-1/(2*K**3*f))>2e-13:
                raise ArithmeticError('Marginal coefficient inconsistently normalized')
            checks={}
            for q in (0.01,0.1):
                analytic=finite_bounds(b,None,q,n,None)
                numerical=independent_enclosure(b,None,q,n)
                lower,upper=numerical['relative_loss_dirichlet'],numerical['relative_loss_neumann']
                if (upper-lower>2e-8 or lower<analytic['relative_loss_lower']-2e-8 or
                        upper>analytic['relative_loss_upper']+2e-8):
                    raise ArithmeticError('Null-lift kernels outside finite-depth enclosure')
                checks[str(q)]=dict(analytic=analytic,independent=numerical)
            ends[event]=dict(phi_b=entry['phi_b'],global_A=ab,norm=n,gravity_norm=f,
                             marginal_relative_loss_coefficient=coefficient,
                             marginal=marginal_check(b,None,n,coefficient),checks=checks)
        delta=.5*math.log(ends['first_minimum']['norm']/ends['bounce']['norm'])
        branches[name]=dict(endpoints=ends,log_I_contrast=delta)
    return dict(global_reference='A(0)=0; G(y)=exp(2A(y)), independent of shell',
                canonical_asymptotic_power=1.5,continuum_threshold=0.,branches=branches)


def calculate(archive):
    history=json.loads(archive.read_bytes())['branches']
    endpoints={key:(data['bounce']['phi_b'],data['first_minimum']['phi_b']) for key,data in history.items()}
    positive='phi_m=0.5'
    negative='phi_m=-0.5'
    r3=brentq(lambda r:contrast(r,*endpoints[negative])+3,2.,90.,xtol=2e-13)
    examples={}
    for r in (*EXPONENTS,8.,r3,32.):
        row={key:contrast(r,*pair) for key,pair in endpoints.items()}
        for key,pair in endpoints.items():
            qc=[integrate(lambda y:math.exp(log_weight(y,r)),-math.atanh(phi),np.inf) for phi in pair]
            if abs(.5*math.log(qc[1]/qc[0])-row[key])>3e-10:
                raise ArithmeticError('Independent contrasts differ')
        if not 0<row[positive]<LOG2/2 or not row[negative]<0:
            raise ArithmeticError('Sign or half-line capacity failed')
        if r==2. and max(abs(row[key]-history[key]['response']['log_I_change']) for key in endpoints)>1e-8:
            raise ArithmeticError('Failed to reproduce prior rolling contrast')
        examples[format(r,'.14g')]=row
    return dict(curvature_k=K,units='Q in inverse L0, warp factor set to one on each fixed leaf',
                family='G=(1-phi^2)^r=sech(y)^(2r)',launch_and_minimum=endpoints,
                positive_contrast_upper_bound=LOG2/2,
                negative_magnitude_three_r=r3,
                positive_contrast_at_negative_magnitude_three=contrast(r3,*endpoints[positive]),
                contrasts=examples,null_lift=null_lift(history),
                leaves={format(r,'.14g'):{str(phi):leaf(phi,r) for phi in (0.,.5,-.5)}
                        for r in EXPONENTS})


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior',type=Path,default=ROOT/'prior')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.prior/'origin'/'ORIGIN.json'),indent=2,allow_nan=False)+'\n')
