"""Exact massive Maxwell–Jüttner momentum transfer between two equal counterstreams."""
import math

import numpy as np
from scipy.special import kve, roots_legendre

POLE=8300.
FRACTION=.30


def m2sigma(eta, pole=POLE, width=FRACTION):
    """Dimensionless m² sigma_T for one unpolarized elastic vector channel."""
    E=np.cosh(eta)
    s=4*E**2
    g=width*(s+2.)
    return 3.*math.pi*(g/E)**2/((s-pole**2)**2+(g*np.tanh(eta))**2)


def scaled_k(n, x):
    x=np.asarray(x,float)
    small=x<1e5
    out=np.empty(x.shape)
    out[small]=kve(n,x[small])
    xx=x[~small]
    a=4.*n*n
    h=1/(8*xx)
    out[~small]=np.sqrt(np.pi/(2*xx))*(1+(a-1)*h*(1+(a-9)*h/2*(1+(a-25)*h/3)))
    return out


def normalized_response(theta, r, eta):
    """e^(2/theta)*I/v in the pair-rapidity representation, including the v=0 limit."""
    E=2*np.cosh(eta)
    gamma=math.cosh(r); v=math.tanh(r)
    x=E*gamma/theta
    y_over_v=2*np.sinh(eta)*gamma/theta
    y=y_over_v*v
    minus=2*np.cosh(eta-r)/theta
    plus=2*np.cosh(eta+r)/theta
    log_minus=-4*np.sinh((eta-r)/2)**2/theta
    log_plus=-4*np.sinh((eta+r)/2)**2/theta
    use_series=y<.045*np.minimum(x,1.)
    result=np.empty_like(eta)
    j=~use_series
    if np.any(j):
        ma,pl=minus[j],plus[j]
        em,ep=np.exp(log_minus[j]),np.exp(log_plus[j])
        bracket=(scaled_k(1,ma)*em+scaled_k(1,pl)*ep
                 -(scaled_k(0,ma)*em-scaled_k(0,pl)*ep)/y[j])
        result[j]=E[j]**2*bracket/(2*y[j]*v)
    if np.any(use_series):
        j=use_series
        xx=x[j]; yy=y[j]
        output=np.zeros_like(xx)
        for i in range(5):
            p=2*i+3
            f=(2*i+2)/math.factorial(2*i+3)
            avg=sum(math.comb(p,k)*scaled_k(abs(p-2*k),xx) for k in range((p+1)//2))*2**(1-p)
            output+=f*(yy**(2*i))*avg
        result[j]=E[j]**2*y_over_v[j]*np.exp(log_minus[j]-yy)*output
    return result


def rate(theta, r, *, order=20, cross_section=m2sigma, tm=True):
    """Return the wall-time momentum drag nu*m²/n. r=abs(asinh(gamma*v))."""
    r=abs(r)
    if theta < 0: raise ValueError("Negative rest temperature")
    if theta==0:
        return float(2*math.sinh(r)*cross_section(np.array(r)))
    width=2*math.asinh(math.sqrt(60*theta)/2)
    a,b=max(0.,r-width),r+width
    if not a<b: return float(2*math.sinh(r)*cross_section(np.array(r)))
    center=math.acosh(POLE/2)
    special=[a,b,r,center-1,center-.25,center,center+.25,center+1]
    special=sorted(set(z for z in special if a<=z<=b))
    edges=[special[0]]
    step=min(1.5,3*math.sqrt(theta))
    for left,right in zip(special,special[1:]):
        pieces=max(1,math.ceil((right-left)/step))
        edges.extend(np.linspace(left,right,pieces+1)[1:])
    x,w=_nodes(order)
    edges=np.array(edges)
    half=np.diff(edges)/2
    eta=(edges[:-1,None]+(x+1)*half[:,None]).ravel()
    weights=(half[:,None]*w).ravel()
    k2=scaled_k(2,1/theta)
    h=(2.5*theta+math.hypot(1.,1.5*theta)) if tm else scaled_k(3,1/theta)/k2
    integrand=np.sinh(eta)**4*cross_section(eta)*normalized_response(theta,r,eta)
    return float((2/(theta*math.cosh(r)**3*h*k2**2))*(weights @ integrand))


from functools import lru_cache


@lru_cache(maxsize=None)
def _nodes(n):
    return roots_legendre(n)


def _table_row(args):
    theta,dy,ny,order=args
    u_scale=math.sqrt(theta/(1+theta))
    return [math.log(rate(theta,math.asinh(u_scale*math.sinh(j*dy)),order=order))
            for j in range(ny)]


def write_table(path,log_step=.125,dy=.05,order=16,jobs=1):
    from concurrent.futures import ProcessPoolExecutor
    t0=-5.
    nt=1+round(11/log_step); ny=1+round(20/dy)
    log_step=11/(nt-1); dy=20/(ny-1)
    tasks=[(10**(t0+i*log_step),dy,ny,order) for i in range(nt)]
    with open(path,'w') as f:
        f.write(f'{nt} {ny} {t0:.15g} {log_step:.15g} {dy:.15g}\n')
        if jobs==1:
            np.savetxt(f,list(map(_table_row,tasks)),fmt='%.12g')
        else:
            with ProcessPoolExecutor(max_workers=jobs) as processes:
                for row in processes.map(_table_row,tasks):
                    np.savetxt(f,[row],fmt='%.12g')


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True)
    parser.add_argument('--jobs',type=int,default=2)
    parser.add_argument('--log-step',type=float,default=.125)
    parser.add_argument('--dy',type=float,default=.05)
    parser.add_argument('--order',type=int,default=16)
    args=parser.parse_args()
    write_table(args.output,args.log_step,args.dy,args.order,args.jobs)
