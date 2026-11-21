"""Wall-time single-particle event hazard between two thermal counterstreams."""
import math
from functools import lru_cache
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.special import erf, kve, roots_legendre


POLE = 8300.
WIDTH = .30


def sigma(eta):
    """m² times the total cross section for the isolated dedicated vector pole."""
    c = np.cosh(eta)
    s = 4*c*c
    g = WIDTH*(s+2)
    return 3*np.pi*(g/c)**2 / ((s-POLE**2)**2 + (g*np.tanh(eta))**2)


def scaled_bessel(order,x):
    """exp(x) K_order(x), including cold arguments beyond SciPy's AMOS limit."""
    x = np.asarray(x, dtype=float)
    out = np.empty_like(x)
    valid = x < 1e5
    out[valid] = kve(order,x[valid])
    large = x[~valid]
    a = 4.*order**2
    u = 1./(8.*large)
    out[~valid] = np.sqrt(np.pi/(2*large))*(1+(a-1)*u*(1+(a-9)*u/2*(1+(a-25)*u/3)))
    return out


def rapidity_integral(theta,r,eta):
    E = 2*np.cosh(eta)
    x = E*math.cosh(r)/theta
    y = 2*np.sinh(eta)*math.sinh(r)/theta
    minus = 2*np.cosh(eta-r)/theta
    plus = 2*np.cosh(eta+r)/theta
    log_minus = -4*np.sinh((eta-r)/2)**2/theta
    log_plus = -4*np.sinh((eta+r)/2)**2/theta
    small = y < .045*np.minimum(x,1.)
    result = np.empty_like(eta)
    if np.any(~small):
        j = ~small
        result[j] = E[j]/(2*y[j])*(scaled_bessel(0,minus[j])*np.exp(log_minus[j])
                                  -scaled_bessel(0,plus[j])*np.exp(log_plus[j]))
    if np.any(small):
        j = small
        xx, yy = x[j], y[j]
        terms = np.zeros_like(xx)
        for i in range(4):
            p = 2*i+1
            cosh_moment = sum(math.comb(p,k)*scaled_bessel(p-2*k,xx)
                              for k in range((p+1)//2))/2**(p-1)
            terms += cosh_moment*yy**(2*i)/math.factorial(2*i+1)
        result[j] = E[j]*np.exp(log_minus[j]-yy)*terms
    return result


@lru_cache(maxsize=None)
def quadrature(order):
    return roots_legendre(order)


def rate(theta,r,order=16,cross_section=sigma):
    """m² Gamma/n_proper in wall coordinate time, averaged over a lab stream."""
    r = abs(r)
    if theta<0:
        raise ValueError('Temperature must be nonnegative')
    if theta==0:
        return 2*math.sinh(r)*cross_section(r)
    edge = 2*math.asinh(math.sqrt(60*theta)/2)
    left, right = max(0.,r-edge), r+edge
    center = math.acosh(POLE/2)
    knots = sorted(set(k for k in (left,right,r,center-1,center-.25,center,center+.25,center+1)
                       if left<=k<=right))
    cuts = [left]
    for a,b in zip(knots,knots[1:]):
        n = max(1,math.ceil((b-a)/min(1.5,3*math.sqrt(theta))))
        cuts.extend(np.linspace(a,b,n+1)[1:])
    cuts = np.array(cuts)
    h = np.diff(cuts)*.5
    x,w = quadrature(order)
    eta = (cuts[:-1,None]+h[:,None]*(1+x)).ravel()
    w = (h[:,None]*w).ravel()
    sh = np.sinh(eta)
    g = np.cosh(eta)*sh**3*cross_section(eta)*rapidity_integral(theta,r,eta)
    return 4*(w @ g)/(theta*(math.cosh(r)*scaled_bessel(2,1/theta))**2)


def cold(theta,u):
    """Nonrelativistic Gaussian interpolation with the exact cold-beam limits."""
    theta,u = np.broadcast_arrays(theta,np.abs(u))
    r = np.arcsinh(u)
    z = np.divide(u,np.sqrt(theta),out=np.zeros_like(u,dtype=float),where=theta>0)
    speed = np.empty_like(z)
    j = z<1e-4
    speed[j] = 2/np.sqrt(np.pi)*(1+z[j]**2/3)
    v = z[~j]
    speed[~j] = (v+.5/v)*erf(v)+np.exp(-v*v)/np.sqrt(np.pi)
    return 2*sigma(r)*np.where(theta==0,u,np.sqrt(theta)*speed)


def table_row(task):
    logT,ys,order = task
    theta = 10**logT
    scale = math.sqrt(theta/(1+theta))
    return [math.log(rate(theta, math.asinh(scale*math.sinh(y)),order)) for y in ys]


def make_table(step,dy,order,jobs):
    ts = np.linspace(-5,6,1+round(11/step))
    ys = np.linspace(0,20,1+round(20/dy))
    tasks = [(t,ys,order) for t in ts]
    if jobs>1:
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            coeff = list(ex.map(table_row,tasks))
    else:
        coeff = list(map(table_row,tasks))
    return ts,ys,np.asarray(coeff)


class Lookup:
    def __init__(self,ts,ys,coeff):
        self.spline = RectBivariateSpline(ts,ys,coeff,s=0)
        self.max_t = ts[-1]
        self.max_y = ys[-1]

    def evaluate(self,theta,u):
        theta,u = np.broadcast_arrays(theta,np.abs(u))
        j = theta>=1e-5
        out = np.empty_like(theta)
        out[~j] = cold(theta[~j],u[~j])
        if np.any(j):
            lt = np.log10(theta[j])
            y = np.arcsinh(u[j]/np.sqrt(theta[j]/(1+theta[j])))
            if np.any(lt>self.max_t) or np.any(y>self.max_y):
                raise ValueError('Event trajectory leaves lookup domain')
            out[j] = np.exp(self.spline.ev(lt,y))
        return out
