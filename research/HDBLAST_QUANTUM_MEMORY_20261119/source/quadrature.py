"""Positive high-band oscillator quadrature with exact static/kinetic moments.

Use x=omega^-2 and dmu=g²(domega)/omega². Gaussian polynomial quadrature
in x matches C=integral dmu and D=integral x dmu for every N>=1.
Construction uses a positive reference measure and reorthogonalized Lanczos
on the scaled interval, avoiding ill-conditioned Hankel moment matrices.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh_tridiagonal
from scipy.special import roots_legendre

LOW,HIGH = 60.,120.
ORDERS = (1,2,4,8,16,32,64)


def high_moments(q,low=LOW,high=HIGH):
    c = (high**2-low**2-q*q*np.log(high**2/low**2))/(8*np.pi**2)
    d = (np.log(high**2/low**2)+q*q*(1/high**2-1/low**2))/(8*np.pi**2)
    return c,d


@lru_cache(maxsize=None)
def gauss_x(q,order,base_order=512,low=LOW,high=HIGH):
    """Return omega,g²,mu,x for a positive Gaussian x rule."""
    if not(0<q<low<high) or not(1<=order<base_order):
        raise ValueError('Require positive q below band, and order < base order')
    c,d = high_moments(q,low,high)
    left,right = high**-2,low**-2
    midpoint,half = (left+right)/2,(right-left)/2
    if order==1:
        x = np.array([d/c]);mu = np.array([c])
    else:
        y,w = roots_legendre(base_order)
        xbase = midpoint+half*y
        mass = half*w*(1-q*q*xbase)/(8*np.pi**2*xbase**2)
        mass *= c/np.sum(mass)
        vector = np.sqrt(mass/c)
        basis = np.empty((order,base_order))
        alpha = np.empty(order)
        beta = np.empty(order-1)
        previous = np.zeros_like(vector)
        previous_beta = 0.
        for k in range(order):
            basis[k] = vector
            alpha[k] = np.dot(vector,y*vector)
            remainder = (y-alpha[k])*vector-previous_beta*previous
            # Two-pass modified Gram-Schmidt protects the high-order rule.
            for _ in range(2):
                remainder -= (basis[:k+1]@remainder)@basis[:k+1]
            if k+1<order:
                beta[k] = np.linalg.norm(remainder)
                if not beta[k]>0:
                    raise ArithmeticError('Lanczos measure lost rank')
                previous,vector = vector,remainder/beta[k]
                previous_beta = beta[k]
        node,eigenvectors = eigh_tridiagonal(alpha,beta,check_finite=False)
        x = midpoint+half*node
        mu = c*eigenvectors[0]**2
    omega,g2 = 1/np.sqrt(x),mu/x
    # Sort by omega for straightforward comparison with direct frequency rules.
    index = np.argsort(omega)
    omega,g2,mu,x = (v[index] for v in (omega,g2,mu,x))
    assert np.all((omega>low)&(omega<high)) and np.all(g2>0)
    np.testing.assert_allclose([np.sum(mu),np.sum(mu*x)],[c,d],rtol=3e-13,atol=1e-15)
    return omega,g2,mu,x


@lru_cache(maxsize=None)
def direct_high(q,order,low=LOW,high=HIGH):
    node,weight = leggauss(order)
    omega = low+(node+1)*(high-low)/2
    dw = weight*(high-low)/2
    g2 = omega*(omega*omega-q*q)*dw/(4*np.pi**2)
    return omega,g2


def make_bath(q,high_order,low_order=64):
    """Original low GL64 plus compressed high; N32/64 are direct references."""
    q = np.asarray(q)
    node,weight = leggauss(low_order)
    omega_low = q[:,None]+(LOW-q[:,None])*(node+1)/2
    dw = (LOW-q[:,None])*weight/2
    g2_low = omega_low*(omega_low**2-q[:,None]**2)*dw/(4*np.pi**2)
    if high_order==0:
        return omega_low,g2_low,'original_low_GL64'
    if high_order in (32,64):
        high = [direct_high(float(k),high_order) for k in q]
        method = 'direct_frequency_GaussLegendre'
    else:
        high = [gauss_x(float(k),high_order)[:2] for k in q]
        method = 'positive_Gaussian_in_inverse_frequency_squared'
    omega_high = np.array([v[0] for v in high])
    g2_high = np.array([v[1] for v in high])
    return np.concatenate([omega_low,omega_high],axis=1),np.concatenate([g2_low,g2_high],axis=1),method


def calculate(q,weight,output):
    output = Path(output)
    output.mkdir(parents=True,exist_ok=True)
    arrays = dict(q=np.asarray(q),weight=np.asarray(weight),orders=np.array(ORDERS))
    report = dict(schema=1,band=[LOW,HIGH],variable='x=omega^-2',
        measure='dmu=(1-q^2*x)/(8*pi^2*x^2)dx = g²(domega)/omega²',
        reconstruction='omega=x^-1/2,g²=mu/x',
        source_rule='N1..16 x-Gauss; N32/64 independent direct frequency references',
        mode_count=len(q),orders={})
    for order in ORDERS:
        values = [gauss_x(float(k),order) for k in q]
        refined = [gauss_x(float(k),order,base_order=1024) for k in q]
        for index,name in enumerate(('omega','g2','mu','x')):
            arrays[f'{name}_N{order}'] = np.array([v[index] for v in values])
        c,d = np.array([high_moments(k) for k in q]).T
        mu,x = arrays[f'mu_N{order}'],arrays[f'x_N{order}']
        moment_error = np.max(abs(np.array([np.sum(mu,axis=1)/c-1,np.sum(mu*x,axis=1)/d-1])))
        # Gaussian degree 2N-1 exactness, checked in the scaled variable so
        # high inverse-frequency powers do not underflow or destroy conditioning.
        base_y,base_w = roots_legendre(1024)
        midpoint,half = (HIGH**-2+LOW**-2)/2,(LOW**-2-HIGH**-2)/2
        base_x = midpoint+half*base_y
        polynomial_error = 0.
        for qi,k in enumerate(q):
            base_mass = half*base_w*(1-k*k*base_x)/(8*np.pi**2*base_x**2)
            base_mass *= c[qi]/np.sum(base_mass)
            gauss_y = (x[qi]-midpoint)/half
            for degree in range(2*order):
                difference = abs(np.sum(mu[qi]*gauss_y**degree)-np.sum(base_mass*base_y**degree))/c[qi]
                polynomial_error = max(polynomial_error,float(difference))
        # Compare the actual finite response and causal kernels, not node
        # eigenvector signs, when checking base-rule stability.
        probes = np.array([0.,.25,.5,1.,2.])*HIGH
        error = 0.
        for values0,values1 in zip(values,refined):
            o0,g0 = values0[:2];o1,g1 = values1[:2]
            r0 = np.sum(g0/(o0*o0+probes[:,None]**2),axis=1)
            r1 = np.sum(g1/(o1*o1+probes[:,None]**2),axis=1)
            error = max(error,float(np.max(abs(r0-r1)/r1)))
        report['orders'][str(order)] = dict(min_strength=float(np.min(arrays[f'g2_N{order}'])),
            max_relative_static_kinetic_moment_error=float(moment_error),
            max_abs_normalized_scaled_polynomial_moment_error=polynomial_error,
            base512_vs1024_max_relative_response_change=error)
        assert moment_error<3e-13 and polynomial_error<3e-12 and error<3e-12
        if order in (32,64):
            direct = [direct_high(float(k),order) for k in q]
            arrays[f'direct_omega_N{order}'] = np.array([v[0] for v in direct])
            arrays[f'direct_g2_N{order}'] = np.array([v[1] for v in direct])
    np.savez_compressed(output/'RULES.npz',**arrays)
    (output/'RULES.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    return report


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--modes',type=Path,required=True,help='JSON with native q and weight arrays')
    parser.add_argument('--output',type=Path,required=True)
    args = parser.parse_args()
    modes = json.loads(args.modes.read_bytes())
    calculate(modes['q'],modes['weight'],args.output)
