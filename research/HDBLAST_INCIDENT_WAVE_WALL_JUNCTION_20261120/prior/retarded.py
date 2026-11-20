"""Transverse photon continuum on a single static tuned half-leaf."""
import argparse
import json
import importlib.util
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import gamma, hankel1

K = 5 / 9
S = (1 + 4/K)/2
NU = S + .5


def primitive(y):
    """Minus three times the unshifted warp logarithm."""
    return y + (2/3)*math.log(math.cosh(y)) + math.tanh(y)**2/6


def geometry(y, primitive_b, weight):
    t = math.tanh(y)
    warp = -(primitive(y) - primitive_b)/3
    warp_y = -(1 + t - t**3/3)/3
    return warp, (1.5*warp_y if weight == "selected" else -2*t + warp_y/2)


def norm(phi):
    return (1 + phi)**2 * (2 - phi)/3


def threshold(phi):
    yb = -math.atanh(phi)
    a0 = (primitive(yb) + (2/3)*math.log(2) - 1/6)/3
    j_inf = 4 * math.exp(-2*a0/K) / K**S
    hankel_coefficient = 2**(NU-.5) * gamma(NU) / math.sqrt(math.pi)
    return 2*j_inf**2 / (math.pi*hankel_coefficient**2*norm(phi)**2)


def continuum(phi, masses, *, outer=8., rtol=1e-10, forward=False, weight="wall"):
    """Return density per unit positive mass, full Green function and flux diagnostic.

    Local brane warp factor is one, all masses in inverse wall length. Each
    positive supplied mass is integrated along the same half-leaf. Use the
    public density for a brane source coupled to the one-copy boundary field.
    """
    m = np.atleast_1d(np.asarray(masses, dtype=float))
    if not -1 < phi < 1 or np.any(m <= 0) or not np.all(np.isfinite(m)):
        raise ValueError('Use interior shell values and positive finite masses.')
    yb = -math.atanh(phi)
    if outer <= yb:
        raise ValueError('Matching surface must be beyond the shell.')
    pb = primitive(yb)
    if weight not in ("wall","selected"):
        raise ValueError("Unknown Maxwell weight")
    s = 1.5 if weight == "selected" else S
    nu = s + .5
    jb = 1. if weight == "selected" else 1 - phi**2
    a0 = (pb + (2/3)*math.log(2) - 1/6)/3
    radius = math.exp(K*outer - a0)/K
    arg = m*radius
    ay, hy = geometry(outer, pb, weight)
    drdz = K * radius * math.exp(ay)
    f = np.sqrt(math.pi*arg/2) * hankel1(nu,arg)
    # d = (d/dz - j_z/j)u. Use an adjacent Hankel order to avoid
    # subtracting two leading powers when m is small.
    d = (drdz*m*np.sqrt(math.pi*arg/2)*hankel1(nu-1,arg)
         -(s*K+hy)*math.exp(ay)*f)
    # Wronskian of this pair is m dr/dz; divide both by sqrt(dr/dz)
    # so the transmitted current is m at each finite matching location.
    f /= math.sqrt(drdz)
    d /= math.sqrt(drdz)
    def rhs(y,state):
        warp, jy = geometry(y,pb,weight)
        dzdy = math.exp(-warp)
        u,v = state.reshape(2,-1)
        return np.r_[jy*u + dzdy*v, -jy*v - dzdy*m*m*u]
    if forward:
        # A Neumann regular solution with psi_b=1 has u_b=j_b, d_b=0.
        # Its real Wronskian against the outgoing solution gives the Jost
        # denominator without any subtraction at the brane.
        init = np.r_[np.full(m.shape,jb),np.zeros_like(m)]
        sol = solve_ivp(rhs,(yb,outer),init,rtol=rtol,atol=rtol*1e-4,method='DOP853')
        if not sol.success: raise ArithmeticError(sol.message)
        regular, regular_d = sol.y[:,-1].reshape(2,-1)
        wronskian = regular*d - regular_d*f
        return 2*m*m/(math.pi*abs(wronskian)**2)
    sol = solve_ivp(rhs,(outer,yb),np.r_[f,d],rtol=rtol,atol=rtol*1e-4,method='DOP853')
    if not sol.success: raise ArithmeticError(sol.message)
    u,v = sol.y[:,-1].reshape(2,-1)
    green = -u / (jb**2*v)
    # The square-amplitude formula avoids cancellation in Im G near m=0.
    density = 2*m*m/(math.pi*jb**2*abs(v)**2)
    current = np.imag(u.conj()*v)/m
    return density, green, current


def spectral_density(phi_b, m_grid, **solver):
    """Continuum measure only, with its continuous extension zero at m=0."""
    masses = np.atleast_1d(np.asarray(m_grid,dtype=float))
    if np.any(masses<0) or not np.all(np.isfinite(masses)):
        raise ValueError('Expected finite nonnegative masses')
    result = np.zeros_like(masses)
    active = masses>0
    if active.any():
        result[active] = continuum(phi_b,masses[active],**solver)[0]
    return result


def spectral_integrals(phi, order):
    """Positive Gauss-Legendre sums on dyadic sub-bands for continuum observables."""
    nodes, weights = np.polynomial.legendre.leggauss(order)
    edges = np.array([0.,.125,.25,.5,1.,2.,4.,8.,16.,32.])
    all_m, all_weights = [], []
    for lo,hi in zip(edges[:-1],edges[1:]):
        m = (hi+lo)/2 + (hi-lo)*nodes/2
        # At high mass the oscillation period shrinks. The inverse-square
        # matching is already accurate at a shallower surface in that band.
        rho = spectral_density(phi,m,outer=8. if hi<=8 else 6.)
        all_m.extend(m)
        all_weights.extend((hi-lo)*weights*rho/2)
    return np.array(all_m),np.array(all_weights)


def cross_checks(phi, m, weights, other_m, other_weights, leaf):
    """Compare independent spectral and Euclidean boundary responses."""
    n = norm(phi)
    high = 32.
    bare_density = 2/(math.pi*(1-phi**2)**2)
    yb = -math.atanh(phi)
    t = math.tanh(yb)
    ay = -(1+t-t**3/3)/3
    hy = -2*t+ay/2
    hyy = -2*(1-t*t)-(1-t*t)**2/6
    correction = (hyy + ay*hy - hy*hy)/2
    spacelike = {}
    moment = leaf.low_moment()[0]/n**2
    for q in (0.,.3,1.,3.):
        direct = (1/leaf.response(q)-1/(n*q*q)) if q else moment
        def finite_tail(cutoff):
            # Boundary high-frequency symbol through second inverse power.
            t0 = math.atan(q/cutoff)/q if q else 1/cutoff
            t2 = (1/cutoff-t0)/(q*q) if q else 1/(3*cutoff**3)
            return bare_density*(t0+correction*t2)
        spectral = float(weights @ (1/(m*m+q*q)) + finite_tail(high))
        coarse = float(other_weights @ (1/(other_m**2+q*q)) + finite_tail(high))
        inner = m<16
        trunc = float(weights[inner] @ (1/(m[inner]**2+q*q)) + finite_tail(16.))
        err = spectral/direct-1
        if abs(err)>4e-6:
            raise ArithmeticError(f'Stieltjes continuum disagrees with elliptic response: {err}')
        spacelike[f'{q:g}'] = dict(euclidean_continuum=direct,timelike_transform=spectral,
                                    relative_difference=err,quadrature_relative=coarse/spectral-1,
                                    upper_matching_relative=trunc/spectral-1)
    pulses={}
    for sigma in (.25,.5,1.,2.,4.,8.,16.):
        ratio = n*float(weights @ np.exp(-(sigma*m)**2))
        coarse = n*float(other_weights @ np.exp(-(sigma*other_m)**2))
        bound = n*moment*high**2*math.exp(-(high*sigma)**2)
        asym = n*threshold(phi)*gamma((2*S-1)/2)/(2*sigma**(2*S-1))
        pulses[f'{sigma:g}'] = dict(continuum_over_pole=ratio,continuum_fraction=ratio/(1+ratio),
                                     quadrature_relative=coarse/ratio-1,omitted_band_upper=bound,
                                     ratio_to_long_pulse=float(ratio/asym))
    return dict(continuum_stieltjes=spacelike,gaussian_current=pulses,
                high_mass_density=bare_density,high_mass_relative_inverse_square=correction)

def selected_weight_result(phi, masses):
    """Test the null-uplift-selected weight against its independent threshold."""
    from scipy.integrate import quad
    b = -math.atanh(phi)
    pb = primitive(b)
    # Remainder beyond 48 is analytic to roundoff at this scale.
    a0 = (pb+(2/3)*math.log(2)-1/6)/3
    end=48.
    n=quad(lambda y:math.exp(-2*(primitive(y)-pb)/3),b,end,
           epsabs=2e-14,epsrel=2e-13)[0]+math.exp(2*(a0-K*end))/(2*K)
    c=1/(4*K**3*n*n)
    density, green, flux=continuum(phi,masses,weight="selected")
    fine=spectral_density(phi,masses,weight="selected",outer=10.,rtol=3e-11)
    indep=continuum(phi,masses,weight="selected",outer=9.,rtol=6e-11,forward=True)
    ratios=density/(masses*c)
    if (abs(ratios[0]-1)>4e-6 or max(abs(fine/density-1))>2e-6 or
            max(abs(indep/density-1))>2e-6):
        raise ArithmeticError("Alternate-weight threshold or integration failed")
    return dict(weight="G=exp(2A); warp normalized at shell",half_norm=n,
                threshold_coefficient=c,threshold_power=1,
                samples=[dict(mass=float(m),density=float(v),threshold_ratio=float(r))
                         for m,v,r in zip(masses,density,ratios)],
                independent_max_relative=float(max(abs(indep/density-1))),
                refinement_max_relative=float(max(abs(fine/density-1))),
                current_maximum_error_above_half=float(max(abs(flux[masses>=.5]-1))))


def calculate():
    other = Path(__file__).resolve().parents[1]/'euclidean/euclidean.py'
    spec = importlib.util.spec_from_file_location('elliptic', other)
    elliptic = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(elliptic)
    m = np.array([.005,.01,.02,.04,.08,.125,.25,.5,.75,1.,1.5,2.,3.,4.,6.,8.])
    results = {}
    alternate = {}
    for phi in (0.,.5,-.5):
        density, green, flux = continuum(phi,m)
        finer, _, _ = continuum(phi,m,outer=10.,rtol=3e-11)
        independent = continuum(phi,m,outer=9.,rtol=6e-11,forward=True)
        alternate[f'{phi:g}']=selected_weight_result(phi,np.r_[.0003,.001,m])
        cases = []
        for i,mass in enumerate(m):
            asym = threshold(phi)*mass**(2*S-2)
            cases.append(dict(mass=mass,density=float(density[i]),green_real=float(green[i].real),
                              green_imaginary_from_flux=float(math.pi*density[i]/(2*mass)),
                              threshold_ratio=float(density[i]/asym),
                              refinement_relative=float(finer[i]/density[i]-1),
                              independent_relative=float(independent[i]/density[i]-1)))
        keep = m >= .5
        mg,w = spectral_integrals(phi,24)
        mg2,w2 = spectral_integrals(phi,16)
        observables=cross_checks(phi,mg,w,mg2,w2,elliptic.Leaf(phi))
        if (max(abs(row['refinement_relative']) for row in cases)>2e-6
                or max(abs(row['independent_relative']) for row in cases)>2e-6
                or max(abs(flux[keep]-1))>5e-7):
            raise ArithmeticError('Outgoing normalization or convergence failed.')
        results[f'{phi:g}'] = dict(**observables,half_norm=norm(phi), pole_residue=1/norm(phi),
              threshold_power=2*S-2,threshold_coefficient=threshold(phi),
              current_maximum_error_above_half=float(max(abs(flux[keep]-1))),
              spectral_current_maximum_relative_above_half=float(max(abs(
                  2*m[keep]*green.imag[keep]/(math.pi*density[keep])-1))),
              samples=cases)
    return dict(mass_unit='inverse proper wall width; warp equals one on the shell',
                action='one retained copy; g5 squared set to one',branches=results,alternate_weight=alternate)


if __name__ == '__main__':
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument('--output',type=Path,required=True)
    options=arg.parse_args()
    options.output.write_text(json.dumps(calculate(),indent=2,allow_nan=False)+'\n')
