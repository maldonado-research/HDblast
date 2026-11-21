"""Equal-current and equal-stress angular and radial counterexamples for the isolated vector pole."""
import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import eval_legendre, kve, roots_legendre


def flux_cross_section(s_minus_four):
    """Invariant Møller flux times m² sigma_T, with m=1, beginning at threshold."""
    s = 4. + s_minus_four
    beta = np.sqrt(s_minus_four / s)
    g2 = (.30 * (s + 2.))**2
    return 6. * math.pi * beta * g2 / ((s - 8300.**2)**2 + g2 * beta**2)


def angular_laws(order):
    mu, weights = roots_legendre(order)
    weights = weights * .5
    p4 = eval_legendre(4, mu)
    smooth = dict(isotropic=weights, plus_P4=weights*(1+p4), minus_P4=weights*(1-p4))
    return [(mu, smooth),
            (np.array([-1, 1])/math.sqrt(3.), dict(cones=np.array([.5,.5]))),
            (np.array([-1., 0., 1.]), dict(poles_equator=np.array([1.,4.,1.])/6.))]


def angular_force(E1, E2, r, h, mu, weights, nphi):
    """Average signed impulse times invariant flux, after dividing by h gamma² v."""
    p1, p2 = math.sqrt((E1-1)*(E1+1)), math.sqrt((E2-1)*(E2+1))
    azimuth = (np.arange(nphi)+.5)*math.pi/nphi
    azimuth_part = 2*np.sin(.5*azimuth)**2
    transverse = np.sqrt((1.-mu)*(1.+mu))
    pt2 = (p2*transverse)[:,None]
    et2 = np.hypot(1.,pt2)
    z2 = (p2*mu)[:,None]
    y2 = np.arcsinh(z2/et2)
    gamma = math.cosh(r)
    sh = math.sinh(r)
    if sh == 0:
        raise ValueError("The rate as evaluated requires nonzero stream drift")
    force = np.empty((len(mu),len(mu)))
    for i in range(len(mu)):
        pt1 = p1*transverse[i]
        et1 = math.hypot(1.,pt1)
        z1 = p1*mu[i]
        dy = math.asinh(z1/et1)-y2
        # dot(p1,p2) - 1 = E_T1 E_T2 (cosh(Delta y) - 1)
        # + E_T1 E_T2 - p_T1 p_T2 - 1 + p_T1 p_T2 (1-cos(phi)).
        transverse_part = ((pt1-pt2)**2/(et1*et2+pt1*pt2+1.) +
                           pt1*pt2*azimuth_part)
        gplus = flux_cross_section(2.*(transverse_part+2.*et1*et2*np.sinh(.5*dy+r)**2))
        gminus = flux_cross_section(2.*(transverse_part+2.*et1*et2*np.sinh(.5*dy-r)**2))
        # Simultaneous reversal of both rest directions pairs equal-weight cells.
        # The expression eliminates cancellation of opposite O(1) signed forces.
        force[i] = np.mean((E1+E2)*(gplus+gminus)/gamma +
                           (z1-z2)*(gplus-gminus)/sh, axis=1)/(4.*E1*E2*h)
    return {name:float(w @ force @ w) for name,w in weights.items()}


def finite_design(E, r, h, directions):
    """Direct two-beam finite sum in the wall frame, with no angular quadrature."""
    directions=np.asarray(directions,float)
    if not (np.allclose(directions.mean(axis=0),0.) and
            np.allclose(directions.T @ directions/len(directions),np.eye(3)/3.)):
        raise ValueError("Finite angular design has incorrect first or second moments")
    spatial=math.sqrt((E-1.)*(E+1.))*directions
    gamma,sh=math.cosh(r),math.sinh(r)
    lab1=spatial.copy()
    lab2=spatial.copy()
    z=spatial[:,2]
    lab1[:,2]=gamma*z+sh*E
    lab2[:,2]=gamma*z-sh*E
    e1=gamma*E+sh*z
    e2=gamma*E-sh*z
    # The two equal-mass four-vectors have s-4 = |p1-p2|² - (E1-E2)².
    delta=lab1[:,None,:]-lab2[None,:,:]
    sm4=np.maximum(0.,np.sum(delta*delta,axis=2)-(e1[:,None]-e2[None,:])**2)
    transfers=delta * flux_cross_section(sm4)[:,:,None] / (2*E*E)
    vector=transfers.mean(axis=(0,1))/(h*gamma*sh)
    if np.max(abs(vector[:2]))>1e-6*abs(vector[2]):
        raise ArithmeticError("Unexpected net transverse transfer")
    return float(vector[2])


def zero_drift_isotropic(E1, E2, h):
    """One-dimensional check from differentiating the invariant impulse at zero drift."""
    p1,p2=math.sqrt((E1-1.)*(E1+1.)),math.sqrt((E2-1.)*(E2+1.))
    denom0=E1*E2+p1*p2+1.
    base=2.*(p1-p2)**2/denom0
    def integrand(q):
        sm4=base+2.*p1*p2*(1.-q)
        if sm4==0.:
            return 0.
        s=4.+sm4
        b2=sm4/s
        g2=(.30*(s+2.))**2
        den=(s-8300.**2)**2+g2*b2
        d_den=2.*(s-8300.**2)+g2*(2.*b2/(s+2.)+4./s**2)
        g=flux_cross_section(sm4)
        d_g=g*(2./(s+2.)+2./(s*sm4)-d_den/den)
        moment=(E2*p1**2+E1*p2**2-(E1+E2)*p1*p2*q)/3.
        return ((E1+E2)*g+4.*moment*d_g)/(4.*E1*E2*h)
    points=[]
    if p1*p2>0:
        pole_q=1.-(8300.**2-4.-base)/(2.*p1*p2)
        if -1.<pole_q<1.:
            points.append(pole_q)
    return quad(integrand, -1., 1., points=points, epsabs=1.e-24, epsrel=1.e-10)[0]


def radial_witness(mean_energy, mean_inverse, fraction):
    """Two-atom solution, with lower energy interpolating from 1 toward 1/<1/E>."""
    lower = 1. + fraction*(1./mean_inverse-1.)
    upper = (mean_energy-lower)/(1.-mean_inverse*lower)
    weight = (upper-mean_energy)/(upper-lower)
    energies = np.array([lower,upper])
    masses = np.array([weight,1.-weight])
    if not (0 < weight < 1 and lower>=1 and upper>lower):
        raise ValueError("Invalid positive radial distribution")
    if not np.allclose([masses@energies,masses@(1./energies)],
                       [mean_energy,mean_inverse],rtol=2e-11,atol=0.):
        raise ArithmeticError("Radial moment check failed")
    return energies,masses


def evaluate(theta, u, order):
    r = math.asinh(abs(u))
    etm = 1.5*theta+math.hypot(1.,1.5*theta)
    htm = etm+theta
    tm = {}
    laws = angular_laws(order)
    for mu,weights in laws:
        for w in weights.values():
            if not np.allclose([np.sum(w),w@mu,w@(mu**2)], [1.,0.,1./3.],atol=2e-14):
                raise ArithmeticError("Angular moment check failed")
        tm.update(angular_force(etm,etm,r,htm,mu,weights,max(96,2*order)))
    tm['octahedron']=finite_design(etm,r,htm,np.r_[np.eye(3),-np.eye(3)])
    cube=np.array([[x,y,z] for x in (-1,1) for y in (-1,1) for z in (-1,1)])/math.sqrt(3.)
    tm['cube']=finite_design(etm,r,htm,cube)
    # K1/K2 is <1/E>; the recurrence gives <E> = 3 Theta + K1/K2.
    B = float(kve(1,1./theta)/kve(2,1./theta))
    A = 3.*theta+B
    mu = laws[0][0]
    weights = dict(isotropic=laws[0][1]['isotropic'])
    synge = {}
    for name,fraction in [('lower_at_rest',0.),('lower_near_harmonic_mean',.9)]:
        energy,w = radial_witness(A,B,fraction)
        total=0.
        zero_drift=0.
        for i in range(2):
            for j in range(i+1):
                coefficient=1. if i==j else 2.
                integral=angular_force(energy[i],energy[j],r,A+theta,mu,weights,max(96,2*order))
                total+=coefficient*w[i]*w[j]*integral['isotropic']
                if abs(u)<.001:
                    zero_drift+=coefficient*w[i]*w[j]*zero_drift_isotropic(energy[i],energy[j],A+theta)
        synge[name]=dict(energies=energy.tolist(),weights=w.tolist(),C=total)
        if abs(u)<.001:
            synge[name]['zero_drift_check']=zero_drift
    return dict(TM=dict(mean_energy=etm,mean_inverse_energy=1./etm,C=tm),
                Synge=dict(mean_energy=A,mean_inverse_energy=B,witnesses=synge))


def continuous_radial_witness(A,B,position):
    """Two compact, finite-width energy bands with exactly adjusted mean inverse energy."""
    lower=1.+position*(1./B-1.)
    def band(center):
        width=.04*(center-1.)
        ratio=(width/center)**2
        inverse=sum(3.*ratio**k/((2*k+1)*(2*k+3)) for k in range(8))/center
        return width,inverse
    def mixture(upper):
        weight=(upper-A)/(upper-lower)
        return weight*band(lower)[1]+(1.-weight)*band(upper)[1]
    high=2.*A
    while mixture(high)<B:
        high*=2.
    upper=brentq(lambda z:mixture(z)-B,A,high,xtol=1e-13)
    w=(upper-A)/(upper-lower)
    return [(lower,band(lower)[0],w),(upper,band(upper)[0],1.-w)]


def smooth_rates(theta,u,number_density,order,radial_order):
    """Finite-phase-density Synge witnesses for the startup example."""
    B=float(kve(1,1./theta)/kve(2,1./theta))
    A=3.*theta+B
    r=math.asinh(abs(u))
    mu,weights=angular_laws(order)[0]
    x,wx=roots_legendre(radial_order)
    weights={'isotropic':weights['isotropic']}
    result={}
    for name,position in [('low_band',.02),('high_band',.9)]:
        bands=continuous_radial_witness(A,B,position)
        points=np.concatenate([energy+width*x for energy,width,w in bands])
        prob=np.concatenate([.75*w*wx*(1.-x*x) for energy,width,w in bands])
        occ_bound=0.
        for energy,width,w in bands:
            emin=energy-width
            occ_bound+=math.pi**2*number_density*.75*w/(width*emin*math.sqrt(emin**2-1.))
        if not np.allclose([sum(prob),prob@points,prob@(1./points)],[1.,A,B],rtol=3e-9,atol=0.):
            raise ArithmeticError("Finite-band moments failed")
        total=0.
        for i in range(len(points)):
            for j in range(i+1):
                c=angular_force(points[i],points[j],r,A+theta,mu,weights,2*order)['isotropic']
                total+=(1. if i==j else 2.)*prob[i]*prob[j]*c
        result[name]=dict(C=total,bands=bands,per_spin_f_upper=occ_bound)
    return result


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--history',type=Path,required=True)
    cli.add_argument('--output',type=Path,required=True)
    args=cli.parse_args()
    meta=json.loads((args.history/'metadata.json').read_bytes())
    cols={v:i for i,v in enumerate(meta['columns'])}
    entry=meta['runs']['dedicated_native']
    rows=np.memmap(args.history/'dedicated_native.bin',dtype='=f8',mode='r').reshape(-1,len(cols))
    if len(rows)!=entry['rows']:
        raise ValueError("Trajectory length differs from metadata")
    time=rows[:,cols['t']]
    if np.any(np.diff(time)<=0):
        raise ValueError("Time knots must increase")
    out={}
    for key,t in [('startup',6.255),('during_release',6.30),('first_hot',entry['first_hot_time'])]:
        theta,u=[float(np.interp(t,time,rows[:,cols[name]])) for name in ('theta','u')]
        fine=evaluate(theta,u,192)
        coarse=evaluate(theta,u,128)
        tm={k:abs(v/coarse['TM']['C'][k]-1) for k,v in fine['TM']['C'].items()}
        sy={k:abs(v['C']/coarse['Synge']['witnesses'][k]['C']-1)
            for k,v in fine['Synge']['witnesses'].items()}
        out[key]=dict(t=t,theta=theta,u=u,**fine,
                      quadrature_comparison=dict(TM=tm,Synge=sy))
    t=out['startup']['t']
    a,D=[float(np.interp(t,time,rows[:,cols[x]])) for x in ('a','D')]
    sample=out['startup']
    density=meta['action']*D/(6.*a**3*meta['daughter_mass_L']**4*math.hypot(1.,sample['u']))
    high=smooth_rates(sample['theta'],sample['u'],density,64,7)
    low=smooth_rates(sample['theta'],sample['u'],density,40,5)
    if max(abs(v['C']/low[k]['C']-1.) for k,v in high.items())>1e-3:
        raise ArithmeticError("Finite-band rate quadrature under-resolved")
    sample['finite_band_Synge']=dict(n_over_m_cubed=density,witnesses=high,
        quadrature_relative_change={k:abs(v['C']/low[k]['C']-1.) for k,v in high.items()})
    for case in out.values():
        forces=case['TM']['C']
        radial=case['Synge']['witnesses']
        case['ratios']=dict(polar_equatorial_over_cones=forces['poles_equator']/forces['cones'],
            smooth_positive_over_negative_P4=forces['plus_P4']/forces['minus_P4'],
            isotropic_Synge_energy_pair=radial['lower_near_harmonic_mean']['C']/radial['lower_at_rest']['C'])
    args.output.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
