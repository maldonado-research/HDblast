"""Leading external-leg soft photons radiated by the recorded equal-mass elastic events."""
import argparse
import json
from pathlib import Path

import numpy as np


def relative_kernel(excess):
    """L(1+x)-1, with L(g)=acosh(g)*g/sqrt(g^2-1)."""
    x=np.asarray(excess,dtype=float)
    if np.any(x<0.) or not np.isfinite(x).all():
        raise ValueError("Invalid relative four-velocity")
    result=np.empty_like(x)
    small=x<1.e-3
    y=x[small]
    result[small]=y*(2./3.+y*(-1./5.+y*(8./105.-y*2./63.)))
    y=x[~small]
    result[~small]=(1.+y)*np.arcsinh(np.sqrt(y*(y+2.)))/np.sqrt(y*(y+2.))-1.
    return result


def relative_derivative(x):
    y=np.asarray(x,dtype=float)
    result=np.empty_like(y)
    small=y<1.e-3
    a=y[small]
    result[small]=2./3.+a*(-2./5.+a*(24./105.-a*8./63.))
    a=y[~small];squared=a*(a+2.)
    result[~small]=((1.+a)-np.arcsinh(np.sqrt(squared))/np.sqrt(squared))/squared
    return result


def charge_kernel(s,mu):
    """Exact leading soft coefficient divided by alpha/pi in the equal-mass case."""
    if np.any(s<4.) or np.any(np.abs(mu)>1.+2.e-12):
        raise ValueError("Events outside physical range")
    mu=np.clip(mu,-1.,1.)
    k=(np.asarray(s)-4.)/4.
    large=2.*k
    diff=k*(1.-mu)
    other=k*(1.+mu)
    delta=relative_kernel(large)-relative_kernel(other)
    # Integrate L' over a short interval. Direct subtraction here loses the
    # forward limit when t is tiny relative to the beam energy.
    near=diff<1.e-4*(1.+other)
    if np.any(near):
        y=other[near]; dx=diff[near]
        delta[near]=dx*(relative_derivative(y)+4*relative_derivative(y+dx/2)+relative_derivative(y+dx))/6.
    return 4.*(relative_kernel(diff)+delta)


def coefficient(s,mu,alpha):
    """Dimensionless B=dN/dlog(omega), including interference and alpha, for any broadcast input."""
    s,mu=np.broadcast_arrays(np.asarray(s,dtype=float),np.asarray(mu,dtype=float))
    result=(alpha/np.pi)*charge_kernel(s.ravel(),mu.ravel()).reshape(s.shape)
    return float(result) if result.ndim==0 else result


def angular_integral(s,mu,nz=64,nphi=96):
    """Independent angular integral of the physical polarization sum with four-pole quadrature."""
    p=np.sqrt((s-4.)/4.)
    if p==0.: return 0.,0.
    gamma=np.hypot(1.,p);beta=p/gamma
    dmin=1./(gamma*(gamma+p))
    dmax=1.+beta
    nodes,gauss=np.polynomial.legendre.leggauss(nz)
    half=.5*np.log(dmax/dmin)
    d=dmin*np.exp(half*(nodes+1.))
    polar_w=.5*half*(d/beta)*gauss
    one_minus_z=(d-dmin)/beta
    z=1.-one_minus_z
    transverse=np.sqrt(np.maximum(0.,one_minus_z*(2.-one_minus_z)))
    phi=(np.arange(nphi)+.37)*2.*np.pi/nphi
    ex=np.array([0.,0.,1.]);ey=np.array([np.sqrt((1.-mu)*(1.+mu)),0.,mu])
    poles=np.array([ex,-ex,ey,-ey])
    total=0.;pulse=np.zeros(3)
    for j,axis in enumerate(poles):
        y=np.array([0.,1.,0.]);x=np.cross(y,axis)
        ray=z[:,None,None]*axis+transverse[:,None,None]*(np.cos(phi)[None,:,None]*x+np.sin(phi)[None,:,None]*y)
        rays=ray.reshape(-1,3)
        den=np.stack([dmin+.5*beta*((rays-pole)**2).sum(axis=1) for pole in poles],axis=1)
        # Partition the angular integral smoothly across four logarithmic polar grids.
        partition=(1./den[:,j]**2)/(1./den**2).sum(axis=1)
        current=2.*beta*(ey/(den[:,2]*den[:,3])[:,None]-ex/(den[:,0]*den[:,1])[:,None])
        power=(np.cross(rays,current)**2).sum(axis=1)
        weight=np.repeat(polar_w/nphi,nphi)*partition*power
        total+=weight.sum()
        pulse+=weight@rays
    return float(total),float(np.linalg.norm(pulse))


def controls():
    tests=[(4.000004,.1),(5.,-.83),(100.,.43),(8300.**2,-.6),(8300.**2,.99)]
    rows=[]
    for s,mu in tests:
        analytic=charge_kernel(np.array([s]),np.array([mu]))[0]
        quad,odd=angular_integral(s,mu)
        rows.append(dict(s=s,cosine=mu,closed_form=float(analytic),angle_integral=quad,
                         relative_error=float(abs(quad/analytic-1)),odd_moment_fraction=float(odd/quad)))
        if abs(quad/analytic-1)>3e-6 or odd/quad>3e-6:
            raise ValueError(("Polarized angular integral disagrees",rows[-1]))
    k=np.array([0.,1.e-7,1.,1.e8]);ss=4.*(1.+k)
    if np.any(charge_kernel(ss,np.ones_like(ss))!=0.) or np.any(charge_kernel(np.full(4,4.),np.linspace(-1.,1.,4))):
        raise ValueError("Missing exact forward or threshold cancellation")
    back=charge_kernel(ss,-np.ones_like(ss))
    if not np.allclose(back,8.*relative_kernel((ss-4.)/2.),rtol=1e-13):
        raise ValueError("Incorrect backward limit")
    low_s=np.array([4.+4.e-7]);low_mu=np.array([.3])
    low=float(charge_kernel(low_s,low_mu)[0]/((16./3.)*((low_s[0]-4.)/4.)*.7)-1.)
    if abs(low)>1e-6: raise ValueError("Incorrect low-speed electric-dipole limit")
    return dict(angular_tests=rows,relative_low_speed_error=low,forward_threshold_backward_passed=True)


def from_shell_legs(first,second):
    """B/(alpha/pi) from all six Minkowski products in the shell frame."""
    legs=(first[:,0],first[:,1],second[:,0],second[:,1])
    sign=(-1,1,1,-1)
    norms=[np.linalg.norm(p,axis=1) for p in legs]
    unit=[np.divide(p,n[:,None],out=np.zeros_like(p),where=n[:,None]>0) for p,n in zip(legs,norms)]
    result=np.zeros(len(first))
    for i in range(4):
        for j in range(i+1,4):
            x=norms[i];y=norms[j]
            rel=(1.+x*x+y*y)/(np.hypot(1.,x)*np.hypot(1.,y)+x*y)+.5*x*y*((unit[i]-unit[j])**2).sum(axis=1)
            result-=2.*sign[i]*sign[j]*relative_kernel(np.maximum(0.,rel-1.))
    return result


def analyze(path):
    with np.load(path) as h:
        s=h["s"];mu=h["mu"];a=h["a"];weight=float(h["weight"]);gas=float(h["gas"])
        qi=float(h["charge_to_mass"]);K=float(h["K"]);mass=float(h["daughter_mass_L"]);S0=float(h["action"])
        alpha=3*qi**2*mass**2/(4*np.pi*K*S0)
        if len(s)!=len(mu) or not np.isfinite(s).all() or not np.isfinite(mu).all() or np.any(a<=0.):
            raise ValueError("Invalid event tape")
        if float(h["channel_fraction"])!=1.:raise ValueError("This formula assumes the same photon field before and after the collision")
        p=h["p_in"];E=h["E_in"].sum(axis=1);p_out=h["p_out"]
        e_out=h["E_out"].sum(axis=1)
        n_total=round(float(h["D"])/weight)
        if abs(n_total*weight/float(h["D"])-1)>1e-12:raise ValueError("Charge weights do not recover total D")
        if len(s) and (np.max(abs(E-e_out)/E)>2e-12 or np.max(np.linalg.norm(p.sum(axis=1)-p_out.sum(axis=1),axis=1)/E)>2e-12):
            raise ValueError("Invalid four-momentum ledger")
        b=coefficient(s,mu,alpha)
        direct=(alpha/np.pi)*from_shell_legs(p,p_out)
        interference_scale=(4.*alpha/np.pi)*relative_kernel((s-4.)/2.)
        if len(b):
            rel_error=float(np.max(abs(direct-b)/np.maximum(abs(b),1e-10*alpha)))
            reflected=coefficient(s,-mu,alpha)
            reflection_error=float(np.max(abs(.5*(b+reflected)-interference_scale)/interference_scale))
            if rel_error>2.e-6 or reflection_error>1.e-12:
                raise ValueError(("Frame or symmetric-angle check failed",rel_error,reflection_error))
        else:
            rel_error=reflection_error=0.
        scale=weight*a*E*.5*(1.-4./s)
        # Width one in xi=omega/omega_kin. Every number is linear in the chosen band width.
        coeff=float(np.sum(scale*b));prop=float(np.sum(b)*np.log(10.)/n_total)
        if np.any(b<0.) or not np.isfinite(b).all():
            raise ValueError("Invalid interference coefficient")
        cm_omega=.5*np.sqrt(s)*(1.-4./s)
        transfer=-.5*(s-4.)*(1.-mu)  # t in daughter mass units
        u_transfer=-.5*(s-4.)*(1.+mu)
        # Positive endpoint test on the Born amplitude's complex s-channel denominator.
        # Fixed coupling and equal masses as in collide.hpp. This is one scale check only.
        pole_dist=np.hypot(s-8300.**2,np.sqrt(1.-4./s)*.30*(s+2.))
        pole_recoil=(s-4.)/pole_dist if len(s) else np.array([])
        count=len(s)
        band=lambda lo,hi:dict(xi_low=lo,xi_high=hi,conformal_emitted=coeff*(hi-lo),
                             fraction_of_endpoint_gas=coeff*(hi-lo)/gas,
                             photons_per_charged_particle=float(np.sum(b)*np.log(hi/lo)/n_total),
                             max_s_shift_over_Born_pole_denominator=float(np.max(pole_recoil,initial=0.)*hi))
        out=dict(tape=path.name,time=float(h["time"]),events=count,charged_particles=n_total,gas=gas,alpha=float(alpha),
                 expected_photons_per_particle_per_decade=prop,
                 energy_over_gas_per_unit_fractional_bandwidth=coeff/gas,
                 soft_bands=[band(lo,10*lo) for lo in [1e-5,1e-4,1e-3]],
                 extrapolations_for_Born_comparison=[band(.01,.1),band(.1,.25)],
                 four_leg_shell_vs_cm_max_relative_error=rel_error,reflection_identity_max_relative_error=reflection_error,
                 integrated_interference_fraction_vs_independent_dipoles=float(coeff/np.sum(scale*interference_scale)-1.) if count else 0.,
                 maximum_relative_error_of_energy_ledger=float(np.max(abs(E-e_out)/E,initial=0.)))
        if count:
            prob=np.array([0.,.1,.5,.9,.99,1.])
            out.update(quantile_probabilities=prob.tolist(),s_quantiles=np.quantile(s,prob).tolist(),
                       root_s_over_pole_mass_quantiles=np.quantile(np.sqrt(s)/8300.,prob).tolist(),
                       cm_omega_kin_quantiles=np.quantile(cm_omega,prob).tolist(),
                       abs_t_quantiles=np.quantile(-transfer,prob).tolist(),abs_u_quantiles=np.quantile(-u_transfer,prob).tolist(),
                       coefficient_quantiles=np.quantile(b,prob).tolist(),cm_boost_quantiles=np.quantile(E/np.sqrt(s),prob).tolist(),
                       soft_energy_fraction_from_top_decile_of_b=float(np.sum((scale*b)[b>=np.quantile(b,.9)])/coeff))
        return out


def main(args):
    control=controls()
    files=sorted(args.artifacts.glob("*.npz"))
    if not files:raise ValueError("No event tapes found")
    results=[]
    for path in files:
        results.append(analyze(path))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(dict(controls=control,runs=results),sort_keys=True,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifacts",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    main(parser.parse_args())
