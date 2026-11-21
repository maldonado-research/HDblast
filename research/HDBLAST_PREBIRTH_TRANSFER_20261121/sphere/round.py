"""Charged test particles in a smooth coexact Maxwell harmonic on the unit three-sphere."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path

import numpy as np
from scipy.special import kv, kve
from scipy.stats import qmc


Q = 1.e7
SEED = 91731
AMPLITUDES = (0., .01, .02, .04, .08, .16)


def directions(c, d, cap, prob):
    """Inverse of acceptance proportional to min(1, |c-d mu|/cap)."""
    mu = 2.*prob - 1.
    use = (d != 0.) & (cap > 0.)
    c, d, cap, u = c[use], d[use], cap[use], prob[use]
    if not c.size:
        return mu
    cuts = np.sort(np.clip(np.stack((-np.ones_like(c), (c-cap)/d, c/d, (c+cap)/d, np.ones_like(c)), axis=1), -1., 1.), axis=1)
    w = np.minimum(1., np.abs(c[:,None]-d[:,None]*cuts)/cap[:,None])
    widths = np.diff(cuts, axis=1)
    areas = (w[:,1:] + w[:,:-1])*widths/2.
    sums = np.column_stack((np.zeros_like(c), np.cumsum(areas,axis=1)))
    target = u*sums[:,-1]
    j = np.minimum(3,(sums[:,1:] <= target[:,None]).sum(axis=1))
    i = np.arange(len(c))
    dx = widths[i,j]
    area = (target-sums[i,j])/dx
    w0, w1 = w[i,j], w[i,j+1]
    coord = 2.*area/(w0+np.sqrt(np.maximum(0.,w0*w0+2*(w1-w0)*area)))
    mu[use] = cuts[i,j]+dx*coord
    return mu


def move(z, u, a, dt):
    """Exact left-invariant great-circle drift at fixed comoving momentum."""
    length=np.linalg.norm(u,axis=1)
    angle=dt*length/np.hypot(a,length)
    vx,vy,vz=u.T
    horizontal=vx+1j*vy
    tangent=np.column_stack((-horizontal*z[:,1].conj()-1j*vz*z[:,0],
                              horizontal*z[:,0].conj()-1j*vz*z[:,1]))
    return np.cos(angle)[:,None]*z + (np.sinc(angle/np.pi)*dt/np.hypot(a,length))[:,None]*tangent


def check_frame():
    """Finite derivatives along three great circles check the chosen curl orientation."""
    sample=qmc.Sobol(3,scramble=True,seed=SEED).random_base2(5)
    z=np.column_stack((np.sqrt(sample[:,0])*np.exp(2j*np.pi*sample[:,1]),
                       np.sqrt(1.-sample[:,0])*np.exp(2j*np.pi*sample[:,2])))
    worst=0.
    for k in (2,5,7,10):
        def harmonic(q):
            v=np.sqrt(k-1.)*q[:,0]**(k-2)
            return np.column_stack((v.real,-v.imag,np.zeros(len(q))))
        deriv=[]
        delta=2.e-5
        for axis in range(3):
            u=np.zeros((len(z),3)); u[:,axis]=1.
            # The frame curves are unit-speed geodesics: a=0 gives speed one.
            deriv.append((harmonic(move(z,u,0.,delta))-harmonic(move(z,u,0.,-delta)))/(2.*delta))
        y=harmonic(z)
        curl=np.column_stack([deriv[(i+1)%3][:,(i+2)%3]-deriv[(i+2)%3][:,(i+1)%3] for i in range(3)])-2.*y
        div=sum(deriv[j][:,j] for j in range(3))
        worst=max(worst,np.max(abs(curl+k*y))/(k*np.sqrt(k-1.)),np.max(abs(div))/np.sqrt(k-1.))
    if worst>1.e-7: raise ValueError('Spherical harmonic curl/divergence check failed')
    return float(worst)


def simulate(history, fields, marker, kappa, alpha, power, steps, seed):
    h=history
    sample=qmc.Sobol(6,scramble=True,seed=seed).random_base2(power)
    # The squared first complex coordinate is uniform under Haar volume on S³.
    born=np.interp(sample[:,0]*marker['D'],h['D'],h['t'])
    s=sample[:,1]
    z0=np.column_stack((np.sqrt(s)*np.exp(2j*np.pi*sample[:,2]),
                        np.sqrt(1.-s)*np.exp(2j*np.pi*sample[:,3])))
    birth=lambda key:np.interp(born,h['t'],h[key])
    pstar=np.sqrt(np.maximum(0.,birth('R')**2-1.))
    mu=directions(birth('c'),birth('d'),birth('target'),sample[:,4])
    pt=pstar*np.sqrt(np.maximum(0.,1.-mu*mu))
    azimuth=2.*np.pi*sample[:,5]
    u0=birth('a')[:,None]*np.column_stack((pt*np.cos(azimuth),pt*np.sin(azimuth),pstar*mu))
    # Follow both members explicitly. Magnetic forces need not send them on opposite paths.
    order=np.argsort(np.r_[born,born])
    born=np.r_[born,born][order]
    z=np.vstack((z0,z0))[order]
    u=np.vstack((u0,-u0))[order]
    sign=np.r_[np.ones(len(z0)),-np.ones(len(z0))][order]
    aT=marker['a']
    t=h['t']; edges=np.linspace(t[0],marker['time'],steps+1)
    mode=int(np.where(fields['kappa']==kappa)[0][0])
    if not np.array_equal(t,fields['t']): raise ValueError('Bulk and birth traces differ')
    E_values=fields['E_cos'][mode]; B_values=fields['B_sin'][mode]
    profile=np.sqrt(kappa-1.)
    initial_pt=np.linalg.norm(u[:,:2],axis=1)
    for left,right in zip(edges[:-1],edges[1:]):
        n=born.searchsorted(right)
        if not n: continue
        start=np.maximum(left,born[:n]); mid=.5*(right+start)
        a=np.interp(mid,t,h['a']); deta=(right-start)/a
        z_mid=move(z[:n],u[:n],a,.5*deta)
        # In the chosen right-handed particle frame, curl Y=-κY and curl e3=-2e3.
        coefficient=profile*z_mid[:,0]**(kappa-2)
        Y=np.column_stack((coefficient.real,-coefficient.imag,np.zeros(n)))
        E=(alpha*np.interp(mid,t,E_values))[:,None]*Y
        B=(alpha*np.interp(mid,t,B_values))[:,None]*Y
        E[:,2]=np.interp(mid,t,h['X'])
        B[:,2]=np.interp(mid,t,h['B'])
        half=.5*Q*sign[:n]*a*a*deta
        p=u[:n]+half[:,None]*E
        rot=(half/np.hypot(a,np.linalg.norm(p,axis=1)))[:,None]*B
        cross=p+np.cross(p,rot)
        p+=np.cross(cross,(2./(1.+np.einsum('ij,ij->i',rot,rot)))[:,None]*rot)
        u[:n]=p+half[:,None]*E
        z[:n]=move(z_mid,u[:n],a,.5*deta)
    sphere_defect=float(np.max(np.abs(np.sum(abs(z)**2,axis=1)-1.)))
    if sphere_defect>2.e-11: raise ValueError('Three-sphere drift lost normalization')
    if not alpha and np.max(abs(np.linalg.norm(u[:,:2],axis=1)-initial_pt))>1.e-10:
        raise ValueError('Axial guide changed the transverse magnitude')
    p=u/aT
    gamma=np.hypot(1.,np.linalg.norm(p,axis=1))
    P=(p.T/gamma)@p
    axial=float(P[2,2]/P.trace())
    pt=np.sort(np.linalg.norm(p[:,:2],axis=1)); theta=marker['theta']
    w=np.hypot(1.,pt)
    F=1.-w*w*kv(2,w/theta)/kv(2,1./theta)
    n=len(pt); rank=np.arange(n)/n
    gap=float(np.maximum(abs(F-rank),abs(F-(rank+1./n))).max())
    thermal=kve(1,1./theta)/kve(2,1./theta)+3.*theta
    # Directions follow the local added harmonic at the observation point on S³.
    # A Lorentz boost along the Hopf axis leaves these two stress components invariant.
    angle=-(kappa-2)*np.angle(z[:,0]); ca,sa=np.cos(angle),np.sin(angle)
    parallel=p[:,0]*ca+p[:,1]*sa; cross=-p[:,0]*sa+p[:,1]*ca
    drift=float(np.interp(marker['time'],h['t'],h['u'])); boost=np.hypot(1.,drift)
    expected_wall=(thermal+theta)*boost-theta/boost
    rest_energy=boost*gamma-sign*drift*p[:,2]
    rest_z=boost*p[:,2]-sign*drift*gamma
    rest_density=np.mean(rest_energy/gamma)
    stress=np.array([np.mean(parallel**2/gamma),np.mean(cross**2/gamma),np.mean(rest_z**2/gamma)])
    target_pressure=theta/boost
    expected_zz=((thermal+theta)*drift*drift+theta)/boost
    matched=dict(signed_positive_stream_u=drift,target_energy_per_particle_wall=expected_wall,
                 mean_energy_over_boosted_thermal=float(gamma.mean()/expected_wall),
                 mean_energy_in_reference_stream_frames_over_thermal=float(np.mean(rest_energy**2/gamma)/rest_density/thermal),
                 perpendicular_pressure_along_pattern_over_target=float(stress[0]/target_pressure),
                 perpendicular_pressure_across_pattern_over_target=float(stress[1]/target_pressure),
                 reference_frame_directional_pressure_fractions=(stress/stress.sum()).tolist(),
                 axial_stress_fraction_target=float(expected_zz/(expected_zz+2.*target_pressure)))
    return dict(kappa=kappa, launch_rms=alpha,boosted_reference=matched,sample_pairs=len(z0),time_steps=steps,scramble=seed,
                transverse_cdf_gap=gap, axial_stress_fraction=axial,mean_gamma_over_thermal=float(gamma.mean()/thermal),
                max_unit_sphere_norm_defect=sphere_defect)


def analytic_bound(h, fields, marker, k, alpha):
    """Post-birth collisionless axial momentum inequality for any initial positions."""
    j=h['t'].searchsorted(marker['time'])
    t=np.r_[h['t'][:j],marker['time']]
    a=np.interp(t,h['t'],h['a']); d=np.interp(t,h['t'],h['D'])
    x=np.interp(t,h['t'],h['X'])
    i=int(np.where(fields['kappa']==k)[0][0])
    b=np.interp(t,fields['t'],fields['B_sin'][i])
    if np.any(x>=0) or np.any(np.diff(a)<0) or np.any(np.diff(d)<0):
        raise ValueError('Reference signs/monotonicity changed')
    dt=np.diff(t)
    # Endpoint rectangles bound each product of linearly interpolated factors, even across zeros of B.
    main=np.sum(dt*a[:-1]*np.minimum(-x[:-1],-x[1:])*d[:-1])
    magnetic=np.sum(dt*a[1:]*np.maximum(abs(b[:-1]),abs(b[1:]))*d[1:])
    factor=Q/(marker['a']*marker['D']); floor=float(np.sqrt(max(h['R'][:j+1].max()**2-1.,0.)))
    cost=factor*max(magnetic*np.sqrt(k-1.),0.)
    v0=factor*main-floor
    thermal=kve(1,1./marker['theta'])/kve(2,1./marker['theta'])+3.*marker['theta']
    return dict(kappa=k, max_harmonic_norm=float(np.sqrt(k-1.)), axial_zero_launch_lower_in_mc=float(v0),
                max_loss_per_unit_launch_in_mc=float(cost),
                unit_mean_thermal_gamma=float(thermal), amplitude_below_which_excess_is_guaranteed=float((v0-thermal)/cost),
                mean_gamma_over_thermal_lower_at_max_scanned_amplitude=float(max(1.,v0-alpha*cost)/thermal))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--history',type=Path,required=True); p.add_argument('--geometry',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True); p.add_argument('--jobs',type=int,default=1)
    args=p.parse_args()
    meta=json.loads((args.history/'markers.json').read_bytes())
    inputs={}
    for label in ('dedicated_half','shared_half','dedicated_native','shared_native'):
        with np.load(args.history/(label+'.npz')) as file: h=dict(file)
        with np.load(args.geometry/(label+'.npz')) as file: fields=dict(file)
        inputs[label]=(h,fields,meta['runs'][label]['markers']['0.01'])
    curl_error=check_frame()
    cases=[]
    for channel in ('dedicated','shared'):
        name=channel+'_half'
        cases.extend((name,k,alpha,10,80,SEED) for k in (2,5,7,10) for alpha in AMPLITUDES if alpha or k==2)
        k2_amp=.08 if channel=='dedicated' else .04
        for k,alpha in ((2,k2_amp),(7,.04 if channel=='dedicated' else .02)):
            cases.extend([(name,k,alpha,12,160,SEED+1),(name,k,alpha,12,320,SEED+1)])
        for alpha in (.02,.04,.08):
            cases.extend([(name,5,alpha,11,160,SEED),(name,5,alpha,11,160,SEED+1),
                          (channel+'_native',5,alpha,10,80,SEED)])
    def one(params):
        label,k,alpha,power,steps,seed=params
        return dict(history=label,**simulate(*inputs[label],k,alpha,power,steps,seed))
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        rows=list(pool.map(one,cases))
    bounds={label:[analytic_bound(*data,k,max(AMPLITUDES)) for k in (2,5,7,10)] for label,data in inputs.items()}
    for row in rows:
        target=next(b for b in bounds[row['history']] if b['kappa']==row['kappa'])
        lower=(target['axial_zero_launch_lower_in_mc']-row['launch_rms']*target['max_loss_per_unit_launch_in_mc'])/target['unit_mean_thermal_gamma']
        if row['mean_gamma_over_thermal'] < lower-.03:
            raise ValueError('Sampled energy inconsistent with population lower bound')
    args.output.write_text(json.dumps(dict(runs=rows,continuum_bounds=bounds,check_normalized_curl_div_error=curl_error),indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
