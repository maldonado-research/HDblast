"""Continuous-time comparison bound for Poisson birth fluctuations in the transverse model."""
import argparse
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.linalg import expm
from scipy.integrate import solve_ivp

ROOT=Path(__file__).resolve().parent
K=2.54


def positive_exponential(beta,zeta,momentum_factor=1.):
    """Scaled positive matrix with constant unit source, solved by independent methods."""
    matrix=np.array([[0.,beta,1.,1.],[momentum_factor,0.,0.,0.],[0.,zeta,0.,0.],[0.,0.,0.,0.]])
    initial=np.array([0.,0.,0.,1.])
    exact=expm(matrix)@initial
    term=initial.copy()
    series=term.copy()
    norm=float(np.linalg.norm(matrix,np.inf))
    order=max(100,math.ceil(4*norm)+100)
    for n in range(1,order+1):
        term=matrix@term/n
        series+=term
    log_tail=(order+1)*math.log(norm)-math.lgamma(order+2)-math.log1p(-norm/(order+2))
    tail=math.exp(log_tail) if log_tail>-740 else 1.e-320
    ode=solve_ivp(lambda t,x:matrix@x,(0.,1.),initial,method="DOP853",rtol=2.e-12,atol=1.e-14).y[:,-1]
    err=float(max(np.max(abs(ode-exact)/(1+abs(exact))),np.max(abs(series-exact)/(1+abs(exact)))))
    if err>2e-10 or tail>1e-11:raise ValueError("Comparison system integration failed")
    # Round upward well beyond the truncation and linear-algebra differences.
    return (series+tail)*(1+1e-10),err,tail


@lru_cache(maxsize=8)
def spatial_wavevector_pairs(radius):
    """Count nonzero cubic integer lattice points inside a ball, modulo k -> -k."""
    total=0
    for i in range(radius+1):
        yz=radius*radius-i*i
        plane=0
        for j in range(math.isqrt(yz)+1):
            plane+=(1 if j==0 else 2)*(2*math.isqrt(yz-j*j)+1)
        total+=(1 if i==0 else 2)*plane
    return (total-1)//2


def bound(table, marker, action, mass, Q, kmax, fundamental, margin, general=False):
    mask=table['t']<marker['time']
    def take(key):return np.append(table[key][mask],marker['time'] if key=='t' else marker[key])
    t,a,X,D=[take(key) for key in ('t','a','X','D')]
    if D[0]!=0. or np.any(np.diff(D)<0.) or np.any(np.diff(a)<0.) or np.any(X>=0.):
        raise ValueError("Bound assumes empty start, expanding wall and negative axial field")
    T=float(t[-1]-t[0]);a0=float(a[0]);b=math.sqrt(1.05**2-1.)
    rate=(1+margin)*float(max(np.diff(D)/np.diff(t)))
    total=(1+margin)*float(D[-1])
    force=(1-margin)*Q*float(min(abs(X)))*a0/max(a)
    integral=min(T,b/force)+math.asinh(max(0.,force*T-b))/force
    age_weight=rate*integral
    S=Q**2/(K*a0**3)
    beta=S*age_weight*T**2
    zeta=beta*kmax*T/a0
    if general:
        guide=take("B")
        gyro_growth=Q*(1+margin)*max(abs(guide))*b*(min(T,b/force)+math.atan(max(0.,force*T-b))/force)
        momentum_factor=math.sqrt(2)*math.exp(gyro_growth)
    else:
        gyro_growth=0.;momentum_factor=1.
    sol,error,tail=positive_exponential(beta,zeta,momentum_factor)
    V=(2*math.pi/fundamental)**3
    w=6*mass/(action*V)
    # The 3D estimate uses one complex Fourier coefficient for both real quadratures.
    eta=((2.*math.sqrt(w*total) if general else math.sqrt(2*w*total))+2*w)/Q
    scale=S*T*eta
    rms=(1. if general else math.sqrt(2))*Q*T*scale*sol[1]
    # One representative of each +/- integer planar wavevector.
    N=math.floor(kmax/fundamental)
    modes=sum(2*math.isqrt(N*N-i*i)+1 for i in range(-N,N+1))
    modes=spatial_wavevector_pairs(N) if general else (modes-1)//2
    return dict(k_max=kmax,fundamental_comoving_k=fundamental,pairs_in_volume=float(D[-1]/w),
                number_of_independent_wavevector_pairs=modes,largest_source_density_slope=rate,
                wavevector_dimensions=(3 if general else 2),transverse_gyro_log_norm_upper=float(gyro_growth),
                minimum_axial_acceleration_coefficient=force,age_weighted_density_upper=age_weight,
                time_window=T,dimensionless_feedback=beta,dimensionless_density_coupling=zeta,
                unit_source_scaled_field=sol[0],unit_source_scaled_momentum=sol[1],
                single_mode_absolute_impulse_rms_upper=float(rms),all_modes_triangle_impulse_rms_upper=float(modes*rms),
                quadrature_relative_check=error,taylor_tail_absolute_bound=tail,
                input_margin=margin)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--history',type=Path,default=ROOT.parent/'history')
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    j=json.loads((args.history/'markers.json').read_text())
    out={}
    for name in j['runs']:
        marker=j['runs'][name]['markers']['0.01']
        with np.load(args.history/(name+'.npz')) as table:
            out[name]={label:bound(table,marker,j['action'],j['daughter_mass_L'],j['charge_to_mass'],32768.,fund,margin)
                       for label,fund,margin in [('one_wavelength',32768.,0.),('common_patch',16.,0.),
                                                 ('common_patch_padded',16.,.01)]}
            for label,fund,margin in [('all_directions_one_wavelength',32768.,0.),
                    ('all_directions_patch',16.,0.),('all_directions_patch_padded',16.,.01)]:
                out[name][label]=bound(table,marker,j['action'],j['daughter_mass_L'],j['charge_to_mass'],
                                       32768.,fund,margin,general=True)
    if any(spatial_wavevector_pairs(n)!=(sum(1 for x in range(-n,n+1) for y in range(-n,n+1)
            for z in range(-n,n+1) if 0 < x*x+y*y+z*z <= n*n)//2) for n in range(1,7)):
        raise ValueError('Lattice count failed its direct enumeration check')
    args.output.write_text(json.dumps({'runs':out},indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':main()
