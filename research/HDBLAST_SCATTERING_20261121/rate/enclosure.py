"""Enclose collision exposure on the interpolated homogeneous field histories by birth ages."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import runpy

import numpy as np


POLE2 = 8300.**2
WIDTH = .30
Q = 1.e7
# Both constants round the analytic maxima upward.
KCAP = 1.118569e-6
PHICAP = (1.-8./POLE2)**(-.5)*(1+1.e-14)
BIRTH_CAP = np.sqrt(1.05**2-1.)*(1+1.e-14)


def kernel_bound(zlo, pmax):
    """Age-by-age bound on the invariant kernel with arbitrary injection directions."""
    elo = np.hypot(1.,zlo)
    ehi = np.hypot(1.,pmax)
    ep = elo[:,None]*elo
    # With both streams forced forward in opposite hemispheres, transverse
    # relative rapidities can cancel, leaving a one-dimensional lower bound.
    both_forward = (zlo[:,None]>0)&(zlo>0)
    lo = np.where(both_forward,2*(1+ep+zlo[:,None]*zlo),4.)
    hi = 2*(1+ehi[:,None]*ehi+pmax[:,None]*pmax)
    # y=(s-M^2)/[width*(s+2)] increases with s. Bound its distance to zero.
    ymin = np.maximum(0.,np.maximum((POLE2-hi)/(WIDTH*(hi+2.)),(lo-POLE2)/(WIDTH*(lo+2.))))
    beta_lo = np.sqrt(np.maximum(0.,1.-4./lo))
    beta_hi = np.sqrt(np.maximum(0.,1.-4./hi))
    beta = np.clip(ymin,beta_lo,beta_hi)
    denom = ymin**2+beta**2
    phi = np.divide(beta,denom,out=np.full_like(beta,PHICAP),where=denom>0.)
    return np.minimum(KCAP,6*np.pi*np.minimum(phi,PHICAP)/ep)*(1+1.e-12)


def coefficient_limits(h, f, marker, case, fit, parent):
    t = np.r_[h['t'][h['t']<marker['time']],marker['time']]
    def at(values):
        return np.interp(t,h['t'],values)
    a = at(h['a'])
    if np.any(np.diff(a)<0) or np.max(at(h['R']))>1.05+1e-13:
        raise ValueError('Recorded radius or birth cap changed')
    early,bearly = parent['mode'](f,2,'early')
    late,blate = parent['mode'](f,2,'late_birth')
    x,u,v = fit if case=='two_pressure_fit' else (.04/np.sqrt(2),0.,.04/np.sqrt(2))
    electric = np.column_stack((at(x*early+u*late),at(v*late),at(h['X'])))
    bperp = np.column_stack((at(x*bearly+u*blate),at(v*blate)))
    axial_force = -np.max(electric[:,2])-np.max(np.linalg.norm(bperp,axis=1))
    force_norm = np.max(np.linalg.norm(electric,axis=1))
    if axial_force<=0 or np.any(np.diff(h['D'])<0):
        raise ValueError('Expected a single axial force sign and monotone injection')
    # Integrals of piecewise lower and upper envelopes between original knots.
    # The norm of a linearly interpolated vector is convex.
    max_b = np.maximum(np.linalg.norm(bperp[:-1],axis=1),np.linalg.norm(bperp[1:],axis=1))
    axial_lo = np.minimum(-electric[:-1,2],-electric[1:,2])-max_b
    total_hi = np.maximum(np.linalg.norm(electric[:-1],axis=1),np.linalg.norm(electric[1:],axis=1))
    proper_lo = (1-1.e-12)*a[:-1]*axial_lo
    proper_hi = (1+1.e-12)*a[1:]*total_hi
    acc_lo = np.r_[0.,np.cumsum(np.diff(t)*proper_lo)]
    acc_hi = np.r_[0.,np.cumsum(np.diff(t)*proper_hi)]
    if np.any(proper_lo<0):
        raise ValueError('Axial force changed sign')
    return dict(t=t, acc_lo=acc_lo, acc_hi=acc_hi, axial_min=float(axial_force*(1-1.e-12)),
                electric_max=float(force_norm*(1+1.e-12)),amin=float(a.min()),amax=float(a.max()))


def estimate(h, marker, limits, factor, cells):
    t = np.linspace(h['t'][0],marker['time'],cells+1)
    a = np.interp(t,h['t'],h['a'])
    z = np.interp(t,h['t'],h['D'])/marker['D']
    weights = np.diff(z)
    lower_integral = np.interp(t,limits['t'],limits['acc_lo'])
    upper_integral = np.interp(t,limits['t'],limits['acc_hi'])
    depth = initial = 0.
    for j,dt in enumerate(np.diff(t)):
        births = weights[:j+1]
        # Accumulating lower/upper positive envelopes permits uncertain
        # birth and observation times anywhere within their respective cells.
        axial_impulse = np.maximum(0.,lower_integral[j]-lower_integral[1:j+2])
        total_impulse = upper_integral[j+1]-upper_integral[:j+1]
        birth_cap = a[1:j+2]*BIRTH_CAP
        zlo = np.maximum(0.,(Q*axial_impulse-birth_cap)/a[j+1]-1.e-7)
        pmax = (birth_cap+Q*total_impulse)/a[j]+1.e-7
        prefactor = factor*marker['D']/a[j]**3
        r0 = prefactor*np.einsum('i,ij,j',births,kernel_bound(zlo,pmax),births,optimize=False)
        lam = prefactor*KCAP*z[j+1]
        gain = dt*lam
        exprel = np.expm1(gain)/gain if gain>0 else 1.
        depth = depth*np.exp(gain)+dt*r0*exprel
        initial += dt*r0
    # With R and lambda bounded above cellwise, cumulative count <=
    # integral(Rupper + 2*lambda_upper*Yupper) = 2Yupper - integral Rupper.
    return dict(cells=cells,unscattered_pair_source_upper=float(initial),
                ever_hit_fraction_upper=float(depth),mean_encounter_upper=float(2*depth-initial))


def test_kernel_enclosure():
    """Check against direct kernels with prescribed lower and upper momentum supports."""
    rng = np.random.default_rng(448389)
    # Include weakly oriented newborns, slow charges, and particles around the pole.
    lo = np.r_[np.zeros(12),10**rng.uniform(-1,5,116)]
    hi = np.maximum(lo,1.e-6)+rng.uniform(0,1.2,len(lo))*(lo+2500.)
    cap = kernel_bound(lo,hi)
    positive = rng.uniform(lo,hi)
    negative = rng.uniform(lo,hi)
    positive[:12] *= rng.choice([-1.,1.],12)
    negative[:12] *= rng.choice([-1.,1.],12)
    def directions(z):
        lateral = rng.uniform(0.,1.,len(z))*np.sqrt(np.maximum(0.,hi**2-z**2))
        azimuth = rng.uniform(0.,2*np.pi,len(z))
        return np.column_stack((lateral*np.cos(azimuth),lateral*np.sin(azimuth),z))
    p = directions(positive)
    q = directions(-negative)
    ep = np.hypot(1,np.linalg.norm(p,axis=1))[:,None]*np.hypot(1,np.linalg.norm(q,axis=1))
    s = np.maximum(4.,2*(1+ep-p@q.T))
    beta = np.sqrt(1-4/s)
    g = WIDTH*(s+2)
    k = 6*np.pi*beta*g**2/(ep*((s-POLE2)**2+(beta*g)**2))
    if np.any(k>cap):
        raise ValueError('Direct invariant kernel exceeded its enclosure')



def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--parent',type=Path,required=True)
    p.add_argument('--fields',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--jobs',type=int,default=2)
    args = p.parse_args()
    test_kernel_enclosure()
    base = args.parent.resolve()
    api = runpy.run_path(str(base/'orbits/pairs.py'))
    hist = base/'parent/parent/history'
    refs = api['read_inputs'](hist,args.fields)
    meta = json.loads((hist/'markers.json').read_text())
    source = json.loads((base/'feedback/FEEDBACK.json').read_text())
    distances = {}
    for r in source['runs']:
        if r['steps']==512 and r['pairs']==4096 and not r['feedback']:
            distances[(r['history'],r['case'])]=r['coarse_cell_total_variation']
    def one(key):
        name,case = key
        h,f,m = refs[name]
        coeff = coefficient_limits(h,f,m,case,api['FIT'][name.split('_')[0]+'_half'][1],api)
        factor = meta['runs'][name]['channel_fraction']*meta['action']/(6*meta['daughter_mass_L']**3)
        rows = [estimate(h,m,coeff,factor,n) for n in (64,128,256)]
        rows.append(estimate(h,m,coeff,factor,512))
        if name.endswith('half'):
            rows.append(estimate(h,m,coeff,factor,768))
            distance = distances[name,'two_pressure_fit' if case=='two_pressure_fit' else 'equal_early_and_late_4pct']
            for r in rows:
                r.update(collisionless_tv=distance,
                         minimum_20_cell_tv_estimate=max(0.,distance-r['ever_hit_fraction_upper']))
        if np.any(np.diff([r['ever_hit_fraction_upper'] for r in rows[:4]])>1.e-12):
            raise ValueError('Failure of nested upper enclosures')
        return dict(history=name,case=case,radius_min=coeff['amin'],radius_max=coeff['amax'],
                    minimum_signed_axial_force_per_charge=coeff['axial_min'],peak_electric_norm=coeff['electric_max'],bounds=rows)
    cases=[(key,case) for key in refs for case in ('two_pressure_fit','early_late_4pct')]
    with ThreadPoolExecutor(args.jobs) as pool:
        rows=list(pool.map(one,cases))
    args.output.write_text(json.dumps(dict(global_kernel_cap=KCAP,runs=rows),sort_keys=True,indent=2,
                                     allow_nan=False)+'\n')


if __name__=='__main__':
    main()
