"""Guaranteed weak local transverse pressure for two driven spherical Maxwell modes."""
import argparse
import json
import math
from itertools import product
from pathlib import Path

import numpy as np


CHARGE = 1.e7


def suffix(v):
    return np.cumsum(v[::-1], axis=0)[::-1]


def source_nodes(h, f, marker):
    stop = marker['time']
    t = np.r_[h['t'][h['t'] < stop], stop]
    if not np.array_equal(h['t'], f['t']):
        raise ValueError('Need fields on benchmark knots')
    def interp(x):
        return np.interp(t, h['t'], x)
    d = interp(h['D'])
    if d[0] != 0. or np.any(np.diff(d) < 0.) or not np.isclose(d[-1], marker['D']):
        raise ValueError('Inconsistent birth measure')
    a, x = interp(h['a']), interp(h['X'])
    if np.any(np.diff(a) < 0.) or np.max(x) >= 0.:
        raise ValueError('Bound assumes increasing radius and one axial electric sign here')
    return dict(t=t, a=a, x=x, guide=interp(h['B']), w=np.diff(d)/d[-1],
                pstar=np.sqrt(max(float(np.max(h['R'][:h['t'].searchsorted(stop)+1]))**2-1., 0.)),
                thermal=marker['theta']/math.hypot(1., float(interp(h['u'])[-1])),
                interp=interp)


def cross_coefficient(q, integrals, pointwise, a, dt):
    """Maximum of |integral from birth to T| throughout each possible birth cell."""
    projected = pointwise @ q
    left = integrals @ q
    right = np.r_[left[1:], 0.]
    largest = np.maximum(abs(left), abs(right))
    cross = projected[:-1]*projected[1:] < 0.
    if np.any(cross):
        value = projected[:-1][cross]
        frac = value/(value-projected[1:][cross])
        # Integral from the left endpoint to the zero of the linear force.
        partial = dt[cross]*frac*value*(a[:-1][cross]/2.+np.diff(a)[cross]*frac/6.)
        largest[cross] = np.maximum(largest[cross], abs(left[cross]-partial))
    return largest


def pressure(nodes, E, B, kappas, rms):
    t, a, x, guide, w, pstar = (nodes[k] for k in ('t','a','x','guide','w','pstar'))
    dt = np.diff(t); n = len(dt)
    sup = np.sqrt(np.array(kappas)-1.)
    norm = rms*sup
    e = (np.asarray(E)*norm[:,None]).T
    b = (np.asarray(B)*norm[:,None]).T
    maxima = lambda v: np.maximum(abs(v[:-1]), abs(v[1:]))
    me, mb = maxima(e), maxima(b)
    eq = kappas[0] == kappas[1]  # Same-κ quadratures are orthogonal everywhere and have equal norm.
    combine = (lambda v: np.max(v,axis=1)) if eq else (lambda v: np.linalg.norm(v,axis=1))
    bh = combine(mb)
    eh = combine(me)
    dz = np.minimum(-x[:-1], -x[1:]) - bh
    r = float(np.max((eh+bh)/dz)) if np.all(dz>0.) else None
    if r is None:
        vmax = 1.
    else:
        vratio = math.hypot(r, pstar*(1.+r))
        vmax = vratio/math.hypot(1.,vratio)
    # Source sampling at the final point, evaluated in the left-invariant frame.
    travel = suffix(dt/a[:-1])
    changes = np.minimum(2.,travel[:,None]*(np.array(kappas)-2.))
    # Equal-κ axes rotate together under the same complex multiplier, including along a path.
    motion = CHARGE/a[-1]*suffix(dt*a[1:]*combine(me*changes))
    magnetic = CHARGE/a[-1]*suffix(dt*a[1:]*(bh+maxima(guide)*vmax))
    longitudinal = np.where(dz>=0.,a[:-1],a[1:])*dz*dt
    complete = suffix(longitudinal)
    lower_z = (CHARGE/a[-1])*(np.r_[complete[1:],0.]+np.minimum(0.,longitudinal))-pstar
    gamma_lower = np.hypot(1.,np.maximum(0.,lower_z))
    force_integral = dt[:,None]*(2.*a[:-1,None]*e[:-1]+a[:-1,None]*e[1:]+
                                 a[1:,None]*e[:-1]+2.*a[1:,None]*e[1:])/6.
    integrals = suffix(force_integral)
    # Choose a common reference color for the two patterns. At each observation
    # point the tested direction is perpendicular to that local vector.
    def upper(theta, diagnostics=False):
        q=np.array([-math.sin(theta),math.cos(theta)])
        projected = cross_coefficient(q,integrals,e,a,dt)
        electric = projected*(CHARGE/a[-1])
        pk=pstar+electric+motion+magnetic
        raw=float(np.dot(w,pk**2/gamma_lower))
        if not diagnostics: return raw
        components=[pstar*np.ones(n), electric, motion, magnetic]
        return dict(weak_pressure_upper_per_particle=raw, upper_over_drifting_thermal=raw/nodes['thermal'],
                    transverse_speed_upper=vmax,
                    normalized_rms_ledger={k:float(np.sqrt(np.dot(w,v*v/gamma_lower)/nodes['thermal']))
                           for k,v in zip(('birth','changing_electric_direction','travel','magnetic'),components)})
    # Add directions from the signed integrated source to a fixed angular grid.
    # Using a finite list of witnesses avoids fitting an ill-determined optimum for tiny fields.
    cov=(integrals.T*(w/gamma_lower))@integrals
    theta=.5*math.atan2(2.*cov[0,1],cov[0,0]-cov[1,1])
    trial=np.r_[np.linspace(0.,np.pi,19,endpoint=False), theta, np.arctan2(integrals[::max(1,n//9),1],integrals[::max(1,n//9),0])]
    return upper(trial[np.argmin([upper(z) for z in trial])],True)


def columns(f):
    if 'profile' in f:
        return {(str(p),int(k)):(e,b) for p,k,e,b in zip(f['profile'],f['kappa'],f['E'],f['B'])}
    return {('early',int(k)):(e,b) for k,e,b in zip(f['kappa'],f['E_cos'],f['B_sin'])}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--history',type=Path,required=True); p.add_argument('--fields',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    meta=json.loads((args.history/'markers.json').read_bytes())
    result={}
    for label in ('dedicated_half','shared_half','dedicated_native','shared_native'):
        with np.load(args.history/(label+'.npz')) as file: h=dict(file)
        with np.load(args.fields/(label+'.npz')) as file: f=dict(file)
        mk=meta['runs'][label]['markers']['0.01']; nodes=source_nodes(h,f,mk)
        src={key:tuple(nodes['interp'](v) for v in pair) for key,pair in columns(f).items()}
        rows=[]
        profiles=sorted(set(p for p,k in src))
        for profile in profiles:
            for ka,kb in product(sorted({k for p,k in src}),repeat=2):
                if ('early',ka) not in src or (profile,kb) not in src: continue
                ea,ba=src['early',ka]; eb,bb=src[profile,kb]
                caps={}
                for rms in (.01,.02,.04,.08,.16):
                    row=pressure(nodes,(ea,eb),(ba,bb),(ka,kb),rms)
                    caps[str(rms)]=row['upper_over_drifting_thermal']
                rows.append(dict(first_mode=ka,second_mode=kb,second_profile=profile,pressure_upper_over_thermal=caps))
        result[label]=dict(target_lateral_pressure=nodes['thermal'],cases=rows)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
