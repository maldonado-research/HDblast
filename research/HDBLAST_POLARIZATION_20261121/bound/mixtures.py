"""Bound weak sideways pressure for arbitrary mixtures in the finite source space."""
import argparse
import json
import math
from pathlib import Path
import runpy

import numpy as np

HELPERS = runpy.run_path(str(Path(__file__).with_name('weak_pressure.py')))
source_nodes = HELPERS['source_nodes']
suffix = HELPERS['suffix']
columns = HELPERS['columns']
CHARGE = HELPERS['CHARGE']
KAPPAS = (2, 5, 7, 40, 258)
BUDGETS = (.005, .01, .02, .04, .08, .16)


def whiten(f, info, nodes):
    profiles = info['profile_order']
    gram = np.array(info['emitted_energy_gram'])
    if len(profiles)!=7 or not np.allclose(gram.T,gram) or not np.allclose(gram.diagonal(),1.,atol=3.e-6):
        raise ValueError('Unexpected source Gram')
    g = np.linalg.cholesky(gram)
    source = columns(f)
    tables = []
    for idx in range(2):
        tables.append(np.array([np.linalg.solve(g,np.array([nodes['interp'](source[p,k][idx]) for p in profiles]))
                                for k in KAPPAS]))
    # Verify E dot a is invariant under source whitening for a test coefficient vector.
    probe = np.cos(np.arange(len(profiles))*.83)
    direct = sum(c*nodes['interp'](source[p,7][0]) for c,p in zip(probe,profiles))
    whitened = (g.T @ probe) @ tables[0][KAPPAS.index(7)]
    if not np.allclose(direct,whitened,rtol=1.e-10,atol=1.e-10):
        raise ValueError('Source whitening changed the field')
    return *tables, float(np.linalg.eigvalsh(gram)[0])


def bound(nodes, ell, window, klist, E, B, epsilon):
    t,a,x,w,pstar=(nodes[k] for k in ('t','a','x','w','pstar'))
    dt=np.diff(t)
    r=math.sqrt(window)
    exterior=min(1.,r+ell); interior=max(0.,r-ell)
    # Whole sphere needs no density assumption. Otherwise Haar volume of |z1|<=r is r².
    mass_factor=1. if window==1. else (exterior/interior)**2
    powers=klist-2
    peaks=np.sqrt(klist-1.)*exterior**powers
    slopes=np.sqrt(klist-1.)*powers*exterior**(powers-1)
    # One square norm per (degree,profile); the two quadratures at a point are orthogonal.
    ntime=len(a)
    e=(E*peaks[:,None,None]).reshape(-1,ntime).T
    b=(B*peaks[:,None,None]).reshape(-1,ntime).T
    size=lambda z: np.linalg.norm(z,axis=1)
    upper=lambda y: np.maximum(y[:-1],y[1:])
    eh=epsilon*upper(size(e)); bh=epsilon*upper(size(b))
    dz=np.minimum(-x[:-1],-x[1:]) - bh
    v_bound=1.
    if np.all(dz>0):
        rate=np.max((eh+bh)/dz)
        ratio=math.hypot(rate,pstar*(1.+rate))
        v_bound=ratio/math.hypot(1.,ratio)
    lower=dt*dz*np.where(dz>0.,a[:-1],a[1:])
    tails=suffix(lower)
    axial=(CHARGE/a[-1])*(np.r_[tails[1:],0.]+np.minimum(0.,lower))-pstar
    gamma=np.hypot(1.,np.maximum(0.,axial))
    guide=upper(abs(nodes['guide']))
    magnetic=CHARGE/a[-1]*suffix(dt*a[1:]*(bh+guide*v_bound))
    travel=suffix(dt/a[:-1])
    angular_variation=np.minimum(2.*peaks[:,None],slopes[:,None]*travel[None,:])
    Eleft = E[:,:,:-1]*angular_variation[:,None,:]
    Eright = E[:,:,1:]*angular_variation[:,None,:]
    motion=epsilon*CHARGE/a[-1]*suffix(dt*a[1:]*np.sqrt(np.maximum(
        np.sum(Eleft**2,axis=(0,1)),np.sum(Eright**2,axis=(0,1)))))
    # I(s) is the signed integral of every whitened E coefficient. Between cells, a and E are linear.
    step=dt[:,None]*(2*a[:-1,None]*e[:-1]+a[:-1,None]*e[1:]+a[1:,None]*e[:-1]+2*a[1:,None]*e[1:])/6.
    I=suffix(step)
    # A direction of the 35-dimensional transfer response, common to all cohorts.
    # Any such choice gives a valid cap. The leading eigenvector keeps most of I.
    cov=(I.T*(w/gamma))@I
    direction=np.linalg.eigh(cov)[1][:,-1]
    def projected(tab):
        return tab-(tab@direction)[:,None]*direction
    diff=size(projected(I)); gradient=upper(size(projected(e)))
    # In each allowed birth cell use either neighboring endpoint, adding the whole cell.
    uncertain=np.minimum(diff,np.r_[diff[1:],0.])+dt*a[1:]*gradient
    electric=epsilon*CHARGE/a[-1]*uncertain
    pk=pstar+electric+motion+magnetic
    thermal=nodes['thermal']
    cap=mass_factor*np.dot(w,pk**2/gamma)/thermal
    return dict(pressure_upper_over_thermal=float(cap), maximum_transverse_speed=v_bound,
        conditional_mass_factor=mass_factor, reference_singular_fraction=float(direction@cov@direction/np.trace(cov)),
        normalized_rms_ledger={name:float(np.sqrt(mass_factor*np.dot(w,v*v/gamma)/thermal))
                               for name,v in [('birth',pstar*np.ones_like(pk)),('electric_direction_change',electric),
                                              ('motion',motion),('magnetic',magnetic)]})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--history',type=Path,required=True);p.add_argument('--fields',type=Path,required=True)
    p.add_argument('--gram',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    meta=json.loads((args.history/'markers.json').read_bytes()); infos=json.loads(args.gram.read_bytes())
    results={}
    for label in ('dedicated_native','dedicated_half','shared_native','shared_half','dedicated_spatial'):
        parent='dedicated_native' if label=='dedicated_spatial' else label
        with np.load(args.history/(parent+'.npz')) as file: h=dict(file)
        with np.load(args.fields/(label+'.npz')) as file: f=dict(file)
        nodes=source_nodes(h,f,meta['runs'][parent]['markers']['0.01'])
        ell=float(np.sum(np.diff(nodes['t'])/nodes['a'][:-1]))
        E,B,smallest=whiten(f,infos[label],nodes)
        test={}
        for low_only in (True,False):
            subset=np.array(KAPPAS)<=7 if low_only else np.ones(len(KAPPAS),bool)
            mode=','.join(str(k) for k in np.array(KAPPAS)[subset])
            test[mode]={str(window):{str(eps):bound(nodes,ell,window,np.array(KAPPAS)[subset],E[subset],B[subset],eps)
                        for eps in BUDGETS} for window in (.25,.50,.75,1.)}
        # For each observation region apply Markov to the nonnegative local
        # smaller horizontal stress per tagged particle, weighted by the local number.
        entire=test[','.join(map(str,KAPPAS))]
        fraction={}
        for eps in BUDGETS:
            limits={}
            for target_fraction in (.1,.5,.9):
                lower=[]
                for window, row in entire.items():
                    zone=float(window)
                    retained=1. if zone==1. else max(0.,math.sqrt(zone)-ell)**2
                    cap=row[str(eps)]['pressure_upper_over_thermal']
                    lower.append(retained*max(0.,1.-cap/target_fraction))
                limits[str(target_fraction)]=max(lower)
            fraction[str(eps)]=limits
        results[label]=dict(particle_fraction_at_sites_below_lateral_target=fraction,
                            pressure_target=nodes['thermal'],max_travel=ell,profiles=infos[label]['profile_order'],
                            smallest_gram_eigenvalue=smallest,tests=test)
    args.output.write_text(json.dumps(results,indent=2,sort_keys=True,allow_nan=False)+'\n')

if __name__=='__main__':
    main()
