"""Follow a pair of polarized outer-source harmonics through collisionless S³ particle trajectories."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import math
from pathlib import Path
import subprocess
import tempfile

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import kv, kve
from scipy.stats import qmc


Q=1.e7
SCRAMBLE=120835
HERE=Path(__file__).resolve().parent


def directions(c,d,cap,prob):
    """Sample the inherited angular acceptance min(1, |c-d mu|/cap)."""
    mu=2.*prob-1.
    use=(d!=0.) & (cap>0.)
    c,d,cap,u=c[use],d[use],cap[use],prob[use]
    if not c.size: return mu
    cuts=np.sort(np.clip(np.stack((-np.ones_like(c),(c-cap)/d,c/d,(c+cap)/d,np.ones_like(c)),axis=1),-1.,1.),axis=1)
    rate=np.minimum(1.,np.abs(c[:,None]-d[:,None]*cuts)/cap[:,None])
    widths=np.diff(cuts,axis=1)
    total=np.column_stack((np.zeros_like(c),np.cumsum(.5*(rate[:,:-1]+rate[:,1:])*widths,axis=1)))
    target=u*total[:,-1]
    j=np.minimum(3,(total[:,1:]<=target[:,None]).sum(axis=1)); i=np.arange(len(j))
    area=(target-total[i,j])/widths[i,j]
    left,right=rate[i,j],rate[i,j+1]
    dx=2*area/(left+np.sqrt(np.maximum(0.,left*left+2*(right-left)*area)))
    mu[use]=cuts[i,j]+dx*widths[i,j]
    return mu


def seeds(h,mark,power,seed):
    r=qmc.Sobol(6,scramble=True,seed=seed).random_base2(power)
    born=np.interp(r[:,0]*mark['D'],h['D'],h['t'])
    at=lambda key: np.interp(born,h['t'],h[key])
    mu=directions(at('c'),at('d'),at('target'),r[:,4])
    pstar=np.sqrt(np.maximum(0.,at('R')**2-1.))
    pt=pstar*np.sqrt(np.maximum(0.,1-mu*mu))
    z1=np.sqrt(r[:,1])*np.exp(2j*np.pi*r[:,2])
    z2=np.sqrt(1.-r[:,1])*np.exp(2j*np.pi*r[:,3])
    az=2.*np.pi*r[:,5]
    return np.column_stack((born,z1.real,z1.imag,z2.real,z2.imag,
                            at('a')*pt*np.cos(az),at('a')*pt*np.sin(az),at('a')*pstar*mu))


def mode(f,k,profile):
    select=np.flatnonzero((f['kappa']==k)&(f['profile']==profile))
    if len(select)!=1: raise ValueError(f'Missing or repeated mode {k}, {profile}')
    return f['E'][select[0]], f['B'][select[0]]


def field_table(h,f,mark,k1,k2,p1,p2,steps):
    if not np.array_equal(h['t'],f['t']): raise ValueError('Unmatched field trace')
    t=np.linspace(h['t'][0],mark['time'],steps+1)
    e1,b1=mode(f,k1,p1);e2,b2=mode(f,k2,p2)
    return np.column_stack([t]+[np.interp(t,h['t'],v) for v in (h['a'],h['X'],h['B'],e1,b1,e2,b2)])


def shape(k,z1,z2,style):
    if style=='highest':
        return np.sqrt(k-1.)*np.conj(z1**(k-2))
    if style!='distributed':
        raise ValueError(style)
    n=k-2
    return np.conj(sum(np.exp(1j*np.pi*np.sqrt(2.)*j*j)*np.sqrt(float(math.comb(n,j))) * z1**j * z2**(n-j) for j in range(n+1)))


def gram_projection(h,tab,sign,y,j,k1,k2):
    """A field-defined pressure comparison basis at the particle's observation point."""
    t,a,x,b,e1,b1,e2,b2=tab.T
    weights=np.diff(np.interp(t,h['t'],h['D']))
    launch=lambda v: -cumulative_trapezoid((a*v)[::-1],t[::-1],initial=0.)[::-1]
    center=lambda v:.5*(v[:-1]+v[1:])
    axial=center(launch(x))*Q/a[-1]
    weight=weights/np.hypot(1.,axial)
    # Predicted transverse impulse of a relativistic particle moving along the axial drive.
    vz=sign*np.sign(np.mean(x))
    distance=-cumulative_trapezoid((1./a)[::-1],t[::-1],initial=0.)[::-1]
    u=center(launch((e1+1j*vz*b1)*np.exp(-1j*(k1-2)*vz*distance)))
    v=center(launch((e2+1j*vz*b2)*np.exp(-1j*(k2-2)*vz*distance)))
    spin=y*y*np.sum(weight*u*u)+2.*y*j*np.sum(weight*u*v)+j*j*np.sum(weight*v*v)
    return .5*np.angle(spin)


def pressure(h,mark,tab,spec,raw):
    k1,k2,p1,p2,a1,a2,phase=spec[:7]; style=spec[7] if len(spec)>7 else "highest"
    sign=raw[:,0]
    z=raw[:,1]+1j*raw[:,2]
    norm=abs(z)**2+raw[:,3]**2+raw[:,4]**2
    z2=raw[:,3]+1j*raw[:,4]
    y1=a1*shape(k1,z,z2,"highest");y2=a2*np.exp(1j*phase)*shape(k2,z,z2,style)
    p=raw[:,5:];g=np.hypot(1.,np.linalg.norm(p,axis=1));pt=np.sort(np.linalg.norm(p[:,:2],axis=1))
    theta=mark['theta']; w=np.hypot(1.,pt)
    cdf=1.-w*w*kv(2,w/theta)/kv(2,1./theta);rank=np.arange(len(pt))/len(pt)
    ks=np.maximum(np.abs(cdf-rank),np.abs(cdf-(rank+1./len(pt))))
    drift=float(np.interp(mark['time'],h['t'],h['u']));boost=np.hypot(1.,drift)
    thermal=kve(1,1./theta)/kve(2,1./theta)+3.*theta
    target_energy=(thermal+theta)*boost-theta/boost
    target_p=theta/boost
    def covariance(matrix,energy):
        return (matrix.T/energy)@matrix/len(energy)/target_p
    def tensor(arr):
        eig=np.linalg.eigvalsh(arr)
        return [float(v) for v in eig]
    aligned=[]; pooled=[]; rest=[]; axial=[]; raw_cov=[]
    for s in (1,-1):
        use=sign==s;mom=p[use];energy=g[use]
        angle=gram_projection(h,tab,s,y1[use],y2[use],k1,k2)
        ca=np.cos(angle);sa=np.sin(angle)
        oriented=np.column_stack((mom[:,0]*ca+mom[:,1]*sa,-mom[:,0]*sa+mom[:,1]*ca))
        A=covariance(oriented,energy);B=covariance(mom[:,:2],energy)
        aligned.append(tensor(A)); pooled.append(tensor(B)); raw_cov.append(B)
        tr_rest=np.column_stack((oriented,boost*mom[:,2]-s*drift*energy))
        rest.append(tensor(covariance(tr_rest,energy)))
        wall=covariance(mom,energy)
        axial.append(wall[2,2]/wall.trace())
    if np.max(abs(norm-1.))>2.e-11: raise ValueError('Orbit left the unit sphere')
    return dict(per_sign_pooled_lateral_covariance_over_target=np.array(raw_cov).tolist(),
                ks_transverse=float(np.max(ks)),ks_cutoff_in_mc=float(pt[np.argmax(ks)]),
                energy_over_drifting_thermal=float(np.mean(g)/target_energy),
                aligned_lateral_pressure_over_thermal=np.mean(aligned,axis=0).tolist(),
                pooled_lateral_pressure_over_thermal=np.mean(pooled,axis=0).tolist(),
                reference_frame_stress_eigenvalues_over_thermal=np.mean(rest,axis=0).tolist(),
                axial_pressure_share=float(np.mean(axial)),spherical_norm_error=float(np.max(abs(norm-1.))))


def one(h,f,mark,spec,power,steps,seed,exe,tmp):
    k1,k2,p1,p2,a1,a2,phase=spec[:7]; style=spec[7] if len(spec)>7 else "highest"
    with tempfile.TemporaryDirectory(dir=tmp) as wd:
        path=Path(wd);tab=field_table(h,f,mark,k1,k2,p1,p2,steps)
        tab.tofile(path/'fields');seeds(h,mark,power,seed).tofile(path/'births')
        subprocess.run([str(exe),str(path/'fields'),str(path/'births'),str(k1),str(k2),str(a1),str(a2),str(phase),style,str(path/'orbits')],check=True)
        raw=np.fromfile(path/'orbits',dtype='=f8').reshape(-1,8)
        if len(raw)!=2*2**power or not np.isfinite(raw).all():raise ValueError('Incomplete trajectory set')
        return dict(k1=k1,k2=k2,profile1=p1,profile2=p2,amplitude1=a1,amplitude2=a2,phase2=phase,second_shape=style,
                    pairs=2**power,steps=steps,seed=seed,metrics=pressure(h,mark,tab,spec,raw))


def conditional(h,f,mark,spec,positions_power,birth_power,steps,seed,exe,tmp):
    """Evaluate actual local tangent stress at a set of common final physical positions."""
    k1,k2,p1,p2,a1,a2,phase=spec[:7];style=spec[7] if len(spec)>7 else 'highest'
    tags=seeds(h,mark,birth_power,seed+8153)
    r=qmc.Sobol(3,scramble=True,seed=seed).random_base2(positions_power)
    z1=np.sqrt(r[:,0])*np.exp(2j*np.pi*r[:,1]);z2=np.sqrt(1.-r[:,0])*np.exp(2j*np.pi*r[:,2])
    target=np.column_stack((z1.real,z1.imag,z2.real,z2.imag))
    tiled=np.tile(tags,(len(r),1))
    tiled[:,1:5]=np.repeat(target,len(tags),axis=0)
    with tempfile.TemporaryDirectory(dir=tmp) as wd:
        path=Path(wd); tab=field_table(h,f,mark,k1,k2,p1,p2,steps)
        tab.tofile(path/'fields');tiled.tofile(path/'tags')
        subprocess.run([str(exe),str(path/'fields'),str(path/'tags'),str(k1),str(k2),str(a1),str(a2),str(phase),style,str(path/'points'),'conditional'],check=True)
        raw=np.fromfile(path/'points',dtype='=f8').reshape(len(r),len(tags),2,11)
    # The complete local tensor at each endpoint, then its eigenvalues. Conditional number
    # weighting cancels in the final number-weighted volume average by homogeneity of lambda_min.
    mom=raw[:,:,:,5:8];det=raw[:,:,:,8];weights=1./det
    g=np.hypot(1.,np.linalg.norm(mom,axis=-1))
    theta=mark['theta'];drift=float(np.interp(mark['time'],h['t'],h['u']));boost=np.hypot(1.,drift)
    target_pressure=theta/boost
    T=np.einsum('xbsi,xbsj,xbs->xsij',mom[:,:,:,:2],mom[:,:,:,:2],weights/g)/len(tags)/target_pressure
    densities=weights.mean(axis=1); local=np.linalg.eigvalsh(T)/densities[:,:,None]
    avg=(local*densities[:,:,None]).sum(axis=(0,1))/densities.sum()
    thermal=kve(1,1./theta)/kve(2,1./theta)+3.*theta
    pt=np.linalg.norm(mom[:,:,:,:2],axis=-1).ravel();idx=np.argsort(pt)
    w=np.hypot(1.,pt[idx]);F=1.-w*w*kv(2,w/theta)/kv(2,1./theta)
    ordered=weights.ravel()[idx]/weights.sum();cdf=np.cumsum(ordered)
    return dict(k1=k1,k2=k2,profile1=p1,profile2=p2,amplitude1=a1,amplitude2=a2,phase2=phase,second_shape=style,
                endpoint_count=len(r),birth_samples_per_endpoint_per_sign=len(tags),steps=steps,seed=seed,
                local_lateral_pressure_over_thermal=avg.tolist(),
                local_weak_pressure_volume_quantiles_10_50_90=np.quantile(local[:,:,0],[.1,.5,.9]).tolist(),
                ks_transverse=float(np.maximum(abs(cdf-F),abs(cdf-ordered-F)).max()),
                energy_over_drifting_thermal=float(np.sum(weights*g)/weights.sum()/((thermal+theta)*boost-theta/boost)),
                mean_relative_density=float(densities.mean()),min_spatial_map_det=float(det.min()),max_spatial_map_det=float(det.max()),
                max_frobenius_map_jacobian_minus_identity=float(raw[:,:,:,10].max()),
                max_endpoint_error=float(raw[:,:,:,9].max()))


def read_inputs(hdir,fdir):
    meta=json.loads((hdir/'markers.json').read_bytes())
    inputs={}
    for label in ('dedicated_half','shared_half','dedicated_native','shared_native'):
        with np.load(hdir/(label+'.npz')) as z:h=dict(z)
        with np.load(fdir/(label+'.npz')) as z:f=dict(z)
        inputs[label]=(h,f,meta['runs'][label]['markers']['0.01'])
    return inputs


# These fixed coefficients are repeated on independently scrambled birth samples.
# The optional refit repeats the coarse moment fits on the pilot scramble.
FIT={
    'dedicated_half':(.01692449608886012, [.23687382734432305,34.602604569781604,12.576440954103287]),
    'shared_half':(.009276723212172279, [.17001443590833457,21.413907155127227,5.818149361182037]),
}
STORED=[(72,'t5.085',.0015),(160,'t5.155',.0035),(200,'t5.205',.0065),
        (258,'t5.205',.014),(360,'t5.205',.13),(515,'t5.225',.52)]


def uniform_pair(v):
    return (2,2,'early','late_birth',v[0],float(np.hypot(*v[1:])),float(np.arctan2(v[2],v[1])))


def refit(inputs,exe,tmp):
    from scipy.optimize import brentq,least_squares
    summary={}
    for label,estimate in [('dedicated_half',(.236,36.62,.349)),('shared_half',(.170,22.167,.265))]:
        def moments(spec,power):
            return one(*inputs[label],spec,power,128,SCRAMBLE,exe,tmp)['metrics']
        def residual(v):
            m=np.mean(moments(uniform_pair(v),9)['per_sign_pooled_lateral_covariance_over_target'],axis=0)
            return [m[0,0]-1.,m[0,1],m[1,1]-1.]
        def trace(a):
            result=moments((7,7,'early','early',a,a,np.pi/2),10)
            return sum(result['pooled_lateral_pressure_over_thermal'])-2.
        a,b,phase=estimate
        v0=np.array([a,b*np.cos(phase),b*np.sin(phase)])
        fitted=least_squares(residual,v0,x_scale=np.maximum(abs(v0),.001),gtol=1e-8,ftol=1e-9,xtol=1e-9)
        if not fitted.success or max(abs(fitted.fun))>1e-6: raise RuntimeError('Moment solve failed')
        summary[label]=[brentq(trace,.004,.04,xtol=2.e-10),fitted.x.tolist()]
    return summary


def active_inputs(path):
    with np.load(path) as z: h=dict(z)
    if not np.isfinite(h['marker_time']): raise ValueError('Absent work marker')
    f=dict(t=h['t'],kappa=h['kappa'],profile=np.array(['first','second']),E=h['mode_E'],B=h['mode_B'])
    t=float(h['marker_time'])
    mark={key:float(np.interp(t,h['t'],h[key])) for key in ('a','theta','D')}
    mark['time']=t
    return h,f,mark


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--history',type=Path,required=True)
    p.add_argument('--fields',type=Path,required=True,help='Primary trace directory')
    p.add_argument('--scan',type=Path,required=True,help='Storage scan trace directory')
    p.add_argument('--source',type=Path,required=True,help='Live source trace directory')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--jobs',type=int,default=1)
    p.add_argument('--refit',action='store_true',help='Also repeat the pilot fits of the pressure controls')
    args=p.parse_args()
    if args.jobs<1: p.error('--jobs must be positive')
    inp=read_inputs(args.history,args.fields)
    runs=[]; conditioned=[]
    def track(name,label,spec,pow=11,steps=192,seed=SCRAMBLE+1):
        runs.append((name,label,spec,(pow,steps,seed)))
    def local(name,label,spec,locpow=7,birthpow=6,steps=144,seed=SCRAMBLE+1):
        conditioned.append((name,label,spec,(locpow,birthpow,steps,seed)))
    for label in ['dedicated_native','dedicated_half','shared_half','dedicated_spatial']:
        key='dedicated_native' if label=='dedicated_spatial' else label
        h,_,mark=inp[key]
        with np.load(args.scan/(label+'.npz')) as z: inp['scan_'+label]=(h,dict(z),mark)

    for label,(rot,v) in FIT.items():
        scale=1. if label=='dedicated_half' else .59
        same=(7,7,'early','early',rot,rot,np.pi/2)
        inverse=uniform_pair(v)
        track('zero_transverse',label,(2,2,'early','early',0.,0.,0.))
        track('single_uniform',label,(2,2,'early','early',.10*scale,0.,0.))
        track('uniform_rotation',label,(2,2,'early','early',.08*scale,.06*scale,np.pi/2))
        track('pooled_pressure_match',label,same,12,256)
        local('pooled_pressure_match',label,same,6,7,160)
        track('uniform_local_pressure_match',label,inverse,12,512)
        track('uniform_local_pressure_match_second_scramble',label,inverse,11,256,SCRAMBLE+17)
        low=(2,40,'early','prepare',.08*scale,.004 if scale==1. else .001,np.pi/2)
        track('low_cdf_discrepancy_pilot_selection',label,low,12,256)
        for k,delay,amp in STORED:
            spec=(2,k,'early',delay,.05*scale,amp*scale,np.pi/2,'distributed')
            key='scan_'+label
            track('distributed_storage',key,spec,11,192)
            if k in (258,515):
                local('distributed_storage',key,spec)
            if k==515:
                local('distributed_storage_scramble',key,spec,6,7,192,SCRAMBLE+41)
                track('distributed_storage_step_check',key,spec,11,384)
                if scale==1.:
                    for other in ('scan_dedicated_native','scan_dedicated_spatial'):
                        track('distributed_storage',other,spec,11,192)
                    for other in ('scan_dedicated_native','scan_dedicated_spatial'):
                        local('distributed_storage_radial_grid',other,spec,6,6,144)
    active=['pair_0.04_active','pair_0.08_active','pair_0.08_half','pair_0.08_radial',
            'pair_0.16_long','delayed_0.0100','simultaneous_time_control',
            'prepare_0.005','prepare_0.010','prepare_0.020']
    for label in active:
        inp[label]=active_inputs(args.source/(label+'.npz'))
        k1,k2=map(int,inp[label][1]['kappa'])
        spec=(k1,k2,'first','second',1.,1.,np.pi/2)
        track('evolving_mean_source',label,spec,11,196,SCRAMBLE)
        if label in ('pair_0.08_half','prepare_0.010'):
            local('evolving_mean_source',label,spec,6,7,160)
    with tempfile.TemporaryDirectory(prefix='polar-orbits-') as wd:
        tmp=Path(wd);exe=tmp/'sphere'
        subprocess.run(['g++','-std=c++17','-O3',str(HERE/'sphere.cpp'),'-o',str(exe)],check=True)
        def evaluate(item):
            fn,name,label,spec,settings=item
            return dict(case=name,history=label,**fn(*inp[label],spec,*settings,exe,tmp))
        jobs=[(one,*case) for case in runs]+[(conditional,*case) for case in conditioned]
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            rows=list(pool.map(evaluate,jobs))
        report=dict(runs=rows[:len(runs)],conditioned=rows[len(runs):],
                    control_coefficients=FIT,source_control_emitted_energy_over_reference={
                    label:sum(x*x for x in v) for label,(_,v) in FIT.items()})
        if args.refit:report['repeated_fits']=refit(inp,exe,tmp)
    if any(abs(x['mean_relative_density']-1.)>.01 or x['max_endpoint_error']>2.e-10
           for x in report['conditioned']):raise RuntimeError('Local volume or endpoint control failed')
    for x in report['runs']:
        if x['case']=='uniform_local_pressure_match':
            cov=x['metrics']['per_sign_pooled_lateral_covariance_over_target']
            if np.max(abs(np.linalg.eigvalsh(cov)-1.))>.025:raise RuntimeError('Pressure calibration did not transfer')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2,allow_nan=False)+'\n')


if __name__=='__main__': main()
