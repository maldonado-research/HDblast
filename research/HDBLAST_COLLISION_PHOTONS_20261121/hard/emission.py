"""Massive charged-leg radiation from a neutral s-channel vector, in the pair CM."""
import argparse
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss


MASS = 8300.
WIDTH_FACTOR = .30
IDENTITY = np.eye(4, dtype=complex)
TWO = np.eye(2, dtype=complex)
ZERO = np.zeros((2, 2), dtype=complex)
SIGMA = np.array([[[0., 1.], [1., 0.]], [[0., -1j], [1j, 0.]], [[1., 0.], [0., -1.]]])
GAMMA = np.array([np.block([[TWO, ZERO], [ZERO, -TWO]])] +
                 [np.block([[ZERO, sig], [-sig, ZERO]]) for sig in SIGMA])
METRIC = np.array([1., -1., -1., -1.])
BASIS = np.eye(4)
POLARIZATIONS = BASIS[1:3]


def dot(a,b):
    return (a*METRIC*b).sum(axis=-1)


def slash(v):
    return np.einsum('...i,ijk->...jk',v*METRIC,GAMMA)


def spinor(p,particle):
    n = np.sqrt(p[...,0]+1)[...,None,None]
    large = n*TWO
    small = np.einsum('...i,ijk->...jk',p[...,1:],SIGMA)/n
    return np.concatenate((large,small) if particle else (small,large),axis=-2)


def adjoint(u):
    return np.swapaxes(u.conj(),-1,-2) @ GAMMA[0]


def scaled(coeff,matrix):
    return np.asarray(coeff)[...,None,None]*matrix


def prop(s):
    return 1./(s-MASS**2+1j*WIDTH_FACTOR*(s+2)*np.sqrt(1.-4./s))


def born_average(s):
    b2=1.-4./s
    return s*s*(3.-2.*b2+b2*b2/3.)*abs(prop(s))**2


def cutoff(s):
    return 1.-4./s


def kinematics(s,x,c,fi,angle=0.):
    """Photon along +z. c is inclination of the incoming or outgoing fermion in the pair rest frame."""
    w=np.sqrt(s)
    e=.5*w
    c=np.asarray(c)
    gam=e*x
    gam=np.asarray(gam)
    k=np.stack((gam,np.zeros_like(gam),np.zeros_like(gam),gam),axis=-1)
    sp=s*(1.-x)
    pair_energy=np.sqrt(sp) if fi else w
    beta=np.sqrt(1.-4./pair_energy**2)
    rest_e=.5*pair_energy
    along=rest_e*beta*c
    side=rest_e*beta*np.sqrt(np.maximum(0.,(1.-c)*(1.+c)))
    sn,cs=np.sin(angle),np.cos(angle)
    transverse_x=side*cs
    transverse_y=side*sn
    if fi:
        # Light-cone components avoid a nearly cancelling boost at high recoil.
        pm=e*(1.-beta*c);pp=e*(1.-x)*(1.+beta*c)
        qm=e*(1.+beta*c);qp=e*(1.-x)*(1.-beta*c)
        p=np.stack([.5*(pp+pm),transverse_x,transverse_y,.5*(pp-pm)],axis=-1)
        q=np.stack([.5*(qp+qm),-transverse_x,-transverse_y,.5*(qp-qm)],axis=-1)
    else:
        p=np.stack([np.full_like(c,e),transverse_x,transverse_y,along],axis=-1)
        q=p*np.array([1,-1,-1,-1])
    pk=gam*e*(1.-beta*c)
    qk=gam*e*(1.+beta*c)
    return p,q,k,pk,qk


def emit_current(p,q,k,pk,qk,eps,fi,vector,raw=False):
    """Dirac spin matrix after emitting along the indicated incoming or outgoing pair."""
    g=slash(vector)
    pol,ks=slash(eps),slash(k)
    if raw:
        if fi:
            operator=scaled(.5/pk,pol@(slash(p+k)+IDENTITY)@g)+scaled(.5/qk,g@(-slash(q+k)+IDENTITY)@pol)
        else:
            operator=scaled(-.5/pk,g@(slash(p-k)+IDENTITY)@pol)+scaled(-.5/qk,pol@(slash(k-q)+IDENTITY)@g)
    else:
        # Apply each external Dirac equation before evaluating collinear denominators.
        eikonal=dot(p,eps)/pk-dot(q,eps)/qk
        if fi:
            operator=scaled(eikonal,g)+scaled(.5/pk,pol@ks@g)-scaled(.5/qk,g@ks@pol)
        else:
            operator=scaled(-eikonal,g)+scaled(.5/pk,g@ks@pol)-scaled(.5/qk,pol@ks@g)
    if fi:
        return adjoint(spinor(p,True))@operator@spinor(q,False)
    return adjoint(spinor(q,False))@operator@spinor(p,True)


def pair_rate(s,x,y,wts,fi):
    """Exact angle average using the opposite pair's isotropic integrated spin tensor."""
    sp=s*(1.-x)
    b=np.sqrt(1.-4./(sp if fi else s))
    ymax=np.arctanh(b)
    rapidity=ymax*y
    cos=np.tanh(rapidity)/b
    jac=ymax/(b*np.cosh(rapidity)**2)
    p,q,k,pk,qk=kinematics(s,x,cos,fi)
    w=np.sqrt(s)
    gam=k[...,0]
    # Three unit spacelike polarizations transverse to the mediator. In the CM of
    # the final pair this is the usual rest-frame spatial sum.
    # q.J=0 allows eps_L=(-gamma,0,0,w-gamma)/sqrt(sp) to be
    # replaced by a purely spatial projector even in the radiative-return tail.
    longitudinal=BASIS[3] if fi else (np.sqrt(sp)/(w-gam))[...,None]*BASIS[3]
    squares=np.zeros(np.broadcast_shapes(np.shape(x),np.shape(cos)))
    for eps in POLARIZATIONS:
        for med in (*POLARIZATIONS,longitudinal):
            j=emit_current(p,q,k,pk,qk,eps,fi,med)
            squares+=np.sum(np.abs(j)**2,axis=(-1,-2))
    # Spin-sum final rest tensor = 4/3(s'+2); initial spin average = 1/3(s+2).
    opposite=(s+2.)/3. if fi else (sp+2.)/3.
    prefactor=s*x*np.sqrt(1.-4./sp)*opposite*abs(prop(s if fi else sp))**2/(4.*np.pi*np.sqrt(1.-4./s)*born_average(s))
    return np.sum(wts*jac*squares*prefactor,axis=-1)


def soft_coefficient(s):
    b=np.sqrt(1.-4./s)
    return 2./np.pi*((1.+b*b)*np.arctanh(b)/b-1.)


def spectra(s,lo,hi,unit,wts):
    """Integrate dx directly at zero, logarithmically in x at low energy and in s' at hard recoil."""
    endpoint=cutoff(s)
    near_pole=np.zeros(0,dtype=bool)
    if hi<=.25:
        if lo==0.:
            x=endpoint*hi*unit
            dx=endpoint*hi*wts
        else:
            span=np.log(hi/lo)
            x=endpoint*lo*np.exp(span*unit)
            dx=x*span*wts
    else:
        # sqrt(log(s'/4)) resolves beta' at the exact massive threshold.
        # Break in log(s') on both shoulders and center of the running-width pole.
        qmin=4. if hi==1. else s*(1.-hi*endpoint)
        qmax=s*(1.-lo*endpoint)
        shoulders=MASS*MASS*np.r_[np.exp(np.array([-1.,-.3,0.,.3,1.])),.70,1.30]
        # A long recoil tail also needs resolution near its upper end even when
        # the pole lies outside the allowed range.
        segments=max(1,int(np.ceil(np.log(qmax/qmin)/3.)))
        internal=np.geomspace(qmin,qmax,segments+1)[1:-1]
        breaks=np.unique(np.r_[qmin,shoulders[(shoulders>qmin)&(shoulders<qmax)],internal,qmax])
        eta=np.sqrt(np.log(breaks/4.))
        xx=[];ww=[];in_pole=[]
        for a,b in zip(eta[:-1],eta[1:]):
            y=a+(b-a)*unit
            q_over_s=4.*np.exp(y*y)/s
            xx.append(1.-q_over_s)
            in_pole.append((q_over_s*s>=.70*MASS**2)&(q_over_s*s<=1.30*MASS**2))
            ww.append(2.*y*q_over_s*(b-a)*wts)
        x=np.concatenate(xx);dx=np.concatenate(ww);near_pole=np.concatenate(in_pole)
    rad=np.column_stack([pair_rate(s,x[:,None],unit[None,:],wts,f) for f in (False,True)])
    recoil_tag=float((dx*x/2.)[near_pole]@rad[near_pole,0]) if hi>.25 else 0.
    return dx@rad,(dx*x/2.)@rad,recoil_tag


def band(s,lo,hi,n):
    nodes,weights=leggauss(n)
    unit=(nodes+1.)/2.;wts=weights/2.
    prob,energy,recoil_tag=spectra(s,lo,hi,unit,wts)
    soft_energy=soft_coefficient(s)*cutoff(s)*(hi-lo)
    row=dict(band_over_endpoint=[lo,hi],mean_fraction_of_pair_energy_over_alpha=float(energy.sum()),
             soft_mean_fraction_over_alpha=float(soft_energy),isr_cm_energy_over_alpha=float(energy[0]),
             fsr_cm_energy_over_alpha=float(energy[1]))
    if hi>.25:row["isr_cm_energy_near_pole_over_alpha"]=recoil_tag
    if lo>0.:
        soft=2.*soft_coefficient(s)*np.log(hi/lo)
        row.update(isr_over_alpha=float(prob[0]),fsr_over_alpha=float(prob[1]),total_over_alpha=float(prob.sum()),
                   soft_over_alpha=float(soft),exact_over_soft=float(prob.sum()/soft))
    return row


def spin_checks(s,x,c,d,azimuth):
    """Test primitive Feynman numerators, spin and charge exchange before any angular averaging."""
    a,b,k,ak,bk=kinematics(s,x,c,False)
    l,r,_,lk,rk=kinematics(s,x,d,True,azimuth)
    a,b,l,r=a.reshape(4),b.reshape(4),l.reshape(4),r.reshape(4)
    ini=[adjoint(spinor(b,False))@GAMMA[i]@spinor(a,True) for i in range(4)]
    def amplitudes(final,eps,raw):
        l,r,lk,rk=final
        jj=[adjoint(spinor(l,True))@GAMMA[i]@spinor(r,False) for i in range(4)]
        mi=np.zeros((2,2,2,2),complex);mf=np.zeros_like(mi)
        for i in range(4):
            v=BASIS[i]*METRIC[i]  # gives gamma^i after slashing
            isr=emit_current(a,b,k,ak,bk,eps,False,v,raw)
            fsr=emit_current(l,r,k,lk,rk,eps,True,v,raw)
            outer=lambda left,right: left[...,None,None]*right
            mi+=METRIC[i]*outer(isr,jj[i])
            mf+=METRIC[i]*outer(ini[i],fsr)
        return mi*prop(s*(1.-x)),mf*prop(s)
    direct=(l,r,lk,rk);swapped=(r,l,rk,lk)
    ward=0.;mismatch=0.;terms=[]
    for final in (direct,swapped):
        a2=b2=inter=0.
        for eps in POLARIZATIONS:
            ii,ff=amplitudes(final,eps,False)
            di,df=amplitudes(final,eps,True)
            wi,wf=amplitudes(final,k/k[0],True)
            denom=np.hypot(np.linalg.norm(ii),np.linalg.norm(ff))
            ward=max(ward,np.linalg.norm(wi)/denom,np.linalg.norm(wf)/denom)
            mismatch=max(mismatch,np.linalg.norm(di-ii)/denom,np.linalg.norm(df-ff)/denom)
            a2+=(np.abs(ii)**2).sum()/4.;b2+=(np.abs(ff)**2).sum()/4.
            inter+=np.vdot(ii,ff).real/2.
        terms.append((a2,b2,inter))
    t,u=np.asarray(terms)
    closure=float(np.max(abs(a+b-l-r-k))/np.sqrt(s))
    shell=max(abs(dot(v,v)-1.) for v in (a,b,l,r))/s
    controls=dict(ward_raw_relative=float(ward),reduced_vs_raw_relative=float(mismatch),
                  charge_swapped_odd_interference_residual=float(abs(t[2]+u[2])/(t[0]+t[1])),
                  interference_fraction=float(t[2]/(t[0]+t[1])),charge_even_change=float(max(abs(t[:2]-u[:2]))/(t[0]+t[1])),
                  momentum_residual=closure,mass_shell_residual=float(shell))
    return controls


def born_check(s,c):
    a,b,_,_,_=kinematics(s,0.,1.,False);l,r,_,_,_=kinematics(s,0.,c,False)
    tree=np.zeros((2,2,2,2),complex)
    for i in range(4):
        ini=adjoint(spinor(b,False))@GAMMA[i]@spinor(a,True)
        fin=adjoint(spinor(l,True))@GAMMA[i]@spinor(r,False)
        tree+=METRIC[i]*ini[...,None,None]*fin
    b2=1-4./s
    analytic=s*s*(3.-2.*b2+b2*b2*c*c)
    spin=float((np.abs(tree)**2).sum()/4)
    beta=np.sqrt(b2);coupling=WIDTH_FACTOR*(s+2.)
    inherited=12.*np.pi*coupling*coupling/(s*((s-MASS**2)**2+b2*coupling**2))
    g2=12.*np.pi*WIDTH_FACTOR
    sigma=g2*g2*born_average(s)/(16.*np.pi*s)
    return max(abs(spin/analytic-1),abs(sigma/inherited-1))


def full_angular_control():
    # Low boosts allow a direct five-angle check independent of the isotropic tensor reduction.
    s,x=12.,.17
    nodes,w=leggauss(24)
    phi=(np.arange(12)+.27)*2*np.pi/12
    ci,cf,az=np.meshgrid(nodes,nodes,phi,indexing="ij")
    weights=np.broadcast_to((w[:,None]*w[None,:])[:,:,None]/(4*len(phi)),ci.shape).ravel()
    a,b,k,ak,bk=kinematics(s,x,ci.ravel(),False)
    l,r,_,lk,rk=kinematics(s,x,cf.ravel(),True,az.ravel())
    def outer(left,right):
        return left[...,None,None]*right[:,None,None,:,:]
    terms=np.zeros((3,len(weights)))
    for eps in POLARIZATIONS:
        ii=np.zeros((len(weights),2,2,2,2),complex)
        ff=np.zeros_like(ii)
        for mu in range(4):
            v=BASIS[mu]*METRIC[mu]
            j=adjoint(spinor(b,False))@GAMMA[mu]@spinor(a,True)
            f=adjoint(spinor(l,True))@GAMMA[mu]@spinor(r,False)
            iemit=emit_current(a,b,k,ak,bk,eps,False,v)
            femit=emit_current(l,r,k,lk,rk,eps,True,v)
            ii+=METRIC[mu]*outer(iemit,f)*prop(s*(1.-x))
            ff+=METRIC[mu]*outer(j,femit)*prop(s)
        dim=(1,2,3,4)
        terms[0]+=.25*np.sum(np.abs(ii)**2,axis=dim)
        terms[1]+=.25*np.sum(np.abs(ff)**2,axis=dim)
        terms[2]+=.5*np.sum((ii.conj()*ff).real,axis=dim)
    phase=s*x*np.sqrt(1.-4./(s*(1.-x)))/(4*np.pi*np.sqrt(1.-4./s)*born_average(s))
    avg=(terms@weights)*phase
    red=np.array([pair_rate(s,x,(nodes+1.)/2.,w/2.,f) for f in (False,True)])
    return dict(s=s,x=x,integrated_isr_fsr_interference_over_alpha=avg.tolist(),
                relative_component_difference=(avg[:2]/red-1.).tolist(),relative_interference=float(avg[2]/red.sum()))


BANDS=((0.,.001),(.001,.01),(.01,.10),(.10,.25),(.25,1.))


def rate_table(spec):
    s,n=spec
    return [band(s,lo,hi,n) for lo,hi in BANDS]


def table_values(row):
    # Integrate photon energy down to zero while keeping only infrared-finite counts.
    return np.array([[r.get("total_over_alpha",0.),r["mean_fraction_of_pair_energy_over_alpha"]] for r in row])


def interpolation(grid,values,s):
    ls=np.log(s)
    lo=np.clip(np.searchsorted(grid,ls)-2,0,len(grid)-4)
    result=np.zeros((len(s),*values.shape[1:]))
    for i in range(4):
        fac=np.ones_like(ls)
        for j in range(4):
            if i!=j:
                fac*=(ls-grid[lo+j])/(grid[lo+i]-grid[lo+j])
        result+=fac[:,None,None]*values[lo+i]
    if np.any(result<0.):raise ValueError("Negative interpolated band")
    return result


def measured(path,grid,table,quarter):
    with np.load(path) as h:
        s=h['s']
        alpha=3.*float(h['charge_to_mass'])**2*float(h['daughter_mass_L'])**2/(4*np.pi*float(h['action'])*float(h['K']))
        if float(h['channel_fraction'])!=1. or np.any(s<=4.):raise ValueError('Unexpected interaction or threshold')
        count=round(float(h['D'])/float(h['weight']))
        weight=float(h['weight'])*h['a']*h['E_in'].sum(axis=1)/float(h['gas'])
    runs=[]
    for skip in (2,1):
        v=interpolation(grid[::skip],table[::skip],s)
        energy=(alpha*weight)@v[:,:,1]
        occupancy=(alpha/count)*v[:,:,0].sum(axis=0)
        runs.append((energy,occupancy))
    fine,coarse=runs[1],runs[0]
    total=float(fine[0].sum())
    ls=np.log(s)
    # Piecewise linear interpolation of the nonnegative recoil tag across its acceptance threshold.
    returned=np.interp(ls,grid,quarter[:,2])
    selected=(s>1.30*MASS**2)
    return_fraction=float((alpha*weight[selected])@returned[selected]/total)
    soft=float(np.sum(weight*alpha*cutoff(s)*soft_coefficient(s)))
    hard_count=float(fine[1][2:].sum())
    def stats(band,number,energy):
        row=dict(band_over_endpoint=list(band),fraction_of_endpoint_gas=float(energy))
        if band[0]>0.:row["photons_per_charged_particle"]=float(number)
        return row
    return dict(tape=path.name,events=len(s),charged_particles=count,alpha=alpha,
                bands=[stats(b,n,e) for b,n,e in zip(BANDS,fine[1],fine[0])],
                hard_isr_return_to_near_pole_share_of_total_energy=return_fraction,
                total_energy_over_gas=total,energy_above_one_percent_over_gas=float(fine[0][2:].sum()),
                photons_above_one_percent_per_charged_particle=hard_count,
                soft_extrapolation_full_energy_over_gas=soft,exact_to_soft_extrapolated_total=total/soft,
                coarse_grid_relative_change_total_energy=float(coarse[0].sum()/total-1.),
                coarse_grid_relative_change_hard_count=float(coarse[1][2:].sum()/hard_count-1.))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--artifacts',type=Path)
    p.add_argument('--nodes',type=int,default=40)
    p.add_argument('--jobs',type=int,default=2)
    opt=p.parse_args()
    representatives=[5.,6000.,8300.,np.sqrt(8.63903e7),np.sqrt(4.81732e8),np.sqrt(3.61012e9),np.sqrt(1.46782e11)]
    result=[]
    for w in representatives:
        s=w*w
        unit,weights=leggauss(opt.nodes)
        def check(f):
            x=1.e-7*cutoff(s)
            return x*pair_rate(s,x,(unit+1)/2,weights/2,f)/soft_coefficient(s)
        fine=rate_table((s,opt.nodes))
        coarse=rate_table((s,opt.nodes//2))
        rec=dict(sqrt_s=float(w),bands=fine,half_order=coarse,soft_limit_ratio_ISR=float(check(False)),soft_limit_ratio_FSR=float(check(True)),
                 born_relative_residual=born_check(s,.37),amplitude_control=spin_checks(s,.13,.32,-.51,.95))
        if s>=MASS**2:
            x=.43;sp=s*(1.-x)
            approx=(1+(1.-x)**2)/x*(np.log(s)-1)/np.pi*(born_average(sp)/sp)/(born_average(s)/s)
            rec["massless_isr_factorization_ratio_at_x043"]=float(pair_rate(s,x,(unit+1)/2,weights/2,False)/approx)
        if (rec['born_relative_residual']>1.e-11 or rec['amplitude_control']['ward_raw_relative']>1.e-10
                or rec['amplitude_control']['charge_swapped_odd_interference_residual']>1.e-10
                or max(abs(check(f)-1.) for f in (False,True))>5.e-6):
            raise ValueError("Amplitude or soft limit failed")
        result.append(rec)
    out=dict(pole_mass=MASS,width_factor=WIDTH_FACTOR,quadrature_nodes=opt.nodes,runs=result,independent_angle_check=full_angular_control())
    if max(abs(v) for v in out['independent_angle_check']['relative_component_difference'])>1.e-9:
        raise ValueError("Inconsistent angular normalization")
    if opt.artifacts:
        paths=sorted(opt.artifacts.glob('*.npz'))
        if not paths:raise ValueError('No collision records found')
        # Use shared grid chosen by the complete tape range and resolve the pole neighborhood.
        bounds=[]
        for path in paths:
            with np.load(path) as z:bounds.append((np.min(z['s']),np.max(z['s'])))
        low=np.log(min(x[0] for x in bounds))-.002
        high=np.log(max(x[1] for x in bounds))+.002
        mid=2*np.log(MASS)
        grid=np.concatenate((np.linspace(low,mid-1.,81)[:-1],np.linspace(mid-1.,mid+1.,81)[:-1],np.linspace(mid+1.,high,97)))
        if not (low<mid-1.<mid+1.<high):raise ValueError('Expected archived tapes to straddle the vector pole')
        with ProcessPoolExecutor(opt.jobs) as pool:
            tables=list(pool.map(rate_table,[(np.exp(v),opt.nodes) for v in grid],chunksize=2))
        values=np.array([table_values(r) for r in tables])
        qkeys=("isr_cm_energy_over_alpha","fsr_cm_energy_over_alpha","isr_cm_energy_near_pole_over_alpha")
        quarter=np.asarray([[r[-1][key] for key in qkeys] for r in tables])
        out["energy_grid"]=dict(nodes=len(grid),min_s=float(np.exp(grid[0])),max_s=float(np.exp(grid[-1])),
            log_s=grid.tolist(),above_quarter_isr_energy_over_alpha=quarter[:,0].tolist(),
            above_quarter_fsr_energy_over_alpha=quarter[:,1].tolist(),above_quarter_near_pole_isr_energy_over_alpha=quarter[:,2].tolist())
        out["collision_histories"]=[measured(path,grid,values,quarter) for path in paths]
    opt.output.parent.mkdir(parents=True,exist_ok=True)
    opt.output.write_text(json.dumps(out,indent=2,sort_keys=True,allow_nan=False)+'\n')


if __name__=='__main__':
    main()
