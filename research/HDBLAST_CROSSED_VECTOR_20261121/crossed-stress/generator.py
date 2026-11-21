"""Instantaneous collision contribution to wall-frame stress for equal-mass unpolarized charges."""
import argparse
import json
from pathlib import Path
import runpy

import numpy as np
from scipy.stats import qmc


def scalar_shear(p, energy, axis):
    return (np.sum(p*axis,axis=-1)**2-np.sum(p*p,axis=-1)/3.)/energy


def integrate_pairs(first, second, random, nodes, weights, azimuths, same, kernel):
    # Last two random coordinates rotate the observer isotropically for the independent null.
    size=len(first)
    p=first[(random[:,0]*size).astype(int)]; q=second[(random[:,1]*size).astype(int)]
    rp=np.linalg.norm(p,axis=1);rq=np.linalg.norm(q,axis=1)
    ep=np.hypot(1.,rp);eq=np.hypot(1.,rq); total=ep+eq
    P=p+q
    npdir=np.divide(p,rp[:,None],out=np.zeros_like(p),where=rp[:,None]>0)
    nqdir=np.divide(q,rq[:,None],out=np.zeros_like(q),where=rq[:,None]>0)
    delta=(rp-rq)**2/(ep*eq+rp*rq+1.)+.5*rp*rq*np.sum((npdir-nqdir)**2,axis=1)
    s=4.+2*delta; mass=np.sqrt(s); kstar=np.sqrt(.5*delta)
    incoming=p-((ep+.5*mass)/(total+mass))[:,None]*P
    norm=np.linalg.norm(incoming,axis=1)
    a=np.divide(incoming,norm[:,None],out=np.zeros_like(p),where=norm[:,None]>0)
    a[norm==0]=[0.,0.,1.]
    ref=np.zeros_like(a); zaxis=np.array([0.,0.,1.])
    vertical=abs(a[:,2])<.8
    ref[vertical]=zaxis;ref[~vertical]=[1.,0.,0.]
    b=np.cross(a,ref);b/=np.linalg.norm(b,axis=1)[:,None]
    c=np.cross(a,b)
    h=2*random[:,3]-1.; ang=2*np.pi*random[:,4];rho=np.sqrt(1.-h*h)
    direction=np.column_stack((rho*np.cos(ang),rho*np.sin(ang),h))
    before=scalar_shear(p,ep,zaxis)+scalar_shear(q,eq,zaxis)
    null_before=scalar_shear(p,ep,direction)+scalar_shear(q,eq,direction)
    # Stratified azimuth with randomized phase; opposite azimuths cancel odd terms.
    phi=2*np.pi*(random[:,2,None]+np.arange(azimuths))/azimuths
    transverse=b[:,None,:]*np.cos(phi)[:,:,None]+c[:,None,:]*np.sin(phi)[:,:,None]
    r=(s-4)/(2*kernel['M']**2); log=np.log1p(r)
    jac0=np.divide(log,r,out=np.ones_like(r),where=r>0)
    flux=np.sqrt(delta*(delta+2.))/(ep*eq)
    prefactor=kernel['G2']**2*flux/(32*np.pi*s)
    output=np.zeros((len(s),2,3))  # pole/full, rate/wall shear/null shear
    max_conservation=0.
    for y,weight in zip(nodes,weights):
        x=np.divide(np.expm1(y*log),r,out=np.full_like(r,y),where=r>0)
        z=1.-x
        jac=jac0*np.exp(y*log)
        star=kstar[:,None,None]*(z[:,None,None]*a[:,None,:]+np.sqrt(x*(2.-x))[:,None,None]*transverse)
        blast=np.sum(star*P[:,None,:],axis=-1)/(mass*(total+mass))[:,None]
        out=star+(.5+blast[:,:,None])*P[:,None,:]
        other=P[:,None,:]-out
        eo=np.sqrt(1.+np.sum(out*out,axis=-1)); eother=np.sqrt(1.+np.sum(other*other,axis=-1))
        max_conservation=max(max_conservation,float(np.max(abs(eo+eother-total[:,None])/total[:,None])))
        diff=(scalar_shear(out,eo,zaxis)+scalar_shear(other,eother,zaxis)).mean(axis=1)-before
        null=(scalar_shear(out,eo,direction[:,None,:])+scalar_shear(other,eother,direction[:,None,:])).mean(axis=1)-null_before
        amp=kernel['terms'](s,z)+kernel['terms'](s,-z)
        # Summing over z and -z is valid for the pair-symmetric outgoing stress.
        pole=np.zeros_like(s) if same else amp[0]
        full=.5*amp[3:].sum(axis=0) if same else amp[:3].sum(axis=0)
        if np.any(full<0.): raise ValueError('Negative spin-averaged differential rate')
        for k,w in enumerate((pole,full)):
            dens=weight*jac*prefactor*w
            output[:,k,0]+=dens;output[:,k,1]+=dens*diff;output[:,k,2]+=dens*null
    if max_conservation>2e-11: raise ValueError('Outgoing energy conservation')
    # Separate the diagnostic region in which the stated helicity partial wave exceeds one half.
    high=mass/kernel['M']>=kernel['unitarity_cut']
    return np.stack((output[~high].sum(axis=0),output[high].sum(axis=0))),max_conservation


def measurement(path,history,kernel,nodes,weights,azimuths,power,batches):
    with np.load(path) as d:
        p=d['p'];sign=d['sign']
    plus=p[sign==1];minus=p[sign==-1]
    if len(plus)!=len(minus) or len(plus)+len(minus)!=len(p): raise ValueError('Unequal populations')
    e=np.sqrt(1.+np.sum(p*p,axis=1)); energy=e.mean()
    qfrac=scalar_shear(p,e,np.array([0.,0.,1.])).mean()/energy
    time=float(path.name[1:8])
    D=np.interp(time,history['t'],history['D']);a=np.interp(time,history['t'],history['a'])
    target_density=float(history['action']/(6*history['daughter_mass_L']**3)*D/a**3)
    batch_results=[]; peak=0.
    for n in range(batches):
        ens=[]
        for j,(p1,p2) in enumerate(((plus,minus),(plus,plus),(minus,minus))):
            u=qmc.Sobol(d=5,scramble=True,seed=93827+409*n+2719*j).random_base2(power)
            summed=np.zeros((2,2,3))
            for start in range(0,len(u),512):
                v,err=integrate_pairs(p1,p2,u[start:start+512],nodes,weights,azimuths,j!=0,kernel)
                summed+=v;peak=max(peak,err)
            ens.append(summed/len(u))
        batch_results.append(np.stack(ens))
    data=np.stack(batch_results)
    by_region=data[:,0]+.5*data[:,1]+.5*data[:,2]
    total=by_region.sum(axis=1)
    data=data.sum(axis=2)  # unlike + half of each sign's like-pair generator
    def stats(x):
        return dict(mean=float(x.mean()),batch_se=float(x.std(ddof=1)/np.sqrt(batches)),batches=x.tolist())
    result={}
    for i,key in enumerate(('isolated_pole','crossed_born')):
        result[key]=dict(encounters_per_particle_per_wall_time=stats(total[:,i,0]*target_density),
                         derivative_shear_over_energy_per_wall_time=stats(total[:,i,1]*target_density/(2*energy)),
                         rotated_null_derivative=stats(total[:,i,2]*target_density/(2*energy)))
    for j,region in enumerate(('below_cut','above_or_equal_cut')):
        result[region]=dict(pole_derivative=stats(by_region[:,j,0,1]*target_density/(2*energy)),
                            crossed_derivative=stats(by_region[:,j,1,1]*target_density/(2*energy)),
                            crossed_minus_pole_derivative=stats((by_region[:,j,1,1]-by_region[:,j,0,1])*target_density/(2*energy)))
    result['unlike_crossed_derivative']=stats(data[:,0,1,1]*target_density/(2*energy))
    result['like_crossed_derivative']=stats((data[:,1,1,1]+data[:,2,1,1])*target_density/(4*energy))
    result['full_over_pole_derivative']=stats(total[:,1,1]/total[:,0,1])
    return dict(snapshot=path.stem,initial_wall_shear_over_energy=float(qfrac),mean_gamma=float(energy),
                normalized_target_density=target_density,max_relative_energy_error=peak,**result)



def analytic_check(kernel):
    """For equal head-on beams in their CM, Delta Q = -Q_in (1-P2)."""
    x,w=np.polynomial.legendre.leggauss(12)
    sample=np.array([[.1,.2,.3,.4,.5]])
    for multiple in (1.,20.,100.):
        s=multiple*kernel['M']**2; r=np.sqrt(s/4.-1); e=np.sqrt(1.+r*r)
        p=np.array([[0.,0.,r]]); flux=2.*r/e; qin=4.*r*r/(3*e)
        reference=kernel['summary'](s)
        for same in (False,True):
            integrated,_=integrate_pairs(p,-p,sample,(x+1.)/2,w/2,4,same,kernel)
            measured=integrated.sum(axis=0)
            labels=((1,'like'),) if same else ((0,'pole'),(1,'unlike'))
            for i,key in labels:
                rate=flux*reference['total'][key]
                stress=-flux*qin*reference['quadrupole'][key]
                if not np.allclose(measured[i,:2],(rate,stress),rtol=1.e-7,atol=0.):
                    raise ValueError('Head-on invariant-frame rate and stress control')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--snapshots',type=Path,nargs='+',required=True)
    ap.add_argument('--kernel',type=Path,required=True)
    ap.add_argument('--history',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--power',type=int,default=13)
    ap.add_argument('--batches',type=int,default=6)
    ap.add_argument('--polar',type=int,default=10)
    ap.add_argument('--azimuths',type=int,default=8)
    ap.add_argument('--cut',type=float,default=np.sqrt(140379347.50895345)/8300.)
    opt=ap.parse_args()
    if opt.batches<2 or opt.azimuths<2 or opt.azimuths%2:raise ValueError('Need independent batches and paired azimuths')
    kernel=runpy.run_path(str(opt.kernel))
    kernel["unitarity_cut"]=opt.cut
    analytic_check(kernel)
    with np.load(opt.history) as hist: history=dict(hist)
    x,w=np.polynomial.legendre.leggauss(opt.polar)
    rows=[measurement(path,history,kernel,(x+1)/2,w/2,opt.azimuths,opt.power,opt.batches) for path in opt.snapshots]
    opt.output.write_text(json.dumps(dict(samples_per_sign_pair_per_batch=2**opt.power,batches=opt.batches,
            polar_order=opt.polar,azimuth_order=opt.azimuths,root_s_over_M_cut=opt.cut,runs=rows),indent=2,allow_nan=False)+'\n')
