"""Independent geometry cache and direct-F trajectory; imported without evaluation."""
from __future__ import annotations
import time,resource
import numpy as np
from kernel import LD,CD,DenseKernel,baseline,source_jet,forcing
from moment_contacts import weighted_subtractions


def cut_prefix(values,counts):
    return np.array([np.sum(values[:n],axis=0,dtype=LD) for n in counts],dtype=LD)


def geometry(k,weights,times,counts,pi,check=lambda:None):
    times=np.asarray(times,dtype=LD)
    out={key:np.zeros((len(times),3),dtype=LD) for key in ('r0','p0')}
    out['moments']=np.zeros((8,len(times),3),dtype=LD)
    measure_weights=weights*k*k/(2*pi*pi)
    out['inverse_k']=cut_prefix(measure_weights/k,counts)
    # Bound every temporary: 256 momenta x32 time points.
    for tstart in range(0,len(times),32):
        check()
        ids=slice(tstart,tstart+32)
        L=-LD(1)/times[ids][None,:]
        for mstart in range(0,len(k),256):
            mstop=min(mstart+256,len(k))
            mb=slice(mstart,mstop)
            kk=k[mb,None]
            omega=np.sqrt(kk*kk+2*L*L)
            iw=1/omega; iw2=iw*iw
            r0,p0=baseline(kk,L)
            weighted=measure_weights[mb,None]
            powers=[]
            power=iw
            for n in range(8):
                powers.append(power*weighted)
                power=power*iw2
            for j,c in enumerate(counts):
                count=max(0,min(mstop,c)-mstart)
                if count==0:continue
                out['r0'][ids,j]+=np.sum((r0*weighted)[:count],axis=0,dtype=LD)
                out['p0'][ids,j]+=np.sum((p0*weighted)[:count],axis=0,dtype=LD)
                for n,p in enumerate(powers):out['moments'][n,ids,j]+=np.sum(p[:count],axis=0,dtype=LD)
    return out


def contacts(times,jet,geo):
    L=-LD(1)/times[:,None]
    H=jet[:,:,None]
    dr,dp=weighted_subtractions(L,geo['moments'],H)
    r=geo['inverse_k']*(L*L*H[0]+L*H[1]/2)-dr-4*H[0]*geo['r0']
    p=geo['inverse_k']*(-L*L*H[0]+L*H[1]/2)-dp-4*H[0]*geo['p0']
    work=-3*H[1]*(geo['r0']+geo['p0'])
    f=L*(r-3*p)+work
    return {'R':r,'P':p,'F':f,'work':work}


def prepare(k,weights,steps,order,counts,pi,a,b,check=lambda:None):
    dt=(b-a)/steps
    # This creates fixed quadrature constants only. No source is sampled here.
    reference=DenseKernel(k[:1],dt,order,order)
    grid=a+np.arange(steps+1,dtype=LD)*dt
    dense=(grid[:-1,None]+reference.t[None,:]).ravel()
    return {'dt':dt,'grid':grid,'dense':dense,'q':order,'source_q':order,'t':reference.t,'tw':reference.tw,
            's':reference.s,'sw':reference.sw,
            'grid_geometry':geometry(k,weights,grid,counts,pi,check),
            'dense_geometry':geometry(k,weights,dense,counts,pi,check)}


def physical_sources(prep,source):
    grid=prep['grid'];dt=prep['dt'];t=prep['t'];s=prep['s']
    dense_jet=source_jet(prep['dense'],source)
    grid_jet=source_jet(grid,source)
    nested=grid[:-1,None,None]+t[None,:,None]*s[None,None,:]
    end=grid[:-1,None]+dt*s[None,:]
    return {'grid':grid_jet,'dense':dense_jet,
            'nested_g':forcing(nested,source_jet(nested,source)),
            'end_g':forcing(end,source_jet(end,source))}


def trajectory(k,weights,u0,w0,prep,source_arrays,epsilon,counts,pi,check=lambda:None):
    steps=len(prep['grid'])-1;q=prep['q']
    output={key:np.zeros((steps+1,3),dtype=LD) for key in ('R','P','F')}
    cg=contacts(prep['grid'],source_arrays['grid'],prep['grid_geometry'])
    cd=contacts(prep['dense'],source_arrays['dense'],prep['dense_geometry'])
    contact_I=np.sum(cd['F'].reshape(steps,q,3)*prep['tw'][None,:,None],axis=(0,1),dtype=LD)
    work_I=np.sum(cd['work'].reshape(steps,q,3)*prep['tw'][None,:,None],axis=(0,1),dtype=LD)
    modal_I=np.zeros(len(k),dtype=LD)
    ends={'u_mid':np.empty_like(u0),'w_mid':np.empty_like(w0),'u_end':np.empty_like(u0),'w_end':np.empty_like(w0)}
    for mstart in range(0,len(k),256):
        check();mstop=min(mstart+256,len(k));mb=slice(mstart,mstop)
        kk=k[mb];mw=weights[mb]*kk*kk/(2*pi*pi)
        kernel=DenseKernel(kk,prep['dt'],prep['source_q'],q)
        u=u0[mb].copy();w=w0[mb].copy()
        for step in range(steps+1):
            if step%32==0:check()
            L=-LD(1)/prep['grid'][step]
            U,W=u/epsilon,w/epsilon
            common=-L*W.real-kk*W.imag
            r=((2*kk*kk+3*L*L)*U.real+common)/(2*kk)
            p=((2*kk*kk/3-L*L)*U.real+common)/(2*kk)
            f=(3*L*L*L*U.real+L*L*W.real)/kk+L*W.imag
            for j,c in enumerate(counts):
                n=max(0,min(mstop,c)-mstart)
                if n:
                    output['R'][step,j]+=np.sum((mw*r)[:n],dtype=LD)
                    output['P'][step,j]+=np.sum((mw*p)[:n],dtype=LD)
                    output['F'][step,j]+=np.sum((mw*f)[:n],dtype=LD)
            if step==steps//2:
                ends['u_mid'][mb]=u;ends['w_mid'][mb]=w
            if step==steps:
                ends['u_end'][mb]=u;ends['w_end'][mb]=w
                break
            du,dw=kernel.dense(u,w,source_arrays['nested_g'][step],epsilon)
            Ld=-LD(1)/(prep['grid'][step]+kernel.t)
            fd=(3*Ld[None,:]**3*(du.real/epsilon)+Ld[None,:]**2*(dw.real/epsilon))/kk[:,None]+Ld[None,:]*(dw.imag/epsilon)
            modal_I[mb]+=np.sum(fd*kernel.tw[None,:],axis=1,dtype=LD)
            u,w=kernel.step(u,w,source_arrays['end_g'][step],epsilon)
    for key in ('R','P','F'):output[key]+=cg[key]
    gd=forcing(prep['dense'],source_arrays['dense'])
    J_Lg=np.sum((-LD(1)/prep['dense']).reshape(steps,q)*gd.reshape(steps,q)*prep['tw'][None,:],dtype=LD)
    output.update({'J_Lg':J_Lg,'modal_I':modal_I,'contact_I':contact_I,'work_I':work_I,
                   'contact_R':cg['R'],'contact_P':cg['P'],'contact_F':cg['F'],**ends})
    return output
