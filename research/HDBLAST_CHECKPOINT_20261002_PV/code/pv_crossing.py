#!/usr/bin/env python3
"""Matched flat-space quantum-source benchmark; no shell or bulk backreaction."""
import argparse
import hashlib
import json
import gc
from pathlib import Path
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import solve_ivp

M,T,R=2.0,1.0,4.0
C=np.array([1.,-3.,3.,-1.])
JS=np.arange(4.)
COMMON=np.linspace(-6.,6.,2001)
LIMITS=dict(wronskian=1e-7,occupation_abs=1e-6,ward=1e-5,
            numeric_source_scaled=1e-5,regulator_source_scaled=1e-3,trace_scaled=1e-4)
CASES=[
 ("primary",16.,128.,24,10.,1e-10,1e-13),
 ("lambda4",4.,32.,24,10.,1e-10,1e-13),
 ("lambda8",8.,64.,24,10.,1e-10,1e-13),
 ("quadrature48",16.,128.,48,10.,1e-10,1e-13),
 ("cutoff64",16.,64.,24,10.,1e-10,1e-13),
 ("cutoff192",16.,192.,24,10.,1e-10,1e-13),
 ("tight",16.,128.,24,10.,1e-12,1e-15),
 ("prehistory8",16.,128.,24,8.,1e-10,1e-13)]
def clean(v):
    if isinstance(v,dict): return {str(k):clean(x) for k,x in v.items()}
    if isinstance(v,(tuple,list)): return [clean(x) for x in v]
    if isinstance(v,np.ndarray): return clean(v.tolist())
    if isinstance(v,(float,np.floating)): return float(v) if np.isfinite(v) else None
    if isinstance(v,(bool,np.bool_)): return bool(v)
    if isinstance(v,(int,np.integer)): return int(v)
    return v
def pulse(t):
    u=np.tanh(t/T)
    return R*u*u,2*R*u*(1-u*u)/T,2*R*(1-u*u)*(1-3*u*u)/T**2
def panels(K,n):
    edges=[0.,.125,.25,.5,1.,2.,4.,8.,16.,32.,64.,128.]
    edges=[z for z in edges if z<K]+[K]
    xx,ww=leggauss(n)
    k=np.concatenate([(a+b)/2+(b-a)*xx/2 for a,b in zip(edges[:-1],edges[1:])])
    wt=np.concatenate([(b-a)*ww/2 for a,b in zip(edges[:-1],edges[1:])])
    return k,wt*k*k/(2*np.pi**2),edges
def jost(t,k,lam):
    """Small-z exact-incoming series, with irrelevant global phase omitted."""
    om_inf=np.sqrt(k[None,:]**2+R+JS[:,None]*lam**2)
    z=.5*(1+np.tanh(t/T))
    if not (0<z<1e-5): raise ValueError("Jost initial slice must have small positive z")
    cc=1-1j*om_inf*T
    term=np.ones_like(cc); f=term.copy();fz=np.zeros_like(cc)
    for n in range(1,40):
        term=term*(R*T*T+n*(n-1))*z/(n*(cc+n-1))
        f+=term;fz+=n*term/z
        if np.max(np.abs(term))<1e-30: break
    else: raise RuntimeError("Jost series did not converge")
    amp=f/np.sqrt(2*om_inf)
    vel=(-1j*om_inf*f+fz*2*z*(1-z)/T)/np.sqrt(2*om_inf)
    return amp,vel,f,fz,om_inf,z,n
def phase(t,k,lam):
    b=k[None,:]**2+JS[:,None]*lam**2
    a=np.sqrt(b+R);h=t/T
    return T*(a*np.arcsinh(a*np.sinh(h)/np.sqrt(b))
              -M*np.arcsinh(M*np.tanh(h)/np.sqrt(b)))
def init(t,k,lam):
    _,_,f,fz,oinf,z,terms=jost(t,k,lam)
    x,_,_=pulse(t)
    om=np.sqrt(k[None,:]**2+x+JS[:,None]*lam**2)
    dz=2*z*(1-z)/T;omdiff=(x-R)/(om+oinf)
    ai=((om+oinf)*f+1j*fz*dz)/(2*np.sqrt(om*oinf))
    bi=(omdiff*f-1j*fz*dz)/(2*np.sqrt(om*oinf))
    th=phase(t,k,lam)
    return ai*np.exp(1j*th),bi*np.exp(-1j*th),terms
def exact_occupation(k,lam):
    p=T*np.sqrt(k[None,:]**2+R+JS[:,None]*lam**2)
    z=np.pi*np.sqrt(R*T*T-.25)
    logcosh=z+np.log1p(np.exp(-2*z))-np.log(2.)
    logsinh=np.pi*p+np.log1p(-np.exp(-2*np.pi*p))-np.log(2.)
    return np.exp(2*(logcosh-logsinh))
def matching_constants(lam):
    rr=R+JS*lam**2
    vp=float(np.dot(C,rr*np.log(rr))/(32*np.pi**2))
    vpp=float(np.dot(C,np.log(rr))/(32*np.pi**2))
    return vp,vpp,vpp/6
def matched_potential(x,lam):
    rr=R+JS*lam**2;u=(x-R)/rr
    h=np.empty(4);hp=np.empty(4)
    for j,d in enumerate(u):
        if d == -1:
            h[j]=-.5;hp[j]=2.
        elif abs(d)<.125:
            h[j]=sum(2*(-1)**(n-1)*d**n/(n*(n-1)*(n-2)) for n in range(3,28))
            hp[j]=sum(2*(-1)**n*d**n/(n*(n-1)) for n in range(2,28))
        else:
            h[j]=(1+d)**2*np.log1p(d)-d-1.5*d*d
            hp[j]=2*((1+d)*np.log1p(d)-d)
    return float(np.dot(C,rr*rr*h)/(64*np.pi**2)),float(np.dot(C,rr*hp)/(64*np.pi**2))
def static_tails(x,lam,K):
    y=x+JS*lam**2;u=y/K**2
    root=np.sqrt(1+u);v=u/(root+1);ell=np.log1p(v/2)
    et=K**4*np.dot(C,u*u*ell+.25*v**4)/(32*np.pi**2)
    qt=K**2*np.dot(C,u*ell-v+.5*u-.375*u*u)/(8*np.pi**2)
    bound=K**4*np.dot(C,v**3*(4+v))/(96*np.pi**2)
    return float(et),float(qt),float(bound)
def trace_tail(x,xd,xdd,lam,K):
    y=x+JS*lam**2;u=y/K**2;root=np.sqrt(1+u)
    g=np.log1p(u/(2*(root+1)))-np.expm1(-.5*np.log1p(u))
    b1=-np.dot(C,g)/(8*np.pi**2)
    rat=np.empty(4);nonzero=y!=0
    rat[nonzero]=np.expm1(-1.5*np.log1p(u[nonzero]))/y[nonzero]
    rat[~nonzero]=-1.5/K**2
    b2=np.dot(C,rat)/(16*np.pi**2)
    return float(.5*(b2*xd*xd+b1*xdd))
def source(t,aa,bb,k,wt,lam,K,const,ref_tails):
    x,xd,xdd=pulse(t)
    y=x+JS[:,None]*lam**2;om=np.sqrt(k[None,:]**2+y)
    oi=np.sqrt(k[None,:]**2+R+JS[:,None]*lam**2)
    th=phase(t,k,lam);nn=np.abs(bb)**2
    coh=np.real(aa*np.conj(bb)*np.exp(-2j*th))
    invdiff=(R-x)/(om*oi*(om+oi))
    qr=(nn+coh)/om+.5*invdiff
    er=om*nn+(x-R)/(2*(om+oi))
    pr=k[None,:]**2*nn/(3*om)-(2*k[None,:]**2/3+y)*coh/om
    pr+=k[None,:]**2*invdiff/6
    qref=qr@wt;eref=er@wt;pref=pr@wt
    et,qt,bo=static_tails(x,lam,K);e0,q0,b0=ref_tails
    vp,vpp,fp=const;d=x-R
    rho=float(np.dot(C,eref)-vp*d-.5*vpp*d*d+et-e0)
    Q=float(np.dot(C,qref)-2*vpp*d+qt-q0)
    pres=float(np.dot(C,pref)+vp*d+.5*vpp*d*d-(et-e0)-(bo-b0)-2*fp*xdd)
    V,Vx=matched_potential(x,lam)
    q0K=[];q0R=[]
    for z in x+JS*lam**2:
        q0K.append(K*K/(8*np.pi**2) if z==0 else
                   (K*np.sqrt(K*K+z)-z*np.arcsinh(K/np.sqrt(z)))/(8*np.pi**2))
    for z in R+JS*lam**2:
        q0R.append((K*np.sqrt(K*K+z)-z*np.arcsinh(K/np.sqrt(z)))/(8*np.pi**2))
    qdyn=qref-(np.array(q0K)-np.array(q0R))
    anomaly=float(-4*V+2*x*Vx-lam**2*np.dot(C*JS,qdyn))
    return rho,pres,Q/2,anomaly,trace_tail(x,xd,xdd,lam,K),nn

def direct_limit_source(t,aa,bb,k,wt):
    """Same matching convention, physical modes only; finite-K convergent integrands."""
    x,xd,xdd=pulse(t)
    om=np.sqrt(k*k+x);oi=np.sqrt(k*k+R);d=x-R
    th=phase(t,k,0.)[0];a=aa[0];b=bb[0]
    ap=a*np.exp(-1j*th);bp=b*np.exp(1j*th)
    n=np.abs(b)**2;co=np.real(a*np.conj(b)*np.exp(-2j*th))
    wr=np.abs(a)**2-n-1
    e_static=d**3*(om+3*oi)/(16*oi**3*(om+oi)**3)
    rho=float(np.dot(wt,om*n+e_static))
    qraw=np.abs(ap+bp)**2/(2*om)
    qren=qraw-wr/(2*om)-1/(2*oi)+d/(4*oi**3)
    S=float(np.dot(wt,qren)/2)
    v=d/(oi*(om+oi))
    p_static=-k*k*v**3*(20+15*v+3*v*v)/(48*om)
    p_dynamic=k*k*n/(3*om)-(2*k*k/3+x)*co/om
    p_kernel=xdd*(2*k*k/3+R)/(16*oi**5)
    p=float(np.dot(wt,p_dynamic+p_static+p_kernel))
    return rho,p,S

def second(v,h,stencil):
    if stencil==5:
        cc=np.array([-1.,16.,-30.,16.,-1.])/12;width=2
    elif stencil==9:
        cc=np.array([-1/560,8/315,-1/5,8/5,-205/72,8/5,-1/5,8/315,-1/560]);width=4
    else:raise ValueError(stencil)
    out=np.full(len(v),np.nan)
    out[width:-width]=sum(c*v[width+i:len(v)-width+i] for i,c in zip(range(-width,width+1),cc))/h**2
    return out
def trace_check(curves):
    h=COMMON[1]-COMMON[0];Q=2*curves[:,4];tar=-curves[:,2]+3*curves[:,3]
    x=np.array([pulse(t)[0] for t in COMMON]);base=-x*Q+curves[:,5]+curves[:,6]
    d5=.5*second(Q,h,5)+base;d9=.5*second(Q,h,9)+base
    coarse=.5*second(Q[::2],2*h,9)+base[::2]
    pick=np.arange(len(Q))[::2];mask=(pick>=8)&(pick<len(Q)-8)
    err5=np.nanmax(np.abs(d5[8:-8]-tar[8:-8]))/M**4
    err9=np.nanmax(np.abs(d9[8:-8]-tar[8:-8]))/M**4
    errc=np.nanmax(np.abs(coarse[mask]-tar[pick[mask]]))/M**4
    refin=np.nanmax(np.abs(d9[pick[mask]]-coarse[mask]))/M**4
    return dict(five_point_max_scaled=float(err5),nine_point_max_scaled=float(err9),
                coarse_nine_point_max_scaled=float(errc),refinement_max_scaled=float(refin),
                threshold=LIMITS["trace_scaled"],pass_nine=bool(err9<=LIMITS["trace_scaled"]),
                finite_K_correction="D_K included; 8 grid indices excluded at boundaries")
def run(spec,out):
    name,lam,K,n,end,rtol,atol=spec
    k,wt,edges=panels(K,n);aa,bb,terms=init(-end,k,lam);nm=aa.size
    const=matching_constants(lam);tails=static_tails(R,lam,K)
    y0=np.concatenate([aa.ravel(),bb.ravel(),[0j]])
    r0=source(-end,aa,bb,k,wt,lam,K,const,tails)[0]
    max_step=.25/np.sqrt(K*K+R+3*lam*lam)
    def rhs(t,y):
        aa=y[:nm].reshape(4,-1);bb=y[nm:2*nm].reshape(4,-1)
        x,xd,xdd=pulse(t);om=np.sqrt(k[None,:]**2+x+JS[:,None]*lam**2)
        h=xd/(4*om*om);e=np.exp(2j*phase(t,k,lam))
        da=h*bb*e;db=h*aa*np.conj(e)
        oi=np.sqrt(k[None,:]**2+R+JS[:,None]*lam**2)
        nn=np.abs(bb)**2;co=np.real(aa*np.conj(bb)*np.conj(e))
        qref=((nn+co)/om+(R-x)/(2*om*oi*(om+oi)))@wt
        qt=static_tails(x,lam,K)[1]
        Q=float(np.dot(C,qref)-2*const[1]*(x-R)+qt-tails[1])
        return np.concatenate([da.ravel(),db.ravel(),[complex(xd*Q/2)]])
    grid=np.concatenate([[-end],COMMON,[end]])
    sol=solve_ivp(rhs,(-end,end),y0,method="DOP853",t_eval=grid,
                  dense_output=False,rtol=rtol,atol=atol,max_step=max_step)
    if not sol.success or len(sol.t)!=len(grid):raise RuntimeError(sol.message)
    arr=np.empty((len(grid),16));wr=0.
    for i,t in enumerate(grid):
        aa=sol.y[:nm,i].reshape(4,-1);bb=sol.y[nm:2*nm,i].reshape(4,-1)
        wr=max(wr,float(np.max(np.abs(np.abs(aa)**2-np.abs(bb)**2-1))))
        rho,p,S,A,tail,_=source(t,aa,bb,k,wt,lam,K,const,tails)
        work=sol.y[-1,i].real;balance=rho-r0-work
        norm=abs(rho)+abs(r0)+abs(work)+M**4
        dr,dp,dS=direct_limit_source(t,aa,bb,k,wt)
        x,xd,xdd=pulse(t)
        om=np.sqrt(k[None,:]**2+x+JS[:,None]*lam**2)
        wd=np.abs(np.abs(aa)**2-np.abs(bb)**2-1)
        absw=np.abs(C)[:,None]*wd
        nb_rho=float(np.sum((absw*om/2)@wt))
        nb_p=float(np.sum((absw*k[None,:]**2/(6*om))@wt))
        nb_S=float(np.sum((absw/(4*om))@wt))
        p_without_curvature=p+2*const[2]*xdd
        arr[i]=[t,x,rho,p,S,A,tail,work,abs(balance)/norm,
                dr,dp,dS,p_without_curvature,nb_rho,nb_p,nb_S]
    aa=sol.y[:nm,-1].reshape(4,-1);bb=sol.y[nm:2*nm,-1].reshape(4,-1)
    th=phase(end,k,lam);om=np.sqrt(k[None,:]**2+pulse(end)[0]+JS[:,None]*lam**2)
    ap=aa*np.exp(-1j*th);bp=bb*np.exp(1j*th)
    ff=(ap+bp)/np.sqrt(2*om);fd=-1j*np.sqrt(om/2)*(ap-bp)
    jout,joutd,*_=jost(-end,k,lam);jout=np.conj(jout);joutd=-np.conj(joutd)
    bet=-1j*(jout*fd-joutd*ff);nout=np.abs(bet)**2;target=exact_occupation(k,lam)
    nerr=float(np.max(np.abs(nout-target)));common=arr[1:-1];tr=trace_check(common)
    wi=float(np.max(arr[:,8]))
    controls=dict(wronskian=bool(wr<=LIMITS["wronskian"]),
                  exact_occupation=bool(nerr<=LIMITS["occupation_abs"]),
                  ward=bool(wi<=LIMITS["ward"]),pressure_trace=bool(tr["pass_nine"]))
    np.savez_compressed(out/(name+"_curves.npz"),curves=common,k=k,weights=wt,n_out=nout,n_exact=target)
    report=dict(name=name,lambda_=lam,K=K,nodes_per_panel=n,nk=len(k),edges=edges,
                t_initial=-end,t_final=end,rtol=rtol,atol=atol,max_step=max_step,nfev=sol.nfev,
                jost_series_terms=terms,wronskian_max=wr,exact_occupation_abs=nerr,
                ward_path_max=wi,ward_endpoint=float(arr[-1,8]),
                ward_absolute_over_m4=float(np.max(np.abs(arr[:,2]-r0-sol.y[-1].real))/M**4),
                trace=tr,controls=controls,
                crossing=dict(rho=float(common[1000,2]),p=float(common[1000,3]),
                              S=float(common[1000,4]),p_potential_only=float(common[1000,12])),
                direct_limit_difference_scaled={key:float(np.max(np.abs(common[:,c]-common[:,d]))/scale)
                    for key,c,d,scale in [("rho",2,9,M**4),("p",3,10,M**4),("S",4,11,M**2)]},
                normalization_reconstruction_absolute_bounds={
                    key:float(np.max(arr[:,col])) for key,col in [("rho",13),("p",14),("S",15)]},
                endpoint_rho=float(arr[-1,2]),
                endpoint_pressure=float(arr[-1,3]),endpoint_S=float(arr[-1,4]),
                common_peak_abs={key:float(np.max(np.abs(common[:,col]))) for key,col in
                                 [("rho",2),("p",3),("S",4),("anomaly",5)]},
                spectrum=dict(k=k,n_physical=nout[0],exact_physical=target[0]),
                matching=dict(x_ref=R,V_value=0,V_x=0,V_xx=0,F_value=0,F_x=0),
                status="finite PV regulator and finite dynamical momentum cutoff",
                ledger="rho-rho_initial-integral xdot*S dt; normalizer |rho|+|rho_initial|+|work|+m_inf^4")
    del sol;gc.collect()
    return report,common
def comparison(a,b):
    keys=[("rho",2,M**4),("p",3,M**4),("S",4,M**2)]
    d={key:float(np.max(np.abs(a[:,col]-b[:,col]))/scale) for key,col,scale in keys}
    d["numeric_target_pass"]=all(d[key]<=LIMITS["numeric_source_scaled"] for key,_,_ in keys)
    d["regulator_target_pass"]=all(d[key]<=LIMITS["regulator_source_scaled"] for key,_,_ in keys)
    return d
def null_control():
    residual=[]
    for lam in [4.,8.,16.]:
        for K in [8*lam,12*lam]:
            v,vx=matched_potential(R,lam);et,qt,bo=static_tails(R,lam,K)
            residual.append(abs(v)+abs(vx)+abs(et-et)+abs(qt-qt)+abs(bo-bo))
    return dict(max_source_absolute=float(max(residual)),pass_=bool(max(residual)<=1e-12),
                scope="Static algebraic prescription null; no nontrivial ODE")
def main():
    pa=argparse.ArgumentParser()
    pa.add_argument("--out",type=Path,default=Path("outputs/pv_actual"))
    pa.add_argument("--registration",type=Path,required=True)
    pa.add_argument("--registration-sha256",required=True)
    pa.add_argument("--source-ref",default="unspecified")
    ar=pa.parse_args()
    regsha=hashlib.sha256(ar.registration.read_bytes()).hexdigest()
    if regsha!=ar.registration_sha256:raise ValueError("Registration hash mismatch")
    ar.out.mkdir(parents=True,exist_ok=True)
    report=dict(model="Smooth Minkowski external-source benchmark, x=4 tanh²(t)",
                source_ref=ar.source_ref,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                protocol_addendum_sha256=hashlib.sha256(Path("PROTOCOL_ADDENDUM.md").read_bytes()).hexdigest(),
                registration_sha256=regsha,registration_commit="32f1fc617b0691bbc10783e7e12b337afb4a3af9",
                limits=LIMITS,cases=[],errors={},
                source_units="rho,p comparison scale m_inf^4; S scale m_inf²",
                state="Exact incoming Jost state for all physical and regulator fields",
                prescription="PV(1,-3,3,-1); V,Vx,Vxx,F,Fx zero at x_ref=4",
                scope="Flat a=1 benchmark; not original shell/bulk evolution")
    curves={}
    for spec in CASES:
        try:
            rr,cv=run(spec,ar.out);report["cases"].append(rr);curves[spec[0]]=cv
            print("HDBLAST_PV_CASE="+json.dumps(clean(rr),separators=(",",":"),allow_nan=False),flush=True)
        except Exception as e:
            report["errors"][spec[0]]=type(e).__name__+": "+str(e)
            print("HDBLAST_PV_ERROR="+spec[0]+":"+str(e),flush=True)
    report["null"]=null_control()
    if "primary" in curves:
        prim=curves["primary"]
        report["comparisons"]={name:comparison(cv,prim) for name,cv in curves.items() if name!="primary"}
        report["regulator_successive"]={
            "4_to_8":comparison(curves["lambda4"],curves["lambda8"]) if "lambda4" in curves and "lambda8" in curves else None,
            "8_to_16":comparison(curves["lambda8"],prim) if "lambda8" in curves else None}
        inf=np.array([(pulse(t)[0]-R)**2/(32*np.pi**2)+pulse(t)[2]/(96*np.pi**2) for t in COMMON])
        report["anomaly_limit_max_scaled"]={name:float(np.max(np.abs(cv[:,5]-inf))/M**4) for name,cv in curves.items()}
        np.savetxt(ar.out/"primary_curves.csv",prim,delimiter=",",
                   header="t,x,rho,p,S,A_Lambda,D_K,scalar_work,ward_normalized,rho_direct,p_direct,S_direct,p_potential_only,norm_rho_bound,norm_p_bound,norm_S_bound",comments="")
        print("HDBLAST_PV_PRIMARY_CURVES="+json.dumps(clean(dict(columns=["t","x","rho","p","S","A_Lambda","D_K","scalar_work","ward_normalized","rho_direct","p_direct","S_direct","p_potential_only","norm_rho_bound","norm_p_bound","norm_S_bound"],rows=prim)),separators=(",",":"),allow_nan=False),flush=True)
    report["registered_controls_pass"]=(len(report["cases"])==len(CASES) and not report["errors"]
        and report["null"]["pass_"] and all(all(z["controls"].values()) for z in report["cases"]))
    report["source_convergence_certified"]=False
    report["convergence_note"]="Finite tested differences only; fixed K/Lambda does not prove joint regulator limit."
    (ar.out/"pv_actual.json").write_text(json.dumps(clean(report),indent=2,allow_nan=False)+"\n")
    print("HDBLAST_PV_ACTUAL_JSON="+json.dumps(clean(report),separators=(",",":"),allow_nan=False),flush=True)
    return int(not report["registered_controls_pass"])
if __name__=="__main__":
    raise SystemExit(main())
