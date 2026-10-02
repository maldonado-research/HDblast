"""Regularity audit of archived FRW reconstructions; no absolute-source claim."""
import argparse,hashlib,importlib.util,json
from pathlib import Path
import numpy as np
from scipy.interpolate import PPoly,make_interp_spline
from scipy.integrate import quad
from scipy.special import sici
from scipy.optimize import newton
import mpmath as mp
END=6.9
GRID=np.linspace(0.,END,2001)
LIMITS=dict(continuity_scaled=1e-9,step_asymptotic=1e-3,interference_log_slope=1e-3)
KINDS=("hermite","cubic","quintic_C4","septic_C6")
DEGREES=dict(hermite=3,cubic=3,quintic_C4=5,septic_C6=7)
CONTINUITY=dict(hermite=1,cubic=2,quintic_C4=4,septic_C6=6)
COLS=["t","ln_a","phi","a","H","phi_dot","Hdot","Hddot","Hdddot","U","R","Rdot","Rddot","rho_R2_regular","p_R2_regular"]
def clean(x):
    if isinstance(x,dict):return {str(k):clean(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [clean(v) for v in x]
    if isinstance(x,np.ndarray):return clean(x.tolist())
    if isinstance(x,(bool,np.bool_)):return bool(x)
    if isinstance(x,(int,np.integer)):return int(x)
    if isinstance(x,(float,np.floating)):return float(x) if np.isfinite(x) else None
    return x
def read_prior(script):
    spec=importlib.util.spec_from_file_location("hdblast_legacy_modes",script)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod
class NativeSmoothSpline:
    """Native B-spline fields; high-precision one-sided derivative limits."""
    def __init__(self,bs):
        self.bs=bs
        self.x=bs.t.copy()
        self.pp_roots=PPoly.from_spline(bs)
        self._mp=[]
        with mp.workdps(80):
            tt=[mp.mpf(float(z)) for z in bs.t]
            count=len(bs.t)-bs.k-1
            cc=[mp.mpf(float(z)) for z in bs.c[:count]]
            degree=bs.k
            for order in range(bs.k+1):
                self._mp.append((tt[:],cc[:],degree,np.array([float(z) for z in tt])))
                if degree:
                    dd=[]
                    for i in range(len(cc)-1):
                        denom=tt[i+degree+1]-tt[i+1]
                        if denom==0: raise ValueError("Unsupported repeated interior derivative knot")
                        dd.append(degree*(cc[i+1]-cc[i])/denom)
                    cc=dd;tt=tt[1:-1];degree-=1
    def __call__(self,t,nu=0):
        if nu>self.bs.k:
            return np.zeros_like(np.asarray(t,dtype=float))
        ans=self.bs(t,nu=nu)
        if not np.all(np.isfinite(ans)): raise ValueError("Nonfinite native B-spline evaluation")
        return ans
    def one_sided_mp(self,t,left,order):
        if order>self.bs.k:return mp.mpf("0")
        knots,coeffs,q,floats=self._mp[order]
        span=int(np.searchsorted(floats,t,side="left" if left else "right")-1)
        span=min(max(span,q),len(coeffs)-1)
        with mp.workdps(80):
            xx=mp.mpf(float(t))
            dd=[coeffs[span-q+j] for j in range(q+1)]
            for rr in range(1,q+1):
                for j in range(q,rr-1,-1):
                    idx=span-q+j
                    denom=knots[idx+q-rr+1]-knots[idx]
                    if denom==0:raise ValueError("Zero deBoor span")
                    alpha=(xx-knots[idx])/denom
                    dd[j]=(1-alpha)*dd[j-1]+alpha*dd[j]
            return +dd[q]
    def solve(self,target,extrapolate=False):
        candidates=self.pp_roots.solve(target,extrapolate=False)
        lo=float(self.bs.t[self.bs.k]);hi=float(self.bs.t[-self.bs.k-1])
        roots=[]
        for z in candidates:
            if not np.isfinite(z) or z<lo or z>hi:continue
            val=float(newton(lambda t:float(self.bs(t)-target),float(z),
                fprime=lambda t:float(self.bs(t,nu=1)),tol=5e-13,maxiter=30))
            if lo<=val<=hi and abs(float(self.bs(val)-target))<=5e-11:
                roots.append(val)
            else:raise ValueError("Native root refinement failed")
        return np.unique(np.array(roots,dtype=float))
def pp_side(pp,t,left,order):
    if isinstance(pp,NativeSmoothSpline):return float(pp.one_sided_mp(t,left,order))
    j=np.searchsorted(pp.x,t,side="left" if left else "right")-1
    j=int(np.clip(j,0,len(pp.x)-2))
    return float(np.polyval(np.polyder(pp.c[:,j],order),t-pp.x[j]))
def build(prior,data,kind):
    if kind in ("hermite","cubic"):
        b=prior.Geometry(data,kind);return b.la,b.ph
    k=DEGREES[kind];tt=data["H0tau"];ll=data["ln_a"]-data["ln_a"][0]
    return NativeSmoothSpline(make_interp_spline(tt,ll,k=k)),NativeSmoothSpline(make_interp_spline(tt,data["phi_b"],k=k))
def fields(la,ph,t,G,star):
    l=la(t);a=np.exp(l);h,hd,hdd,hddd=[la(t,nu=i) for i in range(1,5)]
    phi=ph(t);vel=ph(t,nu=1);x=G*G*(phi-star)**2
    U=a*a*(x-hd-2*h*h);r=6*(hd+2*h*h);rd=6*(hdd+4*h*hd)
    rdd=6*(hddd+4*hd*hd+4*h*hdd)
    er2=6*r*hd-12*h*rd;pr2=2*r*hd+4*rdd+8*h*rd
    ans=np.column_stack([t,l,phi,a,h,vel,hd,hdd,hddd,U,r,rd,rdd,er2,pr2])
    if not np.isfinite(ans).all():raise ValueError('Nonfinite field or high derivative')
    return ans
def side_geometry(la,ph,t,left,G,star):
    L=np.array([pp_side(la,t,left,i) for i in range(8)])
    F=np.array([pp_side(ph,t,left,i) for i in range(3)])
    a=np.exp(L[0]);h,hd,hdd=L[1:4];x=G*G*(F[0]-star)**2
    xd=2*G*G*(F[0]-star)*F[1]
    ans=dict(L=L,F=F,U=a*a*(x-hd-2*h*h),
       Up=a**3*(xd+2*h*x-4*h**3-6*h*hd-hdd),
       R=6*(hd+2*h*h),Rd=6*(hdd+4*h*hd),a=a)
    if not (np.isfinite(L).all() and np.isfinite(F).all() and np.isfinite([ans[k] for k in ("a","U","Up","R","Rd")]).all()):
        raise ValueError("Nonfinite one-sided geometry")
    return ans
def audit_one(prior,data,kind,out,tag):
    la,ph=build(prior,data,kind);tt=data["H0tau"];saved=tt[(tt>=0)&(tt<=END)]
    grid=np.unique(np.concatenate([GRID,saved]))
    vals=fields(la,ph,grid,prior.G,prior.PHISTAR)
    uniform=fields(la,ph,GRID,prior.G,prior.PHISTAR)
    if not np.isfinite(vals).all():raise ValueError("Nonfinite reconstructed field")
    np.savetxt(out/(tag+"_"+kind+"_geometry.csv"),vals,delimiter=",",header=",".join(COLS),comments="")
    knots=np.unique(np.concatenate([la.x,ph.x]));knots=knots[(knots>0)&(knots<END)]
    bounds=np.concatenate([[0.],knots,[END]]);eta=np.zeros(len(bounds))
    for i,(l,r) in enumerate(zip(bounds[:-1],bounds[1:])):
        eta[i+1]=eta[i]+quad(lambda t:float(np.exp(-la(t))),l,r,epsabs=1e-11,epsrel=1e-11,limit=100)[0]
    if not (np.isfinite(eta).all() and np.all(np.diff(eta)>0)):
        raise ValueError("Invalid conformal-time quadrature")
    af=float(np.exp(la(END)));jumps=[];max_cont=0.
    derivative_scale=[max(float(np.max(np.abs(la(grid,nu=i)))),1.) for i in range(8)]
    fscale=[max(float(np.max(np.abs(ph(grid,nu=i)))),1.) for i in range(3)]
    for i,t in enumerate(knots):
        le=side_geometry(la,ph,t,True,prior.G,prior.PHISTAR)
        ri=side_geometry(la,ph,t,False,prior.G,prior.PHISTAR)
        dl=ri["L"]-le["L"];df=ri["F"]-le["F"]
        if isinstance(la,NativeSmoothSpline):
            with mp.workdps(80):
                for poly,scales,maxorder in [(la,derivative_scale,CONTINUITY[kind]),(ph,fscale,min(CONTINUITY[kind],2))]:
                    for q in range(maxorder+1):
                        residual=abs(poly.one_sided_mp(t,False,q)-poly.one_sided_mp(t,True,q))
                        max_cont=max(max_cont,float(residual/mp.mpf(scales[q])))
        for q in range(CONTINUITY[kind]+1):max_cont=max(max_cont,abs(dl[q])/derivative_scale[q])
        for q in range(min(CONTINUITY[kind],2)+1):max_cont=max(max_cont,abs(df[q])/fscale[q])
        if kind=="hermite":n=0;A=-float(np.exp(2*la(t)))*dl[2]
        elif kind=="cubic":n=1;A=-float(np.exp(3*la(t)))*dl[3]
        else:n=DEGREES[kind]-2;A=-float(np.exp(DEGREES[kind]*la(t)))*dl[DEGREES[kind]]
        jumps.append(dict(t=float(t),eta=float(eta[i+1]),delta_L=dl,delta_phi=df,
            delta_U=ri["U"]-le["U"],delta_U_eta_prime=ri["Up"]-le["Up"],
            delta_R=ri["R"]-le["R"],delta_Rdot=ri["Rd"]-le["Rd"],
            first_potential_derivative_order=n,generic_first_jump=float(A)))
    A=np.array([z["generic_first_jump"] for z in jumps]);n=(0 if kind=="hermite" else 1 if kind=="cubic" else DEGREES[kind]-2)
    sumsq=float(np.dot(A,A))
    if not (np.isfinite(A).all() and np.isfinite(sumsq) and np.isfinite(af) and af>0):
        raise ValueError("Nonfinite jump energy coefficient")
    if n==0:energy=dict(type="logarithmic divergent formal high-k excitation energy",diagonal_log_coefficient=sumsq/(32*np.pi**2*af**4))
    else:energy=dict(type="finite formal dephased high-k energy tail",tail_power=2*n,tail_coefficient=sumsq/(2**(2*n+5)*np.pi**2*af**4*(2*n)))
    energy["sum_squared_first_jumps"]=sumsq
    energy["interpretation"]="Interior-knot WKB scattering only; endpoints and particle-basis contributions excluded."
    h0=float(la(0.,nu=1));hd0=float(la(0.,nu=2));a0=float(np.exp(la(0.)));hf=float(la(END,nu=1))
    curv0=a0*a0*(hd0+2*h0*h0)
    state=dict(initial_conformal_curvature=curv0,
        initial_physical_WKB_log_coefficient=curv0**2/(32*np.pi**2*af**4),
        initial_physical_Hamiltonian_K2_coefficient=(a0*h0)**2/(16*np.pi**2*af**4),
        out_physical_Hamiltonian_K2_coefficient=(af*hf)**2/(16*np.pi**2*af**4),
        interpretation="Formal UV state/reference extension versus curvature-aware admissible data. The final physical-Hamiltonian particle basis is basis-dependent; its tail alone is not a state certificate.")
    crossings=ph.solve(prior.PHISTAR,extrapolate=False)
    crossings=crossings[np.isfinite(crossings)&(crossings>=0)&(crossings<=END)]
    hdis=float(np.max(np.abs(la(saved,nu=1)-np.interp(saved,tt,data["H_over_H0"]))))
    vdis=float(np.max(np.abs(ph(saved,nu=1)-np.interp(saved,tt,data["v_over_H0"]))))
    np.savez_compressed(out/(tag+"_"+kind+"_geometry.npz"),grid=grid,geometry=vals,uniform=uniform)
    (out/(tag+"_"+kind+"_jumps.json")).write_text(json.dumps(clean(jumps),indent=2,allow_nan=False)+"\n")
    report=dict(tag=tag,kind=kind,continuity=CONTINUITY[kind],degree=DEGREES[kind],n_saved_in_window=len(saved),
        n_knots=len(knots),n_union_grid=len(grid),one_sided_evaluator=("native BSpline, cached 80-digit derivative coefficients/deBoor" if isinstance(la,NativeSmoothSpline) else "original PPoly"),a_final=af,eta_final=float(eta[-1]),crossings=crossings,
        saved_derivative_disagreement=dict(H_absolute=hdis,phi_dot_absolute=vdis,
          H_scaled=hdis/max(float(np.max(np.abs(data["H_over_H0"]))),1e-30),
          phi_dot_scaled=vdis/max(float(np.max(np.abs(data["v_over_H0"]))),1e-30)),
        continuity_max_scaled=float(max_cont),continuity_pass=bool(max_cont<=LIMITS["continuity_scaled"]),
        maxima={name:float(np.max(np.abs(vals[:,i]))) for i,name in enumerate(COLS) if i>0},
        knot_jump_maxima={name:max([abs(z[name]) for z in jumps],default=0.) for name in
          ("delta_U","delta_U_eta_prime","delta_R","delta_Rdot","generic_first_jump")},
        first_potential_derivative_jump_order=n,beta_power=n+2,energy_tail=energy,state_boundaries=state,
        regularity_scope="Finite-order interpolation, not a physical solution or all-orders Hadamard completion",
        local_R2_scope="Unit-action functional-derivative regular pieces only; knot distributions omitted and diagnosed by curvature jumps.")
    if kind=="hermite":
        y=data["ln_a"]-data["ln_a"][0];d=data["H_over_H0"];dt=np.diff(tt)
        c3=(2*y[:-1]-2*y[1:]+dt*(d[:-1]+d[1:]))/dt**3
        c2=(-3*y[:-1]+3*y[1:]-dt*(2*d[:-1]+d[1:]))/dt**2
        direct=[]
        for t in knots:
            j=int(np.searchsorted(tt,t))
            direct.append(-np.exp(2*la(t))*(2*c2[j]-2*c2[j-1]-6*c3[j-1]*dt[j-1]))
        err=float(np.max(np.abs(A-direct)))
        scale=max(float(np.max(np.abs(A))),1.)
        report["independent_Hermite_formula_error_scaled"]=err/scale
        report["independent_Hermite_formula_pass"]=err/scale<=1e-9
    return report,uniform,jumps
def abrupt_control():
    kk=np.array([8.,16.,32.,64.,128.,256.]);wm=np.sqrt(kk*kk+1);wp=np.sqrt(kk*kk+4)
    n=9/(4*wm*wp*(wm+wp)**2);ratios=16*kk**4*n/9
    beta=(wp-wm)/(2*np.sqrt(wm*wp));absdiff=float(np.max(np.abs(n-beta*beta)))
    return dict(k=kk,n=n,leading_ratio=ratios,exact_formula_difference=absdiff,
        asymptotic_final_relative=abs(float(ratios[-1]-1)),
        pass_=bool(abs(ratios[-1]-1)<=LIMITS["step_asymptotic"] and absdiff<1e-14),
        scope="Exact single-step scattering; not an archived source")
def opposite_jump_control():
    kl,ku=1e3,1e6
    diagonal=2*np.log(ku/kl);interference=-2*(float(sici(2*ku)[1])-float(sici(2*kl)[1]))
    ratio=(diagonal+interference)/diagonal
    return dict(jumps=[1.,-1.],eta=[0.,1.],K=[kl,ku],diagonal_log_term=float(diagonal),
        interference_term=float(interference),coefficient_ratio=float(ratio),
        diagonal_energy_log_coefficient=2/(32*np.pi**2),
        pass_=bool(abs(ratio-1)<=LIMITS["interference_log_slope"]),
        scope="Exact integral of leading opposite-jump amplitudes; no archive-mode fit")
def differences(a,b):
    out={}
    for i,name in enumerate(COLS):
        if i==0:continue
        err=float(np.max(np.abs(a[:,i]-b[:,i])))
        scale=max(float(np.max(np.abs(b[:,i]))),1e-30)
        out[name]=dict(max_abs=err,scaled_by_primary_max=err/scale)
    return out
def main():
    p=argparse.ArgumentParser();base=Path(__file__).resolve().parent.parent
    p.add_argument("--inputs",type=Path,default=base/"data/A1")
    p.add_argument("--prior-script",type=Path,default=base/"code/legacy_quantum_source.py")
    p.add_argument("--out",type=Path,default=base/"outputs")
    p.add_argument("--registration",type=Path,required=True);p.add_argument("--registration-sha256",required=True)
    p.add_argument("--source-ref",default="unspecified");ar=p.parse_args()
    regsha=hashlib.sha256(ar.registration.read_bytes()).hexdigest()
    if regsha!=ar.registration_sha256:raise ValueError("Registration hash mismatch")
    ar.out.mkdir(parents=True,exist_ok=True);prior=read_prior(ar.prior_script)
    report=dict(source_ref=ar.source_ref,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      registration_sha256=regsha,limits=LIMITS,t_window=[0,END],uniform_grid_n=len(GRID),
      reconstructions=[],sources=[],errors={},toy_controls=dict(abrupt_step=abrupt_control(),opposite_jumps=opposite_jump_control()),
      assumption="Minimal scalar,G=100,phi_star=.5; no decay/backreaction. Geometry audit only.")
    curves={}
    for tag in (prior.FINE,prior.COARSE):
        try:data,meta=prior.load(ar.inputs,tag);report["sources"].append(meta)
        except Exception as e:report["errors"][tag]=type(e).__name__+": "+str(e);continue
        for kind in KINDS:
            key=tag+":"+kind
            try:
                rr,cv,jj=audit_one(prior,data,kind,ar.out,tag);report["reconstructions"].append(rr);curves[key]=cv
                print("HDBLAST_FRW_GEOMETRY_CASE="+json.dumps(clean(rr),separators=(",",":"),allow_nan=False),flush=True)
                if tag==prior.FINE and kind=="hermite":
                    print("HDBLAST_FRW_PRIMARY_JUMPS="+json.dumps(clean(jj),separators=(",",":"),allow_nan=False),flush=True)
            except Exception as e:report["errors"][key]=type(e).__name__+": "+str(e)
    ref=prior.FINE+":hermite"
    if ref in curves:
        report["sensitivity_to_fine_hermite"]={key:differences(cv,curves[ref]) for key,cv in curves.items() if key!=ref}
        keys=[ref]+[prior.FINE+":"+k for k in KINDS[1:] if prior.FINE+":"+k in curves]
        np.savez_compressed(ar.out/"plot_reconstructions.npz",columns=np.array(COLS),labels=np.array(keys),
                            curves=np.array([curves[key] for key in keys]))
        sel=[0,2,4,6,9,10,12,14]
        report["selected_curves"]=dict(columns=[COLS[j] for j in sel],samples={key:curves[key][::10][:,sel] for key in keys},
                                     scope="Every tenth uniform sample; full curves retained in CSV/NPZ.")
    report["controls_pass"]=(not report["errors"] and len(report["reconstructions"])==8
       and all(z["continuity_pass"] and z.get("independent_Hermite_formula_pass",True) for z in report["reconstructions"])
       and all(z["pass_"] for z in report["toy_controls"].values()))
    report["absolute_FRW_source_computed"]=False;report["Hadamard_state_certified"]=False
    report["scientific_scope"]="Knot/state UV diagnostics and high-derivative sensitivity. Smooth reconstruction,state,prehistory andcurved finite matching still required."
    (ar.out/"geometry_audit.json").write_text(json.dumps(clean(report),indent=2,allow_nan=False)+"\n")
    print("HDBLAST_FRW_GEOMETRY_JSON="+json.dumps(clean(report),separators=(",",":"),allow_nan=False),flush=True)
    return int(not report["controls_pass"])
if __name__=="__main__":raise SystemExit(main())

