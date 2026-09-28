"""Exploratory compensated compact initial seed; no evolution is performed."""
from pathlib import Path
import json,math
import numpy as np
from scipy.optimize import root
from scipy.integrate import solve_ivp
import constraint_seed as C

def local_bump(psi,center,width):
    x=(psi-center)/width
    return np.exp(1-1/(1-x*x)) if abs(x)<1 else 0.

def run(eps,rtol=3e-13,transport=False):
    def shape(psi,k):return local_bump(psi,.35,.15)-k*local_bump(psi,.70,.15)
    def compute(x,dense=False):
        if transport:
            y0=5e-5;ph=10**x[0];u,u1,u2=C.potential(ph)
            a=-u/36;b=u*u/4320-u1*u1/750;cc=u1/10;dd=u1*(u2/280+u/630)
            ini=[y0+a*y0**3+b*y0**5,ph+cc*y0*y0+dd*y0**4,2*cc*y0+4*dd*y0**3,0.,0.]
            def rhs(y,V):
                r,p,s,z,I=V;uu,up,_=C.potential(p);v=np.sqrt(1+r*r*(s*s/12-uu/6)-I/r**2)
                return [v,s,up-4*v*s/r+eps*shape(p,x[2]),1/r,eps*r**4*s*shape(p,x[2])/6]
            ss=solve_ivp(rhs,(y0,x[1]),ini,method='DOP853',rtol=rtol,
                         atol=[1e-13,1e-22,1e-22,1e-13,1e-20],max_step=.08,dense_output=dense)
            if not ss.success:raise RuntimeError(ss.message)
            def adapt(V):
                r,p,s,z,I=V;uu,_,_=C.potential(p);v=np.sqrt(1+r*r*(s*s/12-uu/6)-I/r**2)
                return np.array([r,v,p,s,z,I])
            ss.y=adapt(ss.y)
            if dense:
                old=ss.sol;ss.sol=lambda y:adapt(old(y))
            return ss
        return C.integrate(x[:2],eps,rtol,dense,y0=5e-5,shape=lambda psi:shape(psi,x[2]))
    def residual(x):
        sol=compute(x);rho,ry,psi,sy,_,I=sol.y[:,-1];sig,sig1,_=C.tension(psi)
        return [(ry/rho-sig/6)/.001,(sy+sig1/2)/.001,I/(eps*rho**4) if eps else x[2]-.40932878]
    found=root(residual,[math.log10(8.7855146e-7),8.22820155,.5],tol=1e-9,options={'eps':1e-8})
    residuals=np.array(residual(found.x));sol=compute(found.x,True)
    if max(abs(residuals))>2e-7:raise RuntimeError('Compensated solve failed: '+str(residuals))
    rho,ry,psi,sy,_,I=sol.y[:,-1];u,up,_=C.potential(psi);D=1-ry*ry+rho*rho*(sy*sy/12-u/6)
    result={'epsilon':eps,'rtol':rtol,'transport_formulation':transport,'root_success':bool(found.success),'root_message':str(found.message),
      'psi_h':float(10**found.x[0]),'yb':float(found.x[1]),'compensation':float(found.x[2]),
      'rho_b':float(rho),'phi_b':float(psi-1),'junction_residual':(residuals[:2]*.001).tolist(),
      'scaled_compensation_residual':float(residuals[2]),'mass_defect_b':float(D),
      'mass_defect_integral':float(I/rho**2),'second_corner_B':float(-12*ry*D),
      'second_corner_phi':float(-4*rho*sy*D),'bump_at_brane':float(shape(psi,found.x[2])),
      'scope':'FLOATING_POINT_COMPENSATED_INITIAL_DATA_ONLY; no evolution; no interval proof'}
    return sol,result

if __name__=='__main__':
    rows=[];profiles={}
    for eps in [0.,1e-6,-1e-6,1e-4,-1e-4]:
        sol,row=run(eps);sol2,row2=run(eps,rtol=8e-14);row['refined']=row2
        sol3,row3=run(eps,rtol=8e-14,transport=True);row['independent_mass_transport']=row3
        rows.append(row);print(json.dumps(row),flush=True)
        ys=np.linspace(5e-5,row3['yb'],4097);v=sol3.sol(ys);v[4]-=v[4,-1]
        profiles[str(eps)]=np.vstack([ys,v])
    out=Path(__file__).resolve().parent
    (out/'BALANCED_CONSTRAINT_SEED_RESULTS.json').write_text(json.dumps({'status':'EXPLORATORY_COMPENSATED_DATA','rows':rows},indent=2)+'\n')
    np.savez_compressed(out/'BALANCED_CONSTRAINT_SEED_PROFILES.npz',**profiles)
