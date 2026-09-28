"""Floating-point initial constraint data, not a modified evolution model.

Requires NumPy and SciPy. No files are read from the user's research folders.
The positive-bump construction satisfies first but generically not second time
compatibility at the brane; report this failure explicitly, never as evolution.
"""
from pathlib import Path
import argparse, json, math
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

IP=1.0357712571566784
C=2/IP-4/3

def potential(psi):
    p=psi-1; w=5/3-psi*psi+psi**3/3; w1=psi*(psi-2)
    return (.5*w1*w1-2*w*w/3, w1*(2*p-4*w/3),
            4*p*p+2*w1-(4/3)*(w1*w1+2*p*w))

def tension(psi,tdet=.001):
    p=psi-1;w=5/3-psi*psi+psi**3/3
    return 2*w+tdet*(1+C*p),2*psi*(psi-2)+tdet*C,4*p

def bump(psi):
    x=(psi-.5)/.3
    return np.exp(1-1/(1-x*x)) if abs(x)<1 else 0.

def integrate(x,eps,rtol=2e-11,dense=False,y0=1e-4,shape=bump):
    psi_h=10**x[0];yb=x[1];u,u1,u2=potential(psi_h)
    a=-u/36;b=u*u/4320-u1*u1/750;cc=u1/10;d=u1*(u2/280+u/630)
    init=[y0+a*y0**3+b*y0**5,1+3*a*y0*y0+5*b*y0**4,
          psi_h+cc*y0*y0+d*y0**4,2*cc*y0+4*d*y0**3,0.,0.]
    # Last two entries accumulate conformal time and integral for rho^2 D.
    def rhs(y,v):
        rho,ry,psi,sy,_,_=v;uu,up,_=potential(psi);bb=shape(psi)
        ryy=(1-ry*ry)/rho-rho*sy*sy/6-rho*uu/3
        syy=up-4*ry*sy/rho+eps*bb
        return [ry,ryy,sy,syy,1/rho,eps*rho**4*sy*bb/6]
    sol=solve_ivp(rhs,(y0,yb),init,method='DOP853',rtol=rtol,
                  atol=[1e-13,1e-13,1e-22,1e-22,1e-13,1e-20],
                  max_step=.08,dense_output=dense)
    if not sol.success:raise RuntimeError(sol.message)
    return sol

def solve_seed(eps,guess=None,rtol=2e-11,y0=1e-4):
    guess=guess if guess is not None else [math.log10(8.7855146e-7),8.22820155]
    def residual(x):
        v=integrate(x,eps,rtol,y0=y0).y[:,-1];rho,ry,psi,sy=v[:4];sig,sig1,_=tension(psi)
        return np.array([ry/rho-sig/6,sy+sig1/2])
    found=root(residual,guess,tol=1e-9,options={'eps':1e-8})
    residuals=residual(found.x)
    if max(abs(residuals))>2e-9:raise RuntimeError('Junction residual too large: '+str(residuals))
    sol=integrate(found.x,eps,rtol,True,y0)
    ys=np.linspace(y0,found.x[1],3001);v=sol.sol(ys);rho,ry,psi,sy=v[:4]
    uu,up,_=potential(psi);ryy=(1-ry*ry)/rho-rho*sy*sy/6-rho*uu/3
    az=ry;bz=ry;azz=rho*ryy;phiz=rho*sy
    ham=-2*uu*rho*rho+6-12*az*az+6*az*bz-6*azz-phiz*phiz
    momentum=np.zeros_like(ham)
    defect=1-ry*ry+rho*rho*(sy*sy/12-uu/6)
    # Check the separate integral identity as a nontrivial evolution-ODE crosscheck.
    ledger=rho*rho*defect-v[5]
    rb,vb,psib,sb=sol.y[:4,-1];ub,upb,_=potential(psib);sig,sig1,sig2=tension(psib)
    db=1-vb*vb+rb*rb*(sb*sb/12-ub/6)
    # b(phi) vanishes on an open neighborhood of the brane for these roots.
    if bump(psib)!=0.:raise RuntimeError('Bump reaches brane; corner formulas do not apply')
    scale=1+6*ry*ry+phiz*phiz
    return sol,found.x,dict(epsilon=eps,rtol=rtol,y0=y0,root_success=bool(found.success),
      root_message=str(found.message),psi_h=float(10**found.x[0]),yb=float(found.x[1]),
      rho_b=float(rb),phi_b=float(psib-1),phi_y_b=float(sb),rho_y_b=float(vb),
      junction_residual=residuals.tolist(),hamiltonian_algebraic_max=float(max(abs(ham))),
      hamiltonian_algebraic_normalized_max=float(max(abs(ham)/scale)),
      momentum_algebraic_max=0.,first_time_corner_residual=[0.,0.,0.],
      mass_defect_b=float(db),mass_defect_integral=float(sol.y[5,-1]/rb**2),
      mass_defect_ledger_absolute_max=float(max(abs(ledger))),
      mass_defect_ledger_scaled_max=float(max(abs(ledger)/(1+rho**2+abs(v[5])))),
      second_time_corner_A=float(4*vb*db-(rb*sig*4*db/6)),
      second_time_corner_B=float(-8*vb*db-(rb*sig*4*db/6)),
      second_time_corner_phi=float(2*rb*sig1*db),
      second_time_corner_status='NOT_COMPATIBLE_GENERALLY' if eps else 'STATIC_NUMERICAL_CONTROL',
      not_claimed=['interval existence certificate','complete corner compatibility','nonlinear evolution','hot-universe mechanism'])

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default=str(Path(__file__).resolve().parent))
    args=ap.parse_args();out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    rows=[];allprofiles={}
    for eps in [0.,1e-8,-1e-8,1e-6,-1e-6,1e-4,-1e-4]:
        sol,x,row=solve_seed(eps)
        fine,xf,rowf=solve_seed(eps,x,rtol=3e-13,y0=5e-5)
        row['refined']=rowf
        row['refinement_delta']={k:row[k]-rowf[k] for k in ['psi_h','yb','rho_b','phi_b','mass_defect_b']}
        # Monotone conformal coordinate sampled on uniform proper distance.
        ys=np.linspace(5e-5,xf[1],4097);values=fine.sol(ys);values[4]-=values[4,-1]
        allprofiles['eps_'+str(eps)]=np.vstack([ys,values])
        rows.append(row);print(json.dumps({k:row[k] for k in ['epsilon','phi_b','junction_residual','mass_defect_b','second_time_corner_B','second_time_corner_phi']}),flush=True)
    np.savez_compressed(out/'CONSTRAINT_SEED_PROFILES.npz',**allprofiles)
    result={'status':'FIRST_ORDER_COMPATIBLE_INITIAL_DATA_ONLY','tdet':.001,'d':0.,'c':C,
            'scalar_bump':'exp(1-1/(1-x^2)) on |x|<1, x=(phi+0.5)/0.3; zero elsewhere',
            'equations_changed_in_evolution':False,'rows':rows}
    (out/'CONSTRAINT_SEED_RESULTS.json').write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
