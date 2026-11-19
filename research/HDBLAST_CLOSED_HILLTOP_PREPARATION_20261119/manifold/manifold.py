"""Nonzero-velocity homogeneous closed-bounce separatrices and re-centered launch thresholds."""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from numpy.polynomial import Polynomial
from scipy.integrate import solve_bvp, solve_ivp
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent


class Field:
    def __init__(self, row):
        self.row = row
        self.knots = np.array(row["spline_knots"])
        self.coef = np.array(row["phi_f_coefficients"])
        self.pm = row["minimum"]
        self.b = row["b"]
        self.v0 = row["v0"]
        self.top = row["top"]
        self.hstar = math.sqrt(self.values(0)[0])
        self.mu2 = -3*self.values(0)[2]
        self.decay = (-3*self.hstar - math.sqrt(9*self.hstar**2 + 4*self.mu2))/2
        self.central = self.index(0)
        # Re-expand the rational force around the stationary point in the two
        # adjacent cells. This evaluates the small force without cancellation.
        self.force = {}
        for k in [self.central-1, self.central]:
            p, f, num = self.polynomials(k)
            dernum = num.deriv()*f - 2*num*f.deriv()
            self.force[k] = (dernum, f)

    def index(self, q):
        k = np.searchsorted(self.knots, self.top+q, side="right")-1
        if np.any((k<0) | (k>=len(self.knots)-1)):
            raise ValueError("Outside extracted spline neighborhood")
        return k

    def polynomials(self, k):
        center = Polynomial([self.top-self.knots[k], 1.])
        p, f = [Polynomial(c[::-1])(center) for c in self.coef[:, k].T]
        num = (1-p/self.pm)**2*(1+self.b*p)/self.v0
        return p, f, num

    def third(self, k):
        _, f, num = self.polynomials(k)
        den = (f*f).coef
        series = []
        for j in range(4):
            series.append((num.coef[j] - sum(den[i]*series[j-i] for i in range(1,j+1)))/den[0])
        return 6*series[3]

    def values(self, q):
        q = np.asarray(q)
        k = self.index(q)
        x = self.top+q-self.knots[k]
        c3,c2,c1,c0 = self.coef[:, k] if q.ndim==0 else self.coef[:, k].transpose(0,2,1)
        p,f = ((c3*x+c2)*x+c1)*x+c0
        dp,df = (3*c3*x+2*c2)*x+c1
        ddp,ddf = 6*c3*x+2*c2
        z=1-p/self.pm
        j=1+self.b*p
        u=z*z*j/(f*f*self.v0)
        wp=-2/(self.pm*z)+self.b/j
        wpp=-2/(self.pm*z)**2-(self.b/j)**2
        l1=wp*dp-2*df/f
        l2=wpp*dp*dp+wp*ddp-2*(ddf/f-(df/f)**2)
        du=u*l1
        for i,(numerator,denominator) in getattr(self,"force",{}).items():
            # q is relative to the declared stationary point. Set the constant
            # on each adjoining polynomial to zero (both at roundoff level).
            du=np.where(k==i, q*Polynomial(numerator.coef[1:])(q)/denominator(q)**3, du)
        return u,du,u*(l1*l1+l2)

    def increments(self, q, dq):
        # Two-point Gauss-Legendre integrals of U' and U'' remain accurate
        # even when q+dq rounds to q.
        half=dq/2
        middle=q+half
        step=half/math.sqrt(3)
        _, a, b=self.values(middle-step)
        _, c, d=self.values(middle+step)
        return half*(a+c),half*(b+d)


def stable(field, velocity, end=5., tol=5e-10):
    hstar=field.hstar
    grid=np.linspace(0,end,201)
    amp=-velocity/math.sqrt(field.mu2)
    shape=amp*np.exp(-math.sqrt(field.mu2)*grid)
    guess=np.vstack([shape,-math.sqrt(field.mu2)*shape,hstar*np.tanh(hstar*grid),
                     np.log(np.cosh(hstar*grid)/hstar)])
    def ode(t,Y):
        q,p,h,l=Y
        u,du,_=field.values(q)
        return np.vstack([p,-3*(h*p+du),u-p*p/3-h*h,h])
    def bc(start,stop):
        u=field.values(start[0])[0]
        return np.array([start[1]-velocity,start[2],start[3]+.5*np.log(u+velocity**2/6),
                         stop[1]-field.decay*stop[0]])
    sol=solve_bvp(ode,bc,grid,guess,tol=tol,max_nodes=10000)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol


def taylor(field):
    end=4.
    x=np.linspace(0,end,161)
    mu2=field.mu2
    third=field.third(field.central)
    def ode(t,Y):
        s,sp,k,kp=Y
        h=field.hstar*np.tanh(field.hstar*t)
        return np.vstack([sp,mu2*s-3*h*sp,kp,mu2*k-3*h*kp-1.5*third*s*s])
    def bc(a,b):
        return [a[1]-1,a[3],b[1]-field.decay*b[0],b[3]-field.decay*b[2]]
    guess=np.zeros((4,len(x)))
    sol=solve_bvp(ode,bc,x,guess,tol=2e-11,max_nodes=10000)
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[0,0],sol.y[2,0]


def perturbed(field, base, seed, end, target_l=None):
    # Normalized differences carry the exact requested initial field separation,
    # even when that separation is much smaller than a binary64 ulp at the center.
    q0,p0,h0,l0=base.sol(0)
    du0,_=field.increments(q0,seed)
    init=[1.,0.,0.,-.5*np.log1p(du0/(field.values(q0)[0]+p0*p0/6))/seed]
    def ode(t,z):
        q,p,h,l=base.sol(t)
        dq,dp,dh,dl=z
        du,ddu=field.increments(q,seed*dq)
        return [dp,-3*(h*dp+dh*(p+seed*dp)+ddu/seed),
                du/seed-(2*p+seed*dp)*dp/3-(2*h+seed*dh)*dh,dh]
    def at_l(t,z):
        return base.sol(t)[3]+seed*z[3]-target_l
    at_l.terminal=True
    at_l.direction=1
    sol=solve_ivp(ode,(0,end),init,method="DOP853",rtol=1e-10,atol=1e-12,
                  events=at_l if target_l is not None else None)
    if not sol.success or (target_l is not None and not len(sol.t_events[0])):
        raise RuntimeError("Deviation integration or radius event failed")
    x=sol.t[-1]
    return x, base.sol(x)+seed*sol.y[:,-1]


def direct_overlap(field,base,seed):
    # At a resolvable finite offset an ordinary IVP also tests the re-centering.
    q,p,h,l=base.sol(0)
    q+=seed
    initial=[q,p,0.,-.5*math.log(field.values(q)[0]+p*p/6)]
    def ode(t,y):
        q,p,h,l=y
        u,du,_=field.values(q)
        return [p,-3*(h*p+du),u-p*p/3-h*h,h]
    end=.55
    sol=solve_ivp(ode,(0,end),initial,method="DOP853",rtol=2e-12,atol=2e-15)
    if not sol.success:raise RuntimeError(sol.message)
    _,yp=perturbed(field,base,seed,end)
    return max(abs((sol.y[:2,-1]-yp[:2])/(yp[:2]-base.sol(end)[:2])))


def evaluate(row):
    field=Field(row)
    if abs(field.mu2-row["source_mu_squared"])>1e-8:
        raise RuntimeError("Spline curvature reconstruction mismatch")
    b1,b2=taylor(field)
    sign=-np.sign(field.pm)
    log_d=row["critical_log_d"]
    seed=sign*math.exp(log_d)
    x0=row["critical_x"]
    l_target=math.log(math.cosh(field.hstar*x0)/field.hstar)
    base0=stable(field,0.)
    _,control=perturbed(field,base0,seed,x0+.1,l_target)
    samples=[]
    for v in [-math.sqrt(.06),-.10,-.01,-.001,0.,.001,.01,.10,math.sqrt(.06)]:
        sol=stable(field,v)
        q0,p0,h0,l0=sol.sol(0.)
        q,p,h,l=sol.sol(np.linspace(0.,5.,1001))
        max_constraint=np.max(abs(h*h-field.values(q)[0]-p*p/6+np.exp(-2*l)))
        def match(log_ratio):
            _,point=perturbed(field,sol,seed*math.exp(log_ratio),x0+.25,l_target)
            return point[0]/control[0]-1.
        if v==0:log_ratio=0.
        else:log_ratio=brentq(match,-.65,.65,xtol=2e-10)
        ratio=math.exp(log_ratio)
        xm,point=perturbed(field,sol,seed*ratio,x0+.25,l_target)
        rec=dict(v=v,center=q0,linear_bend=q0-b1*v,quadratic_residual=q0-b1*v-b2*v*v,
                 a_b=math.exp(l0),bounce_hprime=field.values(q0)[0]-v*v/3,
                 max_constraint=max_constraint,conditional_width_ratio=ratio,
                 conditional_width=abs(seed)*ratio,match_time=xm,
                 relative_match_p_error=point[1]/control[1]-1,match_h_error=point[2]-control[2])
        if abs(v)>=.2 or abs(v)==.001:
            alt=stable(field,v,end=6.,tol=7e-11)
            rec["center_refinement_difference"]=float(alt.sol(0)[0]-q0)
            rec["direct_ivp_relative_difference"]=direct_overlap(field,alt,sign*2e-7)
        samples.append(rec)
    # An independent origin check at resolvable offsets compares odd/even parts.
    neg,pos=samples[3],samples[5]
    odd=(pos["center"]-neg["center"])/.002
    even=(pos["center"]+neg["center"])/(2e-6)
    if abs(odd/b1-1)>2e-5 or abs(even/b2-1)>2e-3:
        raise ArithmeticError("Origin Taylor check failed")
    for r in samples:
        if (r["max_constraint"]>8e-9 or abs(r["relative_match_p_error"])>1e-5 or
            abs(r.get("center_refinement_difference",0.))>2e-11 or
            r.get("direct_ivp_relative_difference",0.)>2e-5):
            raise ArithmeticError(r)
    return dict(critical_zero_v_distance=abs(seed),match_loga=l_target,
                actual_zero_v_match_q=control[0],linear_slope=b1,quadratic_coefficient=b2,
                Uthird_right=field.third(field.central),Uthird_left=field.third(field.central-1),
                velocity_at_quadratic_bend_equal_to_zero_v_distance=math.sqrt(abs(seed/b2)),
                finite_difference_linear=odd,finite_difference_quadratic=even,samples=samples)


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--prior",type=Path,default=HERE.parent/"prior")
    cli.add_argument("--output",type=Path,required=True)
    args=cli.parse_args()
    spec=importlib.util.spec_from_file_location("extract_local_spline",HERE/"extract.py")
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source=module.extract(args.prior)
    result=dict(prior_sha256=source["prior_sha256"],branches={})
    for k,row in source["branches"].items():
        result["branches"][k]=evaluate(row)
        print("Completed "+k,flush=True)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    main()
