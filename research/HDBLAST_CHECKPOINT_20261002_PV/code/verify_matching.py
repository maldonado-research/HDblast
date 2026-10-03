"""Independent 70-digit analytic checks of subtraction functions and Jost data."""
import json
import numpy as np
import mpmath as mp
import pv_crossing as p
mp.mp.dps=70
c=[1,-3,3,-1]
pi=mp.pi
def ylog(y):
    return mp.mpf("0") if y==0 else y*mp.log(y)
def potential(x,L):
    return sum(c[j]*(x+j*L**2)*ylog(x+j*L**2) for j in range(4))/(64*pi**2)
def dpotential(x,L):
    return sum(c[j]*ylog(x+j*L**2) for j in range(4))/(32*pi**2)
def ddpotential(x,L):
    return sum(c[j]*mp.log(x+j*L**2) for j in range(4))/(32*pi**2)
def qtail(x,L,K):
    return sum(c[j]*(-K*mp.sqrt(K*K+x+j*L*L)+(x+j*L*L)*mp.log(K+mp.sqrt(K*K+x+j*L*L))) for j in range(4))/(8*pi**2)
def etail(x,L,K):
    return sum(c[j]*((x+j*L*L)**2*mp.log(K+mp.sqrt(K*K+x+j*L*L))-K*mp.sqrt(K*K+x+j*L*L)*(2*K*K+x+j*L*L)) for j in range(4))/(32*pi**2)
errors=dict(potential=0.,potential_derivative=0.,static_energy=0.,static_Q=0.,static_boundary=0.,trace_tail=0.,jost=0.,jost_derivative=0.,phase=0.)
count=0
for L0 in (4.,8.,16.):
    L=mp.mpf(L0);r=mp.mpf(4)
    for x0 in (0.,.01,.5,2.,3.99,4.):
        x=mp.mpf(x0);d=x-r
        refv=potential(x,L)-potential(r,L)-d*dpotential(r,L)-d*d*ddpotential(r,L)/2
        refvx=dpotential(x,L)-dpotential(r,L)-d*ddpotential(r,L)
        a,b=p.matched_potential(x0,L0)
        errors["potential"]=max(errors["potential"],abs(a-float(refv)))
        errors["potential_derivative"]=max(errors["potential_derivative"],abs(b-float(refvx)))
        for K0 in (4*L0,8*L0,12*L0):
            K=mp.mpf(K0)
            ee,qq,bb=p.static_tails(x0,L0,K0)
            bref=K**3*sum(c[j]*mp.sqrt(K*K+x+j*L*L) for j in range(4))/(12*pi**2)
            errors["static_energy"]=max(errors["static_energy"],abs(ee-float(etail(x,L,K))))
            errors["static_Q"]=max(errors["static_Q"],abs(qq-float(qtail(x,L,K))))
            errors["static_boundary"]=max(errors["static_boundary"],abs(bb-float(bref)))
            xd=mp.mpf("1.3");xdd=mp.mpf("-2.1")
            ref=-(mp.diff(lambda xx:qtail(xx,L,K),x,2)*xd**2+mp.diff(lambda xx:qtail(xx,L,K),x)*xdd)/2
            errors["trace_tail"]=max(errors["trace_tail"],abs(p.trace_tail(x0,float(xd),float(xdd),L0,K0)-float(ref)))
            count+=1
for t0 in (-10.,-8.):
    for L0 in (4.,16.):
        ks=np.array([.0001,.25,2.,8.,64.,128.])
        ff,fd,*_=p.jost(t0,ks,L0)
        for j in range(4):
            for i,k0 in enumerate(ks):
                t=mp.mpf(t0);k=mp.mpf(float(k0));L=mp.mpf(L0)
                om=mp.sqrt(k*k+4+j*L*L);z=(1+mp.tanh(t))/2
                a=(1+mp.sqrt(mp.mpc(-15)))/2;b=1-a;cc=1-1j*om
                F=mp.hyp2f1(a,b,cc,z);Fz=4/cc*mp.hyp2f1(a+1,b+1,cc+1,z)
                ref=F/mp.sqrt(2*om);refd=(-1j*om*F+2*z*(1-z)*Fz)/mp.sqrt(2*om)
                errors["jost"]=max(errors["jost"],abs(ff[j,i]-complex(ref)))
                errors["jost_derivative"]=max(errors["jost_derivative"],abs(fd[j,i]-complex(refd)))
for t0 in (-6.,-.1,0.,.1,6.):
    for L0 in (4.,16.):
        ks=np.array([.0001,.25,2.,8.,128.]);actual=p.phase(t0,ks,L0)
        for j in range(4):
            for i,k0 in enumerate(ks):
                b=mp.mpf(float(k0))**2+j*mp.mpf(L0)**2
                # Independent quadrature of the physical instantaneous frequency.
                ref=mp.quad(lambda tt:mp.sqrt(b+4*mp.tanh(tt)**2),[0,mp.mpf(t0)])
                errors["phase"]=max(errors["phase"],abs(actual[j,i]-float(ref)))
assert max(errors.values())<1e-8,errors
out=dict(status="PASS",precision_decimal_digits=70,static_cases=count,
         max_absolute_errors=errors,
         scope="Independent high-precision static integrals, differentiated trace tail, hypergeometric Jost state and phase quadrature")
print("HDBLAST_MATCHING_CHECK="+json.dumps(out))
