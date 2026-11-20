"""Displaced-vacuum response of the conformal zero mode on the #189 reference history."""
import argparse
import json
import math
from pathlib import Path
import runpy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.special import expit

HERE = Path(__file__).resolve().parent
PRIOR = runpy.run_path(str(HERE.parent/"prior/parent_modes.py"))
TARGET = 2**3*.0001/1.9

def step(z):
    z=np.asarray(z)
    active=(z>0)&(z<1)
    u=np.where(active,z,.5)
    f=expit(1/(1-u)-1/u)
    df=expit(-np.abs(1/(1-u)-1/u))*expit(np.abs(1/(1-u)-1/u))*(u**-2+(1-u)**-2)
    return np.where(active,f,(z>=1).astype(float)),np.where(active,df,0.)

def source(psi,start,width):
    x=(psi-start)/width
    s,ds=step(x)
    t,dt=step(x-1)
    return s-t,(ds-dt)/width

def mass(psi):
    f,df=step((psi-.77)/.03+.5)
    return 1.9+.2*f,.2*df/.03

def exact_cells(a,psi,end,start,width,g,spacing):
    n=math.ceil(end/spacing)
    dt=end/n
    mid=(np.arange(n)+.5)*dt
    radii=a(mid); fields=psi(mid)
    m=mass(fields)[0]; v=g*source(fields,start,width)[0]
    w=np.hypot(1.,radii*m)
    angle=w*dt/radii
    center=radii**3*m**2*v/w**2
    co,si=np.cos(angle),np.sin(angle)
    y=p=0.
    for c,s,w_i,q in zip(co,si,w,center):
        y,p=q+c*(y-q)+s*p/w_i,c*p-w_i*s*(y-q)
    return y,p

def evaluate(a,psi,end,start,width,g):
    da=a.derivative(); dpsi=psi.derivative()
    readouts=[brentq(lambda t:float(psi(t))-p,0.,end) for p in (start,start+width,start+2*width)]+[end]
    # Outside the drive and the mass ramp all work due to their field dependence vanishes.
    def state(t,vec):
        radius=float(a(t)); rdot=float(da(t)); field=float(psi(t)); field_dot=float(dpsi(t))
        m,dm=mass(field); shape,dshape=source(field,start,width)
        v=g*shape; vdot=g*dshape*field_dot
        y,p=vec[:2]
        delta=y-radius*v
        return (p/radius, (radius**3*m*m*v-(1+(m*radius)**2)*y)/radius,
                -radius**3*m*m*delta*vdot,
                radius**2*m*dm*field_dot*delta**2,
                radius*rdot*m*m*delta*(y-2*radius*v))
    sol=solve_ivp(state,(0.,end),np.zeros(5),method="DOP853",atol=1e-13,rtol=2e-10,
                  max_step=.035,t_eval=readouts)
    if not sol.success:raise RuntimeError(sol.message)
    y,p,wv,wm,wa=sol.y[:,-1]
    radius=float(a(end)); m=float(mass(float(psi(end)))[0]); omega=math.hypot(1,radius*m)
    energy=.5*(p*p+omega**2*y*y)
    work_error=abs(energy-wv-wm-wa)
    def action(values):
        y,p=values
        return .5*(p*p+omega*omega*y*y)/omega
    coarse=action(exact_cells(a,psi,end,start,width,g,.00024))
    fine=action(exact_cells(a,psi,end,start,width,g,.00012))
    n=energy/omega
    if work_error>1e-10 or abs(n-fine)>2e-9 or abs(fine-coarse)>1e-8:
        raise ArithmeticError((work_error,n,fine,coarse))
    return dict(g=g,drive_start=start,half_width=width,proper_times=readouts,
                mean_displacement_n_over_S0_volume=n,ratio_to_fluid_target=n/TARGET,
                source_amplitude_for_target=g*math.sqrt(TARGET/n) if n>0 else None,
                conformal_energy=energy,work_from_minimum=wv,work_from_mass=wm,
                work_from_scale_factor=wa,work_balance_error=work_error,
                independent_midpoint_fine_action=fine,midpoint_refinement=abs(fine-coarse))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source",type=Path,default=HERE.parent/"prior/CHARGED.json")
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args(); a,psi,times=PRIOR["history"](args.source)
    cases={"nominal":(.03,.075,.01),"half_amplitude":(.03,.075,.005),
           "narrow":(.03,.035,.01),"late":(.15,.075,.01),"no_drive":(.03,.075,0.)}
    result={k:evaluate(a,psi,times[0],*v) for k,v in cases.items()}
    # On this reference interpolant aM has one minimum before opening. Use
    # deliberately loose bounds from the all-shell vacuum calculation and
    # the observed support endpoints of the selected nominal drive.
    omega_min=math.hypot(1.,3.5)
    amplification=math.sqrt(math.hypot(1.,3.8)*math.hypot(1.,19.2))/omega_min
    support=result["nominal"]["proper_times"]
    if support[2]-support[0]>=1.2 or max(a(support[0]),a(support[2]))>=2:
        raise ValueError("Nominal drive is outside the bounded contraction/expansion interval")
    amplitude_bound=amplification*2**2*1.9**2*1.2/math.sqrt(2*omega_min)
    necessary=math.sqrt(TARGET)/amplitude_bound
    args.output.write_text(json.dumps({"nominal_necessary_amplitude_bound":necessary,
        "nominal_sqrt_action_per_amplitude_upper":amplitude_bound,"history_sha256":PRIOR["DIGEST"],"mL":1,
        "curvature_coupling":"1/6","comoving_number_target_over_S0_volume":TARGET,"cases":result},
        indent=2,sort_keys=True,allow_nan=False)+"\n")

if __name__=="__main__": main()
