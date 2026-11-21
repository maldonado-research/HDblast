"""Bound the change in a homogeneous axial electric field for any subluminal particle current."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


K = 2.54
Q = 1.e7
TARGET_TIMES = (6.27854515, 6.30944143, 6.37, 6.45, 6.5, 6.57, 6.67)


def extrema(t, a, d, v):
    """Integrate the two extremizing histories with a midpoint step, in physical wall time."""
    # In each cell the two forcings use the midpoint of the endpoint charge currents.
    signed_source = Q*d[:,None] * (v[:,None]-np.array([-1., 1.]))/(K*a[:,None])
    dt = np.diff(t)
    w = 2 / (0.5*(a[:-1]+a[1:]))
    damping = 2/K
    force = (signed_source[1:] + signed_source[:-1]) / 2
    # Implicit midpoint for the forced damped oscillator. Solve the two-by-two system.
    h = dt/2
    wheel = h*w
    denominator = 1+h*damping+wheel**2
    fs, gs = np.zeros((2,len(t),2))
    for i in range(len(dt)):
        rhs_f = (1-h[i]*damping)*fs[i] - wheel[i]*gs[i] + dt[i]*force[i]
        rhs_g = wheel[i]*fs[i]+gs[i]
        fs[i+1] = (rhs_f-wheel[i]*rhs_g)/denominator[i]
        gs[i+1] = rhs_g+wheel[i]*fs[i+1]
    return fs, gs


def cell_enclosure(t, a, d, u):
    """Use upper positive quadratures over each interpolated cell to enclose the electrical change."""
    vel = u/np.hypot(1,u)
    max_d = np.maximum(d[:-1],d[1:])
    min_a = np.minimum(a[:-1],a[1:])
    vhi = np.maximum(vel[:-1],vel[1:])
    vlo = np.minimum(vel[:-1],vel[1:])
    drive = (Q*max_d/(K*min_a))[:,None]*np.column_stack((1+vhi,1-vlo))
    dt = np.diff(t)
    damping = 2/K
    decays = np.exp(-damping*dt)
    contribution = (-np.expm1(-damping*dt)/damping)[:,None] * drive
    bound = np.zeros((len(t),2))
    for i in range(len(dt)):
        bound[i+1] = decays[i]*bound[i]+contribution[i]
    return bound


def crossing(t, data):
    hits=np.flatnonzero((data[:-1]<0)&(data[1:]>=0))
    if not hits.size:
        return None
    k=int(hits[0])
    return float(t[k]+(t[k+1]-t[k])*(-data[k])/(data[k+1]-data[k]))


def solve_one(src):
    with np.load(src) as data:
        t, a, d, u, imposed, gauge = (data[key] for key in ("t","a","D","u","X","K"))
    tmax=TARGET_TIMES[-1]
    if t[-1]<tmax or t[0]>TARGET_TIMES[0] or np.any(np.diff(t)<=0) or d[0]!=0 or min(a)<=0:
        raise ValueError("Field history does not span the birth and measurement times")
    keep=np.searchsorted(t,tmax,side="right")
    times=np.unique(np.r_[t[:keep], TARGET_TIMES])
    a,d,u,imposed,gauge=(np.interp(times,t,x) for x in (a,d,u,imposed,gauge))
    t=times
    if np.max(abs(gauge-K))>1.e-10:
        raise ValueError("Changing gauge coefficient in the optical interval")
    age=t[-1]-t[0]
    minimum_kernel=np.exp(-2*age/K) - 2*(age/min(a))**2
    if minimum_kernel<=0.:
        raise ValueError("Sign of the retarded electric response is undetermined")
    v=u/np.hypot(1,u)
    bound=cell_enclosure(t,a,d,u)
    fs, gs=extrema(t,a,d,v)
    best=imposed+fs[:,0]/a**2
    worst=imposed+fs[:,1]/a**2
    outside=np.max(np.column_stack((fs[:,0]-bound[:,0],-fs[:,1]-bound[:,1])))
    if outside>5e-8 or np.min(fs[:,0])<0 or np.max(fs[:,1])>0:
        raise ValueError("Extremal histories left their analytical enclosure")
    # A directly integrated 4-component ODE checks sign, splitting, and both field components.
    def rhs(x,state):
        radius=np.interp(x,t,a)
        velocity=np.interp(x,t,v)
        rest=np.interp(x,t,d)
        f=state[:2];g=state[2:]
        df=-2*f/K-2*g/radius+Q*rest*(velocity-np.array([-1.,1.]))/(K*radius)
        return np.r_[df,2*f/radius]
    ode=solve_ivp(rhs,(t[0],t[-1]),np.zeros(4),method="DOP853",rtol=5.e-10,atol=1.e-11,
                  t_eval=TARGET_TIMES,max_step=.0008)
    if not ode.success:
        raise ValueError(ode.message)
    idx=np.searchsorted(t,TARGET_TIMES)
    numerical=np.vstack((fs[idx].T,gs[idx].T))
    residual=float(np.max(abs(numerical-ode.y))/max(np.abs(numerical).max(),1e-8))
    if residual>5.e-6:
        raise ValueError("Independent integration failed to match the extremal currents")
    # The maximum of a constant-input scalar damped envelope lies at a cell endpoint.
    whole_cell=np.maximum(imposed[:-1],imposed[1:])+np.maximum(bound[:-1,0],bound[1:,0])/np.minimum(a[:-1],a[1:])**2
    first_possible=np.flatnonzero(whole_cell>=0.)
    rows=[]
    for i in idx:
        rows.append(dict(time=float(t[i]), reference_E=float(imposed[i]), reference_v=float(v[i]),
                     field_change_upper=float(bound[i,0]/a[i]**2),
                     field_change_lower=float(-bound[i,1]/a[i]**2),
                     E_at_negative_saturated_current=float(best[i]), E_at_positive_saturated_current=float(worst[i])))
    return dict(rows=rows,first_birth_left_knot=float(t[0]),kernel_lower=float(minimum_kernel),
                certified_negative_through=float(t[first_possible[0]]) if first_possible.size else float(t[-1]),
                reference_first_zero=crossing(t,imposed),earliest_kinematic_zero=crossing(t,best),
                first_zero_with_positive_saturated_current=crossing(t,worst),
                independent_solver_scaled_error=residual)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--inputs",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    opt=p.parse_args()
    results={key:solve_one(opt.inputs/(key+".npz")) for key in ("dedicated_half","shared_half")}
    opt.output.parent.mkdir(parents=True,exist_ok=True)
    opt.output.write_text(json.dumps(results,sort_keys=True,indent=2,allow_nan=False)+"\n")


if __name__=="__main__":
    main()
