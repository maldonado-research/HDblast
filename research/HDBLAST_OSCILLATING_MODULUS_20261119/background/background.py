"""Early scalar, radiation and spatial curvature dynamics around the archived minima."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp, simpson
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

PINS={"blast_common.py":"548827e51f7739aa0511a114123c19b9769785d9bddad8936c72e7735057565a",
      "BLAST_REHEATING_RESULTS.json":"44f30d1cd987e4df3ca9889338e1dbd5b0e35dcf42ed7228857dbaeab41c7945"}
DETUNING=1e-3
FLAT_END=65.
FLAT_ANCHOR=50.
MODES=("off","always","harmonic")

def source(archive):
    for name,sha in PINS.items():
        if hashlib.sha256((archive/name).read_bytes()).hexdigest()!=sha:
            raise ValueError("Unexpected archive file: "+name)
    spec=importlib.util.spec_from_file_location("modulus_archive_common",archive/"blast_common.py")
    mod=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod,json.loads((archive/"BLAST_REHEATING_RESULTS.json").read_text())["B"]

class Potential:
    """Use the archived phi and Planck factor tables with a consistent differentiable potential."""
    def __init__(self,mod,tables,minimum,record):
        self.pm=minimum
        self.b=mod.c_star+2/minimum
        self.record=record
        V=lambda p:(1-p/minimum)**2*(1+self.b*p)
        dV=lambda p:(-2/minimum)*(1-p/minimum)*(1+self.b*p)+self.b*(1-p/minimum)**2
        L=mod.Landscape(V,dV,*tables)
        self.v0=float(L.v0)
        # The registered trajectories stay on the smooth tanh branch, far from the endpoints.
        use=(L.branch==0)&(abs(L.phi)<.98)
        self.theta=L.Th[use]
        self.spline=CubicSpline(self.theta,np.column_stack((L.phi[use],L.f[use])),axis=0,extrapolate=False)
        self.coeff=self.spline.c
        self.tmin=float(brentq(lambda th:self.spline(th)[0]-minimum,self.theta[0],self.theta[-1]))
        self.sign=float(-np.sign(minimum))
        ph,f=self.spline(self.tmin)
        dph=self.spline(self.tmin,1)[0]
        self.mharm=float(np.sqrt(3*2*(1+self.b*ph)*dph**2/(minimum**2*f**2*self.v0)))
        self.Hpl=float(np.sqrt(DETUNING*self.v0/3))
        self.gamma=float(record["alpha"]**2*record["m_over_H_top"]**3*self.Hpl**2/(32*np.pi))
        self.growing=record['k0_dynamics']['s_linear']
        self.dq=self.Hpl/(2*np.pi)
        self.barriers=[float(brentq(lambda th:self.spline(th)[0]-p,self.theta[0],self.theta[-1]))
                       for p in (minimum-.1,minimum+.1)]
        self.capture_energy=min(self.at(t)[0] for t in self.barriers)

    def at(self,th):
        ix=int(np.searchsorted(self.theta,th,side="right")-1)
        if ix<0 or ix>=len(self.theta)-1:raise ValueError("Potential table boundary crossed")
        delta=th-self.theta[ix]
        c=self.coeff[:,ix]
        (p,f)=((c[0]*delta+c[1])*delta+c[2])*delta+c[3]
        (dp,df)=(3*c[0]*delta+2*c[1])*delta+c[2]
        z=1-p/self.pm
        weight=1+self.b*p
        u=z*z*weight/(f*f*self.v0)
        du=(z*(-2*weight/self.pm+z*self.b)*dp/(f*f*self.v0)
            -2*u*df/f)
        return u,du,p

def rhs(pot,curvature,gamma):
    def evolve(x,state):
        th,v,ln_a,h,r,source4,work=state
        a=np.exp(ln_a)
        u,du,_=pot.at(th)
        k=v*v/6
        transfer=2*gamma*k
        return [v,-(3*h+gamma)*v-3*du,h,curvature/a**2-3*k-2*r,
                transfer-4*h*r, a**4*transfer,3*h*a**3*(k-u+r/3)]
    return evolve

def seed(pot,curvature,displacement=None):
    d=pot.dq if displacement is None else displacement
    th=pot.sign*d
    v=pot.growing*th if curvature==0 else 0.
    energy=pot.at(th)[0]+v*v/6
    na=0. if curvature==0 else -np.log(energy)/2
    h=np.sqrt(energy) if curvature==0 else 0.
    return np.array([th,v,na,h,0.,0.,0.])

def segment(pot,curvature,gamma,state,t0,t1,rtol):
    def cross(x,y):return y[0]-pot.tmin
    def turn(x,y):return y[1]
    def max_a(x,y):return y[3]
    max_a.terminal=True
    max_a.direction=-1
    return solve_ivp(rhs(pot,curvature,gamma),(t0,t1),state,events=(cross,turn,max_a),
                     dense_output=True,method='DOP853',rtol=rtol,
                     atol=[rtol*1e-2,rtol*1e-2,rtol*1e-2,rtol*1e-2,
                           rtol*1e-6,rtol*1e-6,rtol*1e-2],
                     max_step=2*np.pi/pot.mharm/14)

def cycles(pot,sol):
    returns=[float(t) for t in sol.t_events[0] if pot.sign*sol.sol(t)[1]>0]
    data=[]
    for t0,t1 in zip(returns,returns[1:]):
        z=np.linspace(t0,t1,129)
        th,v,loga,h,r,_,_=sol.sol(z)
        u=np.array([pot.at(x)[0] for x in th])
        kin=v*v/6
        dwell=t1-t0
        # Keep endpoint cancellation out of the harmonicity criterion.
        ratio=np.ones(len(z))
        dx=th-pot.tmin
        away=abs(dx)>1e-8
        ratio[away]=u[away]/(pot.mharm**2*dx[away]**2/6)
        turns=[float(t) for t in sol.t_events[1] if t0<t<t1]
        if len(turns)!=2:raise ArithmeticError("Expected two oscillator turns per cycle")
        for t in turns:
            tht=sol.sol(t)[0]
            ratio=np.append(ratio,pot.at(tht)[0]/(pot.mharm**2*(tht-pot.tmin)**2/6))
        max_dev=float(np.max(abs(ratio-1)))
        max_h=float(np.max(abs(h))/pot.mharm)
        # Ratio of integrated pressures and integrated energy weights dense phases physically.
        pressure=float(simpson(kin-u,x=z))
        energy=float(simpson(kin+u,x=z))
        data.append(dict(start=t0,end=t1,period=dwell,w_time_weighted=pressure/energy,
                         H_over_m=max_h, max_fractional_potential_deviation=max_dev,
                         harmonic=(max_dev<=.1 and max_h<=.01),
                         min_field=float(min(pot.at(sol.sol(t)[0])[2] for t in turns)),
                         max_field=float(max(pot.at(sol.sol(t)[0])[2] for t in turns))))
    return data

def activation(cyc):
    for a,b in zip(cyc,cyc[1:]):
        if a['harmonic'] and b['harmonic']:return b['end']
    return None

def point(pot,k,x,s):
    th,v,ln_a,h,r,src,w=map(float,s)
    a=float(np.exp(ln_a))
    u,_,p=pot.at(th)
    ec=v*v/6+u
    return dict(x=float(x), theta=th, phi=float(p), v=v, a=a, h=h,
                scalar=float(ec), kinetic=v*v/6, potential=float(u), radiation=r, curvature=k/a**2,
                h_prime=k/a**2-v*v/2-2*r,
                radiation_fraction=float(r/(r+ec)), friedmann_residual=float(h*h-ec-r+k/a**2),
                scalar_comoving=a**3*ec, radiation_source_comoving4=src,
                pressure_work=w)

def branch(pot,k,mode,base,stage,rtol,displacement=None):
    s0=seed(pot,k,displacement)
    t_end=FLAT_END if not k else 20.
    if mode=="always":
        segs=[segment(pot,k,pot.gamma,s0,0,t_end,rtol)];gamma_start=0.
    elif mode=="harmonic" and stage is not None:
        switched=base.sol(stage)
        # Restart exactly on the recorded crossing so it belongs to both adjoining
        # complete-cycle intervals independently of event-root roundoff.
        switched[0]=pot.tmin
        segs=[base,segment(pot,k,pot.gamma,switched,stage,t_end,rtol)]
        gamma_start=stage
    else:
        segs=[base];gamma_start=None
    bounds=[(0.,stage),(stage,segs[-1].t[-1])] if len(segs)>1 else [(0.,segs[0].t[-1])]
    first,last=segs[0],segs[-1]
    if not all(s.success for s in segs):raise ValueError("ODE failed")
    end=float(last.t[-1]);final=point(pot,k,end,last.y[:,-1])
    records=[]
    cross=[]
    turn=[]
    for i,(s,(lo,hi)) in enumerate(zip(segs,bounds)):
        records.extend(c for c in cycles(pot,s) if c['start']>=lo-1e-9 and c['end']<=hi+1e-9)
        cross.extend((t,s.sol(t)) for t in s.t_events[0] if lo+1e-8<t<=hi+1e-8)
        turn.extend((t,s.sol(t)) for t in s.t_events[1] if lo+1e-8<t<=hi+1e-8)
    cross.sort(key=lambda q:q[0]);turn.sort(key=lambda q:q[0])
    anchor=None
    if not k:
        c=next(c for c in records if c["start"]>=FLAT_ANCHOR)
        stop=c["end"];y=next(s.sol(stop) for s,(lo,hi) in zip(segs,bounds) if lo<=stop<=hi)
        anchor=point(pot,k,stop,y)
        anchor['cycle']=c
        anchor['gamma_over_Htop']=pot.gamma if mode!='off' else 0.
    def state(t):
        return next(s.sol(t) for s,(lo,hi) in zip(segs,bounds) if lo-1e-12<=t<=hi+1e-12)
    # Once energy drops below both sides of an interval surrounding the minimum,
    # an expanding solution cannot escape it while the radiation source is passive.
    capture=None
    for t,sc in cross:
        if sc[1]**2/6<pot.capture_energy and sc[3]>0:
            capture=float(t);break
    audit=[]
    for s,(lo,hi) in zip(segs,bounds):
        times=np.linspace(lo,hi,max(5,int(160*(hi-lo))))
        y=s.sol(times);u=np.array([pot.at(x)[0] for x in y[0]])
        a=np.exp(y[2]);kin=y[1]**2/6
        norm=np.max(abs(y[3]**2-(kin+u+y[4]-k/a**2)))
        init=point(pot,k,0,s0)
        work=np.max(abs(a**3*(kin+u+y[4])+y[6]-init['scalar']*init['a']**3))
        rad=np.max(abs(a**4*y[4]-y[5]))
        audit.append([norm,work,rad])
    short={}
    for t in (12.,):
        if t<=end:short[str(t)]=point(pot,k,t,state(t))
    audit_max=np.max(audit,axis=0)
    if audit_max[0]>2e-8 or audit_max[1]>3e-6 or audit_max[2]>3e-9:
        raise ArithmeticError("Friedmann or integrated energy balance failed")
    if k and (final["h_prime"]>=0 or abs(final["h"])>1e-10):
        raise ArithmeticError("Friedmann branch did not reach a decreasing zero")
    return dict(mode=mode,curvature=k,seed=point(pot,k,0,s0),effective_width_switch_x=gamma_start,
                first_minimum_crossing=(point(pot,k,*cross[0]) if cross else None),
                first_scalar_turn=(point(pot,k,*turn[0]) if turn else None),
                minimum_crossings=len(cross),completed_cycles=len(records),
                first_harmonic_cycle=next((c for c in records if c['harmonic']),None),
                capture_sufficient_at_minimum_crossing_x=capture,
                cycle_first=(records[0] if records else None),
                cycle_last=(records[-1] if records else None),
                flat_anchor=anchor,intermediate=short,endpoint=final,
                turnaround=(final if len(last.t_events[2]) else None),
                checks=dict(zip(("max_abs_friedmann","max_abs_comoving_work","max_abs_radiation_source"),
                                map(float,np.max(audit,axis=0)))))

def experiment(archive,rtol=2e-11):
    mod,ref=source(archive)
    tables=mod.build_tables()
    out={}
    for pm in (.5,-.5):
        key=f"phi_m={pm:g}";record=ref[key]
        p=Potential(mod,tables,pm,record)
        row=dict(pm=pm,gamma=p.gamma,Htop_over_Mpl=p.Hpl,t_minimum=p.tmin,
                 mass_exact_interpolant=p.mharm, mass_archived=record["m_over_H_top"],
                 width_model="one real massless minimally coupled scalar",capture_interval_phi=[pm-.1,pm+.1],
                 capture_threshold=p.capture_energy)
        for k in (0,1):
            seed0=seed(p,k)
            base=segment(p,k,0.,seed0,0.,FLAT_END if not k else 20.,rtol)
            a=activation(cycles(p,base))
            runs={mode:branch(p,k,mode,base,a,rtol) for mode in MODES}
            row['flat' if not k else 'closed']=runs
        demo=segment(p,1,0.,seed(p,1,.08),0,20.,rtol)
        demo_t=demo.t[-1]
        row['archived_large_displacement_demo']=point(p,1,demo_t,demo.y[:,-1])
        out[key]=row
    return dict(detuning=DETUNING,units="x=Htop t; densities divided by Vtop=3 Mpl^2 Htop^2",
                primary_width="harmonic",branches=out)

def independent(pot,k,mode,reference):
    """Integrate the same data with a constrained flat rate or the closed acceleration equation."""
    s0=seed(pot,k)
    a0=float(np.exp(s0[2]))
    initial=np.array([s0[0],s0[1],a0,a0*s0[3],0.])
    trigger=reference["effective_width_switch_x"]
    end=reference["flat_anchor"]["x"] if not k else reference["endpoint"]["x"]+.02
    def run(lo,hi,y,gamma):
        def rhs_alt(x,z):
            th,v,a,ad,r=z
            u,du,_=pot.at(th)
            h=ad/a if k else np.sqrt(max(v*v/6+u+r,0.))
            return (v,-(3*h+gamma)*v-3*du,a*h,a*(u-v*v/3-r),gamma*v*v/3-4*h*r)
        def stop(x,z):return z[3]
        stop.direction=-1
        stop.terminal=True
        return solve_ivp(rhs_alt,(lo,hi),y,dense_output=True,events=stop if k else None,
                         method="DOP853",rtol=2e-12,atol=[3e-14,3e-14,3e-14,3e-14,5e-18],
                         max_step=2*np.pi/pot.mharm/28)
    lo=0.
    gamma=pot.gamma if mode=="always" else 0.
    if mode=="harmonic" and trigger is not None:
        first=run(0.,trigger,initial,0.)
        initial=first.y[:,-1]
        gamma=pot.gamma
        lo=trigger
    out=run(lo,end,initial,gamma)
    if not out.success:raise ArithmeticError("Independent evolution failed")
    t=out.t[-1]
    th,v,a,ad,r=out.y[:,-1]
    E=v*v/6+pot.at(th)[0]
    h=ad/a if k else np.sqrt(E+r)
    baseline=reference["turnaround"] if k else reference["flat_anchor"]
    values=dict(x=float(t),a=float(a),h=float(h),scalar=float(E),radiation=float(r))
    errors={key:val-baseline[key] for key,val in values.items()}
    if max(abs(errors[key]) for key in ("x","a","h","scalar"))>4e-5 or abs(errors["radiation"])>2e-9:
        raise ArithmeticError("Independent expansion branch disagrees")
    return dict(independent=values,difference_from_primary=errors)

def after_maximum(pot,curvature_run,always):
    """Follow the post-maximum EFT only until a flagged finite boundary is met."""
    m=curvature_run["turnaround"]
    initial=np.array([m["theta"],m["v"],np.log(m["a"]),0.,m["radiation"],0.,0.])
    initial_a=curvature_run["seed"]["a"]
    lo_theta,hi_theta=[brentq(lambda t:pot.at(t)[2]-field,
                              pot.theta[1],pot.theta[-2]) for field in (.9,-.9)]
    def small(x,y):return np.exp(y[2])-.5*initial_a
    def edge_plus(x,y):return y[0]-lo_theta
    def edge_minus(x,y):return hi_theta-y[0]
    for event in (small,edge_plus,edge_minus):
        event.direction=-1
        event.terminal=True
    def bounce(x,y):return y[3]
    bounce.direction=1
    def recollapse(x,y):return y[3]
    recollapse.direction=-1
    sol=solve_ivp(rhs(pot,1,pot.gamma if always else 0.),(m["x"],m["x"]+12.),
                  initial,method="DOP853",rtol=1e-11,atol=1e-14,
                  max_step=2*np.pi/pot.mharm/128,dense_output=True,
                  events=(small,edge_plus,edge_minus,bounce,recollapse))
    if not sol.success:raise ArithmeticError(sol.message)
    reason=next((name for name,times in zip(("a_half_initial","phi_plus_0.9","phi_minus_0.9"),
                                             sol.t_events[:3]) if len(times)),"time_limit")
    endpoint=point(pot,1,sol.t[-1],sol.y[:,-1])
    signs={name:[point(pot,1,t,sol.sol(t)) for t in events if t>m["x"]+1e-6]
           for name,events in zip(("reexpansions","additional_maxima"),sol.t_events[3:])}
    return dict(first_maximum_x=m["x"],stop_reason=reason,endpoint=endpoint,**signs)

def post_max_check(pot,first,always,comparison):
    """Recheck extrema using the scale factor and its velocity instead of H."""
    item=first["turnaround"]
    a0=first["seed"]["a"]
    gamma=pot.gamma if always else 0.
    def rate(t,z):
        th,v,a,ad,r=z
        u,du,_=pot.at(th)
        h=ad/a
        return (v,-(3*h+gamma)*v-3*du,ad,a*(u-v*v/3-r),gamma*v*v/3-4*h*r)
    def half(t,z):return z[2]-.5*a0
    def field(t,z):return .9-abs(pot.at(z[0])[2])
    def extremum(t,z):return z[3]
    half.terminal=field.terminal=True
    half.direction=field.direction=-1
    seed0=[item["theta"],item["v"],item["a"],0.,item["radiation"]]
    sol=solve_ivp(rate,(item["x"],item["x"]+12),seed0,dense_output=True,
                  events=(half,field,extremum),method="DOP853",rtol=1e-12,atol=2e-16,
                  max_step=2*np.pi/pot.mharm/200)
    if not sol.success:raise ArithmeticError("Scale-factor continuation failed")
    observed=[(float(t),float(sol.sol(t)[2])) for t in sol.t_events[2] if t>item["x"]+1e-5]
    expected=sorted((p["x"],p["a"]) for p in comparison["reexpansions"]+comparison["additional_maxima"])
    if len(observed)!=len(expected):raise ArithmeticError("Post-maximum sign changes disagree")
    abs_t=max([abs(t0-t1) for (t0,_),(t1,_) in zip(observed,expected)]+[0.])
    abs_a=max([abs(a0-a1) for (_,a0),(_,a1) in zip(observed,expected)]+[0.])
    dif_end=sol.t[-1]-comparison["endpoint"]["x"]
    if max(abs_t,abs_a,abs(dif_end))>2e-6:raise ArithmeticError("Post-maximum extrema disagree")
    return dict(scale_factor_zero_velocity_count=len(observed),max_abs_x=abs_t,max_abs_a=abs_a,
                end_time_difference=float(dif_end))

def checks(archive, results):
    mod,old=source(archive)
    tables=mod.build_tables()
    rows={}
    for key,row in results["branches"].items():
        pot=Potential(mod,tables,row["pm"],old[key])
        orbit=old[key]["k0_dynamics"]
        off=row["flat"]["off"]
        cross=off["first_minimum_crossing"]
        older=old[key]["closed_demo(dTheta0=0.08)"]
        demo=row["archived_large_displacement_demo"]
        comparison=dict(
            flat_first_cross_time=cross["x"]-orbit["tau_first_crossing"],
            flat_first_cross_efolds=np.log(cross["a"])-orbit["N_E_first_crossing_of_minimum"],
            flat_first_cross_kinetic=cross["kinetic"]-orbit["K_over_V0_at_crossing"],
            flat_a_at_12=off["intermediate"]["12.0"]["a"]-orbit["a_end"],
            flat_scalar_comoving_at_12=off["intermediate"]["12.0"]["scalar_comoving"]-orbit["rho_a3_late_over_V0"],
            closed_large_displacement_turn_time=demo["x"]-older["tau_max_expansion"],
            closed_large_displacement_turn_scale=demo["a"]-older["a_max"])
        if abs(comparison["flat_first_cross_time"])>3e-4 or abs(comparison["flat_a_at_12"])>.01:
            raise ValueError("Archived flat result changed")
        if abs(comparison["closed_large_displacement_turn_time"])>.0012:
            raise ValueError("Archived closed large-displacement result changed")
        rows[key]=dict(archive_numerical_differences=comparison,
            archive_dust_estimate=old[key]["closed_quantum_start(t=1e-3)"]["a_max_if_dust"],
            measured_undamped_first_maximum=row["closed"]["off"]["turnaround"]["a"],
            independent={geometry:{mode:independent(pot,k,mode,row[geometry][mode])
                                      for mode in ("always","harmonic")}
                         for geometry,k in (("flat",0),("closed",1))},
            closed_after_first_maximum={mode:after_maximum(pot,row["closed"][mode],mode=="always")
                                         for mode in ("always","off")})
        for mode,run in rows[key]["closed_after_first_maximum"].items():
            run["independent_scale_factor_check"]=post_max_check(
                pot,row["closed"][mode],mode=="always",run)
    return dict(branches=rows)

if __name__=="__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--archive",type=Path,required=True)
    cli.add_argument("--output",type=Path,required=True)
    cli.add_argument("--checks-output",type=Path)
    cli.add_argument("--rtol",type=float,default=2e-11)
    args=cli.parse_args()
    data=experiment(args.archive,args.rtol)
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+"\n")
    if args.checks_output:
        args.checks_output.write_text(json.dumps(checks(args.archive,data),indent=2,allow_nan=False)+"\n")
