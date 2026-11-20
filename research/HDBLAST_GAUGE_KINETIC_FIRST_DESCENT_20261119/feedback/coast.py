"""Source-free continuation of completed gauge switches until re-entry or turnaround."""
import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq


def coast(branch, switch_end_record, first_minimum_record=None):
    """Continue one coupled trajectory in absolute proper x from its outgoing 80% crossing."""
    if switch_end_record is None:
        return dict(status='switch_not_reached', peak=None, stop=None, first_minimum=None,
                    half_intervals=[], total_efolds_at_or_above_half=0.)

    start=switch_end_record
    xs,qs,vs,as_,hs,rs=(float(start[k]) for k in ('x','q','v','a','h','photon_density'))
    qi=branch.fitted.sign*branch.fitted.dq
    distance=branch.fitted.tmin-branch.top-qi
    progress=lambda q:(q-qi)/distance
    if abs(progress(qs)-.8)>1.e-8 or as_<=0. or rs<0.:
        raise ValueError('Expected a physical outgoing switch endpoint')
    u0,force0=branch.potential(qs)
    scalar0=u0+vs*vs/6
    if abs(hs*hs+1/as_**2-rs-scalar0)>5e-8:
        raise ValueError('Switch endpoint does not satisfy the coupled Friedmann constraint')
    hprime=1/as_**2-vs**2/2-2*rs
    domain=lambda q:min(branch.top+q-branch.fitted.theta[0],
                        branch.fitted.theta[-1]-branch.top-q)-.02

    def trajectory(direct):
        def fields(y):
            q,v,third,fourth,fifth,_=y
            if direct:
                return q,v,third,fourth/third,fifth
            return q,v,as_*np.exp(third),fourth,rs*np.exp(-4*third)

        def rhs(t,y):
            q,v,a,h,r=fields(y)
            u,force=branch.potential(q)
            if direct:
                return v,force-3*h*v,h*a,a*(u-v*v/3-r),-4*h*r,h*v*v
            return v,force-3*h*v,h,1/a**2-v*v/2-2*r,0.,h*v*v

        def radius_max(t,y): return y[3]
        radius_max.direction=-1
        radius_max.terminal=True
        def reentry(t,y):
            # The input is the outgoing 0.8 event, with a numerically refined q.
            # Subtract that same q to avoid a fictitious immediately adjacent return.
            return (y[0]-qs)/distance
        reentry.direction=-1
        reentry.terminal=True
        def limit(t,y): return domain(y[0])
        limit.direction=-1
        limit.terminal=True
        def passage(t,y): return progress(y[0])-1.
        passage.direction=1
        def stationary(t,y):
            u,_=branch.potential(y[0])
            return y[1]*y[1]/6-2*u
        initial=[qs,vs,as_ if direct else 0.,as_*hs if direct else hs,rs if direct else 0.,0.]
        labels=['first_radius_maximum','first_switch_reentry','calibrated_range_stop']
        events=(radius_max,reentry,limit,passage,stationary)
        if hs< -1e-12 or (abs(hs)<1e-12 and hprime<0):
            status=labels[0] if abs(hs)<1e-12 else 'already_contracting'
        elif vs/distance<0 or (vs==0 and force0/distance<=0):
            status=labels[1]
        elif domain(qs)<=0:
            status=labels[2]
        else:
            status=None
        start_x=xs
        states=[]
        extrema=[]
        minima=[]
        y=initial
        while status is None:
            stop_x=min(xs+256.,start_x+max(2.,start_x-xs))
            sol=solve_ivp(rhs,(start_x,stop_x),y,method='RK45' if direct else 'DOP853',
                          rtol=2e-12 if direct else 2e-11,atol=2e-14 if direct else 2e-13,
                          max_step=.003 if direct else .012,dense_output=True,events=events)
            if not sol.success:
                raise ArithmeticError(sol.message)
            states.append(sol)
            extrema.extend(sol.t_events[4])
            minima.extend(sol.t_events[3])
            start_x=float(sol.t[-1])
            y=sol.y[:,-1]
            status=next((name for name,t in zip(labels,sol.t_events[:3]) if len(t)),None)
            if status is None and start_x>=xs+256.:
                status='coast_horizon'
        end=start_x
        def values(t):
            for sol in states:
                if t<=sol.t[-1]:
                    return sol.sol(t)
            return np.asarray(y)
        def data(t):
            yy=values(t)
            q,v,a,h,r=fields(yy)
            u,_=branch.potential(q)
            e=u+v*v/6
            return dict(x=float(t),a=float(a),q=float(q),v=float(v),h=float(h),photon_density=float(r),
                        parent_density=float(e),fraction_of_matter=float(r/(r+e)),progress=float(progress(q)))
        def half(t):
            yy=values(t)
            q,v,_,_,r=fields(yy)
            u,_=branch.potential(q)
            return float(r-u-v*v/6)
        # During expansion F' has the sign of K-2u whenever 0<F<1.
        # Each adjacent pair of these extrema therefore contains at most one half crossing.
        knots=sorted(set([xs,end,*extrema]))
        peak=max((data(t) for t in knots),key=lambda s:s['fraction_of_matter'])
        half_events=[float(t) for t in knots if half(t)==0.]
        for lo,hi in zip(knots,knots[1:]):
            if half(lo)*half(hi)<0.:
                half_events.append(float(brentq(half,lo,hi,xtol=2e-13)))
        half_events=sorted(set(half_events))
        boundaries=sorted(set([xs,end,*half_events]))
        intervals=[]
        for lo,hi in zip(boundaries,boundaries[1:]):
            if half((lo+hi)/2)>=0.:
                begin,finish=data(lo),data(hi)
                intervals.append(dict(begin_x=lo,end_x=hi,proper_duration=hi-lo,
                                      begin_a=begin['a'],end_a=finish['a'],
                                      delta_N=math.log(finish['a']/begin['a'])))
        samples=[s for sol in states for s in sol.y.T] or [np.asarray(initial)]
        q,v,a,h,r=fields(np.array(samples).T)
        u,_=branch.potential(q)
        balance=np.abs(h*h+1/a**2-r-u-v*v/6)
        e_balance=np.abs(u+v*v/6-scalar0+np.array(samples)[:,-1])
        comoving=np.abs(r*(a/as_)**4-rs)
        return dict(status=status,start=data(xs),stop=data(end),peak=peak,
                    first_minimum=data(minima[0]) if minima else None,
                    half_crossings=half_events,half_intervals=intervals,
                    total_efolds_at_or_above_half=sum(d['delta_N'] for d in intervals),
                    total_proper_duration_at_or_above_half=sum(d['proper_duration'] for d in intervals),
                    constraint_error=float(max(balance)),parent_balance_error=float(max(e_balance)),
                    normalized_comoving_density_error=float(max(comoving)))
    result=trajectory(False)
    alt=trajectory(True)
    if alt['status']!=result['status'] or len(alt['half_intervals'])!=len(result['half_intervals']):
        raise ArithmeticError('Source-free formulations disagree on stopping or equality events')
    errors=[abs(result['total_efolds_at_or_above_half']-alt['total_efolds_at_or_above_half'])]
    for stage in ('stop','peak','first_minimum'):
        if (result[stage] is None)!=(alt[stage] is None):
            raise ArithmeticError('Source-free formulations disagree on first passage')
        if result[stage] is not None:
            errors.extend(abs(result[stage][key]-alt[stage][key]) for key in
                          ('x','q','v','a','h','fraction_of_matter'))
    for t,s in zip(result['half_crossings'],alt['half_crossings']):
        errors.append(abs(t-s))
    result['independent_maximum_difference']=max(errors)
    other=first_minimum_record
    if other is not None:
        if result['first_minimum'] is None:
            raise ArithmeticError('Source-free continuation missed archived first passage')
        result['feedback_first_minimum_difference']=max(abs(float(other[k])-result['first_minimum'][k])
                    for k in ('x','q','v','a','h','photon_density','fraction_of_matter'))
    if max(errors+[result['constraint_error'],result['parent_balance_error'],
                   result['normalized_comoving_density_error'],
                   result.get('feedback_first_minimum_difference',0.)])>5e-8:
        raise ArithmeticError('Source-free trajectory or feedback matching tolerance exceeded')
    return result
