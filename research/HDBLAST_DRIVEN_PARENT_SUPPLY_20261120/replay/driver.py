"""Energy-depleting finite scalar driver on a short archived interval of expansion."""
import argparse
from bisect import bisect_right
import json
import math
from pathlib import Path
import runpy

HERE=Path(__file__).resolve().parent
trace=runpy.run_path(str(HERE/'trace.py'))
step,cubic=trace['step'],trace['cubic']
MASS=2*125*.95
SOURCE_OVER_REFERENCE=2.


def geometry(source):
    rows=source['rows']; times=[r['t'] for r in rows]
    def at(z):
        t=times[0]+z/MASS
        if t>times[-1]: raise ValueError('Pulse is outside the source geometry interval.')
        j=min(len(rows)-2,max(0,bisect_right(times,t)-1))
        l,r=rows[j:j+2]; h=r['t']-l['t']
        a,adot=cubic(l['a'],r['a'],l['a']*l['H'],r['a']*r['H'],h,(t-l['t'])/h)
        return (a/rows[0]['a'])**3,adot/(MASS*a)
    return at


def simulate(source,width,amplitude,hstep,feedback=True):
    at=geometry(source)
    def rhs(z,state):
        x,y,px,py,w,ws,wp=state
        vol,H=at(z)
        up,dup=step(2*x/width)
        down,ddown=step(2*x/width-1)
        minimum=amplitude*(up-down)
        slope=2*amplitude*(dup-ddown)/width
        delta=y-minimum
        return (px,py,(delta*slope if feedback else 0.)-3*H*px,-delta-3*H*py,
                -2*vol*delta*slope*px,3*H*vol*px**2,3*H*vol*(py**2-delta**2))
    def advance(z,y,h):
        k1=rhs(z,y)
        k2=rhs(z+h/2,[u+h*v/2 for u,v in zip(y,k1)])
        k3=rhs(z+h/2,[u+h*v/2 for u,v in zip(y,k2)])
        k4=rhs(z+h,[u+h*v for u,v in zip(y,k3)])
        return [u+h*(a+2*b+2*c+d)/6 for u,a,b,c,d in zip(y,k1,k2,k3,k4)]
    z=0.;state=[0.,0.,1.,0.,0.,0.,0.]
    while z<24:
        new=advance(z,state,hstep)
        if new[0]>=width or (new[0]<0 and state[0]>0):
            target=width if new[0]>=width else 0.
            direction=1 if target else -1
            lo,hi=0.,hstep
            for _ in range(36):
                mid=(lo+hi)/2
                if direction*(advance(z,state,mid)[0]-target)>0:
                    hi=mid
                else:
                    lo=mid
            state=advance(z,state,(lo+hi)/2);z+=(lo+hi)/2;break
        state=new;z+=hstep
    else:
        raise ValueError('No freely oscillating endpoint within the selected interval.')
    vol,H=at(z)
    x,y,px,py,w,ws,wp=state
    parent=vol*(py**2+y**2)
    driver=vol*px**2
    lo=hi=parent
    after=state
    for i in range(round(2*math.pi/hstep)):
        t=z+i*hstep
        after=advance(t,after,hstep)
        if 0.000001<after[0]<width-.000001:
            raise ValueError('The drive resumed in the free comparison period.')
        comoving=at(t+hstep)[0]*(after[1]**2+after[3]**2)
        lo=min(lo,comoving);hi=max(hi,comoving)
    return dict(width=width,amplitude=amplitude,source_reacts=feedback,forward_exit=x>width/2,
                end_z=z,proper_elapsed=z/MASS,H_over_M_at_exit=H,parent_share_of_initial_source=parent,
                source_share_at_exit=driver,transferred_work=w,source_pressure_work=ws,parent_pressure_work=wp,
                parent_action_over_target=parent*SOURCE_OVER_REFERENCE,
                next_period_action_over_target_range=[lo*SOURCE_OVER_REFERENCE,hi*SOURCE_OVER_REFERENCE],
                parent_ledger=parent-w+wp,total_ledger=parent+driver-1+ws+wp)


def calculate(source):
    results={label:simulate(source,W,G,.005,feed) for label,W,G,feed in [
        ('primary',4.,.5,True),('smaller_displacement',4.,.25,True),('shorter_region',2.,.5,True),
        ('longer_region',8.,.5,True),('high_transfer',2.,.75,True),('zero_displacement',4.,0.,True),
        ('prescribed_source_control',4.,.5,False)]}
    a0=source['rows'][0]['a']
    target_energy=source['initial_radius']**3*source['initial_parent_density']
    v0=math.sqrt(2*target_energy*SOURCE_OVER_REFERENCE/a0**3)
    fine=simulate(source,4.,.5,.0025)
    reference=results['primary']
    err=abs(fine['parent_action_over_target']-reference['parent_action_over_target'])
    if err>2e-7 or results['zero_displacement']['parent_action_over_target']!=0.:
        raise ValueError('Convergence or zero-coupling control failed.')
    if max(abs(r['total_ledger']) for r in results.values() if r['source_reacts'])>1e-7:
        raise ValueError('Closed source and parent energy balance failed.')
    return dict(parameters=dict(source_velocity_initial=v0,initial_scale_factor=a0,
                                source_comoving_energy=SOURCE_OVER_REFERENCE*target_energy,
                                reference_comoving_energy=target_energy,mass=MASS,g=.5*v0/MASS,
                                field_width=4*v0/MASS,action_target=target_energy/MASS),
                primary_action_ratio_step_refinement=err,runs=results)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--geometry',type=Path,default=HERE/'geometry.json')
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    args.output.write_text(json.dumps(calculate(json.loads(args.geometry.read_text())),indent=2,
                                      sort_keys=True,allow_nan=False)+'\n')
