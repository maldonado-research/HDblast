"""Classical spherical Maxwell triplet on two fixed Minkowski annuli, per unit wall volume."""
import argparse
import json
import math
from pathlib import Path
import importlib.util

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.special import expit
from scipy.sparse import coo_matrix, diags


def switch(x):
    if x<=0:
        return 0.
    if x>=1:
        return 1.
    return expit(1/(1-x)-1/x)


def coupling(psi, g, cap):
    """Bound the positive kinetic coefficient on both signs of the field."""
    s = abs(g*psi)
    width = .4
    onset = math.log(cap) - width/2
    x = (s-onset)/width
    if x<=0:
        h, slope = s, 1.
    elif x>=1:
        h, slope = math.log(cap), 0.
    else:
        h = onset + width*(x-2.5*x**4+3*x**5-x**6)
        slope = (1-x)**3*(1+3*x+6*x*x)
    K = math.exp(math.copysign(h,psi))
    return K, K*g*slope


def pulse(width, cap, return_width=None):
    start = -3*width
    end = width*(math.log(cap-1)/4 + .0625)
    def shape(t, main_end=end):
        if t<=start or t>=main_end:
            return 0.
        s = t/width
        return math.sqrt(expit(4*s))*expit(-4*s)*switch((s+3)/.35)*switch((main_end-t)/(.125*width))
    edges=[start+.35*width,end-.125*width]
    if return_width:
        main_end=end
        primary=shape
        area=quad(primary,start,end,points=edges,epsabs=1e-13)[0]
        def lobe(x):
            return math.exp(4-1/(x*(1-x))) if 0<x<1 else 0.
        lobe_area=quad(lobe,0.,1.,epsabs=2e-14)[0]
        def shape(t):
            if t<=main_end:
                return primary(t)
            return -area*lobe((t-main_end)/return_width)/(return_width*lobe_area)
        edges.extend([main_end,main_end+return_width/2])
        end=main_end+return_width
    norm = math.sqrt(2*quad(lambda t:shape(t)**2,start,end,points=edges,epsabs=2e-14,epsrel=3e-13)[0])
    return start, end, lambda t:shape(t)/norm


def mesh(radius, gap, elements):
    """Two oriented copies joined at a single node. x increases left to right."""
    h=gap/elements
    nodes=np.array([-1.,-math.sqrt(3/7),0.,math.sqrt(3/7),1.])
    weights=np.array([1/10,49/90,32/45,49/90,1/10])*h/2
    coords=np.arange(2*elements)[:,None]*h+(nodes+1)*h/2
    r=radius-abs(coords-gap)
    w=r/radius
    ids=4*np.arange(2*elements)[:,None]+np.arange(5)
    n=8*elements+1
    bary=np.array([1/math.prod(x-nodes[j] for j in range(5) if i!=j)
                   for i,x in enumerate(nodes)])
    D=np.zeros((5,5))
    for i,x in enumerate(nodes):
        for j,y in enumerate(nodes):
            if i!=j:
                D[i,j]=bary[j]/(bary[i]*(x-y))
        D[i,i]=-D[i].sum()
    D*=2/h
    mass=np.bincount(ids.ravel(),(w*weights).ravel(),minlength=n)
    angular=np.bincount(ids.ravel(),(4/r**2*w*weights).ravel(),minlength=n)
    blocks=np.einsum('ij,ei,ik->ejk',D,w*weights,D,optimize=True)
    rows=np.broadcast_to(ids[:,:,None],blocks.shape).ravel()
    cols=np.broadcast_to(ids[:,None,:],blocks.shape).ravel()
    grad=coo_matrix((blocks.ravel(),(rows,cols)),shape=(n,n)).tocsr()
    # The separate angular array gives the magnetic energy independently at readout.
    return mass,angular,grad,grad+diags(angular)


def run(radius=16.,width=1.,channel='equal',g=math.sqrt(8),cap=20.,
        gap=1.,elements=12,horizon=18.,tol=3e-10,trace=False,return_width=None):
    if not (radius>gap>0 and elements>=2 and width>0):
        raise ValueError('Need a source in the retained annulus')
    mass,angular,grad,stiff=mesh(radius,gap,elements)
    n=len(mass)
    b=n//2
    w_port=(radius-gap)/radius
    kappa=4/radius**2
    start,end,signal=pulse(width,cap,return_width)
    pair={'equal':(1.,1.),'one':(math.sqrt(2),0.),'odd':(1.,-1.),
          'none':(0.,0.)}[channel]
    factors=np.asarray(pair)/math.sqrt(w_port)
    def dynamics(t,s):
        q=s[:n]; P=s[n:2*n]
        psi, v=s[2*n:2*n+2]
        K,derivative=coupling(psi,g,cap)
        speed=P/mass
        speed[b]=P[b]/(mass[b]+K)
        incident=factors*signal(t)
        port_speed=speed[[0,-1]]
        pdot=-(stiff@q)
        pdot[[0,-1]]+=w_port*(2*incident-port_speed)
        pdot[b]-=K*kappa*q[b]
        return np.r_[speed,pdot,v,.5*derivative*(speed[b]**2-kappa*q[b]**2),
                     w_port*incident**2,w_port*(port_speed-incident)**2,
                     .5*derivative*speed[b]**2,.5*derivative*kappa*q[b]**2]
    def snapshot(t,s):
        q=s[:n];P=s[n:2*n]
        psi,v,inc_l,inc_r,out_l,out_r,j_e,j_m=map(float,s[2*n:])
        K,derivative=coupling(psi,g,cap)
        speeds=P/mass
        speeds[b]=P[b]/(mass[b]+K)
        scalar=.5*v*v
        ge=.5*K*speeds[b]**2
        gm=.5*K*kappa*q[b]**2
        bulk_e=.5*mass@(speeds**2)
        bulk_m=.5*(q@(grad@q))
        bulk_angular=.5*angular@(q**2)
        incoming=inc_l+inc_r
        outgoing=out_l+out_r
        residual=incoming-outgoing-scalar-ge-gm-bulk_e-bulk_m-bulk_angular
        frozen=bool(derivative==0 and ((psi*v>=0 and abs(g*psi)>math.log(cap)+.2) or g==0))
        return dict(t=float(t),psi=psi,scalar_velocity=v,K=K,wall_q=float(q[b]),
                    wall_velocity=float(speeds[b]),scalar=scalar,surface_electric=float(ge),
                    surface_magnetic=float(gm),bulk_electric=float(bulk_e),
                    bulk_radial_gradient=float(bulk_m),bulk_angular=float(bulk_angular),
                    incident_left=inc_l,incident_right=inc_r,outgoing_left=out_l,outgoing_right=out_r,
                    balance=float(residual),electric_impulse=j_e,magnetic_impulse=j_m,scalar_decoupled=frozen)
    late=end+horizon*max(1.,width)+2*gap
    ts=np.unique(np.r_[np.linspace(start,late,181),end,end+gap])
    solution=solve_ivp(dynamics,(start,late),np.zeros(2*n+8),method='DOP853',rtol=tol,atol=tol*1e-3,
                       t_eval=ts,max_step=min(.1*width,1.1*gap/elements))
    if not solution.success:
        raise ArithmeticError(solution.message)
    data=[snapshot(t,s) for t,s in zip(solution.t,solution.y.T)]
    energy=max(abs(row['balance']) for row in data)
    final=data[-1]
    if energy>1e-7 or abs(final['incident_left']+final['incident_right']-(channel!='none'))>1e-8:
        raise ArithmeticError((energy,final))
    pre=[r for r in data if r['t']<start+gap-.1]
    record=dict(radius=radius,source_radius=radius-gap,gap=gap,elements_each_side=elements,
                channel=channel,width=width,g=g,cap=cap,compensating_width=return_width,first_possible_arrival=start+gap,
                max_sampled_prefront_field=max((abs(r['wall_q']) for r in pre),default=0.),
                end_of_port_drive=data[list(ts).index(end)],
                direct_front_after_end=data[list(ts).index(end+gap)],readout=final,
                max_energy_residual=energy,max_surface_magnetic=max(r['surface_magnetic'] for r in data),
                max_sampled_scalar=max(r['scalar'] for r in data),earliest_sampled_decoupling=next(
                    (r['t'] for r in data if r['scalar_decoupled']),None),
                max_impulse_residual=max(abs(r['electric_impulse']-r['magnetic_impulse']
                                             -r['scalar_velocity']) for r in data))
    if trace:
        keys=('t','wall_q','scalar','scalar_velocity','surface_magnetic','electric_impulse',
              'magnetic_impulse','outgoing_left','outgoing_right')
        record['trace']={key:[r[key] for r in data] for key in keys}
    return record


def planar():
    start,end,signal=pulse(1.,20.)
    def rhs(t,state):
        P,psi,v=state
        K,derivative=coupling(psi,math.sqrt(8),20.)
        E=P/K
        return [4*signal(t)-2*E, v, .5*derivative*E*E]
    solution=solve_ivp(rhs,(start,end+18),[0.,0.,0.],method='DOP853',rtol=2e-12,atol=1e-14,
                       max_step=.025)
    P,psi,v=solution.y[:,-1]
    if not solution.success or coupling(psi,math.sqrt(8),20.)[1]!=0:
        raise ArithmeticError('Plane limiting solution failed')
    return dict(scalar=.5*v*v,scalar_velocity=v,K=coupling(psi,math.sqrt(8),20.)[0])


def spectral_validation():
    spec=importlib.util.spec_from_file_location('spherical_bessel',Path(__file__).with_name('bessel.py'))
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.calculate(mesh,run,pulse)


def rounded(obj):
    if isinstance(obj,(float,np.floating)):
        return round(float(obj),12)
    if isinstance(obj,dict):
        return {k:rounded(v) for k,v in obj.items()}
    if isinstance(obj,list):
        return [rounded(v) for v in obj]
    return obj


def calculate():
    radii=(2.,4.,8.,16.,32.,128.,4096.)
    cases={f'R{R:g}':run(radius=R,trace=R in (4,32)) for R in radii}
    cases['short_R8']=run(radius=8.,width=.5)
    cases['broad_R8']=run(radius=8.,width=2.,trace=True)
    cases['one_R32']=run(radius=32.,channel='one')
    cases['one_retuned_R32']=run(radius=32.,channel='one',g=4.)
    cases['odd_R32']=run(radius=32.,channel='odd')
    cases['constant_R32']=run(radius=32.,g=0.)
    cases['empty_R32']=run(radius=32.,channel='none')
    cases['deeper_port_R32']=run(radius=32.,gap=2.,elements=24)
    cases['near_port_R4']=run(radius=4.,gap=.5,elements=6)
    cases['deeper_port_R4']=run(radius=4.,gap=1.5,elements=18)
    cases['bipolar_R4']=run(radius=4.,return_width=40.)
    cases['bipolar_R32']=run(radius=32.,return_width=40.)
    for name,row in cases.items():
        final=row['readout']
        if name.startswith('one') and final['scalar']>.5+1e-8:
            raise ArithmeticError('One-sided parity bound failed')
        if name not in ('odd_R32','empty_R32','constant_R32') and not final['scalar_decoupled']:
            raise ArithmeticError('The asymptotic scalar channel remains coupled')
        if row['max_impulse_residual']>2e-9:
            raise ArithmeticError('Scalar impulse account differs')
    if (cases['R4']['readout']['scalar_velocity']>=0 or cases['broad_R8']['readout']['scalar_velocity']>=0
            or cases['R32']['readout']['scalar_velocity']<=0 or cases['odd_R32']['readout']['scalar']>1e-15):
        raise ArithmeticError('Unexpected scalar direction or parity')
    mesh_checks={}
    for name in ('R2','R4','R8','R32','broad_R8'):
        row=cases[name]
        coarse=run(radius=row['radius'],width=row['width'],elements=6)
        finer=run(radius=row['radius'],width=row['width'],elements=18,tol=5e-11)
        differences=[abs(x['readout']['scalar']-y['readout']['scalar']) for x,y in ((coarse,row),(row,finer))]
        if max(differences)>2e-5:
            raise ArithmeticError((name,differences))
        mesh_checks[name]=dict(coarse_scalar=coarse['readout']['scalar'],archive_scalar=row['readout']['scalar'],
                               finer_scalar=finer['readout']['scalar'],adjacent_differences=differences)
    plane=planar()
    difference=abs(cases['R4096']['readout']['scalar']-plane['scalar'])
    if difference>1e-5:
        raise ArithmeticError('Distant spherical wall failed to approach flat response')
    return dict(units='K0=1, dose=1 per unit proper wall volume, times in initial kinetic length/c',
                cases=rounded(cases),mesh_checks=rounded(mesh_checks),flat_limit=rounded(plane),
                scalar_difference_at_R4096=rounded(difference),helmholtz=spectral_validation())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(),indent=2,allow_nan=False)+'\n')
