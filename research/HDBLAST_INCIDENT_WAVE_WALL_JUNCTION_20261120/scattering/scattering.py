"""Incoming transverse wave packets and normal stress on the fixed central tuned leaf."""
import argparse
from functools import lru_cache
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import quad, quad_vec, simpson, solve_ivp
from scipy.linalg import eigh
from scipy.sparse import coo_matrix, diags

ROOT = Path(__file__).resolve().parents[1]
CENTER = 8.
WIDTH = 2.
SAMPLE_DT = .025


@lru_cache(maxsize=1)
def original():
    path = ROOT / 'prior/box.py'
    spec = importlib.util.spec_from_file_location('archived_central_leaf',path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@lru_cache(maxsize=4)
def geometry(length):
    def rhs(z,y):
        return np.exp(-original().primitive(y)/3)
    sol = solve_ivp(rhs,(0.,length),[0.],method='DOP853',dense_output=True,rtol=8e-14,atol=3e-15)
    if not sol.success:
        raise ArithmeticError(sol.message)
    return sol.sol


def weight(z, law, length):
    y = geometry(length)(z)[0]
    warp = np.exp(-original().primitive(y)/3)
    t = np.tanh(y)
    slope = -(1+t-t**3/3)/3
    if law == 'selected':
        return warp**3, 3*warp*slope
    if law == 'wall':
        return warp / np.cosh(y)**4, warp*(slope-4*t)
    raise ValueError('Unknown gauge weight')


def shapes(z):
    """Rows are B, xB, x²B and their derivatives with respect to z."""
    x = np.atleast_1d(np.asarray(z,dtype=float)-CENTER)/WIDTH
    basis = np.zeros((3,len(x)))
    gradient = np.zeros_like(basis)
    inside = np.abs(x)<1
    u=x[inside]
    scale=1-u*u
    b=np.exp(1.-1./scale)
    bx=-2*u*b/scale**2
    basis[:,inside]=[b,u*b,u*u*b]
    gradient[:,inside]=np.array([bx,b+u*bx,2*u*b+u*u*bx])/WIDTH
    return basis,gradient


def initial_data(z,spec):
    """Return (field, velocity) on arbitrary conformal coordinates, arrival toward z=0."""
    q,v=shapes(z)
    c=np.asarray(spec['coefficients'])*spec['amplitude']
    return c@q,c@v


@lru_cache(maxsize=4)
def prepare_packet(law='selected',packet='balanced'):
    if packet not in ('balanced','raw'):
        raise ValueError('Unknown packet')
    def moments(z):
        b,db=shapes(z)
        w,_=weight(z,law,10.)
        return w*np.stack([b[:,0],db[:,0]])
    constraints = quad_vec(moments,CENTER-WIDTH,CENTER+WIDTH,epsabs=1e-18,epsrel=3e-13)[0]
    if packet=='balanced':
        c=np.r_[np.linalg.solve(constraints[:,:2],-constraints[:,2]),1.]
    else:
        c=np.array([1.,0.,0.])
    def energy(z):
        b,db=shapes(z)
        w,_=weight(z,law,10.)
        return w*((c@db[:,0])**2+.5*(c@b[:,0])**2)
    E=quad_vec(energy,CENTER-WIDTH,CENTER+WIDTH,epsabs=1e-18,epsrel=3e-13)[0]
    amp=1/math.sqrt(E)
    if law=='selected':
        # N=integral G dy. Integrating to infinity avoids a finite-box residue.
        norm=quad(lambda y:math.exp(-2*original().primitive(y)/3),0.,np.inf,epsabs=3e-14)[0]
    else:
        norm=2/3
    cq,cv=amp*(constraints@c)
    return dict(law=law,packet=packet,center=CENTER,half_width=WIDTH,k=1.,
                coefficients=c.tolist(),amplitude=amp,half_norm=norm,
                moment_q=float(cq),moment_v=float(cv),
                pole_energy=float((cq*cq+cv*cv)/(2*norm)))


def phase_data():
    path=Path(__file__).with_name('phase.py')
    spec=importlib.util.spec_from_file_location('static_phase_shift',path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.calculate()


def evolution(law, packet, spacing, end, length, exact=False):
    if end >= length-CENTER-WIDTH:
        raise ValueError('Provide a longer domain so that the far end is outside the pulse cone')
    steps=round(end/SAMPLE_DT)
    cells=round(length/(4*spacing))
    if min(steps,cells)<=0 or abs(steps*SAMPLE_DT-end)>1e-12:
        raise ValueError('Positive duration in multiples of sample interval required')
    times=np.linspace(0.,end,steps+1)
    # Polynomial elements with degree four and Legendre-Gauss-Lobatto nodes.
    h=length/cells
    nodes=np.array([-1.,-math.sqrt(3/7),0.,math.sqrt(3/7),1.])
    quad_weights=np.array([1/10,49/90,32/45,49/90,1/10])*h/2
    mesh=np.arange(cells)[:,None]*h+(nodes+1)*h/2
    z=np.r_[mesh[:,:-1].ravel(),length]
    ids=4*np.arange(cells)[:,None]+np.arange(5)
    bary=np.array([1/math.prod(x-nodes[j] for j in range(5) if i!=j) for i,x in enumerate(nodes)])
    deriv=np.zeros((5,5))
    for i,x in enumerate(nodes):
        for j,y in enumerate(nodes):
            if i!=j:
                deriv[i,j]=bary[j]/(bary[i]*(x-y))
        deriv[i,i]=-deriv[i].sum()
    deriv*=2/h
    w,logs=weight(z,law,length)
    weights=w[ids]*quad_weights
    mass=np.bincount(ids.ravel(),weights.ravel(),minlength=len(z))
    local=np.einsum('ab,ea,ac->ebc',deriv,weights,deriv)
    rows=np.broadcast_to(ids[:,:,None],local.shape)
    cols=np.broadcast_to(ids[:,None,:],local.shape)
    stiff=coo_matrix((local.ravel(),(rows.ravel(),cols.ravel())),shape=(len(z),len(z))).tocsr()
    k=1.
    spec=prepare_packet(law,packet)
    q0,v0=initial_data(z,spec)
    norm=spec['half_norm']
    outputs={key:[] for key in ('field','velocity','pressure','energy','p_radial','bulk_force',
                                'right_boundary_pressure','pole_q','pole_v','near_energy','far_energy',
                                'surface_field_bound')}
    dx=h/4
    if exact:
        scale=diags(1/np.sqrt(mass))
        eigen,basis=eigh((scale@stiff@scale).toarray(),driver='evd')
        freq=np.sqrt(eigen+k*k)
        a=basis.T@(np.sqrt(mass)*q0)
        b=basis.T@(np.sqrt(mass)*v0)
        phase=freq[:,None]*times
        wall=basis[0]/np.sqrt(mass[0])
        field=wall@(a[:,None]*np.cos(phase)+(b/freq)[:,None]*np.sin(phase))
        velocity=wall@(b[:,None]*np.cos(phase)-(a*freq)[:,None]*np.sin(phase))
        return dict(field=field,velocity=velocity,pressure=.5*w[0]*(velocity**2-k*k*field**2))
    minus_omega_squared=-(diags(1/mass)@stiff+diags(np.full(len(z),k*k)))
    def rhs(t,u):
        q,v=u.reshape(2,-1)
        return np.r_[v,minus_omega_squared@q]
    solution=solve_ivp(rhs,(0.,end),np.r_[q0,v0],t_eval=times,method='DOP853',
                       rtol=1e-10,atol=2e-13,max_step=min(.04,dx*1.1))
    if not solution.success:
        raise ArithmeticError(solution.message)
    near=mesh<=1.+1e-12
    far=mesh>=length-1.-1e-12
    for state in solution.y.T:
        q,v=state.reshape(2,-1)
        qfield=q[ids]
        velocity=v[ids]
        grad=qfield@deriv.T
        energy=.5*weights*(velocity**2+grad**2+k*k*qfield**2)
        outputs['field'].append(q[0])
        outputs['velocity'].append(v[0])
        outputs['pressure'].append(.5*w[0]*(v[0]**2-k*k*q[0]**2))
        outputs['right_boundary_pressure'].append(.5*w[-1]*(v[-1]**2-k*k*q[-1]**2))
        outputs['energy'].append(energy.sum())
        outputs['p_radial'].append(-float((weights*velocity*grad).sum()))
        outputs['bulk_force'].append(.5*float((weights*logs[ids]*(velocity**2-grad**2-k*k*qfield**2)).sum()))
        outputs['pole_q'].append(mass@q)
        outputs['pole_v'].append(mass@v)
        outputs['near_energy'].append(energy[near].sum())
        outputs['far_energy'].append(energy[far].sum())
        outputs['surface_field_bound'].append(abs(q[-1])+abs(v[-1]))
    out={name:np.asarray(values) for name,values in outputs.items()}
    energy0=out['energy'][0]
    analytic_q=spec['moment_q']*np.cos(k*times)+spec['moment_v']/k*np.sin(k*times)
    analytic_v=spec['moment_v']*np.cos(k*times)-k*spec['moment_q']*np.sin(k*times)
    pole=(out['pole_v']**2+k*k*out['pole_q']**2)/(2*norm)
    impulse=float(simpson(out['pressure'],x=times))
    bulk=float(simpson(out['bulk_force'],x=times))
    right=float(simpson(out['right_boundary_pressure'],x=times))
    momentum=out['p_radial'][-1]-out['p_radial'][0]
    peak=int(np.argmax(out['pressure']))
    return dict(time=times.tolist(),wall_pressure=out['pressure'].tolist(),
                wall_field=out['field'].tolist(),wall_velocity=out['velocity'].tolist(),
                spec=spec,spacing=dx,sample_dt=SAMPLE_DT,end=end,length=length,
                checks=dict(initial_energy=float(energy0),energy_maximum_drift=float(max(abs(out['energy']/energy0-1))),
                            initial_outward_momentum=float(out['p_radial'][0]),
                            final_outward_momentum=float(out['p_radial'][-1]),
                            integrated_wall_pressure=impulse,integrated_background_force=bulk,
                            integrated_right_pressure=right,momentum_balance_error=float((momentum-impulse-bulk+right)/energy0),
                            max_projection_moment_error=float(max(max(abs(out['pole_q']-analytic_q)),
                                                                  max(abs(out['pole_v']-analytic_v)))),
                            initial_discrete_pole_energy=float(pole[0]),maximum_discrete_pole_energy=float(max(pole)),
                            maximum_near_unit_energy=float(max(out['near_energy'])),
                            peak_pressure=float(out['pressure'][peak]),sampled_peak_time=float(times[peak]),
                            minimum_pressure=float(min(out['pressure'])),maximum_wall_field=float(max(abs(out['field']))),
                            pre_arrival_maximum_wall_field=float(max(abs(out['field'][times<=5.75]))),
                            maximum_energy_in_final_unit=float(max(out['far_energy'])),
                            maximum_far_boundary_field_velocity=float(max(out['surface_field_bound']))))


def evolve(law='selected',packet='balanced',spacing=.02,end=26.,length=None):
    """Time series for a normalized incoming packet. Positive pressure pushes the wall toward -z."""
    if length is None:
        length=max(42.,end+CENTER+WIDTH+4.)
    return evolution(law,packet,spacing,end,length)


def support_bound(law):
    ys = geometry(64.)(np.array([6.,10.,16.,32.,64.]))[0]
    norm = prepare_packet(law,'raw')['half_norm']
    def g(y):
        if law == 'selected':
            return np.exp(-2*original().primitive(y)/3)
        return np.exp(-4*(np.logaddexp(y,-y)-np.log(2.)))
    tails=[quad(g,y,np.inf,epsabs=1e-22,epsrel=4e-13)[0]/norm for y in ys]
    return dict(one_copy_norm=norm,packet_interval_fraction=tails[0]-tails[1],
                exterior_norm_fraction=dict(zip(('6','10','16','32','64'),tails)),
                large_radius_bound_exponent=2. if law=='selected' else 7.2)


def report_precision(x):
    # Timeseries and residuals are quoted in units set by E_in=1 and L0=1.
    # Residual differences below this resolution carry no physical information.
    if isinstance(x,(float,np.floating)):
        return round(float(x),12)
    if isinstance(x,dict):
        return {k:report_precision(v) for k,v in x.items()}
    if isinstance(x,list):
        return [report_precision(v) for v in x]
    return x


def calculate():
    results={}
    for law in ('wall','selected'):
        entries={}
        for packet in ('raw','balanced'):
            meshes=[evolve(law,packet,dx) for dx in (.04,.02,.01)]
            coarse,medium,fine=meshes
            check=evolution(law,packet,.04,26.,42.,exact=True)
            spectral={key:float(max(abs(check[key]-np.asarray(coarse['wall_'+key])))) for key in ('field','velocity','pressure')}
            pressure=[np.array(row['wall_pressure']) for row in meshes]
            scale=max(abs(pressure[-1]))
            shifts=[float(max(abs(f-g))/scale) for f,g in zip(pressure,pressure[1:])]
            stable=float(abs(fine['checks']['integrated_wall_pressure']-medium['checks']['integrated_wall_pressure']))
            if (spectral['pressure']>4e-5 or shifts[1]>.05 or shifts[1]>.4*shifts[0] or stable>.02
                    or abs(fine['checks']['momentum_balance_error'])>.008
                    or fine['checks']['energy_maximum_drift']>1e-6
                    or fine['checks']['max_projection_moment_error']>5e-6):
                raise ArithmeticError((law,packet,shifts,stable,spectral,fine['checks']))
            entries[packet]=dict(**fine,mesh_checks=[r['checks'] for r in meshes[:-1]],
                                 relative_pressure_curve_changes=shifts,medium_to_fine_impulse_change=stable,
                                 exact_semidiscrete_errors=spectral)
        results[law]=entries
    return dict(k=1.,geometry='static central tuned leaf; one copy, L0=1, A(0)=0; E_in=1 per unit tangent volume',
                packet_support=[CENTER-WIDTH,CENTER+WIDTH],branches=report_precision(results),
                support_bounds={law:support_bound(law) for law in results},reflection=phase_data())


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(),indent=2,allow_nan=False)+'\n')
