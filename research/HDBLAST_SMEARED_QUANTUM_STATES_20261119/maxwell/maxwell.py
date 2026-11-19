"""Energy-minimizing Maxwell preparations in a transverse angular sector on the archived metric."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import simpson, solve_ivp
from scipy.linalg import block_diag, eigh
from scipy.optimize import minimize_scalar
from scipy.special import eval_gegenbauer, roots_legendre


ELL = 2
WIDTHS = (.1,.2,.3)
HERE = Path(__file__).resolve().parent


def load_passage(prior):
    path = prior/'prior/passage/passage.py'
    spec = importlib.util.spec_from_file_location('closed_maxwell_passage', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def symmetric(matrix):
    return (matrix+matrix.T)/2


def power(matrix, exponent):
    values, vectors = eigh(symmetric(matrix))
    if min(values) <= 0:
        raise ArithmeticError('Non-positive quadratic form')
    return (vectors*values**exponent) @ vectors.T


class MaxwellModes:
    """Axial vector spherical harmonics, one fixed real (ell,m), even radial parity."""
    def __init__(self, geometry, solution, maximum, nodes, round_control=False):
        self.geometry, self.solution = geometry, solution
        self.round_control = round_control
        self.labels = np.arange(ELL, maximum+1, 2)
        self.size = len(self.labels)
        angle = np.pi*(np.arange(nodes)+.5)/nodes
        sine, cosine = np.sin(angle), np.cos(angle)
        self.sine2 = sine*sine
        self.sine = sine
        self.weight = np.pi/nodes
        self.interpolation = geometry.matrices(cosine)[0]
        degree = self.labels-ELL
        poly = eval_gegenbauer(degree, ELL+1, cosine[:,None])
        diff = np.zeros_like(poly)
        diff[:,1:] = 2*(ELL+1)*eval_gegenbauer(degree[1:]-1,ELL+2,cosine[:,None])
        u = sine[:,None]**(ELL+1)*poly
        du = sine[:,None]**ELL*((ELL+1)*cosine[:,None]*poly-self.sine2[:,None]*diff)
        norm = np.sqrt(self.weight*np.sum(u*u,axis=0))
        self.u,self.du = u/norm,du/norm

    def forms(self,t):
        a = self.geometry.background(t)[0]
        if self.round_control:
            m=d=np.zeros_like(self.sine2)
        else:
            m,d = (self.interpolation @ self.solution.sol(t).reshape(6,-1)[:2].T).T
        radial = m+self.sine2*d
        def gram(field,weight):
            return symmetric(field.T @ ((self.weight*weight)[:,None]*field))
        mass = gram(self.u,np.exp(radial))
        stiffness = gram(self.du,np.exp(-radial))
        stiffness += ELL*(ELL+1)*gram(self.u,np.exp(radial-2*m)/self.sine2)
        volume=4*np.pi*a**3*self.weight*np.sum(self.sine2*np.exp(radial+2*m))
        return mass,stiffness,a,volume


def propagate(ham, end, step):
    size = ham.size
    def rhs(t,state):
        mass,stiffness,a,_=ham.forms(t)
        q,p=state.reshape(2,size,2*size)
        return np.vstack((np.linalg.solve(mass,p),-stiffness@q)).ravel()/a
    sol=solve_ivp(rhs,(0,end),np.eye(2*size).ravel(),method='DOP853',
                  dense_output=True,rtol=2e-12,atol=3e-14,max_step=step)
    if not sol.success:
        raise ArithmeticError(sol.message)
    return lambda t:sol.sol(t).reshape(2*size,2*size)


def preparation(ham, flow, width, points):
    size=ham.size
    y,w=roots_legendre(points)
    # Reflection around the time-symmetric bounce cancels the mixed q,p block.
    t=width*(1+y)/2
    weights=w*np.exp(-1/(1-(t/width)**2))
    weights/=weights.sum()
    qform=np.zeros((size,size))
    pform=np.zeros_like(qform)
    for time,weight in zip(t,weights):
        mass,k,a,volume=ham.forms(time)
        T=flow(time)
        q,p=T[:size],T[size:]
        form=symmetric(q.T@k@q+p.T@np.linalg.solve(mass,p))*(weight/(a*volume))
        qform+=form[:size,:size]
        pform+=form[size:,size:]
    half=power(pform,.5)
    invhalf=power(pform,-.5)
    omega=symmetric(invhalf@power(half@qform@half,.5)@invhalf)
    residual=np.linalg.norm(omega@pform@omega-qform)/np.linalg.norm(qform)
    launch=np.vstack((power(omega,-.5),-1j*power(omega,.5)))/np.sqrt(2.)
    return launch, residual


def readout(ham, flow, time, launch):
    mass,k,a,_=ham.forms(time)
    w2,v=eigh(k,mass)
    q,p=np.split(flow(time)@launch,2)
    omega=np.sqrt(w2)
    # Project onto the instantaneous physical-time Hamiltonian at the same endpoint.
    qmode=v.T@mass@q
    pmode=v.T@p
    beta=np.sqrt(omega[:,None]/2)*qmode-1j*pmode/np.sqrt(2*omega[:,None])
    occupation=np.sum(abs(beta)**2,axis=1)
    # A one-quantum addition in the first specified launch column has the
    # positive stress difference of one normalized complex classical solution.
    q0,p0=q[:,0],p[:,0]
    added_energy=float((p0.conj()@np.linalg.solve(mass,p0)+q0.conj()@k@q0).real/a)
    return dict(added_first_launch_mode_energy=added_energy, low_occupations={str(n):float(occupation[i]) for i,n in enumerate(ham.labels) if n<=8},
                excitation_energy=float(omega@occupation/a),
                sector_occupation=float(occupation.sum()))


def run(geometry,solution,events,maximum=24,nodes=160,samples=48,step=.005,round_control=False):
    ham=MaxwellModes(geometry,solution,maximum,nodes,round_control)
    flow=propagate(ham,events['after'],step)
    checks=[]
    results={}
    for width in WIDTHS:
        launch,residual=preparation(ham,flow,width,samples)
        checks.append(residual)
        results[str(width)]={name:readout(ham,flow,t,launch) for name,t in events.items()}
    mass,k,_,_=ham.forms(0)
    form=block_diag(k,np.linalg.inv(mass))
    T=flow(events['after'])
    size=ham.size
    j=np.block([[np.zeros((size,size)),np.eye(size)],[-np.eye(size),np.zeros((size,size))]])
    result=dict(maximum=maximum, states=results,
                maximum_initial_riccati_relative_residual=float(max(checks)),
                symplectic_max_abs_error=float(np.max(abs(T.T@j@T-j))),
                launch_hamiltonian_minimum_eigenvalue=float(eigh(form,eigvals_only=True)[0]))
    return result


def distances(left,right):
    low=energy=high=0.
    for width in left['states']:
        for label,record in left['states'][width].items():
            other=right['states'][width][label]
            low=max(low,max(abs(v-other['low_occupations'][n]) for n,v in record['low_occupations'].items()))
            energy=max(energy,abs(record['excitation_energy']-other['excitation_energy']))
            high=max(high,abs(record['sector_occupation']-other['sector_occupation']))
    return dict(max_low_occupation_shift=float(low),max_sector_occupation_shift=float(high),
                max_sector_excitation_energy_shift=float(energy))


def redshift_bounds(geometry, solution, start, stop, angular_points, time_points):
    """Numerical extrema in the positivity bound for total electromagnetic energy."""
    angle=np.linspace(0.,np.pi/2,angular_points)
    # Only values enter this bound. The frozen spacetime basis is even in cos(chi).
    basis=geometry.matrices(np.cos(angle))[0]
    grid=np.linspace(start,stop,time_points)
    rates=[]
    for t in grid:
        H=geometry.background(t)[1]
        b,e=solution.sol(t).reshape(6,-1)[[2,3]]
        def local(chi, u):
            ht=H+u@b
            w=np.sin(chi)**2*(u@e)
            return np.array([-ht+w,-ht-w])
        samples=local(angle,basis)
        def maximum(sign):
            extremes=[]
            for j in (0,1):
                vals=sign*samples[j]
                extremes.append(float(np.max(vals)))
                peaks=np.flatnonzero((vals[1:-1]>vals[:-2])&(vals[1:-1]>=vals[2:]))+1
                for i in peaks:
                    def objective(theta):
                        u=geometry.matrices(np.array([np.cos(theta)]))[0][0]
                        return -sign*local(theta,u)[j]
                    best=minimize_scalar(objective,bounds=(angle[i-1],angle[i+1]),method='bounded')
                    extremes.append(-best.fun)
            return max(extremes)
        rates.append([-maximum(-1),maximum(1)])
    return np.exp(simpson(np.array(rates),x=grid,axis=0)).tolist()


def parent_energy(geometry,solution,time):
    record=geometry.record(time,solution.sol(time))
    # The archive reports kinetic + gradient + (potential minus hilltop), in Vstar units.
    mean=1+record['energy_over_Vstar']['signed_sum']
    return 3*geometry.branch.u_top*mean*record['volume']


def bounce_containment(geometry,solution,crossing,angular_points,time_points):
    bounds=redshift_bounds(geometry,solution,0.,crossing,angular_points,time_points)
    initial=parent_energy(geometry,solution,0.)
    final=parent_energy(geometry,solution,crossing)
    if min(initial,final)<=0:
        raise ArithmeticError('Integrated parent energy must be positive at both comparison slices')
    parent_ratio=final/initial
    return dict(classical_maxwell_energy_ratio=bounds,
                parent_total_energy_bounce=initial,parent_total_energy_crossing=final,
                parent_total_energy_ratio=parent_ratio,
                photon_to_parent_ratio_change=[x/parent_ratio for x in bounds],
                homogeneous_round_a_bounce_over_a_crossing=(
                    geometry.background(0.)[0]/geometry.background(crossing)[0]))


def tidal_invariant(geometry,solution,t):
    """Cross-check the spherical electric Weyl eigenvalue gap in two decompositions."""
    state=solution.sol(t)
    fields=geometry.fields(state,geometry.check)
    y,w=geometry.check_y,geometry.check_w
    _,H,inv,_,_,dk,dj,*_=geometry.geometry(t,fields,y,w)
    (m,d,b,e,_,_),(_,_,_,qy),_=fields
    ht=H+b
    hr=ht+w*e
    gradient2=inv*w*qy*qy
    space=inv*(dk-dj)+ht*w*e-gradient2/2
    e_dot=geometry.check[0]@geometry.equation(t,state).reshape(6,-1)[3]
    time=-w*e_dot-w*e*(hr+ht)+gradient2/2
    weight=geometry.check_weight*np.exp(3*m+w*d)
    difference=float(np.max(abs(time-space)))/geometry.branch.h**2
    # Spherical symmetry sets the magnetic Weyl tensor to zero.
    return dict(weyl_squared_spatial_average_over_h4=float(
                    16/3*(weight@(space**2))/weight.sum()/geometry.branch.h**4),
                maximum_3plus1_vs_time_gap_difference_over_h2=difference)


def calculate(prior):
    passage=load_passage(prior)
    gravity=prior/'prior/gravity'
    solver=passage.dynamic.intrinsic.gravity_module(gravity)
    branch=solver.constraint.load_branches(gravity/'previous')['phi_m=-0.5']
    runs={}
    for amp in (.001,-.001):
        geometry,solution,info=passage.geometry_through_passage(solver,branch,amp,64,192)
        times={'bounce':0.,**info['sample_times']}
        base=run(geometry,solution,times)
        # Increasing n changes the measured sector total. Compare low channels separately.
        wide=run(geometry,solution,times,maximum=36)
        coarse=run(geometry,solution,times,maximum=12)
        fine_geometry,fine_solution,fine_info=passage.geometry_through_passage(solver,branch,amp,96,256,
                                                                           tolerance=5e-13,max_step=.002)
        fine_times={'bounce':0.,**fine_info['sample_times']}
        fine=run(fine_geometry,fine_solution,fine_times,nodes=224,samples=80,step=.0025)
        uniform=run(geometry,solution,times,round_control=True)
        round_null=max(r['sector_occupation'] for entry in uniform['states'].values() for r in entry.values())
        resolution=distances(base,fine)
        radial=distances(base,wide)
        if base['maximum_initial_riccati_relative_residual']>2e-11 or base['symplectic_max_abs_error']>2e-10:
            raise ArithmeticError('Variational or canonical control failed')
        if round_null>1e-21 or radial['max_low_occupation_shift']>5e-7 or resolution['max_low_occupation_shift']>5e-11:
            raise ArithmeticError(f'Maxwell convergence failure: {round_null}, {radial}, {resolution}')
        budget = redshift_bounds(fine_geometry,fine_solution,fine_times['before'],fine_times['after'],257,49)
        budget_coarse = redshift_bounds(geometry,solution,times['before'],times['after'],129,25)
        containment=bounce_containment(fine_geometry,fine_solution,fine_times['crossing'],257,801)
        coarse_containment=bounce_containment(geometry,solution,times['crossing'],129,401)
        containment_change=max(abs(x-y) for key in containment for x,y in
                               zip(np.atleast_1d(containment[key]),np.atleast_1d(coarse_containment[key])))
        changes = [abs(x-y) for x,y in zip(budget,budget_coarse)]
        state = base['states']['0.2']
        positive_ratio = state['after']['added_first_launch_mode_energy']/state['before']['added_first_launch_mode_energy']
        positive_bounce_ratio=state['crossing']['added_first_launch_mode_energy']/state['bounce']['added_first_launch_mode_energy']
        lower,upper=containment['classical_maxwell_energy_ratio']
        if containment_change>1e-6 or not lower<positive_bounce_ratio<upper:
            raise ArithmeticError(f'Bounce-to-crossing budget check failed: {containment_change}')
        containment['maximum_refinement_shift']=float(containment_change)
        containment['window_0p2_added_mode_energy_ratio']=positive_bounce_ratio
        if max(changes)>1e-6 or not budget[0]<positive_ratio<budget[1]:
            raise ArithmeticError('Maxwell redshift inequality or convergence failed')
        mass,stiff,_,_=MaxwellModes(geometry,solution,24,160,True).forms(0.)
        uniform_spectrum_error=float(np.max(abs(eigh(stiff,mass,eigvals_only=True)-np.arange(3,26,2)**2)))
        if uniform_spectrum_error>1e-9:
            raise ArithmeticError('Vacuum cylinder spectral control failed')
        tidal={label:tidal_invariant(fine_geometry,fine_solution,time) for label,time in fine_times.items()}
        if max(v['maximum_3plus1_vs_time_gap_difference_over_h2'] for v in tidal.values())>2e-8:
            raise ArithmeticError(f'Tidal decomposition check failed: {tidal}')
        runs[f'{amp:g}']=dict(times=times,primary=base,weyl=tidal,bounce_to_first_crossing=containment,
            shape_round_cylinder_frequency_squared_error=uniform_spectrum_error,
            maxwell_stress_energy_ratio_bounds=budget, bounds_refinement=max(changes),
            window_0p2_added_mode_energy_ratio=positive_ratio,expanded=wide,restricted=coarse,
            equal_band_refinement=resolution,radial_24_vs_36=radial,radial_12_vs_24=distances(base,coarse),
            round_control_max_sector_occupation=round_null)
    return dict(ell=ELL,one_real_m=True,polarization='axial',radial_parity='even',
                normalization='one real (ell,m) block; physical-time energy in inverse model time units',
                weights='time-averaged proper-volume mean electromagnetic energy density',
                windows=list(WIDTHS),runs=runs)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior',type=Path,default=HERE.parent/'prior')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.prior.resolve()),indent=2,allow_nan=False)+'\n')
