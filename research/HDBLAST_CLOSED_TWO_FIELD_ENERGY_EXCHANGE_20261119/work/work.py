"""State-difference work on the two archived first-passage backgrounds."""
import argparse
from copy import copy
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre


INITIAL_MODES = (2, 4, 8)
POWERS = ("metric", "coupling")
HERE = Path(__file__).resolve().parent


def passage_module(prior):
    spec = importlib.util.spec_from_file_location("archived_passage", prior/"passage/passage.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def evaluate_energy(ham, t, field, momentum):
    mass, stiffness = ham.matrices(t)
    velocity = np.linalg.solve(mass, momentum)
    return np.real(np.sum(momentum.conj()*velocity + field.conj()*(stiffness @ field), axis=0))


class MixedTimes:
    """Freeze the matter coefficient while varying the metric, or conversely."""
    def __init__(self, geometry_solution, matter_time):
        self.solution, self.matter_time = geometry_solution, matter_time

    def sol(self, metric_time):
        a = self.solution.sol(metric_time).reshape(6, -1)
        b = self.solution.sol(self.matter_time).reshape(6, -1)
        return np.concatenate((a[:4], b[4:])).ravel()


def partial_difference(ham, t, field, momentum, dx, vary_metric):
    varied = copy(ham)
    results = []
    for h in (-dx, dx):
        tg, tz = (t+h, t) if vary_metric else (t, t+h)
        varied.solution = MixedTimes(ham.solution, tz)
        results.append(evaluate_energy(varied, tg, field, momentum))
    return (results[1]-results[0])/(2*dx)


def power_and_energy(ham, t, field, momentum):
    geom = ham.geometry
    local = ham.interpolation @ ham.solution.sol(t).reshape(6,geom.n).T
    m,d,b,e,q,pi = local.T
    a, hubble, _ = geom.background(t)
    w = ham.sine2
    radial_metric = m+w*d
    tangential_expansion = hubble+b
    radial_expansion = tangential_expansion+w*e
    if ham.kinetic == "matter":
        theta = geom.branch.top+q
        f = geom.branch.fitted.spline(theta)[:,1]
        z = 1/f
        log_rate = -geom.branch.fitted.spline(theta, 1)[:,1]*pi/f
    else:
        z = np.ones_like(q)
        log_rate = np.zeros_like(q)
    weight = ham.weights*z
    mass_density = weight*a**3*np.exp(radial_metric+2*m)
    angular_density = weight*a*6*np.exp(radial_metric)/w
    radial_density = weight*a*np.exp(2*m-radial_metric)
    mass, _ = ham.matrices(t)
    velocity = np.linalg.solve(mass, momentum)
    kinetic = mass_density[:,None]*abs(ham.value @ velocity)**2
    angular = angular_density[:,None]*abs(ham.value @ field)**2
    radial = radial_density[:,None]*abs(ham.derivative @ field)**2
    metric = (-kinetic*(radial_expansion+2*tangential_expansion)[:,None]
              + radial*(2*tangential_expansion-radial_expansion)[:,None]
              + angular*radial_expansion[:,None])
    coupling = log_rate[:,None]*(angular+radial-kinetic)
    return (kinetic+angular+radial).sum(axis=0), np.stack((metric.sum(axis=0),coupling.sum(axis=0)))


def trajectory(passage, geom, solution, events, kinetic, cutoff=36, angles=192, step=.005,
               tolerance=1e-11, time_nodes=40, independent_checks=False):
    spec = passage.dynamic
    ham = spec.FieldHamiltonian(geom,solution,cutoff,angles,kinetic)
    initial_field,initial_momentum = spec.initial_modes(ham)
    columns = [ham.labels.tolist().index(n) for n in INITIAL_MODES]
    nrad, states = len(ham.labels),len(INITIAL_MODES)
    def unpack(z):
        return z.reshape(2,nrad,states)
    def rhs(t,z):
        field,momentum = unpack(z)
        mass,stiffness = ham.matrices(t)
        return np.stack((np.linalg.solve(mass,momentum),-stiffness@field)).ravel()

    points = dict(bounce=0.,**events)
    z0 = np.stack((initial_field[:,columns],initial_momentum[:,columns])).ravel()
    modes = solve_ivp(rhs,(0.,points['after']),z0,method='DOP853',dense_output=True,
                      rtol=tolerance,atol=2e-14,max_step=step)
    if not modes.success:
        raise ArithmeticError(modes.message)
    samples = {}
    gram_error = matrix_error = split_error = global_error = 0.
    for name,t in points.items():
        field,momentum = unpack(modes.sol(t))
        energy,power = power_and_energy(ham,t,field,momentum)
        gram = 1j*(field.conj().T@momentum - momentum.conj().T@field)
        gram_error = max(gram_error,float(np.max(abs(gram-np.eye(states)))))
        original = evaluate_energy(ham,t,field,momentum)
        omega,alpha,beta = spec.project(ham,t,field,momentum)
        endpoint = omega@(abs(alpha)**2+abs(beta)**2)
        matrix_error = max(matrix_error,float(np.max(abs(original-energy))),float(np.max(abs(endpoint-energy))))
        samples[name] = (energy,power)
        if independent_checks and name=='crossing':
            # Vary the explicit metric and matter times independently at fixed canonical phase point.
            for j in range(len(POWERS)):
                rough = partial_difference(ham,t,field,momentum,2e-4,j==0)
                fine = partial_difference(ham,t,field,momentum,1e-4,j==0)
                split_error = max(split_error,float(np.max(abs((4*fine-rough)/3-power[j]))))

    intervals = dict(to_before=('bounce','before'),first_half=('before','crossing'),
                     second_half=('crossing','after'),bracket=('before','after'),
                     whole=('bounce','after'))
    def integral(start,stop,nodes):
        y,w = roots_legendre(nodes)
        times = (start+stop)/2 + y*(stop-start)/2
        values = [power_and_energy(ham,t,*unpack(modes.sol(t)))[1] for t in times]
        return (stop-start)/2*np.einsum('i,ijk->jk',w,values)

    mode_rows = {str(n):dict(initial_energy=float(samples["bounce"][0][i]),samples={},intervals={})
                 for i,n in enumerate(INITIAL_MODES)}
    for name,(energy,power) in samples.items():
        for i,n in enumerate(INITIAL_MODES):
            mode_rows[str(n)]['samples'][name] = dict(x=float(points[name]),energy=float(energy[i]),
                metric_power=float(power[0,i]),coupling_power=float(power[1,i]))
    for key,(first,last) in intervals.items():
        work = integral(points[first],points[last],time_nodes)
        if independent_checks:
            refined = integral(points[first],points[last],time_nodes+20)
            global_error=max(global_error,float(np.max(abs(refined-work))))
        for i,n in enumerate(INITIAL_MODES):
            delta=float(samples[last][0][i]-samples[first][0][i])
            residual=float(delta-work[:,i].sum())
            mode_rows[str(n)]['intervals'][key] = dict(
                energy_change=delta,metric_work=float(work[0,i]),coupling_work=float(work[1,i]),
                balance_residual=residual)
    return dict(modes=mode_rows,checks=dict(canonical_error=gram_error,
        direct_and_spectral_energy_error=matrix_error,
        **(dict(independent_partial_power_error=split_error,
                temporal_quadrature_difference=global_error) if independent_checks else {})))


def vacuum_balance(passage, geom, solution, events, kinetic, cutoff, *, fine=False):
    """Excitation energy above the instantaneous ground state of the whole finite band."""
    spec = passage.dynamic
    ham = spec.FieldHamiltonian(geom,solution,cutoff,272 if fine else 192,kinetic)
    field0, momentum0 = spec.initial_modes(ham)
    n = len(ham.labels)
    def rhs(t,z):
        field, momentum = z.reshape(2,n,n)
        mass, stiffness = ham.matrices(t)
        return np.stack((np.linalg.solve(mass,momentum),-stiffness @ field)).ravel()
    trajectory = solve_ivp(rhs,(0.,events["after"]),np.stack((field0,momentum0)).ravel(),
                           method="DOP853",rtol=2e-12 if fine else 2e-11,atol=2e-14,
                           max_step=.0025 if fine else .005,dense_output=True)
    if not trajectory.success:
        raise ArithmeticError(trajectory.message)
    def evaluate(t):
        field,momentum=trajectory.sol(t).reshape(2,n,n)
        mass,omega,vectors=ham.modes(t)
        ref_field=vectors/np.sqrt(2*omega)
        ref_momentum=-1j*(mass@vectors)*np.sqrt(omega/2)
        energy,power=power_and_energy(ham,t,field,momentum)
        ref_energy,ref_power=power_and_energy(ham,t,ref_field,ref_momentum)
        # The Feynman-Hellmann derivative of the finite ground-state energy
        # is its expectation of the same partial Hamiltonian derivative.
        return np.sum(energy-ref_energy)/2, np.sum(power-ref_power,axis=1)/2
    bounds={"bounce":0.,**events}
    checks=[]
    energies={}
    for key,t in bounds.items():
        energy,_=evaluate(t)
        u,p=trajectory.sol(t).reshape(2,n,n)
        omega,_,beta=spec.project(ham,t,u,p)
        modal_energy=np.sum(omega[:,None]*abs(beta)**2)
        checks.append(abs(energy-modal_energy))
        energies[key]=float(modal_energy)
    rows={}
    pairs={"to_before":("bounce","before"),"bracket":("before","after"),"whole":("bounce","after")}
    quad_nodes,weights=roots_legendre(80 if fine else 56)
    for lab,(first,last) in pairs.items():
        a,b=bounds[first],bounds[last]
        powers=[evaluate((a+b)/2+(b-a)*y/2)[1] for y in quad_nodes]
        wg,wz=(b-a)/2*(weights@powers)
        delta=energies[last]-energies[first]
        checks.append(abs(delta-wg-wz))
        rows[lab]=dict(energy_change=delta,metric_work=float(wg),coupling_work=float(wz),
                       balance_residual=float(delta-wg-wz))
    if max(checks)>2e-8:
        raise ArithmeticError(f"Finite-vacuum balance or modal energy mismatch: {max(checks)}")
    return dict(maximum_n=cutoff,sample_excitation_energies=energies,intervals=rows,
                max_balance_or_modal_energy_error=float(max(checks)))


def shifted(a,b):
    maximum = 0.
    for label,mode in a['modes'].items():
        for key,interval in mode['intervals'].items():
            for term in ('energy_change','metric_work','coupling_work'):
                maximum=max(maximum,abs(interval[term]-b['modes'][label]['intervals'][key][term]))
    return maximum


def closed_de_sitter_sum(passage, h, t=.47):
    """Exact homogeneous control for the UV behavior of instantaneous subtraction."""
    a = np.cosh(h*t)/h
    s = np.sinh(h*t)
    eta = np.arctan(s)
    def occupations(labels):
        k = labels+1
        sin_term = np.sin(k*eta)/k
        return ((sin_term-s*np.cos(k*eta))**2/(k*k-1)+(s*sin_term)**2)/4

    # Independent solution with arbitrary normalized initial Cauchy data.
    labels = np.arange(2,50,2)
    omega0 = h*np.sqrt(labels*(labels+2))
    mass0 = h**-3
    q0 = np.diag(1/np.sqrt(2*mass0*omega0))
    p0 = np.diag(-1j*np.sqrt(mass0*omega0/2))
    field,momentum = passage.dynamic.de_sitter_modes(h,labels,t,q0,p0)
    mw = a*a*np.sqrt(labels*(labels+2))
    beta = np.sqrt(mw/2)*field.diagonal()-1j*momentum.diagonal()/np.sqrt(2*mw)
    mode_error=float(np.max(abs(abs(beta)**2-occupations(labels))))

    coefficient = s*s/(8*a)
    cutoffs = [512,2048,8192,32768,131072,524288]
    rows=[]
    for cutoff in cutoffs:
        n = np.arange(2,cutoff+1,2, dtype=float)
        energy = float(np.sqrt(n*(n+2))@occupations(n)/a)
        row = dict(maximum_n=cutoff,excitation_energy=energy)
        if rows:
            last=rows[-1]
            observed=(energy-last["excitation_energy"])/np.log(cutoff/last["maximum_n"])
            row["logarithmic_slope_over_prediction"]=float(observed/coefficient)
        rows.append(row)
    if mode_error>1e-12 or abs(rows[-1]["logarithmic_slope_over_prediction"]-1)>5e-5:
        raise ArithmeticError("Independent homogeneous logarithmic check failed")
    return dict(x=t,h=h,a=float(a),aH=float(s),log_coefficient=float(coefficient),
                exact_mode_check=mode_error,partial_sums=rows)


def launch_state_check(h):
    """Bounce energy difference between two explicit nonzero-mode states on the round sphere."""
    labels = np.arange(1,49)
    ell = labels+1.
    lam = ell*ell-1
    freq = h*np.sqrt(lam)
    mass0 = h**-3
    normalization = h/np.sqrt(2*ell*lam)
    q = normalization*ell
    momentum = -1j*normalization*lam/h**2
    comm = 1j*(q.conj()*momentum - momentum.conj()*q)
    beta = np.sqrt(mass0*freq/2)*q-1j*momentum/np.sqrt(2*mass0*freq)
    exact = 1/(4*ell*np.sqrt(lam)*(ell+np.sqrt(lam))**2)
    number_error = float(np.max(abs(abs(beta)**2-exact)))
    bare_energy = (abs(momentum)**2/mass0+mass0*freq**2*abs(q)**2)/2
    energy_error = float(np.max(abs(bare_energy-freq/2-freq*exact)))
    comm_error = float(np.max(abs(comm-1)))
    # Sum the manifestly positive, cancellation-free expression including all L² modes.
    def shell(n):
        k = n+1
        return h*k/(4*(k+np.sqrt(k*k-1))**2)
    rows = []
    for cutoff in (100,1000,10000,100000,1000000):
        total = np.sum(shell(np.arange(1,cutoff+1,dtype=float)))
        doubled = np.sum(shell(np.arange(cutoff+1,2*cutoff+1,dtype=float)))
        rows.append(dict(maximum_n=cutoff,reference_minus_instantaneous_energy=float(total),
                         doubling_log_slope_over_prediction=float(doubled/(np.log(2)*h/16))))
    if (number_error>1e-14 or energy_error>1e-12 or comm_error>1e-13
            or abs(rows[-1]['doubling_log_slope_over_prediction']-1)>3e-6):
        raise ArithmeticError("Exact normalized launch modes or full-sphere logarithmic coefficient disagreed")
    return dict(h=h,excluded_mode="n=0",multiplicity="(n+1)^2 real harmonics",
                total_energy_log_coefficient=h/16,
                physical_density_log_coefficient=h**4/(32*np.pi**2),
                reference_commutator_error=comm_error,bogoliubov_number_error=number_error,
                direct_energy_error=energy_error,partial_sums=rows)


def calculate(prior):
    passage = passage_module(prior)
    module = passage.dynamic.intrinsic.gravity_module(prior/'gravity')
    key,branch = next((k,v) for k,v in module.constraint.load_branches(prior/'gravity/previous').items()
                      if v.fitted.pm<0)
    data = dict(branch=key,units="one extra launch quantum, hbar=1, inverse model-time energy",
                state_difference="one initial normal-mode quantum minus the same launch vacuum",
                angular_sector=dict(ell=2,one_m=True,parity="even n",primary_cutoff=36,
                                    initial_labels=list(INITIAL_MODES)),
                bracket_halfwidth=passage.HALF_WIDTH,runs={})
    data['homogeneous_vacuum_control']=closed_de_sitter_sum(passage,branch.h)
    data['launch_state_control']=launch_state_check(branch.h)
    worst_balance=worst_control=0.
    for amp in (.001,-.001):
        geom,sol,entry=passage.geometry_through_passage(module,branch,amp,64,160)
        high,highsol,highentry=passage.geometry_through_passage(
            module,branch,amp,96,256,tolerance=2e-13,max_step=.0015)
        events=entry['sample_times']
        case=dict(first_location="poles" if amp>0 else "equator",times=events,kinetic={})
        data['runs'][f'{amp:g}']=case
        for kinetic in ('matter','geometric'):
            primary=trajectory(passage,geom,sol,events,kinetic,independent_checks=True)
            case['kinetic'][kinetic]=primary
            changes={}
            for name, geom_i,sol_i,t,options in [
                    ('n24',geom,sol,events,dict(cutoff=24)),
                    ('n48',geom,sol,events,dict(cutoff=48)),
                    ('angular272',geom,sol,events,dict(angles=272)),
                    ('field_time',geom,sol,events,dict(step=.0025,tolerance=1e-12)),
                    ('background',high,highsol,highentry['sample_times'],{})]:
                changes[name]=shifted(primary,trajectory(passage,geom_i,sol_i,t,kinetic,**options))
            primary['max_absolute_integrated_changes']=changes
            worst_control=max(worst_control,max(changes.values()))
            c=primary['checks']
            if (c['canonical_error']>1e-9 or c['direct_and_spectral_energy_error']>5e-11 or
                c['independent_partial_power_error']>1e-8 or c['temporal_quadrature_difference']>3e-8):
                raise ArithmeticError(f'Power, quadrature or canonical control failed: {c}')
            for n,row in primary['modes'].items():
                for interval in row['intervals'].values():
                    worst_balance=max(worst_balance,abs(interval['balance_residual']))
                    if kinetic=='geometric' and interval['coupling_work']!=0.:
                        raise ArithmeticError('Constant matter coefficient has nonzero coupling work')
            vacuum = {}
            for cutoff in (12,24,36,48):
                vacuum[str(cutoff)] = vacuum_balance(passage,geom,sol,events,kinetic,cutoff)
            primary["instantaneous_subtracted_vacuum"]=vacuum
            check=vacuum_balance(passage,high,highsol,highentry['sample_times'],kinetic,36,fine=True)
            worst=max(abs(vacuum['36']['intervals'][k][v]-check['intervals'][k][v])
                      for k in check['intervals'] for v in ('energy_change','metric_work','coupling_work'))
            if worst>1e-7:
                raise ArithmeticError(f'Finite-vacuum spatial and temporal refinement failed: {worst}')
            primary["fixed_n36_vacuum_max_refinement"]=float(worst)
        print('Integrated state-difference work',amp,flush=True)
    data['checks']=dict(max_absolute_balance_residual=worst_balance,
                        max_absolute_work_or_energy_refinement=worst_control)
    if worst_balance>5e-8 or worst_control>5e-7:
        raise ArithmeticError(f'Energy or refinement check failed: {data["checks"]}')
    return data


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior',type=Path,default=HERE.parent/'prior')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.write_text(json.dumps(calculate(args.prior.resolve()),indent=2,allow_nan=False)+'\n')
