"""Time-smeared energy-minimizing finite-sector states on two Einstein--scalar spacetimes."""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import block_diag, eigh
from scipy.special import roots_legendre


HERE = Path(__file__).resolve().parent
WIDTHS = (.1, .2, .3)


def sym(matrix):
    return (matrix + matrix.T)/2


def powers(matrix, *exponents):
    values, vectors = eigh(sym(matrix))
    if values.min() <= 0:
        raise ArithmeticError("Quadratic form lost positivity")
    return [(vectors*values**p) @ vectors.T for p in exponents]


def canonical_ground(form):
    """Full Williamson covariance in original canonical variables, with hbar=1."""
    count = form.shape[0]//2
    identity = np.eye(count)
    zeros = np.zeros_like(identity)
    J = np.block([[zeros, identity], [-identity, zeros]])
    root, inverse = powers(form, .5, -.5)
    antisymmetric = root @ J @ root
    absolute, = powers(antisymmetric.T @ antisymmetric, .5)
    return sym(inverse @ absolute @ inverse/2)


def energy_ground(a, b):
    """Exact solution for a time-reflection invariant sampled quadratic form diag(a,b)."""
    br, bi = powers(b, .5, -.5)
    frequencies, inverse = powers(br @ a @ br, .5, -.5)
    return block_diag(sym(br @ inverse @ br/2), sym(bi @ frequencies @ bi/2))


def expectation(form, covariance):
    return float(np.einsum('ij,ji->', form, covariance)/2)


def integrate_interval(function, start, stop, order):
    x,w = roots_legendre(order)
    t = (stop+start+(stop-start)*x)/2
    return (stop-start)/2 * np.tensordot(w, np.stack([function(x) for x in t]), axes=(0,0))


class Evolution:
    def __init__(self, ham, times, order, step, tolerance):
        self.ham = ham
        self.count = len(ham.labels)
        mass, stiff, *_ = ham.matrices(0.)
        # Fixed symplectic rescaling tames the condition number of the fundamental solution.
        qscale, pscale = powers(energy_ground(stiff,np.linalg.inv(mass))[:self.count,:self.count]*2, .5, -.5)
        self.scaling = block_diag(qscale, pscale)
        self.J = np.block([[np.zeros_like(qscale), np.eye(self.count)],
                           [-np.eye(self.count), np.zeros_like(qscale)]])
        dim = 2*self.count
        def rhs(t, state):
            q, p = np.split(state.reshape(dim,dim), 2)
            m,k,*_ = ham.matrices(t)
            return np.vstack((np.linalg.solve(m,p), -k@q)).ravel()
        # Retain the transfer only where sampled. Memory scales with the number
        # of quadrature nodes instead of the number of internal integration steps.
        gauss = (roots_legendre(order)[0]+1)/2
        start, end = times['before'], times['after']
        bracket = (end+start+(end-start)*roots_legendre(order//2)[0])/2
        query = np.unique(np.r_[list(times.values()),bracket,*[w*gauss for w in WIDTHS]])
        flow = solve_ivp(rhs,(0.,end),self.scaling.ravel(),method='DOP853',t_eval=query,
                         atol=tolerance/300,rtol=tolerance,max_step=step)
        if not flow.success:
            raise ArithmeticError(flow.message)
        self.transfers = {t:state.reshape(dim,dim) for t,state in zip(flow.t,flow.y.T)}

    def state(self,t):
        return self.transfers[t]

    def volume(self,t):
        ham = self.ham
        m,d = (ham.interpolation @ ham.solution.sol(t).reshape(6,-1)[:2].T).T
        a=ham.geometry.background(t)[0]
        return float(2*np.pi**2*a**3*(ham.weights@np.exp(3*m+ham.w*d)))

    def forms(self, t, energy_only=False):
        q,p = np.split(self.state(t), 2)
        mass,stiff,mz,kz,mg,kg = self.ham.matrices(t)
        v = np.linalg.solve(mass,p)
        E = sym(q.T@stiff@q + p.T@v)
        if energy_only:
            return E
        return np.stack((E, sym(q.T@kz@q-v.T@mz@v), sym(q.T@kg@q-v.T@mg@v)))

    def sample_state(self,width,order):
        x,w = roots_legendre(order)
        y = (x+1)/2
        prob = w*np.exp(-1/(1-y*y))
        prob /= sum(prob)
        # Time reflection of the full symmetric support cancels the mixed block exactly.
        matrix = sum(weight*self.forms(t,True)/self.volume(t) for t,weight in zip(width*y,prob))
        count=self.count
        a,b = sym(matrix[:count,:count]), sym(matrix[count:,count:])
        sampled = block_diag(a,b)
        cov = energy_ground(a,b)
        full = canonical_ground(sampled)
        physical = self.scaling @ cov @ self.scaling.T
        cqq, cpp = physical[:count,:count], physical[count:,count:]
        energy = expectation(sampled,cov)
        # Compare to the best product of independent round (q_n,p_n) oscillator states.
        undo = np.linalg.inv(self.scaling)
        raw = undo.T @ sampled @ undo
        product = float(np.sqrt(np.diag(raw)[:count]*np.diag(raw)[count:]).sum()/2)
        qq_corr = cqq/np.sqrt(np.outer(np.diag(cqq),np.diag(cqq)))
        pp_corr = cpp/np.sqrt(np.outer(np.diag(cpp),np.diag(cpp)))
        stats=dict(own_sampled_mean_density=energy, independent_round_oscillator_mean_density=product,
                   correlated_state_improvement=float(product-energy),
                   maximum_offdiagonal_field_correlation=float(np.max(abs(qq_corr-np.eye(count)))),
                   maximum_offdiagonal_momentum_correlation=float(np.max(abs(pp_corr-np.eye(count)))),
                   direct_Williamson_covariance_error=float(np.max(abs(cov-full))),
                   purity_error=float(np.max(abs(cov@self.J@cov-self.J/4))),
                   low_variances={str(n):dict(field=float(cqq[i,i]),momentum=float(cpp[i,i]))
                                 for i,n in enumerate(self.ham.labels) if n in (2,4,8)})
        if product < energy-1e-10 or max(stats['direct_Williamson_covariance_error'],stats['purity_error'])>4e-9:
            raise ArithmeticError("Sampled covariance checks failed")
        return cov, sampled, stats


def analyze(phase, passage, geometry, solution, event, *, cutoff, points=192,
            nodes=64, step=.003, tolerance=2e-12, geometric=False):
    times=dict(bounce=0., **event['sample_times'])
    ham=phase.WorkHamiltonian(passage,geometry,solution,cutoff,points,geometric=geometric)
    ev=Evolution(ham,times,nodes,step,tolerance)
    covs,forms,stats = {},{},{}
    for width in WIDTHS:
        covs[width], forms[width], stats[str(width)] = ev.sample_state(width,nodes)
    end_forms={label:ev.forms(t,True) for label,t in times.items()}
    volumes={label:ev.volume(t) for label,t in times.items()}
    work=integrate_interval(lambda t:ev.forms(t)[1:],times['before'],times['after'],nodes//2)
    full_ledger=float(np.max(abs(end_forms['after']-end_forms['before']-work.sum(axis=0))))
    if full_ledger>3e-9:
        raise ArithmeticError("Full two-part work matrix failed to close")
    differences={}
    for width in (.1,.3):
        delta=covs[width]-covs[.2]
        energy={name:expectation(g,delta) for name,g in end_forms.items()}
        wz,wg=[expectation(g,delta) for g in work]
        residual=energy['after']-energy['before']-wz-wg
        penalties=dict(on_comparison=float(expectation(forms[width],-delta)),
                       on_reference=float(expectation(forms[.2],delta)))
        if min(penalties.values())<-1e-10 or abs(residual)>3e-10:
            raise ArithmeticError("Variational or relative work check failed")
        differences[str(width)]=dict(energy=energy,spatial_mean_density={
            name:E/volumes[name] for name,E in energy.items()},coupling_work=wz,metric_work=wg,
            energy_change=energy['after']-energy['before'],balance_residual=residual,
            objective_penalties=penalties)
    symplectic_error=max(np.max(abs(ev.state(t).T@ev.J@ev.state(t)-ev.J)) for t in times.values())
    if symplectic_error>3e-9:
        raise ArithmeticError("Canonical evolution failed")
    return dict(cutoff=cutoff, crossing=times['crossing'],volumes=volumes,state_statistics=stats,
                differences_from_width_0p2=differences, symplectic_error=float(symplectic_error),
                full_matrix_balance_error=full_ledger)


def round_control(phase, passage, geometry):
    """Compare canonical matrix minimization to exact nonzero de Sitter modes."""
    class RoundSolution:
        def sol(self, t):
            return np.zeros(6*geometry.n)

    times = dict(bounce=0., before=.395, crossing=.4, after=.405)
    ham = phase.WorkHamiltonian(passage,geometry,RoundSolution(),16,96,geometric=True)
    ev = Evolution(ham,times,64,.002,2e-13)
    h = geometry.branch.h
    L = ham.labels + 1.
    lam = L**2-1.
    norm = np.sqrt(2*L*lam)
    result = {}
    for width in WIDTHS:
        cov, _, stats = ev.sample_state(width,64)
        physical = ev.scaling @ cov @ ev.scaling.T
        # Independent exact Klein-Gordon modes with a^3 Wronskian equal to i.
        z, weights = roots_legendre(96)
        t = width*(z+1)/2
        eta = np.arctan(np.sinh(h*t))[:,None]
        c = np.cos(eta)
        phase_factor = np.exp(-1j*L*eta)
        u = h*(L*c+1j*np.sin(eta))*phase_factor/norm
        du = -1j*h**2*lam*c**2*phase_factor/norm
        om2 = lam*(h*c)**2
        weight = weights*np.exp(-1/(1-(t/width)**2))
        weight /= sum(weight)
        diagonal = weight @ ((abs(du)**2 + om2*abs(u)**2)/2)
        anomalous = weight @ ((du**2 + om2*u**2).real/2)
        squeeze = -.5*np.arctanh(anomalous/diagonal)
        qvar = (h*L/norm)**2*np.exp(2*squeeze)
        pvar = (lam/(h*norm))**2*np.exp(-2*squeeze)
        density = np.sqrt((diagonal-anomalous)*(diagonal+anomalous)).sum()/(2*np.pi**2)
        result[str(width)] = dict(maximum_full_covariance_error=float(np.max(abs(
            physical-np.diag(np.r_[qvar,pvar])))), sampled_density_error=float(abs(
            density-stats['own_sampled_mean_density'])))
    if max(x for row in result.values() for x in row.values())>5e-10:
        raise ArithmeticError("Round exact-mode check failed")
    return result


def changes(base,other):
    by_width={}
    for label,d in base['differences_from_width_0p2'].items():
        b=other['differences_from_width_0p2'][label]
        by_width[label]={k:b[k]-d[k] for k in ('coupling_work','metric_work','energy_change')}
        by_width[label]['mean_density_at_crossing'] = (b['spatial_mean_density']['crossing']
                                                     -d['spatial_mean_density']['crossing'])
        by_width[label]['mean_density_at_bounce'] = (b['spatial_mean_density']['bounce']
                                                   -d['spatial_mean_density']['bounce'])
    return by_width


def calculate(prior):
    spec=importlib.util.spec_from_file_location('frozen_phase_work',prior/'checks/phase_work.py')
    phase=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(phase)
    old=prior/'prior'
    passage=phase.load(old/'passage/passage.py')
    gravity=passage.dynamic.intrinsic.gravity_module(old/'gravity')
    branch=gravity.constraint.load_branches(old/'gravity/previous')['phi_m=-0.5']
    result=dict(branch='phi_m=-0.5',angular_ell=2,one_angular_real_block=True,
                widths=list(WIDTHS),reference_width=.2,
                weight="exp[-1/(1-(x/T)^2)] in Einstein geodesic proper time, normalized on (-T,T)",
                objective="expected original canonical Hamiltonian divided by actual whole-slice proper volume",
                finite_sector="reflection-even radial n=2,4,...,nmax in one real ell=2,m block",runs={})
    for amplitude in (.001,-.001):
        g,sol,event=passage.geometry_through_passage(gravity,branch,amplitude,64,160)
        runs={}
        for cutoff in (16,32,48,64,96,128,160):
            runs[str(cutoff)]=analyze(phase,passage,g,sol,event,cutoff=cutoff,points=max(192,cutoff+64),
                                      step=min(.003,.180/cutoff))
        geometry,solution,finer=passage.geometry_through_passage(gravity,branch,amplitude,96,256,
                                                                  tolerance=2e-13,max_step=.0015)
        control=analyze(phase,passage,geometry,solution,finer,cutoff=48,points=272,nodes=96,
                        tolerance=5e-13,step=.0015)
        peak=analyze(phase,passage,geometry,solution,finer,cutoff=160,points=288,nodes=96,
                     tolerance=5e-13,step=.00075)
        best_convergence=changes(runs['160'],peak)
        constant=analyze(phase,passage,g,sol,event,cutoff=48,geometric=True)
        convergence=changes(runs['48'],control)
        max_error=max(abs(v) for ref in (convergence,best_convergence) for row in ref.values() for v in row.values())
        if max_error>2e-9:
            raise ArithmeticError(f"Fixed-band resolution error {max_error}")
        for row in constant['differences_from_width_0p2'].values():
            if row['coupling_work'] != 0:
                raise ArithmeticError("Constant Z has a coupling signal")
        result['runs'][f'{amplitude:g}']=dict(crossing=event['crossing_x'],radial_cutoffs=runs,
            geometry_time_space_refinement_at_48=convergence,constant_Z_at_48=constant,
            geometry_time_space_refinement_at_160=best_convergence,
            cutoff_64_to_96=changes(runs['64'],runs['96']),
            cutoff_96_to_128=changes(runs['96'],runs['128']),
            cutoff_128_to_160=changes(runs['128'],runs['160']))
    result['round_closed_de_Sitter_control_n16']=round_control(phase,passage,g)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior',type=Path,default=HERE.parent/'prior')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    data=calculate(args.prior.resolve())
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
