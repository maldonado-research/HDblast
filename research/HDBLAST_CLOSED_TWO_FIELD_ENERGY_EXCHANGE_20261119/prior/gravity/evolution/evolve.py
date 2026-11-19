"""Local synchronous Einstein--scalar evolution on a complete closed S^3.

Metric: -dx^2 + exp(2L)dchi^2 + exp(2M)sin(chi)^2 dOmega^2.
L=M+(1-y^2)D with y=cos(chi) removes axis divisions analytically.
Metric variables are deviations from the closed de Sitter solution, and all
fields are represented by reflection-even Chebyshev U polynomials.  The
quadrature measure is normalized sin(chi)^2 dchi, on the actual S^3.
"""
import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import eval_chebyu, eval_gegenbauer

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('constraint_slice', HERE.parent/'initial/constraint.py')
constraint = importlib.util.module_from_spec(spec)
spec.loader.exec_module(constraint)
SAMPLE_TIMES = (0., .01, .02, .04, .06, .08, .10, .12)


def exp_remainder(z):
    """(exp(z)-1-z)/z^2, including its smooth value at zero."""
    z = np.asarray(z)
    small = abs(z) < .01
    out = np.empty_like(z)
    v = z[small]
    out[small] = .5+v*(1/6+v*(1/24+v*(1/120+v*(1/720+v/5040))))
    out[~small] = (np.expm1(z[~small])-z[~small])/z[~small]**2
    return out


class Evolution:
    def __init__(self, branch, maximum_degree=24, nodes=96, monitor_nodes=257):
        self.branch = branch
        self.degree = maximum_degree
        self.degrees = np.arange(0, maximum_degree+1, 2)
        self.n = len(self.degrees)
        self.chi, self.weight = constraint.quadrature(nodes)
        self.y = np.cos(self.chi)
        self.w = 1-self.y**2
        self.U, self.Uy, self.Uyy = self.matrices(self.y)
        self.project = self.U.T*self.weight
        chi_check, weight_check = constraint.quadrature(max(2*nodes+1, 193))
        self.monitor_chi = np.linspace(0., math.pi, monitor_nodes)
        self.monitor_all = self.matrices(np.cos(self.monitor_chi))
        self.monitor = self.monitor_all[0]
        self.check_y, self.check_weight = np.cos(chi_check), weight_check
        self.check_w = 1-self.check_y**2
        self.check = self.matrices(self.check_y)
        self.gram_error = float(np.max(abs(self.project@self.U-np.eye(self.n))))

    def matrices(self, y):
        y = np.asarray(y)[:, None]
        U = eval_chebyu(self.degrees, y)
        Uy, Uyy = np.zeros_like(U), np.zeros_like(U)
        Uy[:, 1:] = 2*eval_gegenbauer(self.degrees[1:]-1, 2, y)
        Uyy[:, 1:] = 8*eval_gegenbauer(self.degrees[1:]-2, 3, y)
        return U, Uy, Uyy

    def background(self, x):
        h = self.branch.h
        a = math.cosh(h*x)/h
        H = h*math.tanh(h*x)
        return a, H, 1/a**2

    def fields(self, state, matrices=None):
        U, Uy, Uyy = (self.U, self.Uy, self.Uyy) if matrices is None else matrices
        m, d, b, e, q, pi = state.reshape(6, self.n)
        return ((U@m, U@d, U@b, U@e, U@q, U@pi),
                (Uy@m, Uy@d, Uy@b, Uy@q), (Uyy@m, Uyy@q))

    def geometry(self, x, fields, y, w):
        (m,d,b,e,q,pi), (my,dy,by,qy), (myy,qyy) = fields
        a,H,inv_bg = self.background(x)
        ell = m+w*d
        inv = inv_bg*np.exp(-2*ell)
        ly = my+w*dy-2*y*d
        F = y-w*my
        # K=-Rss/(R*exp(-2L)); J=(exp(2wD)-F^2)/w.
        dK = 2*y*my-w*myy+(my-ly)*F
        dJ = 2*d*np.where(abs(2*w*d)<1e-7, 1+w*d+2*w*w*d*d/3,
                         np.expm1(2*w*d)/np.where(2*w*d != 0, 2*w*d, 1)) + 2*y*my-w*my**2
        delta_ric_t = 2*inv_bg*np.expm1(-2*ell)+inv*(dJ+dK)
        delta_R3 = 6*inv_bg*np.expm1(-2*ell)+2*inv*(dJ+2*dK)
        return a,H,inv,ly,F,dK,dJ,delta_ric_t,delta_R3

    def equation(self, x, state):
        fields = self.fields(state)
        (m,d,b,e,q,pi), (my,dy,by,qy), (myy,qyy) = fields
        a,H,inv,ly,F,dK,dJ,delta_ric_t,delta_R3 = self.geometry(x,fields,self.y,self.w)
        u,force = self.branch.potential(q)
        theta = 3*(H+b)+self.w*e
        rhs_b = -delta_ric_t-6*H*b-3*b*b-self.w*e*(H+b)+3*(u-self.branch.u_top)
        G = 4*d*d*exp_remainder(2*self.w*d)
        rhs_e = inv*(myy-ly*my+2*d+self.y*dy+G+qy*qy)-theta*e
        lap = inv*(self.w*qyy-3*self.y*qy+self.w*(2*my-ly)*qy)
        rhs_pi = lap-theta*pi+force
        z = state.reshape(6,self.n)
        return np.r_[z[2],z[3],self.project@rhs_b,self.project@rhs_e,z[5],self.project@rhs_pi]

    def record(self, x, state, initial_mean=0.):
        fields = self.fields(state,self.check)
        (m,d,b,e,q,pi), (my,dy,by,qy), (myy,qyy) = fields
        y,w,weights = self.check_y,self.check_w,self.check_weight
        a,H,inv,ly,F,dK,dJ,dric,dR3 = self.geometry(x,fields,y,w)
        u,_ = self.branch.potential(q)
        du = u-self.branch.u_top
        logweight = 3*m+w*d
        volume_shift = float(weights@np.expm1(logweight))
        volume_ratio = 1+volume_shift
        norm = weights*np.exp(logweight)/volume_ratio
        mean = float(norm@q)
        variance = float(norm@((q-mean)**2))
        theta_delta = 3*b+w*e
        mean_theta_delta = float(norm@theta_delta)
        theta_var = float(norm@((theta_delta-mean_theta_delta)**2))
        acceleration_shift = float(norm@(du-pi*pi/3-2*(w*e)**2/9)+2*theta_var/9)
        ham = dR3+12*H*b+6*b*b+4*w*(H+b)*e-pi*pi-inv*w*qy*qy-6*du
        mom_regular = 2*(by+F*e)+pi*qy
        mom = np.sqrt(inv*w)*mom_regular
        # Independent trace evolution identity, projected RHS versus Raychaudhuri.
        rhs = self.equation(x,state).reshape(6,self.n)
        Bdot = self.check[0]@rhs[2]+self.branch.h**2/math.cosh(self.branch.h*x)**2
        Edot = self.check[0]@rhs[3]
        trace_dot = 3*Bdot+w*Edot
        theta = 3*H+theta_delta
        alpha = H+b+w*e
        beta = H+b
        ray = trace_dot+alpha*alpha+2*beta*beta+pi*pi-3*u
        monitored = self.monitor@state.reshape(6,self.n).T
        mm,dd,bb,ee,qmonitor,pimonitor = monitored.T
        monitor_w = np.sin(self.monitor_chi)**2
        monitor_inv = np.exp(-2*(mm+monitor_w*dd))/a**2
        monitor_qy = self.monitor_all[1]@state.reshape(6,self.n)[4]
        monitor_u,_ = self.branch.potential(qmonitor)
        kinetic = float(norm@(pi*pi/6))/self.branch.u_top
        gradient = float(norm@(inv*w*qy*qy/6))/self.branch.u_top
        potential = float(norm@du)/self.branch.u_top
        return dict(x=float(x),mean=mean,mean_change=mean-initial_mean,
                    round_mean=float(state.reshape(6,self.n)[4,0]), variance=variance,
                    stddev=math.sqrt(max(variance,0.)), maximum_field=float(np.max(abs(qmonitor))),
                    volume=2*math.pi**2*a**3*volume_ratio,
                    volume_over_background=volume_ratio, volume_fractional_difference=volume_shift,
                    volume_radius_over_background_minus_one=float(np.expm1(math.log1p(volume_shift)/3)),
                    H_volume=H+mean_theta_delta/3, H_volume_minus_background=mean_theta_delta/3,
                    volume_radius_acceleration_over_radius=self.branch.u_top+acceleration_shift,
                    acceleration_difference_from_background=acceleration_shift,
                    maximum_local_expansion_difference=float(np.max(3*bb+monitor_w*ee)/3),
                    maximum_abs_metric_deviation=float(max(np.max(abs(mm)),np.max(abs(mm+monitor_w*dd)))),
                    maximum_abs_shear=float(np.max(abs(monitor_w*ee))),
                    minimum_spatial_volume_ratio=float(np.min(np.exp(3*mm+monitor_w*dd))),
                    largest_local_energy_terms_over_Vstar=dict(
                        kinetic=float(np.max(pimonitor*pimonitor)/(6*self.branch.u_top)),
                        gradient=float(np.max(monitor_inv*monitor_w*monitor_qy*monitor_qy)/(6*self.branch.u_top)),
                        abs_potential_shift=float(np.max(abs(monitor_u-self.branch.u_top))/self.branch.u_top)),
                    energy_over_Vstar=dict(kinetic=kinetic,gradient=gradient,potential_minus_top=potential,signed_sum=kinetic+gradient+potential),
                    constraints=dict(max_abs_hamiltonian_over_6h2=float(np.max(abs(ham))/(6*self.branch.u_top)),
                                     rms_hamiltonian_over_6h2=float(math.sqrt(norm@(ham*ham))/(6*self.branch.u_top)),
                                     max_abs_physical_momentum_over_h2=float(np.max(abs(mom))/self.branch.u_top),
                                     max_abs_regular_momentum_over_h=float(np.max(abs(mom_regular))/self.branch.h),
                                     max_abs_raychaudhuri_over_3h2=float(np.max(abs(ray))/(3*self.branch.u_top))))

    def initial(self, amplitude, offset):
        coeff = constraint.solve_slice(self.branch,amplitude,max_degree=max(24,self.degree),nodes=max(128,len(self.chi)),offset=offset)
        psi = constraint.evaluate(coeff,self.chi)
        state = np.zeros((6,self.n))
        state[0] = self.project@(2*np.log(psi))
        state[4,0], state[4,1] = offset,amplitude
        return state.ravel()

    def evolve(self, amplitude, offset=0., endpoint=.12, field_cap=None, tolerance=2e-11, max_step=.006):
        initial = self.initial(amplitude,offset)
        initial_mean = self.record(0.,initial)['mean']
        def stop(x,z):
            return field_cap-np.max(abs(self.monitor@z.reshape(6,self.n)[4]))
        stop.terminal=True
        stop.direction=-1
        sol = solve_ivp(self.equation,(0.,endpoint),initial,method='DOP853',rtol=tolerance,atol=2e-15,
                        max_step=max_step,dense_output=True,events=stop if field_cap is not None else None)
        if not sol.success:
            raise ArithmeticError(sol.message)
        samples = {f'{x:g}':self.record(x,sol.sol(x),initial_mean) for x in SAMPLE_TIMES if x<=sol.t[-1]+1e-13}
        at_stop = self.record(sol.t[-1],sol.y[:,-1],initial_mean)
        audit = [self.record(x,sol.sol(x),initial_mean) for x in np.linspace(0.,sol.t[-1],17)]
        max_constraints = {key:max(r['constraints'][key] for r in audit) for key in audit[0]['constraints']}
        return dict(amplitude=amplitude,offset=offset,maximum_degree=self.degree,quadrature_nodes=len(self.chi),
                    stop='peak_field' if sol.t_events is not None and len(sol.t_events[0]) else 'endpoint',
                    samples=samples,at_probe=samples.get('0.12'),at_stop=at_stop,
                    max_sampled_constraints=max_constraints,
                    number_rhs_evaluations=sol.nfev,success=sol.success),sol


def fixed_record(branch, amplitude, offset, times, degree=24, nodes=96):
    evol = Evolution(branch,degree,nodes)
    z = np.zeros((2,evol.n)); z[0,0]=offset;z[0,1]=amplitude
    def rhs(x,state):
        c,v=state.reshape(2,evol.n)
        q,pi=evol.U@c,evol.U@v
        _,force=branch.potential(q)
        a,H,_=evol.background(x)
        return np.r_[v,evol.project@force-evol.degrees*(evol.degrees+2)*c/a**2-3*H*v]
    sol=solve_ivp(rhs,(0.,max(times)),z.ravel(),method='DOP853',rtol=2e-12,atol=2e-15,max_step=.006,dense_output=True)
    if not sol.success: raise ArithmeticError(sol.message)
    records={}
    for x in times:
        c,v=sol.sol(x).reshape(2,evol.n)
        records[f'{x:g}']=dict(x=x,mean=float(c[0]),mean_change=float(c[0]-offset),stddev=float(np.linalg.norm(c[1:])),
                               maximum_field=float(np.max(abs(evol.monitor@c))))
    return records


def validate_homogeneous(branch):
    """Independent closed-FRW ODE check at a nonstationary homogeneous field."""
    Q = .001
    solver = Evolution(branch,24,96)
    run,_ = solver.evolve(0.,offset=Q)
    u0,_ = branch.potential(Q)
    def rhs(x,z):
        q,pi,H,loga = z
        u,force = branch.potential(q)
        return [pi,force-3*H*pi,u-pi*pi/3-H*H,H]
    sol=solve_ivp(rhs,(0.,.12),[Q,0.,0.,-.5*math.log(u0)],method='DOP853',
                  rtol=2e-13,atol=2e-15,max_step=.003)
    q,pi,H,loga=sol.y[:,-1]
    last=run['at_probe']
    return dict(offset=Q,mean_difference_from_FRW=last['mean']-q,
                H_volume_difference_from_FRW=last['H_volume']-H,
                volume_relative_difference_from_FRW=last['volume']/(2*math.pi**2*math.exp(3*loga))-1,
                maximum_abs_shear=last['maximum_abs_shear'],constraints=run['max_sampled_constraints'])


def controlled_extensions(branch, primary):
    """Stop at the previous study's early field caps; compare at equal proper time."""
    cap=.04 if branch.fitted.pm>0 else .08
    runs={}
    solutions={}
    for A in (.001,-.001,.0005,-.0005):
        run,sol=primary.evolve(A,endpoint=.65,field_cap=cap)
        if run['stop']!='peak_field':
            raise ArithmeticError('Extension missed early field cap')
        runs[f'{A:g}'],solutions[f'{A:g}']=run,sol
    common_time=min(r['at_stop']['x'] for r in runs.values())
    refine=Evolution(branch,24,96)
    refinements=[]
    temporal=[]
    for textA,run in runs.items():
        A=float(textA)
        initial_mean=run['samples']['0']['mean']
        run['at_common_time']=primary.record(common_time,solutions[textA].sol(common_time),initial_mean)
        times=sorted(set([float(t) for t in run['samples']]+[common_time,run['at_stop']['x']]))
        fixed=fixed_record(branch,A,0.,times)
        run['fixed_geometry_samples']=fixed
        run['fixed_geometry_at_stop']=fixed[f"{run['at_stop']['x']:g}"]
        run['fixed_geometry_at_common_time']=fixed[f'{common_time:g}']
        check,_=refine.evolve(A,endpoint=.65,field_cap=cap)
        delta={name:check['at_stop'][name]-run['at_stop'][name]
               for name in ('x','mean','mean_change','stddev','volume_fractional_difference',
                            'H_volume_minus_background','acceleration_difference_from_background',
                            'maximum_abs_metric_deviation','maximum_abs_shear')}
        refinements.append(dict(amplitude=A,maximum_degree=24,quadrature_nodes=96,
                                stop_differences_from_primary=delta,
                                max_sampled_constraints=check['max_sampled_constraints']))
        check,_=primary.evolve(A,endpoint=.65,field_cap=cap,tolerance=2e-12,max_step=.003)
        temporal.append(dict(amplitude=A,relative_tolerance=2e-12,maximum_step=.003,
                             stop_differences_from_primary={name:check['at_stop'][name]-run['at_stop'][name]
                              for name in ('x','mean','mean_change','stddev','volume_fractional_difference',
                                           'H_volume_minus_background','acceleration_difference_from_background')}))
    return dict(peak_field_cap=cap,common_comparison_time=common_time,runs=runs,
                spatial_refinements=refinements,temporal_refinements=temporal,
                maximum_verified_time=max(r['at_stop']['x'] for r in runs.values()))


def calculate(previous, endpoint=.12, extensions=True):
    data=dict(formulation='closed S3, unit-lapse zero-shift Gaussian proper time, time-symmetric constraint initial data',
              primary_maximum_degree=32,primary_quadrature_nodes=128,probe_x=.12,branches={})
    for key,branch in constraint.load_branches(previous).items():
        primary=Evolution(branch,32,128)
        runs={}
        for A in (0.,.001,-.001,.0005,-.0005):
            run,_=primary.evolve(A,endpoint=endpoint)
            run['fixed_geometry_samples']=fixed_record(branch,A,0.,[float(t) for t in run['samples']])
            runs[f'{A:g}']=run
        centered={}
        for A in (.001,-.001):
            offset=constraint.centered_offset(branch,A)
            run,_=primary.evolve(A,offset=offset,endpoint=endpoint)
            run['fixed_geometry_samples']=fixed_record(branch,A,offset,[float(t) for t in run['samples']])
            centered[f'{A:g}']=run
        refinements=[]
        for degree,nodes in ((16,64),(24,96)):
            other=Evolution(branch,degree,nodes)
            for A in (.001,-.001,.0005,-.0005):
                run,_=other.evolve(A,endpoint=endpoint)
                base=runs[f'{A:g}']
                deltas={observable:run['at_probe'][observable]-base['at_probe'][observable]
                        for observable in ('mean','mean_change','stddev','volume_fractional_difference','H_volume_minus_background','acceleration_difference_from_background')}
                refinements.append(dict(amplitude=A,maximum_degree=degree,quadrature_nodes=nodes,samples=run['samples'],
                                        at_stop=run['at_stop'],probe_differences_from_primary=deltas,
                                        max_sampled_constraints=run['max_sampled_constraints']))
        pairs={}
        for A in (.001,.0005):
            pos,neg=runs[f'{A:g}']['at_probe'],runs[f'{-A:g}']['at_probe']
            initial_pos,initial_neg=runs[f'{A:g}']['samples']['0'],runs[f'{-A:g}']['samples']['0']
            fixed_pos=runs[f'{A:g}']['fixed_geometry_samples']['0.12'];fixed_neg=runs[f'{-A:g}']['fixed_geometry_samples']['0.12']
            pairs[f'{A:g}']=dict(even_mean=(pos['mean']+neg['mean'])/2,
                                even_mean_change=(pos['mean_change']+neg['mean_change'])/2,
                                even_initial_mean=(initial_pos['mean']+initial_neg['mean'])/2,
                                odd_mean=(pos['mean']-neg['mean'])/2,
                                even_mean_change_over_A_squared=(pos['mean_change']+neg['mean_change'])/(2*A*A),
                                fixed_even_mean_change=(fixed_pos['mean_change']+fixed_neg['mean_change'])/2,
                                even_acceleration_difference=(pos['acceleration_difference_from_background']+neg['acceleration_difference_from_background'])/2)
        temporal=[]
        for A in (.001,-.001):
            run,_=primary.evolve(A,endpoint=endpoint,tolerance=2e-12,max_step=.003)
            base=runs[f'{A:g}']['at_probe']
            temporal.append(dict(amplitude=A,relative_tolerance=2e-12,maximum_step=.003,
                                 probe_differences_from_primary={name:run['at_probe'][name]-base[name]
                                  for name in ('mean','mean_change','stddev','volume_fractional_difference',
                                               'H_volume_minus_background','acceleration_difference_from_background')}))
        entry=dict(h=branch.h,mu_squared=branch.mu2,gram_error=primary.gram_error,runs=runs,
                   proper_mean_centered_runs=centered,amplitude_pairs_at_probe=pairs,
                   spatial_refinements=refinements,temporal_refinements=temporal,
                   homogeneous_FRW_check=validate_homogeneous(branch))
        if extensions:
            entry['controlled_extensions']=controlled_extensions(branch,primary)
        data['branches'][key]=entry
    data['verification']=verify(data)
    return data


def verify(data):
    """Fail the replay if the measured numerical checks cease to be controlled."""
    runs=[]
    spatial=[]
    temporal=[]
    centered=[]
    frw=[]
    extension_runs=[]
    for entry in data['branches'].values():
        runs+=list(entry['runs'].values())+list(entry['proper_mean_centered_runs'].values())
        centered += [abs(r['samples']['0']['mean']) for r in entry['proper_mean_centered_runs'].values()]
        spatial += [r['probe_differences_from_primary'] for r in entry['spatial_refinements']]
        temporal += [r['probe_differences_from_primary'] for r in entry['temporal_refinements']]
        frw.append(entry['homogeneous_FRW_check'])
        zero=entry['runs']['0']['at_stop']
        if any(zero[n]!=0. for n in ('mean','stddev','volume_fractional_difference','maximum_abs_metric_deviation')):
            raise ArithmeticError('Stationary de Sitter control moved')
        if 'controlled_extensions' in entry:
            ex=entry['controlled_extensions']
            extension_runs+=list(ex['runs'].values())
            spatial += [r['stop_differences_from_primary'] for r in ex['spatial_refinements']]
            temporal += [r['stop_differences_from_primary'] for r in ex['temporal_refinements']]
    runs+=extension_runs
    keys=runs[0]['max_sampled_constraints'].keys()
    constraints={n:max(r['max_sampled_constraints'][n] for r in runs) for n in keys}
    spatial_max={n:max(abs(r.get(n,0.)) for r in spatial) for n in set().union(*(r.keys() for r in spatial))}
    temporal_max={n:max(abs(r.get(n,0.)) for r in temporal) for n in set().union(*(r.keys() for r in temporal))}
    if max(constraints.values())>1e-7:
        raise ArithmeticError('Einstein constraint/evolution audit failed')
    if max(spatial_max.values())>1e-9 or max(temporal_max.values())>1e-9:
        raise ArithmeticError('Spatial/temporal refinement failed')
    if max(centered)>1e-14:
        raise ArithmeticError('Proper-mean-centered initial condition failed')
    if max(abs(r[n]) for r in frw for n in ('mean_difference_from_FRW','H_volume_difference_from_FRW','volume_relative_difference_from_FRW'))>1e-10:
        raise ArithmeticError('Independent homogeneous FRW comparison failed')
    return dict(passed=True,max_sampled_constraints=constraints,
                maximum_spatial_refinement_differences=spatial_max,
                maximum_temporal_refinement_differences=temporal_max,
                maximum_abs_centered_initial_mean=max(centered),
                maximum_verified_time=max(r['at_stop']['x'] for r in runs),
                maximum_final_abs_metric_deviation=max(r['at_stop']['maximum_abs_metric_deviation'] for r in runs),
                minimum_final_spatial_volume_ratio=min(r['at_stop']['minimum_spatial_volume_ratio'] for r in runs))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--previous',type=Path,default=Path(__file__).resolve().parents[1]/'previous')
    parser.add_argument('--output',type=Path,default=HERE/'results.json')
    parser.add_argument('--endpoint',type=float,default=.12)
    parser.add_argument('--skip-extensions',action='store_true')
    args=parser.parse_args()
    if args.endpoint<.12:
        parser.error('--endpoint must include the fixed comparison time x=0.12')
    result=calculate(args.previous,args.endpoint,not args.skip_extensions)
    args.output.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
