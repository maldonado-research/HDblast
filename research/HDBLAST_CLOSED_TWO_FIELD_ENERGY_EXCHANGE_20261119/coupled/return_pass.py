"""Whole-sphere exchange until the same polar geodesic crosses the minimum again."""
import argparse
import importlib.util
import json
from pathlib import Path

path = Path(__file__).resolve().with_name("coupled.py")
spec = importlib.util.spec_from_file_location("two_fields", path)
coupled = importlib.util.module_from_spec(spec)
spec.loader.exec_module(coupled)

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import roots_legendre


def run(model, seed, phase, step, tolerance):
    e = model.e
    y0 = model.initial(.001, seed, phase)
    qmin = model.branch.fitted.tmin-model.branch.top
    def polar(t, z, i):
        return float(e.monitor[0] @ z.reshape(8,e.n)[i])
    def return_event(t, z):
        return polar(t,z,4)-qmin
    return_event.direction = -1
    return_event.terminal = True
    def entry_event(t, z):
        return polar(t,z,4)-qmin
    entry_event.direction = 1
    def turn_event(t, z):
        return polar(t,z,5)
    turn_event.direction = -1
    # The initial parent velocity vanishes, so select the turning root AFTER entry.
    sol = solve_ivp(model.equation, (0.,.7), y0, method="DOP853", events=[entry_event,turn_event,return_event],
                    max_step=step,rtol=tolerance,atol=2e-15,dense_output=True)
    if not sol.success or len(sol.t_events[0])!=1 or len(sol.t_events[2])!=1:
        raise ArithmeticError("Did not reach both polar crossings")
    entry, back = sol.t_events[0][0],sol.t_events[2][0]
    turns = [x for x in sol.t_events[1] if entry<x<back]
    if len(turns)!=1:
        raise ArithmeticError("Expected one polar reversal between crossings")
    turn = turns[0]
    times = dict(launch=0., first_crossing=entry, turning_point=turn, return_crossing=back)
    snapshots = {}
    for label,t in times.items():
        z=sol.sol(t)
        r=model.diagnostics(t,z)
        _,_,b,shear,q,pi,_,_=(e.monitor @ z.reshape(8,e.n).T).T
        H=e.background(t)[1]
        local=H+b+(np.sin(e.monitor_chi)**2)*shear/3
        snapshots[label]=dict(t=t,polar_q=q[0],polar_parent_velocity=pi[0],polar_H=local[0],
                              minimum_local_H=min(local),H_volume=r[7],proper_volume_over_reference=r[6],
                              E_parent=r[0],E_daughter=r[1],daughter_transfer_rate=r[4],
                              daughter_fraction_of_total=r[1]/(r[0]+r[1]))
    def ledger(lo,hi,order):
        x,w=roots_legendre(order)
        t=lo+.5*(hi-lo)*(x+1)
        rates=np.array([model.diagnostics(v,sol.sol(v))[2:5] for v in t])
        wp,wd,transfer=.5*(hi-lo)*(w @ rates)
        dt=model.diagnostics(hi,sol.sol(hi))[:2]-model.diagnostics(lo,sol.sol(lo))[:2]
        return dict(parent_energy_change=dt[0],daughter_energy_change=dt[1],parent_pressure_work=wp,
                    daughter_pressure_work=wd,transfer_into_daughter=transfer,
                    daughter_residual=dt[1]+wd-transfer,parent_residual=dt[0]+wp+transfer,
                    total_residual=sum(dt)+wp+wd)
    ranges=dict(launch_to_first=(0.,entry),first_to_turn=(entry,turn),turn_to_return=(turn,back),
                first_to_return=(entry,back),launch_to_return=(0.,back))
    ledger_rows={key:ledger(lo,hi,80) for key,(lo,hi) in ranges.items()}
    # Piecewise and one-pass time integrals give separate closure and quadrature tests.
    names=('launch_to_first','first_to_turn','turn_to_return')
    diagscales=("transfer_into_daughter","parent_pressure_work","daughter_pressure_work")
    piece_error=max(abs(sum(ledger_rows[n][k] for n in names)-ledger_rows['launch_to_return'][k])
                    for k in diagscales)
    finer=ledger(0.,back,112)
    quadrature_error=max(abs(finer[k]-ledger_rows['launch_to_return'][k]) for k in diagscales)
    peak=np.zeros(3)
    domain_margin=np.inf
    charge0=model.diagnostics(0.,y0)[5]
    charge_drift=0.
    # Apply the off-collocation constraint and domain tests on every accepted step,
    # and at the three physical events. All trial RK stages already check the full
    # interpolated-domain enclosure via the frozen TwoFields.equation.
    for t in [*sol.t,entry,turn,back]:
        z=sol.sol(t)
        peak=np.maximum(peak,model.diagnostics(t,z,True))
        c=z.reshape(8,e.n)[4]
        q=e.monitor @ c
        err=abs(c) @ model.angular_error_factor
        domain_margin=min(domain_margin,q.min()-err-model.domain[0],model.domain[1]-q.max()-err)
        charge_drift=max(charge_drift,abs(model.diagnostics(t,z)[5]-charge0))
    if domain_margin<=0 or max(peak)>5e-8 or quadrature_error>2e-9 or piece_error>2e-9:
        raise ArithmeticError(f"Uncontrolled return at seed {seed}: {peak}, {domain_margin}")
    residual_max=max(abs(row[key]) for row in ledger_rows.values()
                     for key in ('daughter_residual','parent_residual','total_residual'))
    if residual_max>2e-9 or charge_drift>1e-10:
        raise ArithmeticError("Energy balance or shift charge failed at return")
    if not model.matter and any(row["transfer_into_daughter"]!=0 for row in ledger_rows.values()):
        raise ArithmeticError("Nonzero transfer with constant kinetic coupling")
    initial=snapshots['launch']['E_daughter']
    if seed==0 and any(r['E_daughter']!=0 for r in snapshots.values()):
        raise ArithmeticError("Daughter excited in the classical zero solution")
    summary={key:100*row["transfer_into_daughter"]/initial for key,row in ledger_rows.items()} if seed else {}
    return dict(daughter_seed=seed,launch_phase=phase,matter=model.matter,maximum_degree=e.degree,
                maximum_step=step,snapshots=snapshots,balances=ledger_rows,
                transfer_percent_of_initial_daughter_energy=summary,
                minimum_certified_domain_margin=domain_margin,constraint_peaks=dict(hamiltonian=peak[0],momentum=peak[1],raychaudhuri=peak[2]),
                max_shift_charge_drift=charge_drift,interval_composition_error=piece_error,temporal_quadrature_error=quadrature_error)


def calculate(prior):
    module=coupled.load(prior)
    branch=module.constraint.load_branches(prior/"gravity/previous")["phi_m=-0.5"]
    result=dict(branch="phi_m=-0.5",parent_amplitude=.001,
                boundary="first return crossing on the same polar normal geodesic",
                normalization="integrated energies divided by stationary round launch volume times u_top",
                runs={},refinement={})
    for method,deg,nodes,step,tol in [('primary',128,320,.003,2e-12),('space',160,384,.003,2e-12),('time',128,320,.0015,4e-13)]:
        models={matter:coupled.TwoFields(module,branch,deg,nodes,matter) for matter in (True,False)}
        for phase in ('velocity','displacement'):
            for seed in (.01,.005,0.):
                for matter,model in models.items():
                    if not seed and (phase!='velocity' or not matter):
                        continue
                    key=f"{phase}/{seed:g}/"+('matter' if matter else 'constant')
                    case=run(model,seed,phase,step,tol)
                    if method=='primary':
                        result['runs'][key]=case
                    else:
                        ref=result['runs'][key]
                        event_diff=max(abs(case['snapshots'][k]['t']-r['t']) for k,r in ref['snapshots'].items())
                        rate_diff=max(abs(case['balances'][k]['transfer_into_daughter']-r['transfer_into_daughter']) for k,r in ref['balances'].items())
                        H_diff=max(abs(case['snapshots'][k]['H_volume']-r['H_volume']) for k,r in ref['snapshots'].items())
                        if event_diff>2e-8 or rate_diff>2e-12 or H_diff>2e-8:
                            raise ArithmeticError("Return refinement failed")
                        result['refinement'].setdefault(key,{})[method]=dict(maximum_event_time_difference=event_diff,
                            maximum_transfer_difference=rate_diff,maximum_H_volume_difference=H_diff,
                            constraint_peaks=case['constraint_peaks'])
    comp={}
    for phase in ('velocity','displacement'):
        for seed in (.01,.005):
            key=f"{phase}/{seed:g}"
            matter=result['runs'][key+'/matter']['snapshots']
            control=result['runs'][key+'/constant']['snapshots']
            unseeded=result['runs']['velocity/0/matter']['snapshots']
            comp[key]={k:dict(matter_minus_constant_time=row['t']-control[k]['t'],
                              matter_minus_unseeded_time=row['t']-unseeded[k]['t'])
                       for k,row in matter.items() if k!='launch'}
    result['paired_timings']=comp
    result['passed']=True
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior',type=Path,default=Path(__file__).resolve().parents[1]/'prior')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    data=calculate(args.prior)
    args.output.write_text(json.dumps(data,indent=2,allow_nan=False)+'\n')
